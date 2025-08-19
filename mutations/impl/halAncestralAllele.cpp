#include "hal.h"
#include "halCLParser.h"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace hal;

static void initParser(CLParser &parser) {
    parser.addArgument("halFile", "input hal file");
    parser.addArgument("refGenome", "reference genome name");
    parser.addArgument("targetGenome", "target (ancestral) genome name");
    parser.addArgument("positionsFile", "bed/gff file with reference coordinates");
    parser.addArgument("outputFile", "output tab-delimited file");
}

static bool parseLine(const string &line, string &chrom, hal_index_t &start, hal_index_t &end) {
    if (line.empty() || line[0] == '#') {
        return false;
    }
    stringstream ss(line);
    ss >> chrom >> start >> end;
    return !ss.fail();
}

int main(int argc, char **argv) {
    CLParser parser;
    initParser(parser);

    string halPath, refName, tgtName, posPath, outPath;
    try {
        parser.parseOptions(argc, argv);
        halPath = parser.getArgument<string>("halFile");
        refName = parser.getArgument<string>("refGenome");
        tgtName = parser.getArgument<string>("targetGenome");
        posPath = parser.getArgument<string>("positionsFile");
        outPath = parser.getArgument<string>("outputFile");
    } catch (exception &e) {
        cerr << e.what() << endl;
        parser.printUsage(cerr);
        return 1;
    }

    try {
        AlignmentConstPtr alignment(openHalAlignment(halPath, &parser));
        const Genome *refGenome = alignment->openGenome(refName);
        if (refGenome == NULL) {
            throw hal_exception("Reference genome " + refName + " not found");
        }
        const Genome *tgtGenome = alignment->openGenome(tgtName);
        if (tgtGenome == NULL) {
            throw hal_exception("Target genome " + tgtName + " not found");
        }

        ifstream posStream(posPath.c_str());
        if (!posStream) {
            throw hal_exception("Unable to open positions file: " + posPath);
        }
        ofstream outStream(outPath.c_str());
        if (!outStream) {
            throw hal_exception("Unable to open output file: " + outPath);
        }

        set<const Genome *> targets;
        targets.insert(tgtGenome);

        string line;
        while (getline(posStream, line)) {
            string chrom;
            hal_index_t start = 0, end = 0;
            if (!parseLine(line, chrom, start, end)) {
                continue;
            }
            const Sequence *refSeq = refGenome->getSequence(chrom);
            if (refSeq == NULL) {
                continue;
            }
            hal_index_t absStart = refSeq->getStartPosition() + start;
            DnaIteratorPtr refIt = refGenome->getDnaIterator(absStart);
            char refBase = toupper(refIt->getBase());

            ColumnIteratorPtr colIt = refGenome->getColumnIterator(&targets, 0, absStart, absStart, false, false);
            const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
            vector<char> bases;
            for (ColumnIterator::ColumnMap::const_iterator cmIt = colMap->begin(); cmIt != colMap->end(); ++cmIt) {
                const Sequence *seq = cmIt->first;
                if (seq->getGenome() != tgtGenome) {
                    continue;
                }
                ColumnIterator::DNASet *dnaSet = cmIt->second;
                for (size_t i = 0; i < dnaSet->size(); ++i) {
                    char b = toupper(dnaSet->at(i)->getBase());
                    bases.push_back(b);
                }
            }

            string ancAllele = "N";
            string evidence = "Missing";
            if (!bases.empty()) {
                map<char, int> counts;
                int total = 0;
                for (size_t i = 0; i < bases.size(); ++i) {
                    char b = bases[i];
                    if (b == 'N' || b == '-') {
                        continue;
                    }
                    counts[b]++;
                    total++;
                }
                if (total == 0) {
                    ancAllele = "N";
                    evidence = "Missing";
                } else if (total == 1) {
                    ancAllele = string(1, counts.begin()->first);
                    evidence = "Direct";
                } else {
                    // find majority and detect ties
                    char maj = 'N';
                    int majCount = 0;
                    bool tie = false;
                    for (map<char, int>::const_iterator it = counts.begin(); it != counts.end(); ++it) {
                        if (it->second > majCount) {
                            maj = it->first;
                            majCount = it->second;
                            tie = false;
                        } else if (it->second == majCount) {
                            tie = true;
                        }
                    }
                    std::ostringstream info;
                    for (map<char, int>::const_iterator it = counts.begin(); it != counts.end(); ++it) {
                        if (it != counts.begin()) {
                            info << ",";
                        }
                        info << it->first << "=" << it->second;
                    }
                    if (tie) {
                        ancAllele = "N";
                        evidence = string("Paralogous:") + info.str();
                    } else {
                        ancAllele = string(1, maj);
                        evidence = string("MajorityVote:") + info.str();
                    }
                }
            }

            outStream << chrom << '\t' << start << '\t' << end << '\t' << refBase << '\t' << tgtName << '\t' << ancAllele
                      << '\t' << evidence << '\n';
        }
    } catch (hal_exception &e) {
        cerr << "hal exception: " << e.what() << endl;
        return 1;
    }
    return 0;
}

