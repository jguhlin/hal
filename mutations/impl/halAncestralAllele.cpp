#include "hal.h"
#include "halCLParser.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace hal;

struct Position {
    string chrom;
    hal_index_t start;
    hal_index_t end;
    size_t originalOrder;
    
    bool operator<(const Position& other) const {
        if (chrom != other.chrom) {
            return chrom < other.chrom;
        }
        return start < other.start;
    }
};

static void initParser(CLParser &parser) {
    parser.addArgument("halFile", "input hal file");
    parser.addArgument("refGenome", "reference genome name");
    parser.addArgument("targetGenome", "target genome name (or comma-separated list of genomes to try in order)");
    parser.addArgument("positionsFile", "bed/gff file with reference coordinates");
    parser.addArgument("outputFile", "output tab-delimited file");
    parser.addOptionFlag("noSort", "disable position sorting optimization (preserves input order)", false);
    parser.addOption("progress", "report progress every N positions processed", 0);
}

static bool parseLine(const string &line, string &chrom, hal_index_t &start, hal_index_t &end) {
    if (line.empty() || line[0] == '#') {
        return false;
    }
    stringstream ss(line);
    ss >> chrom >> start >> end;
    return !ss.fail();
}

static vector<string> splitAncestors(const string &ancestorList) {
    vector<string> ancestors;
    stringstream ss(ancestorList);
    string ancestor;
    while (getline(ss, ancestor, ',')) {
        // Remove whitespace
        ancestor.erase(0, ancestor.find_first_not_of(" \t"));
        ancestor.erase(ancestor.find_last_not_of(" \t") + 1);
        if (!ancestor.empty()) {
            ancestors.push_back(ancestor);
        }
    }
    return ancestors;
}

static bool findParalogsInGenome(const Genome *refGenome, hal_index_t absStart, 
                                const Genome *targetGenome, vector<char>& bases) {
    try {
        // Safety check: validate inputs
        if (refGenome == NULL || targetGenome == NULL || absStart < 0) {
            return false;
        }
        
        // Note: Within-species search only when refGenome == targetGenome
        bool isWithinSpecies = (refGenome == targetGenome);
        
        // Safer approach: Use ColumnIterator with duplications enabled
        // This avoids the complex segment mapping that was causing HDF5 errors
        set<const Genome *> targets;
        targets.insert(targetGenome);
        
        // Create column iterator with duplications enabled (dupeMode = true)
        ColumnIteratorPtr colIt = refGenome->getColumnIterator(&targets, 0, absStart, absStart, true, false);
        if (colIt == NULL) {
            return false;
        }
        
        const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
        if (colMap == NULL) {
            return false;
        }
        
        bool foundData = false;
        for (ColumnIterator::ColumnMap::const_iterator cmIt = colMap->begin(); cmIt != colMap->end(); ++cmIt) {
            const Sequence *seq = cmIt->first;
            if (seq == NULL || seq->getGenome() != targetGenome) {
                continue;
            }
            
            ColumnIterator::DNASet *dnaSet = cmIt->second;
            if (dnaSet == NULL) {
                continue;
            }
            
            // Collect all bases from this sequence (including duplications)
            for (size_t j = 0; j < dnaSet->size(); ++j) {
                try {
                    DnaIteratorPtr dnaIt = dnaSet->at(j);
                    if (dnaIt != NULL) {
                        hal_index_t foundPos = dnaIt->getArrayIndex();
                        char base = toupper(dnaIt->getBase());
                        
                        // When searching within-species (refGenome == targetGenome),
                        // exclude the query position itself to avoid circular reasoning
                        if (refGenome == targetGenome && foundPos == absStart) {
                            continue;  // Skip the query position itself
                        }
                        
                        if (base != 'N' && base != '-' && base != 0) {
                            bases.push_back(base);
                            foundData = true;
                        }
                    }
                } catch (...) {
                    // Skip this DNA iterator if it fails
                    continue;
                }
            }
        }
        return foundData;
    } catch (...) {
        // Catch any exceptions and return false to fail gracefully
        return false;
    }
}

static bool findAncestralAllele(const Genome *refGenome, hal_index_t absStart,
                               const vector<set<const Genome *>>& targetSets,
                               const vector<const Genome*>& targetGenomes,
                               const vector<string>& ancestorNames,
                               string& ancAllele, string& evidence, string& usedAncestor) {
    // Try each ancestor in order
    for (size_t ancestorIdx = 0; ancestorIdx < targetSets.size(); ++ancestorIdx) {
        const set<const Genome *>& targets = targetSets[ancestorIdx];
        const Genome* tgtGenome = targetGenomes[ancestorIdx];
        
        // Step 1: Try direct orthologous lookup (without duplications) - main mapping
        ColumnIteratorPtr colIt = refGenome->getColumnIterator(&targets, 0, absStart, absStart, false, false);
        const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
        vector<char> bases;
        
        if (colIt != NULL && colMap != NULL) {
            for (ColumnIterator::ColumnMap::const_iterator cmIt = colMap->begin(); cmIt != colMap->end(); ++cmIt) {
                const Sequence *seq = cmIt->first;
                if (seq->getGenome() != tgtGenome) {
                    continue;
                }
                ColumnIterator::DNASet *dnaSet = cmIt->second;
                if (dnaSet != NULL) {
                    for (size_t j = 0; j < dnaSet->size(); ++j) {
                        char b = toupper(dnaSet->at(j)->getBase());
                        bases.push_back(b);
                    }
                }
            }
        }

        // Step 2: If no direct hit, try finding ancestral paralogs (with duplications enabled)
        bool foundParalogs = false;
        if (bases.empty()) {
            foundParalogs = findParalogsInGenome(refGenome, absStart, tgtGenome, bases);
        }

        if (!bases.empty()) {
            map<char, int> counts;
            int total = 0;
            for (size_t j = 0; j < bases.size(); ++j) {
                char b = bases[j];
                if (b == 'N' || b == '-') {
                    continue;
                }
                counts[b]++;
                total++;
            }
            
            if (total > 0) {
                usedAncestor = ancestorNames[ancestorIdx];
                
                // Add ancestor source info to evidence if using multiple ancestors
                string sourceInfo = "";
                if (ancestorNames.size() > 1) {
                    sourceInfo = string("@") + ancestorNames[ancestorIdx];
                    if (ancestorIdx > 0) {
                        sourceInfo += string("(fallback:") + to_string(ancestorIdx) + ")";
                    }
                }
                
                // Add paralog info if paralogs were used
                string methodInfo = foundParalogs ? "AncestralParalog" : "Direct";
                
                if (total == 1) {
                    ancAllele = string(1, counts.begin()->first);
                    evidence = methodInfo + sourceInfo;
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
                        evidence = string("AncestralParalogTie:") + info.str() + sourceInfo;
                    } else {
                        ancAllele = string(1, maj);
                        string voteType = foundParalogs ? "AncestralParalogVote:" : "MajorityVote:";
                        evidence = voteType + info.str() + sourceInfo;
                    }
                }
                return true; // Found an allele
            }
        }
    }
    
    // Step 3: Last resort - look for within-species paralogs in reference genome itself
    vector<char> refParalogs;
    bool foundRefParalogs = findParalogsInGenome(refGenome, absStart, refGenome, refParalogs);
    
    if (foundRefParalogs && !refParalogs.empty()) {
        map<char, int> counts;
        int total = 0;
        for (size_t j = 0; j < refParalogs.size(); ++j) {
            char b = refParalogs[j];
            if (b == 'N' || b == '-') {
                continue;
            }
            counts[b]++;
            total++;
        }
        
        if (total > 0) {
            usedAncestor = "WithinSpecies";  // Indicate within-species paralogs
            
            if (total == 1) {
                ancAllele = string(1, counts.begin()->first);
                evidence = "WithinSpeciesParalog";
            } else {
                // Find majority from reference paralogs
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
                    evidence = string("WithinSpeciesParalogTie:") + info.str();
                } else {
                    ancAllele = string(1, maj);
                    evidence = string("WithinSpeciesParalogVote:") + info.str();
                }
            }
            return true;
        }
    }
    
    // Truly no allele found anywhere
    ancAllele = "N";
    if (ancestorNames.size() > 1) {
        evidence = string("Missing(tried:") + to_string(ancestorNames.size()) + "+self)";
    } else {
        evidence = "Missing(+self)";
    }
    usedAncestor = ancestorNames.empty() ? "Unknown" : ancestorNames[0];
    return false;
}

int main(int argc, char **argv) {
    CLParser parser;
    initParser(parser);

    string halPath, refName, tgtName, posPath, outPath;
    bool noSort = false;
    size_t progressInterval = 0;
    try {
        parser.parseOptions(argc, argv);
        halPath = parser.getArgument<string>("halFile");
        refName = parser.getArgument<string>("refGenome");
        tgtName = parser.getArgument<string>("targetGenome");
        posPath = parser.getArgument<string>("positionsFile");
        outPath = parser.getArgument<string>("outputFile");
        noSort = parser.getFlag("noSort");
        progressInterval = parser.getOption<size_t>("progress");
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
        
        // Parse target genomes - auto-detect comma-separated lists
        vector<string> genomeNames;
        vector<const Genome*> targetGenomes;
        bool useMultipleGenomes = false;
        
        if (tgtName.find(',') != string::npos) {
            // Multiple genomes detected
            useMultipleGenomes = true;
            genomeNames = splitAncestors(tgtName);
            if (genomeNames.empty()) {
                throw hal_exception("No valid genome names provided in: " + tgtName);
            }
            cerr << "Using multiple genomes: ";
            for (size_t i = 0; i < genomeNames.size(); ++i) {
                if (i > 0) cerr << ", ";
                cerr << genomeNames[i];
                const Genome *genome = alignment->openGenome(genomeNames[i]);
                if (genome == NULL) {
                    throw hal_exception("Target genome " + genomeNames[i] + " not found");
                }
                targetGenomes.push_back(genome);
            }
            cerr << endl;
        } else {
            // Single genome
            genomeNames.push_back(tgtName);
            const Genome *tgtGenome = alignment->openGenome(tgtName);
            if (tgtGenome == NULL) {
                throw hal_exception("Target genome " + tgtName + " not found");
            }
            targetGenomes.push_back(tgtGenome);
        }

        ifstream posStream(posPath.c_str());
        if (!posStream) {
            throw hal_exception("Unable to open positions file: " + posPath);
        }
        ofstream outStream(outPath.c_str());
        if (!outStream) {
            throw hal_exception("Unable to open output file: " + outPath);
        }

        // Create target sets for each ancestor (for ColumnIterator)
        vector<set<const Genome *>> targetSets;
        for (size_t i = 0; i < targetGenomes.size(); ++i) {
            set<const Genome *> targets;
            targets.insert(targetGenomes[i]);
            targetSets.push_back(targets);
        }

        // First pass: read all positions
        vector<Position> positions;
        string line;
        size_t originalOrder = 0;
        while (getline(posStream, line)) {
            string chrom;
            hal_index_t start = 0, end = 0;
            if (!parseLine(line, chrom, start, end)) {
                continue;
            }
            Position pos;
            pos.chrom = chrom;
            pos.start = start;
            pos.end = end;
            pos.originalOrder = originalOrder++;
            positions.push_back(pos);
        }
        posStream.close();

        if (positions.empty()) {
            cerr << "No valid positions found in input file" << endl;
            return 1;
        }

        cerr << "Loaded " << positions.size() << " positions" << endl;

        // Optimization: sort by chromosome and position for better cache locality
        vector<size_t> processingOrder;
        if (noSort) {
            // Preserve original order
            for (size_t i = 0; i < positions.size(); ++i) {
                processingOrder.push_back(i);
            }
        } else {
            // Sort for cache efficiency
            for (size_t i = 0; i < positions.size(); ++i) {
                processingOrder.push_back(i);
            }
            sort(processingOrder.begin(), processingOrder.end(), 
                 [&positions](size_t a, size_t b) {
                     return positions[a] < positions[b];
                 });
            cerr << "Sorted positions for optimal processing" << endl;
        }

        // Store results to restore original order if needed
        vector<string> results(positions.size());

        // Process positions
        for (size_t i = 0; i < processingOrder.size(); ++i) {
            size_t idx = processingOrder[i];
            const Position& pos = positions[idx];

            if (progressInterval > 0 && (i + 1) % progressInterval == 0) {
                cerr << "Processed " << (i + 1) << "/" << positions.size() << " positions" << endl;
            }

            const Sequence *refSeq = refGenome->getSequence(pos.chrom);
            if (refSeq == NULL) {
                string usedGenome = useMultipleGenomes ? genomeNames[0] : tgtName;
                results[idx] = pos.chrom + "\t" + to_string(pos.start) + "\t" + to_string(pos.end) + 
                              "\tN\t" + usedGenome + "\tN\tMissing";
                continue;
            }
            
            hal_index_t absStart = refSeq->getStartPosition() + pos.start;
            DnaIteratorPtr refIt = refGenome->getDnaIterator(absStart);
            char refBase = toupper(refIt->getBase());

            string ancAllele = "N";
            string evidence = "Missing";
            string usedAncestor;
            
            findAncestralAllele(refGenome, absStart, targetSets, targetGenomes, genomeNames,
                               ancAllele, evidence, usedAncestor);

            results[idx] = pos.chrom + "\t" + to_string(pos.start) + "\t" + to_string(pos.end) + "\t" + 
                          refBase + "\t" + usedAncestor + "\t" + ancAllele + "\t" + evidence;
        }

        // Output results in original order
        for (size_t i = 0; i < results.size(); ++i) {
            outStream << results[i] << '\n';
        }
    } catch (hal_exception &e) {
        cerr << "hal exception: " << e.what() << endl;
        return 1;
    }
    return 0;
}

