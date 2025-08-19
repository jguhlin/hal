rootDir = ..
include ${rootDir}/include.mk
modObjDir = ${objDir}/mutations

libHalMutations_srcs = impl/halBranchMutations.cpp impl/halMutationsStats.cpp impl/halSummarizeMutations.cpp
libHalMutations_objs = ${libHalMutations_srcs:%.cpp=${modObjDir}/%.o}
halIndels_srcs = impl/halIndels.cpp
halIndels_objs = ${halIndels_srcs:%.cpp=${modObjDir}/%.o}
halBranchMutations_srcs = impl/halBranchMutationsMain.cpp
halBranchMutations_objs = ${halBranchMutations_srcs:%.cpp=${modObjDir}/%.o}
halSnps_srcs = impl/halSnps.cpp
halSnps_objs = ${halSnps_srcs:%.cpp=${modObjDir}/%.o}
halAncestralAllele_srcs = impl/halAncestralAllele.cpp
halAncestralAllele_objs = ${halAncestralAllele_srcs:%.cpp=${modObjDir}/%.o}
halSummarizeMutations_srcs = impl/halSummarizeMutationsMain.cpp
halSummarizeMutations_objs = ${halSummarizeMutations_srcs:%.cpp=${modObjDir}/%.o}
srcs = ${libHalMutations_srcs} ${halIndels_srcs} ${halBranchMutations_srcs} ${halSnps_srcs} ${halAncestralAllele_srcs} ${halSummarizeMutations_srcs}
objs = ${srcs:%.cpp=${modObjDir}/%.o}
depends = ${srcs:%.cpp=%.depend}
progs = ${binDir}/halIndels ${binDir}/halBranchMutations ${binDir}/halSnps ${binDir}/halAncestralAllele ${binDir}/halSummarizeMutations
otherLibs = ${libHalMutations}

all : libs progs
libs: ${libHalMutations}
progs: ${progs}

clean :
	rm -f ${libHalMutations} ${objs} ${progs} ${depends}

test: ${binDir}/halAncestralAllele ../bin/halRandGen ../bin/halLiftover ../bin/hal2fasta
	@mkdir -p tests/output
	../bin/halRandGen --preset small --seed 0 --testRand --format hdf5 tests/output/test.hdf5.hal
	echo 'Genome_1_seq 0 1' > tests/output/positions.bed
	echo 'Genome_1_seq 1 2' >> tests/output/positions.bed
	${binDir}/halAncestralAllele tests/output/test.hdf5.hal Genome_1 Genome_0 tests/output/positions.bed tests/output/out.tsv
	test -s tests/output/out.tsv
include ${rootDir}/rules.mk

# don't fail on missing dependencies, they are first time the .o is generates
-include ${depends}


# Local Variables:
# mode: makefile-gmake
# End:

