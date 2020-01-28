###################################################
#
# Makefile for genomegaMap
#
###################################################

#
# Macros
#

# Header file locations for gcat-project
INCLUDE = -Isrc -I/usr/include/gcat -I/usr/include/gcat/myutils
# C++ compiler
CC = g++
# C++ compiler for MPICH2
MPICC = mpic++
# C++ linker
LD = g++
# C++ standard compiler options
CXXFLAGS = -Wall -w -O3 -g -D __NOEXTERN_FOR_CINCLUDE -D _NOT_MYUTILS_DEBUG
# C++ compiler options for gcat library code
CC_OPTIONS = $(CXXFLAGS) -fPIC
# C++ linker options
#LNK_OPTIONS = -L/mnt/lustre/home/djw/xerces-c-3.0.1-x86_64-linux-gcc-3.4/lib\
#		-L/usr/lib64/mpich2\
#		-lxerces-c\
#		-lgsl\
#		-lgslcblas
LNK_OPTIONS = -Wl,--no-as-needed -lxerces-c -L./

#
# Build genomegaMap
#

LIBGENOMEGAMAP_OBJECTS = \
		./Distribution_genomegaMap.o\
		./MCMC_genomegaMapMoves.o\
		./RandomVariable_Codon61Count.o\
		./RandomVariable_Codon61Sequence.o\
		./Transformation_NY98_PDRM.o\
		./Transformation_NY98_PDRM_Single.o\
		./Transformation_Omega2GammaVector.o\
		./Transformation_ParentDependentRateMatrix.o\
		./Utility_mutation.o\
		./Utility_genomegaMapUtils.o\
		./Utility_paml.o\
		./genomegaMapLibrary.o\
		./genomegaMapXML.o

all : libgenomegaMap.so

clean : cleanobj
	rm libgenomegaMap.so

cleanobj :
	rm $(LIBGENOMEGAMAP_OBJECTS)

#install : gcat
#		cp gcat gcat

libgenomegaMap.so : $(LIBGENOMEGAMAP_OBJECTS)
	$(LD) $(LNK_OPTIONS) -lgsl -lgcat-core -shared -o libgenomegaMap.so $(LIBGENOMEGAMAP_OBJECTS)

#
# Build the parts of genomegaMap
#

./Distribution_genomegaMap.o : src/genomegaMap/Distributions/genomegaMap.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/Distributions/genomegaMap.cpp -c $(INCLUDE) -o ./Distribution_genomegaMap.o

./MCMC_genomegaMapMoves.o : src/genomegaMap/Inference/MCMC/genomegaMapMoves.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/Inference/MCMC/genomegaMapMoves.cpp -c $(INCLUDE) -o ./MCMC_genomegaMapMoves.o

./RandomVariable_Codon61Count.o : src/genomegaMap/RandomVariables/Codon61Count.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/RandomVariables/Codon61Count.cpp -c $(INCLUDE) -o ./RandomVariable_Codon61Count.o

./RandomVariable_Codon61Sequence.o : src/genomegaMap/RandomVariables/Codon61Sequence.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/RandomVariables/Codon61Sequence.cpp -c $(INCLUDE) -o ./RandomVariable_Codon61Sequence.o

./Transformation_NY98_PDRM.o : src/genomegaMap/Transformations/NY98_PDRM.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/Transformations/NY98_PDRM.cpp -c $(INCLUDE) -o ./Transformation_NY98_PDRM.o

./Transformation_NY98_PDRM_Single.o : src/genomegaMap/Transformations/NY98_PDRM_Single.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/Transformations/NY98_PDRM_Single.cpp -c $(INCLUDE) -o ./Transformation_NY98_PDRM_Single.o

./Transformation_Omega2GammaVector.o : src/genomegaMap/Transformations/Omega2GammaVector.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/Transformations/Omega2GammaVector.cpp -c $(INCLUDE) -o ./Transformation_Omega2GammaVector.o

./Transformation_ParentDependentRateMatrix.o : src/genomegaMap/Transformations/ParentDependentRateMatrix.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/Transformations/ParentDependentRateMatrix.cpp -c $(INCLUDE) -o ./Transformation_ParentDependentRateMatrix.o

./Utility_mutation.o : src/genomegaMap/Utilities/mutation.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/Utilities/mutation.cpp -c $(INCLUDE) -o ./Utility_mutation.o

./Utility_genomegaMapUtils.o : src/genomegaMap/Utilities/genomegaMapUtils.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/Utilities/genomegaMapUtils.cpp -c $(INCLUDE) -o ./Utility_genomegaMapUtils.o

./Utility_paml.o : src/genomegaMap/Utilities/paml.c
	$(CC) $(CC_OPTIONS) src/genomegaMap/Utilities/paml.c -c $(INCLUDE) -o ./Utility_paml.o

./genomegaMapLibrary.o : src/genomegaMap/genomegaMapLibrary.cpp
	$(CC) $(CC_OPTIONS) src/genomegaMap/genomegaMapLibrary.cpp -c $(INCLUDE) -o ./genomegaMapLibrary.o

./genomegaMapXML.o : src/genomegaMap/genomegaMapXML.cpp src/genomegaMap/genomegaMap1.0.xsd.h
	$(CC) $(CC_OPTIONS) src/genomegaMap/genomegaMapXML.cpp -c $(INCLUDE) -o ./genomegaMapXML.o

src/genomegaMap/genomegaMap1.0.xsd.h : src/genomegaMap/genomegaMap1.0.xsd
	(cd src/genomegaMap && xxd -i genomegaMap1.0.xsd > genomegaMap1.0.xsd.h)

##### END RUN ####
