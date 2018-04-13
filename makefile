###################################################
#
# Makefile for gcat-omegaMap
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
CXXFLAGS = -Wall -w -O3 -g -D __NOEXTERN_FOR_CINCLUDE -D _NOT_MYUTILS_DEBUG -arch $(ARCH)
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
# Build gcat-omegaMap
#

LIB_OMEGAMAP_OBJECTS = \
		./Distribution_omegaMapNCD.o\
		./Distribution_omegaMapUnlinked.o\
		./MCMC_omegaMapMoves.o\
		./RandomVariable_Codon61Count.o\
		./RandomVariable_Codon61Sequence.o\
		./Transformation_NY98_PDRM.o\
		./Transformation_NY98_PDRM_Single.o\
		./Transformation_Omega2GammaVector.o\
		./Transformation_ParentDependentRateMatrix.o\
		./Utility_mutation.o\
		./Utility_omegaMapUtils.o\
		./Utility_paml.o\
		./omegaMapLibrary.o\
		./omegaMapXML.o

all : libgcat_omegaMap.so

clean : cleanobj
	rm libgcat_omegaMap.so

cleanobj :
	rm $(LIB_OMEGAMAP_OBJECTS)

#install : gcat
#		cp gcat gcat

libgcat_omegaMap.so : $(LIB_OMEGAMAP_OBJECTS)
	$(LD) $(LNK_OPTIONS) -lgsl -lgcat-core -shared -o libgcat_omegaMap.so $(LIB_OMEGAMAP_OBJECTS)

#
# Build the parts of gcat_omegaMap
#

./Distribution_omegaMapNCD.o : src/omegaMap/Distributions/omegaMapNCD.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Distributions/omegaMapNCD.cpp -c $(INCLUDE) -o ./Distribution_omegaMapNCD.o

./Distribution_omegaMapUnlinked.o : src/omegaMap/Distributions/omegaMapUnlinked.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Distributions/omegaMapUnlinked.cpp -c $(INCLUDE) -o ./Distribution_omegaMapUnlinked.o

./MCMC_omegaMapMoves.o : src/omegaMap/Inference/MCMC/omegaMapMoves.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Inference/MCMC/omegaMapMoves.cpp -c $(INCLUDE) -o ./MCMC_omegaMapMoves.o

./RandomVariable_Codon61Count.o : src/omegaMap/RandomVariables/Codon61Count.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/RandomVariables/Codon61Count.cpp -c $(INCLUDE) -o ./RandomVariable_Codon61Count.o

./RandomVariable_Codon61Sequence.o : src/omegaMap/RandomVariables/Codon61Sequence.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/RandomVariables/Codon61Sequence.cpp -c $(INCLUDE) -o ./RandomVariable_Codon61Sequence.o

./Transformation_NY98_PDRM.o : src/omegaMap/Transformations/NY98_PDRM.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Transformations/NY98_PDRM.cpp -c $(INCLUDE) -o ./Transformation_NY98_PDRM.o

./Transformation_NY98_PDRM_Single.o : src/omegaMap/Transformations/NY98_PDRM_Single.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Transformations/NY98_PDRM_Single.cpp -c $(INCLUDE) -o ./Transformation_NY98_PDRM_Single.o

./Transformation_Omega2GammaVector.o : src/omegaMap/Transformations/Omega2GammaVector.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Transformations/Omega2GammaVector.cpp -c $(INCLUDE) -o ./Transformation_Omega2GammaVector.o

./Transformation_ParentDependentRateMatrix.o : src/omegaMap/Transformations/ParentDependentRateMatrix.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Transformations/ParentDependentRateMatrix.cpp -c $(INCLUDE) -o ./Transformation_ParentDependentRateMatrix.o

./Utility_mutation.o : src/omegaMap/Utilities/mutation.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Utilities/mutation.cpp -c $(INCLUDE) -o ./Utility_mutation.o

./Utility_omegaMapUtils.o : src/omegaMap/Utilities/omegaMapUtils.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/Utilities/omegaMapUtils.cpp -c $(INCLUDE) -o ./Utility_omegaMapUtils.o

./Utility_paml.o : src/omegaMap/Utilities/paml.c
	$(CC) $(CC_OPTIONS) src/omegaMap/Utilities/paml.c -c $(INCLUDE) -o ./Utility_paml.o

./omegaMapLibrary.o : src/omegaMap/omegaMapLibrary.cpp
	$(CC) $(CC_OPTIONS) src/omegaMap/omegaMapLibrary.cpp -c $(INCLUDE) -o ./omegaMapLibrary.o

./omegaMapXML.o : src/omegaMap/omegaMapXML.cpp src/omegaMap/gcat_omegaMap1.0.xsd.h
	$(CC) $(CC_OPTIONS) src/omegaMap/omegaMapXML.cpp -c $(INCLUDE) -o ./omegaMapXML.o

src/omegaMap/gcat_omegaMap1.0.xsd.h : src/omegaMap/omegaMap1.0.xsd
	(cd src/omegaMap && xxd -i omegaMap1.0.xsd > gcat_omegaMap1.0.xsd.h)

##### END RUN ####
