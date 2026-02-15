 
#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory
UPDIR = ../..
# Path to subdirectories
SUBDIR = 
# Path to directory for binaries
BINPATH = ../../bin

#
# Libraries and paths (which are not defined globally)
#
		
LIBS = -lAlign2 -lPhylogenesis -lBiopool -ltools -lPhylogenesis

LIB_PATH = -L.

INC_PATH = -I. -I../../Align2/Sources -I ../../Biopool/Sources -I../../tools/ 

#
# Objects and headers
#

SOURCES = DistanceMatrix.cc ClusterDistanceMatrix.cc UPGMA.cc NewickPhylogeneticTreeExporter.cc NeighborJoining.cc \
	ClustalW.cc SubMatrixFactory.cc \
	MsaNWG.cc CWScoringScheme.cc CWGapFunction.cc Sequence.cc \
	MultipleAlignment.cc Cluster.cc

OBJECTS = DistanceMatrix.o ClusterDistanceMatrix.o UPGMA.o NewickPhylogeneticTreeExporter.o \
	NeighborJoining.o ClustalW.o \
	MsaNWG.o CWScoringScheme.o CWGapFunction.o Sequence.o SubMatrixFactory.o \
	MultipleAlignment.o Cluster.o Sequence.o

TARGETS =  

EXECS =  

LIBRARY = libPhylogenesis.a

#
# Install rule
#

compile: all

install: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(UPDIR)/lib

#
# Call global Makefile to do the job
#
all: install
include ../../Makefile.global
