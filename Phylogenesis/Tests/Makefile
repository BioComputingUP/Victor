#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = ../..
# Path to subdirectories.
SUBDIR=
# Path to directory for binaries:
BINPATH = ../../bin


#
# Libraries and paths (which are not defined globally).
#

LIBS =  -lPhylogenesis -lAlign2 -lBiopool -ltools -L/usr/lib/ -lm -ldl -lcppunit

LIB_PATH = -L.

INC_PATH = -I. -I ../../Biopool/ -I../../tools/ -I../../Align2/Sources -I../../Phylogenesis/Sources 

#
# Objects and headers
#

SOURCES =  TestPhylogenesis.cc TestDistanceMatrix.h \
			TestClusterDistanceMatrix.h TestUPGMA.h TestNJ.h \
			TestNewickPhylogeneticTreeExporter.h TestClustalw.h

OBJECTS =  $(SOURCES:.cpp=.o)

TARGETS = TestPhylogenesis

EXECS = TestPhylogenesis

LIBRARY = TESTlibPhylogenesis.a

#
# Install rule
#

compile: all
	
all: install

install: $(LIBRARY) $(TARGETS)
	mv $(EXECS) $(UPDIR)/bin
	mv $(LIBRARY) $(UPDIR)/lib
	
#
# Call global Makefile to do the job.
#

include ../../Makefile.global
