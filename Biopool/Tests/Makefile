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

LIBS =  -lEnergy -lTorsion -lBiopool  -ltools -L/usr/lib/ -lm -ldl -lcppunit

LIB_PATH = -L.

INC_PATH = -I. -I../../tools/ -I../../Biopool/Sources -I../../Energy/Sources -I../../Energy/Sources/TorsionPotential

#
# Objects and headers
#

SOURCES =  TestBiopool.cc Lin2dSolver.h TestLin2dSolver.h TestAtom.h

OBJECTS =  $(SOURCES:.cpp=.o)

TARGETS = TestBiopool

EXECS = TestBiopool

LIBRARY = TESTlibBiopool.a

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
