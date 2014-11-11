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

LIBS =  -lEnergy -lTorsion -lBiopool  -ltools

LIB_PATH = -L.

INC_PATH = -I. -I../../tools/ -I../../Biopool/Sources -I../../Energy/Sources -I../../Energy/Sources/TorsionPotential

#
# Objects and headers
#

SOURCES =   PdbCorrector.cc PdbSecondary.cc PdbEditor.cc Pdb2Seq.cc pdb2secondary.cc pdbshifter.cc \
	pdbMover.cc

OBJECTS =   PdbCorrector.o PdbSecondary.o PdbEditor.o Pdb2Seq.o pdb2secondary.o pdbshifter.o \
	pdbMover.o

TARGETS = PdbCorrector PdbSecondary PdbEditor Pdb2Seq pdb2secondary pdbshifter \
	pdbMover

EXECS = PdbCorrector PdbSecondary PdbEditor Pdb2Seq pdb2secondary pdbshifter \
	pdbMover

LIBRARY = APPSlibBiopool.a

#
# Install rule
#

compile: all

all: install
install: $(LIBRARY) $(TARGETS)
	mv $(EXECS) $(UPDIR)/bin/
	mv $(LIBRARY) $(UPDIR)/lib/
	
#
# Call global Makefile to do the job.
#

include ../../Makefile.global
