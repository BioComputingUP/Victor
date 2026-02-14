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
	pdbMover.cc CifEditor.cc CifSecondary.cc Cif2Secondary.cc CifMover.cc \
	CifShifter.cc CifCorrector.cc Cif2Seq.cc CifReaderWriter.cc

OBJECTS =   PdbCorrector.o PdbSecondary.o PdbEditor.o Pdb2Seq.o pdb2secondary.o pdbshifter.o \
	pdbMover.o CifEditor.o CifSecondary.o Cif2Secondary.o CifMover.o \
	CifShifter.o CifCorrector.o Cif2Seq.o CifReaderWriter.o

TARGETS = PdbCorrector PdbSecondary PdbEditor Pdb2Seq pdb2secondary pdbshifter \
	pdbMover CifEditor CifSecondary Cif2Secondary CifMover CifShifter \
	CifCorrector Cif2Seq CifReaderWriter

EXECS = PdbCorrector PdbSecondary PdbEditor Pdb2Seq pdb2secondary pdbshifter \
	pdbMover CifEditor CifSecondary Cif2Secondary CifMover CifShifter \
	CifCorrector Cif2Seq CifReaderWriter

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
