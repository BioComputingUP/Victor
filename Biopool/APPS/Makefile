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

SOURCES =  

OBJECTS =  

TARGETS = AtomTest GroupTest AminoAcidTest SpacerTest LigandSetTest \
	ReferenceTest SideChainTest BuildTest IntSpacerTest \
	PdbCorrector PdbSecondary PdbEditor Pdb2Seq pdb2secondary pdbshifter \
	pdbMover  addCbTest ProteinTest ProteinTestLoadAllchains CreateReference InsertTest

EXECS = AtomTest GroupTest AminoAcidTest SpacerTest LigandSetTest \
	ReferenceTest SideChainTest BuildTest IntSpacerTest \
	PdbCorrector PdbSecondary PdbEditor Pdb2Seq pdb2secondary pdbshifter \
	pdbMover  addCbTest ProteinTest ProteinTestLoadAllchains CreateReference InsertTest

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
