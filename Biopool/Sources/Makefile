#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = ../..
# Path to subdirectories.
SUBDIR = ../../tools
# Path to directory for binaries:
BINPATH = ../../bin


#
# Libraries and paths (which are not defined globally).
#

LIBS = -lEnergy -lTorsion -ltools

LIB_PATH = -L.

INC_PATH = -I.   -I../../tools 

#
# Objects and headers
#
 

SOURCES = Identity.cc SimpleBond.cc Bond.cc \
 Component.cc Monomer.cc Polymer.cc Atom.cc Group.cc SideChain.cc \
 AminoAcid.cc Spacer.cc IntSaver.cc IntLoader.cc XyzLoader.cc \
 XyzSaver.cc SeqSaver.cc PdbLoader.cc PdbSaver.cc SeqLoader.cc \
 IntCoordConverter.cc SideChainConstructor.cc RelLoader.cc RelSaver.cc \
 SeqConstructor.cc Ligand.cc LigandSet.cc SolvExpos.cc AminoAcidHydrogen.cc Nucleotide.cc

OBJECTS = Identity.o SimpleBond.o Bond.o \
 Component.o Monomer.o Polymer.o Atom.o Group.o SideChain.o \
 AminoAcid.o Spacer.o IntSaver.o IntLoader.o XyzLoader.o \
 XyzSaver.o SeqSaver.o PdbLoader.o PdbSaver.o SeqLoader.o \
 IntCoordConverter.o SideChainConstructor.o RelLoader.o RelSaver.o \
 SeqConstructor.o Ligand.o LigandSet.o SolvExpos.o Protein.o   AminoAcidHydrogen.o Nucleotide.o
 

TARGETS =   

EXECS =   

LIBRARY = libBiopool.a

#
# Install rule
#

#start: subdirs

compile: all

allinstall: $(LIBRARY) $(TARGETS) 
	    mv $(LIBRARY) $(UPDIR)/lib
	 

install: subinstall allinstall

all: install
#
# Call global Makefile to do the job.
#

include ../../Makefile.global

