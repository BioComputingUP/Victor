#--*-makefile-*-----------------------------------------------------------------
#
#   Standard makefile - Phylo
#
#-------------------------------------------------------------------------------

# Path to project directory
PROJECT_ROOT = ../..
BINPATH      = $(PROJECT_ROOT)/bin/
LIBPATH      = $(PROJECT_ROOT)/lib/
UPDIR        = $(PROJECT_ROOT)
SUBDIR       = 



#
# Libraries and paths (which are not defined globally)
#

LIBS     =
INC_PATH = -I.
LIB_PATH = -L.



#
# Objects and headers
#

SOURCES = PhylogeneticTree.cc RootedTree.cc UnrootedTree.cc NewickParser.cc \
          Sequence.cc AminoAcid.cc SubstitutionMatrix.cc \
          DistanceMatrix.cc DistanceMatrixBuilder.cc \
          IdentityPercentage.cc LevenshteinDistance.cc FengDoolittleDistance.cc \
          FitchMargoliash.cc UPGMA.cc NJ.cc \
          MultipleAlignment.cc \
          FengDoolittle.cc ClustalW.cc
OBJECTS = PhylogeneticTree.o RootedTree.o UnrootedTree.o NewickParser.o \
          Sequence.o AminoAcid.o SubstitutionMatrix.o \
          DistanceMatrix.o DistanceMatrixBuilder.o \
          IdentityPercentage.o LevenshteinDistance.o FengDoolittleDistance.o \
          FitchMargoliash.o UPGMA.o NJ.o \
          MultipleAlignment.o \
          FengDoolittle.o ClustalW.o
TARGETS = 
EXECS   = 
LIBRARY = libPhylo.a



#
# Install rule
#

compile: all

all: install
install: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(LIBPATH)



#
# Call global Makefile to do the job
#
include ../../Makefile.global

# Overriding to use C++11
STANDARDFLAGS = -Wall -Wextra -pedantic -std=c++0x -DNEXCEPTIONS -DLINUX -c
