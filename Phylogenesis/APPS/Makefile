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

LIBS = -lPhylogenesis -lAlign2  -ltools 

LIB_PATH = -L.

INC_PATH = -I. -I../../Phylogenesis/Sources -I../../tools/ -I../../Align2/Sources 

#
# Objects and headers 
#

SOURCES =  PhylogeneticTrees.cc Msa.cc

OBJECTS =  PhylogeneticTrees.o Msa.o

TARGETS =   PhylogeneticTrees Msa
 
EXECS =  PhylogeneticTrees Msa
 
LIBRARY = APPSlibPhylogenesis.a

#
# Install rule
#

compile: all

all: install
install: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(UPDIR)/lib/
	mv $(EXECS) $(BINPATH)/

#
# Call global Makefile to do the job
#

include ../../Makefile.global
