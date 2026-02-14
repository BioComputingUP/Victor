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

LIBS = -lAlign2 -ltools 

LIB_PATH = -L.

INC_PATH = -I. -I../../tools/ -I../../Align2/Sources

#
# Objects and headers
#

SOURCES = PhyloSupport.cc NewickTree.cc ClustalW.cc SeqNodeGraph.cc

OBJECTS = PhyloSupport.o NewickTree.o ClustalW.o SeqNodeGraph.o

TARGETS =  

EXECS =  

LIBRARY = libPhylo.a

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
