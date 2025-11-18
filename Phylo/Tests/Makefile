#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory
PROJECT_ROOT = ../..
BINPATH      = $(PROJECT_ROOT)/bin/
LIBPATH      = $(PROJECT_ROOT)/lib/
UPDIR        = $(PROJECT_ROOT)
SUBDIR       = 



#
# Libraries and paths (which are not defined globally).
#

LIBS     = -lPhylo -lAlign2 -ltools -L/usr/lib/ -lm -ldl -lcppunit
INC_PATH = -I. -I ../../Phylo/Sources -I../../tools/ -I../../Align2/Sources
LIB_PATH = -L.



#
# Objects and headers
#

SOURCES =  TestPhylo.cc TestTree.h TestDistanceMatrix.h TestClustering.h TestMSA.h
OBJECTS =  $(SOURCES:.cpp=.o)
TARGETS = TestPhylo
EXECS   = TestPhylo
LIBRARY = TESTlibPhylo.a



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

# Overriding to use C++11
STANDARDFLAGS = -Wall -Wextra -pedantic -std=c++0x -DNEXCEPTIONS -DLINUX -c