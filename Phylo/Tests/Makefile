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

LIBS =  -lPhylo -lAlign2 -ltools -L/usr/lib/ -lm -ldl -lcppunit

LIB_PATH = -L.

INC_PATH = -I. -I../../tools/ -I../../Align2/Sources -I../../Phylo/Sources

#
# Objects and headers
#

SOURCES =  TestPhylogen.cc TestPhylo.h

OBJECTS =  $(SOURCES:.cpp=.o)

TARGETS = TestPhylogen

EXECS = TestPhylogen

LIBRARY = TESTlibPhylogen.a

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
