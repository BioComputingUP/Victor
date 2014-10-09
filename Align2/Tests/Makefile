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

LIBS =  -lAlign2 -lBiopool -ltools -L/usr/lib/ -lm -ldl -lcppunit

LIB_PATH = -L.

INC_PATH = -I. -I ../../Biopool/ -I../../tools/ -I../../Align2/Sources

#
# Objects and headers
#

SOURCES =  TestAlign2.cc TestAlignmentData.h TestAlign.h TestAlignment.h

OBJECTS =  $(SOURCES:.cpp=.o)

TARGETS = TestAlign2

EXECS = TestAlign2

LIBRARY = TESTlibAlign2.a

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
