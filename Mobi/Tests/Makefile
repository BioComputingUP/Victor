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

LIBS = -lMobi -lBiopool -ltools -L/usr/lib/  -lm -ldl -lcppunit

LIB_PATH = -L.

INC_PATH = -I. -I../../tools/ -I../../Biopool/Sources -I../../Mobi/Sources

#
# Objects and headers
#

SOURCES =  TestMobi.cc TestVectorCollection.h TestMobiProtein.h TestTM.h TestMobiMethods.h

OBJECTS =  $(SOURCES:.cpp=.o)

TARGETS = TestMobi

EXECS = TestMobi

LIBRARY = TESTlibMobi.a



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
