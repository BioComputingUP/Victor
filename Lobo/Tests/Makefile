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

LIBS =  -L$(LIB) -lLobo -lBiopool  -ltools -lEnergy -lTorsion -L/usr/lib/ -lm -ldl -lcppunit

LIB_PATH = -L.

INC_PATH = -I.  -I../../Energy/Sources -I../../Energy/TorsionPotential/Sources -I../../Biopool/Sources  -I../../tools -I../../Lobo/Sources

#
# Objects and headers
#

SOURCES =  TestLobo.cc TestVectorTransformation.h

OBJECTS =  $(SOURCES:.cpp=.o)

TARGETS = TestLobo

EXECS = TestLobo

LIBRARY = TESTlibLobo.a

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
