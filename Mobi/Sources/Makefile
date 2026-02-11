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

LIBS = -lBiopool -ltools

LIB_PATH = -L.

INC_PATH = -I. -I../../tools -I../../Biopool/Sources

#
# Objects and headers
#
 

SOURCES = 	 ProteinModels.cc	TmScore.cc	StandardDeviation.cc	SecondaryStructure.cc	MobiSaver.cc


OBJECTS = 	 ProteinModels.o	TmScore.o	StandardDeviation.o	SecondaryStructure.o	MobiSaver.o


TARGETS =   

EXECS =   

LIBRARY = libMobi.a

#
# Install rule
#

#start: subdirs

compile: all

allinstall: $(LIBRARY) $(TARGETS) 
	    mv $(LIBRARY) $(UPDIR)/lib/
	 

install: subinstall allinstall

all: install
#
# Call global Makefile to do the job.
#

include ../../Makefile.global

