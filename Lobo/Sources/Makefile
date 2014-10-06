#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = ../..
# Path to subdirectories.
SUBDIR=../../tools
# Path to directory for binaries:
BINPATH =  ../../bin


#
# Libraries and paths (which are not defined globally).
#

LIBS = -lLobo -lBiopool  -ltools -lEnergy -lTorsion

LIB_PATH = -L.

INC_PATH = -I.   -I$(UPDIR)/Energy/Sources -I$(UPDIR)/Energy/Sources/TorsionPotential -I$(UPDIR)/Biopool  -I$(UPDIR)/tools -I$(UPDIR)/Lobo
#
# Objects and headers
#

SOURCES = RamachandranData.cc VectorTransformation.cc LoopTableEntry.cc  \
   LoopTable.cc LoopModel.cc LoopExtractor.cc ranking_helper.cc \
   ranking_helper2.cc  globalStatistic.cc RankAnalyzer.cc  

OBJECTS = RamachandranData.o VectorTransformation.o   LoopTableEntry.o \
   LoopTable.o LoopModel.o   LoopExtractor.o ranking_helper.o ranking_helper2.o \
    globalStatistic.o RankAnalyzer.o  

TARGETS =  
    

EXECS =   

LIBRARY = libLobo.a

#
# Install rule
#

compile: all
	
allinstall: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(UPDIR)/lib
	 

install: subinstall allinstall
	
all: install	
#
# Call global Makefile to do the job.
#
export VICTOR_ROOT=$(dir $(patsubst %/,%, $(dir $(patsubst %/,%, $(shell pwd)))))
include $(UPDIR)/Makefile.global