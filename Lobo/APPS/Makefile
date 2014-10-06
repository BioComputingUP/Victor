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
LIB = ../../lib

#
# Libraries and paths (which are not defined globally).
#

LIBS = -L$(LIB) -lLobo -lBiopool  -ltools -lEnergy -lTorsion


INC_PATH = -I.  -I../../Energy/Sources -I../../Energy/TorsionPotential/Sources -I../../Biopool/Sources  -I../../tools -I../../Lobo/Sources
#
# Objects and headers
#

SOURCES = 
    
    

OBJECTS =   

TARGETS = loboLUT LoopTablePlot  ClusterLoopTable  \
     ClusterRama     LoopTableTest LoopModelTest FullModelTest lobo  \
   loop2torsion scatEdit   createLoopTestset backboneAnalyzer  loboFull
# 
EXECS = loboLUT LoopTablePlot   ClusterLoopTable  \
     ClusterRama       LoopTableTest FullModelTest lobo   \
   loop2torsion scatEdit  createLoopTestset backboneAnalyzer  LoopModelTest  \
        loboLUT_all  loboFull

LIBRARY = APPSlibLobo.a

#
# Install rule
#

compile: all

install: $(LIBRARY) $(TARGETS)
	mv $(EXECS) $(UPDIR)/bin
	mv $(LIBRARY) $(UPDIR)/lib 
	cp $(UPDIR)/bin/loboLUT_all ./loboLUT_all
	

#
# Call global Makefile to do the job.
#
export VICTOR_ROOT=$(dir $(patsubst %/,%, $(dir $(patsubst %/,%, $(shell pwd)))))
include $(UPDIR)/Makefile.global
