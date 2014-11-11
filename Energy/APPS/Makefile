#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = ../..
# Path to subdirectories.
SUBDIR = ../Sources/TorsionPotential
# Path to directory for binaries:
BINPATH = ../../bin


#
# Libraries and paths (which are not defined globally).
#

LIBS = -lEnergy -lTorsion -lBiopool  -ltools

LIB_PATH = -L.

INC_PATH = -I. -I../Sources/TorsionPotential/ -I../../tools/ -I../../Biopool/Sources -I../../Energy/Sources

#
# Objects and headers
#
 

SOURCES =  frst.cc  correlation.cc  energy2zscore.cc   frstZscore.cc  mutationGenerator.cc   \
              pdb2tor.cc    tap2plot.cc  pdb2solv.cc  \
          solv2energy.cc  pdb2contact.cc  tapRef.cc  pdb2energy.cc taptable.cc pdb2tap.cc pdb2torenergy.cc 
           
OBJECTS =  frst.o  correlation.o  energy2zscore.o   frstZscore.o  mutationGenerator.o   \
              pdb2tor.o    tap2plot.o  pdb2solv.o  \
          solv2energy.o  pdb2contact.o  tapRef.o  pdb2energy.o taptable.o pdb2tap.o pdb2torenergy.o 
 

TARGETS =  frst  correlation  energy2zscore   frstZscore  mutationGenerator   \
              pdb2tor    tap2plot  pdb2solv  \
          solv2energy  pdb2contact  tapRef  pdb2energy taptable pdb2tap pdb2torenergy 
 

EXECS =  frst  correlation  energy2zscore   frstZscore  mutationGenerator   \
              pdb2tor    tap2plot    pdb2solv \
          solv2energy  pdb2contact  tapRef  pdb2energy taptable pdb2tap pdb2torenergy 
           
LIBRARY = APPSlibEnergy.a



#
# Install rule
#

#start: subdirs

compile: all

all: install
allinstall: $(LIBRARY) $(TARGETS)
	mv $(EXECS) $(UPDIR)/bin/
	mv $(LIBRARY) $(UPDIR)/lib/
	echo $(VICTOR_ROOT)
	export PATH=$(VICTOR_ROOT)bin/:$(PATH)

install: subinstall allinstall

#
# Call global Makefile to do the job.
#

export VICTOR_ROOT=$(dir $(patsubst %/,%, $(dir $(patsubst %/,%, $(shell pwd)))))
include ../../Makefile.global
