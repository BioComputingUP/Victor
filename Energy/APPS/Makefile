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
 

SOURCES =   
           
OBJECTS =    
 

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

allinstall: $(LIBRARY) $(TARGETS)
	mv $(EXECS) $(UPDIR)/bin 
	mv $(LIBRARY) $(UPDIR)/lib 
	echo $(VICTOR_ROOT)
	export PATH=$(VICTOR_ROOT)bin/:$(PATH)

install: subinstall allinstall

#
# Call global Makefile to do the job.
#

export VICTOR_ROOT=$(dir $(patsubst %/,%, $(dir $(patsubst %/,%, $(shell pwd)))))
include ../../Makefile.global
