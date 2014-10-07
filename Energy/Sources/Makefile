#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = ../..
# Path to subdirectories.
SUBDIR = TorsionPotential
# Path to directory for binaries:
BINPATH = ../../bin


#
# Libraries and paths (which are not defined globally).
#

LIBS = -lEnergy -lTorsion -lBiopool -ltools

LIB_PATH = -L.

INC_PATH = -I. -I./TorsionPotential/ -I../../tools/

#
# Objects and headers
#
 

SOURCES =   PolarSolvationPotential.cc SolvationPotential.cc RapdfPotential.cc  \
            EnergyFeatures.cc \
          EffectiveSolvationPotential.cc   
           
           

OBJECTS = PolarSolvationPotential.o SolvationPotential.o RapdfPotential.o  \
            EnergyFeatures.o \
          EffectiveSolvationPotential.o  

 

TARGETS =   

EXECS =   

LIBRARY = libEnergy.a

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
export VICTOR_ROOT=$(dir $(patsubst %/,%, $(dir $(patsubst %/,%, $(shell pwd)))))
include ../../Makefile.global
