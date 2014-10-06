#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = ../../..
# Path to subdirectories.
SUBDIR=../../tools 
# Path to directory for binaries:
BINPATH = ../../../bin


#
# Libraries and paths (which are not defined globally).
#

LIBS = -lEnergy -lTorsion -lBiopool  -ltools

LIB_PATH = -L. 

INC_PATH = -I.  -I../../tools  

#
# Objects and headers
#

SOURCES =  PhiPsi.cc  PhiPsiOmegaChi1.cc PhiPsiOmegaChi1Chi2PreAngle.cc PhiPsiPreAngle.cc PhiPsiOmegaPreAngle.cc\
           PhiPsiOmegaChi1Chi2.cc PhiPsiOmega.cc Chi1Chi2.cc PhiPsiOmegaPreAngle.cc PhiPsiOmegaChi1PreAngle.cc Omega.cc

OBJECTS =  PhiPsi.o    PhiPsiOmega.o PhiPsiOmegaChi1.o  PhiPsiOmegaPreAngle.o PhiPsiOmegaChi1PreAngle.o\
           PhiPsiOmegaChi1Chi2.o  PhiPsiOmegaChi1Chi2PreAngle.o  PhiPsiPreAngle.o Chi1Chi2.o PhiPsiOmegaPreAngle.o Omega.o

TARGETS =

EXECS =

LIBRARY = libTorsion.a

#
# Install rule
#

compile: all
allinstall: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) ../../../lib
	 

install: subinstall allinstall

all: install	
	
#
# Call global Makefile to do the job.
#
export VICTOR_ROOT=$(dir $(patsubst %/,%, $(dir $(patsubst %/,%, $(shell pwd)))))
include ../../../Makefile.global

