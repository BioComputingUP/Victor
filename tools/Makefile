#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = ..
# Path to subdirectories.
SUBDIR=
# Path to directory for binaries:
BINPATH = ../bin

#
# Libraries and paths (which are not defined globally).
#

LIBS =  -ltools

LIB_PATH = -L.

INC_PATH = -I.

#
# Objects and headersscale*
#

SOURCES = vector3.cc matrix3.cc vglStd.cc config.cc GetArg.cc \
 String2Number.cc timer.cc IoTools.cc StatTools.cc  Kabsch.cc
OBJECTS = vector3.o matrix3.o vglStd.o config.o GetArg.o \
 String2Number.o timer.o IoTools.o StatTools.o  Kabsch.o
TARGETS =  

LIBRARY = libtools.a

#
# Install rule
#

compile: all
all: install
install: $(LIBRARY) $(TARGET) 
	mkdir -p ../lib
	mv libtools.a ../lib/
	mkdir -p ../bin	
	mkdir -p ../Biopool/data	
	mkdir -p ../Lobo/data
	mkdir -p ../Energy/data
	export PATH=$(VICTOR_ROOT)bin/:$(PATH)
	ln -fs $(VICTOR_ROOT)data/AminoAcidHydrogenData.txt $(VICTOR_ROOT)Biopool/data/AminoAcidHydrogenData.txt
	ln -fs $(VICTOR_ROOT)data/allaminoacids_new.pdb $(VICTOR_ROOT)Energy/data/allaminoacids_new.pdb
	ln -fs $(VICTOR_ROOT)data/frst.model $(VICTOR_ROOT)Energy/data/frst.model
	ln -fs $(VICTOR_ROOT)data/frst2.range $(VICTOR_ROOT)Energy/data/frst2.range
	ln -fs $(VICTOR_ROOT)data/polar.par $(VICTOR_ROOT)Energy/data/polar.par
	ln -fs $(VICTOR_ROOT)data/ram.par $(VICTOR_ROOT)Energy/data/ram.par
	ln -fs $(VICTOR_ROOT)data/solv.par $(VICTOR_ROOT)Energy/data/solv.par
	ln -fs $(VICTOR_ROOT)data/tor.par $(VICTOR_ROOT)Energy/data/tor.par
	ln -fs $(VICTOR_ROOT)data/tor.par $(VICTOR_ROOT)Lobo/data/tor.par

#
# Call global Makefile to do the job.
# 
export VICTOR_ROOT=$(dir $(patsubst %/,%, $(shell pwd)))
include ../Makefile.global

