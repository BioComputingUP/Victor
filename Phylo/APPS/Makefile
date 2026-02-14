#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory
UPDIR = ../..
# Path to subdirectories
SUBDIR =
# Path to directory for binaries
BINPATH = ../../bin

#
# Libraries and paths (which are not defined globally)
#

LIBS = -lPhylo -lAlign2 -ltools 

LIB_PATH = -L.

INC_PATH = -I. -I../../tools/ -I../../Align2/Sources -I../../Phylo/Sources

#
# Objects and headers
#

SOURCES =  phylogen.cc 

OBJECTS =  phylogen.o 
 
TARGETS =  phylogen
 

EXECS =  phylogen
 

LIBRARY = APPSlibPhylo.a

#
# Install rule
#

compile: all

all: install
install: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(UPDIR)/lib/
	mv $(EXECS) $(BINPATH)/

#
# Call global Makefile to do the job
#

include ../../Makefile.global
