#--*-makefile-*-----------------------------------------------------------------
#
#   Standard makefile - Phylo/APPS
#
#-------------------------------------------------------------------------------

# Path to project directory
PROJECT_ROOT = ../..
BINPATH      = $(PROJECT_ROOT)/bin/
LIBPATH      = $(PROJECT_ROOT)/lib/
UPDIR        = $(PROJECT_ROOT)
SUBDIR       = 



#
# Libraries and paths (which are not defined globally)
#

LIBS     = -lPhylo -lAlign2 -ltools
INC_PATH = -I. -I $(PROJECT_ROOT)/Phylo/Sources
LIB_PATH = -L. 



#
# Objects and headers
#

SOURCES = phyltree.cc clustalw.cc fengdoolittle.cc
OBJECTS = phyltree.o clustalw.o fengdoolittle.o
TARGETS = phyltree clustalw fengdoolittle
EXECS   = phyltree clustalw fengdoolittle
LIBRARY = APPSlibPhylo.a



#
# Install rule
#

compile: all

all: install
install: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(LIBPATH)
	mv $(EXECS) $(BINPATH)



#
# Call global Makefile to do the job
#
include ../../Makefile.global

# Overriding to use C++11
STANDARDFLAGS = -Wall -Wextra -pedantic -std=c++0x -DNEXCEPTIONS -DLINUX -c