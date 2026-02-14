#--*- makefile -*--------------------------------------------------------------
#
#   Standard makefile
#
#------------------------------------------------------------------------------

# Path to project directory.
UPDIR = .
# Path to subdirectories.
SUBDIRS =  tools Energy/Sources   Biopool/Sources  Align2/Sources Phylo/Sources Energy/Sources/TorsionPotential  Lobo/Sources Lobo/APPS Energy/APPS   Biopool/APPS   Align2/APPS Phylo/APPS

#
# Libraries and paths (which are not defined globally).
#

LIBS =

LIB_PATH = lib

INC_PATH =

BINPATH = bin

#
# Objects and headers
#

SOURCES =
OBJECTS =
EXECS =
TARGET  =
LIBRARY =

#
# Install rule
#

start: subdirs

compile: all

clean: subclean

install: subinstall
	
test: all
	

depend: subdepend

#
# Call global Makefile to do the job.
#
export VICTOR_ROOT=$(dir $(patsubst %/,%, $(shell pwd)))
include Makefile.global






