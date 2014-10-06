 
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

LIBS = -lAlign2 -lBiopool  -ltools

LIB_PATH = -L.

INC_PATH = -I. -I ../../Biopool/Sources -I../../tools/

#
# Objects and headers
#

SOURCES = Alignment.cc AlignmentBase.cc \
          Align.cc NWAlign.cc SWAlign.cc FSAlign.cc NWAlignNoTermGaps.cc \
          AlignmentData.cc SequenceData.cc SecSequenceData.cc \
          VGPFunction.cc VGPFunction2.cc \
          Substitution.cc SubMatrix.cc StructuralAlignment.cc\
          ScoringScheme.cc ScoringS2S.cc ScoringP2S.cc ScoringP2P.cc \
          PssmInput.cc Profile.cc HenikoffProfile.cc PSICProfile.cc SeqDivergenceProfile.cc \
          LogAverage.cc CrossProduct.cc DotPFreq.cc DotPOdds.cc Pearson.cc JensenShannon.cc EDistance.cc AtchleyDistance.cc AtchleyCorrelation.cc Panchenko.cc Zhou.cc \
          ThreadingInput.cc Ss2Input.cc ProfInput.cc Sec.cc Threading.cc Ss2.cc Prof.cc ThreadingSs2.cc ThreadingProf.cc  \
          ReverseScore.cc stringtools.cc

OBJECTS = Alignment.o AlignmentBase.o \
          Align.o NWAlign.o SWAlign.o FSAlign.o NWAlignNoTermGaps.o \
          AlignmentData.o SequenceData.o SecSequenceData.o \
          VGPFunction.o VGPFunction2.o \
          Substitution.o SubMatrix.o StructuralAlignment.o\
          ScoringScheme.o ScoringS2S.o ScoringP2S.o ScoringP2P.o \
          PssmInput.o Profile.o HenikoffProfile.o PSICProfile.o SeqDivergenceProfile.o \
          LogAverage.o CrossProduct.o DotPFreq.o DotPOdds.o Pearson.o JensenShannon.o EDistance.o AtchleyDistance.o AtchleyCorrelation.o Panchenko.o Zhou.o \
          ThreadingInput.o Ss2Input.o ProfInput.o Sec.o Threading.o Ss2.o Prof.o ThreadingSs2.o ThreadingProf.o  \
          ReverseScore.o stringtools.o

TARGETS =  

EXECS =  

LIBRARY = libAlign2.a

#
# Install rule
#

compile: all

install: $(LIBRARY) $(TARGETS)
	mv $(LIBRARY) $(UPDIR)/lib

#
# Call global Makefile to do the job
#
all: install
include ../../Makefile.global
