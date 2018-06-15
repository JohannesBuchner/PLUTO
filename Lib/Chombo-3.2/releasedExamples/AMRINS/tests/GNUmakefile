# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_AMRINS_tests

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = testGradient

# the location of the Chombo "lib" directory
CHOMBO_HOME = $(shell pwd)/../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTimeDependent AMRTools  BoxTools 

# the locations of the source code directories
base_dir = .
src_dirs = ../util ../BCutil

# input file for 'run' target
INPUT = inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

NEW_GRADIENT = TRUE
ifeq ($(NEW_GRADIENT),TRUE)
XTRACPPFLAGS += -DNEW_GRADIENT
endif

# Note: tell Make.rules to not use -q for 'run' target
_space = # # NOTE: keep the space after the "=" or bad things happen
VERBOSE = $(_space)

# application-specific targets
