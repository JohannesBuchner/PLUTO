# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_AMRGodunov_execPolytropic

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = amrGodunov

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTimeDependent AMRTools BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../srcPolytropic 

# input file for 'run' target
INPUT = ramp.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# Note: tell Make.rules to not use -q for 'run' target
_space = # # NOTE: keep the space after the "=" or bad things happen
VERBOSE = $(_space)

# application-specific targets
