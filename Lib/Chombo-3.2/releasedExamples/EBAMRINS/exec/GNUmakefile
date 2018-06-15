# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_EBAMRINS_exec

# This is an embedded boundary example
USE_EB=TRUE

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = viscousDriver

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRTimeDependent EBAMRElliptic  EBAMRTools Workshop EBTools AMRTimeDependent AMRElliptic AMRTools BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../src ../srcLib

# input file for 'run' target
INPUT = sphere.inputs

include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
