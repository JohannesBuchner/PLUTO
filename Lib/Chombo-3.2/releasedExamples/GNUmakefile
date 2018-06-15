# -*- Mode: Makefile -*- 

## This makefile builds one or more subdirectories.

# trace the chain of included makefiles
makefiles += releasedExamples

## Define the variables needed by Make.example.subdirs

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../lib

# the subdirectories to build
exec_dirs = AMRPoisson   \
            AMRGodunov   \
            AMRINS       \
            EBMenagerie  \
            EBAMRPoisson \
            EBAMRGodunov \
            EBAMRINS     \
            EBAMRCNS

# the shared code for building in example subdirectories
include $(CHOMBO_HOME)/mk/Make.example.subdirs

# application-specific variables

# application-specific targets
