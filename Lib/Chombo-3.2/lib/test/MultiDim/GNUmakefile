# -*- Mode: Makefile; -*-

# the location of Chombo lib dir
CHOMBO_HOME = ../../

ebase = sliceTest

MINDIM = 1
MAXDIM = 3

# this is the local GNUmakefile which contains this example's multidim
# build info -- libraries and source directory information
MULTIDIM_MAKEFILE = GNUmakefile.multidim

# this is the default target, which just invokes the all-multidim target
all: all-multidim

# this file contains the basic rules used to build multidim codes (using the 
# GNUmakefile.multidim in this directory), including the shell script
# which orchestrates the make process 
include $(CHOMBO_HOME)/mk/Make.multidim.basic

