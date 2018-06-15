# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase = compare

## location of 'lib' directory
CHOMBO_HOME = ../..

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames =  AMRTools AMRTimeDependent BoxTools


# input file for 'run' target
INPUT = inputs.compare
override VERBOSE = # #program doesnt handle -q option


include $(CHOMBO_HOME)/mk/Make.test
