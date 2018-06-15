# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase = removeGhosts

## location of 'lib' directory
CHOMBO_HOME = ../..

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames = AMRTools BoxTools

override VERBOSE = # #program doesnt handle -q option

include $(CHOMBO_HOME)/mk/Make.test
