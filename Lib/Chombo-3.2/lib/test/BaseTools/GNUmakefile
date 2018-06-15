# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_BaseTools

ebase =  clock testTask testCH_Attach testRefCountedPtr testParmParse \
   test_complex test_parstream testRootSolver

# note that BaseTools library should be included by default, even 
# if we don't specify it here
LibNames = 

include $(CHOMBO_HOME)/mk/Make.test
