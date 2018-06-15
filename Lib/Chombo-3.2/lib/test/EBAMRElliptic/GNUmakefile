# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_EBAMRElliptic

ebase := testDirVTEBBC testRelaxEB testBCGEB poissonHeatTest

LibNames := EBAMRElliptic AMRElliptic EBAMRTimeDependent EBAMRTools Workshop EBTools AMRTimeDependent AMRTools BoxTools

include $(CHOMBO_HOME)/mk/Make.test
