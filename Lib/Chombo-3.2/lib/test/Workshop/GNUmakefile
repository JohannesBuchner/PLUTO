# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_Workshop

ebase = divergeTest pointCoarseningTest ldBaseIFFABTest cylinderTest coarseningTest fabTestTwo   \
        impFuncTest iffabExchangeTest linearizationTest normTest \
        rampTest sphereConvTest sphereTest eieioTest irregFABArith ebisWriteAllTest

LibNames = Workshop EBAMRTools EBTools AMRTools BoxTools

include $(CHOMBO_HOME)/mk/Make.test
