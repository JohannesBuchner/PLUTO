# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_EBAMRTools

ebase = mginterpTest aggpwlfpTest correct1d2dTest fastFRTest newIntRedistTest      \
        halfQuadTest allRegFluxRegTest aveConserveTest averageTest    \
        coarsenTest averageFluxTest pwlinterpTest fpExactTest         \
        levelRedistTest fluxRegTest fullRedistTest quadCFITestEBCross \
	halfTensorTest bigHalfQuadTest

LibNames = EBAMRTools Workshop EBTools AMRTools BoxTools

INPUT = ""

include $(CHOMBO_HOME)/mk/Make.test

CPPFLAGS += $(subst FALSE,,$(subst TRUE,-DCH_PROFILE,$(PROFILE)))
