# -*- Mode: Makefile -*- 

export makefiles += releaseExamples_AMRINS_exec

### This makefile produces an executable for each 
### name in the `ebase' variable

ebase := ns

#set compiler and CHOMBO_HOME

MACHINE = $(shell uname)
UNAMEM = $(shell uname -m)

DIM=2

TRAPFPE=FALSE
ifeq ($(TRAPFPE),TRUE)	
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
# XTRACPPFLAGS += -DTRAP_FPE
endif

# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
# XTRACPPFLAGS += #-DPAPI 
PAPILIBS = #/usr/local/lib/libpapi.a
PAPILIBS = #/usr/common/usg/papi/2.1/lib/libpapi.a -L/usr/lpp/pmtoolkit/lib -lpmapi

#  have a default for machines not defined here.
# make this a relative path so it will find the Chombo 
# on same level as dMartin.
CHOMBO_HOME = ../../../lib

#this is for the NERSC SP (gseaborg)
#ifeq ($(MACHINE),AIX)
#CHOMBO_HOME = $(HOME)/Chombo/lib
#CHOMBO_HOME = /u2/dmartin/code/Chombo/lib
#CHOMBO_HOME = /usr/common/homes/n/ndk/Chombo/lib
#endif

include $(CHOMBO_HOME)/mk/Make.defs
include $(CHOMBO_HOME)/mk/Make.defs.config
XTRALIBFLAGS += $(PAPILIBS)

#cases := $(shell echo case*)

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames :=  AMRElliptic AMRTools AMRTimeDependent BoxTools 

base_dir = .

AMRINS_BASE_DIR = ../

EXEC_DIR = .

PROJ_DIR = $(AMRINS_BASE_DIR)/projection
src_dirs = $(PROJ_DIR)

UTIL_DIR = $(AMRINS_BASE_DIR)/util
src_dirs += $(UTIL_DIR)

BCUTIL_DIR = $(AMRINS_BASE_DIR)/BCutil
src_dirs += $(BCUTIL_DIR)

NS_DIR = ../NavierStokes
src_dirs += $(NS_DIR)

ADV_DIR = $(AMRINS_BASE_DIR)/advection
src_dirs += $(ADV_DIR)

INPUT = inputs.ring

include $(CHOMBO_HOME)/mk/Make.example

ifeq ($(DEBUG), TRUE)
# petermc, 12 Oct 2007:  never set XTRACPPFLAGS in makefile
#XTRACPPFLAGS += -DDEBUG
endif
