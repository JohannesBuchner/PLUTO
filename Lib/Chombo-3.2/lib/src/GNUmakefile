# -*- Mode: Makefile -*- 

### This makefile controls building Chombo libraries.
### It is nearly identical to `$CHOMBO_HOME/test/GNUmakefile'.

### The way this makefile works is the target specified by the user
### on the command line is treated as the `action' to perform on all
### the `lib/src' subdirectories.  In the rules, the subdirectories are
### treated like targets and the action is specified by a variable.
### The makefile builds each directory target by invoking make
### recursively on the directory with the action as the recursive target.
### This assumes that the makefiles in the subdirectories implement all
### the actions in the $(_all_actions) variable.

export makefiles+=lib_src

## This variable specifies where we are.  Lower-level make files will need an
## absolute path so we use 'pwd' rather than "."
export CHOMBO_HOME := $(shell pwd)/..

##
## These variables define all the subdirectories that can be built.
## These library names are used by the makefiles in the application directories.
##[NOTE: if you change the names, change the dependencies below.]

lib_ch_targets := BaseTools BoxTools AMRTools AMRElliptic OldAMRElliptic AMRTimeDependent ParticleTools InfiniteDomain HOAMRTools MOLAMRTimeDependent HOMappedGeometry HOMappedConservation
lib_eb_targets := EBTools EBAMRTools Workshop EBAMRElliptic EBAMRTimeDependent
lib_mf_targets := MFElliptic MFTools
lib_ccse_targets := CCSE

## This variable defines the names of the subdirectories, in the order they will be build.
## One or more can be used as a target on the command line.
lib_targets := $(lib_ch_targets)
ifeq ($(USE_EB),TRUE)
  lib_targets += $(lib_eb_targets)
endif
ifeq ($(USE_MF),TRUE)
  lib_targets += $(lib_mf_targets)
endif
ifeq ($(USE_CCSE),TRUE)
  lib_targets += $(lib_ccse_targets)
endif

lib_targets := $(shell \echo $(lib_targets) `\ls` | \tr ' ' '\n' | \sort | \uniq -d)

##
## These are the actions (targets for the recursive make) supported by this
## makefile.  It assumes the makefiles in the subdirectories support them too.
##   all        builds all the libraries
##   lib        synonym for `all'
##   clean      deletes files for this configuration
##   realclean  deletes files for all configurations
##
_all_actions = all lib clean realclean include

.PHONY: $(_all_actions) $(lib_targets) usage vars

##
## Default target
##
usage:
	@echo "This makefile supports targets:"
	@echo "     $(_all_actions) $(lib_targets) usage vars"

##
## For each action, set the action variable to tell the recursive
## makefile what to do.  There should be a separate dependency for
## every value in $(_all_actions).
##[NOTE: the `NODEPENDS=TRUE' variable disables automatic
##       dependency analysis in the recursive makefile, which is only
##       needed for actions that compile or link code.]
##[NOTE: running `include' before running `lib' or `all' will not put 
##       any of the ChomboFortran header files in the include directory.
##       Caveat usor.]
##
all       : _action = all
lib       : _action = lib
clean     : _action = clean     NODEPENDS=TRUE
realclean : _action = realclean NODEPENDS=TRUE
include:    _action = include   NODEPENDS=TRUE# #[NOTE: don't move the ":"]

# this dependency ensures the action is done for every directory.
$(_all_actions) : $(lib_targets)

# this rule performs the action
$(lib_targets):
	$(MAKE) --no-print-directory --directory $@ $(_action)

# this rule prints out various 'make' variables to the screen for debugging
vars:
	@$(MAKE) --no-print-directory -f $(CHOMBO_HOME)/mk/Make.printVariables $@
