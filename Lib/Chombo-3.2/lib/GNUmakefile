
### This is the top-level Chombo GNUmakefile.  It requires GNU make.
### It builds some or all of the libraries and tests directories, 
### depending on the targets specified on the command line.
###
### User-definable variables are in UPPER case.
### Internal variables are in lower case.
### File-local variables are in lower case with a leading "_".
###
### The default target is "all".
###
### The main CONFIGURATION variables are:
###   DEBUG = {TRUE,FALSE} compile with debugging (TRUE) or optimization (FALSE)
###   DIM = {2,3} number of dimensions to build the code with (i.e. CH_SPACEDIM macro)
###   PRECISION = {DOUBLE,FLOAT} floating point data type to build with
###   PROFILE = {TRUE,FALSE} profiling is enabled (TRUE) or disabled (FALSE)
###   MPI = {TRUE,FALSE} build serial code (FALSE) or parallel code with MPI (TRUE)
###                      (defines 'MPI' cpp macro)
###   CXX = C++ compiler (the serial C++ compiler when $MPI=TRUE)
###   MPICXX = C++ compiler to use when $MPI=TRUE
###            (in which case CXX is the serial compiler used by MPICXX)
###   FC = Fortran compiler
###   USE_64 = {FALSE,TRUE} if TRUE, compile for 64bit points and add ".64" to filenames
###   XTRACONFIG = additional string to add to filenames (used to build separate
###                versions of Chombo that have the same config. variables but
###                different values of other variables (see below))
### Other variables that control what gets built but don't affect the filenames: 
###   USE_EB     = {TRUE,FALSE} if TRUE, compile EBChombo code (default: FALSE)
###   USE_HDF    = {TRUE,FALSE} if TRUE, compile HDF5 code (default: TRUE)
###                                     (defines 'HDF5' cpp macro)
###   USE_FFTW   = {TRUE,FALSE} if TRUE, compile InfiniteDomain using FFTW library (default: FALSE)
###                                     (requires locally-set FFTWDIR variable)
###   USE_MF     = {TRUE,FALSE} if TRUE, compile MultiFluid code (default: FALSE)
###   USE_MT     = {TRUE,FALSE} if TRUE, compile Chombo memory tracking code (default: TRUE)
###                                     (defines 'CH_USE_MEMORY_TRACKING' cpp macro)
###   USE_TIMER  = {TRUE,FALSE} if TRUE, compile with Chombo Timers enabled (default: TRUE)
###   USE_PYTHON = {TRUE,FALSE} if TRUE, compile PyParse code (default: FALSE)
###   DIST_EB    = {TRUE,FALSE} if TRUE, put EBChombo code in tar file (default: FALSE)
###   DIST_MF    = {TRUE,FALSE} if TRUE, put MultiFluid code in tar file (default: FALSE)
###
### Other variables that can be set on the gmake command line:
###   CPPFLAGS = options to pass to the C preprocessor (for both C++ and Fortran code)
###   CXXFLAGS = options to pass to the C++ compiler
###   FFLAGS = options to pass to the Fortran compiler
###   LDFLAGS = options to pass to the linker
###   LIBFLAGS = options to pass to the linker to specify libraries
###   ECHO   controls verbosity; default is "@", change to "" (null string)
###          to echo commands that _aren't_ echoed by default
###   QUIET  controls verbosity; default is "", change to "@" to disable
###          echo of commands that _are_ echoed by default
###
###  For each of the above, there is another variable `XTRA*', which instead of
###  replacing the default value of the variable, adds onto the end of it.  
###  I.e: XTRACPPFLAGS XTRACXXFLAGS XTRAFFLAGS XTRALDFLAGS XTRALIBFLAGS
###
### If you dont like the default values, you can reset them on the command line
###  or in the file "mk/Make.defs.local".
###
### Usage:
###   target:    purpose:
###   -------    ----------
###   all        build all libraries and tests for this configuration
###   check      check the user-customizable variables are correct
###   clean      delete all built files for this configuration
###   defs       print values of user-customizable variables
###   dist       build a tar file of the whole dir. tree (does a "realclean" first)
###   dist-noclean same as "dist" except doesn't do a "realclean" first
###   doc        build the documentation directories (config vars are ignored)
###   example    build the example programs in $(CHOMBO_HOME)/example
###   include    copies all the Chombo .H files to $(CHOMBO_HOME)/include
###   lib        build only the libraries for the specified configuration
###   realclean  delete everything that wasn't in the original distribution
###   realclean-example   same as 'realclean' but for 'example' dirs
###   run        build and run the test programs (see run_targets)
###   run-tests  only runs the test programs, doesn't rebuild anything (see run_targets)
###   setup      checks to see if there is a Make.defs.local file
###   tags       build TAGS file for emacs
###   test       build only the test programs (will also build necessary libraries)
###   usage      list the available targets and variables
###   vars       print variable values for all variables for debugging make files
###   <library-name> build only this library and its dependencies (see lib_targets)
###   <test-name>    build only this test and its dependencies (see test_targets)
###
### Notes:
### 1) this file _does_not_ include the 'Make.rules' file.
###
##################################################################################

export makefiles+=lib

#
# This variable specifies where we are.  Lower-level make files will need an
# absolute path so we use 'pwd' rather than "."
#
export CHOMBO_HOME = $(shell pwd)

##
## Default values of configuration variables.
##
include $(CHOMBO_HOME)/mk/Make.defs.defaults

## Name of EBChombo repository to use
EBCHOMBO = EBChomboRel

##
## These variables define all the subdirectories that can be built.
## These library names are used by the makefiles in the application directories.
##[NOTE: if you change the names, change the dependencies below.]

ifeq ($(DIM),1)
lib_ch_targets := BaseTools BoxTools AMRTools AMRElliptic OldAMRElliptic AMRTimeDependent ParticleTools HOAMRTools MOLAMRTimeDependent HOMappedGeometry HOMappedConservation 
extra_test_targets := DotProduct ChomboFortran
endif
ifeq ($(DIM),2)
lib_ch_targets := BaseTools BoxTools AMRTools AMRElliptic OldAMRElliptic AMRTimeDependent ParticleTools HOAMRTools MOLAMRTimeDependent HOMappedGeometry HOMappedConservation
extra_test_targets := DotProduct ChomboFortran
endif
ifeq ($(DIM),3)
lib_ch_targets := BaseTools BoxTools AMRTools AMRElliptic OldAMRElliptic AMRTimeDependent ParticleTools InfiniteDomain HOAMRTools MOLAMRTimeDependent HOMappedGeometry HOMappedConservation
extra_test_targets := DotProduct ChomboFortran
endif
ifeq ($(DIM),4)
lib_ch_targets := BaseTools BoxTools AMRTools AMRTimeDependent HOAMRTools HOMappedGeometry HOMappedConservation
extra_test_targets := ChomboFortran
endif
ifeq ($(DIM),5)
lib_ch_targets := BaseTools BoxTools AMRTools AMRTimeDependent HOAMRTools HOMappedGeometry HOMappedConservation
extra_test_targets := ChomboFortran
endif
ifeq ($(DIM),6)
lib_ch_targets := BaseTools BoxTools AMRTools AMRTimeDependent HOAMRTools HOMappedGeometry HOMappedConservation
extra_test_targets := ChomboFortran
endif

lib_eb_targets := EBTools EBAMRTools Workshop EBAMRElliptic EBAMRTimeDependent
lib_new_eb_targets := NewEBTools NewEBAMRTools NewWorkshop NewEBAMRElliptic NewEBAMRTimeDependent
lib_mf_targets := MFElliptic MFTools
lib_ccse_targets := CCSE

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
ifeq ($(USE_PYTHON),TRUE)
  lib_targets += PyParse
endif

lib_targets := $(shell \echo $(lib_targets) `\ls src` | \tr ' ' '\n' | \sort | \uniq -d)

# some src subdirs don't have tests; ChomboFortran and DotProduct have tests
# but no src subdir
#[NOTE: this has to be "=", not ":=" for it to work correctly]
test_targets = $(addprefix test/,$(filter-out ArrayView EBArrayView InfiniteDomain HOAMRTools MOLAMRTimeDependent,$(lib_targets)) $(extra_test_targets))

#filter out higher-dimension tests which don't work, even though the libraries compile
ifeq ($(DIM),4)
test_targets := $(filter-out test/HOMappedGeometry test/HOMappedConservation ,$(test_targets))
endif
ifeq ($(DIM),5)
test_targets := $(filter-out test/HOMappedGeometry test/HOMappedConservation ,$(test_targets))
endif
ifeq ($(DIM),6)
test_targets := $(filter-out test/HOMappedGeometry test/HOMappedConservation ,$(test_targets))
endif

run_targets  = $(subst test/,run/,$(test_targets))

##
## These variables are needed by the `dist' targets
##
pkg_id := Chombo-3.2
ifeq ($(DIST_EB),TRUE)
  pkg_id := EB$(pkg_id)
endif
today  := $(shell date +%d%b%Y)

###
### End of variables.  Beginning of targets.
###

.PHONY: usage all lib lib-only test defs doc example example-only tags clean realclean realclean-example \
  check dist distclean dist-noclean doxygen setup vars .vars \
  $(test_targets) $(run_targets) $(lib_targets) run run-tests run-only include regression


# This is the default target.  It prints out the supported targets.
usage:
	@echo "This GNUmakefile supports the following targets and variables:"
	@echo "     all lib test run run-tests example example-only clean realclean doc tags defs vars .vars dist"
	@echo "     $(lib_targets)"
	@echo "     $(test_targets)"
	@echo "     $(run_targets)"
	@echo "     DIM={2,3} DEBUG={TRUE,FALSE} OPT={TRUE,FALSE,HIGH} PRECISION={DOUBLE,FLOAT} PROFILE={FALSE,TRUE} MPI={FALSE,TRUE}"
	@echo "     CXX= FC= MPICXX= XTRACONFIG= [XTRA]{CXX,CPP,F,LD,LIB}FLAGS="
	@echo "     USE_EB={FALSE,TRUE} USE_MF={FALSE,TRUE} USE_MT={FALSE,TRUE} USE_64={FALSE,TRUE} USE_HDF={FALSE,TRUE} USE_FFTW={FALSE,TRUE}"


# This is the primary target that most users will use.

all: .check libtest

test: .check libtest

libtest:
	+$(ECHO)$(MAKE) --no-print-directory lib
	+$(ECHO)$(MAKE) --no-print-directory $(test_targets)

# These targets are included for convenience.  They just invoke other targets.

#NOTE: the include is done recursively to allow parallel makes
lib: .check $(lib_targets)
	+@$(MAKE) --no-print-directory include

run: .check testrun

testrun:
	+$(ECHO)$(MAKE) --no-print-directory test
	+$(ECHO)$(MAKE) --no-print-directory $(run_targets)

run-tests run-only: 
	+$(ECHO)$(MAKE) --no-print-directory $(run_targets)

example: .check
	+$(MAKE) --no-print-directory --directory=../example all

example-only: .check
	+$(MAKE) --no-print-directory --directory=../example example-only

MF: .check $(lib_ch_targets) $(lib_eb_targets) $(lib_mf_targets) $(test_mf_targets) include

# This target can be run after unpacking the distribution tar file
setup:
	@if [ ! -e mk/Make.defs.local ]; then cp mk/Make.defs.local.template mk/Make.defs.local; fi

# This target should be run after the 'run' target. It prints a summary of all *.prof files
profile:
	$(ECHO)find . -name \*$(config).prof -print -exec awk '$$1~/%/||$$1+0>=5{print;next}/time/{print;s=1;next}s==1{print;s=2;next}s==2{if(NF>2){t=$$2;next}printf "total: %9.2f\n\n",t;exit}' {} \;

# This target copies all the header files to the `include' directory.
#[NOTE: this target also appears in Make.rules.]

ifeq ($(MULTIDIM),TRUE) 
MULTIDIM_INCLUDE_DIR = $(CHOMBO_HOME)/include/multidim
endif

include: $(CHOMBO_HOME)/include  $(MULTIDIM_INCLUDE_DIR)
ifneq ($(wildcard $(CHOMBO_HOME)/include/*.H),)
	$(QUIET)chmod -R u+w $(CHOMBO_HOME)/include
endif
	$(QUIET)$(foreach d,$(lib_targets),cp -p $(CHOMBO_HOME)/src/$(d)/*.H $(CHOMBO_HOME)/include ;)
#build multidim headers if necessary -- hardwire to build DIM=1->6
ifeq ($(MULTIDIM),TRUE)
	$(QUIET)cat $(foreach d,$(lib_targets), $(CHOMBO_HOME)/src/$(d)/multidim/dim-independent-headers.txt)> $(CHOMBO_HOME)/include/multidim/dim-independent-headers.txt
	$(QUIET)$(CHOMBO_HOME)/util/multidim/make_headers.sh $(CHOMBO_HOME)/include $(CHOMBO_HOME)/include 1 6
	$(QUIET)cp -p $(CHOMBO_HOME)/src/MultiDim/*.H.transdim $(CHOMBO_HOME)/include/multidim
	$(QUIET)find $(CHOMBO_HOME)/include/multidim -name '*.H*' -exec chmod u-w {} \;
endif
	$(QUIET)find $(CHOMBO_HOME)/include -name '*.H' -exec chmod u-w {} \;

$(CHOMBO_HOME)/include:
	$(QUIET)mkdir -p $@

$(MULTIDIM_INCLUDE_DIR) : $(CHOMBO_HOME)/include
	$(QUIET)mkdir -p $@

# This target builds an emacs TAGS file for all the lib,test sources.
#[Note: etags is very fast so just redo the whole file, don't bother with dependencies.]
TAGS tags:
	-etags {src,util,test}/*/*.{H,cpp} -l fortran {src,util,test}/*/*.{ChF,F,f}

# This target removes all non-essential files from the Chombo directory tree for this configuration
clean:
	$(MAKE) --no-print-directory --directory=src clean
	$(MAKE) --no-print-directory --directory=test clean

# This target removes all non-essential files from the Chombo directory tree for all configurations
#[NOTE: the "&&" after `cd include` prevents deleting all *.H files if subdir "include" doesn't exist. <dbs>]
realclean:
	-cd include && chmod -R u+w . && find . -name "*.H" -exec $(RM) {} \;
ifneq ($(wildcard include/multidim),)
	-cd include/multidim && find . -name "*.H.*" -exec $(RM) {} \;	
	-cd include/multidim && find . -name "dim*.txt" -exec $(RM) {} \;
endif
	-find . -depth -type d \( -name d -o -name p -o -name f -o -name o \) -exec $(RM) -r {} \;
	-find . \( -name \*.ex -o -name \*.a -o -name \*_F.H -o -name \*.hdf5 -o -name \*.h5 -o -name pout.\[0-9\]\* -o -name \*PROF.prof \) -a \( ! -path './src/CCSE/*' \) -exec $(RM) {} \;
	-find . \( -name \*.tmp\? -o -name core -o -name \*.stackdump -o -name \*~ -o -name work.pc\* -o -name gmon.out -o -name leak.out\* -o -name .\#\* \) -exec $(RM) {} \;
	-find . -depth  \( -name .check.*  \) -exec $(RM) -r {} \;
	-$(RM)  .check
	-$(RM) installed_multidim.flag

realclean-example:
	$(MAKE) --no-print-directory --directory=../example realclean

distclean: realclean realclean-example
	-$(RM) -r doc/doxygen/latex doc/doxygen/html TAGS

# These targets build tar files of the Chombo directory. `dist' does a `realclean' first.
dist: distclean dist-noclean

# compute the name of the Chombo root directory, in case it isn't "Chombo"
#[NOTE: this uses '=' intead of ':=' so the shell commands are only run when needed.]
_dir_words = $(subst /, ,$(dir $(shell pwd)))
chombo_dir = $(word $(words $(_dir_words)),$(_dir_words))# last word is the one we want

dist-nodoc: distclean tarball

dist-noclean: doc tags tarball

regression:
	util/regression/regression.sh $(shell pwd)/util/regression ${DIM} ${DEBUG} ${FC} ${CXX} ${OPT} $(shell pwd)/../..
	mv util/regression/*log.regression.* regression 
	@echo "Regression log files are put into regression subdirectory"
	@echo "Check regression/status.log.<datestamp> for status information."
	@echo "If a status != 0, then that example is broken."

tarball:
	echo "$(pkg_id)" > ../VERSION ; echo "$(today)" >> ../VERSION
	sed s/Chombo/$(chombo_dir)/ mk/dist.excludes    > /tmp/dist.excludes.tmp
ifeq ($(DIST_EB),TRUE)
	$(RM) /tmp/dist.excludes.tmp
	sed s/Chombo/$(chombo_dir)/ mk/dist.excludes.EB > /tmp/dist.excludes.tmp
endif
ifeq ($(DIST_MF),TRUE)
	$(RM) /tmp/dist.excludes.tmp
	sed s/Chombo/$(chombo_dir)/ mk/dist.excludes.MF > /tmp/dist.excludes.tmp
endif
	cd ../.. ; tar cvfX $(pkg_id)-$(today).tar /tmp/dist.excludes.tmp --exclude=CVS --exclude=\*~ $(chombo_dir)
	cd ../.. ; gzip $(pkg_id)-$(today).tar
	$(RM) /tmp/dist.excludes.tmp

dist-regression: realclean
ifeq ($(DIST_EB),TRUE)
	sed s/Chombo/$(chombo_dir)/ mk/dist.excludes.EB > /tmp/dist.excludes.tmp
else
	sed s/Chombo/$(chombo_dir)/ mk/dist.excludes    > /tmp/dist.excludes.tmp
endif
	cd ../.. ; tar cfX regression.tar /tmp/dist.excludes.tmp $(chombo_dir)
	cd ../.. ; gzip regression.tar
	mv ../../regression.tar.gz .
	$(RM) /tmp/dist.excludes.tmp

# These targets build the HTML documentation from the code.
doc doxygen:
	mkdir -p doc/doxygen
ifeq ($(USE_MF),TRUE)
	$(DOXYGEN) doxygen.config.mf | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)'
else
ifeq ($(USE_EB),TRUE)
	$(DOXYGEN) doxygen.config.eb | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)'
else
	$(DOXYGEN) doxygen.config    | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)'
endif
endif
	@echo "" ; echo "point browser at $(CHOMBO_HOME)/doc/doxygen/html/index.html" ; echo ""

# Special target to run some basic checks that the customization variables are correct.
# The '.check' file is a real file, not a phony target.
# It depends on the Make.defs.local customization file so the checks will be rerun
# whenever the customizations change.
##[NOTE: it would be better if it depended on the Make.defs.local.<hostname> file too,
##       but that doesn't have to exist, so it can't depend on it. <dbs>]
check: .check

.check: mk/Make.defs.local
	@$(MAKE) --no-print-directory --directory $(CHOMBO_HOME) -f mk/Make.check check
	@touch .check

# print out a helpful error message if Make.defs.local isn't found and exit
mk/Make.defs.local:
	@echo ""
	@echo "You do not appear to have a 'Make.defs.local' file in lib/mk"
	@echo "This is where the Chombo build system looks for configuration"
	@echo "information specific to your system. See the documentation in"
	@echo "doc/chomboDesign.pdf."
	@echo ""
	@echo "There is a 'Make.defs.local.template' file available in lib/mk"
	@echo "which can be used as starting point for customization, and there"
	@echo "are tested configurations specific to various commonly used"
	@echo "machines (hopper, jaguar, etc) in lib/mk/local"
	@echo ""
	@exit 1

# This target prints out the configuration and compiler make variables,
#  i.e. the ones an ordinary user is likely to customize.
# It is intended for debugging.
defs:
	@$(MAKE) --no-print-directory -f $(CHOMBO_HOME)/mk/Make.printVariables defs

# These targets print out all the 'make' variables to the screen for debugging.
vars:
	@echo "pkg_id=$(pkg_id)"
	@echo "date=$(date)"
	@echo "today=$(today)"
	@echo "lib_targets=$(lib_targets)"
	@echo "test_targets=$(test_targets)"
	@echo "run_targets=$(run_targets)"
	@$(MAKE) --no-print-directory -f $(CHOMBO_HOME)/mk/Make.printVariables $@

.vars:
	@echo ".VARIABLES=$(.VARIABLES)"
##
## These targets are used by lower-level makefiles to build libraries by
## recusively invoking this makefile.  They just map into the mixed case
## names in $(lib_targets).  By the time they need to do this, the lower
## level makefiles don't remember the mixed case names.
##

amrelliptic:          AMRElliptic
amrnodeelliptic:      AMRNodeElliptic
amrtimedependent:     AMRTimeDependent
amrtools:             AMRTools
arrayview:            ArrayView
basetools:            BaseTools
boxtools:             BoxTools
ccse:                 CCSE
ebamrelliptic:        EBAMRElliptic
ebamrtimedependent:   EBAMRTimeDependent
ebamrtools:           EBAMRTools
ebarrayview:          EBArrayView
ebtools:              EBTools
oldamrelliptic:       OldAMRElliptic
infinitedomain:       InfiniteDomain
#multiblock:           MultiBlock
mfelliptic:           MFElliptic
mftools:              MFTools
molamrtimedependent:  MOLAMRTimeDependent
particletools:        ParticleTools
pic:                  PIC
workshop:             Workshop
homappedgeometry:     HOMappedGeometry
homappedconservation: HOMappedConservation
hoamrtools:           HOAMRTools                           


# These targets do the real work.
$(lib_targets):
	+$(make-Lib)

$(test_targets):
	+$(make-Test)

$(run_targets):
	+$(make-Run)

##
## Macro definitions used by targets.
##
define make-Lib
$(MAKE) --no-print-directory --directory=$(CHOMBO_HOME)/src/$(@F)  lib
endef

define make-Test
$(MAKE) --no-print-directory --directory=$(CHOMBO_HOME)/test/$(@F) test
endef

define make-Run
$(MAKE) --no-print-directory --directory=$(CHOMBO_HOME)/test/$(@F) run
endef
