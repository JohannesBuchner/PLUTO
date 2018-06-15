## This file contains compiler variables for all flavors of the Fortran compiler called "f77"

## The compiler variables are set in Make.defs using these defaults and local overrides.

## The important things to know:
##  1) this file should _only_ set variables that begin with "def".
##  2) the $cxxname and $fname variables should be used to determine
##     the names of the C++ and Fortran compilers, since $CXX and $FC.
##     may contain a directory path and/or options.
##  3) The 'Make.defs' file automatically merges the $def* variables into
##     the variables actually used in the make rules that compile/link code.
##  4) The user can override any of these variables on the command line.
##  5) The $cxx*flags variables (except for $cxxcppflags) are also used in
##     the link rule, so any option that is the same for the linker and
##     compiler should not be repeated in the $ld*flags variables.

## Compiler Variables:
##  defcppcomflags = C-preprocessor options for both C++ and Fortran code
##  defcppdbgflags = C-preprocessor options for both C++ and Fortran code when DEBUG=TRUE
##  defcppoptflags = C-preprocessor options for both C++ and Fortran code when OPT=TRUE or HIGH
##  defcxxcppflags = C-preprocessor options for only C++ code
##  defcxxcomflags = C++ compiler/linker options
##  defcxxdbgflags = C++ compiler/linker options when DEBUG=TRUE
##  defcxxoptflags = C++ compiler/linker options when OPT=TRUE or HIGH
##  defcxxprofflags= C++ compiler/linker options when PROFILE=TRUE
##  deffcppflags   = C-preprocessor options for only Fortran code
##  deffcomflags   = Fortran compiler options
##  deffdbgflags   = Fortran compiler options when DEBUG=TRUE
##  deffoptflags   = Fortran compiler options when OPT=TRUE or HIGH
##  deffprofflags  = Fortran compiler options when PROFILE=TRUE
##  defflibflags   = linker options to specify the Fortran libraries
##  defldcomflags  = linker options
##  deflddbgflags  = linker options when DEBUG=TRUE
##  defldoptflags  = linker options when OPT=TRUE or HIGH
##  defldprofflags = linker options when PROFILE=TRUE

makefiles+=local/Make.defs.f77

# HP Tru64 (DEC) Fortran compiler
ifeq ($(system),OSF1)
  deffoptflags = -fast
  # if using f77 without cxx, find the directory that has the Fortran libraries
  ifneq ($(cxxname),cxx)
    _ifclib := -L$(dir $(shell which $(firstword $(FC))))../lib
  endif
  defflibflags = $(_ifclib) -lfor
endif

# HP-UX Fortran compiler
ifeq ($(system),HPUX)
  deffoptflags = -fast
  cppcallsfort = -DCH_FORT_NOUNDERSCORE
endif

# SGI IRIX (MipsPro) Fortran compiler
ifeq ($(findstring IRIX,$(system)),IRIX)
  deffcomflags := -64
  deffdbgflags := -g -DEBUG:subscript_check=ON
  deffoptflags := -Ofast -LNO -IPA
endif

# Sun Fortran compiler
ifeq ($(findstring $(system),Solaris SunOS),$(system))
  deffdbgflags = -g -C -errtags
  ifeq ($(arch),i86pc)
    deffoptflags = -O4 -fast -pentium
  endif
  ifeq ($(arch),sparc)
    deffoptflags = -O4 -fast -dalign
  endif
  deffprofflags = -xpg
  # if using f77 without CC, find the directory that has the Fortran libraries
  ifneq ($(cxxname),CC)
    _ifclib := -L$(dir $(shell which $(firstword $(FC))))../lib
  endif
  defflibflags = $(_ifclib) -lF77 -lM77 -lsunmath
endif

# in the unlikely event someone is using 'f77' when they mean 'g77'
ifeq ($(system),Linux)
  #[NOTE: 'override' is dangerous because it can create hard-to-find bugs elsewhere,
  #       but Make.defs.GNU computes $fname for itself, so there's no choice.]
  override fname := g77
  # $warning appeared in 3.78 (or 3.79??)
  _makemajorver := $(word 1,$(subst ., ,$(MAKE_VERSION)))
  _makeminorver := $(word 2,$(subst ., ,$(MAKE_VERSION)))
  # TJL - Changed moving from csh to sh...
  # ifeq (3-0,$(_makemajorver)-$(shell test $(_makeminorver) -gt 77 ; echo $$status))
  ifeq (3-0,$(_makemajorver)-$(shell test $(_makeminorver) -gt 77 ; echo $$?))
    $(warning assuming f77 command on Linux is really GNU g77)
  endif
  include $(CHOMBO_HOME)/mk/compiler/Make.defs.GNU
endif
