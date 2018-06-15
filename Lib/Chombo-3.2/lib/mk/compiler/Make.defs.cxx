### This is a compiler-specific definitions file for the HP Tru64 C++ compiler "cxx"
## It documents all the compiler variables.

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
##  cppcallsfort   = preprocessor macro as used in FORT_PROTO.H to specify
##                    how the C++ code should construct names for calling Fortran

makefiles+=compiler/Make.defs.cxx

cxxname := $(notdir $(firstword $(CXX)))
fname   := $(notdir $(firstword $(FC)))

ifeq ($(cxxname),cxx)
  defcxxcppflags = -D__USE_STD_IOSTREAM
  defcxxcomflags = -std gnu
  defcxxdbgflags = -g
  defcxxoptflags = -fast
  #[NOTE: these may be needed on some versions of cxx, but not v6.5 on OSF1 v5.1]
  #defldcomflags = -Wl,-S
  syslibflags   = -lm
  # -C disables removing C++ comments starting with "//" 
  # since this is valid code in Fortran
  CH_CPP = $(CXX) -E -C -x c++
endif
