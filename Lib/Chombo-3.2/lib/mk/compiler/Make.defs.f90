## This file contains compiler variables for all flavors of the Fortran compiler called "f90"

## The compiler variables are set in Make.defs using these defaults and local overrides.

makefiles+=local/Make.defs.f90

# HP Tru64 (DEC) Fortran compiler
ifeq ($(system),OSF1)
  deffoptflags = -fast
  # if using fort/f90 without cxx, find the directory that has the Fortran libraries
  ifneq ($(cxxname),cxx)
    _ifclib := -L$(dir $(shell which $(firstword $(FC))))../lib
  endif
  defflibflags = -lfor
endif

# HP-UX Fortran compiler
ifeq ($(system),HPUX)
  deffoptflags = -fast
  cppcallsfort = -DCH_FORT_NOUNDERSCORE
endif

# SGI IRIX (MipsPro) Fortran compiler
ifeq ($(findstring IRIX,$(system)),IRIX)
  deffdbgflags = -g -DEBUG:subscript_check=ON -64
  deffoptflags = -Ofast -LNO -IPA -64
endif

# Sun Fortran compiler on Linux
ifeq ($(findstring Sun,$(shell $(FC) -V 2>&1)),Sun)
  deffdbgflags = -g
  deffoptflags = -O
  deffprofflags = -xpg
endif
