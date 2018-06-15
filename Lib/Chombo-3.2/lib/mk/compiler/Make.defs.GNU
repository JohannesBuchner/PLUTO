### This is a compiler-specific definitions file for GNU compilers.

## It sets values for the defaults variables ($def*), which will be used
## if the calling makefile doesn't set the corresponding variables.
## The corresponding variables are the ones without the "def" prefix.

## Defaults variables:
##  defcppcomflags = C-preprocessor options for both C++ and Fortran code
##  defcppdbgflags = C-preprocessor options for both C++ and Fortran code when DEBUG=TRUE
##  defcppoptflags = C-preprocessor options for both C++ and Fortran code when OPT=TRUE or HIGH
##  defcxxcppflags = C-preprocessor options for only C++ code
##  defcxxcomflags = C++ compiler/linker options
##  defcxxdbgflags : for C++ compiling and linking when DEBUG=TRUE
##  defcxxoptflags : for C++ compiling and linking when OPT=TRUE or HIGH
##  defcxxprofflags: for C++ compiling and linking when PROFILE=TRUE
##  deffcppflags   = C-preprocessor options for only Fortran code
##  deffcomflags   = Fortran compiler options
##  deffdbgflag    : for Fortran compiling when DEBUG=TRUE
##  deffoptflag    : for Fortran compiling when OPT=TRUE or HIGH
##  deffprofflags  : for Fortran compiling when PROFILE=TRUE
##  defflibflags   : library options needed to link the Fortran code
##  defldcomflags  = linker options
##  deflddbgflags  = linker options when DEBUG=TRUE
##  defldoptflags  = linker options when OPT=TRUE or HIGH
##  defldprofflags = linker options when PROFILE=TRUE
##  cppcallsfort   = preprocessor macro as used in FORT_PROTO.H to specify
##                    how the C++ code should construct names for calling Fortran
##
## Special object models:
##  USE_64 = TRUE   Compile/link with 64bit pointers and longs
##
##  Note: 64bit only works on AMD64 or PowerPC computers with gcc 3.3 or later.  <dbs,Mar04>
##
## Notes:
##  GNU C/C++ compilers set the preprocessor macros: (use -E -dM options to display)
##    __GNUC__       (int) defined for gcc or g++, == major version number
##    __GNUG__       (int) defined for g++, == major version number
##    __GNUC_MINOR__ (int) defined for gcc or g++, == minor version number

makefiles+=compiler/Make.defs.GNU

cxxname := $(notdir $(firstword $(subst -, ,$(CXX))))
fname   := $(notdir $(firstword $(subst -, ,$(FC))))

# special case for g++ version3 installed on a system where version2 is the default
ifeq ($(cxxname),g++3)
  cxxname := g++
endif
ifeq ($(fname),g773)
  fname := g77
endif

################################################################

ifeq ($(cxxname),g++)
  _gppversion := $(subst ., ,$(shell $(CXX) -dumpversion))
  _gppmajorver:= $(word 1,$(_gppversion))
  _gppminorver:= $(word 2,$(_gppversion))

  # `-ftemplate-depth' is needed to keep mangled function names from
  #  getting too long for some linkers to handle.
  # HDF5 uses `long long', so disable warnings about it.
  # `-Wno-sign-compare' turns off warnings about comparing signed and unsigned int -
  #   this is a meaningful warning, but it appears way too often in our code.

  _cxxbaseflags := -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-25 

  # g++ v3.4.0 started objecting to some uses of HOFFSET()
  # TJL - Changed moving from csh to sh...
  ifeq (00,$(shell test $(_gppmajorver) -ge 3 ; echo $$?)$(shell test $(_gppminorver) -ge 4 ; echo $$?))
    _cxxbaseflags += -Wno-invalid-offsetof
  endif

  ifeq (0,$(shell test $(_gppmajorver) -gt 3 ; echo $$?))
    _cxxbaseflags += -Wno-variadic-macros
  endif

  # NOTE: -malign-double causes MPI programs that use pout() to
  # segfault in ~ostream(), so dont use it
  defcxxcomflags := $(_cxxbaseflags)

ifeq ($(OPENMPCC),TRUE)
  defcxxoptflags := -O2 -funroll-loops -fopenmp
  defcxxdbgflags := -g -pedantic -Wall -fopenmp
else
  defcxxoptflags := -O2 -funroll-loops -Wno-unknown-pragmas
  defcxxdbgflags := -g -pedantic -Wall -Wno-unknown-pragmas
endif
  #newer versions of gcc (4.3.2 works) can use this

  # gcc v4 doesn't understand this yet
  ifeq (0,$(shell test $(_gppmajorver) -lt 4 ; echo $$?))
    # gcc v3.3 on MacOS doesn't like it either
    ifneq ($(system),Darwin)
      defcxxdbgflags += -fbounds-checking
    endif
  endif

  # add extra symbol info to debug output for gdb
  ifeq ($(system),OSF1)
    # -gstabs breaks on gcc 3.3 and osf1 5.1 (on halem)
    # its not clear which is the cause, so for now I'll assume gcc is the problem
    ifeq ($(_gppmajorver)$(_gppminorver),33)
      defcxxdbgflags += -ggdb
    else
      defcxxdbgflags += -gstabs
    endif
  else
    defcxxdbgflags += -ggdb
  endif

  # non-default object models are meaningful ONLY on some machines
  ifeq ($(USE_64),TRUE)
    ifeq ($(arch),x86_64)
      # works on an AMD64 machine that defaults to 32bit
      # may also work on a Power{PC,4,5} but I'm not positive
      defcxxcomflags += -m64
    endif
    # on IBM Power (like seaborg)
    ### WARNING: I haven't tested this <dbs,5Mar04> !!!
    ifeq ($(system),AIX)
      defcxxcomflags += -maix64
    endif
    # Darwin Intel CPU's (Macs)
    ifeq ($(system)-$(arch),Darwin-i386)
      defcxxcomflags += -m64
    endif
    # Apple Macintosh
    ifeq ($(system)-$(arch),Darwin-powerpc)
      defcxxcomflags += -mpowerpc64
    endif
  else
    ifeq ($(arch),x86_64)
      # on AMD64 machines with default 64bit, you have to ask
      # for an emulation mode to get 32 bits
      defcxxcomflags += -m32
      defldcomflags = -Wl,-melf_i386
    endif
    # on IBM Power with default 64bit mode
    ### WARNING: I haven't tested this <dbs,5Mar04> !!!
    ifeq ($(system),AIX)
      defcxxcomflags += -maix32
    endif
    # Darwin Intel CPU's (Macs)
    ifeq ($(system)-$(arch),Darwin-i386)
      defcxxcomflags += -m32
    endif
    # Apple Macintosh
    ifeq ($(system)-$(arch),Darwin-powerpc)
      defcxxcomflags += -mno-powerpc64
    endif
  endif

  ##XXX <dbs> this is handled better in Make.defs.
  ##  # workaround for bug on MacOS where g++ and g77 dont play well
  ##  # together so force the preprocessor to use the Fortran compiler
  ifeq ($(system),Darwin)
    CH_CPP = /usr/bin/cpp -E -P
  endif


  # the preprocessor for v4.8+ adds some header gibberish that we 
  # suppress with -nostdinc.
  ifeq (00,$(shell test $(_gppmajorver) -ge 4 ; echo $$?)$(shell test $(_gppminorver) -ge 8 ; echo $$?))
    CH_CPP += -nostdinc
  endif

endif

################################################################

ifeq ($(fname),g77)
  _g77version := $(subst ., ,$(shell $(FC) -dumpspecs | awk '/^\*version:/{getline;print $$1}'))
  _g77majorver:= $(word 1,$(_g77version))
  _g77minorver:= $(word 2,$(_g77version))

  deffcomflags = -fno-second-underscore
  deffoptflags = -O2 -funroll-loops
  # if machine has i686 in the uname -a, then add the -malign-double opt flag
  ifeq ($(findstring i686,$(unamea)),i686)
    deffoptflags += -malign-double
  endif

  deffdbgflags = -g -Wimplicit -Wall -ffortran-bounds-check

  # add extra symbol info to debug output for gdb
  ifeq ($(system),OSF1)
    # -gstabs breaks on gcc 3.3 and osf1 5.1 (on halem)
    # its not clear which is the cause, so for now I'll assume gcc is the problem
    ifeq ($(_g77majorver)$(_g77minorver),33)
      deffdbgflags += -ggdb
    else
      deffdbgflags += -gstabs
    endif
  else
    deffdbgflags += -ggdb
  endif

  # non-default object models are meaningful ONLY on some machines
  ifeq ($(USE_64),TRUE)
    ifeq ($(arch),x86_64)
      # works on an AMD64 machine that defaults to 32bit
      # may also work on a Power{PC,4,5} but I'm not positive
      deffcomflags += -m64
    endif
    # Darwin Intel CPU's (Macs)
    ifeq ($(system)-$(arch),Darwin-i386)
      deffcomflags += -m64
    endif
    # on IBM Power (like seaborg)
    ### WARNING: I haven't tested this <dbs,5Mar04> !!!
    ifeq ($(system),AIX)
      deffcomflags += -maix64
    endif
  else
    ifeq ($(arch),x86_64)
      # on AMD64 machines with default 64bit, you have to ask
      # for an emulation mode to get 32 bits
      deffcomflags += -m32
    endif
    # Darwin Intel CPU's (Macs)
    ifeq ($(system)-$(arch),Darwin-i386)
      deffcomflags += -m32
    endif
    # on IBM Power with default 64bit mode
    ### WARNING: I haven't tested this <dbs,5Mar04> !!!
    ifeq ($(system),AIX)
      deffcomflags += -maix32
    endif
  endif

  # libraries needed by g77
  # g++ v2 doesn't have libfrtbegin but g77 v3 needs it
  _libg2c_dir := $(dir $(shell $(FC) -print-file-name=libg2c.a))
  ifeq ($(_libg2c_dir),)
    $(error cannot find libg2c.a using $(FC) -print-file-name=libg2c.a)
  endif
  _libfrt_dir := $(dir $(shell $(FC) -print-file-name=libfrtbegin.a))
  # set library directory (if not using g++) and library names; use libfrt only if it exists
  ifeq ($(cxxname),g++)
    defflibflags := -lg2c
  else
    defflibflags := -L$(_libg2c_dir) -lg2c
  endif
  ifneq ($(_libfrt_dir),./)  #[NOTE: $dir returns ./ when there is no directory part]
    # in the unlikely event libg2c and libfrtbegin are in different directories
    ifneq ($(_libfrt_dir),$(_libg2c_dir))
      defflibflags += -L$(_libfrt_dir)
    endif
    defflibflags += -lfrtbegin
  endif
  defflibflags += -lm
endif

###########################################################################

ifeq ($(fname),g95)
  deffcomflags = -fno-second-underscore
  deffoptflags = -O2 -funroll-loops
  # if machine has i686 in the uname -a, then add the -malign-double opt flag
  ifeq ($(findstring i686,$(unamea)),i686)
    deffoptflags += -malign-double
  endif

  defflibflags = -L/opt/local/lib -lf95
  deffdbgflags = -g -fimplicit-none -Wall -fbounds-check
endif

################################################################

# GNU Fortran95 compiler in GCC version 4
ifeq ($(fname),gfortran)
  _g77version := $(subst ., ,$(shell $(FC) -dumpspecs | awk '/^\*version:/{getline;print $$1}'))
  _g77majorver:= $(word 1,$(_g77version))
  _g77minorver:= $(word 2,$(_g77version))

  deffcomflags = -Wno-line-truncation -fno-second-underscore
  # these flags get rid of unused variables warnings. (-Wno-unused-parameter -Wno-unused-variable)
  deffcomflags += -Wno-unused-parameter -Wno-unused-variable
  deffoptflags = -O2 -funroll-loops
  deffdbgflags = -g -ggdb -Wall -W

  # Try this for gcc 4.2 and above...
  # newer versions of gcc (4.2.1 works) can use this
  ifeq (00,$(shell test $(_gppmajorver) -ge 4 ; echo $$?)$(shell test $(_gppminorver) -ge 2 ; echo $$?))
     deffdbgflags += -fbounds-check
  endif

  # Removed because this flag doesn't seem to play well with "valgrind"
  # Try this for gcc 4.3 and above...
  # newer versions of gcc (4.3.2 works) can use this
  # ifeq (00,$(shell test $(_gppmajorver) -ge 4 ; echo $$?)$(shell test $(_gppminorver) -ge 3 ; echo $$?))
  #    deffdbgflags += -fstack-check 
  # endif

  # Note this is a copy/paste from the g77 ifeq above
  # non-default object models are meaningful ONLY on some machines
  ifeq ($(USE_64),TRUE)
    ifeq ($(arch),x86_64)
      # works on an AMD64 machine that defaults to 32bit
      # may also work on a Power{PC,4,5} but I'm not positive
      deffcomflags += -m64
    endif
    # Darwin Intel CPU's (Macs)
    ifeq ($(system)-$(arch),Darwin-i386)
      deffcomflags += -m64
      ifneq ($(wildcard /usr/local/lib/x86_64), ) # This only helps if it exists!!! -JNJ
        defflibflags += -L/usr/local/lib/x86_64
      endif
    endif
    # on IBM Power
    ### WARNING: I haven't tested this <dbs,5Mar04> !!!
    ifeq ($(system),AIX)
      deffcomflags += -maix64
    endif
  else
    ifeq ($(arch),x86_64)
      # on AMD64 machines with default 64bit, you have to ask
      # for an emulation mode to get 32 bits
      deffcomflags += -m32
    endif
    # Darwin Intel CPU's (Macs)
    ifeq ($(system)-$(arch),Darwin-i386)
      deffcomflags += -m32
    endif
    # on IBM Power with default 64bit mode
    ### WARNING: I haven't tested this <dbs,5Mar04> !!!
    ifeq ($(system),AIX)
      deffcomflags += -maix32
    endif
  endif

ifeq ($(OPENMPCC),TRUE)
  defflibflags  += -lgfortran -lm -lgomp
else
  defflibflags  += -lgfortran -lm 
endif
endif

################################################################

# Special link flags for DEC/Compaq/HP Tru64/OSF1
ifeq ($(system),OSF1)
  defldcomflags = -Wl,-S
endif
