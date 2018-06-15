## This file contains compiler variables for all flavors of the C++ compiler called "CC"

makefiles+=local/Make.defs.CC

# The IBM compiler can be called CC, but it's usually called xlC
ifeq ($(system),AIX)
  include $(CHOMBO_HOME)/mk/compiler/Make.defs.IBM
endif

# Check for SGI OSs (IRIX or IRIX64)
ifeq ($(findstring IRIX,$(system)),IRIX)
  ## suppress=* turns off msgs:
  ##    (1009) nested comment not allowed (this is Dave M's handiwork)
  ##    (1047) assert macro redefined
  ##    (1110) statement unreachable
  ##    (1116) non-void function should return a value
  ##    (1174) variable declared but never referenced
  ##    (1209) REMARK: The controlling expression is constant. while(1)
  ##    (1424) template parameter "X" not used in declaring argument types of function template
  ##    (1460) function is redeclared "inline" after being used
  ##    (1468) REMARK: Inline function  cannot be explicitly instantiated.
  ##    (1552) variable set but never referenced
  ##    (3201) The parameter "X" was never referenced
  basecxxflags = -64 -J4 -fullwarn -trapuv -LANG:std -LANG:ansi-for-init-scope=ON \
  -DEBUG:suppress=1009,1047,1110,1116,1169,1174,1209,1375,1424,1460,1468,1506,1552,3201,3625,3649 \
   -diag_suppress 1009,1047,1110,1116,1169,1174,1209,1375,1424,1460,1468,1506,1552,3201,3625,3649

  defcxxoptflags = -O2 -OPT:Olimit=3500 -IPA $(basecxxflags)
  defcxxdbgflags = -g $(basecxxflags)

  # Profiling under IRIX using speedshop needs no special compile/link options
  defcxxprofflags =
endif

# Check for Sun OSs
ifeq ($(findstring $(system),SunOS Solaris),$(system))
  # definitions for Sun Workshop compiler
  arch := $(filter i86pc sparc,$(unamea))
  defcxxdbgflags  = -g -ftrap=common
  defcxxprofflags = -xpg
  ifeq ($(arch),i86pc)
    defcxxoptflags = -O3 -fast
  endif
  ifeq ($(arch),sparc)
    defcxxoptflags = -O3 -dalign
  endif
  CH_AR = $(CXX) -xar -o
endif

# Check for Sun compiler on Linux.
#[NOTE: Compiler defines the macro var __SUNPRO_CC to "0x<majorver><minorver><subminver>"
ifeq ($(findstring Sun,$(shell $(CXX) -V 2>&1)),Sun)
  defcxxdbgflags = -g
  defcxxoptflags = -O
  defcxxprofflags = -xpg
endif

#  Check for Cray X1 OS (ie Unicos)
#[NOTE: this needs more work for opt and prof]
ifeq ($(system),CRAY)
  _basecxxflags = -h new_for_init -h instantiate=used
  defcxxdbgflags := -g $(_basecxxflags)
  defcxxoptflags := -O 3 -e o $(_basecxxflags)
  # Or to be even more agressive
  ## defcxxoptflags := -V -F -O inline3,scalar3,vector3,stream3,task0 -rm $(_basecxxflags)
endif
