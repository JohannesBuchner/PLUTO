## This file defines variables for use on the login nodes of the ORNL XT3/4
## machine 'jaguar'.  
##
## NOTE: everything is always in 64bit mode

makefiles+=local/Make.defs.jaguar

CXX=CC
FC=ftn
MPICXX=CC -target=linux
USE_64=TRUE
USE_HDF=TRUE
MPI=TRUE
DEBUG=FALSE 
OPT=HIGH

CH_CPP=$(CXX) -E

# Add -DCH_DISABLE_SIGNALS to compile line for use on Catamount
#cxxcomflags += -DCH_DISABLE_SIGNALS
#cppdbgflags += -DCH_DISABLE_SIGNALS
#cxxdbgflags += -g -DCH_DISABLE_SIGNALS 
cxxoptflags += -march=barcelona -ffast-math -O3
foptflags += -O2
# the pgf libs are needed for linking parallel HDF5
flibflags += -lgfortran 

HDFINCFLAGS=-I$(HDF5_DIR)/include 
HDFMPIINCFLAGS=-I$(HDF5_DIR)/include 
HDFLIBFLAGS=-L$(HDF5_DIR)/lib
HDFMPILIBFLAGS=-L$(HDF5)/lib

