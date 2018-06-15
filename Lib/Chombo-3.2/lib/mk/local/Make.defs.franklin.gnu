## This file defines variables for use on the login nodes of the NERSC Linux
## machine 'franklin'.  
##
## NOTE: everything is always in 64bit mode

makefiles+=local/Make.defs.franklin.gnu

CXX=CC
FC=ftn
MPICXX=CC
#MPICXX=CC -target=linux
USE_64=TRUE

CH_CPP=$(CXX) -march=barcelona -E

cxxoptflags += -march=barcelona -ffast-math -O3
foptflags += -O2
# the pgf libs are needed for linking parallel HDF5
flibflags += -lgfortran -L/opt/pgi/default/linux86-64/default/lib -lpgf90 -lpgf90rtl -lpgftnrtl -lpgc -lpghpf_mpi -lpghpf2

# The appropriate 'module' must be loaded for this to work.
# For serial, do    'module load hdf5'
# For parallel, do  'module load hdf5_par'
# To use the 64bit compilation model of xlC/xlf, you should set USE_64=TRUE and:
#  for serial, do   'module load hdf5_64'
#  for parallel, do 'module load hdf5_par_64'

ifneq ($(USE_HDF),FALSE)
  # The NERSC HDF5 modules use different variables in different HDF versions (not smart)
  #[NOTE: the HDF5C variable in some of the modules has options that don't compile Chombo]
  ifeq ($(HDF5_PAR_DIR),)
    HDFINCFLAGS=-I$(HDF5_DIR)/include 
    HDFMPIINCFLAGS=-I$(HDF5_DIR)/include 
  else
    HDFINCFLAGS=-I$(HDF5_PAR_DIR)/include 
    HDFMPIINCFLAGS=-I$(HDF5_PAR_DIR)/include 
  endif
  ifeq ($(HDF5_LIB),)
    HDFLIBFLAGS=$(HDF5) -lhdf5_hl -lhdf5
    HDFMPILIBFLAGS=$(HDF5) -lhdf5_hl -lhdf5
  else
    HDFLIBFLAGS=$(HDF5_LIB)
    HDFMPILIBFLAGS=$(HDF5_LIB)
  endif
endif

# Check that the right HDF module is loaded.
ifneq ($(USE_HDF),FALSE)
  ifeq ($(MPI),TRUE)
    ifeq ($(findstring parallel,$(HDF5_PAR_DIR)),)
      $(error HDF5 directory [$(HDF5_DIR)] is not parallel but MPI is TRUE.  Did you load the right module?)
    endif
  else
    ifeq ($(findstring serial,$(HDF5_DIR)),)
      $(error HDF5 directory [$(HDF5_DIR)] is not serial but MPI is FALSE.  Did you load the right module?)
    endif
  endif
endif

ifeq ($(USE_64),FALSE)
  $(error Are you sure you want to run non-64bit?)
endif
