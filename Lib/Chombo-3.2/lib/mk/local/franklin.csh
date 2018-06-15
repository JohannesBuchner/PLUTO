#!/bin/csh

module swap PrgEnv-pgi PrgEnv-gnu
module unload hdf5
module load hdf5_par
setenv USE_EB TRUE
setenv CVS_RSH ssh
alias make3 make all  DEBUG=FALSE OPT=HIGH DIM=3 MPI=TRUE
alias make2 make all  DEBUG=FALSE OPT=HIGH DIM=2 MPI=TRUE
