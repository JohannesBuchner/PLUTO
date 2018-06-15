#!/bin/csh

module swap PrgEnv-pgi PrgEnv-gnu
module unload hdf5
module load hdf5-parallel

alias make3 make all DIM=3 DEBUG=FALSE OPT=HIGH DIM=3 MPI=TRUE
alias make2 make all DIM=3 DEBUG=FALSE OPT=HIGH DIM=2 MPI=TRUE
