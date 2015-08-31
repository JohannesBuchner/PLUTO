# Makefile for the Relativistic HydroDynamics (RHD) Module
#
# Note that we add Src/RMHD to VPATH since RHD and RMHD
# module share some files (see below)

VPATH += $(SRC)/RHD/:$(SRC)/RMHD
INCLUDE_DIRS += -I$(SRC)/RHD/
OBJ += eigenv.o energy_solve.o entropy_solve.o \
       fluxes.o mappers.o pressure_fix.o prim_eqn.o \
       hll_speed.o hll.o hllc.o \
       set_solver.o tvdlf.o two_shock.o
 
#  The following files are shared by RHD/ and RMHD/ modules

OBJ += four_vel.o rhs.o vel_limiter.o   
