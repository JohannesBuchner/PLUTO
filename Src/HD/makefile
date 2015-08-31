# Makefile for the HydroDynamics (HD) Module
#
# Note that we add Src/MHD to VPATH since HD and MHD
# share some files (see below)

VPATH += $(SRC)/HD $(SRC)/MHD

INCLUDE_DIRS += -I$(SRC)/HD
OBJ +=  ausm.o eigenv.o fluxes.o mappers.o  prim_eqn.o \
        hll_speed.o hll.o  hllc.o set_solver.o \
        tvdlf.o two_shock.o roe.o   

# The next set of files are shared between the HD and MHD directories

OBJ +=  get_area_flux.o rhs.o rhs_source.o

