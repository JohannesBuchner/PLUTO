# *********************************************************
#                                                     
#                  PLUTO 4.2  Makefile  
#                                                     
# *********************************************************

pluto:                              # Default target

ARCH         =                      # machine architecture flags here
PLUTO_DIR    =                      # define PLUTO main directory here
SRC          = $(PLUTO_DIR)/Src
INCLUDE_DIRS = -I. -I$(SRC)
VPATH        = ./:$(SRC)/New:$(SRC):$(SRC)/Time_Stepping:$(SRC)/States

include $(PLUTO_DIR)/Config/$(ARCH)

# ---------------------------------------------------------
#         Set headers and object files 
# ---------------------------------------------------------

HEADERS = pluto.h prototypes.h structs.h definitions.h macros.h mod_defs.h plm_coeffs.h
OBJ = adv_flux.o arrays.o boundary.o check_states.o  \
      cmd_line_opt.o entropy_switch.o  \
      flag_shock.o flatten.o get_nghost.o   \
      init.o int_bound_reset.o input_data.o \
      mappers3D.o mean_mol_weight.o \
      parse_file.o plm_coeffs.o rbox.o \
      set_indexes.o set_geometry.o set_output.o \
      tools.o var_names.o  

OBJ += bin_io.o colortable.o initialize.o jet_domain.o \
       main.o restart.o runtime_setup.o show_config.o  \
       set_image.o set_grid.o startup.o split_source.o \
       userdef_output.o write_data.o write_tab.o \
       write_img.o write_vtk.o 

include $(SRC)/Math_Tools/makefile

# ---------------------------------------------------------
#  Define macros by adding -D<name> where <name> has been
#  set to TRUE in the system configuration file (.defs) 
# ---------------------------------------------------------

ifeq ($(strip $(PARALLEL)), TRUE)
 CFLAGS += -I$(SRC)/Parallel -DPARALLEL
 include $(SRC)/Parallel/makefile
 ifeq ($(strip $(USE_ASYNC_IO)), TRUE)
  CFLAGS += -DUSE_ASYNC_IO
 endif
endif

ifeq ($(strip $(USE_HDF5)), TRUE)
 CFLAGS += -DUSE_HDF5
 OBJ    += hdf5_io.o
endif
      
ifeq ($($strip $(USE_PNG)), TRUE)
 CFLAGS += -DUSE_PNG
endif

-include local_make

# ---------------------------------------------------------
#   Additional_header_files_here   ! dont change this line
# ---------------------------------------------------------


# ---------------------------------------------------------
#   Additional_object_files_here   ! dont change this line
# ---------------------------------------------------------


# ---------------------------------------------------------
#    PLUTO target rule
# ---------------------------------------------------------

pluto: $(OBJ) 
	$(CC) $(OBJ) $(LDFLAGS) -o $@

# ---------------------------------------------------------
#                    Suffix rule
# ---------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE_DIRS) $<

clean:
	@rm -f	*.o
	@echo make clean: done

# ---------------------------------------------------------
#          Dependencies for object files
# ---------------------------------------------------------

$(OBJ):  $(HEADERS)

