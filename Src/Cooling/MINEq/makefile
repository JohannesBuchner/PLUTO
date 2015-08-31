# Makefile for MINEq cooling   

VPATH        += $(SRC)/Cooling/MINEq/
INCLUDE_DIRS += -I$(SRC)/Cooling/MINEq/

COOL_OBJ  = comp_equil.o ion_init.o jacobian.o maxrate.o radiat.o make_tables.o
OBJ     += $(COOL_OBJ)
HEADERS += cooling.h  

$(COOL_OBJ):  $(HEADERS) cooling_defs.h
 

