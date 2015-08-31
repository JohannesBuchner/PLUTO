# Makefile for SNEq sources   

VPATH        += $(SRC)/Cooling/H2_COOL
INCLUDE_DIRS += -I$(SRC)/Cooling/H2_COOL

COOL_OBJ = jacobian.o maxrate.o radiat.o comp_equil.o
OBJ     += $(COOL_OBJ)
HEADERS += cooling.h

$(COOL_OBJ):  $(HEADERS)
 

