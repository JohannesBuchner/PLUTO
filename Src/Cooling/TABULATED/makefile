# Makefile for Raymond sources   

VPATH        += $(SRC)/Cooling/TABULATED
INCLUDE_DIRS += -I$(SRC)/Cooling/TABULATED

COOL_OBJ = jacobian.o maxrate.o radiat.o
OBJ     += $(COOL_OBJ)
HEADERS += cooling.h

$(COOL_OBJ):  $(HEADERS) 
 

