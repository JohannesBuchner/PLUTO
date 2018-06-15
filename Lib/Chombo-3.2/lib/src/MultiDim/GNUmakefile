# -*- Mode: Makefile -*- 

### This makefile produces a library for the name in the 'lbase' variable
### depending on the libraries named in the 'LibNames' variable

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

lbase := multidim

LibNames := 

include $(CHOMBO_HOME)/mk/Make.lib

newconfig := .$(system)$(_obj64).$(_cxxname).$(_fname)$(_precision)$(_debug)$(_opt)$(_profile)$(_mpi)$(_pic)$(XTRACONFIG)

newcppflags = $(CPPFLAGS) \
$(subst FALSE,,$(subst TRUE,-DCH_USE1D,$(USE_1D)))\
$(subst FALSE,,$(subst TRUE,-DCH_USE2D,$(USE_2D)))\
$(subst FALSE,,$(subst TRUE,-DCH_USE3D,$(USE_3D)))\
$(subst FALSE,,$(subst TRUE,-DCH_USE4D,$(USE_4D)))\
$(subst FALSE,,$(subst TRUE,-DCH_USE5D,$(USE_5D)))\
$(subst FALSE,,$(subst TRUE,-DCH_USE6D,$(USE_6D)))\
-I$(CHOMBO_HOME)/include/multidim


mylib : 
	$(MAKE) lib config="$(newconfig)" MULTIDIM=$(MULTIDIM) NAMESPACE=$(NAMESPACE) CPPFLAGS="$(newcppflags) "
