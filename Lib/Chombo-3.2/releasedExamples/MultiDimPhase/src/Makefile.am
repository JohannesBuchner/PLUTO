include ../Automake.rules
nodist_fort_HEADERS = \
    phaseF_F.H\
    SolidF_F.H\
    PROB_F_F.H
fortdir = $(pkgincludedir)
GENERATED_FORTRAN = \
    phaseF.f\
    SolidF.f\
    PROB_F.f
EXTRA_DIST += *.ChF
example_MultiDimPhase_src_LTLIBRARIES = libexample_MultiDimPhase_src.la
nodist_libexample_MultiDimPhase_src_la_SOURCES = $(GENERATED_FORTRAN)\
    AdvectIBC.cpp\
    1DPotentialSolve.cpp\
    amrPhase.cpp
example_MultiDimPhase_srcdir = $(pkglibdir)
EXTRA_DIST += *.cpp
headers_HEADERS = $(srcdir)/*.H
headersdir = $(pkgincludedir)
AM_LDFLAGS += -L$(CHOMBO_INSTALLDIR)/lib/Chombo \
 -lboxtools1D -lboxtools2D \
 -lamrtools1D -lamrtools2D \
 -lamrelliptic1D -lamrelliptic2D
