include ../Automake.rules
nodist_fort_HEADERS = \
    OldPatchGodunovF_F.H\
    GodunovTraceF_F.H\
    IAMRGodunovF_F.H\
    PatchAdvectionF_F.H
fortdir = $(pkgincludedir)
GENERATED_FORTRAN = \
    OldPatchGodunovF.f\
    GodunovTraceF.f\
    IAMRGodunovF.f\
    PatchAdvectionF.f
EXTRA_DIST += *.ChF
example_MultiDimPhase_advectSrc_LTLIBRARIES = libexample_MultiDimPhase_advectSrc.la
nodist_libexample_MultiDimPhase_advectSrc_la_SOURCES = $(GENERATED_FORTRAN)\
    OldLevelGodunov.cpp\
    OldLoHiCenter.cpp\
    OldPatchGodunov.cpp\
    OldPhysIBC.cpp\
    PatchAdvection.cpp\
    GodunovTrace.cpp
example_MultiDimPhase_advectSrcdir = $(pkglibdir)
EXTRA_DIST += *.cpp
headers_HEADERS = $(srcdir)/*.H
headersdir = $(pkgincludedir)
headers_HEADERS += $(srcdir)/*.fh
AM_LDFLAGS += -L$(CHOMBO_INSTALLDIR)/lib/Chombo \
 -lboxtools1D -lboxtools2D \
 -lamrtools1D -lamrtools2D \
 -lamrelliptic1D -lamrelliptic2D
