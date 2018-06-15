include ../Automake.rules
nodist_fort_HEADERS = \
    CellToEdgeF_F.H\
    EdgeToCellF_F.H\
    DivergenceF_F.H\
    GradientF.1D_F.H GradientF.2D_F.H
## You need to make GradientF.{1,2}D.ChF -- symlinks to or copies of GradientF.ChF.
fortdir = $(pkgincludedir)
GENERATED_FORTRAN = \
    CellToEdgeF.f\
    EdgeToCellF.f\
    DivergenceF.f\
    GradientF.1D.f  GradientF.2D.f
EXTRA_DIST += *.ChF
example_MultiDimPhase_util_LTLIBRARIES = libexample_MultiDimPhase_util.la
nodist_libexample_MultiDimPhase_util_la_SOURCES = $(GENERATED_FORTRAN)\
    ReductionCopier.cpp\
    ReductionOps.cpp\
    Mask.cpp\
    Divergence.cpp\
    Gradient.cpp\
    EdgeToCell.cpp\
    CellToEdge.cpp
example_MultiDimPhase_utildir = $(pkglibdir)
EXTRA_DIST += *.cpp
headers_HEADERS = $(srcdir)/*.H
headersdir = $(pkgincludedir)
AM_LDFLAGS += -L$(CHOMBO_INSTALLDIR)/lib/Chombo \
 -lboxtools1D -lboxtools2D \
 -lamrtools1D -lamrtools2D \
 -lamrelliptic1D -lamrelliptic2D
