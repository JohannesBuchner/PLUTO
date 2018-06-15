# -*- Mode: Makefile -*- 

### This makefile produces a library for the name in the 'lbase' variable
### depending on the libraries named in the 'LibNames' variable

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

lbase := boxtools

LibNames :=

# put this first so 'Cfile' isn't the default target
include $(CHOMBO_HOME)/mk/Make.lib

.PHONY: Cfile

Cfile:
	$(RM) CFile.Cpp CFile.h
	grep -h -v 'include' SPACE.H REAL.H Stencil.H CONSTANTS.H Vector.H LoHiSide.H LoHiSide.cpp parstream.H MayDay.H MayDay.cpp_Ffree Misc.H  Arena.H Arena.cpp RefCountedPtr.H IntVect.H IntVect.cpp RealVect.H RealVect.cpp Interval.H Box.H Box.cpp ProblemDomain.H ProblemDomain.cpp  BoxIterator.H BoxIterator.cpp BaseFab.H BaseFabMacros.H BaseFabImplem.H  FArrayBox.H FArrayBox.cpp Pool.H Pool.cpp BitSet.H DenseIntVectSet.H Tuple.H TreeIntVectSet.H SPMD.H SPMDI.H BitSet.cpp DenseIntVectSet.cpp parstream.cpp TreeIntVectSet.cpp IntVectSet.H IntVectSet.cpp  SPMD.cpp DataIndex.H  BoxLayout.H LayoutIterator.H DataIterator.H  BaseFab.cpp BoxLayout.cpp DisjointBoxLayout.H DisjointBoxLayout.cpp LayoutIterator.cpp  DataIterator.cpp LayoutData.H LayoutDataI.H Copier.H BoxLayoutData.H BoxLayoutDataI.H BoxLayoutData.cpp Copier.cpp LevelData.H LevelDataI.H HDF5Portable.H CH_HDF5.H CH_HDF5.cpp MeshRefine.H MeshRefine.cpp memtrack.H memtrack.cpp BRMeshRefine.H BRMeshRefine.cpp LoadBalance.H LoadBalance.cpp  ../AMRTools/AMRIO.H ../AMRTools/AMRIO.cpp > CFile.Cpp
	grep -h -v 'include' SPACE.H REAL.H Stencil.H CONSTANTS.H Vector.H LoHiSide.H  parstream.H MayDay.H Misc.H  Arena.H RefCountedPtr.H IntVect.H RealVect.H  Interval.H Box.H ProblemDomain.H  BoxIterator.H  BaseFab.H BaseFabMacros.H BaseFabImplem.H  FArrayBox.H Pool.H  BitSet.H  DenseIntVectSet.H  Tuple.H TreeIntVectSet.H SPMD.H SPMDI.H  IntVectSet.H DataIndex.H  BoxLayout.H LayoutIterator.H DataIterator.H  DisjointBoxLayout.H LayoutData.H LayoutDataI.H BoxLayoutData.H  Copier.H  LevelData.H LevelDataI.H HDF5Portable.H CH_HDF5.H  memtrack.H MeshRefine.H BRMeshRefine.H LoadBalance.H  ../AMRTools/AMRIO.H > CFile.h
	cp -p CFile.Cpp ~/AMR-CCA/CFile.cpp
	cp -p CFile.h   ~/AMR-CCA/CFile.h
