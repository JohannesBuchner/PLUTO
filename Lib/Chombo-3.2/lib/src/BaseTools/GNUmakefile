# -*- Mode: Makefile -*- 

### This makefile produces a library for the name in the 'lbase' variable
### depending on the libraries named in the 'LibNames' variable

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

lbase := basetools
LibNames :=

# put this first so 'Cfile' isn't the default target
include $(CHOMBO_HOME)/mk/Make.lib

.PHONY: Cfile

Cfile:
	$(RM) CFile.Cpp CFile.h
	grep -h -v 'include' SPACE.H REAL.H  CONSTANTS.H Vector.H LoHiSide.H LoHiSide.cpp parstream.H MayDay.H MayDay.cpp_Ffree Misc.H  Arena.H Arena.cpp RefCountedPtr.H ../BoxTools/IntVect.H ../BoxTools/IntVect.cpp ../BoxTools/RealVect.H ../BoxTools/RealVect.cpp Interval.H ../BoxTools/Box.H ../BoxTools/Box.cpp ../BoxTools/ProblemDomain.H ../BoxTools/ProblemDomain.cpp ../BoxTools/BoxIterator.H ../BoxTools/BoxIterator.cpp ../BoxTools/BaseFab.H ../BoxTools/BaseFabMacros.H ../BoxTools/BaseFabImplem.H ../BoxTools/FArrayBox.H ../BoxTools/FArrayBox.cpp Pool.H Pool.cpp BitSet.H ../BoxTools/DenseIntVectSet.H Tuple.H ../BoxTools/TreeIntVectSet.H ../BoxTools/SPMD.H ../BoxTools/SPMDI.H BitSet.cpp ../BoxTools/DenseIntVectSet.cpp parstream.cpp ../BoxTools/TreeIntVectSet.cpp ../BoxTools/IntVectSet.H ../BoxTools/IntVectSet.cpp  ../BoxTools/SPMD.cpp ../BoxTools/DataIndex.H  ../BoxTools/BoxLayout.H ../BoxTools/LayoutIterator.H ../BoxTools/DataIterator.H ../BoxTools/BaseFab.cpp ../BoxTools/BoxLayout.cpp ../BoxTools/DisjointBoxLayout.H ../BoxTools/DisjointBoxLayout.cpp ../BoxTools/LayoutIterator.cpp ../BoxTools/DataIterator.cpp ../BoxTools/LayoutData.H ../BoxTools/LayoutDataI.H ../BoxTools/Copier.H ../BoxTools/BoxLayoutData.H ../BoxTools/BoxLayoutDataI.H ../BoxTools/BoxLayoutData.cpp ../BoxTools/Copier.cpp ../BoxTools/LevelData.H ../BoxTools/LevelDataI.H HDF5Portable.H ../BoxTools/CH_HDF5.H ../BoxTools/CH_HDF5.cpp ../BoxTools/MeshRefine.H ../BoxTools/MeshRefine.cpp memtrack.H memtrack.cpp ../BoxTools/BRMeshRefine.H ../BoxTools/BRMeshRefine.cpp LoadBalance.H LoadBalance.cpp ../AMRTools/AMRIO.H ../AMRTools/AMRIO.cpp > CFile.Cpp
	grep -h -v 'include' SPACE.H REAL.H  CONSTANTS.H Vector.H LoHiSide.H  parstream.H MayDay.H Misc.H  Arena.H RefCountedPtr.H ../BoxTools/IntVect.H ../BoxTools/RealVect.H  Interval.H ../BoxTools/Box.H ../BoxTools/ProblemDomain.H ../BoxTools/BoxIterator.H  ../BoxTools/BaseFab.H ../BoxTools/BaseFabMacros.H ../BoxTools/BaseFabImplem.H  ../BoxTools/FArrayBox.H Pool.H  BitSet.H  ../BoxTools/DenseIntVectSet.H  Tuple.H ../BoxTools/TreeIntVectSet.H ../BoxTools/SPMD.H ../BoxTools/SPMDI.H  ../BoxTools/IntVectSet.H ../BoxTools/DataIndex.H  ../BoxTools/BoxLayout.H LayoutIterator.H ../BoxTools/DataIterator.H  ../BoxTools/DisjointBoxLayout.H ../BoxTools/LayoutData.H ../BoxTools/LayoutDataI.H ../BoxTools/BoxLayoutData.H ../BoxTools/Copier.H  ../BoxTools/LevelData.H ../BoxTools/LevelDataI.H HDF5Portable.H ../BoxTools/CH_HDF5.H  memtrack.H ../BoxTools/MeshRefine.H BRMeshRefine.H LoadBalance.H  ../AMRTools/AMRIO.H > CFile.h
	cp -p CFile.Cpp ~/AMR-CCA/CFile.cpp
	cp -p CFile.h   ~/AMR-CCA/CFile.h
