#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Fri, Dec 3, 1999

#ifndef _AMRIO_H_
#define _AMRIO_H_

#include "REAL.H"

#ifndef WRAPPER
#include "LevelData.H"
#include "REAL.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "IntVectSet.H"
#include "IVSFAB.H"
#include "RealVect.H"
#include <string>
#include "NamespaceHeader.H"
using std::string;
using std::istream;
#endif

#ifdef CH_USE_HDF5
///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.

*/
void
WriteAMRHierarchyHDF5(const string& filename,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_vectRatio,
                      const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing in each direction at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio in each direction at all levels
 (ith entry is refinement ratio in each direction between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.

*/
void
WriteAnisotropicAMRHierarchyHDF5(
    const string& filename,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_vectRatios,
    const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
handle  :  handle to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is not blocking.

*/
void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_vectRatio,
                      const int& a_numLevels);

//
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
handle  :  handle to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing in each direction at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio in each direction at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is not blocking.

*/
void
WriteAnisotropicAMRHierarchyHDF5(
    HDF5Handle& handle,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_vectRatios,
    const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_domain :  domain at coarsest level.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.
*/
void
WriteAMRHierarchyHDF5(const string& filename,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Box& a_domain,
                      const Vector<int>& a_vectRatio,
                      const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 handle  :  handle to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_domain :  domain at coarsest level.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is not blocking.
*/
void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Box& a_domain,
                      const Vector<int>& a_vectRatio,
                      const int& a_numLevels);

#ifndef WRAPPER
///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking
*/
int
ReadAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<int>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking
*/
int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<int>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking.
*/
int
ReadAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<int>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking.
*/
int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<int>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking
*/
int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     RealVect& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking
*/
int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     RealVect& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking.
*/
int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking.
*/
int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels);

#endif

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file created uses the Boxes that underlie the
    IntVectSet.  If a_domain is 0, then the domain used in the file is the
    minBox of a_dataPtr.
*/
void viewIVS(const IntVectSet* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default.  If a_domain is 0, then the domain used in the file is the
    minBox of a_dataPtr.
*/
void viewVectorBox(const Vector<Box>* a_dataPtr);

///
/** Writes a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a FArrayBox.  This is useful for debugging.  *a_dataPtr is written
    to a file named fab.hdf5
*/
void writeFAB(const FArrayBox* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single FArrayBox.  This is useful for debugging.
*/
void viewFAB(const FArrayBox* a_dataPtr);

///
/** Make a copy of *a_dataPtr as FArrayBox, then call viewFAB on FArrayBox.
*/
void viewBFI(const BaseFab<int>* a_dataPtr);

///
/** Make a copy of *a_dataPtr as FArrayBox, then call viewFAB on FArrayBox.
*/
void viewBFIV(const BaseFab<IntVect>* a_dataPtr);

///
/** Make a copy of *a_dataPtr as FArrayBox, then call viewFAB on FArrayBox.
*/
void viewBFRV(const BaseFab<RealVect>* a_dataPtr);

///
/** Make a copy of *a_dataPtr as FArrayBox, then call viewFAB on FArrayBox.
*/
void viewIVSFAB(const IVSFAB<Real>* a_dataPtr);

///
/** Make a copy of *a_dataPtr as FArrayBox, then call viewFAB on FArrayBox.
*/
void viewIVSFABI(const IVSFAB<int>* a_dataPtr);

///
/** Make a copy of *a_dataPtr as FArrayBox, then call viewFAB on FArrayBox.
*/
void viewIVSFABIV(const IVSFAB<IntVect>* a_dataPtr);

///
/** Make a copy of *a_dataPtr as FArrayBox, then call viewFAB on FArrayBox.
*/
void viewIVSFABRV(const IVSFAB<RealVect>* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running chombobrowser with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single FArrayBox.  This is useful for debugging.
*/
void browseFAB(const FArrayBox* a_dataPtr);

///
/** Writes a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a BaseFab<Real>.  This is useful for debugging.  *a_dataPtr is written
    to a file named fab.hdf5
*/
void
writeBFR(const BaseFab<Real>* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single BaseFab<Real>.  This is useful for debugging.
*/
void
viewBFR(const BaseFab<Real>* a_dataPtr);

///
/** Writes a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a FArrayBox.  This is useful for debugging.  *a_dataPtr is written
    to the file given by a_filename.  If a non-empty Vector of strings is
    given, they are used as the component names in the HDF5 header; else
    the names are set to "component_0", "component_1", ....
*/
void
writeFABname(const FArrayBox      * a_dataPtr,
             const char           * a_filename,
             const Vector<string> & a_compNames = Vector<string>());

///
/** View a modified version of *a_dataPtr where the data are set to zero
    on the grids of *a_dataFinePtr coarsened by a_refRatio.
    Write the results to an HDF5 plotfile (to a temporary file)
    and then run VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single LevelData<FArrayBox>.  This is useful for debugging.
*/
void
viewLevelNoFine(const LevelData<FArrayBox>* a_dataPtr,
                const LevelData<FArrayBox>* a_dataFinePtr,
                int a_refRatio);

///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a single LevelData<FArrayBox>.  Useful for debugging.  *a_dataPtr is
    written to a file named LDF.hdf5.
*/
void
writeLevel(const LevelData<FArrayBox>* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single LevelData<FArrayBox>.  This is useful for debugging.
*/
void
viewLevel(const LevelData<FArrayBox>* a_dataPtr);

///
/** View fluxBoxes stored by *a_dataPtr by writing them to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single LevelData<FArrayBox>.  This is useful for debugging.
    NOTE: the components of \a a_dataPtr are stored in dimension-major order.
*/
void
viewFluxLevel(const LevelData<FluxBox>* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single LevelData<FArrayBox>.  This is useful for debugging.
*/
void
browseLevel(const LevelData<FArrayBox>* a_dataPtr);

///
/** View *a_dataPtr by writing only the interior values (no ghost cell
    values) to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single LevelData<FArrayBox>.  This is useful for debugging.
*/
void
viewLevelNoGhost(const LevelData<FArrayBox>* a_dataPtr);

///
/** View the multilevel dataset *a_dataPtr by writing a to an (temporary)
    HDF5 plotfile using the same format as WriteAMRHierarchyHDF5,
    for a Vector<LevelData<FArrayBox>*>, and then running VisIt from a
    python script to look at the data. Useful for debugging.
*/
void
viewVectorLevel(const Vector<LevelData<FArrayBox>*>* a_dataPtr,
                const Vector<int>*                   a_refRatios);

///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a single LevelData<FArrayBox>. Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeLevelname(const LevelData<FArrayBox>* a_dataPtr,
               const char*                 a_filename);

///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a Vector<LevelData<FArrayBox>*>. Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeVectorLevelName(const Vector<LevelData<FArrayBox>*>* a_dataPtr,
                     const Vector<int>*                   a_refRatios,
                     const char*                          a_filename);

///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a single DisjointBoxLayout (do this by creating a LevelData<FArrayBox>
    and then calling writeLevel.  Useful for debugging.  *a_dataPtr is
    written to a file named DBL.hdf5.
*/
void
writeDBL(const DisjointBoxLayout* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single DisjointBoxLayout (do this by creating a
    LevelData<FArrayBox> and then calling writeLevel.  This is useful
    for debugging.
*/
void
viewDBL(const DisjointBoxLayout* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt with a python script which brings up a data
    browser by default. The file has the same format as WriteAMRHierarchyHDF5,
    but for a single LevelData<BaseFab<int> >.  This is useful for debugging.
*/
void
viewLevelBFI(const LevelData<BaseFab<int> >* a_dataPtr);

void writeCopier(const Copier* a_copier);

void viewCopier(const Copier* a_copier);

///
/** Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
    for a single DisjointBoxLayout (do this by creating a LevelData<FArrayBox>
    and then calling writeLevel.  Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeDBLname(const DisjointBoxLayout* a_dataPtr,
             const char*              a_filename);

///
/** Write a plotfile with levels that don't necessarily begin with
0 (will translate things so that lBase comes out as level 0, etc)
Levels in Vector arguments should correspond with actual levels.
(and so should be defined for a_levels.begin() through a_levels.end().
   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at base level (given by  a_levels.begin()).\\
 a_dx     :  grid spacing at base level.\\
 a_dt     :  time step at base level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.

*/
void
WritePartialAMRHierarchyHDF5(const string& filename,
                             const Vector<DisjointBoxLayout>& a_vectGrids,
                             const Vector<LevelData<FArrayBox>* > & a_vectData,
                             const Vector<string>& a_vectNames,
                             const Box& a_baseDomain,
                             const Real& a_baseDx,
                             const Real& a_dt,
                             const Real& a_time,
                             const Vector<int>& a_vectRatio,
                             const Interval& a_levels);

#endif // CH_USE_HDF5

#include "NamespaceFooter.H"

#endif
