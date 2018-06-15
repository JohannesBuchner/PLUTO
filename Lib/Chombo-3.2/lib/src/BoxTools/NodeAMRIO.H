#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeAMRIO.H
// petermc, 4 March 2003
// adapted from AMRIO by DTGraves, Fri, Dec 3, 1999

#ifndef _NODEAMRIO_H_
#define _NODEAMRIO_H_
#include "REAL.H"

#include "LevelData.H"
#include "RealVect.H"
#include "Vector.H"
#include "NodeFArrayBox.H"
#include <string>

#include "IntVect.H"
#include "NamespaceHeader.H"
using std::string;
using std::istream;

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   This is blocking.
*/
void
WriteAMRHierarchyHDF5(/// file to send output to
                      const string& filename,
                      /// grids at each level
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      /// data indexed by level
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                      /// names of variables in <i>a_vectData</i>
                      const Vector<string>& a_vectNames,
                      /// domain at coarsest level
                      const Box& a_domain,
                      /// grid spacing at coarsest level
                      const Real& a_dx,
                      /// time step at coarsest level
                      const Real& a_dt,
                      /// time
                      const Real& a_time,
                      /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                      const Vector<int>& a_vectRatio,
                      /// number of levels to write out
                      const int& a_numLevels,
                      /// physical coordinates of lower corner of cell at IntVect 0
                      const RealVect& a_origin = RealVect::Zero,
                      /// Interval of components to be written
                      const Interval& a_comps = Interval() );

/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   This is blocking.
*/
void
WriteAnisotropicAMRHierarchyHDF5(/// file to send output to
                      const string& filename,
                      /// grids at each level
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      /// data indexed by level
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                      /// names of variables in <i>a_vectData</i>
                      const Vector<string>& a_vectNames,
                      /// domain at coarsest level
                      const Box& a_domain,
                      /// grid spacing in each direction at coarsest level
                      const RealVect& a_dx,
                      /// time step at coarsest level
                      const Real& a_dt,
                      /// time
                      const Real& a_time,
                      /// refinement ratios in each direction between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                      const Vector<IntVect>& a_vectRatios,
                      /// number of levels to write out
                      const int& a_numLevels,
                      /// physical coordinates of lower corner of cell at IntVect 0
                      const RealVect& a_origin = RealVect::Zero,
                      /// Interval of components to be written
                      const Interval& a_comps = Interval() );
///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   This function differs from the previous one in that
   <i>a_vectData</i>[<i>ind</i>][<i>lev</i>]
   indexes a LevelData<NodeFArrayBox>* (which may have several components)
   at level <i>lev</i>.

   This is blocking.
*/
void
WriteAMRHierarchyHDF5(/// file to send output to
                      const string& filename,
                      /// grids at each level
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      /// data indexed by set and then by level
                      const Vector<Vector<LevelData<NodeFArrayBox>* > >& a_vectData,
                      /// names of variables in <i>a_vectData</i>:  first the names of the components of <i>a_vectData</i>[0], then the names of the components of <i>a_vectData[1]</i>, etc.
                      const Vector<string>& a_vectNames,
                      /// domain at coarsest level
                      const Box& a_domain,
                      /// grid spacing at coarsest level
                      const Real& a_dx,
                      /// time step at coarsest level
                      const Real& a_dt,
                      /// time
                      const Real& a_time,
                      /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                      const Vector<int>& a_vectRatio,
                      /// number of levels to write out
                      const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   This is not blocking.
*/
void
WriteAMRHierarchyHDF5(/// handle to send output to
                      HDF5Handle& handle,
                      /// grids at each level
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      /// data indexed by level
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                      /// names of variables in <i>a_vectData</i>
                      const Vector<string>& a_vectNames,
                      /// domain at coarsest level
                      const Box& a_domain,
                      /// grid spacing at coarsest level
                      const Real& a_dx,
                      /// time step at coarsest level
                      const Real& a_dt,
                      /// time
                      const Real& a_time,
                      /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                      const Vector<int>& a_vectRatio,
                      /// number of levels to write out
                      const int& a_numLevels,
                      /// physical coordinates of lower corner of cell at IntVect 0
                      const RealVect& a_origin = RealVect::Zero,
                      /// Interval of components to be written
                      const Interval& a_comps = Interval() );

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   This is not blocking.
*/
void
WriteAnisotropicAMRHierarchyHDF5(/// handle to send output to
                      HDF5Handle& handle,
                      /// grids at each level
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      /// data indexed by level
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                      /// names of variables in <i>a_vectData</i>
                      const Vector<string>& a_vectNames,
                      /// domain at coarsest level
                      const Box& a_domain,
                      /// grid spacing at coarsest level
                      const RealVect& a_dx,
                      /// time step at coarsest level
                      const Real& a_dt,
                      /// time
                      const Real& a_time,
                      /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                      const Vector<IntVect>& a_vectRatios,
                      /// number of levels to write out
                      const int& a_numLevels,
                      /// physical coordinates of lower corner of cell at IntVect 0
                      const RealVect& a_origin = RealVect::Zero,
                      /// Interval of components to be written
                      const Interval& a_comps = Interval() );

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   This is blocking.
*/
void
WriteAMRHierarchyHDF5(/// file to send output to
                      const string& filename,
                      /// grids at each level
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      /// data indexed by level
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                      /// domain at coarsest level
                      const Box& a_domain,
                      /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                      const Vector<int>& a_vectRatio,
                      /// number of levels to write out
                      const int& a_numLevels);

///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   This is not blocking.
*/
void
WriteAMRHierarchyHDF5(/// handle to send output to
                      HDF5Handle& handle,
                      /// grids at each level
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      /// data indexed by level
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                      /// domain at coarsest level
                      const Box& a_domain,
                      /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                      const Vector<int>& a_vectRatio,
                      /// number of levels to write out
                      const int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   Returns status code with values:
   - 0: success
   - -1: bogus number of levels
   - -2: bogus number of components
   - -3: error in readlevel
   - -4: file open failed

   This is blocking.
*/
int
ReadAMRHierarchyHDF5(/// file to read input from
                     const string& filename,
                     /// grids at each level
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     /// data indexed by level
                     Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                     /// names of variables in <i>a_vectData</i>
                     Vector<string>& a_vectNames,
                     /// domain at coarsest level
                     Box& a_domain,
                     /// grid spacing at coarsest level
                     Real& a_dx,
                     /// time step at coarsest level
                     Real& a_dt,
                     /// time
                     Real& a_time,
                     /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                     Vector<int>& a_refRatio,
                     /// number of levels to read in
                     int& a_numLevels,
                     bool a_setGhost = false);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   Returns status code with values:
   - 0: success
   - -1: bogus number of levels
   - -2: bogus number of components
   - -3: error in readlevel
   - -4: file open failed

   This is not blocking.
*/
int
ReadAMRHierarchyHDF5(/// handle to read input from
                     HDF5Handle& handle,
                     /// grids at each level
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     /// data indexed by level
                     Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                     /// names of variables in <i>a_vectData</i>
                     Vector<string>& a_vectNames,
                     /// domain at coarsest level
                     Box& a_domain,
                     /// grid spacing at coarsest level
                     Real& a_dx,
                     /// time step at coarsest level
                     Real& a_dt,
                     /// time
                     Real& a_time,
                     /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                     Vector<int>& a_refRatio,
                     /// number of levels to read in
                     int& a_numLevels,
                     bool a_setGhost = false);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   Returns status code with values:
   - 0: success
   - -1: bogus number of levels
   - -2: bogus number of components
   - -3: error in readlevel
   - -4: file open failed

   This is blocking.
*/
int
ReadAMRHierarchyHDF5(/// file to read input from
                     const string& filename,
                     /// grids at each level
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     /// data indexed by level
                     Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                     /// domain at coarsest level
                     Box& a_domain,
                     /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                     Vector<int>& a_refRatio,
                     /// number of levels to read in
                     int& a_numLevels,
                     bool a_setGhost = false);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   Returns status code with values:
   - 0: success
   - -1: bogus number of levels
   - -2: bogus number of components
   - -3: error in readlevel
   - -4: file open failed

   This is not blocking.
*/
int
ReadAMRHierarchyHDF5(/// handle to input from
                     HDF5Handle& handle,
                     /// grids at each level
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     /// data indexed by level
                     Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                     /// domain at coarsest level
                     Box& a_domain,
                     /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                     Vector<int>& a_refRatio,
                     /// number of levels to read in
                     int& a_numLevels,
                     bool a_setGhost = false);

///
/**
   Writes a plotfile using the same format as WriteAMRHierarchyHDF5, but
   for a NodeFArrayBox.  This is useful for debugging.
   *<i>a_dataPtr</i> is written to a file named nfab.hdf5
*/
void
writeNFAB(const NodeFArrayBox* a_dataPtr);

///
/**
   View *<i>a_dataPtr</i> by writing it to an HDF5 plotfile (to a
   temporary file) and then running ChomboVis with a python script
   which brings up a data browser by default. The file has the same
   format as WriteAMRHierarchyHDF5, but for a single NodeFArrayBox.
   This is useful for debugging.
*/
void
viewNFAB(const NodeFArrayBox* a_dataPtr);

///
/**
   Writes a plotfile using the same format as WriteAMRHierarchyHDF5, but
   for a NodeFArrayBox.  This is useful for debugging.
   *<i>a_dataPtr</i> is written to the file given by a_filename.
*/
void
writeNFABname(const NodeFArrayBox* a_dataPtr,
              const char*      a_filename);

///
/**
   Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
   for a single LevelData<NodeFArrayBox>.  Useful for debugging.
   *<i>a_dataPtr</i> is written to a file named LDF.hdf5.
*/
void
writeNodeLevel(const LevelData<NodeFArrayBox>* a_dataPtr);

///
/**
   View *<i>a_dataPtr</i> by writing it to an HDF5 plotfile (to a
   temporary file) and then running ChomboVis with a python script
   which brings up a data browser by default. The file has the same
   format as WriteAMRHierarchyHDF5, but for a single
   LevelData<NodeFArrayBox>.  This is useful for debugging.
*/
void
viewNodeLevel(const LevelData<NodeFArrayBox>* a_dataPtr);

///
/**
   Write a plotfile using the same format as WriteAMRHierarchyHDF5, but
   for a single LevelData<FArrayBox>. Useful for debugging.
   *<i>a_dataPtr</i> is written to the file given by <i>a_filename</i>.
*/
void
writeNodeLevelname(const LevelData<NodeFArrayBox>* a_dataPtr,
                   const char*                 a_filename);

///
/**
   Write a plotfile with levels that don't necessarily begin with
   0 (will translate things so that base level comes out as level 0, etc.).
   Levels in Vector arguments should correspond with actual levels,
   and so should be defined for <i>a_levels</i>.begin() through
   <i>a_levels</i>.end().

   This is blocking.
*/
void
WritePartialAMRHierarchyHDF5(/// file to send output to
                             const string& filename,
                             /// grids at each level
                             const Vector<DisjointBoxLayout>& a_vectGrids,
                             /// data indexed by level
                             const Vector<LevelData<FArrayBox>* > & a_vectData,
                             /// names of variables in <i>a_vectData</i>
                             const Vector<string>& a_vectNames,
                             /// domain at base level given by <i>a_levels</i>.begin()
                             const Box& a_baseDomain,
                             /// grid spacing at base level
                             const Real& a_baseDx,
                             /// time step at base level
                             const Real& a_dt,
                             /// time
                             const Real& a_time,
                             /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                             const Vector<int>& a_vectRatio,
                             /// indices of levels to write out
                             const Interval& a_levels);

///
/**
   Write a plotfile with levels that don't necessarily begin with
   0 (will translate things so that base level comes out as level 0, etc.).
   Levels in Vector arguments should correspond with actual levels,
   and so should be defined for <i>a_levels</i>.begin() through
   <i>a_levels</i>.end().

   This function differs from the previous one in that
   <i>a_vectData</i>[<i>ind</i>][<i>lev</i>]
   indexes a LevelData<NodeFArrayBox>* (which may have several components)
   at level <i>lev</i>.

   This is blocking.
*/
void
WritePartialAMRHierarchyHDF5(/// file to send output to
                             const string& filename,
                             /// grids at each level
                             const Vector<DisjointBoxLayout>& a_vectGrids,
                             /// data indexed by set and then by level
                             const Vector<Vector<LevelData<NodeFArrayBox>* > >& a_vectData,
                             /// names of variables in <i>a_vectData</i>
                             const Vector<string>& a_vectNames,
                             /// domain at base level given by <i>a_levels</i>.begin()
                             const Box& a_baseDomain,
                             /// grid spacing at base level
                             const Real& a_baseDx,
                             /// time step at base level
                             const Real& a_dt,
                             /// time
                             const Real& a_time,
                             /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
                             const Vector<int>& a_vectRatio,
                             /// indices of levels to write out
                             const Interval& a_levels);

#include "NamespaceFooter.H"
#endif // include guard
