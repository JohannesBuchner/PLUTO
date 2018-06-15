#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Thurs, Aug 16, 2001

#ifndef _EBAMRIO_H_
#define _EBAMRIO_H_

#include "LevelData.H"
#include "REAL.H"
#include "Vector.H"
#include "EBCellFAB.H"
#include "IntVectSet.H"
#include "EBISBox.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "EBFluxFAB.H"
#include <string>
#include "NamespaceHeader.H"
class EBIndexSpace;

using std::string;
using std::istream;

#ifdef CH_USE_HDF5
///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Writes
   single-valued data using regular Chombo WriteAMRHierarchyHDF5

   {\bf Arguments:}\\
   a_filename :  file to output to.\\
   a_vectGrids : grids at each level.\\
   a_vectData :  data at each level.\\
   a_vectNames:  names of variables.\\
   a_domain :  domain at coarsest level.\\
   a_dx     :  grid spacing at coarsest level.\\
   a_dt     :  time step at coarsest level.\\
   a_time     :  time.\\
   a_vectRatio :  refinement ratio at all level.\\
   a_numLevels :  number of levels to output.\\
   a_replaceCovered:  flag for replacing covered values with a_coveredVal.\\
   a_coveredVal:  sets Covered values to these if a_replaceCovered is true.\\
   (ith entry is refinement ratio between levels i and i + 1).\\
   This is blocking.
*/
void
writeEBHDF5(const string& filename,
            const Vector<DisjointBoxLayout>& a_vectGrids,
            const Vector<LevelData<EBCellFAB>* > & a_vectData,
            const Vector<string>& a_vectNames,
            const ProblemDomain& a_domain,
            const Real& a_dx,
            const Real& a_dt,
            const Real& a_time,
            const Vector<int>& a_vectRatio,
            const int& a_numLevels,
            const bool& a_replaceCovered,
            const Vector<Real>& a_coveredValues,
            IntVect a_ghostVect= IntVect::Zero);

/**
    Two phases (each a Vector<LevelData<EBCellFAB>>) : most general case.
*/
void
writeEBHDF5(const string& filename,
            const Vector<DisjointBoxLayout>& a_vectGrids,
            const Vector<LevelData<EBCellFAB>* > * a_phase1,
            const Vector<LevelData<EBCellFAB>* > * a_phase2,
            const Vector<LevelData<FArrayBox>* > * a_levelSet,
            const Vector<string>& a_vectNames,
            const ProblemDomain& a_domain,
            const Real& a_dx,
            const Real& a_dt,
            const Real& a_time,
            const Vector<int>& a_vectRatio,
            const int& a_numLevels,
            const bool& a_replaceCovered,
            const Vector<Real>& a_coveredValues,
            IntVect a_ghostVect= IntVect::Zero);

///
/**
    Write a single LevelData<EBCellFAB>.
 */
void
writeEBHDF5(const string& a_filename,
            const DisjointBoxLayout& a_grids,
            const LevelData<EBCellFAB>& a_data,
            const Vector<string>& a_vectNames,
            const ProblemDomain& a_domain,
            const Real& a_dx,
            const Real& a_dt,
            const Real& a_time,
            const bool& a_replaceCovered,
            const Vector<Real>& a_coveredValues,
            IntVect a_ghostVect= IntVect::Zero);
///
/**
    Write a single EBCellFAB.
 */
void
writeEBHDF5(const string& a_filename,
            const Box& a_grid,
            const EBCellFAB& a_data,
            const Vector<string>& a_vectNames,
            const ProblemDomain& a_domain,
            const Real& a_dx,
            const Real& a_dt,
            const Real& a_time,
            const bool& a_replaceCovered,
            const Vector<Real>& a_coveredValues,
            IntVect a_ghostVect= IntVect::Zero);
///
/** Writes a plotfile using the same format as writeEBHDF5, but
    for a EBCellFAB.  This is useful for debugging.  *a_dataPtr is written
    to a file named fab.hdf5
*/
void
writeEBFAB(const EBCellFAB* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running ChomboVis with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single EBCellFAB.  This is useful for debugging.
*/
void
viewEBFAB(const EBCellFAB* a_dataPtr);


///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running browse with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single EBCellFAB.  This is useful for debugging.
*/
void
browseEBFAB(const EBCellFAB* a_dataPtr);

///
/** Writes a plotfile using the same format as writeEBHDF5, but
    for a EBCellFAB.  This is useful for debugging.  *a_dataPtr is written
    to the file given by a_filename.
*/
void
writeEBFABname(const EBCellFAB* a_dataPtr,
               const char*      a_filename);

///
/** View the face centered, regular data (as cell center data) in *a_dataPtr
    by writing data on the a_dir face to an HDF5 plotfile (a temporary file).
    Then, run VisIt with a python script, bringing up a data browser by
    default. The file has the same format as writeEBHDF5, but the components
    of the EBFluxFAB are written out successively.
    This is useful for debugging.
*/
void
viewEBFluxLD(const LevelData<EBFluxFAB>* a_dataPtr, int a_dir);

///
/** View the faced centered, regular data (as cell center data) in *a_dataPtr
    by writing it to an HDF5 plotfile (to a temporary file) and then running
    ChomboVis with a python script which brings up a data browser by default.
    The file has the same format as writeEBHDF5, but for a single EBFaceFAB.
    This is useful for debugging.
*/
void writeEBFluxLDname(const LevelData<EBFluxFAB>* a_fluxLD, int a_dir, 
                       const char* a_filename);

///
/** Write the faced centered, regular data (as cell center data) in *a_dataPtr
    by writing it to a named HDF5 plotfile 
*/
void
viewEBFace(const EBFaceFAB* a_dataPtr);


///
/** View the faced centered, regular data (as cell center data) in *a_dataPtr
    by writing it to an HDF5 plotfile (to a temporary file) and then running
    chombobrowser with a python script which brings up a data browser by default.
    The file has the same format as writeEBHDF5, but for a single EBFaceFAB.
    This is useful for debugging.
*/
void
browseEBFace(const EBFaceFAB* a_dataPtr);

///
/** Write a plotfile using the same format as writeEBHDF5, but
    for a single LevelData<EBCellFAB>.  Useful for debugging.  *a_dataPtr is
    written to a file named EBLDF.hdf5.
*/
void
writeEBLevel(const LevelData<EBCellFAB>* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running ChomboVis with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single LevelData<EBCellFAB>.  This is useful for debugging.
*/
void
viewEBLevel(const LevelData<EBCellFAB>* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running chombobrowser with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single LevelData<EBCellFAB>.  This is useful for debugging.
*/
void
browseEBLevel(const LevelData<EBCellFAB>* a_dataPtr);

///
/** Write a plotfile using the same format as writeEBHDF5, but
    for a single LevelData<EBCellFAB>. Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeEBLevelname(const LevelData<EBCellFAB>* a_dataPtr,
                 const char*                 a_filename);

///
/** Write a plotfile using writeEBHDF5.  Useful for debugging.  *a_dataPtr is
    written to a file named EBAMR.hdf5.
*/
void
writeEBAMR(const Vector<LevelData<EBCellFAB>* >* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running ChomboVis with a python script which brings up a data
    browser by default.  This is useful for debugging.
*/
void
viewEBAMR(const Vector<LevelData<EBCellFAB>* >* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running VisIt.  This is useful for debugging.
*/
void
visEBAMR(const Vector<LevelData<EBCellFAB>* >* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running chombobrowser with a python script which brings up a data
    browser by default.  This is useful for debugging.
*/
void
browseEBAMR(const Vector<LevelData<EBCellFAB>* >* a_dataPtr);

///
/** Write a plotfile using writeEBHDF5. Useful for debugging.  *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeEBAMRname(const Vector<LevelData<EBCellFAB>* >* a_dataPtr,
               const char*                           a_filename);


///
/** Writes a plotfile using the basic Chombovis format, but
    for a BaseIVFAB<Real>.  This is somewhat limited in functionality,
    since all it does is copy the BaseIVFAB to a FArrayBox in the appropriate
    cells, and then call writeFAB.  This is useful for debugging.
    *a_dataPtr is written to a file named ivfab.hdf5

*/
void
writeIVFAB(const BaseIVFAB<Real>* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    using writeIVFAB and then running ChomboVis with a python script which
    brings up a data browser by default. The file has the same format
    as writeEBHDF5, but for a single EBCellFAB.  This is useful for
    debugging.
*/
void
viewIVFAB(const BaseIVFAB<Real>* a_dataPtr);


///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running browse with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single BaseIVFAB<Real>.  This is useful for debugging.
*/
void
browseIVFAB(const BaseIVFAB<Real>* a_dataPtr);

///
/** Writes a plotfile using the same format as writeEBHDF5, but
    for a BaseIVFAB<Real>.  This is useful for debugging.
    *a_dataPtr is written to the file given by a_filename.
*/
void
writeIVFABname(const BaseIVFAB<Real>* a_dataPtr,
               const char*      a_filename);


///
/** Write a plotfile using the same format as writeLevel, but
    for a single LevelData<BaseIVFAB<Real> >.  Useful for debugging.
    *a_dataPtr is written to a file named IVLDF.hdf5.
*/
void
writeIVLevel(const LevelData<BaseIVFAB<Real> >* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running ChomboVis with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single LevelData<BaseIVFAB<Real> >.  This is useful for
    debugging.
*/
void
viewIVLevel(const LevelData<BaseIVFAB<Real> >* a_dataPtr);

///
/** View *a_dataPtr by writing it to an HDF5 plotfile (to a temporary file)
    and then running chombobrowser with a python script which brings up a data
    browser by default. The file has the same format as writeEBHDF5,
    but for a single LevelData<BaseIVFAB<Real> >.  This is useful for
    debugging.
*/
void
browseIVLevel(const LevelData<BaseIVFAB<Real> >* a_dataPtr);

///
/** Write a plotfile using the same format as writeEBHDF5, but
    for a single LevelData<BaseIVFAB<Real> >. Useful for debugging.
    *a_dataPtr is
    written to the file given by a_filename.
*/
void
writeIVLevelname(const LevelData<BaseIVFAB<Real> >* a_dataPtr,
                    const char*                 a_filename);

void
createEBFile(HDF5Handle& a_handle,
             const std::string& a_filename,
             int a_numLevels,
             const Vector<int>& a_refRatios,
             const ProblemDomain& a_coarseDomain,
             const RealVect& a_coarseDx,
             const IntVect& a_ghost);

void
createEBFile(HDF5Handle& a_handle,
             const std::string& a_filename,
             int a_numLevels,
             const Vector<int>& a_refRatios,
             const ProblemDomain& a_coarseDomain,
             const RealVect& a_origin,
             const RealVect& a_coarseDx,
             const RealVect& a_aspect,
             const IntVect& a_ghost);

void
headerEBFile(HDF5Handle& a_handle,
             int a_numLevels,
             const Vector<int>& a_refRatios,
             const ProblemDomain& a_coarseDomain,
             const RealVect& a_coarseDx,
             const IntVect& a_ghost);

void
headerEBFile(HDF5Handle& a_handle,
             int a_numLevels,
             const Vector<int>& a_refRatios,
             const ProblemDomain& a_coarseDomain,
             const RealVect& a_origin,
             const RealVect& a_coarseDx,
             const RealVect& a_aspect,
             const IntVect& a_ghost);

void
writeCellCenteredNames(HDF5Handle& a_handle,
                       const Vector<std::string>& a_names);
void
writeNodeCenteredNames(HDF5Handle& a_handle,
                       const Vector<std::string>& a_names);

void
writeFaceCenteredNames(HDF5Handle& a_handle, int a_direction,
                       const Vector<std::string>& a_names);

/// if no CellCenteredNames have been assigned, then dummy values of CC1..CCn are set

void
writeCellCentered(HDF5Handle& a_handle,
                  int a_level,
                  const LevelData<EBCellFAB>* a_data,
                  Interval interval = Interval());

void
writeCellCentered(HDF5Handle& a_handle,
                  int a_level,
                  const LevelData<FArrayBox>* a_data,
                  Interval interval = Interval());

void
writeNodeCentered(HDF5Handle& a_handle,
                  int a_level,
                  const LevelData<FArrayBox>* a_data,
                  Interval interval = Interval());

void
writeFaceCentered(HDF5Handle& a_handle,
                  int a_level,
                  int a_direction,
                  const LevelData<EBFaceFAB>* a_data,
                  Interval interval = Interval());

void
readCellCentered(HDF5Handle& a_handle,
                 int a_level,
                 const EBIndexSpace* eb,
                 const int ebghost,
                 LevelData<EBCellFAB>* a_data);

void
setWhichCellIndex(int a_whichCellIndex);

int
getWhichCellIndex();

void
setCoveredCellValue(Real a_coveredCellValue);

Real
getCoveredCellValue();

void multiFaceValues(const EBFaceFAB* a_face,
                     const int        a_side,
                     const int        a_iv0,
                     const int        a_iv1,
                     const int        a_iv2);

void multiCellValues(const EBCellFAB* a_face,
                     const int        a_iv0,
                     const int        a_iv1,
                     const int        a_iv2);

#endif  //CH_USE_HDF5

#include "NamespaceFooter.H"
#endif
