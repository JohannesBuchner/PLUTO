#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBDEBUGOUT_H_
#define _EBDEBUGOUT_H_


#include "DisjointBoxLayout.H"
#include "BoxLayout.H"
#include "BaseIFFAB.H"
#include "BaseIVFAB.H"
#include "FaceIterator.H"
#include "Vector.H"
#include "Box.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "VolIndex.H"
#include "Stencils.H"
#include "NamespaceHeader.H"

class EBDebugPoint
{
public:
  static IntVect s_ivd;
};

class EBCellFAB;
class EBFluxFAB;
//extern "C"
//{
///
///
void printMaxMinEBFF(EBFluxFAB* a_data);

///
void printMaxMinEBCF(EBCellFAB* a_data);

///
void printMaxMinLDCell(LevelData<EBCellFAB>* a_data);

///
void printMaxMinLDFlux(LevelData<EBFluxFAB>* a_data);

///
void
printPointEBCF(EBCellFAB* a_dat);

void
printLocalEBCF(EBCellFAB* a_dat);

void
printPointEBFace(EBFaceFAB* a_dat);

///
void printPointLDCell(LevelData<EBCellFAB>* a_dat);

///
void printLocalLDCell(LevelData<EBCellFAB>* a_dat);

///
void printPointEBFF(EBFluxFAB* a_dat);

///
void printPointLDFlux(LevelData<EBFluxFAB>* a_dat);

///
void
dumpLDEBCF(const LevelData<EBCellFAB>*  ldptr);

///
void
getMaxEBLevel(const LevelData<EBCellFAB>*  ldptr);

///
void
getMaxEBFAB(const EBCellFAB*  ldptr);


/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpVVoF(const Vector<VolIndex>* a_vectPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpVFace(const Vector<FaceIndex>* a_vectPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpFace(const FaceIndex* a_vectPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpFaceSten(const FaceStencil* a_stenPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpVoFSten(const   VoFStencil* a_stenPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpIFFAB(const BaseIFFAB<Real>* a_stenPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpIVFAB(const BaseIVFAB<Real>* a_stenPtr);


///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpLDBIVF(const LayoutData< BaseIVFAB<Real> >* a_stenPtr);

///
/**
 */
void dumpEBLDDBL(const LevelData<EBCellFAB>* a_stenPtr);

///
/**
 */
void dumpEBAMRThresh(const Vector<LevelData<EBCellFAB>*>* a_stenPtr, Real a_minDump);

///
/**
 */
void dumpEBAMRIrregThresh(const Vector<LevelData<EBCellFAB>*>* a_stenPtr, Real a_minDump);


///
/**
 */
void dumpEBAMR(const Vector<LevelData<EBCellFAB>*>* a_stenPtr);

///
/**
 */
void dumpEBAMRIrreg(const Vector<LevelData<EBCellFAB>*>* a_stenPtr);

///
/**
   dumps the fab to std out.
*/
void dumpEBFAB(const EBCellFAB* a_fab);

///
/**
   dumps the fab to std out.
*/
void dumpEBFaceIVS(const EBFaceFAB* a_fab, const IntVectSet& a_ivs, Real a_thresh);

///
void dumpEBFaceThresh(const EBFaceFAB* a_fab, Real a_thresh);

///
void dumpEBFace(const EBFaceFAB* a_fab);

///
void dumpEBFlux(const EBFluxFAB* a_fab);

///
void dumpEBFluxThresh(const EBFluxFAB* a_fab, Real a_thresh);

///
void dumpEBLevelFlux(const LevelData<EBFluxFAB>* a_fab);

///
void dumpEBLevelFluxThresh(const LevelData<EBFluxFAB>* a_fab, Real a_thresh);

///
/**
   dumps the irreg bits of the fab to std out.
*/
void dumpEBFABIrreg(const EBCellFAB* a_fab);

///
/**
   dumps the geometry data of the irreg bits of the fab to std out.
*/
void dumpEBFABIrregGeometry(const EBCellFAB* a_fab);

///
/**
   dumps the irreg bits of the level to std out.
*/
void dumpEBLevelIrreg(const LevelData<EBCellFAB>* a_fab);


///
/**
   dumps the irreg bits of the level to std out.
*/
void dumpEBLevelIrregThresh(const LevelData<EBCellFAB>* a_fab, Real a_thresh);

///
/**
*/
void dumpEBFABIVS(const EBCellFAB* a_fab, const IntVectSet* a_ivs, Real a_thresh = 0);

///
/**
*/
void dumpEBLevel(const LevelData<EBCellFAB>* a_level);

///same as above with ghost cells too
void dumpEBLevelAll(const LevelData<EBCellFAB>* a_level);


///
/**
*/
void dumpEBLevelThresh(const LevelData<EBCellFAB>* a_level, Real a_thresh);

//}

#include "NamespaceFooter.H"
#endif
