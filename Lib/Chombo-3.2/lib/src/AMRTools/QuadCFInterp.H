#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _QUADCFINTERP_H_
#define _QUADCFINTERP_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>
#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "QuadCFStencil.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

class CFIVS;
class TensorFineStencilSet;
///Quadratic coarse-fine interpolation utility
/**
    Class to fill ghost cells on coarse-fine
    interface using quadratic interpolation
*/
class QuadCFInterp
{
public:

  ///
  /**
   */
  static void
  homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                      LevelData<FArrayBox>& a_tanGrad,
                      LayoutData<CFIVS> a_loCFIVS[SpaceDim],
                      LayoutData<CFIVS> a_hiCFIVS[SpaceDim],
                      Real a_dxLevel,
                      Real a_dxCrse,
                      int a_ncomp,
                      LayoutData<TensorFineStencilSet> a_loTanStencilSets[SpaceDim],
                      LayoutData<TensorFineStencilSet> a_hiTanStencilSets[SpaceDim]);


  ///
  /**
   */
  static void
  homogeneousCFInterpTanGrad(LevelData<FArrayBox>& a_tanGrad,
                             const LevelData<FArrayBox>& a_phi,
                             const DataIndex& a_DatInd,
                             int a_idir,
                             Side::LoHiSide a_hiorlo,
                             Real a_dxLevel,
                             Real a_dxCrse,
                             int a_ncomp,
                             LayoutData<TensorFineStencilSet> a_loTanStencilSets[SpaceDim],
                             LayoutData<TensorFineStencilSet> a_hiTanStencilSets[SpaceDim]);


  ///
  /**
   */
  static void
  homogeneousCFInterpPhi(LevelData<FArrayBox>& a_phif,
                         const DataIndex& a_datInd,
                         int a_idir,
                         Side::LoHiSide a_hiorlo,
                         LayoutData<CFIVS> a_loCFIVS[SpaceDim],
                         LayoutData<CFIVS> a_hiCFIVS[SpaceDim],
                         Real a_dxLevel,
                         Real a_dxCrse,
                         int a_ncomp);

  static void
  interpPhiOnIVS(LevelData<FArrayBox>& a_phif,
                 const FArrayBox& a_phistar,
                 const DataIndex& a_datInd,
                 const int a_idir,
                 const Side::LoHiSide a_hiorlo,
                 const IntVectSet& a_interpIVS,
                 Real a_dxLevel,
                 Real a_dxCrse,
                 int a_ncomp);

  ///constructors, destructors, defines

  ///
  /**
     Full Constructor.  makes all coarse-fine
     information and sets internal variables
     calls full define
  */
  QuadCFInterp(
               const DisjointBoxLayout& a_fineBoxes,
               const DisjointBoxLayout* a_coarBoxes,
               Real  a_dxFine,
               int a_refRatio,
               int a_nComp,
               const Box& a_domf);

  ///
  /**
     Full Constructor.  makes all coarse-fine
     information and sets internal variables
     calls full define
  */
  QuadCFInterp(
               const DisjointBoxLayout& a_fineBoxes,
               const DisjointBoxLayout* a_coarBoxes,
               Real  a_dxFine,
               int a_refRatio,
               int a_nComp,
               const ProblemDomain& a_domf);

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables
  */
  void define(
              const DisjointBoxLayout& a_fineBoxes,
              const DisjointBoxLayout* a_coarBoxes,
              Real  a_dxFine,
              int a_refRatio,
              int a_nComp,
              const ProblemDomain& a_domf);

  ///
  /**
     Default constructor leaves QCFI undefined.
  */
  QuadCFInterp();

  ///
  virtual ~QuadCFInterp();

  ///
  /**
     return QCFI to undefined state
  */
  void clear();

  ///
  /**
     Coarse / Fine interpolation operator.
  */
  void coarseFineInterp(
                        LevelData<FArrayBox>& a_phif,
                        const LevelData<FArrayBox>& a_phic
                        );

  ///
  /**
     has full define function been called?  return true if so
  */
  bool isDefined() const;

  static bool newCFInterMode;

protected:
  //boxes for coarse buffers
  BoxLayout m_coarBoxes;

  //fake level identifier.
  int m_level;

  //number of components
  int m_nComp;

  //refinement ratio between levels
  int m_refRatio;

  //mesh spacing at fine level
  Real m_dxFine;

  //cfstencils on low-side faces
  //lives on fine input grid layout
  LayoutData<QuadCFStencil> m_loQCFS[SpaceDim];

  //cfstencils on high-side faces
  //lives on fine input grid layout
  LayoutData<QuadCFStencil> m_hiQCFS[SpaceDim];


  BoxLayoutData<FArrayBox> m_coarBuffer;
  Copier m_copier;

  DisjointBoxLayout m_inputFineLayout;
  DisjointBoxLayout m_inputCoarLayout;
  ///has full define function been called?
  bool m_isDefined;
  bool m_fineCoversCoarse;
  ProblemDomain m_domainFine;
protected:

  /* Coarse / Fine interpolation operator.
     using one sided differencing
     only at a particular face
  */
  void coarseFineInterp(BaseFab<Real> & a_phif,
                        const BaseFab<Real> & a_phic,
                        const QuadCFStencil& a_qcfs,
                        const Side::LoHiSide a_hiorlo,
                        const int a_idir,
                        const Interval& a_variables) const;

  //get extended phi (lives next to interpivs)
  void getPhiStar(BaseFab<Real> & a_phistar,
                  const BaseFab<Real> & a_phic,
                  const QuadCFStencil& a_qcfs,
                  const Side::LoHiSide a_hiorlo,
                  const int a_idir,
                  const Interval& a_variables) const;

  //interpolate over correct intvectset
  void interpOnIVS(BaseFab<Real> & a_phif,
                   const BaseFab<Real> & a_phiStar,
                   const QuadCFStencil& a_qcfs,
                   const Side::LoHiSide a_hiorlo,
                   const int a_idir,
                   const Interval& a_variables) const;
};

#include "NamespaceFooter.H"
#endif
