#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TENSORCFINTERP_H_
#define _TENSORCFINTERP_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "BaseFab.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "QuadCFStencil.H"
#include "ProblemDomain.H"
#include "TensorFineStencilSet.H"
#include "NamespaceHeader.H"

/// Quadratic coarse-fine interpolation utility for tensors
/**
    Class to fill ghost cells on coarse-fine interface using
    quadratic interpolation.  Also computes tangential
    derivatives
*/
class TensorCFInterp
{
public:

  ///constructors, destructors, defines

  ///
  /**
     Full Constructor.  makes all coarse-fine
     information and sets internal variables
     calls full define
  */
  TensorCFInterp(
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
  TensorCFInterp(
               const DisjointBoxLayout& a_fineBoxes,
               const DisjointBoxLayout* a_coarBoxes,
               Real  a_dxFine,
               int a_refRatio,
               int a_nComp,
               const ProblemDomain& a_domf);

  ///
  /**
     Full define function.  Makes all coarse-fine
     information and sets internal variables
  */
  void define(
              const DisjointBoxLayout& a_fineBoxes,
              const DisjointBoxLayout* a_coarBoxes,
              Real  a_dxFine,
              int a_refRatio,
              int a_nComp,
              const Box& a_domf);

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
     Default constructor leaves TCFI undefined.
  */
  TensorCFInterp();

  ///
  virtual ~TensorCFInterp();

  ///
  /**
     return TCFI to undefined state
  */
  void clear();

  ///
  /**
     Coarse / Fine interpolation operator.
  */
  void coarseFineInterp(
                        LevelData<FArrayBox>& a_phif,
                        LevelData<FArrayBox>& a_gradf,
                        const LevelData<FArrayBox>& a_phic
                        );

  ///
  /**
     Coarse-fine interpolation operator with homogeneous BCs.
     Does same interpolation as coarseFineInterp function,
     but with all coarse-level data set to 0.
  */
  void coarseFineInterpH(
                         LevelData<FArrayBox>& a_phif,
                         LevelData<FArrayBox>& a_gradf
                         );

  ///
  /**
     has full define function been called?  return true if so
  */
  bool isDefined() const;

  /// returns component in gradient for variable ivar in direction gradDir
  /** returns ivar*SpaceDim + gradDir
   */
  static inline int gradIndex(int ivar, int gradDir)
  {
    return ivar*SpaceDim + gradDir;
  }

protected:
  DisjointBoxLayout m_coarsenedFineBoxes;

  //fake level identifier.
  int m_level;

  //number of components
  int m_nComp;

  //refinement ratio between levels
  int m_refRatio;

  //mesh spacing at fine level
  Real m_dxFine;

  // problem domain at fine level
  ProblemDomain m_probDomain;

  //cfstencils on low-side faces
  //lives on fine input grid layout
  LayoutData<QuadCFStencil> m_loQCFS[SpaceDim];

  //cfstencils on high-side faces
  //lives on fine input grid layout
  LayoutData<QuadCFStencil> m_hiQCFS[SpaceDim];


  /// this contains info for the fine-level differencing
  LayoutData<TensorFineStencilSet> m_loFineStencilSets[SpaceDim];

  ///
  LayoutData<TensorFineStencilSet> m_hiFineStencilSets[SpaceDim];

  LevelData<FArrayBox> m_coarsenedFineBuffer;

  DisjointBoxLayout m_inputFineLayout;
  DisjointBoxLayout m_inputCoarLayout;
  ///has full define function been called?
  bool m_isDefined;
  bool m_fineCoversCoarse;
protected:

  /* Coarse / Fine interpolation operator.
     using one sided differencing
     only at a particular face
  */
  void coarseFineInterp(BaseFab<Real> & a_phif,
                        BaseFab<Real> & a_gradf,
                        BaseFab<Real> & a_tanGradStar,
                        const BaseFab<Real> & a_phic,
                        const QuadCFStencil& a_qcfs,
                        const Side::LoHiSide a_hiorlo,
                        const int a_idir,
                        const Interval& a_variables) const;

  //get extended phi (lives next to interpivs)
  void getPhiStar(BaseFab<Real> & a_phistar,
                  BaseFab<Real> & a_tanGradStar,
                  const BaseFab<Real> & a_phic,
                  const QuadCFStencil& a_qcfs,
                  const Side::LoHiSide a_hiorlo,
                  const int a_idir,
                  const Interval& a_variables) const;

  //interpolate over correct intvectset
  void interpOnIVS(BaseFab<Real> & a_phif,
                   BaseFab<Real> & a_gradf,
                   const BaseFab<Real> & a_phiStar,
                   const QuadCFStencil& a_qcfs,
                   const Side::LoHiSide a_hiorlo,
                   const int a_idir,
                   const Interval& a_variables) const;

  // compute tangential gradients
  void computeTanGrad(BaseFab<Real>& a_gradf,
                      const BaseFab<Real>& a_phiFine,
                      const BaseFab<Real>& a_tanGradStar,
                      const TensorFineStencilSet& a_fineStencils,
                      const QuadCFStencil& a_qcfs,
                      const Side::LoHiSide& a_hiorlo,
                      const int a_idir,
                      const Interval& a_variables) const;
};

#include "NamespaceFooter.H"
#endif
