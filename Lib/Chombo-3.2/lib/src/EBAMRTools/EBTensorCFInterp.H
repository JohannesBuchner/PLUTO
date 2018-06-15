#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBTENSORCFINTERP_H_
#define _EBTENSORCFINTERP_H_

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
#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "QuadCFInterp.H"
#include "BaseIVFAB.H"
#include "TensorCFInterp.H"
#include "CornerCopier.H"
#include "EBCFData.H"
#include "EBQuadCFInterp.H"
#include "NamespaceHeader.H"

/// Quadratic coarse-fine interpolation utility for tensors
/**
    Class to fill ghost cells on coarse-fine interface using
    quadratic interpolation.  Also computes tangential
    derivatives
*/
class EBTensorCFInterp: public TensorCFInterp
{
public:

  ///
  EBTensorCFInterp(const DisjointBoxLayout& a_fineBoxes,
                   const DisjointBoxLayout& a_coarBoxes,
                   const EBISLayout&        a_ebislFine,
                   const EBISLayout&        a_ebislCoar,
                   const ProblemDomain&     a_domainCoar,
                   const int&               a_nref,
                   const int&               a_nvar,
                   const Real&              a_dxFine,
                   const LayoutData<IntVectSet>&  a_cfivs,
                   const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance(),
                   bool a_doEBCFCrossing = true);

  ///
  ~EBTensorCFInterp();

  ///
  /**
     Coarse / Fine interpolation operator.
  */
  void coarseFineInterp(LevelData<EBCellFAB>&       a_phif,
                        LevelData<EBCellFAB>&       a_tanGradf,
                        const LevelData<EBCellFAB>& a_phic);


  ///
  /**
     Coarse-fine interpolation operator with homogeneous BCs.
     Does same interpolation as coarseFineInterp function,
     but with all coarse-level data set to 0.
  */
  void coarseFineInterpH(LevelData<EBCellFAB>& a_phif,
                         LevelData<EBCellFAB>& a_tanGradf);


protected:

  void buildEBCFCornerStencils(  const LayoutData<IntVectSet>& a_cfivs);
  void buildEBCFCrossingStencils(const LayoutData<IntVectSet>& a_cfivs);


  void interpEBCFCrossing(LevelData<EBCellFAB>&       a_fineData,
                          LevelData<EBCellFAB>&       a_tanGradF,
                          const LevelData<EBCellFAB>& a_coarData);

  void interpEBCFCorners(LevelData<EBCellFAB>&       a_fineData,
                         LevelData<EBCellFAB>&       a_tanGradF,
                         const LevelData<EBCellFAB>& a_coarData);



  RefCountedPtr<EBCFData>        m_ebcfdata;
  RefCountedPtr<EBQuadCFInterp>  m_ebquadcfi;

  LevelData<EBCellFAB> m_ebBufferCoarsenedFine;

  LayoutData<BaseIVFAB<VoFStencil> > m_coarStencilLo[SpaceDim];
  LayoutData<BaseIVFAB<VoFStencil> > m_coarStencilHi[SpaceDim];
  LayoutData<BaseIVFAB<VoFStencil> > m_stencilCorners;
  LayoutData<BaseIVFAB<VoFStencil> > m_stencilEdges;

  //fine index space fine layout
  LayoutData<IntVectSet> m_ebcfivsLo[SpaceDim];
  LayoutData<IntVectSet> m_ebcfivsHi[SpaceDim];

  // Copier object to handle corner cells
  CornerCopier m_cornerCopier;

  bool m_doEBCFCrossing;
  int m_refRat;
  static IntVect s_ivDebFine;
  static IntVect s_ivDebCoar;
  ProblemDomain m_domainFine, m_domainCoar;


private:
  //weak construction is bad
  EBTensorCFInterp()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  EBTensorCFInterp(EBTensorCFInterp& a_input)
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(EBTensorCFInterp& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
