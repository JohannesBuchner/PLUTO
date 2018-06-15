#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBGRADDIVFILTER_H_
#define _EBGRADDIVFILTER_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>
#include <REAL.H>
#include <IntVect.H>
#include <Box.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <BaseFab.H>
#include "EBCellFAB.H"
#include "EBTensorCFInterp.H"
#include "CFIVS.H"
#include "EBFluxFAB.H"
#include "EBFaceFAB.H"
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"

///
/**
   Class to filter out divergence modes by imposing
   u += lambda * grad div u
 */
class EBGradDivFilter
{
public:


  ///
  /**
     Constructor.   Send coarse objects as NULL if they do not
     exist.
   */
  EBGradDivFilter(const DisjointBoxLayout&  a_fineGrids,
                  const DisjointBoxLayout*  a_coarGridsPtr,
                  const EBISLayout&         a_ebislFine,
                  const EBISLayout*         a_ebislCoarPtr,
                  const ProblemDomain&      a_domainFine,
                  const RealVect&           a_dxFine,
                  const int&                a_refRat,
                  const Real&               a_lambdaScale = 0.125,
                  const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());


  ///
  ~EBGradDivFilter();


  ///
  /**
     makes vel := vel + lambda*(grad div  u).
     velCoar = NULL if no coarser
     If lowOrderOneSidedGrad, the one sided differences
     in the gradient of the velocity are simply 2 point.
     Otherwise we do the higher-order 3 point difference.
     Default true here because it seems to add stability.
  */
  void filter(LevelData<EBCellFAB>&       a_velFine,
              const LevelData<EBFluxFAB>& a_fluxVelFine,
              const LevelData<EBCellFAB>* a_velCoar,
              bool  a_lowOrderOneSidedGrad = false,
              bool  a_noExtrapToCovered    = false);

  ///
  /**
     Returns grad(div(u)).
     Velfine is nonconst because of coarse-fine interpolation and
     exchange.      For convergence tests.
     NOT pre multiplied by lambda
     If lowOrderOneSidedGrad, the one sided differences
     in the gradient of the velocity are simply 2 point.
     Otherwise we do the higher-order 3 point difference.
     Default false here because this is used in convergence tests.
  */
  void
  gradDiv(LevelData<EBCellFAB>&       a_gradDivVel,
          LevelData<EBCellFAB>&       a_velFine,
          const LevelData<EBFluxFAB>& a_fluxVelFine,
          const LevelData<EBCellFAB>* a_velCoar,
          bool a_lowOrderOneSidedGrad = false,
          bool a_noExtrapToCovered    = false);

  ///
  Real getDomainDivergence(const  EBCellFAB&      a_gradVel,
                           const  EBCellFAB&      a_vel,
                           const  EBFluxFAB&      a_fluxVel,
                           const Box&             a_grid,
                           const EBISBox&         a_ebisBox,
                           const int&             a_faceDir,
                           const FaceIndex&       a_face,
                           const Side::LoHiSide&  a_side);

  ///
  /**
     computes  lambda*(grad div  u).
     Does one-side differences at boundary.
     pre multiplied by lambda if a_multiplyByLambda = true
  */
  void gradDiv(EBCellFAB&       a_gradDivVel,
               const EBCellFAB& a_gradVel,
               const EBCellFAB& a_velFine,
               const EBFluxFAB& a_fluxVelFine,
               const EBFluxFAB& a_divUFaceCent,
               const Box&       a_gridFine,
               const EBISBox&   a_ebisBoxFine,
               const DataIndex& a_dit,
               bool a_multiplyByLambda,
               bool a_noExtrapToCovered);

  ///
  /**
     Determines which component  of spacedim*spacedim fab
     the velocity gradient of the velocity goes into.
     velDir is the component of the velocity.
     derivdir is the compoent of the gradient.
   */
  int getGradComp(int a_velDir, int a_derivDir);

  ///
  /**
     computes cell centered gradient of cell centered velocity
     and puts it into gradvel fab of size spacedim*spacedim.
     getGradComp determines which component of fab.
     If lowOrderOneSidedGrad, the one sided differences
     in the gradient of the velocity are simply 2 point.
     Otherwise we do the higher-order 3 point difference.
     Default false here because this is used in convergence tests.
  */
  void
  gradVel(EBCellFAB&             a_gradVel,
          const  EBCellFAB&      a_velFine,
          const Box&             a_gridFine,
          const EBISBox&         a_ebisBoxFine,
          bool a_lowOrderOneSidedGrad);
  void
  faceDivergence(EBFaceFAB&             a_divVel,
                 const  EBCellFAB&      a_gradVel,
                 const  EBCellFAB&      a_velFine,
                 const EBFluxFAB&       a_fluxVelFine,
                 const Box&             a_gridFine,
                 const EBISBox&         a_ebisBoxFine,
                 const int&             a_faceDir);

  void
  cellGradient(EBCellFAB&             a_gradDivVel,
               const EBCellFAB&       a_gradVel,
               const  EBCellFAB&      a_velFine,
               const EBFluxFAB&       a_fluxVelFine,
               const EBFluxFAB&       a_divVel,
               const Box&             a_gridFine,
               const EBISBox&         a_ebisBoxFine,
               const DataIndex&       a_dit,
               bool a_multiplyByLambda,
               bool a_noExtrapToCovered);
protected:

  void getAreaFracs(Vector<FaceIndex> a_facesLo[SpaceDim],
                    Vector<FaceIndex> a_facesHi[SpaceDim],
                    bool              a_hasFacesLo[SpaceDim],
                    bool              a_hasFacesHi[SpaceDim],
                    RealVect&         a_areaFracLo,
                    RealVect&         a_areaFracHi,
                    const VolIndex&   a_vof,
                    const EBISBox&    a_ebisBox);

  void fillLambda();

  DisjointBoxLayout  m_gridsFine;
  EBISLayout         m_ebislFine;

  const DisjointBoxLayout* m_gridsCoarPtr;
  const EBISLayout*        m_ebislCoarPtr;
  EBLevelGrid m_eblgFine;
  ProblemDomain      m_domainCoar;
  ProblemDomain      m_domainFine;
  int                m_refRat;
  bool               m_hasCoarser;
  EBTensorCFInterp*  m_tensorCFI;
  RealVect           m_dxFine;

  //contains tangential gradient info
  LevelData<EBCellFAB> m_gradVel;
  LevelData<EBFluxFAB> m_faceDivCell;
  LevelData<EBFluxFAB> m_faceDivCent;
  LayoutData<BaseIVFAB<VoFStencil> > m_johanStencil;
  LayoutData<BaseIVFAB<char>  > m_dropOrder;
  LayoutData<BaseIVFAB<Real>  > m_distanceAlongLine;
  LevelData<EBCellFAB>          m_lambda;

  Real m_lambdaScale;

private:

  //weak construction is bad
  EBGradDivFilter()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  EBGradDivFilter(const EBGradDivFilter& a_input)
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBGradDivFilter& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
