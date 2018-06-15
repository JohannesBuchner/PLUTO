#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Qinghai Zhang, May 2009

#ifndef _NEWPOISSONOP4_H_
#define _NEWPOISSONOP4_H_

#include "AMRMultiGrid.H"
#include "REAL.H"
#include "Box.H"
#include "LevelDataOps.H"
#include "BCFunc.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

class NewPoissonOp4 : public MGLevelOp<FArrayBox>
{
public:
  NewPoissonOp4()
  {
  }

  virtual ~NewPoissonOp4()
  {
  }

  static const int m_nGhost = 2;

  /**
     \name LinearOp functions */
  /*@{*/

  void define(const RealVect& a_dx,
              const ProblemDomain& a_domain,
              BCFunc a_bc);

  virtual void createCoarsened(    FArrayBox& a_lhs,
                                   const FArrayBox& a_rhs,
                                   const int& a_refRat);

  virtual void residual(  FArrayBox& a_lhs,
                          const FArrayBox& a_phi,
                          const FArrayBox& a_rhs,
                          bool a_homogeneous = false);
  virtual void preCond(   FArrayBox& a_correction,
                          const FArrayBox& a_residual);
  virtual void applyOp(   FArrayBox& a_lhs,
                          const FArrayBox& a_phi,
                          bool a_homogeneous = false);
  virtual void create(    FArrayBox& a_lhs,
                          const FArrayBox& a_rhs);
  virtual void assign(    FArrayBox& a_lhs,
                          const FArrayBox& a_rhs) ;
  virtual Real dotProduct(const FArrayBox& a_1,
                          const FArrayBox& a_2) ;
  virtual void incr( FArrayBox& a_lhs,
                     const FArrayBox& a_x,
                     Real a_scale) ;
  virtual void axby( FArrayBox& a_lhs, const FArrayBox& a_x,
                     const FArrayBox& a_y,
                     Real a, Real b) ;
  virtual void scale(FArrayBox& a_lhs, const Real& a_scale) ;

  virtual Real norm(const FArrayBox& a_x, int a_ord);

  virtual void setToZero(FArrayBox& a_x);

  /*@}*/

  /**
     \name MGLevelOp functions */
  /*@{*/

  virtual void relax(FArrayBox& a_e,
                     const FArrayBox& a_residual,
                     int iterations);

  virtual void createCoarser(FArrayBox& a_coarse,
                             const FArrayBox& a_fine,
                             bool ghost);
  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(FArrayBox& a_resCoarse,
                                FArrayBox& a_phiFine,
                                const FArrayBox& a_rhsFine);

  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(FArrayBox& a_phiThisLevel,
                                const FArrayBox& a_correctCoarse);

  /*@}*/


protected:

  RealVect                m_dx;
  RealVect                m_dxCrse;
  ProblemDomain           m_domain;
  BCFunc                  m_bc;

  void levelGSRB(FArrayBox& a_e,  const FArrayBox& a_residual);

};

class NewPoissonOp4Factory : public MGLevelOpFactory<FArrayBox>
{
public:

  NewPoissonOp4Factory();

  NewPoissonOp4Factory(RealVect& a_dx, BCFunc a_bc);

  void define(const RealVect& m_dx, BCFunc m_bc);

  virtual ~NewPoissonOp4Factory()
  {
  }

  virtual NewPoissonOp4* MGnewOp(const ProblemDomain& a_FineindexSpace,
                                   int   a_depth,
                                   bool  a_homoOnly = true);

  virtual void MGreclaim(NewPoissonOp4* a_reclaim);

  RealVect                m_dx;
  BCFunc                  m_bc;

};

#include "NamespaceFooter.H"
#endif
