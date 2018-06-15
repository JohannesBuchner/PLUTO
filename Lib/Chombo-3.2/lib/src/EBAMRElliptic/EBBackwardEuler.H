#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBBACKWARDEULER_H_
#define _EBBACKWARDEULER_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AMRMultiGrid.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "AMRTGA.H"
#include "NamespaceHeader.H"

///
/**
   Solves dphi/dt = L phi + rho
   using backward euler.
   Uses TGAHelmOp functionaltity.
 **/
class EBBackwardEuler
{
public:


  ///
  ~EBBackwardEuler();

  ///
  /**
   **/
  EBBackwardEuler(const RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > > & a_solver,
                  const AMRLevelOpFactory<LevelData<EBCellFAB> >&              a_factory,
                  const ProblemDomain&                                         a_level0Domain,
                  const Vector<int>&                                           a_refRat,
                  int a_numLevels = -1, int a_verbosity = 3);


  ///
  /**
     This advances a parabolic pde from a_phiOld to a_phiNew using TGA on a non-moving domain with source term a_source
  **/
  void oneStep(Vector<LevelData<EBCellFAB>*>&             a_phiNew,
               Vector<LevelData<EBCellFAB>*>&             a_phiOld,
               Vector<LevelData<EBCellFAB>*>&             a_source,
               const Real&             a_dt,
               int                     a_lbase,
               int                     a_lmax,
               bool                    a_zeroPhi = true,
               const bool              a_kappaWeighted = false);

  ///
  void resetAlphaAndBeta(const Real& a_alpha,
                         const Real& a_beta);

  void
  computeDiffusion(Vector< LevelData<EBCellFAB>* >&       a_diffusiveTerm,
                   Vector< LevelData<EBCellFAB>* >&       a_phiOld,
                   Vector< LevelData<EBCellFAB>* >&       a_src,
                   Real a_oldTime,  Real a_dt,
                   int a_lbase, int a_lmax, bool a_zeroPhi);


protected:
  void solveHelm(Vector<LevelData<EBCellFAB>* >&       a_ans,
                 Vector<LevelData<EBCellFAB>* >&       a_rhs,
                 int               a_lbase,
                 int               a_lmax,
                 Real              a_dt,
                 bool              a_zeroPhi);

  //fills a_ans = dt*kappa*a_source + dt*phiOld
  void createEulerRHS(Vector<LevelData<EBCellFAB>* >&       a_ans,
                      Vector<LevelData<EBCellFAB>* >&       a_source,
                      Vector<LevelData<EBCellFAB>* >&       a_phiOld,
                      int               a_lbase,
                      int               a_lmax,
                      Real              a_dt,
                      const bool        a_kappaWeighted = false);
  void createData(Vector<LevelData<EBCellFAB>* >&       a_source,
                  int                a_lbase,
                  int                a_lmax);

  TGAHelmOp<LevelData<EBCellFAB> >*
  newOp(const ProblemDomain&                             a_indexSpace,
        const AMRLevelOpFactory<LevelData<EBCellFAB> >&  a_opFact);

private:
  //You do not own these operators!!  don't delete it.   the memory is
  //owned by the solver
  Vector<TGAHelmOp<LevelData<EBCellFAB> > * >                m_ops;
  Vector< LevelData<EBCellFAB>* >                            m_rhst;
  ProblemDomain                                              m_level0Domain;
  Vector<int>                                                m_refRat;
  RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >       m_solver;
  int m_verbosity, m_numLevels;
  bool                                           m_dataCreated;

  //copy constructor and operator= disallowed for all the usual reasons
  EBBackwardEuler(const EBBackwardEuler& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBBackwardEuler& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  /// weak construction is bad.   Ref Counted pointers are your friends.
  EBBackwardEuler()
  {
    MayDay::Error("invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
