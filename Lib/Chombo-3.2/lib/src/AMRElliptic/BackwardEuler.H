#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BACKWARDEULER_H_
#define _BACKWARDEULER_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AMRMultiGrid.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "AMRTGA.H"
#include "NamespaceHeader.H"

///
/**
   Solves dphi/dt = L phi + rho
   using backward euler.
   Uses TGAHelmOp functionaltity.
 **/
class BackwardEuler
{
public:
  ///
  ~BackwardEuler();

  ///
  /**
   **/
  BackwardEuler(const RefCountedPtr<AMRMultiGrid< LevelData<FArrayBox> > > & a_solver,
                const AMRLevelOpFactory<LevelData<FArrayBox> >&              a_factory,
                const ProblemDomain&                                         a_level0Domain,
                const Vector<int>&                                           a_refRat,
                int a_numLevels = -1, int a_verbosity = 3);

  ///
  /**
     This advances a parabolic pde from a_phiOld to a_phiNew using TGA on a non-moving domain with source term a_source
  **/
  void oneStep(Vector<LevelData<FArrayBox>*>&             a_phiNew,
               Vector<LevelData<FArrayBox>*>&             a_phiOld,
               Vector<LevelData<FArrayBox>*>&             a_source,
               const Real&             a_dt,
               int                     a_lbase,
               int                     a_lmax,
               bool                    a_zeroPhi = true);

  ///
  void resetAlphaAndBeta(const Real& a_alpha,
                         const Real& a_beta);

  void
  computeDiffusion(Vector< LevelData<FArrayBox>* >&       a_diffusiveTerm,
                   Vector< LevelData<FArrayBox>* >&       a_phiOld,
                   Vector< LevelData<FArrayBox>* >&       a_src,
                   Real a_oldTime,  Real a_dt,
                   int a_lbase, int a_lmax, bool a_zeroPhi);


protected:
  void solveHelm(Vector<LevelData<FArrayBox>* >&       a_ans,
                 Vector<LevelData<FArrayBox>* >&       a_rhs,
                 int               a_lbase,
                 int               a_lmax,
                 Real              a_dt,
                 bool              a_zeroPhi);

  //fills a_ans = dt*kappa*a_source + dt*phiOld
  void createEulerRHS(Vector<LevelData<FArrayBox>* >&       a_ans,
                      Vector<LevelData<FArrayBox>* >&       a_source,
                      Vector<LevelData<FArrayBox>* >&       a_phiOld,
                      int               a_lbase,
                      int               a_lmax,
                      Real              a_dt);
  void createData(Vector<LevelData<FArrayBox>* >&       a_source,
                  int                a_lbase,
                  int                a_lmax);

  TGAHelmOp<LevelData<FArrayBox> >*
  newOp(const ProblemDomain&                             a_indexSpace,
        const AMRLevelOpFactory<LevelData<FArrayBox> >&  a_opFact);

private:
  //You do not own these operators!!  don't delete it.   the memory is
  //owned by the solver
  Vector<TGAHelmOp<LevelData<FArrayBox> > * >                m_ops;
  Vector< LevelData<FArrayBox>* >                            m_rhst;
  ProblemDomain                                              m_level0Domain;
  Vector<int>                                                m_refRat;
  RefCountedPtr<AMRMultiGrid< LevelData<FArrayBox> > >       m_solver;
  int m_verbosity, m_numLevels;
  bool                                           m_dataCreated;

  //copy constructor and operator= disallowed for all the usual reasons
  BackwardEuler(const BackwardEuler& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const BackwardEuler& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  /// weak construction is bad.   Ref Counted pointers are your friends.
  BackwardEuler()
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"

#endif
