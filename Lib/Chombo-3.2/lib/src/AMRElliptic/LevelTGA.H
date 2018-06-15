#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELTGA_H_
#define _LEVELTGA_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include <REAL.H>
#include <Box.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <FArrayBox.H>
#include "FluxBox.H"
#include "LevelFluxRegister.H"
#include <ProblemDomain.H>
#include "BaseLevelTGA.H"
#include "BaseLevelBackwardEuler.H"
#include "NamespaceHeader.H"

/// Implements TGA algorithm to update solution to heat equation
/** The LevelTGA class implements the Runge-Kutta-based approach
    to solving the heat equation due to Twizell, Gumel, and Arigu.
*/
class LevelTGA : public BaseLevelTGA<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>
{

public:

  /// full constructor
  /**
   */
  LevelTGA(const Vector<DisjointBoxLayout>&                               a_grids,
           const Vector<int>&                                             a_refRat,
           const ProblemDomain&                                           a_level0Domain,
           RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >&      a_opFact,
           const RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >&     a_solver)
    :BaseLevelTGA<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver)
  {
    ;
  }

  /// destructor
  virtual ~LevelTGA()
  {
  }

  /// computes time-centered diffusive term for explicit updating
  /** compute time-centered L(phi) for use in subsequent update
      operations.  In this case, we do solve for phiNew, then
      subtract source and old phi back out to get L(phi).
  */
  void computeDiffusion(LevelData<FArrayBox>& a_DiffusiveTerm,
                        LevelData<FArrayBox>& a_phiOld,
                        LevelData<FArrayBox>& a_src,
                        LevelFluxRegister* a_FineFluxRegPtr,
                        LevelFluxRegister* a_crseFluxRegPtr,
                        const LevelData<FArrayBox>* a_crsePhiOldPtr,
                        const LevelData<FArrayBox>* a_crsePhiNewPtr,
                        Real a_oldTime,
                        Real a_crseOldTime,
                        Real a_crseNewTime,
                        Real a_dt,
                        int a_level);

  /// do solve
  /** update phi to a_phiNew
   */
  void updateSoln(LevelData<FArrayBox>& a_phiNew,
                  LevelData<FArrayBox>& a_phiOld,
                  LevelData<FArrayBox>& a_src,
                  LevelFluxRegister* a_fineFluxRegPtr,
                  LevelFluxRegister* a_crseFluxRegPtr,
                  const LevelData<FArrayBox>* a_crsePhiOldPtr,
                  const LevelData<FArrayBox>* a_crsePhiNewPtr,
                  Real a_oldTime,
                  Real a_crseOldTime,
                  Real a_crseNewTime,
                  Real a_dt,
                  int a_level,
                  bool a_zeroPhi = true,
                  int a_fluxStartComponent = 0);

  void setSourceGhostCells(LevelData<FArrayBox>&    a_src,
                           const DisjointBoxLayout& a_grids,
                           int a_level);
protected:

};

/// Implements BackwardEuler algorithm to update solution to heat equation
/** The LevelBackwardEuler class implements the Runge-Kutta-based approach
    to solving the heat equation due to Twizell, Gumel, and Arigu.
*/
class LevelBackwardEuler : public BaseLevelBackwardEuler<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>
{

public:

  /// full constructor
  /**
   */
  LevelBackwardEuler(const Vector<DisjointBoxLayout>&                               a_grids,
           const Vector<int>&                                             a_refRat,
           const ProblemDomain&                                           a_level0Domain,
           RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >&      a_opFact,
           const RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >&     a_solver)
    :BaseLevelBackwardEuler<LevelData<FArrayBox>, FluxBox, LevelFluxRegister>(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver)
  {
    ;
  }

  /// destructor
  virtual ~LevelBackwardEuler()
  {
  }

  /// computes time-centered diffusive term for explicit updating
  /** compute time-centered L(phi) for use in subsequent update
      operations.  In this case, we do solve for phiNew, then
      subtract source and old phi back out to get L(phi).
  */
  void computeDiffusion(LevelData<FArrayBox>& a_DiffusiveTerm,
                        LevelData<FArrayBox>& a_phiOld,
                        LevelData<FArrayBox>& a_src,
                        LevelFluxRegister* a_FineFluxRegPtr,
                        LevelFluxRegister* a_crseFluxRegPtr,
                        const LevelData<FArrayBox>* a_crsePhiOldPtr,
                        const LevelData<FArrayBox>* a_crsePhiNewPtr,
                        Real a_oldTime,
                        Real a_crseOldTime,
                        Real a_crseNewTime,
                        Real a_dt,
                        int a_level);

  /// do solve
  /** update phi to a_phiNew
   */
  void updateSoln(LevelData<FArrayBox>& a_phiNew,
                  LevelData<FArrayBox>& a_phiOld,
                  LevelData<FArrayBox>& a_src,
                  LevelFluxRegister* a_fineFluxRegPtr,
                  LevelFluxRegister* a_crseFluxRegPtr,
                  const LevelData<FArrayBox>* a_crsePhiOldPtr,
                  const LevelData<FArrayBox>* a_crsePhiNewPtr,
                  Real a_oldTime,
                  Real a_crseOldTime,
                  Real a_crseNewTime,
                  Real a_dt,
                  int a_level,
                  bool a_zeroPhi = true,
                  int a_fluxStartComponent = 0);

  void setSourceGhostCells(LevelData<FArrayBox>&    a_src,
                           const DisjointBoxLayout& a_grids,
                           int a_level);
protected:

};

#include "NamespaceFooter.H"

#endif
