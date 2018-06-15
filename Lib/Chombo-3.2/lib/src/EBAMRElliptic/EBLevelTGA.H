#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLEVELTGA_H_
#define _EBLEVELTGA_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "EBFluxFAB.H"
#include "EBFluxRegister.H"
#include "ProblemDomain.H"
#include "BaseLevelTGA.H"
#include "BaseLevelBackwardEuler.H"
#include "BaseLevelCrankNicolson.H"
#include "NamespaceHeader.H"


//! \class EBLevelTGA
//! This class implements an imlplicit integrator for a diffusion
//! equation using an approach devised by Twizzell, Gumel, and Arigu.
class EBLevelTGA : public BaseLevelTGA<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>
{

public:

  //! Create a new TGA level integrator.
  //! \param a_grids The disjoint box layout on which the level integrator is defined.
  //! \param a_refRat The refinement ratios for the boxes.
  //! \param a_level0Domain The coarsest grid level defining the problem domain.
  //! \param a_opFact A factory that generates level diffusion operators.
  //! \param a_solver A multigrid solver.
  EBLevelTGA(const Vector<DisjointBoxLayout>&                               a_grids,
             const Vector<int>&                                             a_refRat,
             const ProblemDomain&                                           a_level0Domain,
             RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&      a_opFact,
             const RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >&     a_solver)
    :BaseLevelTGA<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver)
  {
    m_isEBLGSet = false;
  }

  //! Destructor
  virtual ~EBLevelTGA()
  {
  }

  //! Computes the time-centered diffusive term for explicit updating.
  //! This integrates the diffusion equation for a new value of a quantity phi
  //! and then computes the diffusion term using a finite difference stencil
  //! in time.
  //! \param a_diffusionTerm The term computed by this method.
  //! \param a_phiOld The value of phi at the beginning of the time step,
  //!                 at the current grid level.
  //! \param a_src The source term in the diffusion equation at the current
  //!              grid level.
  //! \param a_fineFluxRegPtr A pointer to the flux register at the finer
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crseFluxRegPtr A pointer to the flux register at the coarser
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crsePhiOldPtr A pointer to the LevelData object holding the
  //!                        old value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_crsePhiNewPtr A pointer to the LevelData object holding the
  //!                        updated value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_oldTime The time at the beginning of the integration step on
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_dt The time step on the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If this flag is true, the initial estimate of phi will
  //!                  be set to zero. Otherwise, \a a_phiOld will be used.
  void computeDiffusion(LevelData<EBCellFAB>& a_DiffusiveTerm,
                        LevelData<EBCellFAB>& a_phiOld,
                        LevelData<EBCellFAB>& a_src,
                        EBFluxRegister* a_FineFluxRegPtr,
                        EBFluxRegister* a_crseFluxRegPtr,
                        const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                        const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                        Real a_oldTime,
                        Real a_crseOldTime,
                        Real a_crseNewTime,
                        Real a_dt,
                        int a_level,
                        bool a_zeroPhi = true,
                        bool a_rhsAlreadyKappaWeighted = false
                        );

  //! Integrates the diffusion equation, storing the result in \a a_phiNew.
  //! \param a_phiNew The new value of phi at the current grid level.
  //! \param a_phiOld The value of phi at the beginning of the time step,
  //!                 at the current grid level.
  //! \param a_src The source term in the diffusion equation at the current
  //!              grid level.
  //! \param a_fineFluxRegPtr A pointer to the flux register at the finer
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crseFluxRegPtr A pointer to the flux register at the coarser
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crsePhiOldPtr A pointer to the LevelData object holding the
  //!                        old value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_crsePhiNewPtr A pointer to the LevelData object holding the
  //!                        updated value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_oldTime The time at the beginning of the integration step on
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_dt The time step on the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If this flag is true, the initial estimate of phi will
  //!                  be set to zero. Otherwise, \a a_phiOld will be used.
  void updateSoln(LevelData<EBCellFAB>& a_phiNew,
                  LevelData<EBCellFAB>& a_phiOld,
                  LevelData<EBCellFAB>& a_src,
                  EBFluxRegister* a_fineFluxRegPtr,
                  EBFluxRegister* a_crseFluxRegPtr,
                  const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                  const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                  Real a_oldTime,
                  Real a_crseOldTime,
                  Real a_crseNewTime,
                  Real a_dt,
                  int a_level,
                  bool a_zeroPhi = true,
                  bool a_rhsAlreadyKappaWeighted = false,
                  int a_fluxStartComponent = 0);

  void setSourceGhostCells(LevelData<EBCellFAB>&    a_src,
                           const DisjointBoxLayout& a_grids,
                           int a_lev);

  void setEBLG(Vector<EBLevelGrid>& a_eblg)
  {
    m_isEBLGSet = true;
    m_eblg      = a_eblg;
  }

protected:
  bool m_isEBLGSet;
  Vector<EBLevelGrid> m_eblg;


};


//! \class EBLevelBackwardEuler
//! This class implements an imlplicit integrator for a diffusion
//! equation using an approach devised by Twizzell, Gumel, and Arigu.
class EBLevelBackwardEuler : public BaseLevelBackwardEuler<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>
{

public:

  //! Create a new BackwardEuler level integrator.
  //! \param a_grids The disjoint box layout on which the level integrator is defined.
  //! \param a_refRat The refinement ratios for the boxes.
  //! \param a_level0Domain The coarsest grid level defining the problem domain.
  //! \param a_opFact A factory that generates level diffusion operators.
  //! \param a_solver A multigrid solver.
  EBLevelBackwardEuler(const Vector<DisjointBoxLayout>&                               a_grids,
                       const Vector<int>&                                             a_refRat,
                       const ProblemDomain&                                           a_level0Domain,
                       RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&      a_opFact,
                       const RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >&     a_solver)
    :BaseLevelBackwardEuler<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver)
  {
    m_isEBLGSet = false;
  }

  //! Destructor
  virtual ~EBLevelBackwardEuler()
  {
  }

  //! Computes the time-centered diffusive term for explicit updating.
  //! This integrates the diffusion equation for a new value of a quantity phi
  //! and then computes the diffusion term using a finite difference stencil
  //! in time.
  //! \param a_diffusionTerm The term computed by this method.
  //! \param a_phiOld The value of phi at the beginning of the time step,
  //!                 at the current grid level.
  //! \param a_src The source term in the diffusion equation at the current
  //!              grid level.
  //! \param a_fineFluxRegPtr A pointer to the flux register at the finer
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crseFluxRegPtr A pointer to the flux register at the coarser
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crsePhiOldPtr A pointer to the LevelData object holding the
  //!                        old value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_crsePhiNewPtr A pointer to the LevelData object holding the
  //!                        updated value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_oldTime The time at the beginning of the integration step on
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_dt The time step on the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If this flag is true, the initial estimate of phi will
  //!                  be set to zero. Otherwise, \a a_phiOld will be used.
  void computeDiffusion(LevelData<EBCellFAB>& a_DiffusiveTerm,
                        LevelData<EBCellFAB>& a_phiOld,
                        LevelData<EBCellFAB>& a_src,
                        EBFluxRegister* a_FineFluxRegPtr,
                        EBFluxRegister* a_crseFluxRegPtr,
                        const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                        const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                        Real a_oldTime,
                        Real a_crseOldTime,
                        Real a_crseNewTime,
                        Real a_dt,
                        int a_level,
                        bool a_zeroPhi = true,
                        bool a_rhsAlreadyKappaWeighted = false
                        );

  //! Integrates the diffusion equation, storing the result in \a a_phiNew.
  //! \param a_phiNew The new value of phi at the current grid level.
  //! \param a_phiOld The value of phi at the beginning of the time step,
  //!                 at the current grid level.
  //! \param a_src The source term in the diffusion equation at the current
  //!              grid level.
  //! \param a_fineFluxRegPtr A pointer to the flux register at the finer
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crseFluxRegPtr A pointer to the flux register at the coarser
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crsePhiOldPtr A pointer to the LevelData object holding the
  //!                        old value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_crsePhiNewPtr A pointer to the LevelData object holding the
  //!                        updated value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_oldTime The time at the beginning of the integration step on
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_dt The time step on the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If this flag is true, the initial estimate of phi will
  //!                  be set to zero. Otherwise, \a a_phiOld will be used.
  void updateSoln(LevelData<EBCellFAB>& a_phiNew,
                  LevelData<EBCellFAB>& a_phiOld,
                  LevelData<EBCellFAB>& a_src,
                  EBFluxRegister* a_fineFluxRegPtr,
                  EBFluxRegister* a_crseFluxRegPtr,
                  const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                  const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                  Real a_oldTime,
                  Real a_crseOldTime,
                  Real a_crseNewTime,
                  Real a_dt,
                  int a_level,
                  bool a_zeroPhi = true,
                  bool a_rhsAlreadyKappaWeighted = false,
                  int a_fluxStartComponent = 0);


  void setEBLG(Vector<EBLevelGrid>& a_eblg)
  {
    m_isEBLGSet = true;
    m_eblg      = a_eblg;
  }

protected:
  bool m_isEBLGSet;
  Vector<EBLevelGrid> m_eblg;


};

//! \class EBLevelCrankNicolson
//! This class implements an imlplicit integrator for a diffusion
//! equation using an approach devised by Twizzell, Gumel, and Arigu.
class EBLevelCrankNicolson : public BaseLevelCrankNicolson<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>
{

public:

  //! Create a new CrankNicolson level integrator.
  //! \param a_grids The disjoint box layout on which the level integrator is defined.
  //! \param a_refRat The refinement ratios for the boxes.
  //! \param a_level0Domain The coarsest grid level defining the problem domain.
  //! \param a_opFact A factory that generates level diffusion operators.
  //! \param a_solver A multigrid solver.
  EBLevelCrankNicolson(const Vector<DisjointBoxLayout>&                               a_grids,
                       const Vector<int>&                                             a_refRat,
                       const ProblemDomain&                                           a_level0Domain,
                       RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&      a_opFact,
                       const RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >&     a_solver)
    :BaseLevelCrankNicolson<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver)
  {
    m_isEBLGSet = false;
  }

  //! Destructor
  virtual ~EBLevelCrankNicolson()
  {
  }

  //! Computes the time-centered diffusive term for explicit updating.
  //! This integrates the diffusion equation for a new value of a quantity phi
  //! and then computes the diffusion term using a finite difference stencil
  //! in time.
  //! \param a_diffusionTerm The term computed by this method.
  //! \param a_phiOld The value of phi at the beginning of the time step,
  //!                 at the current grid level.
  //! \param a_src The source term in the diffusion equation at the current
  //!              grid level.
  //! \param a_fineFluxRegPtr A pointer to the flux register at the finer
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crseFluxRegPtr A pointer to the flux register at the coarser
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crsePhiOldPtr A pointer to the LevelData object holding the
  //!                        old value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_crsePhiNewPtr A pointer to the LevelData object holding the
  //!                        updated value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_oldTime The time at the beginning of the integration step on
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_dt The time step on the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If this flag is true, the initial estimate of phi will
  //!                  be set to zero. Otherwise, \a a_phiOld will be used.
  void computeDiffusion(LevelData<EBCellFAB>& a_DiffusiveTerm,
                        LevelData<EBCellFAB>& a_phiOld,
                        LevelData<EBCellFAB>& a_src,
                        EBFluxRegister* a_FineFluxRegPtr,
                        EBFluxRegister* a_crseFluxRegPtr,
                        const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                        const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                        Real a_oldTime,
                        Real a_crseOldTime,
                        Real a_crseNewTime,
                        Real a_dt,
                        int a_level,
                        bool a_zeroPhi = true,
                        bool a_rhsAlreadyKappaWeighted = false
                        );

  //! Integrates the diffusion equation, storing the result in \a a_phiNew.
  //! \param a_phiNew The new value of phi at the current grid level.
  //! \param a_phiOld The value of phi at the beginning of the time step,
  //!                 at the current grid level.
  //! \param a_src The source term in the diffusion equation at the current
  //!              grid level.
  //! \param a_fineFluxRegPtr A pointer to the flux register at the finer
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crseFluxRegPtr A pointer to the flux register at the coarser
  //!                         adjacent grid level (or NULL if there is none).
  //! \param a_crsePhiOldPtr A pointer to the LevelData object holding the
  //!                        old value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_crsePhiNewPtr A pointer to the LevelData object holding the
  //!                        updated value of phi on the coarser adjacent grid
  //!                        level (or NULL if there is none).
  //! \param a_oldTime The time at the beginning of the integration step on
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      on the coarser adjacent grid level.
  //! \param a_dt The time step on the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If this flag is true, the initial estimate of phi will
  //!                  be set to zero. Otherwise, \a a_phiOld will be used.
  void updateSoln(LevelData<EBCellFAB>& a_phiNew,
                  LevelData<EBCellFAB>& a_phiOld,
                  LevelData<EBCellFAB>& a_src,
                  EBFluxRegister* a_fineFluxRegPtr,
                  EBFluxRegister* a_crseFluxRegPtr,
                  const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                  const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                  Real a_oldTime,
                  Real a_crseOldTime,
                  Real a_crseNewTime,
                  Real a_dt,
                  int a_level,
                  bool a_zeroPhi = true,
                  bool a_rhsAlreadyKappaWeighted = false,
                  int a_fluxStartComponent = 0);


  void setEBLG(Vector<EBLevelGrid>& a_eblg)
  {
    m_isEBLGSet = true;
    m_eblg      = a_eblg;
  }

protected:
  bool m_isEBLGSet;
  Vector<EBLevelGrid> m_eblg;


};

#include "NamespaceFooter.H"
#endif


