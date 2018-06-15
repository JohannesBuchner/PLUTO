#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASELEVELCRANKNICOLSON_H_
#define _BASELEVELCRANKNICOLSON_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include <REAL.H>
#include <Box.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <ProblemDomain.H>
#include "AMRTGA.H"
#include "BaseLevelHeatSolver.H"
#include "NamespaceHeader.H"

//! \class BaseLevelCrankNicolson
//! This base class implements the 1st-order implicit L0-stable time
//! integration algorithm for
//! solving elliptic equations. It relies upon linear algebraic operations
//! defined in the underlying Helmholtz operators.
//! \tparam LevelDataType The type used to store data at a grid level.
//!                       This is usually LevelData<T>, where T is some
//!                       cell-centered FArrayBox type.
//! \tparam FluxDataType The type used to store flux data at a grid
//!                      level. This is usually an array box clas that stores
//!                      fluxes.
//! \tparam FluxRegisterType The type used to store flux register data for
//!                          interactions between different grid levels.
//!                          This is usually a flux register class.
template <class LevelDataType,
          class FluxDataType,
          class FluxRegisterType>
class BaseLevelCrankNicolson : public BaseLevelHeatSolver<LevelDataType, FluxDataType, FluxRegisterType>
{

  public:

  //! Initializes the base class of a CrankNicolson time integrator. This must be called
  //! by any subclass of BaseLevelCrankNicolson.
  //! \param a_grids The DisjointBoxLayout on which the CrankNicolson scheme is to operate.
  //! \param a_refRatio An array containing the refinement ratios between the
  //!                   hierarchical AMR grid levels for the domain.
  //! \param a_level0Domain The domain at the coarsest AMR grid level.
  //! \param a_opFact A factory typename LevelDataTypehat is used to generate Helmholtz
  //!                 operators to be used by the scheme.
  //! \param a_solver An AMR Multigrid solver for solving the linear systems
  //!                 at each stage of the CrankNicolson integration scheme.
  BaseLevelCrankNicolson(const Vector<DisjointBoxLayout>&            a_grids,
               const Vector<int>&                          a_refRat,
               const ProblemDomain&                        a_level0Domain,
               RefCountedPtr<AMRLevelOpFactory<LevelDataType> >&     a_opFact,
               const RefCountedPtr<AMRMultiGrid<LevelDataType> >&     a_solver)
    :BaseLevelHeatSolver<LevelDataType,FluxDataType,FluxRegisterType>(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver)
  {
  }

  //! Destructor, called after destructors of BaseLevelCrankNicolson subclasses.
  virtual ~BaseLevelCrankNicolson()
  {
  }

  //! Integrates the helmholtz equation represented by this object, placing
  //! the new solution in \a a_phiNew.
  //! \param a_phiNew The new solution (the value of phi at time n + 1) will
  //!                 be stored here.
  //! \param a_phiOld The old solution (the value of phi at time n).
  //! \param a_src The source term on the right hand side of the Helmholtz
  //!              equation.
  //! \param a_flux This will store the flux computed at the current grid
  //!               level during the solution of the Helmholtz equation.
  //! \param a_fineFluxRegPtr A pointer to the flux register representing the
  //!                         finer grid level adjacent to this one, or NULL
  //!                         if there is no finer grid level.
  //! \param a_crseFluxRegPtr A pointer to the flux register representing the
  //!                         coarser grid level adjacent to this one, or NULL
  //!                         if there is no coarser grid level.
  //! \param a_oldTime The time at the beginning of the integration step at
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_dt The size of the integration step at the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If set to true, \a a_phiNew will be set to zero before
  //!                  the integration takes place. Otherwise, a_phiNew is
  //!                  assumed to be an initial estimate for the solution in
  //!                  the iterative linear solve.
  //! \param a_fluxStartComponent An index identifying the component at which
  //!                             flux data begins within \a a_fineFluxRegPtr
  //!                             and \a a_crseFluxRegPtr.
  void updateSoln(LevelDataType&           a_phiNew,
                  LevelDataType&           a_phiOld,
                  LevelDataType&           a_src,
                  LevelData<FluxDataType>& a_flux,
                  FluxRegisterType*        a_fineFluxRegPtr,
                  FluxRegisterType*        a_crseFluxRegPtr,
                  const LevelDataType*     a_crsePhiOldPtr,
                  const LevelDataType*     a_crsePhiNewPtr,
                  Real                     a_oldTime,
                  Real                     a_crseOldTime,
                  Real                     a_crseNewTime,
                  Real                     a_dt,
                  int                      a_level,
                  bool                     a_zeroPhi = true,
                  bool                     a_rhsAlreadyKappaWeighted = false,
                  int                      a_fluxStartComponent = 0)
  {
    CH_assert(!this->m_ops[a_level]->isTimeDependent());
    int ncomp = a_phiNew.nComp();
    Interval intervBase(0, ncomp-1);
    Interval intervFlux(a_fluxStartComponent, a_fluxStartComponent + ncomp-1);

    CH_assert(a_level >= 0);
    CH_assert(a_level <  this->m_grids.size());
    CH_assert((a_level == 0) || (a_crsePhiOldPtr != NULL));
    CH_assert((a_level == 0) || (a_crsePhiNewPtr != NULL));
    CH_assert(a_crseNewTime >= a_crseOldTime);
    CH_assert(a_dt >= 0.);

    LevelDataType rhst, phit;
    this->m_ops[a_level]->create(rhst, a_src);
    this->m_ops[a_level]->create(phit, a_phiNew);

    this->m_ops[a_level]->setToZero(phit);
    this->m_ops[a_level]->setToZero(rhst);
    if (a_zeroPhi)
      {
        this->m_ops[a_level]->setToZero(a_phiNew);
      }

    LevelDataType coarseData;
    if ((a_crsePhiOldPtr != NULL) && (a_level > 0))
      {
        this->m_ops[a_level-1]->create(coarseData, *a_crsePhiOldPtr);
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                   a_oldTime, a_crseOldTime, a_crseNewTime, a_level-1);
      }

    //this makes phit = (k a I + dt/2 k L) phiOld
    //inhomogeneous BC for solution
    this->applyHelm(phit, a_phiOld, &coarseData, a_level, 0.5, a_dt, false);

    //set  rhst = dt*src
    this->m_ops[a_level]->incr(rhst, a_src   , a_dt);

    //multiply phi old by kappa*acoef
    this->m_ops[a_level]->diagonalScale(phit, true);

    //multiply rhs by kappa (but NOT by a) rhs = kappa dt rhs
    if (!a_rhsAlreadyKappaWeighted)
      this->m_ops[a_level]->kappaScale(rhst);

    //add both together to make rhs for crank solve = (kappa a I + kappa  dt/2 L) phiold  + kappa rhs dt
    this->m_ops[a_level]->incr(rhst, phit, 1.0);

    this->solveHelm(a_phiNew, coarseData, rhst, a_level, 0.5, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew,       a_level, 0.5, a_dt, -1.0, false);

    // now increment flux registers -- note that because of the way
    // we defined the fluxes, the dt multiplier is already in the
    // flux
    if ((a_fineFluxRegPtr != NULL) && (a_level < this->m_grids.size()-1))
      {
        Real fluxMult = 1.0;
        for (DataIterator dit = this->m_grids[a_level].dataIterator(); dit.ok(); ++dit)
          {
            FluxDataType& thisFlux = a_flux[dit];
            for (int dir=0; dir<SpaceDim; ++dir)
              {
                a_fineFluxRegPtr->incrementCoarse(thisFlux[dir],
                                                  fluxMult, dit(),
                                                  intervBase, // source
                                                  intervFlux, // dest
                                                  dir);
              }
          }
      } // end if there is a finer-level

    if ((a_crseFluxRegPtr != NULL) && (a_level > 0))
      {
        Real fluxMult = 1.0;

        for (DataIterator dit = this->m_grids[a_level].dataIterator(); dit.ok(); ++dit)
          {
            FluxDataType& thisFlux = a_flux[dit];
            for (int dir=0; dir<SpaceDim; ++dir)
              {
                a_crseFluxRegPtr->incrementFine(thisFlux[dir], fluxMult, dit(),
                                                intervBase, // source
                                                intervFlux, // dest
                                                dir);
              }
          }
      } // end if there is a coarser level
  }

private:
  // Disallowed operators.
  BaseLevelCrankNicolson& operator=(const BaseLevelCrankNicolson&);
  BaseLevelCrankNicolson(const BaseLevelCrankNicolson& a_opin);
  BaseLevelCrankNicolson();
};

#include "NamespaceFooter.H"
#endif
