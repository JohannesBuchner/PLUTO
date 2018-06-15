#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASELEVELTGA_H__
#define _BASELEVELTGA_H__

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

//! \class BaseLevelTGA
//! This base class implements the 2nd-order implicit L0-stable time
//! integration algorithm developed by Twizell, Gumel, and Arigu for
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
class BaseLevelTGA : public BaseLevelHeatSolver<LevelDataType, FluxDataType, FluxRegisterType>
{

  public:

  //! Initializes the base class of a TGA time integrator. This must be called
  //! by any subclass of BaseLevelTGA.
  //! \param a_grids The DisjointBoxLayout on which the TGA scheme is to operate.
  //! \param a_refRatio An array containing the refinement ratios between the
  //!                   hierarchical AMR grid levels for the domain.
  //! \param a_level0Domain The domain at the coarsest AMR grid level.
  //! \param a_opFact A factory typename LevelDataTypehat is used to generate Helmholtz
  //!                 operators to be used by the scheme.
  //! \param a_solver An AMR Multigrid solver for solving the linear systems
  //!                 at each stage of the TGA integration scheme.
  BaseLevelTGA(const Vector<DisjointBoxLayout>&            a_grids,
               const Vector<int>&                          a_refRat,
               const ProblemDomain&                        a_level0Domain,
               RefCountedPtr<AMRLevelOpFactory<LevelDataType> >&     a_opFact,
               const RefCountedPtr<AMRMultiGrid<LevelDataType> >&     a_solver)
    :BaseLevelHeatSolver<LevelDataType,FluxDataType,FluxRegisterType>(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver),
     m_mu1(0.0),
     m_mu2(0.0),
     m_mu3(0.0),
     m_mu4(0.0),
     m_r1(0.0)
  {
    Real tgaEpsilon = 1.e-12;
#ifdef CH_USE_FLOAT
    tgaEpsilon = sqrt(tgaEpsilon);
#endif
    Real a = 2.0 - sqrt(2.0) - tgaEpsilon;
    m_mu1 = (a - sqrt( a*a - 4.0*a + 2.0))/2.0 ;
    m_mu2 = (a + sqrt( a*a - 4.0*a + 2.0))/2.0 ;
    m_mu3 = (1.0 - a);
    m_mu4 = 0.5 - a;

    Real discr = sqrt(a*a - 4.0*a + 2.0);
    m_r1 = (2.0*a - 1.0)/(a + discr);
  }

  //! Destructor, called after destructors of BaseLevelTGA subclasses.
  virtual ~BaseLevelTGA()
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
    // If our operators are time-independent, do the "easy" thing.
    if (!this->m_ops[a_level]->isTimeDependent())
      {
        this->updateSolnWithTimeIndependentOp(a_phiNew, a_phiOld, a_src, a_flux,
                                              a_fineFluxRegPtr, a_crseFluxRegPtr,
                                              a_crsePhiOldPtr, a_crsePhiNewPtr,
                                              a_oldTime, a_crseOldTime, a_crseNewTime,
                                              a_dt, a_level, a_zeroPhi,
                                              a_fluxStartComponent);
        return;
      }
    else
      {
        // We have to assume that the operator and its coefficients are
        // time-dependent.
        this->updateSolnWithTimeDependentOp(a_phiNew, a_phiOld, a_src, a_flux,
                                            a_fineFluxRegPtr, a_crseFluxRegPtr,
                                            a_crsePhiOldPtr, a_crsePhiNewPtr,
                                            a_oldTime, a_crseOldTime, a_crseNewTime,
                                            a_dt, a_level, a_zeroPhi,
                                            a_fluxStartComponent);
        return;
      }
  }

  //! Computes the time-centered diffusion term L(phi). This can be used to
  //! find contributions to the solution from diffusion. The diffusion term
  //! is computed by computing a finite difference approximation for d phi/dt
  //! using the updated and original values of phi and the time step. Most of
  //! the arguments given here are passed along to updateSoln and retain their
  //! significance therein.
  //! \param a_diffusiveTerm The diffusion term L(phi) will be stored here.
  //! \param a_phiOld The old solution (the value of phi at time n).
  //! \param a_src The source term on the right hand side of the Helmholtz
  //!              equation (used to fine the new value of phi).
  //! \param a_fineFluxRegPtr A pointer to the flux register representing the
  //!                         finer grid level adjacent to this one, or NULL
  //!                         if there is no finer grid level.
  //! \param a_crseFluxRegPtr A pointer to the flux register representing the
  //!                         coarser grid level adjacent to this one, or NULL
  //!                         if there is no coarser grid level.
  //! \param a_crsePhiOldTime A pointer to the value of phi at the beginning
  //!                         of the integration step at the coarser grid
  //!                         level adjacent to this one, or NULL if there
  //!                         is no coarser grid level.
  //! \param a_crsePhiNewTime A pointer to the value of phi at the end
  //!                         of the integration step at the coarser grid
  //!                         level adjacent to this one, or NULL if there
  //!                         is no coarser grid level.
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
  //! \param a_zeroPhi If set to true, the new value of phi will be set to
  //!                  zero before the integration takes place. Otherwise, it
  //!                  will be set to the value in \a a_diffusiveTerm.
  void computeDiffusion(LevelDataType&           a_diffusiveTerm,
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
                        bool                     a_rhsAlreadyKappaWeighted = false)
  {
    // in this function, we call the updateSoln function, compute
    // the update, then subtract off the extra pieces to return the
    // diffusive part of the update

    if (!this->m_ops[a_level]->isTimeDependent())
      {
        // The operator has no time-dependent parameters. Life is easier.

        // first compute updated solution
        LevelDataType phiNew;

        this->m_ops[a_level]->create(phiNew, a_phiOld);
        this->m_ops[a_level]->setToZero(phiNew);
        if (!a_zeroPhi)
          {
            this->m_ops[a_level]->assign(phiNew, a_phiOld);
          }

        updateSoln(phiNew, a_phiOld, a_src, a_flux,
                   a_fineFluxRegPtr, a_crseFluxRegPtr,
                   a_crsePhiOldPtr, a_crsePhiNewPtr,
                   a_oldTime, a_crseOldTime,
                   a_crseNewTime, a_dt, a_level, a_zeroPhi, a_rhsAlreadyKappaWeighted);

        // now subtract everything off to leave us with diffusive term
        this->m_ops[a_level]->incr(phiNew, a_phiOld, -1.0);
        this->m_ops[a_level]->scale(phiNew, 1.0/a_dt);

        //now multiply by a if there is an a
        this->m_ops[a_level]->diagonalScale(phiNew, false);

        // and finally, subtract off a_src
        this->m_ops[a_level]->incr(phiNew, a_src, -1.0);

        // what's left should be the time-centered diffusive part of the update
        this->m_ops[a_level]->assign(a_diffusiveTerm, phiNew);
      }
    else
      {
        // The operator has time-dependent coefficients. We must be more careful!

        //                                                      n+1    n
        //                                                n   (a    - a )
        // There's an extra source term (phi (da/dt) = phi  * -----------
        //                                                       dt
        // from the time-changing density that we need to subtract from the RHS.
        LevelDataType phidadt, aOld, aNew, rhs;
        this->m_ops[a_level]->create(phidadt, a_phiOld);
        this->m_ops[a_level]->create(aOld, a_phiOld);
        this->m_ops[a_level]->create(aNew, a_phiOld);
        this->m_ops[a_level]->create(rhs, a_phiOld);
        for (DataIterator dit = a_phiOld.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
          {
            // Set the old and new a coefficients.
            this->m_ops[a_level]->setTime(a_oldTime, 0.0, a_dt);
            aOld[dit()].setVal(1.);
            this->m_ops[a_level]->diagonalScale(aOld);
            this->m_ops[a_level]->setTime(a_oldTime, 1.0, a_dt);
            aNew[dit()].setVal(1.);
            this->m_ops[a_level]->diagonalScale(aNew);

            // Compute the above term.
            phidadt[dit()].axby(aNew[dit()], aOld[dit()], 1.0/a_dt, -1.0/a_dt);
            phidadt[dit()] *= a_phiOld[dit()];

            // Make a new right hand side out of the difference between the
            // source term and phidadt.
            rhs[dit()].axby(a_src[dit()], phidadt[dit()], 1.0, -1.0);
          }

        // Compute the updated solution.
        LevelDataType phiNew;
        this->m_ops[a_level]->create(phiNew, a_phiOld);
        this->m_ops[a_level]->setToZero(phiNew);
        if (!a_zeroPhi)
          {
            this->m_ops[a_level]->assign(phiNew, a_phiOld);
          }
        updateSoln(phiNew, a_phiOld, rhs, a_flux,
                   a_fineFluxRegPtr, a_crseFluxRegPtr,
                   a_crsePhiOldPtr, a_crsePhiNewPtr,
                   a_oldTime, a_crseOldTime,
                   a_crseNewTime, a_dt, a_level, a_zeroPhi);

        //         n+1         n
        // ([a phi]   - [a phi] )
        // ---------------------- - a_src  -> a_diffusiveTerm.
        //          dt
        for (DataIterator dit = a_phiOld.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
          {
            aNew[dit()] *= phiNew[dit()];
            aOld[dit()] *= a_phiOld[dit()];
            a_diffusiveTerm[dit()].axby(aNew[dit()], aOld[dit()], 1.0/a_dt, -1.0/a_dt);
            a_diffusiveTerm[dit()] -= a_src[dit()];
          }
      }
  }

protected:
  //! Sets the value of the source term in the Helmholtz equation on ghost
  //! cells. This method should be overridden by subclasses of BaseLevelTGA.
  //! \param a_src The source term in the Helmholtz equation whose ghost cell
  //!              values are to be set.
  //! \param a_dbl The disjoint box layout that indexes \a a_src.
  //! \param a_level The grid level at which the ghost cell values are to be
  //!                set.
  virtual void setSourceGhostCells(LevelDataType&        a_src,
                                   const DisjointBoxLayout& a_dbl,
                                   int a_level) = 0;

  //! The times within the integration step at which the operators are
  //! evaluated.
  Real m_mu1, m_mu2, m_mu3, m_mu4, m_r1;

private:

  // Update the solution assuming that the operator's coefficients are
  // independent of time. Same arguments as updateSoln.
  void updateSolnWithTimeIndependentOp(LevelDataType&           a_phiNew,
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

    LevelDataType rhst, srct, phis;
    this->m_ops[a_level]->create(rhst, a_src);
    this->m_ops[a_level]->create(srct, a_phiNew);
    this->m_ops[a_level]->create(phis, a_phiNew);

    this->m_ops[a_level]->setToZero(srct);
    this->m_ops[a_level]->setToZero(rhst);
    //this makes srct = a_src*dt
    this->m_ops[a_level]->incr(srct, a_src, a_dt);

    //save input guess if we are not using zero
    if (!a_zeroPhi)
      {
        this->m_ops[a_level]->setToZero(phis);
        this->m_ops[a_level]->incr(phis, a_phiNew, 1.0);
      }

    // Divide the source S by the identity coefficient a. srct = rhs*dt/a
    this->m_ops[a_level]->divideByIdentityCoef(srct);

    LevelDataType coarseData;

    if ((a_crsePhiOldPtr != NULL) && (a_level > 0))
      {
        this->m_ops[a_level-1]->create(coarseData, *a_crsePhiOldPtr);
        setSourceGhostCells(srct, this->m_grids[a_level], a_level);
      }

    //from here on k is kappa and L is kappa L
    //this makes rhs hold       (k*a I + mu4 k L) (S/a) dt
    this->applyHelm(rhst, srct, NULL, a_level, m_mu4, a_dt, true); //true means the application is homogeneous
    this->incrementFlux(a_flux, srct, a_level, m_mu4, a_dt, -1.0, true);

    // first need to compute coarse-level BC at this level's old time
    if (a_level > 0)
      {
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         a_oldTime, a_crseOldTime, a_crseNewTime, a_level-1);
      }

    //this makes a_phiNew hold (k*a I + mu3 k L) phi^n
    //'true' apply CF and domain BC
    this->applyHelm(a_phiNew,   a_phiOld, &coarseData, a_level, m_mu3, a_dt, true);
    this->incrementFlux(a_flux, a_phiOld,              a_level, m_mu3, a_dt, -1., false);

    //this makes rhs hold (k*a I + mu3 L) phi^n + dt(k*a I +mu4 L) S/a
    this->m_ops[a_level]->incr(rhst, a_phiNew, 1.0);

    // now construct coarse-level BC for intermediate solve
    // intermediate solution will be at time = oldTime + (1-r1)dt
    if (a_level > 0)
      {
        Real intermediateTime = a_oldTime + (1.0-m_r1)*a_dt;

        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         intermediateTime, a_crseOldTime, a_crseNewTime, a_level-1);
      }

    //if user thought her original guess was good, use it.
    if (!a_zeroPhi)
      {
        this->m_ops[a_level]->setToZero(a_phiNew);
        this->m_ops[a_level]->incr(a_phiNew, phis,  1.0);
      }

    //by now rhs =  ((k*a I + mu3 L) phi^n + dt(k*a I +mu4 L) S/a )
    //this makes phinew = (k*a I - mu2 k L)^-1 (rhs)
    this->solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu2, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew,       a_level, m_mu2, a_dt, -1.0, false);

    //this puts the answer into rhst so we can do the final solve
    this->m_ops[a_level]->assign(rhst, a_phiNew);

    // now construct CF-BC for final solve
    if (a_level > 0)
      {
        Real newTime = a_oldTime + a_dt;
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         newTime, a_crseOldTime, a_crseNewTime, a_level-1);
      }

    //by now rhs =  ((k*a I + mu3 L) phi^n + dt(k*a I +mu4 L) S/a )
    //this makes rhs hold k*a[(k*a I - mu2 L)^-1 (rhs)]
    this->m_ops[a_level]->diagonalScale(rhst, true);

    //if user thought her original guess was good, use it again.
    if (!a_zeroPhi)
      {
        this->m_ops[a_level]->setToZero(a_phiNew);
        this->m_ops[a_level]->incr(a_phiNew, phis,  1.0);
      }

    //this makes phinew = (k*a I - mu1 L)^-1 [ka ((k*a I - mu2 L)^-1 (orig rhs))]
    this->solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu1, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew, a_level, m_mu1, a_dt, -1.0, false);

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

  // Update the solution assuming that the operator's coefficients change
  // with time. Same arguments as updateSoln.
  void updateSolnWithTimeDependentOp(LevelDataType&           a_phiNew,
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
                                     int                      a_fluxStartComponent = 0)
  {
    int ncomp = a_phiNew.nComp();
    Interval intervBase(0, ncomp-1);
    Interval intervFlux(a_fluxStartComponent, a_fluxStartComponent + ncomp-1);

    CH_assert(a_level >= 0);
    CH_assert(a_level <  this->m_grids.size());
    CH_assert((a_level == 0) || (a_crsePhiOldPtr != NULL));
    CH_assert((a_level == 0) || (a_crsePhiNewPtr != NULL));
    CH_assert(a_crseNewTime >= a_crseOldTime);
    CH_assert(a_dt >= 0.);

    LevelDataType rhst, srct, phis;
    this->m_ops[a_level]->create(rhst, a_src);
    this->m_ops[a_level]->create(srct, a_phiNew);
    this->m_ops[a_level]->create(phis, a_phiNew);

    this->m_ops[a_level]->setToZero(srct);
    this->m_ops[a_level]->setToZero(rhst);
    this->m_ops[a_level]->incr(srct, a_src, 1.0);

    //save input guess if we are not using zero
    if (!a_zeroPhi)
      {
        this->m_ops[a_level]->setToZero(phis);
        this->m_ops[a_level]->incr(phis, a_phiNew, 1.0);
      }

    // In the comments, we use superscripts to denote different time centerings
    // in the operator L, the identity coefficient a, and the solution phi.
    // A '0' superscript means the quantity at time n.
    // A '1/2' superscript means the quantity at time n + 1/2.
    // A '*' superscript means the quantity at the "intermediate" time used
    //       by the TGA scheme.
    // A '1' superscript means the quantity at time n + 1.
    //
    // In embedded boundary problems we use k for kappa and L means kappa L.

    //-----------------------------------------------------
    //                                         1/2
    //          dt       1/2          1/2     S
    // Compute ----  x (a    I + mu4 L   ) x ---    -> rhst
    //           1/2                           1/2
    //          a                             a
    //-----------------------------------------------------

    // First, feed the half-step time to the operator to set its coefficients.
    this->m_ops[a_level]->setTime(a_oldTime, 0.5, a_dt);

    //                                                  1/2
    // Divide the source S by the identity coefficient a   .
    this->m_ops[a_level]->divideByIdentityCoef(srct);

    // Set the ghost cells in the source at the coarser level.
    LevelDataType coarseData;
    if ((a_crsePhiOldPtr != NULL) && (a_level > 0))
      {
        this->m_ops[a_level-1]->create(coarseData, *a_crsePhiOldPtr);
        setSourceGhostCells(srct, this->m_grids[a_level], a_level);
      }

    //     1/2                    1/2
    // (k a    I + mu4 k L) (S / a   ) -> rhst
    // (The false parameter tells the operator to use the source's ghost data
    //  instead of applying boundary conditions).
    //homogeneous BC for L(rhs)
    this->applyHelm(rhst, srct, NULL, a_level, m_mu4, a_dt, true);
    this->incrementFlux(a_flux, srct, a_level, a_dt*m_mu4, a_dt, -1.0, true);

    //--------------------------------------------------------
    //          1        0           0       0
    // Compute --- x (k a I + mu4 k L ) x phi  and add it to rhst.
    //           0
    //          a
    //--------------------------------------------------------

    // First we compute the coarse-level boundary condition data at time n.
    if (a_level > 0)
      {
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         a_oldTime, a_crseOldTime, a_crseNewTime, a_level-1);
      }

    // Feed the old time to the operator.
    this->m_ops[a_level]->setTime(a_oldTime, 0.0, a_dt);

    //             0         0     0
    // Compute (k a I + mu3 L ) phi and store it in a_phiNew for now.
    // (The true parameter means apply the coarse-fine and domain boundary
    //  conditions).
    //inhommogeneous bc for solution
    this->applyHelm(a_phiNew, a_phiOld, &coarseData, a_level, m_mu3, a_dt, false);
    this->incrementFlux(a_flux, a_phiOld, a_level, m_mu3, a_dt, -1., false);

    //                     n
    // Divide a_phiNew by a .
    this->m_ops[a_level]->divideByIdentityCoef(a_phiNew);

    // Add a_phiNew to rhst to obtain the right hand side of our
    // TGA update equation.
    this->m_ops[a_level]->incr(rhst, a_phiNew, 1.0);

    //-----------------------------------------------------------------
    //           *           *     n    *
    // Solve (k a I + mu4 k L ) phi  = a rhst and place the answer in rhst,
    //-----------------------------------------------------------------

    // Construct the coarse-level boundary condition data for our TGA
    //                                               *   n
    // intermediate solve. The intermediate time is t = t + (1-r1)dt.
    if (a_level > 0)
      {
        Real intermediateTime = a_oldTime + (1.0-m_r1)*a_dt;

        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         intermediateTime, a_crseOldTime, a_crseNewTime, a_level-1);
      }

    // If the user thought her original guess was good, use it.
    if (!a_zeroPhi)
      {
        this->m_ops[a_level]->setToZero(a_phiNew);
        this->m_ops[a_level]->incr(a_phiNew, phis,  1.0);
      }

    // Feed the intermediate time to the operator to get the time-centered
    // coefficients.
    this->m_ops[a_level]->setTime(a_oldTime, 1.0 - m_r1, a_dt);

    //                   *
    // Multiply rhst by a .
    this->m_ops[a_level]->diagonalScale(rhst, false);

    //                   *         * -1
    // This computes (k a I - mu2 L )   (rhst) and stores it in a_phiNew.
    this->solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu2, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew, a_level, m_mu2, a_dt, -1.0, false);

    // This puts the answer into rhst.
    this->m_ops[a_level]->assign(rhst, a_phiNew);

    //--------------------------------------------------
    //           1           1     1                      1
    // Solve (k a I + k mu4 L ) phi  = rhst and scale by a .
    //--------------------------------------------------

    // Compute the coarse-fine boundary condition data for the final solve.
    if (a_level > 0)
      {
        Real newTime = a_oldTime + a_dt;
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         newTime, a_crseOldTime, a_crseNewTime, a_level-1);
      }

    // Feed the new time to the operator.
    this->m_ops[a_level]->setTime(a_oldTime, 1.0, a_dt);

    //                  1
    // Scale rhst by k a .
    this->m_ops[a_level]->diagonalScale(rhst);

    // If the user thought her original guess was good, use it again.
    if (!a_zeroPhi)
      {
        this->m_ops[a_level]->setToZero(a_phiNew);
        this->m_ops[a_level]->incr(a_phiNew, phis, 1.0);
      }

    //             1           1  -1    1     *           * -1  *
    // Compute [k a I - mu1 k L )]   k a [(k a I - mu2 k L )   a (orig rhs)].
    this->solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu1, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew, a_level, m_mu1, a_dt, -1.0, false);

    // Now increment the flux registers -- note that because of the way
    // we defined the fluxes, the dt multiplier is already in the
    // flux.
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
  BaseLevelTGA& operator=(const BaseLevelTGA&);
  BaseLevelTGA(const BaseLevelTGA& a_opin);
  BaseLevelTGA();
};

#include "NamespaceFooter.H"

#endif
