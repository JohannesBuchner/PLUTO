#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASELEVELHEATSOLVER_H__
#define _BASELEVELHEATSOLVER_H__

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
#include "NamespaceHeader.H"

//! \class BaseLevelHeatSolver
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
class BaseLevelHeatSolver
{

  public:

  ///
  BaseLevelHeatSolver(const Vector<DisjointBoxLayout>&            a_grids,
                      const Vector<int>&                          a_refRat,
                      const ProblemDomain&                        a_level0Domain,
                      RefCountedPtr<AMRLevelOpFactory<LevelDataType> >&     a_opFact,
                      const RefCountedPtr<AMRMultiGrid<LevelDataType> >&     a_solver):
  m_grids(a_grids),
  m_refRat(a_refRat),
  m_level0Domain(a_level0Domain),
  m_ops(),
  m_solver(a_solver)
  {
    m_ops.resize(a_grids.size());
    Vector< AMRLevelOp<LevelDataType> * >& amrops =  m_solver->getAMROperators();

    for (int ilev = 0; ilev < m_ops.size(); ilev++)
      {
        m_ops[ilev] = dynamic_cast<LevelTGAHelmOp<LevelDataType,FluxDataType>* >(amrops[ilev]);
        if (m_ops[ilev]==NULL)
          {
            MayDay::Error("dynamic cast failed---is that operator really a TGAHelmOp?");
          }
      }
  }

  //! Destructor, called after destructors of BaseLevelHeatSolver subclasses.
  virtual ~BaseLevelHeatSolver()
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
  virtual void updateSoln(LevelDataType&           a_phiNew,
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
                          bool                     a_rhsAlreadyKappaWeighted= false,
                          int                      a_fluxStartComponent = 0) = 0;

  ///
  virtual
  void updateSoln(LevelDataType& a_phiNew,
                  LevelDataType& a_phiOld,
                  LevelDataType& a_src,
                  FluxRegisterType* a_fineFluxRegPtr,
                  FluxRegisterType* a_crseFluxRegPtr,
                  const LevelDataType* a_crsePhiOldPtr,
                  const LevelDataType* a_crsePhiNewPtr,
                  Real a_oldTime,
                  Real a_crseOldTime,
                  Real a_crseNewTime,
                  Real a_dt,
                  int a_level,
                  bool a_zeroPhi = true,
                  bool a_rhsAlreadyKappaWeighted = false,
                  int a_fluxStartComponent = 0)
  {
    MayDay::Error("need to overwrite this function with one that calls with your flux or call more general one");
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
  virtual void computeDiffusion(LevelDataType&           a_diffusiveTerm,
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
                                bool                     a_rhsAlreadyKappaWeighted = false
                                )
  {
    // The operator has no time-dependent parameters. Life is easier.

    // first compute updated solution
    LevelDataType phiNew;

    m_ops[a_level]->create(phiNew, a_phiOld);
    m_ops[a_level]->setToZero(phiNew);
    if (!a_zeroPhi)
      {
        m_ops[a_level]->assign(phiNew, a_phiOld);
      }

    updateSoln(phiNew, a_phiOld, a_src, a_flux,
               a_fineFluxRegPtr, a_crseFluxRegPtr,
               a_crsePhiOldPtr, a_crsePhiNewPtr,
               a_oldTime, a_crseOldTime,
               a_crseNewTime, a_dt, a_level, a_zeroPhi, a_rhsAlreadyKappaWeighted);

    // now subtract everything off to leave us with diffusive term
    m_ops[a_level]->incr(phiNew, a_phiOld, -1.0);
    m_ops[a_level]->scale(phiNew, 1.0/a_dt);

    //now multiply by a if there is an a
    m_ops[a_level]->diagonalScale(phiNew, false);

    // and finally, subtract off a_src
    m_ops[a_level]->incr(phiNew, a_src, -1.0);

    // what's left should be the time-centered diffusive part of the update
    m_ops[a_level]->assign(a_diffusiveTerm, phiNew);
  }

  ///
  /**
     calls set time and calls operator with given alpha and beta
   */
  virtual void applyOperator(LevelDataType&          a_ans,
                             const LevelDataType&    a_phi,
                             const LevelDataType*    a_phiC,
                             int                     a_level,
                             Real                    a_alpha,
                             Real                    a_beta,
                             bool                    a_applyBC)
  {
    m_ops[a_level]->setAlphaAndBeta(a_alpha, a_beta);

    LevelDataType zero;
    m_ops[a_level]->create(zero, a_ans);
    m_ops[a_level]->setToZero(zero);

    // set a_ans = helm(a_phi)
    //           = (I + factor*laplacian)(a_phi)
    if (a_applyBC)
      {
        if ( (a_phiC == NULL) || (a_level==0))
          {
            m_ops[a_level]->applyOp(a_ans, a_phi, false);
          }
        else
          {
            m_ops[a_level]->AMROperatorNF(a_ans, a_phi, *a_phiC, false);
          }
      }
    else
      {
        m_ops[a_level]->applyOpNoBoundary(a_ans, a_phi);
      }
  }

  //! Applies the Helmholtz operator to the solution \a a_phi at the given
  //! grid level. This will set \a a_ans to
  //! (I + \a a_mu * \a a_dt * L(\a a_phi).
  //! \param a_ans This will store the value of the Helmholtz operator applied
  //!              to \a a_phi.
  //! \param a_phi The value of the solution to which the Helmholtz operator
  //!              is applied.
  //! \param a_phiC A pointer to the value of the solution at the coarser
  //!               adjacent grid level, or NULL if there is no coarser grid
  //!               level.
  //! \param a_level The grid level at which the Helmholtz operator is applied.
  //! \param a_mu A number between 0 and 1 that defines the time within the
  //!             integration step at which the Helmholtz operator is evaluated.
  //! \param a_dt The size of the integration step at the current grid level.
  //! \param a_applyBC If set to true, inhomogeneous boundary conditions will
  //!                  be applied to the Helmholtz operator (both at the domain
  //!                  boundary and at the coarse-fine grid interfaces). If
  //!                  set to false, homogeneous boundary conditions will
  //!                  applied to the Helmholtz operator at these locations.
  //! protected because not flexible (alpha == 1)
  void applyHelm(LevelDataType&          a_ans,
                 const LevelDataType&    a_phi,
                 const LevelDataType*    a_phiC,
                 int                     a_level,
                 Real                    a_mu,
                 Real                    a_dt,
                 bool                    a_homogeneousBC)
  {
    // Set the operator's alpha and beta coefficients.
    Real factor  = a_mu*a_dt;
    m_ops[a_level]->setAlphaAndBeta(1.0, factor);

    LevelDataType zero;
    m_ops[a_level]->create(zero, a_ans);
    m_ops[a_level]->setToZero(zero);

    // set a_ans = helm(a_phi)
    //           = (I + factor*laplacian)(a_phi)
    if ( (a_phiC == NULL) || (a_level==0))
      {
        m_ops[a_level]->applyOp(a_ans, a_phi, a_homogeneousBC);
      }
    else
      {
        m_ops[a_level]->AMROperatorNF(a_ans, a_phi, *a_phiC, a_homogeneousBC);
      }
  }

  //! Adds flux contributions from the Helmholtz operator at the current
  //! grid level.
  //! \param a_diffusiveFlux The flux to which Helmholtz operator flux
  //!                        contributions will be accumulated.
  //! \param a_phi The solution at the current grid level to which the
  //!              Helmholtz operator is applied in order to evaluate the flux.
  //! \param a_level The grid level at which the flux is computed.
  //! \param a_mu A number between 0 and 1 that defines the time within the
  //!             integration step at which the Helmholtz operator is evaluated.
  //! \param a_dt The size of the integration step at the current grid level.
  //! \param a_sign A factor applied to the derivative term in the Helmholtz
  //!               equation. This allows the sign to be made positive or
  //!               negative as is necessary for the flux calculation.
  //! \param a_setToZero If true, \a a_diffusiveFlux will be set to zero before
  //!                    the fluxes are accumulated to it. Otherwise, fluxes
  //!                    will be added to the existing value.
  void incrementFlux(LevelData<FluxDataType>& a_diffusiveFlux,
                     LevelDataType&           a_phi,
                     int                      a_level,
                     Real                     a_mu,
                     Real                     a_dt,
                     Real                     a_sign,
                     bool                     a_setToZero)
  {
    Real factor  = a_sign*a_dt*a_mu;
    m_ops[a_level]->setAlphaAndBeta(1.0, factor);

    // increment flux
    m_ops[a_level]->fillGrad(a_phi);
    for (DataIterator  dit = a_phi.dataIterator(); dit.ok(); ++dit)
      {
        FluxDataType& thisFlux = a_diffusiveFlux[dit];
        FluxDataType tempFlux;
        tempFlux.define(thisFlux);

        tempFlux.setVal(0.0);
        if (a_setToZero)
          {
            thisFlux.setVal(0.0);
          }

        m_ops[a_level]->getFlux(tempFlux, a_phi, m_grids[a_level][dit], dit(), 1.0);
        thisFlux += tempFlux;
      }
  }

  //! Solves the Helmholtz equation
  //! (I - \a a_mu * \a a_dt * L(\a a_phi) = \a a_rhs
  //! for \a a_phi. Here it is assumed that a solution \a a_phiC exists
  //! on a coarser grid level.
  //! \param a_phi This will store the solution to the Helmholtz equation on
  //!              the given grid level.
  //! \param a_phiC The value of the solution given at the coarser adjacent
  //!               grid level. If there is no coarser grid level, this
  //!               parameter is ignored.
  //! \param a_rhs The right hand side of the Helmholtz equation at the given
  //!              grid level.
  //! \param a_level The grid level at which the Helmholtz equation is solved.
  //! \param a_mu A number between 0 and 1 that defines the time within the
  //!             integration step at which the Helmholtz operator is evaluated.
  //! \param a_dt The size of the integration step at the current grid level.
  //! \param a_zeroPhi If set to true, \a a_phi will be set to zero before
  //!                  the (iterative) solve takes place. Otherwise, a_phi is
  //!                  assumed to be an initial estimate for the solution.
  void solveHelm(LevelDataType& a_phi,
                 LevelDataType& a_phiC,
                 LevelDataType& a_rhs,
                 int            a_level,
                 Real           a_mu,
                 Real           a_dt,
                 bool           a_zeroPhi = true)
  {
    if (a_zeroPhi)
      {
        m_ops[a_level]->setToZero(a_phi);
      }
    Vector<LevelDataType* > phi(m_grids.size(), NULL);
    Vector<LevelDataType* > rhs(m_grids.size(), NULL);
    phi[a_level] = &a_phi;
    rhs[a_level] = &a_rhs;
    if (a_level > 0)
      {
        phi[a_level-1] = &a_phiC;
      }

    Real factor  = -a_dt*a_mu;
    resetSolverAlphaAndBeta(1.0, factor);

    m_solver->solve(phi, rhs, a_level, a_level, a_zeroPhi);
    int solverExitStatus = m_solver->m_exitStatus;
    if (solverExitStatus==2 || solverExitStatus==4 || solverExitStatus==6)
      {
        // These status codes correspond to the cases in which
        // norm is neither small enough nor reduced enough.
        // Either we've reached the maximum number of iterations,
        // or we've hung.
        pout() << "BaseLevelTGA:: WARNING: solver exitStatus == "
               << solverExitStatus << std::endl;
      }
  }

  //! Sets the alpha and beta parameters in each Helmholtz operator to the
  //! given values.
  //! \param a_alpha The identity term coefficient in the Helmholtz operator.
  //! \param a_beta The coefficient in front of the discrete derivative term
  //!               in the Helmholtz operator.
  void resetSolverAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta)
  {
    Vector<MGLevelOp<LevelDataType>* > ops = m_solver->getAllOperators();
    for (int iop = 0; iop < ops.size(); iop++)
      {
        LevelTGAHelmOp<LevelDataType, FluxDataType>* helmop =
          dynamic_cast<LevelTGAHelmOp<LevelDataType, FluxDataType>*>(ops[iop]);
        helmop->setAlphaAndBeta(a_alpha, a_beta);
      }
    for (int iop = 0; iop < m_ops.size(); iop++)
      {
        m_ops[iop]->setAlphaAndBeta(a_alpha, a_beta);
      }
  }

  //! Creates a new Helmholtz operator for use by the TGA integrator.
  //! \param a_indexSpace The ProblemDomain on which the new operator is
  //!                     defined.
  //! \param a_opFact A factory typename LevelDataTypehat generates new Helmholtz operators.
  //! \returns A pointer to a newly-allocated LevelTGAHelmOp instance.
  virtual LevelTGAHelmOp<LevelDataType, FluxDataType>*
  newOp(const ProblemDomain& a_indexSpace,
        RefCountedPtr<AMRLevelOpFactory<LevelDataType> >&   a_opFact)
  {
    LevelTGAHelmOp<LevelDataType, FluxDataType>* retval =
      dynamic_cast<LevelTGAHelmOp<LevelDataType, FluxDataType>*>(a_opFact->AMRnewOp(a_indexSpace));

    return retval;
  }

  //! Returns the number of grid levels on which this integrator operates.
  int size() const
  {
    return m_grids.size();
  }

protected:

  //! Interpolates a given quantity linearly in time using its beginning- and
  //! end-of-step values and placing the result in \a a_data.
  //! \param a_oldData The value of the quantity at the beginning of the time
  //!                  step at the grid level identified by \a a_level.
  //! \param a_newData The value of the quantity at the end of the time
  //!                  step at the grid level identified by \a a_level.
  //! \param a_time The time at which the quantity is to be interpolated
  //!               within the integration step.
  //! \param a_oldTime The beginning of the integration step at the grid level
  //!                  identified by \a a_level.
  //! \param a_newTime The end of the integration step at the grid level
  //!                  identified by \a a_level.
  //! \param a_level An index identifying the grid level at which the
  //!                interpolation takes place.
  void timeInterp(LevelDataType&       a_data,
                  const LevelDataType& a_oldData,
                  const LevelDataType& a_newData,
                  Real                 a_time,
                  Real                 a_oldTime,
                  Real                 a_newTime,
                  int                  a_level)
  {
    Real eps = 1.0e-10;
    CH_assert(a_newTime >= a_oldTime);
    Real diff = (a_newTime - a_oldTime);
    this->m_ops[a_level]->setToZero(a_data);
    if (diff < eps)
      {
        //no real time advance and don't want to divide by zero
        this->m_ops[a_level]->incr(a_data, a_oldData, 1.0);
      }
    else
      {
        Real factor = (a_time-a_oldTime)/(a_newTime - a_oldTime);
        this->m_ops[a_level]->incr(a_data, a_oldData, 1.0-factor);
        this->m_ops[a_level]->incr(a_data, a_newData, factor);
      }
  }

  //! The disjoint box layouts at every AMR grid level.
  Vector<DisjointBoxLayout>                             m_grids ;

  //! The refinement ratios between AMR grid levels.
  Vector<int>                                           m_refRat;

  //! The coarsest domain on which the Helmholtz equation is integrated.
  ProblemDomain                                         m_level0Domain;

  //! An array of the solver's Helmholtz operators at each grid level,
  //! casted to LevelTGAHelmOp instances. These are owned by the solver,
  //! so we shouldn't delete these.
  Vector< LevelTGAHelmOp<LevelDataType, FluxDataType>*  >  m_ops;

  //! The multigrid solver used to solve the Helmholtz equation.
  RefCountedPtr<AMRMultiGrid<LevelDataType> >       m_solver;
};

#include "NamespaceFooter.H"
#endif
