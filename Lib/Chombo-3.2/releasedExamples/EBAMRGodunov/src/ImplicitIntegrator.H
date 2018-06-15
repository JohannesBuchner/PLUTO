#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IMPLICITINTEGRATOR_H__
#define _IMPLICITINTEGRATOR_H__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "SPACE.H"
#include "REAL.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRTGA.H"
#include "BaseImplicitIntegrationStrategy.H"

#include "NamespaceHeader.H"

//! \class BasicImplicitIntegrator
//! This class provides an interface for the implicit time integration of
//! various elliptic equations. An instance of this class uses strategy
//! objects to perform the actual integration. The motivation for this
//! approach is to provide a single interface that can be used to provide
//! various time integration schemes for an algorithm without changing its
//! structure.
//!
//! DO NOT CREATE SUBCLASSES OF THIS CLASS IF YOU ARE IMPLEMENTING A NEW TIME
//! INTEGRATION SCHEME. Implement an ImplicitIntegrationStrategy and use it to
//! create an ImplicitIntegrator. This class uses a traits class to decide
//! what types to use. For an analogous example of this kind of thing, check
//! out the basic_string class in the C++ Standard Template Library
//! (did you know that std::string is actually a typedef?)
//! \tparam ImplicitIntegratorTraits A traits class that contains type
//!                                  information for the integrator. This
//!                                  should be RegularImplicitIntegratorTraits
//!                                  or EBImplicitIntegratorTraits.
//! \tparam LevelDataType The type used to store data at a grid level.
//!                       This is usually LevelData<T>, where T is some
//!                       cell-centered FArrayBox type.
//! \tparam FluxDataType The type used to store flux data at a grid
//!                      level. This is usually an array box clas that stores
//!                      fluxes.
//! \tparam FluxRegisterType The type used to store flux register data for
//!                          interactions between different grid levels.
//!                          This is usually a flux register class.
template <typename ImplicitIntegratorTraits>
class BasicImplicitIntegrator
{
  public:

  //! The type of LevelData used by this Integrator.
  typedef typename ImplicitIntegratorTraits::LevelDataType LevelDataType;
  typedef typename ImplicitIntegratorTraits::FluxDataType FluxDataType;
  typedef typename ImplicitIntegratorTraits::FluxRegisterType FluxRegisterType;
  typedef typename ImplicitIntegratorTraits::GridPolicy GridPolicy;
  typedef typename ImplicitIntegratorTraits::GridType GridType;
  typedef typename ImplicitIntegratorTraits::OpFactoryType OpFactoryType;
  typedef typename ImplicitIntegratorTraits::SolverType SolverType;

  //! Constructs an implicit time integrator with the given strategy. This
  //! method shouldn't be called directly unless the caller knows what she is
  //! doing--use a factory method instead.
  //! \param updateStrategy A pointer to an ImplicitIntegrationStrategy
  //!                       subclass that implements the algorithm for
  //!                       implicitly computing a solution at time n + 1
  //!                       from that at time n.
  //! \param a_opFact A pointer to a factory that is used to generate the
  //!                 underlying Helmholtz operators that implement the time
  //!                 integration scheme.
  //! \param a_solver A pointer to an AMRMultiGrid solver that solves the
  //!                 underlying equations within the integration scheme.
  //! \param a_grids The box layout on which the domain is defined.
  BasicImplicitIntegrator(RefCountedPtr<BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits> >& a_updateStrategy,
                          RefCountedPtr<OpFactoryType>& a_opFact,
                          const RefCountedPtr<SolverType>& a_solver,
                          const Vector<GridType>& a_grids);

  //! Destructor.
  ~BasicImplicitIntegrator();

  //! Integrates the equation
  //! \f$a(\vec{x})\frac{\partial \phi}{\partial t} = L(\phi, t) + S(\vec{x}, t)$,
  //! producing a value for \f$\phi^{n+1}$ from \f$\phi^n$.
  //! \param a_phiNew The new solution \f$\phi^{n+1} will be stored here.
  //! \param a_phiOld The old solution \f$\phi^n$.
  //! \param a_src The source term $S$ on the right hand side of the Helmholtz
  //!              equation.
  //! \param a_fineFluxRegPtr A pointer to the flux register representing the
  //!                         finer grid level adjacent to this one, or NULL
  //!                         if there is no finer grid level.
  //! \param a_crseFluxRegPtr A pointer to the flux register representing the
  //!                         coarser grid level adjacent to this one, or NULL
  //!                         if there is no coarser grid level.
  //! \param a_crsePhiOldPtr A pointer to the coarser-level solution at the
  //!                        beginning of its time step, or NULL if no coarser
  //!                        level exists.
  //! \param a_crsePhiNewPtr A pointer to the coarser-level solution at the
  //!                        end of its time step, or NULL if no coarser
  //!                        level exists.
  //! \param a_oldTime The time \f$t^n$ at the beginning of the integration step at
  //!                  the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_dt The size \f$\Delta t$ of the integration step at the current
  //!             grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If set to true, \a a_phiNew will be set to zero before
  //!                  the integration takes place. Otherwise, a_phiNew is
  //!                  assumed to be an initial estimate for the solution in
  //!                  the iterative linear solve.
  //! \param a_fluxStartComponent An index identifying the component at which
  //!                             flux data begins within \a a_fineFluxRegPtr
  //!                             and \a a_crseFluxRegPtr.
  void updateSolution(LevelDataType& a_phiNew,
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
                      int a_fluxStartComponent = 0);

  //! Computes the time-centered diffusion term
  //! \f$L(\phi) = a(\vec{x})\frac{\partial \phi}{\partial t} - S$.
  //! \param a_diffusiveTerm The diffusion term \f$L(\phi)$ will be stored here.
  //! \param a_phiOld The solution \f$\phi^n$ at time \f$n$.
  //! \param a_src The source term \f$S$ on the right hand side of the
  //!              diffusion equation.
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
  //! \param a_oldTime The time \f$t^n$ at the beginning of the integration step
  //!                  at the current grid level.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_dt The size \f$\Delta t$ of the integration step at the current grid level.
  //! \param a_level The current grid level.
  //! \param a_zeroPhi If set to true, the new value of phi will be set to
  //!                  zero before the integration takes place. Otherwise, it
  //!                  will be set to the value in \a a_diffusiveTerm.
  void computeDiffusion(LevelDataType& a_diffusiveTerm,
                        LevelDataType& a_phiOld,
                        LevelDataType& a_src,
                        FluxRegisterType* a_FineFluxRegPtr,
                        FluxRegisterType* a_crseFluxRegPtr,
                        const LevelDataType* a_crsePhiOldPtr,
                        const LevelDataType* a_crsePhiNewPtr,
                        Real a_oldTime,
                        Real a_crseOldTime,
                        Real a_crseNewTime,
                        Real a_dt,
                        int a_level,
                        bool a_zeroPhi = true,
                        int a_fluxStartComponent = 0);

  //! Explicitly computes \f$L(\phi)$, applying the diffusion operator \f$L$
  //! to the solution \f$\phi$ at the given time \f$t$.
  //! \param a_Lphi The computed value of \f$L(\phi)$ is stored here.
  //! \param a_phi The solution $\phi$ to which the diffusion operator
  //!              is applied.
  //! \param a_crsePhiOldTime A pointer to the value of phi at the beginning
  //!                         of the integration step at the coarser grid
  //!                         level adjacent to this one, or NULL if there
  //!                         is no coarser grid level.
  //! \param a_crsePhiNewTime A pointer to the value of phi at the end
  //!                         of the integration step at the coarser grid
  //!                         level adjacent to this one, or NULL if there
  //!                         is no coarser grid level.
  //! \param a_time The time \f$t$ at which the boundary conditions are applied.
  //!               If \a a_applyBC is false, this argument is ignored.
  //! \param a_crseOldTime The time at the beginning of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_crseNewTime The time at the end of the integration step
  //!                      at the coarser adjacent grid level. This parameter
  //!                      is ignored if there is no coarser grid level.
  //! \param a_dt The size of the integration step at the current grid level.
  //! \param a_level The current grid level.
  //! \param a_applyBC If set to true, inhomogeneous boundary conditions will
  //!                  be applied to the diffusion operator (both at the domain
  //!                  boundary and at the coarse-fine grid interfaces). If
  //!                  set to false, homogeneous boundary conditions will
  //!                  applied to the operator at these locations.
  void applyDiffusionOperator(LevelDataType& a_Lphi,
                              LevelDataType& a_phi,
                              const LevelDataType* a_crsePhiOldPtr,
                              const LevelDataType* a_crsePhiNewPtr,
                              Real a_time,
                              Real a_crseOldTime,
                              Real a_crseNewTime,
                              int a_level,
                              bool a_applyBC);

  //! Returns the number of grid levels on which this integrator operates.
  int numLevels() const
  {
    return m_grids.size();
  }

  //! Returns a reference to the Helmholtz operator at the given
  //! AMR level.
  LevelTGAHelmOp<LevelDataType, FluxDataType>& HelmholtzOperator(int a_level)
  {
    CH_assert(a_level >= 0);
    CH_assert(a_level < m_grids.size());
    return *(m_ops[a_level]);
  }

  private:

  // Increments the coarser- and finer-level flux registers using the
  // flux data in a_flux.
  void
  incrementFluxRegisters(FluxRegisterType*        a_fineFluxRegPtr,
                         FluxRegisterType*        a_crseFluxRegPtr,
                         const LevelData<FluxDataType>& a_flux,
                         int a_level,
                         int a_fluxStartComponent,
                         int a_numSolutionComponents) const;

  //! The update strategy used to perform the implicit time integration.
  RefCountedPtr<BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits> > m_strategy;

  //! The operator factory used to construct Helmholtz operators.
  RefCountedPtr<OpFactoryType> m_opFact;

  //! The multigrid solver used to solve the Helmholtz equation.
  RefCountedPtr<SolverType> m_solver;

  //! The disjoint box layouts at every AMR level.
  Vector<GridType> m_grids ;

  //! An array of the solver's AMR Helmholtz operators at each AMR
  //! level, casted to LevelTGAHelmOp instances. These are owned by the
  //! solver, so we shouldn't delete these.
  Vector< LevelTGAHelmOp<LevelDataType, FluxDataType>* > m_ops;

  // Disallowed operators.
  BasicImplicitIntegrator();
  BasicImplicitIntegrator(const BasicImplicitIntegrator&);
  BasicImplicitIntegrator& operator=(const BasicImplicitIntegrator&);
};


//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
BasicImplicitIntegrator<ImplicitIntegratorTraits>::
BasicImplicitIntegrator(RefCountedPtr<BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits> >& a_updateStrategy,
                        RefCountedPtr<typename ImplicitIntegratorTraits::OpFactoryType>& a_opFact,
                        const RefCountedPtr<typename ImplicitIntegratorTraits::SolverType>& a_solver,
                        const Vector<typename ImplicitIntegratorTraits::GridType>& a_grids):
  m_strategy(a_updateStrategy),
  m_opFact(a_opFact),
  m_solver(a_solver),
  m_grids(a_grids),
  m_ops(a_grids.size())
{
  // Get our own list of operators for AMR refinement levels from the
  // multigrid solver.
  Vector< AMRLevelOp<LevelDataType> * >& amrops =  m_solver->getAMROperators();
  for (int ilev = 0; ilev < m_ops.size(); ilev++)
  {
    m_ops[ilev] = dynamic_cast<LevelTGAHelmOp<LevelDataType,FluxDataType>* >(amrops[ilev]);
    if (m_ops[ilev] == NULL)
    {
      MayDay::Error("ImplicitIntegrator: dynamic cast failed---is that operator really a TGAHelmOp?");
    }
  }

  // Set up the strategy.
  m_strategy->setGrids(m_grids);
  m_strategy->setMultiGridSolver(m_solver);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
BasicImplicitIntegrator<ImplicitIntegratorTraits>::
~BasicImplicitIntegrator()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BasicImplicitIntegrator<ImplicitIntegratorTraits>::
updateSolution(LevelDataType& a_phiNew,
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
               bool a_zeroPhi,
               int a_fluxStartComponent)
{
  // Check the inputs.
  CH_assert(a_phiOld.nComp() == a_phiNew.nComp());
  CH_assert(a_level >= 0);
  CH_assert(a_level <  m_grids.size());
  CH_assert((a_level == 0) || (a_crsePhiOldPtr != NULL));
  CH_assert((a_level == 0) || (a_crsePhiNewPtr != NULL));
  CH_assert((a_level > 0) && (a_crseNewTime >= a_crseOldTime) || (a_level == 0));
  CH_assert(a_dt >= 0.);

  // Set the underlying Helmholtz operator to be used by the strategy.
  m_strategy->setAMRLevel(a_level);
  m_strategy->setHelmholtzOp(m_ops[a_level]);

  // Initialize a flux for the strategy to keep track of flux contributions.
  m_strategy->initFlux(a_phiNew.nComp());

  // Zero the initial solution if it is desired, and have the strategy
  // adopt the same practice.
  m_strategy->zeroInitialSolution(a_phiNew, a_zeroPhi);

  // Set coarse data for boundary condition interpolation.
  m_strategy->setCoarseData(a_crsePhiOldPtr, a_crsePhiNewPtr,
                            a_crseOldTime, a_crseNewTime);

  // Call the strategy.
  m_strategy->updateSolution(a_phiNew, a_phiOld, a_src, a_oldTime, a_dt);

  // Increment the flux registers.
  incrementFluxRegisters(a_fineFluxRegPtr, a_crseFluxRegPtr, m_strategy->flux(),
                         a_level, a_fluxStartComponent, a_phiNew.ncomp());

  // That's it.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BasicImplicitIntegrator<ImplicitIntegratorTraits>::
computeDiffusion(LevelDataType& a_diffusiveTerm,
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
                 bool a_zeroPhi,
                 int a_fluxStartComponent)
{
  // Check the inputs.
  CH_assert(a_phiOld.nComp() == a_diffusiveTerm.nComp());
  CH_assert(a_level >= 0);
  CH_assert(a_level <  m_grids.size());
  CH_assert((a_level == 0) || (a_crsePhiOldPtr != NULL));
  CH_assert((a_level == 0) || (a_crsePhiNewPtr != NULL));
  CH_assert(((a_level > 0) && (a_crseNewTime >= a_crseOldTime)) || (a_level == 0));
  CH_assert(a_dt >= 0.);

  // Set the underlying Helmholtz operator to be used by the strategy.
  m_strategy->setAMRLevel(a_level);
  m_strategy->setHelmholtzOp(m_ops[a_level]);

  // Initialize a flux for the strategy to keep track of flux contributions.
  m_strategy->initFlux(a_phiOld.nComp());

  // Set coarse data for boundary condition interpolation.
  m_strategy->setCoarseData(a_crsePhiOldPtr, a_crsePhiNewPtr,
                            a_crseOldTime, a_crseNewTime);

  // Zero the initial solution if it is desired, and have the strategy
  // adopt the same practice.
  m_strategy->zeroInitialSolution(a_diffusiveTerm, a_zeroPhi);

  // Call the strategy.
  m_strategy->computeDiffusion(a_diffusiveTerm, a_phiOld, a_src,
                               a_oldTime, a_dt);

  // Increment the flux registers.
  incrementFluxRegisters(a_fineFluxRegPtr, a_crseFluxRegPtr, m_strategy->flux(),
                         a_level, a_fluxStartComponent, a_phiOld.nComp());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BasicImplicitIntegrator<ImplicitIntegratorTraits>::
incrementFluxRegisters(FluxRegisterType* a_fineFluxRegPtr,
                       FluxRegisterType* a_crseFluxRegPtr,
                       const LevelData<FluxDataType>& a_flux,
                       int a_level,
                       int a_fluxStartComponent,
                       int a_numSolutionComponents) const
{
  // Increment the coarser-level flux register.
  Interval intervSoln(0, a_numSolutionComponents-1);
  Interval intervFlux(a_fluxStartComponent,
                      a_fluxStartComponent + a_numSolutionComponents - 1);
  if ((a_fineFluxRegPtr != NULL) && (a_level < m_grids.size()))
  {
    Real fluxMult = 1.0;
    const DisjointBoxLayout& dbl = GridPolicy::disjointBoxLayout(m_grids[a_level]);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const FluxDataType& thisFlux = a_flux[dit];
      for (int dir = 0; dir< SpaceDim; ++dir)
      {
        a_fineFluxRegPtr->incrementCoarse(thisFlux[dir],
                                          fluxMult, dit(),
                                          intervSoln, // source
                                          intervFlux, // dest
                                          dir);
      }
    }
  }

  // Now increment the finer-level flux register.
  if ((a_crseFluxRegPtr != NULL) && (a_level > 0))
  {
    Real fluxMult = 1.0;
    const DisjointBoxLayout& dbl = GridPolicy::disjointBoxLayout(m_grids[a_level]);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const FluxDataType& thisFlux = a_flux[dit];
      for (int dir = 0; dir < SpaceDim; ++dir)
      {
        a_crseFluxRegPtr->incrementFine(thisFlux[dir], fluxMult, dit(),
                                        intervSoln, // source
                                        intervFlux, // dest
                                        dir);
      }
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BasicImplicitIntegrator<ImplicitIntegratorTraits>::
applyDiffusionOperator(LevelDataType& a_Lphi,
                       LevelDataType& a_phi,
                       const LevelDataType* a_crsePhiOldPtr,
                       const LevelDataType* a_crsePhiNewPtr,
                       Real a_time,
                       Real a_crseOldTime,
                       Real a_crseNewTime,
                       int a_level,
                       bool a_applyBC)
{
  // Check the inputs.
  CH_assert(a_phi.nComp() == a_Lphi.nComp());
  CH_assert(a_level >= 0);
  CH_assert(a_level <  m_grids.size());
  CH_assert((a_level == 0) || (a_crsePhiOldPtr != NULL));
  CH_assert((a_level == 0) || (a_crsePhiNewPtr != NULL));
  CH_assert(((a_level > 0) && (a_crseNewTime >= a_crseOldTime)) || (a_level == 0));

  // Set the underlying Helmholtz operator to be used by the strategy.
  m_strategy->setAMRLevel(a_level);
  m_strategy->setHelmholtzOp(m_ops[a_level]);

  // Initialize a flux for the strategy to keep track of flux contributions.
  m_strategy->initFlux(a_phi.nComp());

  // Set coarse data for boundary condition interpolation.
  m_strategy->setCoarseData(a_crsePhiOldPtr, a_crsePhiNewPtr,
                            a_crseOldTime, a_crseNewTime);

  // Call the strategy.
  m_strategy->applyDiffusionOperator(a_Lphi, a_phi, a_time, a_applyBC);

  // No need to increment flux registers here since this is an explicit
  // application of the operator. We're done.
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

// Include regular geometry stuff.
#include "ImplicitIntegrationStrategy.H"

#include "NamespaceHeader.H"

typedef BasicImplicitIntegrator<RegularImplicitIntegratorTraits> ImplicitIntegrator;

#include "NamespaceFooter.H"
#endif
