#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEIMPLICITINTEGRATIONSTRATEGY_H_
#define _BASEIMPLICITINTEGRATIONSTRATEGY_H_

#include <iostream>
#include "SPACE.H"
#include "REAL.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRMultiGrid.H"
#include "AMRTGA.H"

#include "NamespaceHeader.H"

//! \class BaseImplicitIntegrationStrategy
//! This base class is a strategy pattern that represents the process of
//! implicitly integrating an elliptic partial differential equation
//! \f$a(\vec{x})\frac{\partial \phi}{\partial t} = L(\phi, t) + S(\vec{x}, t)$
//! on a single AMR grid level using Helmholtz operators of the form
//! \f$H(\phi) = \alpha a(\vec{x})\phi + \beta L(\phi)$ to represent the time
//! discretization. If you wish to implement your favorite implicit time integration
//! scheme, you can do so by creating a subclass of BaseImplicitIntegrationStrategy.
//! This class contains several helpful methods you can use in your implementation.
//! NOTE: These integration strategies make use of LevelTGAHelmOps, but these
//! NOTE: Helmholtz operators \f${H}$ have nothing to do with the TGA integration
//! NOTE: scheme. This remains a naming weirdness in Chombo.
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
class BaseImplicitIntegrationStrategy
{

  typedef typename ImplicitIntegratorTraits::LevelDataType LevelDataType;
  typedef typename ImplicitIntegratorTraits::FluxDataType FluxDataType;
  typedef typename ImplicitIntegratorTraits::FluxRegisterType FluxRegisterType;
  typedef typename ImplicitIntegratorTraits::GridPolicy GridPolicy;
  typedef typename ImplicitIntegratorTraits::GridType GridType;
  typedef typename ImplicitIntegratorTraits::HelmholtzOpType HelmholtzOpType;

  public:

  //! Base class constructor for implicit integration strategies. This must
  //! be called first and foremost in the constructor of any subclass of
  //! BaseImplicitIntegrationStrategy.
  BaseImplicitIntegrationStrategy();

  //! Destructor, called after the destructor of an BaseImplicitIntegrationStrategy
  //! subclass.
  virtual ~BaseImplicitIntegrationStrategy();

  //! Sets the AMR hierarchy of disjoint box layouts to be used by the
  //! strategy.
  //! \param a_grids A pointer to a vector of disjoint box layouts.
  void setGrids(const Vector<GridType>& a_grids);

  //! Sets the multigrid solver for the strategy.
  //! \param a_op A pointer to a LevelTGAHelmOp to be used as the Helmholtz
  //!             operator.
  void setMultiGridSolver(RefCountedPtr<AMRMultiGrid<LevelDataType> >& a_solver);

  //! Sets the current AMR refinement level for the strategy.
  //! \param a_level The level on which this strategy will operate.
  void setAMRLevel(int a_level);

  //! Sets the Helmholtz operator that this class uses to perform its
  //! time integration.
  //! \param a_op A pointer to a LevelTGAHelmOp to be used as the Helmholtz
  //!             operator.
  void setHelmholtzOp(HelmholtzOpType* a_op);

  //! Returns the Helmholtz operator that this class uses to perform its
  //! time integration.
  HelmholtzOpType& helmholtzOp();

  //! Initializes the flux for keeping track of contributions to flux registers.
  //! \param a_numComps The number of components in the solution for which the
  //!                   flux is to be accumulated.
  void initFlux(int a_numComps);

  //! Returns a reference to the object containing the contributions of diffusive
  //! fluxes used in the integration of the solution.
  const LevelData<FluxDataType>& flux() const;

  //! Sets a flag that determines whether the solution will be zeroed when
  //! the Helmholtz equation is solved, or whether it will be used as an initial
  //! guess.
  //! \param a_initialGuess If \a a_zeroSolution is false, this data will be
  //!                       saved as an initial estimate of the solution.
  //! \param a_zeroSolution If set to true, the solution will be set to zero
  //!                       before all solves of the Helmholtz equation.
  //!                       Otherwise, it will be used as initial guess.
  void zeroInitialSolution(const LevelDataType& a_initialGuess,
                           bool a_zeroSolution);

  //! Sets pointers for coarser-level solution data in order to evaluate
  //! boundary conditions at times within the current step.
  void setCoarseData(const LevelDataType* a_oldCoarseDataPtr,
                     const LevelDataType* a_newCoarseDataPtr,
                     Real a_oldCoarseTime,
                     Real a_newCoarseTime);

  //! Integrates the equation
  //! \f$a(\vec{x})\frac{\partial \phi}{\partial t} = L(\phi, t) + S(\vec{x}, t)$,
  //! producing a value for \f$\phi^{n+1}$ from \f$\phi^n$.
  //! \param a_phiNew The new solution \f$\phi^{n+1} will be stored here.
  //! \param a_phiOld The old solution \f$\phi^n$.
  //! \param a_src The source term $S$ on the right hand side of the Helmholtz
  //!              equation.
  //! \param a_oldTime The time \f$t^n$ at the beginning of the integration
  //!                  step \f$n$ at the current grid level.
  //! \param a_dt The size \f$\Delta t$ of the integration step at the current
  //!             grid level.
  virtual void updateSolution(LevelDataType&           a_phiNew,
                              LevelDataType&           a_phiOld,
                              LevelDataType&           a_src,
                              Real                     a_oldTime,
                              Real                     a_dt) = 0;

  //! Computes the time-centered diffusion term
  //! \f$L(\phi) = a(\vec{x})\frac{\partial \phi}{\partial t} - S$.
  //! \param a_diffusiveTerm The diffusion term \f$L(\phi)$ will be stored here.
  //! \param a_phiOld The solution \f$\phi^n$ at time \f$n$.
  //! \param a_src The source term \f$S$ on the right hand side of the
  //!              diffusion equation.
  //! \param a_oldTime The time \f$t^n$ at the beginning of the integration step
  //!                  at the current grid level.
  //! \param a_dt The size \f$\Delta t$ of the integration step at the current grid level.
  virtual void computeDiffusion(LevelDataType&           a_diffusiveTerm,
                                LevelDataType&           a_phiOld,
                                LevelDataType&           a_src,
                                Real                     a_oldTime,
                                Real                     a_dt) = 0;

  //! Explicitly computes \f$L(\phi)$, applying the diffusion operator \f$L$
  //! to the solution \f$\phi$ at the given time \f$t$.
  //! \param a_Lphi The computed value of \f$L(\phi)$ is stored here.
  //! \param a_phi The solution $\phi$ to which the diffusion operator
  //!              is applied.
  //! \param a_applyBC If set to true, inhomogeneous boundary conditions will
  //!                  be applied to the Helmholtz operator (both at the domain
  //!                  boundary and at the coarse-fine grid interfaces). If
  //!                  set to false, homogeneous boundary conditions will
  //!                  applied to the Helmholtz operator at these locations.
  //! \param a_time The time \f$t$ at which the boundary conditions are applied.
  //!               If \a a_applyBC is false, this argument is ignored.
  void applyDiffusionOperator(LevelDataType&          a_Lphi,
                              const LevelDataType&    a_phi,
                              Real                    a_time,
                              bool                    a_applyBC);

  protected:

  //! Returns the current AMR refinement level.
  int level() const;

  //! Returns the current DisjointBoxLayout object.
  DisjointBoxLayout disjointBoxLayout() const;

  //! Returns the current grid object, which is a DisjointBoxLayout
  //! for regular geometry and an EBLevelGrid for EB geometry.
  const GridType& grid() const;

  //! Sets the value of the source term in the Helmholtz equation on ghost
  //! cells. This method should be overridden by subclasses of
  //! BaseImplicitIntegrationStrategy.
  //! \param a_src The source term in the Helmholtz equation whose ghost cell
  //!              values are to be set.
  virtual void setSourceGhostCells(LevelDataType& a_src) = 0;

  //! Sets the parameters \f$\alpha$ and \f$\beta$ in the Helmholtz operator
  //! \f$\alpha a(\vec{x}) I + \beta L$ to the given values.
  //! \param a_alpha The identity term coefficient in the Helmholtz operator.
  //! \param a_beta The coefficient in front of the discrete derivative term
  //!               in the Helmholtz operator.
  void setAlphaAndBeta(Real a_alpha,
                       Real a_beta);

  //! Solves the Helmholtz equation
  //! \f$[\alpha a(\vec{x}) - \mu \Delta t * L]\phi = S$
  //! for \f$\phi$. Boundary conditions are evaluated at the time \f$t$.
  //! \a a_time + \a a_stepFrac * \a a_dt. The diffusive flux is accounted for
  //! in this operation.
  //! \param a_phi The solution $\phi$ to the Helmholtz equation.
  //! \param a_src The source term \f$S$ in the Helmholtz equation.
  //! \param a_mu A number between 0 and 1 that defines the time within the
  //!             integration step at which the Helmholtz operator is evaluated.
  //! \param a_dt The size of the integration step.
  //! \param a_time The time \f$t$ at which boundary conditions are applied.
  void solveHelmholtz(LevelDataType& a_phi,
                      LevelDataType& a_src,
                      Real           a_mu,
                      Real           a_dt,
                      Real           a_time);

  //! Applies the Helmholtz operator to the solution \phi, setting
  //! \a a_ans to \f$I + \mu \Delta t L(\phi)$ and applying boundary
  //! conditions at time \f$t$. NOTE that the operator's time is not set by
  //! this method. You must set the time on the operator separately.
  //! \param a_ans The value of the Helmholtz operator applied to
  //!              to \f$\phi$: \f$I + \mu L(\phi)$.
  //! \param a_phi The solution \f$\phi$ to which the Helmholtz operator
  //!              is applied.
  //! \param a_mu \f$\mu$, the factor multiplying the operator L in the
  //!             Helmholtz expression.
  //! \param a_dt \f$\Delta t$, the size of the integration step.
  //! \param a_applyBC If set to true, inhomogeneous boundary conditions will
  //!                  be applied to the Helmholtz operator (both at the domain
  //!                  boundary and at the coarse-fine grid interfaces). If
  //!                  set to false, homogeneous boundary conditions will
  //!                  applied to the Helmholtz operator at these locations.
  //! \param a_time \f$t$, the time at which boundary conditions are applied.
  //!               If \a a_applyBC is false, this argument is ignored.
  void applyHelmholtz(LevelDataType&          a_ans,
                      const LevelDataType&    a_phi,
                      Real                    a_mu,
                      Real                    a_dt,
                      bool                    a_applyBC,
                      Real                    a_time);

  private:

  //! Adds a flux contribution from the Helmholtz operator to m_flux.
  //! \param a_phi The solution to which the Helmholtz operator is applied
  //!              in order to evaluate the flux.
  //! \param a_mu A number between 0 and 1 that defines the time within the
  //!             integration step at which the Helmholtz operator is evaluated.
  //! \param a_dt The size of the integration step.
  void incrementFlux(const LevelDataType& a_phi,
                     Real a_mu,
                     Real a_dt);

  //! Sets the time at which boundary conditions are evaluated.
  //! \param a_time The time to which BCs are set.
  void setBoundaryConditionTime(Real a_time);

  //! The disjoint box layouts at every AMR grid level.
  Vector<GridType> m_grids;

  //! The current AMR refinement level. Needed to use the multigrid solver.
  int m_level;

  //! The Helmholtz operator we use to solve the Helmholtz equation
  //! at a particular AMR level.
  LevelTGAHelmOp<LevelDataType, FluxDataType>* m_op;

  //! The multigrid solver used to solve the Helmholtz equation.
  RefCountedPtr<AMRMultiGrid<LevelDataType> > m_solver;

  //! A flag that dictates whether the initial solution to Helmholtz solves
  //! should be zeroed.
  bool m_zeroInitialSolution;

  //! The initial guess for the solution (if it exists).
  LevelDataType m_initialGuess;

  // Coarse data for boundary conditions.
  LevelDataType m_coarseData, *m_oldCoarseDataPtr, *m_newCoarseDataPtr;
  Real m_oldCoarseTime, m_newCoarseTime;

  // This keeps track of flux contributions.
  LevelData<FluxDataType> m_flux;

  // Disallowed operators.
  BaseImplicitIntegrationStrategy(const BaseImplicitIntegrationStrategy&);
  BaseImplicitIntegrationStrategy& operator=(const BaseImplicitIntegrationStrategy&);
};


//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
BaseImplicitIntegrationStrategy():
  m_grids(),
  m_level(0),
  m_op(NULL),
  m_solver(),
  m_zeroInitialSolution(false),
  m_initialGuess(),
  m_coarseData(),
  m_oldCoarseDataPtr(NULL),
  m_newCoarseDataPtr(NULL),
  m_oldCoarseTime(0.0),
  m_newCoarseTime(0.0),
  m_flux()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
~BaseImplicitIntegrationStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
setGrids(const Vector<GridType>& a_grids)
{
  m_grids = a_grids;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
setMultiGridSolver(RefCountedPtr<AMRMultiGrid<LevelDataType> >& a_solver)
{
  m_solver = a_solver;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
setAMRLevel(int a_level)
{
  CH_assert((a_level >= 0) && (a_level < m_grids.size()));
  m_level = a_level;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
setHelmholtzOp(typename ImplicitIntegratorTraits::HelmholtzOpType* a_op)
{
  CH_assert(a_op != NULL);
  m_op = a_op;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
typename ImplicitIntegratorTraits::HelmholtzOpType&
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
helmholtzOp()
{
  CH_assert(m_op != NULL);
  return *m_op;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
initFlux(int a_numComps)
{
  // Define a flux for this AMR level and initialize it to zero.
  GridPolicy::defineFlux(m_flux, m_grids[m_level], a_numComps, IntVect::Zero);
  const DisjointBoxLayout& dbl = GridPolicy::disjointBoxLayout(m_grids[m_level]);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    m_flux[dit()].setVal(0.0);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
const LevelData<typename BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::FluxDataType>&
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
flux() const
{
  return m_flux;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
int
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
level() const
{
  CH_assert(m_level >= 0);
  CH_assert(m_level < m_grids.size());
  return m_level;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
DisjointBoxLayout
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
disjointBoxLayout() const
{
  CH_assert(m_level >= 0);
  CH_assert(m_level < m_grids.size());
  return GridPolicy::disjointBoxLayout(m_grids[m_level]);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
const typename ImplicitIntegratorTraits::GridType&
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
grid() const
{
  CH_assert(m_level >= 0);
  CH_assert(m_level < m_grids.size());
  return m_grids[m_level];
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
zeroInitialSolution(const LevelDataType& a_initialGuess,
                    bool a_zeroSolution)
{
  m_zeroInitialSolution = a_zeroSolution;
  if (!a_zeroSolution)
  {
    CH_assert(m_op != NULL);
    m_op->create(m_initialGuess, a_initialGuess);
    m_op->assign(m_initialGuess, a_initialGuess);
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
setCoarseData(const LevelDataType* a_oldCoarseDataPtr,
              const LevelDataType* a_newCoarseDataPtr,
              Real a_oldCoarseTime,
              Real a_newCoarseTime)
{
  CH_assert(((a_oldCoarseDataPtr != NULL) && (a_newCoarseDataPtr != NULL)) ||
            ((a_oldCoarseDataPtr == NULL) && (a_newCoarseDataPtr == NULL)));
  CH_assert(a_oldCoarseTime <= a_newCoarseTime);
  m_oldCoarseDataPtr = const_cast<LevelDataType*>(a_oldCoarseDataPtr);
  m_newCoarseDataPtr = const_cast<LevelDataType*>(a_newCoarseDataPtr);
  m_oldCoarseTime = a_oldCoarseTime;
  m_newCoarseTime = a_newCoarseTime;

  if (m_oldCoarseDataPtr != 0)
  {
    // Make sure coarse data is defined at the right grid level.
    CH_assert(m_level > 0);
    GridPolicy::defineData(m_coarseData, m_grids[m_level-1],
                           m_oldCoarseDataPtr->nComp(), 1*IntVect::Unit);
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
setAlphaAndBeta(Real a_alpha,
                Real a_beta)
{
  CH_assert(m_op != NULL);
  CH_assert(!m_solver.isNull());
  Vector<MGLevelOp<LevelDataType>* > ops = m_solver->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
  {
    LevelTGAHelmOp<LevelDataType, FluxDataType>* helmop =
      dynamic_cast<LevelTGAHelmOp<LevelDataType, FluxDataType>*>(ops[iop]);
    helmop->setAlphaAndBeta(a_alpha, a_beta);
  }
  m_op->setAlphaAndBeta(a_alpha, a_beta);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
applyHelmholtz(LevelDataType&       a_ans,
               const LevelDataType& a_phi,
               Real                 a_mu,
               Real                 a_dt,
               bool                 a_applyBC,
               Real                 a_time)
{
  // Set the operator's alpha and beta coefficients.
  Real factor  = a_mu*a_dt;
  m_op->setAlphaAndBeta(1.0, factor);

  // Zero the answer.
  m_op->setToZero(a_ans);

  // set a_ans = helm(a_phi)
  //           = (I + factor*laplacian)(a_phi)
  if (a_applyBC)
  {
    setBoundaryConditionTime(a_time); // Set the proper time on BCs.
    if (m_level == 0)
    {
      m_op->applyOp(a_ans, a_phi, false);
    }
    else
    {
      m_op->AMROperatorNF(a_ans, a_phi, m_coarseData, false);
    }
  }
  else
  {
    m_op->applyOpNoBoundary(a_ans, a_phi);
  }

  // Add the diffusive flux contribution.
  incrementFlux(a_phi, a_mu, a_dt);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
setBoundaryConditionTime(Real a_time)
{
  if ((m_oldCoarseDataPtr != NULL) && (m_newCoarseDataPtr != NULL))
  {
    Real eps = 1.0e-10;
    Real diff = (m_newCoarseTime - m_oldCoarseTime);
    m_op->setToZero(m_coarseData);
    if (diff < eps)
    {
      //no real time advance and don't want to divide by zero
      m_op->incr(m_coarseData, *m_oldCoarseDataPtr, 1.0);
    }
    else
    {
      Real factor = (a_time - m_oldCoarseTime)/(m_newCoarseTime - m_oldCoarseTime);
      for (DataIterator dit = m_coarseData.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      {
        m_coarseData[dit()].axby((*m_oldCoarseDataPtr)[dit()],
                                 (*m_newCoarseDataPtr)[dit()],
                                 1.0-factor, factor);
      }
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
solveHelmholtz(LevelDataType& a_phi,
               LevelDataType& a_src,
               Real           a_mu,
               Real           a_dt,
               Real           a_time)
{
  // Zero a_phi before we start, or copy the initial guess to it.
  if (m_zeroInitialSolution)
    m_op->setToZero(a_phi);
  else
    m_op->assign(a_phi, m_initialGuess);

  // Make sure the proper time is set for boundary conditions.
  setBoundaryConditionTime(a_time);

  Vector<LevelDataType* > phi(m_grids.size(), NULL);
  Vector<LevelDataType* > rhs(m_grids.size(), NULL);
  phi[m_level] = &a_phi;
  rhs[m_level] = &a_src;
  if (m_level > 0)
  {
    phi[m_level-1] = &m_coarseData;
  }

  Real factor  = -a_dt*a_mu;
  setAlphaAndBeta(1.0, factor);

  m_solver->solve(phi, rhs, m_level, m_level, m_zeroInitialSolution);
  int solverExitStatus = m_solver->m_exitStatus;
  if (solverExitStatus==2 || solverExitStatus==4 || solverExitStatus==6)
  {
    // These status codes correspond to the cases in which
    // norm is neither small enough nor reduced enough.
    // Either we've reached the maximum number of iterations,
    // or we've hung.
    pout() << "BaseImplicitIntegrationStrategy:: WARNING: solver exitStatus == "
      << solverExitStatus << std::endl;
  }

  // Add the diffusive flux contribution.
  incrementFlux(a_phi, a_mu, a_dt);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
incrementFlux(const LevelDataType& a_phi,
              Real a_mu,
              Real a_dt)
{
  // NOTE that elliptical operators all have negative flux contributions
  // since they are brought to the LHS from the RHS.
  Real factor = -a_dt*a_mu;
  m_op->setAlphaAndBeta(1.0, factor);

  // increment flux
  FluxDataType tempFlux;
  m_op->fillGrad(a_phi);
  DisjointBoxLayout dbl = GridPolicy::disjointBoxLayout(m_grids[m_level]);
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  {
    FluxDataType& thisFlux = m_flux[dit];
    m_op->getFlux(tempFlux, a_phi, dbl[dit], dit(), 1.0);
    thisFlux += tempFlux;
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename ImplicitIntegratorTraits>
void
BaseImplicitIntegrationStrategy<ImplicitIntegratorTraits>::
applyDiffusionOperator(LevelDataType& a_Lphi,
                       const LevelDataType& a_phi,
                       Real           a_time,
                       bool           a_applyBC)
{
  // Set the coefficients.
  setAlphaAndBeta(0.0, 1.0);

  // Zero the answer.
  m_op->setToZero(a_Lphi);

  // Make sure the proper time is set.
  setBoundaryConditionTime(a_time);
  m_op->setTime(a_time, 0.0, 0.0);

  // set a_Lphi = helm(a_phi)
  //            = laplacian(a_phi)
  if (a_applyBC)
    {
      if (m_level == 0)
        {
          pout() << "apply diffusion operator for level" << m_level << endl;
          m_op->applyOp(a_Lphi, a_phi, false);
        }
      else
        {
          pout() << "apply diffusion operator (NF) for level" << m_level << endl;
          m_op->AMROperatorNF(a_Lphi, a_phi, m_coarseData, false);
        }
    }
  else
    {
      m_op->applyOpNoBoundary(a_Lphi, a_phi);
    }
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
#endif
