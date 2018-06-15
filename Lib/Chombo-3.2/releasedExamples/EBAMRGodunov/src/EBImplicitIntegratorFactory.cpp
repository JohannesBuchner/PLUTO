#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBImplicitIntegratorFactory.H"
#include "EBCrankNicholsonStrategy.H"
#include "EBTGAStrategy.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
RefCountedPtr<EBImplicitIntegrator>
EBImplicitIntegratorFactory::
CrankNicholson(RefCountedPtr<EBImplicitIntegrator::OpFactoryType>& a_opFactory,
               RefCountedPtr<EBImplicitIntegrator::SolverType>& a_solver,
               const Vector<EBImplicitIntegrator::GridType>& a_grids)
{
  RefCountedPtr<BaseImplicitIntegrationStrategy<EBImplicitIntegratorTraits> >
    strategy(new EBCrankNicholsonStrategy(0.5));
  RefCountedPtr<EBImplicitIntegrator> integrator(
    new EBImplicitIntegrator(strategy, a_opFactory, a_solver, a_grids));
  return integrator;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RefCountedPtr<EBImplicitIntegrator>
EBImplicitIntegratorFactory::
backwardEuler(RefCountedPtr<EBImplicitIntegrator::OpFactoryType>& a_opFactory,
              RefCountedPtr<EBImplicitIntegrator::SolverType>& a_solver,
              const Vector<EBImplicitIntegrator::GridType>& a_grids)
{
  RefCountedPtr<BaseImplicitIntegrationStrategy<EBImplicitIntegratorTraits> >
    strategy(new EBCrankNicholsonStrategy(1.0));
  RefCountedPtr<EBImplicitIntegrator> integrator(
    new EBImplicitIntegrator(strategy, a_opFactory, a_solver, a_grids));
  return integrator;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RefCountedPtr<EBImplicitIntegrator>
EBImplicitIntegratorFactory::
TGA(RefCountedPtr<EBImplicitIntegrator::OpFactoryType>& a_opFactory,
    RefCountedPtr<EBImplicitIntegrator::SolverType>& a_solver,
    const Vector<EBImplicitIntegrator::GridType>& a_grids)
{
  RefCountedPtr<BaseImplicitIntegrationStrategy<EBImplicitIntegratorTraits> >
    strategy(new EBTGAStrategy());
  RefCountedPtr<EBImplicitIntegrator> integrator(
    new EBImplicitIntegrator(strategy, a_opFactory, a_solver, a_grids));
  return integrator;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
