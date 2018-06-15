#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBIMPLICITINTEGRATORFACTORY_H_
#define _EBIMPLICITINTEGRATORFACTORY_H_

#include "EBImplicitIntegrator.H"

#include "NamespaceHeader.H"

//! \struct EBImplicitIntegratorFactory
//! This factory class creates reference-counted pointer to fresh
//! instances of implicit integrators for use with Embedded Boundary
//! methods.
class EBImplicitIntegratorFactory
{

  public:

  //! Returns a reference-counted pointer to a new Crank-Nicholson integrator.
  //! \param a_opFactory The factory used to construct Helmholtz operators.
  //! \param a_solver The Multigrid solver used to solve the elliptic system.
  //! \param a_grids A vector of EBLevelGrids describing the system.
  static RefCountedPtr<EBImplicitIntegrator>
  CrankNicholson(RefCountedPtr<EBImplicitIntegrator::OpFactoryType>& a_opFactory,
                 RefCountedPtr<EBImplicitIntegrator::SolverType>& a_solver,
                 const Vector<EBImplicitIntegrator::GridType>& a_grids);

  //! Returns a reference-counted pointer to a new Backward Euler integrator.
  //! \param a_opFactory The factory used to construct Helmholtz operators.
  //! \param a_solver The Multigrid solver used to solve the elliptic system.
  //! \param a_grids A vector of EBLevelGrids describing the system.
  static RefCountedPtr<EBImplicitIntegrator>
  backwardEuler(RefCountedPtr<EBImplicitIntegrator::OpFactoryType>& a_opFactory,
                RefCountedPtr<EBImplicitIntegrator::SolverType>& a_solver,
                const Vector<EBImplicitIntegrator::GridType>& a_grids);

  //! Returns a reference-counted pointer to a new Twizell-Gumel-Arigu integrator.
  //! \param a_opFactory The factory used to construct Helmholtz operators.
  //! \param a_solver The Multigrid solver used to solve the elliptic system.
  //! \param a_grids A vector of EBLevelGrids describing the system.
  static RefCountedPtr<EBImplicitIntegrator>
  TGA(RefCountedPtr<EBImplicitIntegrator::OpFactoryType>& a_opFactory,
      RefCountedPtr<EBImplicitIntegrator::SolverType>& a_solver,
      const Vector<EBImplicitIntegrator::GridType>& a_grids);

  private:

  // Not instantiable!
  EBImplicitIntegratorFactory();
  EBImplicitIntegratorFactory(const EBImplicitIntegratorFactory&);
  EBImplicitIntegratorFactory& operator=(const EBImplicitIntegratorFactory&);

  // Nor destructible.
  ~EBImplicitIntegratorFactory();
};

#include "NamespaceFooter.H"
#endif
