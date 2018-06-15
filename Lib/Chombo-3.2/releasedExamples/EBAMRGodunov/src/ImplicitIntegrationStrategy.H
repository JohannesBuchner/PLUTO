#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IMPLICITINTEGRATIONSTRATEGY_H_
#define _IMPLICITINTEGRATIONSTRATEGY_H_

#include "BaseImplicitIntegrationStrategy.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelFluxRegister.H"

#include "NamespaceHeader.H"

//! \struct RegularGridPolicy
//! This class regularizes grid objects, giving callers access to
//! DisjointBoxLayout objects regardless of the grid type.
struct RegularGridPolicy
{
  //! Given a "grid" object, return the corresponding disjoint box layout.
  //! \param a_grid A grid object.
  static const DisjointBoxLayout& disjointBoxLayout(const DisjointBoxLayout& a_grid)
  {
    return a_grid;
  }

  //! Define the given flux FAB.
  static void defineData(LevelData<FArrayBox>& a_data,
                         const DisjointBoxLayout& a_grid,
                         int a_ncomp,
                         const IntVect& a_ghostStencil)
  {
    a_data.define(a_grid, a_ncomp, a_ghostStencil);
  }

  //! Define the given flux FAB.
  static void defineFlux(LevelData<FluxBox>& a_flux,
                         const DisjointBoxLayout& a_grid,
                         int a_ncomp,
                         const IntVect& a_ghostStencil)
  {
    a_flux.define(a_grid, a_ncomp, a_ghostStencil);
  }
};

//! \struct RegularImplicitIntegratorTraits
//! This class spells out the data types used by regular-geometry
//! implicit integrators.
struct RegularImplicitIntegratorTraits
{
  typedef LevelData<FArrayBox> LevelDataType;
  typedef FluxBox FluxDataType;
  typedef LevelFluxRegister FluxRegisterType;
  typedef RegularGridPolicy GridPolicy;
  typedef DisjointBoxLayout GridType;
  typedef LevelTGAHelmOp<LevelDataType, FluxDataType> HelmholtzOpType;
  typedef AMRLevelOpFactory<LevelDataType> OpFactoryType;
  typedef AMRMultiGrid<LevelDataType> SolverType;
};

//! \class ImplicitIntegrationStrategy
//! This base class implements the first-order Backward Euler method for
//! embedded boundary algorithms.
class ImplicitIntegrationStrategy:
  public BaseImplicitIntegrationStrategy<RegularImplicitIntegratorTraits>
{

  public:

  typedef RegularImplicitIntegratorTraits::LevelDataType LevelDataType;
  typedef RegularImplicitIntegratorTraits::FluxDataType FluxDataType;
  typedef RegularImplicitIntegratorTraits::FluxRegisterType FluxRegisterType;
  typedef RegularImplicitIntegratorTraits::GridPolicy GridPolicy;
  typedef RegularImplicitIntegratorTraits::GridType GridType;
  typedef RegularImplicitIntegratorTraits::HelmholtzOpType HelmholtzOpType;

  //! Base constructor for embedded-boundary implicit time integration
  //! strategies. This must be called by all subclasses.
  ImplicitIntegrationStrategy();

  //! Destructor.
  ~ImplicitIntegrationStrategy();

  void setSourceGhostCells(LevelDataType& a_src);

};

#include "NamespaceFooter.H"
#endif
