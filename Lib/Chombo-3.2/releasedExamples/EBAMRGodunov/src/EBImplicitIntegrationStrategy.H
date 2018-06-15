#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBIMPLICITINTEGRATIONSTRATEGY_H_
#define _EBIMPLICITINTEGRATIONSTRATEGY_H_

#include "BaseImplicitIntegrationStrategy.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBFluxFAB.H"
#include "EBFluxRegister.H"
#include "EBFluxFactory.H"
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"

//! \struct EBGridPolicy
//! This policy class regularizes grid objects, giving callers access to
//! to DisjointBoxLayout objects regardless of the grid type.
struct EBGridPolicy
{
  //! Given a "grid" object, return the corresponding disjoint box layout.
  //! \param a_grid A grid object.
  static DisjointBoxLayout disjointBoxLayout(const EBLevelGrid& a_grid)
  {
    return a_grid.getDBL();
  }

  //! Define the given level data.
  static void defineData(LevelData<EBCellFAB>& a_data,
                         const EBLevelGrid& a_grid,
                         int a_ncomp,
                         const IntVect& a_ghostStencil)
  {
    EBCellFactory fact(a_grid.getEBISL());
    a_data.define(a_grid.getDBL(), a_ncomp, a_ghostStencil, fact);
  }

  //! Define the given flux FAB.
  static void defineFlux(LevelData<EBFluxFAB>& a_flux,
                         const EBLevelGrid& a_grid,
                         int a_ncomp,
                         const IntVect& a_ghostStencil)
  {
    EBFluxFactory fact(a_grid.getEBISL());
    a_flux.define(a_grid.getDBL(), a_ncomp, a_ghostStencil, fact);
  }
};

//! \struct EBImplicitIntegratorTraits
//! This class spells out the data types used by EB implicit integrators.
struct EBImplicitIntegratorTraits
{
  typedef LevelData<EBCellFAB> LevelDataType;
  typedef EBFluxFAB FluxDataType;
  typedef EBFluxRegister FluxRegisterType;
  typedef EBGridPolicy GridPolicy;
  typedef EBLevelGrid GridType;
  typedef LevelTGAHelmOp<LevelDataType, FluxDataType> HelmholtzOpType;
  typedef AMRLevelOpFactory<LevelDataType> OpFactoryType;
  typedef AMRMultiGrid<LevelDataType> SolverType;
};

//! \class EBImplicitIntegrationStrategy
//! This base class implements the first-order Backward Euler method for
//! embedded boundary algorithms.
class EBImplicitIntegrationStrategy:
  public BaseImplicitIntegrationStrategy<EBImplicitIntegratorTraits>
{

  public:

  typedef EBImplicitIntegratorTraits::LevelDataType LevelDataType;
  typedef EBImplicitIntegratorTraits::FluxDataType FluxDataType;
  typedef EBImplicitIntegratorTraits::FluxRegisterType FluxRegisterType;
  typedef EBImplicitIntegratorTraits::GridPolicy GridPolicy;
  typedef EBImplicitIntegratorTraits::GridType GridType;
  typedef EBImplicitIntegratorTraits::HelmholtzOpType HelmholtzOpType;

  //! Base constructor for embedded-boundary implicit time integration
  //! strategies. This must be called by all subclasses.
  EBImplicitIntegrationStrategy();

  //! Destructor.
  ~EBImplicitIntegrationStrategy();

  void setSourceGhostCells(LevelData<EBCellFAB>& a_src);

};

#include "NamespaceFooter.H"
#endif
