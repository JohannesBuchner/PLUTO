#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cerr;

#include "EBISLayout.H"
#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "EBAMRPoissonOp.H"
#include "AMRTGA.H"
#include "BaseBCValue.H"
#include "BaseDomainBC.H"
#include "slowEBCOFactory.H"
#include "slowEBCO.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "NeumannPoissonEBBC.H"
#include "ParmParse.H"
#include "MultilevelLinearOp.H"

#include "functionsF_F.H"
#include "TrigBCValue.H"
#include "TrigBCFlux.H"
#include "SphericalHarmonicBCValue.H"
#include "SphericalHarmonicBCFlux.H"
#include "MarshaValue.H"
#include "MarshaFlux.H"
#include    "EBAMRPoissonOpFactory.H"
#include "EBViscousTensorOpFactory.H"
#include  "EBConductivityOpFactory.H"

#ifndef _POISSONUTILITIES_H_
#define _POISSONUTILITIES_H_

///
/**
   In an attempt to avoid the extensive cut-and-paste reuse of
   problem definition stuff that happened in ebamrgodunov,
   here is a central repository of the stuff required to
   do an ebamrpoisson problem.  The single-level stuff was
   done earlier and will be left alone.
 */

class PoissonParameters
{
public:
  IntVect       nCells;
  int           maxGridSize;
  int           whichGeom;
  int           blockFactor;
  int           bufferSize;
  Real          fillRatio;
  int           ebBcType;
  int           domBcType;
  int           maxLevel;
  int           numLevels;
  Vector<int>   refRatio;
  ProblemDomain coarsestDomain;
  RealVect      coarsestDx;
  RealVect      domainLength;
  RealVect      probLo;
  RealVect      probHi;
  IntVect       ghostPhi;
  IntVect       ghostRHS;
  bool          noRefCorners;
  void coarsen(int a_factor);
  void  refine(int a_factor);

};
extern void
defineViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&                    a_solver,
                          const Vector<DisjointBoxLayout>&                        a_grids,
                          const Vector<EBISLayout>&                               a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&                    a_bottomSolver,
                          const PoissonParameters&                                a_params,
                          Vector<RefCountedPtr<LevelData<EBCellFAB> > >       &   aco,
                          Vector<RefCountedPtr<LevelData<EBFluxFAB> > >       &   eta,
                          Vector<RefCountedPtr<LevelData<EBFluxFAB> > >       &   lambda,
                          Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&   etaIrreg,
                          Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&   lambdaIrreg,
                          RefCountedPtr<EBViscousTensorOpFactory>&                a_opFactory);



extern void
newAndFillSmoothCellData(Vector<LevelData<EBCellFAB>* >&  a_scalar,
                         const Vector<DisjointBoxLayout>& a_grids,
                         const Vector<EBISLayout>&        a_ebisl,
                         const PoissonParameters&         a_params,
                         int ncomp,
                         bool a_useGhostPhi);

extern void
defineSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
             const Vector<DisjointBoxLayout>&             a_grids,
             const Vector<EBISLayout>&                    a_ebisl,
             LinearSolver<LevelData<EBCellFAB > >&        a_bottomSolver,
             const PoissonParameters&                     a_params,
             const Real&                                  a_time,
             Real                                         a_alpha,
             Real                                         a_beta,
             RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >& operatorFactory);

/*****/
extern void
coarsenBoxes(Vector< Vector<Box>      >&    a_boxesCoar,
             const Vector<Vector<Box> >&    a_boxesFine,
             int a_refToCoar);
extern void
getSimpleBoxes(Vector<Vector<Box> >&   a_boxes,
               Vector<int>&            a_refRat,
               const Box&              a_domain);

extern void
defineGrids(Vector< DisjointBoxLayout     >& a_grids,
            Vector< EBISLayout            >& a_ebisl,
            const Vector<Vector<Box>      >& a_boxes,
            const Vector<int              >& a_refRat,
            const PoissonParameters        & a_params);

extern int getOrderEB();

extern RealVect getTrigRV();
extern RealVect getSphericalHarmonicRV();

///get stuff from input file
extern void getPoissonParameters (      PoissonParameters&  a_params,
                                        bool a_forceSingleLevel = false);

///define ebindexspace from parameters (and parmparse)
extern void definePoissonGeometry(const PoissonParameters&  a_params);


///define the tga ubersolver.
/**
   Memory for solver gets allocated in this routine.  must be deleted later
 */
extern RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > > newTGASolver(const Vector<DisjointBoxLayout>&       a_grids,
                                                                  const Vector<EBISLayout>&              a_ebisl,
                                                                  const PoissonParameters&               a_params,
                                                                  BiCGStabSolver<LevelData<EBCellFAB> >& a_bottomSolver,
                                                                  const int&                             a_coarsestLevel,
                                                                  const Real&                            a_diffConst,
                                                                  const IntVect& a_phiGhost,
                                                                  const IntVect& a_rhsGhost);


///
/**
   Define layouts that in which all irregular cells are refined.
   No other refinement criteria are used.
 */
extern void getAllIrregRefinedLayouts(Vector<DisjointBoxLayout>& a_grids,
                                      Vector<EBISLayout>&        a_ebisl,
                                      const PoissonParameters&   a_params);

///
/**
*/
extern void compareError(const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
                         const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
                         const Vector< DisjointBoxLayout >&       a_gridsFine,
                         const Vector< DisjointBoxLayout >&       a_gridsCoar,
                         const Vector< EBISLayout >&              a_ebislFine,
                         const Vector< EBISLayout >&              a_ebislCoar,
                         const PoissonParameters&                 a_paramsFine,
                         const PoissonParameters&                 a_paramsCoar);

///
/**
*/
extern void getCoarseLayoutsFromFine(Vector<DisjointBoxLayout>&       a_gridsCoar,
                                     Vector< EBISLayout >&            a_ebislCoar,
                                     const Vector<DisjointBoxLayout>& a_gridsFine,
                                     const PoissonParameters&         a_params);

extern void getBCFactories(RefCountedPtr<BaseDomainBCFactory>& a_baseDomainBCFactory,
                           RefCountedPtr<BaseEBBCFactory>&     a_baseEBBCFactory,
                           const PoissonParameters&            a_params);

///
/**
*/
/**************/
extern void
getEBAMRPFactory(RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&  a_levelOpFactory,
                 const Vector<DisjointBoxLayout>&        a_grids,
                 const Vector< EBISLayout >&             a_ebisl,
                 const PoissonParameters&                a_params,
                 const int&                              a_numPreCondIters,
                 const int&                              a_relaxType,
                 const Real&                             a_time,
                 const Real&                             a_alpha=0.0,
                 const Real&                             a_beta = 1.0);

/********/
extern void
defineSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
             const Vector<DisjointBoxLayout>&             a_grids,
             const Vector<EBISLayout>&                    a_ebisl,
             LinearSolver<LevelData<EBCellFAB > >&        a_bottomSolver,
             const PoissonParameters&                     a_params,
             const Real&                                  a_time,
             Real                                         a_alpha=0.0,
             Real                                         a_beta=1.0);

/********/
extern void
defineMGBCGSolver(BiCGStabSolver<Vector<LevelData<EBCellFAB>* > >& a_solver,
                  MultilevelLinearOp<EBCellFAB>&                   a_mlOp,
                  const Vector<DisjointBoxLayout>&                 a_grids,
                  const Vector<EBISLayout>&                        a_ebisl,
                  const PoissonParameters&                         a_params,
                  const Real&                                      a_time,
                  Real                                             a_alpha,
                  Real                                             a_beta);
///
/**
*/
extern void setTrigPhi(LevelData<EBCellFAB>&     a_phi,
                       const  RealVect&          a_dx,
                       const  PoissonParameters& a_params,
                       const  Real& a_time=0.0,
                       int ivar = 0);

extern void setSphericalHarmonicPhi(LevelData<EBCellFAB>&     a_phi,
                                    const  RealVect&          a_dx,
                                    const  PoissonParameters& a_params,
                                    const  Real& a_time=0.0);

extern void setMarshaPhi(LevelData<EBCellFAB>&     a_phi,
                         const  RealVect&          a_dx,
                         const  PoissonParameters& a_params);

extern void setTrigKappaLOfPhi(LevelData<EBCellFAB>&    a_kappaLOfPhi,
                               const RealVect&          a_dx,
                               const PoissonParameters& a_params,
                               const Real&              a_alpha=0.0,
                               const Real&              a_beta=1.0,
                               const Real&              a_time=0.0,
                               int ivar = 0);

extern void setSphericalHarmonicKappaLOfPhi(LevelData<EBCellFAB>&    a_kappaLOfPhi,
                                            const RealVect&          a_dx,
                                            const PoissonParameters& a_params,
                                            const Real&              a_alpha=0.0,
                                            const Real&              a_beta=1.0,
                                            const Real&              a_time=0.0);

extern void setMarshaKappaLOfPhi(LevelData<EBCellFAB>&    a_kappaLOfPhi,
                                 const RealVect&          a_dx,
                                 const PoissonParameters& a_params);

extern void setTrigSource(LevelData<EBCellFAB>&    a_src,
                          const RealVect&          a_dx,
                          const PoissonParameters& a_params,
                          const Real&              a_diffConst,
                          const Real&              a_time=0.0);

extern void setSphericalHarmonicSource(LevelData<EBCellFAB>&    a_src,
                                       const RealVect&          a_dx,
                                       const PoissonParameters& a_params,
                                       const Real&              a_diffConst,
                                       const Real&              a_time=0.0);

extern void setRZPhi(LevelData<EBCellFAB>&    a_phi,
                     const RealVect&          a_dx,
                     const PoissonParameters& a_params);

extern void setRZKappaLOfPhi(LevelData<EBCellFAB>&    a_kappaLOfPhi,
                             const RealVect&          a_dx,
                             const PoissonParameters& a_params,
                             const Real&              a_alpha,
                             const Real&              a_beta);


extern void setKLVViscous  (LevelData<EBCellFAB>&    a_kappaLOfMag,
                            const Real&              a_dx,
                            const PoissonParameters& a_params);

extern void setVelViscous  (LevelData<EBCellFAB>&    a_kappaLOfMag,
                            const Real&          a_dx,
                            const PoissonParameters& a_params);
extern void
defineViscousTensorCoef(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&         a_aco,
                        Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&         a_eta,
                        Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&         a_lambda,
                        Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&  a_etaIrreg,
                        Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&  a_lambdaIrreg,
                        const Vector<DisjointBoxLayout>&                       a_grids,
                        const Vector<EBISLayout>&                              a_ebisl,
                        const PoissonParameters&                               a_params);
extern void
defineViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
                          const Vector<DisjointBoxLayout>&             a_grids,
                          const Vector<EBISLayout>&                    a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&         a_bottomSolver,
                          const PoissonParameters&                     a_params);


extern void
defineViscousTensorSolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
                          RefCountedPtr<EBViscousTensorOpFactory>&     a_factory,
                          const Vector<DisjointBoxLayout>&             a_grids,
                          const Vector<EBISLayout>&                    a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&         a_bottomSolver,
                          const PoissonParameters&                     a_params);

extern void
getEBVTOFactory(RefCountedPtr<EBViscousTensorOpFactory>&                       a_factory,
                const Vector<RefCountedPtr<LevelData<EBCellFAB > > >&          a_acoeffi,
                const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_eta,
                const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_lambda,
                const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_etaIrreg,
                const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_lambdaIrreg,
                const Vector<DisjointBoxLayout>&                               a_grids,
                const Vector<EBISLayout>&                                      a_ebisl,
                const PoissonParameters&                                       a_params);

///
extern
void getViscousBCFactories(RefCountedPtr<BaseDomainBCFactory>&                     a_domBC,
                           RefCountedPtr<BaseEBBCFactory>&                         a_ebBC,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&    a_eta,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&    a_lambda,
                           const Vector<DisjointBoxLayout>&                        a_grids,
                           const Vector<EBISLayout>&                               a_ebisl,
                           const PoissonParameters&                                a_params);

extern void setViscousCoefs(LevelData<EBCellFAB >&             a_acoeffi,
                            LevelData<EBFluxFAB>&              a_eta,
                            LevelData<EBFluxFAB>&              a_lambda,
                            LevelData<BaseIVFAB<Real> >   &    a_etaIrreg,
                            LevelData<BaseIVFAB<Real> >   &    a_lambdaIrreg,
                            const  Real&                       a_dx,
                            const  PoissonParameters&          a_params);

extern void setMagResistive(LevelData<EBCellFAB>&     a_mag,
                            const  Real&              a_dx,
                            const  PoissonParameters& a_params,
                            int a_whichMag = -1);

extern void setEtaResistive(LevelData<EBFluxFAB>&              a_eta,
                            LevelData<BaseIVFAB<Real> >   &    a_etaIrreg,
                            const  Real&                       a_dx,
                            const  PoissonParameters&          a_params,
                            int a_whichEta = -1);

extern void
defineConductivitySolver( AMRMultiGrid<LevelData<EBCellFAB> >&                     a_solver,
                          const Vector<DisjointBoxLayout>&                         a_grids,
                          const Vector<EBISLayout>&                                a_ebisl,
                          LinearSolver<LevelData<EBCellFAB> >&                     a_bottomSolver,
                          const PoissonParameters&                                 a_params,
                          Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                          Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                          Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                          RefCountedPtr<EBConductivityOpFactory>&                  a_opFactory);

extern void
defineSlowConductivitySolver( AMRMultiGrid<LevelData<EBCellFAB> >&                     a_solver,
                              const Vector<DisjointBoxLayout>&                         a_grids,
                              const Vector<EBISLayout>&                                a_ebisl,
                              LinearSolver<LevelData<EBCellFAB> >&                     a_bottomSolver,
                              const PoissonParameters&                                 a_params,
                              Vector<RefCountedPtr<LevelData<EBCellFAB> > >&           a_aco,
                              Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bco,
                              Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                              RefCountedPtr<slowEBCOFactory>&                          a_opFactory);


extern void
defineConductivityCoef(Vector<RefCountedPtr<LevelData<EBCellFAB> > >&         a_aco,
                       Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&         a_bco,
                       Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&  a_bcoIrreg,
                       const Vector<DisjointBoxLayout>&                       a_grids,
                       const Vector<EBISLayout>&                              a_ebisl,
                       const PoissonParameters&                               a_params);

extern void setConductivityCoefs(LevelData<EBCellFAB >&             a_aco,
                                 LevelData<EBFluxFAB>&              a_bco,
                                 LevelData<BaseIVFAB<Real> >   &    a_bcoIrreg,
                                 const  Real&                       a_dx,
                                 const  PoissonParameters&          a_params);

extern void
getConductivityFactory(RefCountedPtr<EBConductivityOpFactory>&                        a_factory,
                       const Vector<RefCountedPtr<LevelData<EBCellFAB > > >&          a_acoeffi,
                       const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bcoeffi,
                       const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                       const Vector<DisjointBoxLayout>&                               a_grids,
                       const Vector<EBISLayout>&                                      a_ebisl,
                       const PoissonParameters&                                       a_params);


extern void
getSlowConductivityFactory(RefCountedPtr<slowEBCOFactory>&                                a_factory,
                           const Vector<RefCountedPtr<LevelData<EBCellFAB > > >&          a_acoeffi,
                           const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&           a_bcoeffi,
                           const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&    a_bcoIrreg,
                           const Vector<DisjointBoxLayout>&                               a_grids,
                           const Vector<EBISLayout>&                                      a_ebisl,
                           const PoissonParameters&                                       a_params);

extern void
defineConductivitySolver(AMRMultiGrid<LevelData<EBCellFAB> >&         a_solver,
                         const Vector<DisjointBoxLayout>&             a_grids,
                         const Vector<EBISLayout>&                    a_ebisl,
                         LinearSolver<LevelData<EBCellFAB> >&         a_bottomSolver,
                         const PoissonParameters&                     a_params);

extern
void getConductivityBCFactories(RefCountedPtr<BaseDomainBCFactory>&                     a_domBC,
                                RefCountedPtr<BaseEBBCFactory>&                         a_ebBC,
                                const Vector<RefCountedPtr<LevelData<EBFluxFAB> > >&    a_bcoef,
                                const Vector<DisjointBoxLayout>&                        a_grids,
                                const Vector<EBISLayout>&                               a_ebisl,
                                const PoissonParameters&                                a_params);

Real getDistance(const RealVect& a_point,
                 const RealVect& a_lineDirection,
                 const RealVect& a_linePoint);

#endif
