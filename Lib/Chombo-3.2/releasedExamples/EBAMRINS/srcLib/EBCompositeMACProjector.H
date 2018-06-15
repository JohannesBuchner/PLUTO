#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCOMPOSITEMACPROJECTOR_H_
#define _EBCOMPOSITEMACPROJECTOR_H_

#include "EBLevelMACProjector.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "Box.H"
#include "REAL.H"
#include "LevelData.H"
#include "EBFluxFAB.H"
#include "EBCellFAB.H"
#include "EBSimpleSolver.H"
#include "GMRESSolver.H"
#include "EBAMRPoissonOp.H"
#include "MultiGrid.H"
#include "EBPWLFillPatch.H"
#include "BiCGStabSolver.H"
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"

///
/**
   Projection to take out pure gradient component of a face-centered vector field u
   which lives on an AMR Hierarchy. \\
   u-= Grad(Lapl^-1 (Div(u)))
 */
class EBCompositeMACProjector
{
public:

  ///
  ~EBCompositeMACProjector();

  ///
  /**
     a_eblg:               AMR hierarchy of grids \\
     a_refRat:            Refinement ratios between levels. a_refRat[i] is between levels i and i+1 \\
     a_coarsestDomain:    Computational domain at level 0\\
     a_coarsestDx:        Grid spacing at level 0 \\
     a_origin:            Physical location of the lower corner of the domain \\
     a_baseDomainBCVel :  Boundary conditions for velocity \\
     a_baseDomainBCPhi :  Boundary conditions of phi (for Lapl(phi) = div(u) solve \\
     a_subtractOffMean :  Set this to be true if you want the mean of m_phi = zero \\
     a_numLevels       :  If data is defined on a set of levels less than the vector lengths, this is the number of defined levels.
     **********************This must be less than or equal to vector length.  Set to -1 if you want numlevels = vector length.***\\
     a_verbosity       :  3 is the normal amount of verbosity.  Salt to taste.     \\
     a_preCondIters    :  number of iterations to do for pre-conditioning \\
     a_time            :  time for boundary conditions \\
     a_relaxType       :  0 for Jacobi, 1 for gs, 2 for line relax.     \\
     a_bottomSolverType:  0 for BiCGStab, 1 for EBSimpleSolver     \\
     The embedded boundary's boundary conditions are always no-flow.
  */
  EBCompositeMACProjector(const Vector<EBLevelGrid>&                      a_eblg,
                          const Vector<int>&                              a_refRat,
                          const Vector<RefCountedPtr<EBQuadCFInterp> >&   a_quadCFI,
                          const RealVect&                                 a_coarsestDx,
                          const RealVect&                                 a_origin,
                          const RefCountedPtr<BaseDomainBCFactory>&       a_baseDomainBCVel,
                          const RefCountedPtr<BaseDomainBCFactory>&       a_baseDomainBCPhi,
                          const RefCountedPtr<BaseEBBCFactory>&           a_ebbcPhi,
                          const bool&                                     a_subtractOffMean,
                          const int&                                      a_numLevels,
                          const int&                                      a_verbosity,
                          const int&                                      a_numPreCondIters,
                          const Real&                                     a_time,
                          const int&                                      a_relaxType,
                          const int&                                      a_bottomSolverType);

  ///
  /**
     velocity--input and output as divergence-free
     gradient--output-pure gradient component of input velocity.
     velocity-= Grad(Lapl^-1 (Div(velocity)))
     gradient = Grad(Lapl^-1 (Div(velocity))).
     Velocity must be single-component.
     If the boundary velocity is NULL,
     returns m_exitStatus of AMRMultiGrid.H
   */
  int
  project(Vector<LevelData<EBFluxFAB>* >&              a_velocity,
          Vector<LevelData<EBFluxFAB>* >&              a_gradient,
          const Real&                                  a_gradCoef=1.0,
          const Real&                                  a_divuCoef=1.0,
          const Vector<LevelData<BaseIVFAB<Real> >* >* a_boundaryVelo=NULL);

  ///
  void
  setSolverParams(int  a_numSmooth,
                  int  a_itermax,
                  int  a_mgcycle,
                  Real a_hang,
                  Real a_tolerance,
                  int  a_verbosity=3,
                  Real a_normThresh=1.e-12);

  ///
  void setTime(Real a_time);

  ///
  void enforceVelocityBC(Vector<LevelData<EBFluxFAB>* >&  a_velo,
                         const bool&                      a_doDivFreeOutflow=false);

  ///
  void enforceVelocityBC(Vector<LevelData<EBFluxFAB>* >&  a_velo,
                                const int&                a_comp,
                                const bool&               a_doDivFreeOutflow=false);

  ///
  void enforceGradientBC(Vector<LevelData<EBFluxFAB>* >&       a_grad,
                         const Vector<LevelData<EBCellFAB>* >& a_phi);

  ///
  /**
     null boundary velocity means no embedded boundary contribution to divergence.
     This is always true in the context of the projection.  The boundary vel is there
     for testing purposes.
  */
  void kappaDivergence(Vector<LevelData<EBCellFAB>* >&              a_divu,
                       Vector<LevelData<EBFluxFAB>* >&              a_velo,
                       const Vector<LevelData<BaseIVFAB<Real> >* >* a_boundaryVelo=NULL);

  ///
  void gradient(Vector<LevelData<EBFluxFAB>* >&       a_grad,
                const Vector<LevelData<EBCellFAB>* >& a_phi);

  ///
  /**
     This assumes that each face of the velocity fabs have the same component of the
     velocity and must be corrected by the graident.  The normal face velocity gets corrected by the
     gradient on its own faces.    The velocity on faces tangential to the velocity direction
     get corrected by an average of neighboring normal faces.  on faces that overlap the coarse-fine
     interface, we only use faces that are in the valid region of the level.
   */
  void correctVelocityComponent(Vector<LevelData<EBFluxFAB>* >& a_velocity,
                                Vector<LevelData<EBFluxFAB>* >& a_gradient,
                                int                             a_icomp);


  void correctTangentialVelocity(EBFaceFAB&        a_velocity,
                                 const EBFaceFAB&  a_gradient,
                                 const Box&        a_grid,
                                 const EBISBox&    a_ebisBox,
                                 const IntVectSet& a_cfivs);
  void
  correctVelocityComponent(BaseIVFAB<Real>       &  a_coveredVel,
                           const Vector<VolIndex>&  a_coveredFace,
                           const IntVectSet      &  a_coveredSets,
                           const EBFluxFAB       &  a_macGradient,
                           const EBISBox         &  a_ebisBox,
                           int    a_idir,  Side::LoHiSide a_sd,    int a_icomp);

  void
  correctVelocityComponent(Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >      &  a_coveredVelLo,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >      &  a_coveredVelHi,
                           const Vector< LayoutData< Vector< Vector<VolIndex> > >* >&  a_coveredFaceLo,
                           const Vector< LayoutData< Vector< Vector<VolIndex> > >* >&  a_coveredFaceHi,
                           const Vector< LayoutData< Vector< IntVectSet > >* >      &  a_coveredSetsLo,
                           const Vector< LayoutData< Vector< IntVectSet > >* >      &  a_coveredSetsHi,
                           const Vector<LevelData<EBFluxFAB>* >                     &  a_macGradient,
                           int                                                         a_faceDir,
                           int                                                         a_velComp);

  ///
  /**
     Provide an initial guess for phi with this function, otherwise phi=0 for init.
   */
  void
  setInitialPhi(Vector<LevelData<EBCellFAB>* >&  a_phi);

  AMRMultiGrid<LevelData<EBCellFAB> >& getSolver()
  {return m_solver;}
  Vector<LevelData<EBCellFAB>* >& getPhi()
  {return m_phi;}
protected:

  /*****/
  void
  initialize(Vector<LevelData<EBCellFAB>* >&  a_phi);


  //coarse fine interface with next COARSER level
  Vector<LevelData<EBCellFAB>* >        m_divu;
  Vector<LevelData<EBCellFAB>* >        m_phi;
  Vector<EBFluxRegister*>               m_fluxReg;
  AMRMultiGrid<LevelData<EBCellFAB> >   m_solver;
  BiCGStabSolver<LevelData<EBCellFAB> > m_bottomSolverBiCG;
  EBSimpleSolver                        m_bottomSolverSimp;
  GMRESSolver<LevelData<EBCellFAB> >    m_bottomSolverGMRes;
  EBAMRPoissonOpFactory*                m_opFactory;
  RefCountedPtr<BaseDomainBCFactory>    m_baseDomainBCFactVel;
  RefCountedPtr<BaseDomainBCFactory>    m_baseDomainBCFactPhi;
  RefCountedPtr<BaseEBBCFactory>        m_baseEBBCVel;

  Vector<EBLevelGrid>                      m_eblg;
  Vector<RefCountedPtr<EBQuadCFInterp> >   m_quadCFI;
  Vector<RealVect>                 m_dx;
  RealVect                         m_origin;
  Vector<int>                      m_refRat;
  int                              m_numLevels;
  int                              m_bottomSolverType;
  bool                             m_subtractOffMean;
  bool                             m_useInitialGuess;
  Real                             m_time;
private:

  EBCompositeMACProjector()
  {
    MayDay::Error("Weak construction is bad where there is this much memory management");
  }

  EBCompositeMACProjector(const EBCompositeMACProjector& a_input)
  {
    MayDay::Error("We disallow copy construction for objects with pointered data");
  }

  void operator=(const EBCompositeMACProjector& a_input)
  {
    MayDay::Error("We disallow assignment for objects with pointered data");
  }

};

/** Make reusable functions standalone
 */
/// utility function for tangential gradients
Real
getAverageFaceGrad(const VolIndex&   a_vof,
                   const EBFaceFAB&  a_macGradient,
                   const EBISBox&    a_ebisBox,
                   int               a_faceDir);

#include "NamespaceFooter.H"
#endif
