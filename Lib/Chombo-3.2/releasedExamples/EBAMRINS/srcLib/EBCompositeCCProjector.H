#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCOMPOSITECCPROJECTOR_H_
#define _EBCOMPOSITECCPROJECTOR_H_

#include "EBLevelCCProjector.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "Box.H"
#include "REAL.H"
#include "LevelData.H"
#include "EBFluxFAB.H"
#include "EBCellFAB.H"
#include "EBSimpleSolver.H"
#include "EBAMRPoissonOp.H"
#include "MultiGrid.H"
#include "LevelFluxRegister.H"
#include "EBQuadCFInterp.H"
#include "EBCompositeMACProjector.H"

#include "NamespaceHeader.H"

///
/**
   Projection to take out pure gradient component of a face-centered vector field u
   which lives on an AMR Hierarchy. \\
   u-= Grad(Lapl^-1 (Div(u)))
   This is done as interface to the composite mac projector.  Velocity is averaged to faces
   then a mac projection is executed.  The gradient is then averaged back to cells
   and subtracted off the velocity.  So it really looks more like \\
   u -=     AveFtoC(G(DG^-1)(D(AveCToF(u))))
 */
class EBCompositeCCProjector
{
public:

  ///
  ~EBCompositeCCProjector();

  ///
  /**
     a_dbl:               AMR hierarchy of grids \\
     a_ebisl:             Hierarchy of EBISLayouts                  \\
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
     a_relaxType       :  0 means point Jacobi, 1 is Gauss-Seidel. \\
     a_bottomSolverType:  0 for BiCGStab, 1 for EBSimpleSolver     \\
     The embedded boundary's boundary conditions are always no-flow. \\
     Will define a new mac projector  if you do not send in one.
  */
  EBCompositeCCProjector(const Vector<EBLevelGrid>&                      a_dbl,
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
                         const int&                                      a_bottomSolverType,
                         EBCompositeMACProjector*                        a_inputMAC= NULL);

  ///
  /**
     velocity--input and output as divergence-free
     gradient--output-pure gradient component of input velocity.
     velocity-= Grad(Lapl^-1 (Div(velocity)))
     gradient = Grad(Lapl^-1 (Div(velocity))).
     Velocity must be have spacedim components.
     returns m_exitStatus of AMRMultiGrid.H
   */
  int
  project(Vector<LevelData<EBCellFAB>* >&               a_velocity,
          Vector<LevelData<EBCellFAB>* >&               a_gradient,
          const Real&                                   a_gradCoef=1.0,
          const Real&                                   a_divuCoef=1.0,
          const Vector<LevelData<BaseIVFAB<Real> >* >*  a_boundaryVelo=NULL);

  ///
  Vector<LevelData<EBFluxFAB>* >& getMacVelocity()
  {
    return m_macVelo;
  }

  ///
  const Vector<LevelData<EBFluxFAB>* >& getMacVelocity() const
  {
    return m_macVelo;
  }

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

  void
  gradient(Vector<LevelData<EBCellFAB>* >&  a_gradient,
           Vector<LevelData<EBCellFAB>* >&  a_phi);

  /// boundary velo = NULL means boundary velocity set to zero for divergence
  void kappaDivergence(Vector<LevelData<EBCellFAB>* >&              a_divu,
                       Vector<LevelData<EBCellFAB>* >&              a_velo,
                       const Vector<LevelData<BaseIVFAB<Real> >* >* a_boundaryVelo = NULL);

  void
  averageVelocityToFaces(Vector<LevelData<EBFluxFAB>* >&        a_macVelo,
                         Vector<LevelData<EBCellFAB>* >&        a_velocity);

  void
  averageVelocityToFaces(Vector<LevelData<EBFluxFAB>* >&        a_macVelo,
                         Vector<LevelData<EBCellFAB>* >&        a_velocity,
                         const int &                            a_comp);

  void
  averageStressToFaces(Vector<LevelData<EBFluxFAB>* >&        a_macVelo,
                       Vector<LevelData<EBCellFAB>* >&        a_velocity);

  void
  averageFaceToCells(Vector<LevelData<EBCellFAB>* >&         a_cellData,
                     const Vector<LevelData<EBFluxFAB>* >&   a_faceData);

  Vector<LevelData<EBCellFAB>* >& getPhi()
  {return m_macProjector->getPhi();}

  ///
  /**
     Provide an initial guess for phi with this function, otherwise phi=0 for init.
   */
  void
  setInitialPhi(Vector<LevelData<EBCellFAB>* >&  a_phi)
  {m_macProjector->setInitialPhi(a_phi);}

  protected:

  //mac projection objects
  EBCompositeMACProjector*         m_macProjector;
  //whether mac projection was sent in or not
  bool                             m_externMAC;
  bool                             m_subtractOffMean;
  Vector<LevelData<EBFluxFAB>* >   m_macVelo;
  Vector<LevelData<EBFluxFAB>* >   m_macGrad;

  //grid information
  Vector<EBLevelGrid>                      m_eblg;
  Vector<RefCountedPtr<EBQuadCFInterp> >   m_quadCFI;
  Vector<RefCountedPtr<EBPWLFillPatch> >   m_pwlCFI;
  Vector<RealVect>                         m_dx;
  Vector<int>                              m_refRat;
  int                                      m_numLevels;

private:

  EBCompositeCCProjector()
  {
    MayDay::Error("Weak construction is bad where there is this much memory management");
  }

  EBCompositeCCProjector(const EBCompositeCCProjector& a_input)
  {
    MayDay::Error("We disallow copy construction for objects with pointered data");
  }

  void operator=(const EBCompositeCCProjector& a_input)
  {
    MayDay::Error("We disallow assignment for objects with pointered data");
  }

};

#include "NamespaceFooter.H"
#endif
