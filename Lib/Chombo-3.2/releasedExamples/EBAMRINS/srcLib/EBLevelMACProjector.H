#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLEVELMACPROJECTOR_H_
#define _EBLEVELMACPROJECTOR_H_
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "Box.H"
#include "REAL.H"
#include "LevelData.H"
#include "EBFluxFAB.H"
#include "EBCellFAB.H"
#include "EBSimpleSolver.H"
#include "EBAMRPoissonOp.H"
#include "EBAMRPoissonOpFactory.H"
#include "MultiGrid.H"
#include "BiCGStabSolver.H"

#include "NamespaceHeader.H"

///
/**
   Class to project a face-centered. velocity field on a level.
   u -= G (DG^-1)(D u)
   This class does not assume that the boundary conditions at the embedded boundary are no-flow.
   Faces over coarse-fine interfaces are not treated any differently than interior faces.
 */
class EBLevelMACProjector
{
public:

  ///
  ~EBLevelMACProjector();

  ///
  /**
     a_grids:          boxes on the level                                      \\
     a_dx:             grid spacing                                            \\
     a_origin:         physical location of the lower corner of the domain     \\
     a_bottomSolver:   bottom solver for multigrid                             \\
     a_ebbcPhi:        boundary conditions for phi at embedded boundary        \\
     a_ebbcVel:        boundary conditions for velocity at embedded boundary   \\
     a_domainbcPhi:    boundary conditons for phi at domain boundary           \\
     a_domainbcVel:    boundary conditons for velocity at domain boundary      \\
     a_numSmooths:     number of smoothings that multigrid uses                \\
     a_mgCycle:        1 for v cycle, 2 for w cycle                            \\
     a_maxIterations:  max number of multigrid iterations                      \\
     a_tolerance:      factor to reduce residual by in solve                   \\
     a_maxDepth:       maximum multigrid depth.  -1 for all the way down.      \\
     a_time:           time for boundary conditions                            \\
   */
  EBLevelMACProjector(const DisjointBoxLayout&                        a_grids,
                      const EBISLayout&                               a_ebisl,
                      const ProblemDomain&                            a_domain,
                      const RealVect&                                 a_dx,
                      const RealVect&                                 a_origin,
                      const LinearSolver<LevelData<EBCellFAB> >&      a_bottomSolver,
                      const RefCountedPtr<BaseEBBCFactory>&           a_ebbcPhi,
                      const RefCountedPtr<BaseDomainBCFactory>&       a_domainbcPhi,
                      const RefCountedPtr<BaseDomainBCFactory>&       a_domainbcVel,
                      const int  &                                    a_numSmooths,
                      const int  &                                    a_mgCycle,
                      const int  &                                    a_maxIterations,
                      const Real &                                    a_tolerance,
                      const int  &                                    a_maxDepth,
                      const Real&                                     a_time,
                      const IntVect&                                  a_nghostPhi,
                      const IntVect&                                  a_nghostRhs);


  ///
  /**
     velocity--input and output as divergence-free
     gradient--pure gradient component of input velocity.
     Setting boundary velocity to null makes it zero.
  */
  void
  project(LevelData<EBFluxFAB>&  a_velocity,
          LevelData<EBFluxFAB>&  a_gradient,
          const LevelData< BaseIVFAB<Real> >* const a_boundaryVelocity = NULL);

  static void setVerbose(bool a_verbose);
  static void setCurLevel(int a_curLevel);
  static void setDebugString(string a_debugString);
  static string getDebugString();
  static int getCurLevel();
  static bool getVerbose();

protected:

  void
  solve(LevelData<EBCellFAB>&       a_phi,
        const LevelData<EBCellFAB>& a_rhs);

  int m_maxIterations;
  Real m_tolerance;
  DisjointBoxLayout                  m_grids;
  EBISLayout                         m_ebisl;
  ProblemDomain                      m_domain;
  RealVect                           m_dx;
  RealVect                           m_origin;
  EBAMRPoissonOpFactory*             m_levelOpFactory;
  RefCountedPtr<BaseEBBCFactory>     m_ebbcPhi;
  RefCountedPtr<BaseDomainBCFactory> m_domainBCFactPhi;
  RefCountedPtr<BaseDomainBCFactory> m_domainBCFactVel;

  IntVect m_nghostPhi;
  IntVect m_nghostRhs;

  MultiGrid<LevelData<EBCellFAB> >   m_solver;

  //stuff for debugging output
  static int  s_curLevel;
  static bool s_verbose;
  static string s_debugString;
private:
  EBLevelMACProjector(const EBLevelMACProjector& a_input)
  {
    MayDay::Error("We disallow copy construction for objects with pointered data");
  }

  void operator=(const EBLevelMACProjector& a_input)
  {
    MayDay::Error("We disallow assignment for objects with pointered data");
  }

  EBLevelMACProjector()
  {
    MayDay::Error("weak construction is bad");
  }
};

//make reusable functions standalone
extern void
macGradient(LevelData<EBFluxFAB>&        a_gradPhi,
            const LevelData<EBCellFAB>&  a_phi,
            const DisjointBoxLayout&     a_dbl,
            const EBISLayout&            a_ebisl,
            const ProblemDomain&         a_domain,
            const RealVect&              a_dx);

extern void
macGradient(EBFluxFAB&           a_gradPhi,
            const EBCellFAB&     a_phi,
            const EBISBox&       a_ebisBox,
            const Box&           a_box,
            const ProblemDomain& a_domain,
            const RealVect&      a_dx);

///NULL boundary velocity signals no contribution from boundary face to divergence
extern void
macKappaDivergence(LevelData<EBCellFAB>&              a_divVel,
                   const LevelData<EBFluxFAB>&        a_velocity,
                   const DisjointBoxLayout&           a_dbl,
                   const EBISLayout&                  a_ebisl,
                   const ProblemDomain&               a_domain,
                   const RealVect&                    a_dx,
                   const LevelData<BaseIVFAB<Real> >* a_boundaryVelo = NULL);

///NULL boundary velocity signals no contribution from boundary face to divergence
extern void
macKappaDivergence(EBCellFAB&             a_divVel,
                   const EBFluxFAB&       a_velocity,
                   const EBISBox&         a_ebisBox,
                   const Box&             a_box,
                   const ProblemDomain&   a_domain,
                   const RealVect&        a_dx,
                   const BaseIVFAB<Real>* a_boundaryVel = NULL);

extern void
macEnforceVelocityBC(LevelData<EBFluxFAB>&              a_velocity,
                     const DisjointBoxLayout&           a_grids,
                     const EBISLayout&                  a_ebisl,
                     const ProblemDomain&               a_domain,
                     const RealVect&                    a_dx,
                     const RealVect&                    a_origin,
                     RefCountedPtr<BaseDomainBCFactory> a_domainbcVel,
                     const Real&                        a_time,
                     const bool&                        a_doDivFreeOutflow,
                     const LevelData<BaseIVFAB<Real> >* const a_boundaryVelocity);

extern void
macEnforceVelocityBC(LevelData<EBFluxFAB>&              a_velocity,
                     const DisjointBoxLayout&           a_grids,
                     const EBISLayout&                  a_ebisl,
                     const ProblemDomain&               a_domain,
                     const RealVect&                    a_dx,
                     const RealVect&                    a_origin,
                     RefCountedPtr<BaseDomainBCFactory> a_domainbcVel,
                     const Real&                        a_time,
                     const bool&                        a_doDivFreeOutflow,
                     const int&                         a_comp,
                     const LevelData<BaseIVFAB<Real> >* const a_boundaryVelocity);

extern void
macEnforceGradientBC(LevelData<EBFluxFAB>&              a_gradPhi,
                     const LevelData<EBCellFAB>&        a_phi,
                     const DisjointBoxLayout&           a_grids,
                     const EBISLayout&                  a_ebisl,
                     const ProblemDomain&               a_domain,
                     const RealVect&                    a_dx,
                     const RealVect&                    a_origin,
                     const Real&                        a_time,
                     RefCountedPtr<BaseDomainBCFactory> a_domainBCFactPhi);

#include "NamespaceFooter.H"
#endif
