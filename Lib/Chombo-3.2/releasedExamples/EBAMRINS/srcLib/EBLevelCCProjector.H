#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLEVELCCPROJECTOR_H_
#define _EBLEVELCCPROJECTOR_H_
#include "EBLevelMACProjector.H"
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
#include "EBQuadCFInterp.H"

#include "NamespaceHeader.H"

///
/**
   Class to project a face-centered. velocity field on a level.  \\
   u -= G (DG^-1)(D u)
   This is done as interface to the level mac projector.  Velocity is averaged to faces
   then a mac projection is executed.  The gradient is then averaged back to cells
   and subtracted off the velocity.  So it really looks more like \\
   u -=     AveFtoC(G(DG^-1)(D(AveCToF(u))))
   Ghost cells over the coarse-fine interface are filled
   by quadradic coarse-fine interpolation.
   This class does not assume that the boundary conditions at the embedded boundary are no-flow.
 */
class EBLevelCCProjector
{
public:

  ///
  ~EBLevelCCProjector();

  ///
  /**
     a_grids:          boxes on the level                                      \\
     a_ebisl:          ebislayout at this level                                \\
     a_domain:         domain of compuation at this refinement                 \\
     a_dx:             grid spacing                                            \\
     a_origin:         physical location of the lower corner of the domain     \\\
     a_bottomSolver:   bottom solver for multigrid                             \\
     a_ebbcPhi:        boundary conditions for phi at embedded boundary        \\
     a_refToCoar:      refinement ratio to next coarser AMR level              \\
     a_hasCoaser       true if there is a coarser AMR level                    \\
     a_domainbcPhi:    boundary conditons for phi at domain boundary           \\
     a_domainbcVel:    boundary conditons for velocity at domain boundary      \\
     a_numSmooths:     number of smoothings that multigrid uses                \\
     a_mgCycle:        1 for v cycle, 2 for w cycle                            \\
     a_maxIterations:  max number of multigrid iterations                      \\
     a_tolerance:      factor to reduce residual by in solve                   \\
     a_maxDepth:       maximum multigrid depth.  -1 for all the way down.      \\
     a_time:           time for boundary conditions                            \\
     a_doEBCFCrossing: true to enable EB Coarse-Fine interface code.           \\
  */
  EBLevelCCProjector(const DisjointBoxLayout &                        a_grids,
                     const DisjointBoxLayout &                        a_gridsCoar,
                     const EBISLayout &                               a_ebisl,
                     const EBISLayout &                               a_ebislCoar,
                     const ProblemDomain &                            a_domain,
                     const RealVect &                                 a_dx,
                     const RealVect &                                 a_origin,
                     const int &                                      a_refToCoar,
                     const bool &                                     a_hasCoarser,
                     const LinearSolver<LevelData<EBCellFAB> > &      a_bottomSolver,
                     const RefCountedPtr<BaseEBBCFactory> &           a_ebbcPhi,
                     const RefCountedPtr<BaseDomainBCFactory> &       a_domainbcPhi,
                     const RefCountedPtr<BaseDomainBCFactory> &       a_domainbcVel,
                     const int   &                                    a_numSmooths,
                     const int   &                                    a_mgCycle,
                     const int   &                                    a_maxIterations,
                     const Real  &                                    a_tolerance,
                     const int   &                                    a_maxDepth,
                     const Real&                                      a_time,
                     const bool&                                      a_doEBCFCrossing,
                     const IntVect&                                   a_nghostPhi,
                     const IntVect&                                   a_nghostRhs)  ;


  ///
  /**
     velocity--input and output as divergence-free
     gradient--pure gradient component of input velocity.
     Coarse-Fine Boundary conditions for the velocity are interpolated
     from the coarser level.  For no-flow through eb, leave boundaryvel = NULL.
   */
  void
  project(LevelData<EBCellFAB>&         a_velocity,
          LevelData<EBCellFAB>&         a_gradient,
          const LevelData<EBCellFAB>&   a_velCoar,
          const LevelData<BaseIVFAB<Real> >* const a_boundaryVel = NULL);

  void
  kappaDivergence(LevelData<EBCellFAB>&         a_divergence,
                  LevelData<EBCellFAB>&         a_velocity,
                  const LevelData<EBCellFAB>&   a_velCoar,
                  const LevelData<BaseIVFAB<Real> >* const a_boundaryVel = NULL);

protected:

  DisjointBoxLayout                  m_grids;
  EBISLayout                         m_ebisl;
  DisjointBoxLayout                  m_gridsCoar;
  EBISLayout                         m_ebislCoar;
  ProblemDomain                      m_domain;
  RealVect                           m_dx;
  EBLevelMACProjector*               m_macProjector;
  EBQuadCFInterp                     m_patcher;
  bool                               m_doEBCFCrossing;
  int                                m_refToCoar;
  bool                               m_hasCoarser;
  IntVect                            m_nghostPhi;
  IntVect                            m_nghostRhs;
  LayoutData<IntVectSet>             m_cfivs;

private:

  EBLevelCCProjector()
  {
    MayDay::Error("Weak construction is bad");
  }
  EBLevelCCProjector(const EBLevelCCProjector& a_input)
  {
    MayDay::Error("We disallow copy construction for objects with pointered data");
  }
  void operator=(const EBLevelCCProjector& a_input)
  {
    MayDay::Error("We disallow assignment for objects with pointered data");
  }

};

//reusable functions made external
extern void
ccpAverageVelocityToFaces(LevelData<EBFluxFAB>&         a_macVeloc,
                          const LevelData<EBCellFAB>&   a_cellVeloc,
                          const DisjointBoxLayout&      a_grids,
                          const EBISLayout&             a_ebisl,
                          const ProblemDomain&          a_domain,
                          const RealVect&               a_dx,
                          const LayoutData<IntVectSet>& a_cfivs);

extern void
ccpAverageVelocityToFaces(EBFaceFAB&             a_faceVel,
                          const EBCellFAB&       a_cellVel,
                          const EBGraph&         a_ebGraph,
                          const Box&             a_grid,
                          const int&             a_idir,
                          const ProblemDomain&   a_domain,
                          const RealVect&        a_dx);

extern void
ccpAverageVelocityToFaces(LevelData<EBFluxFAB>&         a_macVeloc,
                          const LevelData<EBCellFAB>&   a_cellVeloc,
                          const DisjointBoxLayout&      a_grids,
                          const EBISLayout&             a_ebisl,
                          const ProblemDomain&          a_domain,
                          const RealVect&               a_dx,
                          const LayoutData<IntVectSet>& a_cfivs,
                          const int&                    a_comp);

extern void
ccpAverageVelocityToFaces(EBFaceFAB&             a_faceVel,
                          const EBCellFAB&       a_cellVel,
                          const EBGraph&         a_ebGraph,
                          const Box&             a_grid,
                          const int&             a_idir,
                          const ProblemDomain&   a_domain,
                          const RealVect&        a_dx,
                          const int&             a_comp);

extern void
ccpAverageStressToFaces(LevelData<EBFluxFAB>&         a_macVeloc,
                          const LevelData<EBCellFAB>&   a_cellVeloc,
                          const DisjointBoxLayout&      a_grids,
                          const EBISLayout&             a_ebisl,
                          const ProblemDomain&          a_domain,
                          const RealVect&               a_dx,
                          const LayoutData<IntVectSet>& a_cfivs);

extern void
ccpAverageStressToFaces(EBFaceFAB&             a_faceVel,
                        const EBCellFAB&       a_cellVel,
                        const EBGraph&         a_ebGraph,
                        const Box&             a_grid,
                        const int&             a_idir,
                        const ProblemDomain&   a_domain,
                        const RealVect&        a_dx);

extern void
ccpAverageFaceToCells(LevelData<EBCellFAB>&         a_cellData,
                      const LevelData<EBFluxFAB>&   a_macData,
                      const DisjointBoxLayout&      a_grids,
                      const EBISLayout&             a_ebisl,
                      const ProblemDomain&          a_domain,
                      const RealVect&               a_dx);

extern void
ccpAverageFaceToCells(EBCellFAB&             a_cellData,
                      const EBFluxFAB&       a_fluxData,
                      const EBGraph&         a_ebGraph,
                      const Box&             a_grid,
                      const ProblemDomain&   a_domain,
                      const RealVect&        a_dx);

extern void
ccpAverageFaceToCellsScalar(LevelData<EBCellFAB>&         a_cellData,
                            const LevelData<EBFluxFAB>&   a_macData,
                            const DisjointBoxLayout&      a_grids,
                            const EBISLayout&             a_ebisl,
                            const ProblemDomain&          a_domain,
                            const RealVect&               a_dx,
                            const int &                   a_dir);

extern void
ccpAverageFaceToCellsScalar(EBCellFAB&             a_cellData,
                            const EBFluxFAB&       a_fluxData,
                            const EBGraph&         a_ebGraph,
                            const Box&             a_grid,
                            const ProblemDomain&   a_domain,
                            const RealVect&        a_dx,
                            const int &            a_dir);

extern void
ccpCellGradientFromFaceData(LevelData<EBCellFAB>&         a_cellData,
                            const LevelData<EBFluxFAB>&   a_macData,
                            const DisjointBoxLayout&      a_grids,
                            const EBISLayout&             a_ebisl,
                            const ProblemDomain&          a_domain,
                            const RealVect&               a_dx,
                            const int &                   a_dir);

extern void
ccpCellGradientFromFaceData(EBCellFAB&             a_cellData,
                            const EBFluxFAB&       a_fluxData,
                            const EBGraph&         a_ebGraph,
                            const Box&             a_grid,
                            const ProblemDomain&   a_domain,
                            const RealVect&        a_dx,
                            const int &            a_dir);

extern void
ccpExtrapolateToDomainBoundaries(LevelData<EBFluxFAB>&    a_macData,
                                 const DisjointBoxLayout& a_grids,
                                 const EBISLayout&        a_ebisl,
                                 const ProblemDomain&     a_domain,
                                 const RealVect&          a_dx);

extern void
ccpExtrapolateToDomainBoundaries(EBFaceFAB&             a_faceData,
                                 const EBGraph&         a_ebGraph,
                                 const Box&             a_grid,
                                 const int&             a_idir,
                                 const ProblemDomain&   a_domain,
                                 const RealVect&        a_dx);

extern void
ccpLinearInterp(Real &                    a_dataOnLine,
                const Vector<RealVect>&   a_faceLoc,
                const RealVect&           a_intersectLoc,
                const Vector<Real>&       a_interpolationData,
                const int&                a_planeDir);

extern void
ccpBilinearInterp(Real &                             a_dataOnLine,
                  const Vector<RealVect>&            a_faceLoc,
                  const RealVect&                    a_intersectLoc,
                  const Vector<Real>&                a_interpolationData,
                  const Tuple<int, SpaceDim-1>&      a_planeDir);


extern void
ccpExtrapFaceToCovered(bool&                   a_dropOrder,
                           Real&                   a_extrapVal,
                           const EBFaceFAB&        a_primFace,
                           const VolIndex&         a_vof,
                           const int&              a_faceDir,
                           const Side::LoHiSide&   a_sd,
                           const RealVect&         a_normal,
                           const RealVect&         a_dx,
                           const int&              a_icomp);

///
/**
   Multi-D extrapolation to covered face
*/
extern void
ccpJohansenExtrapFaceToCovered(bool&                   a_dropOrder,
                               Real&                   a_extrapVal,
                               const EBFaceFAB&        a_primFace,
                               const EBISBox&          a_ebisBox,
                               const VolIndex&         a_vof,
                               const int&              a_faceDir,
                               const Side::LoHiSide&   a_sd,
                               const RealVect&         a_normal,
                               const RealVect&         a_dx,
                               const int&              a_icomp) ;

///
/**
   do one dimensional extrapolation to covered face
   preferring face direction
*/
extern Real
ccpOneDCoveredExtrapValue(const VolIndex&         a_vof,
                          const int&              a_idir,
                          const Side::LoHiSide    a_side,
                          const EBFaceFAB &       a_faceGrad,
                          const EBISBox &         a_ebisBox,
                          const Box &             a_grid,
                          const ProblemDomain &   a_domain,
                          const RealVect &        a_dx,
                          const int&              a_icomp);

///
/**
   do multi-D  extrapolation to covered face
   if possible.  if not, drop to 1d
   preferring face direction
*/
extern Real
ccpGetCoveredExtrapValue(const VolIndex&         a_vof,
                         const int&              a_idir,
                         const Side::LoHiSide    a_side,
                         const EBFaceFAB &       a_faceGrad,
                         const EBISBox &         a_ebisBox,
                         const Box &             a_grid,
                         const ProblemDomain &   a_domain,
                         const RealVect &        a_dx,
                         const int&              a_icomp);

#include "NamespaceFooter.H"
#endif
