#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CCProjectorF_F.H"
#include "EBLevelCCProjector.H"
#include "NeumannPoissonEBBC.H"
#include "NeumannPoissonDomainBC.H"
#include "EBFluxFactory.H"
#include "EBAMRPoissonOp.H"
#include "EBArith.H"
#include "FaceIterator.H"

#include "NamespaceHeader.H"

/*****/
EBLevelCCProjector::
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
                   const IntVect&                                   a_nghostRhs)
{
  m_grids          = a_grids;
  m_gridsCoar      = a_gridsCoar;
  m_ebislCoar      = a_ebislCoar;
  m_ebisl          = a_ebisl;
  m_domain         = a_domain;
  m_dx             = a_dx;
  m_hasCoarser     = a_hasCoarser;
  m_refToCoar      = a_refToCoar;
  m_nghostPhi      = a_nghostPhi;
  m_nghostRhs      = a_nghostRhs;
  m_doEBCFCrossing = a_doEBCFCrossing;

  //define mac projector here.
  m_macProjector = new EBLevelMACProjector(m_grids, m_ebisl, m_domain, m_dx, a_origin, a_bottomSolver,
                                           a_ebbcPhi, a_domainbcPhi, a_domainbcVel,
                                           a_numSmooths, a_mgCycle, a_maxIterations, a_tolerance,
                                           a_maxDepth, a_time,
                                           m_nghostPhi, m_nghostRhs);

  if (m_hasCoarser)
    {
      ProblemDomain domainCrse = coarsen(m_domain, m_refToCoar);
      //patcher is defined with SpaceDim variables
      EBArith::defineCFIVS(m_cfivs, m_grids, m_domain);

      m_patcher.define(m_grids, m_gridsCoar,
                       m_ebisl, m_ebislCoar,
                       domainCrse, a_refToCoar, SpaceDim, m_cfivs);

    }
}
/*****/
void
EBLevelCCProjector::
project(LevelData<EBCellFAB> &         a_velocity,
        LevelData<EBCellFAB> &         a_gradient,
        const LevelData<EBCellFAB> &   a_velCoar,
        const LevelData<BaseIVFAB<Real> >* const a_boundaryVel)
{
  //interpolate at coarse-fine interface
  Interval interv(0, SpaceDim-1);
  if (m_hasCoarser)
    {
      m_patcher.interpolate(a_velocity,
                            a_velCoar,
                            interv);
    }
  a_velocity.exchange(interv);

  EBFluxFactory fluxfact(m_ebisl);
  LevelData<EBFluxFAB>  macVel(m_grids, 1, IntVect::Unit, fluxfact);
  LevelData<EBFluxFAB> macGrad(m_grids, 1, IntVect::Zero, fluxfact);

  //average velocity to faces
  ccpAverageVelocityToFaces(macVel, a_velocity, m_grids, m_ebisl, m_domain, m_dx,
                            m_cfivs);

  //use macprojector on the result
  m_macProjector->project(macVel, macGrad, a_boundaryVel);

  //extrapolate gradient to domain faces so that
  //we get reasonable answers near domain boundary
  ccpExtrapolateToDomainBoundaries(macGrad, m_grids, m_ebisl, m_domain, m_dx);

  //put the gradient into cell centers
  ccpAverageFaceToCells(a_gradient, macGrad, m_grids, m_ebisl, m_domain, m_dx);

  //subtract cell-centered gradient off input vel
  EBLevelDataOps::incr(a_velocity, a_gradient, -1.0);
}
void
EBLevelCCProjector::
kappaDivergence(LevelData<EBCellFAB> &         a_divergence,
                LevelData<EBCellFAB> &         a_velocity,
                const LevelData<EBCellFAB> &   a_velCoarOld,
                const LevelData<BaseIVFAB<Real> >* const a_boundaryVel)
{
  //interpolate at coarse-fine interface
  Interval interv(0, SpaceDim-1);
  if (m_hasCoarser)
    {
      m_patcher.interpolate(a_velocity,
                            a_velCoarOld,
                            interv);
    }
  a_velocity.exchange(interv);

  EBFluxFactory fluxfact(m_ebisl);
  LevelData<EBFluxFAB>  macVel(m_grids, 1, IntVect::Unit, fluxfact);

  //average velocity to faces
  ccpAverageVelocityToFaces(macVel, a_velocity, m_grids, m_ebisl, m_domain, m_dx,
                            m_cfivs);

  //use macprojector on the result
  macKappaDivergence(a_divergence, macVel, m_grids, m_ebisl, m_domain, m_dx, a_boundaryVel);
}
/*****/
EBLevelCCProjector::
~EBLevelCCProjector()
{
  delete m_macProjector;
}
/*****/
void
ccpExtrapolateToDomainBoundaries(LevelData<EBFluxFAB> &        a_macData,
                                 const DisjointBoxLayout &     a_grids,
                                 const EBISLayout &            a_ebisl,
                                 const ProblemDomain &         a_domain,
                                 const RealVect &              a_dx)
{

  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      const Box grid = a_grids.get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (!a_domain.isPeriodic(idir))
            {
              EBFaceFAB &       faceGrad = a_macData[dit()][idir];
              ccpExtrapolateToDomainBoundaries(faceGrad, ebgraph, grid, idir, a_domain, a_dx);
            }
        }
    }
}
/*****/
void
ccpExtrapolateToDomainBoundaries(EBFaceFAB &             a_faceData,
                                 const EBGraph &         a_ebGraph,
                                 const Box &             a_grid,
                                 const int &             a_idir,
                                 const ProblemDomain &   a_domain,
                                 const RealVect &        a_dx)
{
  CH_assert(!a_domain.isPeriodic(a_idir));

  //iterate over boundary faces
  for (SideIterator sit; sit.ok(); ++sit)
    {
      IntVect ivSideGrid =   a_grid.sideEnd(sit());
      IntVect ivSideDom  = a_domain.domainBox().sideEnd(sit());

      //check to see if we actually are at domain boundary
      if (ivSideGrid[a_idir] == ivSideDom[a_idir])
        {
          //create box and ivs adjacent to domain boundary
          //and interior to domain
          Box sideBox = adjCellBox(a_grid, a_idir, sit(), 1);
          int ishift = -sign(sit());
          sideBox.shift(a_idir, ishift);
          IntVectSet sideIVS(sideBox);
          FaceStop::WhichFaces stopCritSide;
          if (sit() == Side::Lo)
            {
              stopCritSide = FaceStop::LoBoundaryOnly;
            }
          else
            {
              stopCritSide = FaceStop::HiBoundaryOnly;
            }

          //do simple extrapolation to boundary faces from neighboring faces in same direction
          for (FaceIterator faceit(sideIVS, a_ebGraph, a_idir, stopCritSide);  faceit.ok(); ++faceit)
            {
              const FaceIndex & bndryFace = faceit();
              Real extrapValue = EBArith::extrapFaceGradToOutflow(bndryFace,sit(),a_idir,a_ebGraph,a_faceData,0);
              a_faceData(bndryFace,0) = extrapValue;
            }
        }
    }
}
/*****/
void
ccpAverageVelocityToFaces(EBFaceFAB &             a_faceVel,
                          const EBCellFAB &       a_cellVel,
                          const EBGraph &         a_ebGraph,
                          const Box &             a_grid,
                          const int &             a_idir,
                          const ProblemDomain &   a_domain,
                          const RealVect &        a_dx)
{
  CH_TIME("EBLevelCCProjector::ccpAverageVelocityToFaces");
  FaceStop::WhichFaces stopCrit;
  if (a_domain.isPeriodic(a_idir))
    {
      stopCrit = FaceStop::SurroundingWithBoundary;
    }
  else
    {
      stopCrit = FaceStop::SurroundingNoBoundary;
    }
  //initially set to zero so the boundary faces will all be zero
  a_faceVel.setVal(0.);

  BaseFab<Real> &       regFaceVel = a_faceVel.getSingleValuedFAB();
  const BaseFab<Real> & regCellVel = a_cellVel.getSingleValuedFAB();

  //only want non-domain boundary faces
  Box faceBox = a_grid;
  faceBox.grow(a_idir, 1);
  faceBox  &= a_domain;
  faceBox.grow(a_idir, -1);
  faceBox.surroundingNodes(a_idir);
  //do regular cells in fortran
  FORT_CCPAVECELLTOFACE(CHF_FRA1(regFaceVel, 0),
                        CHF_CONST_FRA1(regCellVel, a_idir),
                        CHF_CONST_INT(a_idir),
                        CHF_BOX(faceBox));

  //fix up irregular cells.  Only multiValued cells should and their neighbors
  //need fixing
  IntVectSet ivsIrreg = a_ebGraph.getMultiCells(a_grid);
  ivsIrreg.grow(1);
  ivsIrreg  &= a_grid;
  for (FaceIterator faceit(ivsIrreg, a_ebGraph, a_idir, stopCrit); faceit.ok(); ++faceit)
    {
      const FaceIndex & face = faceit();

      a_faceVel(face, 0) = 0.5*(a_cellVel(face.getVoF(Side::Hi), a_idir) +
                                a_cellVel(face.getVoF(Side::Lo), a_idir));
    }

}
/*****/
void
ccpAverageVelocityToFaces(LevelData<EBFluxFAB> &         a_macVeloc,
                          const LevelData<EBCellFAB> &   a_cellVeloc,
                          const DisjointBoxLayout &      a_grids,
                          const EBISLayout &             a_ebisl,
                          const ProblemDomain &          a_domain,
                          const RealVect &               a_dx,
                          const LayoutData<IntVectSet>&  a_cfivs)
{
  CH_TIME("EBLevelCCProjector::ccpAverageVelocityToFaces(level)");
  int ibox = 0;
  Interval interv(0,SpaceDim-1);
  LevelData<EBCellFAB>& velo = (LevelData<EBCellFAB>&) a_cellVeloc;
  velo.exchange(interv);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      const Box grid = a_grids.get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB &       faceVel = a_macVeloc[dit()][idir];
          const EBCellFAB & cellVel = a_cellVeloc[dit()];

          ccpAverageVelocityToFaces(faceVel, cellVel, ebgraph, grid, idir, a_domain, a_dx);

          ibox++;
        }
    }
}
/*****/
void
ccpAverageVelocityToFaces(EBFaceFAB &             a_faceVel,
                          const EBCellFAB &       a_cellVel,
                          const EBGraph &         a_ebGraph,
                          const Box &             a_grid,
                          const int &             a_idir,
                          const ProblemDomain &   a_domain,
                          const RealVect &        a_dx,
                          const int &             a_comp)
{
  CH_TIME("EBLevelCCProjector::ccpAverageVelocityToFaces");
  FaceStop::WhichFaces stopCrit;
  if (a_domain.isPeriodic(a_idir))
    {
      stopCrit = FaceStop::SurroundingWithBoundary;
    }
  else
    {
      stopCrit = FaceStop::SurroundingNoBoundary;
    }
  //initially set to zero so the boundary faces will all be zero
  a_faceVel.setVal(0.);

  BaseFab<Real> &       regFaceVel = a_faceVel.getSingleValuedFAB();
  const BaseFab<Real> & regCellVel = a_cellVel.getSingleValuedFAB();

  //only want non-domain boundary faces
  Box faceBox = a_grid;
  faceBox.grow(a_idir, 1);
  faceBox  &= a_domain;
  faceBox.grow(a_idir, -1);
  faceBox.surroundingNodes(a_idir);
  //do regular cells in fortran
  FORT_CCPAVECELLTOFACE(CHF_FRA1(regFaceVel, 0),
                        CHF_CONST_FRA1(regCellVel, a_comp),
                        CHF_CONST_INT(a_idir),
                        CHF_BOX(faceBox));

  //fix up irregular cells.  Only multiValued cells should and their neighbors
  //need fixing
  IntVectSet ivsIrreg = a_ebGraph.getMultiCells(a_grid);
  ivsIrreg.grow(1);
  ivsIrreg  &= a_grid;
  for (FaceIterator faceit(ivsIrreg, a_ebGraph, a_idir, stopCrit); faceit.ok(); ++faceit)
    {
      const FaceIndex & face = faceit();

      a_faceVel(face, 0) = 0.5*(a_cellVel(face.getVoF(Side::Hi), a_comp) +
                                a_cellVel(face.getVoF(Side::Lo), a_comp));
    }

}
/*****/
void
ccpAverageVelocityToFaces(LevelData<EBFluxFAB> &         a_macVeloc,
                          const LevelData<EBCellFAB> &   a_cellVeloc,
                          const DisjointBoxLayout &      a_grids,
                          const EBISLayout &             a_ebisl,
                          const ProblemDomain &          a_domain,
                          const RealVect &               a_dx,
                          const LayoutData<IntVectSet>&  a_cfivs,
                          const int &                    a_comp)
{
  CH_TIME("EBLevelCCProjector::ccpAverageVelocityToFaces(level)");
  int ibox = 0;
  Interval interv(a_comp,a_comp);
  LevelData<EBCellFAB>& velo = (LevelData<EBCellFAB>&) a_cellVeloc;
  velo.exchange(interv);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      const Box grid = a_grids.get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB &       faceVel = a_macVeloc[dit()][idir];
          const EBCellFAB & cellVel = a_cellVeloc[dit()];

          ccpAverageVelocityToFaces(faceVel, cellVel, ebgraph, grid, idir, a_domain, a_dx, a_comp);

          ibox++;
        }
    }
}
/*****/
void
ccpAverageStressToFaces(EBFaceFAB &             a_faceVel,
                        const EBCellFAB &       a_cellVel,
                        const EBGraph &         a_ebGraph,
                        const Box &             a_grid,
                        const int &             a_idir,
                        const ProblemDomain &   a_domain,
                        const RealVect &        a_dx)
{
  CH_TIME("EBLevelCCProjector::ccpAverageVelocityToFaces");
  FaceStop::WhichFaces stopCrit;
  if (a_domain.isPeriodic(a_idir))
    {
      stopCrit = FaceStop::SurroundingWithBoundary;
    }
  else
    {
      stopCrit = FaceStop::SurroundingNoBoundary;
    }
  //initially set to zero so the boundary faces will all be zero
  a_faceVel.setVal(0.);

  BaseFab<Real> &       regFaceVel = a_faceVel.getSingleValuedFAB();
  const BaseFab<Real> & regCellVel = a_cellVel.getSingleValuedFAB();

  //only want non-domain boundary faces
  Box faceBox = a_grid;
  faceBox.grow(a_idir, 1);
  faceBox  &= a_domain;
  faceBox.grow(a_idir, -1);
  faceBox.surroundingNodes(a_idir);
  //do regular cells in fortran
  FORT_CCPAVECELLTOFACE(CHF_FRA1(regFaceVel, 0),
                        CHF_CONST_FRA1(regCellVel, 0),
                        CHF_CONST_INT(a_idir),
                        CHF_BOX(faceBox));

  //fix up irregular cells.  Only multiValued cells should and their neighbors
  //need fixing
  IntVectSet ivsIrreg = a_ebGraph.getMultiCells(a_grid);
  ivsIrreg.grow(1);
  ivsIrreg  &= a_grid;
  for (FaceIterator faceit(ivsIrreg, a_ebGraph, a_idir, stopCrit); faceit.ok(); ++faceit)
    {
      const FaceIndex & face = faceit();

      a_faceVel(face, 0) = 0.5*(a_cellVel(face.getVoF(Side::Hi), 0) +
                                a_cellVel(face.getVoF(Side::Lo), 0));
    }

}
/*****/
void
ccpAverageStressToFaces(LevelData<EBFluxFAB> &         a_macVeloc,
                        const LevelData<EBCellFAB> &   a_cellVeloc,
                        const DisjointBoxLayout &      a_grids,
                        const EBISLayout &             a_ebisl,
                        const ProblemDomain &          a_domain,
                        const RealVect &               a_dx,
                        const LayoutData<IntVectSet>&  a_cfivs)
{
  CH_TIME("EBLevelCCProjector::ccpAverageVelocityToFaces(level)");
  int ibox = 0;
  CH_assert(a_cellVeloc.nComp() == 1);
  Interval interv(0,0);
  LevelData<EBCellFAB>& velo = (LevelData<EBCellFAB>&) a_cellVeloc;
  velo.exchange(interv);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      const Box grid = a_grids.get(dit());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB &       faceVel = a_macVeloc[dit()][idir];
          const EBCellFAB & cellVel = a_cellVeloc[dit()];

          ccpAverageStressToFaces(faceVel, cellVel, ebgraph, grid, idir, a_domain, a_dx);

          Real temp = 1.;
          temp += 1.;
        }
      ibox++;
    }
}
/*****/
void
ccpAverageFaceToCells(EBCellFAB &             a_cellData,
                      const EBFluxFAB &       a_fluxData,
                      const EBGraph &         a_ebGraph,
                      const Box &             a_grid,
                      const ProblemDomain &   a_domain,
                      const RealVect &        a_dx)
{
  IntVectSet ivsIrreg = a_ebGraph.getIrregCells(a_grid);
  int icomp = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const EBFaceFAB & faceData = a_fluxData[idir];
      const BaseFab<Real> & regFaceData =   faceData.getSingleValuedFAB();
      BaseFab<Real> &       regCellData = a_cellData.getSingleValuedFAB();
      FORT_CCPAVEFACETOCELL( CHF_FRA1(regCellData, idir),
                             CHF_CONST_FRA1(regFaceData, icomp),
                             CHF_CONST_INT(idir),
                             CHF_BOX(a_grid));

      for (VoFIterator vofit(ivsIrreg, a_ebGraph); vofit.ok(); ++vofit)
        {
          const VolIndex & vof = vofit();
          int numFaces = 0;
          Real cellVal = 0.;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex> faces = a_ebGraph.getFaces(vof, idir, sit());
              //if we have faces, then use them.  otherwise, need to extrapolate to covered
              if (faces.size() > 0)
                {
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      cellVal += faceData(faces[iface], icomp);
                      numFaces++;
                    }
                }
              else
                {
                  const EBISBox&  ebisBox = a_cellData.getEBISBox();
                  Real extrapValCov = ccpGetCoveredExtrapValue(vof, idir, sit(),
                                                               faceData,
                                                               ebisBox,
                                                               a_grid,
                                                               a_domain,
                                                               a_dx,
                                                               icomp);


                  //add in covered face value as if it were a normal value
                  cellVal += extrapValCov;
                  numFaces++;
                }
            }
          if (numFaces > 1)
            {
              cellVal /= Real(numFaces);
            }
          a_cellData(vof, idir) = cellVal;
        }
    } //end loop over directions
}
/*****/
void
ccpAverageFaceToCellsScalar(EBCellFAB &             a_cellData,
                            const EBFluxFAB &       a_fluxData,
                            const EBGraph &         a_ebGraph,
                            const Box &             a_grid,
                            const ProblemDomain &   a_domain,
                            const RealVect &        a_dx,
                            const int &             a_dir)
{
  IntVectSet ivsIrreg = a_ebGraph.getIrregCells(a_grid);
  int icomp = 0;
  const EBFaceFAB & faceData = a_fluxData[a_dir];
  const BaseFab<Real> & regFaceData =   faceData.getSingleValuedFAB();
  BaseFab<Real> &       regCellData = a_cellData.getSingleValuedFAB();
  FORT_CCPAVEFACETOCELL( CHF_FRA1(regCellData, icomp),
                         CHF_CONST_FRA1(regFaceData, icomp),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(a_grid));

  for (VoFIterator vofit(ivsIrreg, a_ebGraph); vofit.ok(); ++vofit)
    {
      const VolIndex & vof = vofit();
      int numFaces = 0;
      Real cellVal = 0.;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Vector<FaceIndex> faces = a_ebGraph.getFaces(vof, a_dir, sit());
          //if we have faces, then use them.  otherwise, need to extrapolate to covered
          if (faces.size() > 0)
            {
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  cellVal += faceData(faces[iface], icomp);
                  numFaces++;
                }
            }
          else
            {
              const EBISBox&  ebisBox = a_cellData.getEBISBox();
              Real extrapValCov = ccpGetCoveredExtrapValue(vof, a_dir, sit(),
                                                           faceData,
                                                           ebisBox,
                                                           a_grid,
                                                           a_domain,
                                                           a_dx,
                                                           icomp);

              //add in covered face value as if it were a normal value
              cellVal += extrapValCov;
              numFaces++;
            }
        }
      if (numFaces > 1)
        {
          cellVal /= Real(numFaces);
        }
      a_cellData(vof, icomp) = cellVal;
    }
}

/*****/
void
ccpCellGradientFromFaceData(EBCellFAB &             a_cellGrad,
                            const EBFluxFAB &       a_fluxData,
                            const EBGraph &         a_ebGraph,
                            const Box &             a_grid,
                            const ProblemDomain &   a_domain,
                            const RealVect &        a_dx,
                            const int &             a_dir)
{
  IntVectSet ivsIrreg = a_ebGraph.getIrregCells(a_grid);
  int icomp = 0;
  const EBFaceFAB & faceData = a_fluxData[a_dir];
  const BaseFab<Real> & regFaceData =   faceData.getSingleValuedFAB();
  BaseFab<Real> &       regCellData = a_cellGrad.getSingleValuedFAB();
  FORT_CCPCELLGRADFROMFACEDATA( CHF_FRA1(regCellData, icomp),
                                CHF_CONST_FRA1(regFaceData, icomp),
                                CHF_CONST_INT(a_dir),
                                CHF_CONST_REAL(a_dx[0]),
                                CHF_BOX(a_grid));

  for (VoFIterator vofit(ivsIrreg, a_ebGraph); vofit.ok(); ++vofit)
    {
      const VolIndex & vof = vofit();
      Real cellVal = 0.;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Vector<FaceIndex> faces = a_ebGraph.getFaces(vof, a_dir, sit());
          //if we have faces, then use them.  otherwise, need to extrapolate to covered
          if (faces.size() > 0)
            {
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  cellVal += isign*faceData(faces[iface], icomp);
                }
            }
          else
            {
              const EBISBox&  ebisBox = a_cellGrad.getEBISBox();
              Real extrapValCov = ccpGetCoveredExtrapValue(vof, a_dir, sit(),
                                                           faceData,
                                                           ebisBox,
                                                           a_grid,
                                                           a_domain,
                                                           a_dx,
                                                           icomp);

              //add in covered face value as if it were a normal value
              cellVal += isign*extrapValCov;
            }
        }
      a_cellGrad(vof, icomp) = cellVal/a_dx[0];
    }
}

//do multi-d extrapolation a-la colella, graves modiano keen
//if cannot do that do extrapolation in face direction
//if cannot do that do it in tangential directions
Real
ccpGetCoveredExtrapValue(const VolIndex&         a_vof,
                         const int&              a_idir,
                         const Side::LoHiSide    a_side,
                         const EBFaceFAB &       a_faceData,
                         const EBISBox &         a_ebisBox,
                         const Box &             a_grid,
                         const ProblemDomain &   a_domain,
                         const RealVect &        a_dx,
                         const int&              a_icomp)
{

  Real extrapValCov = 0.;

  bool dropOrder = false;
  RealVect normal = a_ebisBox.normal(a_vof);
  /**/
  ccpJohansenExtrapFaceToCovered(dropOrder,
                                 extrapValCov,
                                 a_faceData,
                                 a_ebisBox,
                                 a_vof,
                                 a_idir,
                                 a_side,
                                 normal,
                                 a_dx,
                                 a_icomp) ;


  /**/
  if (dropOrder)
    {
      extrapValCov = ccpOneDCoveredExtrapValue(a_vof, a_idir, a_side,
                                               a_faceData,
                                               a_ebisBox,
                                               a_grid,
                                               a_domain,
                                               a_dx,
                                               a_icomp);
    }
  return extrapValCov;
}
/**/
//do one dimensional extrapolation
//prefer face direction
Real
ccpOneDCoveredExtrapValue(const VolIndex&         a_vof,
                          const int&              a_idir,
                          const Side::LoHiSide    a_side,
                          const EBFaceFAB &       a_faceData,
                          const EBISBox &         a_ebisBox,
                          const Box &             a_grid,
                          const ProblemDomain &   a_domain,
                          const RealVect &        a_dx,
                          const int&              a_icomp)
{
  const bool doNormalDir = true;
  const bool doTangenDir = true;

  Real extrapValCovNormal = -1.e99;
  Real extrapValCovTangen = 0.0;
  int numExtrapCoveredTangen = 0;
  int orderNormal = 0;
  int orderTangen = 0;

  //extrap normal to idir
  if (doNormalDir)
    {
      Real nearVal    = 1.e99;
      Real farVal     = 1.e99;
      Real farFarVal  = 1.e99;
      Side::LoHiSide flipSide = flip(a_side);
      Vector<FaceIndex> nearFaces  = a_ebisBox.getFaces(a_vof, a_idir, flipSide);
      bool hasNearFace    = false;
      bool hasFarFace     = false;
      bool hasFarFarFace  = false;
      hasNearFace = ((nearFaces.size() == 1) && !nearFaces[0].isBoundary());
      if (hasNearFace)
        {
          nearVal = a_faceData(nearFaces[0], 0);
          if (!nearFaces[0].isBoundary())
            {
              VolIndex farVoF = nearFaces[0].getVoF(flipSide);
              Vector<FaceIndex> farFaces = a_ebisBox.getFaces(farVoF, a_idir, flipSide);
              hasFarFace = ((farFaces.size() == 1) && !farFaces[0].isBoundary());
              if (hasFarFace)
                {
                  farVal = a_faceData(farFaces[0], 0);
                  const VolIndex & nextNextVoF = farFaces[0].getVoF(flipSide);
                  Vector<FaceIndex> farFarFaces = a_ebisBox.getFaces(nextNextVoF, a_idir, flipSide);
                  hasFarFarFace = ((farFarFaces.size() == 1) && !farFarFaces[0].isBoundary());
                  if (hasFarFarFace)
                    {
                      farFarVal = a_faceData(farFarFaces[0], 0);
                    }
                }
            }
        }
      if (hasNearFace && hasFarFace && hasFarFarFace)
        {
          extrapValCovNormal = (4.*nearVal - 3.*farVal + farFarVal)/2.;
          orderNormal = 3;
        }
      else if (hasNearFace && hasFarFace)
        {
          extrapValCovNormal = 2.*nearVal - farVal;
          orderNormal = 2;
        }
      else if (hasNearFace)
        {
          extrapValCovNormal = nearVal;
          orderNormal = 1;
        }
      else
        {
          extrapValCovNormal = 0.0;
          orderNormal = 0;
        }
    }

  //extrap tangential to a_idir
  if (doTangenDir)
    {
      for (int jdir = 0; jdir < SpaceDim; jdir++)
        {
          if (a_idir != jdir)
            {
              for (SideIterator tanSit; tanSit.ok(); ++tanSit)
                {

                  Real nearVal    = 1.e99;
                  Real farVal     = 1.e99;
                  Vector<FaceIndex> tanFaces  = a_ebisBox.getFaces(a_vof, jdir, tanSit());
                  bool hasNearFace     = false;
                  bool hasFarFace      = false;
                  if (tanFaces.size() == 1 && !tanFaces[0].isBoundary())
                    {
                      VolIndex tanVoF = tanFaces[0].getVoF(tanSit());
                      Vector<FaceIndex> nearFaces  = a_ebisBox.getFaces(tanVoF, a_idir, a_side);
                      if (nearFaces.size() == 1 && !nearFaces[0].isBoundary())
                        {
                          const IntVect& faceIV = nearFaces[0].gridIndex(a_side);
                          const Box & faceRegion = a_faceData.getRegion();
                          if (faceRegion.contains(faceIV))
                            {
                              hasNearFace = true;
                              nearVal = a_faceData(nearFaces[0], 0);
                              Vector<FaceIndex> farTanFaces  = a_ebisBox.getFaces(tanVoF, jdir, tanSit());
                              if (farTanFaces.size() == 1 && !farTanFaces[0].isBoundary())
                                {
                                  VolIndex farTanVoF = farTanFaces[0].getVoF(tanSit());
                                  Vector<FaceIndex> farFaces = a_ebisBox.getFaces(farTanVoF, a_idir, a_side);
                                  if (farFaces.size() == 1 && !farFaces[0].isBoundary())
                                    {
                                      const IntVect& farFaceIV = farFaces[0].gridIndex(a_side);
                                      if (faceRegion.contains(farFaceIV))
                                        {
                                          hasFarFace = true;
                                          farVal = a_faceData(farFaces[0], 0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                  if (hasNearFace && hasFarFace)
                    {
                      if (orderTangen==2)
                        {//add to previus order 2 extraps
                          extrapValCovTangen += 2.*nearVal - farVal;
                          numExtrapCoveredTangen++;
                          orderTangen = 2;
                        }
                      else if (orderTangen==0 || orderTangen==1)
                        {//this is the first time with this order
                          extrapValCovTangen = 2.*nearVal - farVal;
                          numExtrapCoveredTangen = 1;
                          orderTangen = 2;
                        }
                    }
                  else if (hasNearFace)
                    {
                      if (orderTangen==1)
                        {//add to previus order 1 extraps
                          extrapValCovTangen += nearVal;
                          numExtrapCoveredTangen++;
                          orderTangen = 1;
                        }
                      else if (orderTangen==0)
                        {//this is the first time with this order
                          extrapValCovTangen = nearVal;
                          numExtrapCoveredTangen = 1;
                          orderTangen = 1;
                        }
                    }
                }
            }
        }
    } //end tangential extrapolation

  //now we compute the extrapValCov with the highest order possible
  //   based on the different attempts (normal and tangential from above)
  Real extrapValCov;
  extrapValCov = extrapValCovNormal;
  if (orderNormal > orderTangen)
    {
      extrapValCov = extrapValCovNormal;
    }
  else if (orderNormal < orderTangen)
    {
      CH_assert(numExtrapCoveredTangen > 0);
      extrapValCov = extrapValCovTangen/Real(numExtrapCoveredTangen);
    }
  else
    {
      if (orderNormal > 0)
        {//average all extrap values since they are of the same order
          CH_assert(numExtrapCoveredTangen > 0);
          extrapValCov = (extrapValCovTangen/Real(numExtrapCoveredTangen) + extrapValCovNormal)/2.0;
        }
      else
        {
          CH_assert(orderNormal==0);
          CH_assert(orderTangen==0);
          extrapValCov = 0.0;
//           MayDay::Error("EBLevelCCProjector::ccpOneDCoveredExtrapValue:: can't extrap to covered...");
        }
    }
  return extrapValCov;
}
/*****/
void
ccpAverageFaceToCells(LevelData<EBCellFAB> &        a_cellData,
                      const LevelData<EBFluxFAB> &  a_macData,
                      const DisjointBoxLayout &     a_grids,
                      const EBISLayout &            a_ebisl,
                      const ProblemDomain &         a_domain,
                      const RealVect &              a_dx)
{
  LevelData<EBFluxFAB>& nonConstFlux = (LevelData<EBFluxFAB>&) a_macData;
  nonConstFlux.exchange();
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      const Box grid = a_grids.get(dit());
      const EBFluxFAB &  macData =  a_macData[dit()];
      EBCellFAB &       cellData = a_cellData[dit()];

      ccpAverageFaceToCells(cellData, macData, ebgraph, grid, a_domain, a_dx);
    }
}
/*****/
void
ccpAverageFaceToCellsScalar(LevelData<EBCellFAB> &        a_cellData,
                            const LevelData<EBFluxFAB> &  a_macData,
                            const DisjointBoxLayout &     a_grids,
                            const EBISLayout &            a_ebisl,
                            const ProblemDomain &         a_domain,
                            const RealVect &              a_dx,
                            const int &                   a_dir)
{
  LevelData<EBFluxFAB>& nonConstFlux = (LevelData<EBFluxFAB>&) a_macData;
  nonConstFlux.exchange();
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      const Box grid = a_grids.get(dit());
      const EBFluxFAB &  macData =  a_macData[dit()];
      EBCellFAB &       cellData = a_cellData[dit()];

      ccpAverageFaceToCellsScalar(cellData, macData, ebgraph, grid, a_domain, a_dx, a_dir);
    }
}
/*****/
void
ccpCellGradientFromFaceData(LevelData<EBCellFAB> &        a_cellData,
                            const LevelData<EBFluxFAB> &  a_macData,
                            const DisjointBoxLayout &     a_grids,
                            const EBISLayout &            a_ebisl,
                            const ProblemDomain &         a_domain,
                            const RealVect &              a_dx,
                            const int &                   a_dir)
{
  LevelData<EBFluxFAB>& nonConstFlux = (LevelData<EBFluxFAB>&) a_macData;
  nonConstFlux.exchange();
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      const Box grid = a_grids.get(dit());
      const EBFluxFAB &  macData =  a_macData[dit()];
      EBCellFAB &       cellData = a_cellData[dit()];

      ccpCellGradientFromFaceData(cellData, macData, ebgraph, grid, a_domain, a_dx, a_dir);
    }
}
/***/
void
ccpGetExtrapInformation( Tuple<int,SpaceDim-1>&   a_planeDirs,  // plane of interpolation
                         int&                     a_interpPlane, //plane (line in 2d) normal  in which interpolation occurs
                         RealVect&                a_startingPt, //location of center of the covered cell.
                         Vector<VolIndex>&        a_vofsStencil, //all vofs in monotone path
                         IntVect&                 a_signNormal,
                         const EBISBox&           a_ebisBox,
                         const RealVect&          a_normal,     //
                         const RealVect&          a_dx,         //
                         const VolIndex&          a_vof,        //vof that has covered face
                         const int&               a_faceDir,    //direction of both covered face and data
                         const Side::LoHiSide&    a_sd)         //whether low or hi side is covered
{

  //find the interpolation plane---plane with biggest normal
  a_interpPlane =0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if ((Abs(a_normal[idir]/a_dx[idir])) > (Abs(a_normal[a_interpPlane]/a_dx[a_interpPlane])))
        {
          a_interpPlane = idir;
        }
    }
  //compute tangential directions
  a_planeDirs = PolyGeom::computeTanDirs(a_interpPlane);

  a_signNormal = IntVect::Unit;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_normal[idir] < 0.0)
        {
          a_signNormal[idir] = -1;
        }
    }

  //compute the location in space of the center of the covered face
  const IntVect& iv = a_vof.gridIndex();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != a_faceDir)
        {
          a_startingPt[idir] = a_dx[idir]*(iv[idir] + 0.5);
        }
      else
        {
          //sign(a_side)      goes between -1 and 1
          //sign(a_side) + 1  goes between  0 and 2
          int faceOffset = (sign(a_sd) + 1)/2;
          a_startingPt[idir] = a_dx[idir]*(iv[idir] + faceOffset);
        }
    }
  int radius = 2;
  EBArith::getAllVoFsInMonotonePath(a_vofsStencil, a_vof, a_ebisBox, radius);
}
/**/
void
ccpGetFacesInInterpolationStencil(bool&                          a_dropOrder,     //whether all faces were able to be found
                                  Vector<FaceIndex>&             a_faces,         //faces in stencil.
                                  Vector<RealVect>&              a_faceLoc,       //location of face centers in interpolation
                                                                                  //relative to center of intersection plane
                                  RealVect&                      a_intersectLoc,
                                  const RealVect&                a_startingPt,
                                  const Vector<VolIndex>&        a_vofsStencil,   //all vofs in monotone path
                                  const IntVect&                 a_signNormal,
                                  const Tuple<int,SpaceDim-1>&   a_planeDirs,     // plane of interpolation
                                  const int&                     a_interpPlane,
                                  const EBISBox&                 a_ebisBox,
                                  const RealVect&                a_normal,        //
                                  const RealVect&                a_dx,            //
                                  const VolIndex&                a_vof,           //vof that has covered face
                                  const int&                     a_faceDir,       //direction of both covered face and data
                                  const Side::LoHiSide&          a_sd,            //whether low or hi side is covered
                                  const int&                     a_istepOut )     //number of steps away from covered face
{
  //find which normal is biggest
  //calculate the location where the line intersects the interpolation plane
  IntVect cellsFromStart;
  cellsFromStart[a_interpPlane] = a_signNormal[a_interpPlane]*(a_istepOut + 1);
  Real distToPlane = Abs(a_dx[a_interpPlane]*cellsFromStart[a_interpPlane]);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //wacked normal condition
      if (Abs(a_normal[a_interpPlane]) < 1.0e-6)
        {
          a_dropOrder = true;
          return;
        }
      Real normalRat = Abs(a_normal[idir]/a_normal[a_interpPlane]);

      Real distMoved =  distToPlane*normalRat;
      a_intersectLoc[idir]= a_startingPt[idir] + a_signNormal[idir]*distMoved;

      if (idir != a_interpPlane)
        {
          cellsFromStart[idir] = a_signNormal[idir]*int(distMoved/a_dx[idir]);
        }
    }
  Vector<IntVect> ivStenLo;
  //if the covered face were actually a face, this would be its low cell
  IntVect ivLoCoveredFace = a_vof.gridIndex();
  if (a_sd == Side::Lo)
    {
      ivLoCoveredFace -= BASISV(a_faceDir);
    }

  IntVect ivLoInter = ivLoCoveredFace + cellsFromStart;
  IntVect ivHiInter = ivLoInter + BASISV(a_faceDir);
  FaceIndex faceCenter(VolIndex(ivLoInter, 0), VolIndex(ivHiInter, 0), a_faceDir);
  RealVect loCenterLoc = EBArith::getFaceLocation(faceCenter, a_dx, RealVect::Zero);

  IntVect  planeSign = IntVect::Unit;
  Real eps = 1.0e-10;
  for (int itan = 0; itan < SpaceDim-1; itan++)
    {
      //if it hits on a corner point i do not want it to drop order just out of spite
      if (Abs(a_intersectLoc[a_planeDirs[itan]] - loCenterLoc[a_planeDirs[itan]]) < eps)
        {
          if (cellsFromStart[a_planeDirs[itan]] > 0)
            {
              planeSign[a_planeDirs[itan]] = -1;
            }
          else
            {
              planeSign[a_planeDirs[itan]] =  1;
            }
        }
      else if (a_intersectLoc[a_planeDirs[itan]] > loCenterLoc[a_planeDirs[itan]])
        {
          planeSign[a_planeDirs[itan]] = 1;
        }
      else if (a_intersectLoc[a_planeDirs[itan]] < loCenterLoc[a_planeDirs[itan]])
        {
          planeSign[a_planeDirs[itan]] = -1;
        }
    }

  // this prevents a drop order case next to domain boundaires
  // when the normal in the tangential direction equals zero
  IntVect flipSign = IntVect::Unit;
  {
    IntVect tangenDir = 9999*IntVect::Unit;
    tangenDir[0] = a_planeDirs[0];
#if CH_SPACEDIM==3
    tangenDir[1] = a_planeDirs[1];
#endif
    const ProblemDomain& domain = a_ebisBox.getDomain();
    for (int idir = 0; idir < SpaceDim-1; idir++)
      {
        if (a_normal[tangenDir[idir]] == 0.0)
          {
            //check if we are next to the high domain edge
            IntVect ivHi = a_vof.gridIndex();
            ivHi[tangenDir[idir]] += 1;
            if (!domain.contains(ivHi))
              {
                flipSign[idir] = -1;
              }
          }
      }
  }

#if CH_SPACEDIM==2
  ivStenLo.resize(2);

  ivStenLo[0] = ivLoInter;
  ivStenLo[1] = ivLoInter + flipSign[0]*planeSign[a_planeDirs[0]]*BASISV(a_planeDirs[0]);

  a_faceLoc.resize(2, RealVect::Zero);
  a_faceLoc[1][a_planeDirs[0]] = a_dx[a_planeDirs[0]];

#elif CH_SPACEDIM==3
  ivStenLo.resize(4);
  ivStenLo[0] = ivLoInter;
  ivStenLo[1] = ivLoInter + (flipSign[0]*planeSign[a_planeDirs[0]]*BASISV(a_planeDirs[0]));
  ivStenLo[2] = ivLoInter + (flipSign[1]*planeSign[a_planeDirs[1]]*BASISV(a_planeDirs[1]));
  ivStenLo[3] = ivLoInter + (flipSign[1]*planeSign[a_planeDirs[1]]*BASISV(a_planeDirs[1]) +
                             flipSign[0]*planeSign[a_planeDirs[0]]*BASISV(a_planeDirs[0]));

  //if x y are planeDirs
  //first point is at the origin.   The other points are
  //dx or dy away (can ignore signs here)
  // point  0 is at  (0  , 0 )
  // point  1 is at  (dx , 0 )
  // point  2 is at  (0  , dy)
  // point  3 is at  (dx , dy)
  a_faceLoc.resize(4, RealVect::Zero);
  a_faceLoc[1][a_planeDirs[0]] = a_dx[a_planeDirs[0]];
  a_faceLoc[2][a_planeDirs[1]] = a_dx[a_planeDirs[1]];
  a_faceLoc[3][a_planeDirs[0]] = a_dx[a_planeDirs[0]];
  a_faceLoc[3][a_planeDirs[1]] = a_dx[a_planeDirs[1]];

#else
  bogus_ch_spacedim();
#endif

  Vector<IntVect> ivStenHi(ivStenLo.size());
  for (int isten = 0; isten < ivStenLo.size(); isten++)
    {
      //no sign of normal here.  just doing lo-hi sides of the face.
      ivStenHi[isten]  = ivStenLo[isten]  + BASISV(a_faceDir);
    }

  //run away if we are dealing  with multivalued cells
  //the only way to do those properly is to do monotone path stuff
  //from the original cell and that is too expensive
  for (int isten = 0; isten < ivStenLo.size(); isten++)
    {
      //first check to see everything is in domain.
      const ProblemDomain& domain = a_ebisBox.getDomain();
      if ((!domain.contains(ivStenLo[isten])) ||
         (!domain.contains(ivStenHi[isten])))
        {
          a_dropOrder = true;
          return;
        }

      if ((a_ebisBox.numVoFs(ivStenLo[isten]) > 1) ||
         (a_ebisBox.numVoFs(ivStenHi[isten]) > 1))
        {
          a_dropOrder = true;
          return;
        }
    }

  //check to see if all the cells we will need are in monotone path from
  // the center  cell.  If not, return with dropOrder = true
  a_faces.resize(ivStenLo.size());
  for (int isten = 0; isten < ivStenLo.size(); isten++)
    {
      VolIndex whichVoFLo, whichVoFHi;  //gets the vof in cell that is in stencil
      bool isVoFLo = EBArith::isVoFHere(whichVoFLo,a_vofsStencil, ivStenLo[isten]);
      bool isVoFHi = EBArith::isVoFHere(whichVoFHi,a_vofsStencil, ivStenHi[isten]);
      if ((!isVoFLo) || (!isVoFHi))
        {
          a_dropOrder = true;
          return;
        }
      a_faces[isten] = FaceIndex(whichVoFLo, whichVoFHi, a_faceDir);
    }

  //now normalize intersectLoc so that it starts at loCenterLoc and goes
  // between 0 and the appropriate dx.  this is so it can be reused in the interpolation
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_intersectLoc[idir] = Abs(a_intersectLoc[idir] - loCenterLoc[idir]);
    }
}
/***/
void
ccpLinearInterp(Real &                    a_dataOnLine,
                const Vector<RealVect>&   a_faceLoc,
                const RealVect&           a_intersectLoc,
                const Vector<Real>&       a_interpolationData,
                const int&                a_planeDir)
{
  CH_assert(a_faceLoc.size()           >= 2);
  CH_assert(a_interpolationData.size() >= 2);
  Real x0 = a_faceLoc[0][a_planeDir];
  Real x1 = a_faceLoc[1][a_planeDir];
  Real f0 = a_interpolationData[0];
  Real f1 = a_interpolationData[1];

  //linear function f = Ax + B
  Real A = (f0    - f1   )/(x0-x1);
  Real B = (x0*f1 - x1*f0)/(x0-x1);
  Real x = a_intersectLoc[a_planeDir];

  a_dataOnLine = A*x + B;
}

/***/
void
ccpBilinearInterp(Real &                             a_dataOnLine,
                  const Vector<RealVect>&            a_faceLoc,
                  const RealVect&                    a_intersectLoc,
                  const Vector<Real>&                a_interpolationData,
                  const Tuple<int, SpaceDim-1>&      a_planeDirs)
{
  CH_assert(a_faceLoc.size()           >= 4);
  CH_assert(a_interpolationData.size() >= 4);

  //bilinear function f = Ax + By + Cxy + D
  //I do not want to solve the general 4x4 system so I
  //will use the fact that I know which point is which.

  // point  0 is at  (0  , 0 )
  // point  1 is at  (dx , 0 )
  // point  2 is at  (0  , dy)
  // point  3 is at  (dx , dy)
  Real dx = a_faceLoc[3][a_planeDirs[0]];
  Real dy = a_faceLoc[3][a_planeDirs[1]];
  Real f0 = a_interpolationData[0];
  Real f1 = a_interpolationData[1];
  Real f2 = a_interpolationData[2];
  Real f3 = a_interpolationData[3];

  Real D = f0;
  Real A = (f1-D)/dx;
  Real B = (f2-D)/dy;
  Real C = (f3 - (A*dx + B*dy + D))/(dx*dy);

  Real x = a_intersectLoc[a_planeDirs[0]];
  Real y = a_intersectLoc[a_planeDirs[1]];

  a_dataOnLine = A*x + B*y + C*x*y + D;
}

/***/
void
ccpJohansenExtrapFaceToCovered(bool&                   a_dropOrder,
                               Real&                   a_extrapVal,
                               const EBFaceFAB&        a_primFace,
                               const EBISBox&          a_ebisBox,
                               const VolIndex&         a_vof,
                               const int&              a_faceDir,
                               const Side::LoHiSide&   a_sd,
                               const RealVect&         a_normal,
                               const RealVect&         a_dx,
                               const int&              a_icomp)
{

  Tuple<int,SpaceDim-1>   planeDirs;    //arranged in order of increasing size of normal
  RealVect                startingPt;   //location of center of the covered cell.
  IntVect                 signNormal;  //sign of normals components
  int                     interpPlane;  //plane (line in 2d) in which interpolation occurs
  Vector<VolIndex>        vofsStencil; //all vofs in monotone path

  ccpGetExtrapInformation( planeDirs,
                           interpPlane,
                           startingPt,
                           vofsStencil,
                           signNormal,
                           a_ebisBox,
                           a_normal,
                           a_dx,
                           a_vof,
                           a_faceDir,
                           a_sd) ;

  Real dataOnLine[2];
  a_dropOrder = false;
  //get data at each intersection point tin the stencil
  for (int iStepOut = 0; iStepOut < 2; iStepOut++)
    {
      //faces in stencil.
      Vector<FaceIndex>             faces;
      //location of face centers in interpolation
      Vector<RealVect>              faceLoc;
      RealVect intersectLoc;

      ccpGetFacesInInterpolationStencil( a_dropOrder,
                                         faces,
                                         faceLoc,
                                         intersectLoc,
                                         startingPt,
                                         vofsStencil,
                                         signNormal,
                                         planeDirs,
                                         interpPlane,
                                         a_ebisBox,
                                         a_normal,
                                         a_dx,
                                         a_vof,
                                         a_faceDir,
                                         a_sd,
                                         iStepOut );
      if (a_dropOrder)
        {
          //just get out if we do not have the stencil
          return;
        }

      Vector<Real> interpolationData(faces.size());
      for (int iface = 0; iface < faces.size(); iface++)
        {
          interpolationData[iface] = a_primFace(faces[iface], a_icomp);
        }

#if CH_SPACEDIM==2

      ccpLinearInterp(dataOnLine[iStepOut],
                      faceLoc,
                      intersectLoc,
                      interpolationData,
                      planeDirs[0]);

#elif CH_SPACEDIM==3

      ccpBilinearInterp(dataOnLine[iStepOut],
                        faceLoc,
                        intersectLoc,
                        interpolationData,
                        planeDirs);

#else
      //bogus ch_spacedim_macro
      this_is_gonna_be_ugly_in_4d();

#endif

    }

  if (!a_dropOrder)
    {
      a_extrapVal = 2.*dataOnLine[0] - dataOnLine[1];
    }
}
/***/

#include "NamespaceFooter.H"
