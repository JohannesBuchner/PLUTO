#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DebugOut.H"
#include "EBLevelAdvect.H"
#include "BaseIVFactory.H"
#include "BaseIFFactory.H"
#include "BaseIFFAB.H"
#include "EBFluxFAB.H"
#include "FaceIterator.H"
#include "REAL.H"
#include "EBCellFactory.H"
#include "FaceIterator.H"
#include "EBPatchAdvectF_F.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include "EBAMRIO.H"
#include "NamespaceHeader.H"
const IntVect    f_ivLo(D_DECL(143,31,0));
const IntVect    f_ivHi(D_DECL(143,32,0));
const IntVect f_ivLoGho(D_DECL(144,31,0));
const IntVect f_ivHiGho(D_DECL(144,32,0));

void
EBLevelAdvect::
resetBCs(const RefCountedPtr<EBPhysIBCFactory>&  a_advectBC)
{
  if (m_isDefined)
    {
      for (DataIterator dit = m_ebPatchAdvect.dataIterator(); dit.ok(); ++dit)
        {
          m_ebPatchAdvect[dit()]->setEBPhysIBC(*a_advectBC);
        }
    }
}
/*****************************/
EBLevelAdvect::
EBLevelAdvect()
{
  m_isDefined = false;
}
/*****************************/
bool
EBLevelAdvect::
isDefined() const
{
  return m_isDefined;
}
/*****************************/
EBLevelAdvect::
~EBLevelAdvect()
{
  if (m_isDefined)
    {
      for (DataIterator dit = m_ebPatchAdvect.dataIterator(); dit.ok(); ++dit)
        {
          if (m_ebPatchAdvect[dit()] != NULL)
            {
              delete m_ebPatchAdvect[dit()];
              m_ebPatchAdvect[dit()] = NULL;
            }
        }
    }
}
/*****************************/
EBLevelAdvect::
EBLevelAdvect(const DisjointBoxLayout&           a_thisDBL,
              const DisjointBoxLayout&           a_coarDBL,
              const EBISLayout&                  a_thisEBISL,
              const EBISLayout&                  a_coarEBISL,
              const ProblemDomain&               a_domain,
              const int&                         a_nRefine,
              const RealVect&                    a_dx,
              const bool&                        a_hasCoarser,
              const bool&                        a_hasFiner,
              const EBPatchGodunovFactory* const a_patchGodunov,
              const EBIndexSpace*          const a_eb)
{
  CH_TIME("EBLevelAdvect::EBLevelAdvect");
  m_isDefined = false;
  define(a_thisDBL,
         a_coarDBL,
         a_thisEBISL,
         a_coarEBISL,
         a_domain,
         a_nRefine,
         a_dx,
         a_hasCoarser,
         a_hasFiner,
         a_patchGodunov,
         a_eb);
}

/*****************************/
void
EBLevelAdvect::
define(const DisjointBoxLayout&           a_thisDBL,
       const DisjointBoxLayout&           a_coarDBL,
       const EBISLayout&                  a_thisEBISL,
       const EBISLayout&                  a_coarEBISL,
       const ProblemDomain&               a_domain,
       const int&                         a_nRefine,
       const RealVect&                    a_dx,
       const bool&                        a_hasCoarser,
       const bool&                        a_hasFiner,
       const EBPatchGodunovFactory* const a_patchGodunov,
       const EBIndexSpace*          const a_eb)
{
  CH_TIME("EBLevelAdvect::define");
  CH_assert(a_dx[0] > 0.0);
  CH_assert(a_nRefine > 0);

  if (m_isDefined)
    {
      for (DataIterator dit = m_ebPatchAdvect.dataIterator(); dit.ok(); ++dit)
        {
          if (m_ebPatchAdvect[dit()] != NULL)
            {
              delete m_ebPatchAdvect[dit()];
              m_ebPatchAdvect[dit()] = NULL;
            }
        }
    }

  m_isDefined = true;
  m_thisGrids    = a_thisDBL;
  m_thisEBISL    = a_thisEBISL;
  m_refRatCrse   = a_nRefine;
  m_dx           = a_dx;
  m_domain       = a_domain;
  m_hasCoarser   = a_hasCoarser;
  m_hasFiner     = a_hasFiner;
  if (m_hasCoarser)
    {
      m_coarGrids = a_coarDBL;
      m_coarEBISL = a_coarEBISL;
    }

  m_ebPatchAdvect.define(m_thisGrids);
  m_nVar         =  1;

  for (DataIterator dit = m_ebPatchAdvect.dataIterator(); dit.ok(); ++dit)
    {
      EBPatchGodunov* ebPatchGodunov = a_patchGodunov->create();
      m_ebPatchAdvect[dit()] = dynamic_cast<EBPatchAdvect*>(ebPatchGodunov);
      if (m_ebPatchAdvect[dit()] == NULL)
        {
          MayDay::Error("problem in casting to patch advection class");
        }

      m_ebPatchAdvect[dit()]->define(m_domain, m_dx);
      IntVectSet cfivs; //not used here.  only used in flux interpolation
      Real time = 0; //needs to get reset later
      Real dt   = 0; //needs to get reset later
      m_ebPatchAdvect[dit()]->setValidBox(m_thisGrids[dit()], m_thisEBISL[dit()], cfivs, time, dt);
    }
  m_nGhost = 4;

  if (m_hasCoarser)
    {
      CH_TIME("EBLevelAdvect::define::fillPatchDefine");
      ProblemDomain domainCrse = coarsen(m_domain, m_refRatCrse);

      //patcher is defined with the number of conserved vars.
      m_fillPatch.define(m_thisGrids, m_coarGrids,
                         m_thisEBISL, m_coarEBISL,
                         domainCrse, m_refRatCrse, m_nVar,
                         m_nGhost, a_eb);

      m_fillPatchVel.define(m_thisGrids, m_coarGrids,
                            m_thisEBISL, m_coarEBISL,
                            domainCrse, m_refRatCrse, SpaceDim,
                            m_nGhost, a_eb);

    }

}
/*****************************/
void
EBLevelAdvect::
advectToFacesCol(LevelData< EBFluxFAB >&                         a_extrapState,
                 LayoutData< Vector< BaseIVFAB<Real>* > >&       a_coveredPrimLo,
                 LayoutData< Vector< BaseIVFAB<Real>* > >&       a_coveredPrimHi,
                 const LayoutData< Vector< Vector<VolIndex> > >& a_coveredFaceLo,
                 const LayoutData< Vector< Vector<VolIndex> > >& a_coveredFaceHi,
                 const LayoutData< Vector< IntVectSet> >&        a_coveredSetsLo,
                 const LayoutData< Vector< IntVectSet> >&        a_coveredSetsHi,
                 const LevelData< EBCellFAB >&                   a_consState,
                 const LevelData< EBCellFAB >&                   a_normalVel,
                 const LevelData< EBFluxFAB >&                   a_advectionVel,
                 const LevelData< EBCellFAB >*                   a_consStateCoarseOld,
                 const LevelData< EBCellFAB >*                   a_consStateCoarseNew,
                 const LevelData< EBCellFAB >*                   a_normalVelCoarseOld,
                 const LevelData< EBCellFAB >*                   a_normalVelCoarseNew,
                 const Real&                                     a_timeCoarseOld,
                 const Real&                                     a_timeCoarseNew,
                 const Real&                                     a_timeFine,
                 const Real&                                     a_dt,
                 const LevelData<EBCellFAB>* const               a_source,
                 const LevelData<EBCellFAB>* const               a_sourceCoarOld,
                 const LevelData<EBCellFAB>* const               a_sourceCoarNew)
{
  CH_TIME("EBLevelAdvect::advectToFacesCol");
  CH_assert(isDefined());

  //create temp data with the correct number of ghost cells
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  Interval consInterv(0, m_nVar-1);
  Interval intervSD(0, SpaceDim-1);

  // LevelData<EBCellFAB>& consTemp = (LevelData<EBCellFAB>&) a_consState;
  // LevelData<EBCellFAB>& veloTemp = (LevelData<EBCellFAB>&) a_normalVel;

  EBCellFactory factory(m_thisEBISL);
  LevelData<EBCellFAB> consTemp(m_thisGrids, m_nVar, ivGhost, factory);
  LevelData<EBCellFAB> veloTemp(m_thisGrids, SpaceDim, ivGhost, factory);
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      consTemp[dit()].setVal(0.);
    }

  a_consState.copyTo(consInterv, consTemp, consInterv);
  a_normalVel.copyTo(intervSD, veloTemp, intervSD);
  // Fill ghost cells using fillInterp, and copyTo.
  if (m_hasCoarser)
    {
      CH_TIME("fillPatch");
      m_fillPatch.interpolate(consTemp,
                              *a_consStateCoarseOld,
                              *a_consStateCoarseNew,
                              a_timeCoarseOld,
                              a_timeCoarseNew,
                              a_timeFine,
                              consInterv);

      m_fillPatchVel.interpolate(veloTemp,
                                 *a_normalVelCoarseOld,
                                 *a_normalVelCoarseNew,
                                 a_timeCoarseOld,
                                 a_timeCoarseNew,
                                 a_timeFine,
                                 intervSD);
    }
  // Exchange all the data between grids
  {
    CH_TIME("initial_exchange");
    consTemp.exchange(consInterv);
    veloTemp.exchange(intervSD);
  }

  LevelData<EBCellFAB>* srcTmpPtr = NULL;
  if (a_source != NULL)
    {
      // srcTmpPtr = (LevelData<EBCellFAB>*) a_source;

      srcTmpPtr = new LevelData<EBCellFAB>(m_thisGrids, m_nVar, ivGhost, factory);
      for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
        {
          (*srcTmpPtr)[dit()].setVal(0.);
        }
      a_source->copyTo(consInterv, *srcTmpPtr, consInterv);

      if ( (a_sourceCoarOld != NULL) &&
           (a_sourceCoarNew != NULL) &&
           (m_hasCoarser) )
        {
          CH_TIME("fillPatch");

          m_fillPatch.interpolate(*srcTmpPtr,
                                  *a_sourceCoarOld,
                                  *a_sourceCoarNew,
                                  a_timeCoarseOld,
                                  a_timeCoarseNew,
                                  a_timeFine,
                                  consInterv);
          {
            CH_TIME("initial_exchange");
            srcTmpPtr->exchange(consInterv);
          }
        }
    }

  {
    CH_TIME("advectToFaces");
    int ibox = 0;
    for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& cellBox = m_thisGrids.get(dit());
        const EBISBox& ebisBox = m_thisEBISL[dit()];
        if (!ebisBox.isAllCovered())
          {
            //the viscous term goes into here
            EBCellFAB dummy;
            EBCellFAB* srcPtr = &dummy;
            if (srcTmpPtr != NULL)
              {
                srcPtr = (EBCellFAB*)(&((*srcTmpPtr)[dit()]));
              }

            const EBCellFAB& source = *srcPtr;
            //unused in this case
            BaseIVFAB<Real> boundaryPrim;
            bool doBoundaryPrim = false;

            EBFluxFAB& extrapFAB  = a_extrapState[dit()];
            advectToFaces(extrapFAB,
                          boundaryPrim,
                          a_coveredPrimLo[dit()],
                          a_coveredPrimHi[dit()],
                          a_coveredFaceLo[dit()],
                          a_coveredFaceHi[dit()],
                          a_coveredSetsLo[dit()],
                          a_coveredSetsHi[dit()],
                          consTemp[dit()],
                          veloTemp[dit()],
                          a_advectionVel[dit()],
                          cellBox, ebisBox,
                          a_dt,   a_timeFine,
                          source,  dit(),doBoundaryPrim);

            ibox++;
          }
      }
  }
  if (srcTmpPtr != NULL)
    {
      delete srcTmpPtr;
    }
}
/*****************************/
void
EBLevelAdvect::
advectToFaces(EBFluxFAB&                         a_extrapState,
              BaseIVFAB<Real>&                   a_boundaryPrim,
              Vector< BaseIVFAB<Real>* >&        a_coveredPrimLo,
              Vector< BaseIVFAB<Real>* >&        a_coveredPrimHi,
              const Vector< Vector<VolIndex> >&  a_coveredFaceLo,
              const Vector< Vector<VolIndex> >&  a_coveredFaceHi,
              const Vector< IntVectSet  >&       a_coveredSetsLo,
              const Vector< IntVectSet  >&       a_coveredSetsHi,
              const EBCellFAB &                  a_consState,
              const EBCellFAB &                  a_normalVel,
              const EBFluxFAB &                  a_advectionVel,
              const Box&                         a_cellBox,
              const EBISBox&                     a_ebisBox,
              const Real&                        a_dt,
              const Real&                        a_time,
              const EBCellFAB &                  a_source,
              const DataIndex&                   a_dit,
              bool   a_doBoundaryPrim)
{
  CH_TIME("EBLevelAdvect::advectToFaces");
  IntVectSet cfivs; //not used here.  only used in flux interpolation
  m_ebPatchAdvect[a_dit]->setTimeAndDt(a_time, a_dt);
  EBCellFAB& primState = m_ebPatchAdvect[a_dit]->getPrimState();
  m_ebPatchAdvect[a_dit]->setVelocities(a_normalVel, a_advectionVel);

  //placeholder (no flattening used here)
  EBCellFAB flattening;
  //not reused
  EBCellFAB  slopesPrim[SpaceDim];
  EBCellFAB  slopesSeco[SpaceDim];
  bool verbose = false;


  m_ebPatchAdvect[a_dit]->extrapolatePrim(a_extrapState,
                                          primState,
                                          slopesPrim,
                                          slopesSeco,
                                          a_coveredPrimLo,
                                          a_coveredPrimHi,
                                          a_coveredSetsLo,
                                          a_coveredSetsHi,
                                          a_coveredFaceLo,
                                          a_coveredFaceHi,
                                          flattening,
                                          a_consState,
                                          a_source,
                                          a_cellBox,
                                          a_dit,
                                          verbose);


  if (a_doBoundaryPrim)
    {
      IntVectSet irregIVS = a_ebisBox.getIrregIVS(a_cellBox);
      m_ebPatchAdvect[a_dit]->computeEBIrregFlux(a_boundaryPrim,
                                                 primState,
                                                 slopesPrim,
                                                 irregIVS,
                                                 a_source);
    }
}
/*****************************/
void
EBLevelAdvect::
advectToFacesBCG(LevelData< EBFluxFAB >&                         a_extrapState,
                 const LevelData< EBCellFAB >&                   a_consState,
                 const LevelData< EBCellFAB >&                   a_normalVel,
                 const LevelData< EBFluxFAB >&                   a_advectionVel,
                 const LevelData< EBCellFAB >*                   a_consStateCoarseOld,
                 const LevelData< EBCellFAB >*                   a_consStateCoarseNew,
                 const LevelData< EBCellFAB >*                   a_normalVelCoarseOld,
                 const LevelData< EBCellFAB >*                   a_normalVelCoarseNew,
                 const Real&                                     a_timeCoarseOld,
                 const Real&                                     a_timeCoarseNew,
                 const Real&                                     a_timeFine,
                 const Real&                                     a_dt,
                 const LevelData<EBCellFAB>* const               a_source,
                 const LevelData<EBCellFAB>* const               a_sourceCoarOld,
                 const LevelData<EBCellFAB>* const               a_sourceCoarNew)
{
  CH_TIME("EBLevelAdvect::advectToFacesBCG (level)");
  CH_assert(isDefined());

  //create temp data with the correct number of ghost cells
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  Interval consInterv(0, m_nVar-1);
  Interval intervSD(0, SpaceDim-1);

  // LevelData<EBCellFAB>& consTemp = (LevelData<EBCellFAB>&) a_consState;
  // LevelData<EBCellFAB>& veloTemp = (LevelData<EBCellFAB>&) a_normalVel;

  EBCellFactory factory(m_thisEBISL);
  LevelData<EBCellFAB> consTemp(m_thisGrids, m_nVar, ivGhost, factory);
  LevelData<EBCellFAB> veloTemp(m_thisGrids, SpaceDim, ivGhost, factory);
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      consTemp[dit()].setVal(0.);
    }

  a_consState.copyTo(consInterv, consTemp, consInterv);
  a_normalVel.copyTo(intervSD, veloTemp, intervSD);
  // Fill ghost cells using fillInterp, and copyTo.
  if (m_hasCoarser)
    {
      CH_TIME("fillPatch");
      m_fillPatch.interpolate(consTemp,
                              *a_consStateCoarseOld,
                              *a_consStateCoarseNew,
                              a_timeCoarseOld,
                              a_timeCoarseNew,
                              a_timeFine,
                              consInterv);

      m_fillPatchVel.interpolate(veloTemp,
                                 *a_normalVelCoarseOld,
                                 *a_normalVelCoarseNew,
                                 a_timeCoarseOld,
                                 a_timeCoarseNew,
                                 a_timeFine,
                                 intervSD);
    }
  // Exchange all the data between grids
  {
    CH_TIME("initial_exchange");
    consTemp.exchange(consInterv);
    veloTemp.exchange(intervSD);
  }

  LevelData<EBCellFAB>* srcTmpPtr = NULL;
  if (a_source != NULL)
    {
      // srcTmpPtr = (LevelData<EBCellFAB>*) a_source;

      srcTmpPtr = new LevelData<EBCellFAB>(m_thisGrids, m_nVar, ivGhost, factory);
      for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
        {
          (*srcTmpPtr)[dit()].setVal(0.);
        }
      a_source->copyTo(consInterv, *srcTmpPtr, consInterv);

      if ( (a_sourceCoarOld != NULL) &&
           (a_sourceCoarNew != NULL) &&
           (m_hasCoarser) )
        {
          CH_TIME("fillPatch");

          m_fillPatch.interpolate(*srcTmpPtr,
                                  *a_sourceCoarOld,
                                  *a_sourceCoarNew,
                                  a_timeCoarseOld,
                                  a_timeCoarseNew,
                                  a_timeFine,
                                  consInterv);
          {
            CH_TIME("initial_exchange");
            srcTmpPtr->exchange(consInterv);
          }
        }
    }

  {
    CH_TIME("advectToFaces");
    int ibox = 0;
    for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& cellBox = m_thisGrids.get(dit());
        const EBISBox& ebisBox = m_thisEBISL[dit()];
        if (!ebisBox.isAllCovered())
          {
            //the viscous term goes into here
            EBCellFAB dummy;
            EBCellFAB* srcPtr = &dummy;
            if (srcTmpPtr != NULL)
              {
                srcPtr = (EBCellFAB*)(&((*srcTmpPtr)[dit()]));
              }

            const EBCellFAB& source = *srcPtr;
            //unused in this case
            BaseIVFAB<Real> boundaryPrim;
            bool doBoundaryPrim = false;

            EBFluxFAB& extrapFAB  = a_extrapState[dit()];
            advectToFacesBCG(extrapFAB,
                             boundaryPrim,
                             consTemp[dit()],
                             veloTemp[dit()],
                             a_advectionVel[dit()],
                             cellBox,
                             ebisBox,
                             a_dt,
                             a_timeFine,
                             source,
                             dit(),
                             doBoundaryPrim);



            ibox++;
          }
      }
  }
  if (srcTmpPtr != NULL)
    {
      delete srcTmpPtr;
    }
}
/*****************************/
void
EBLevelAdvect::
advectToFacesBCG(EBFluxFAB&                         a_extrapState,
                 BaseIVFAB<Real>&                   a_boundaryPrim,
                 const EBCellFAB &                  a_consState,
                 const EBCellFAB &                  a_normalVel,
                 const EBFluxFAB &                  a_advectionVel,
                 const Box&                         a_cellBox,
                 const EBISBox&                     a_ebisBox,
                 const Real&                        a_dt,
                 const Real&                        a_time,
                 const EBCellFAB &                  a_source,
                 const DataIndex&                   a_dit,
                 bool   a_doBoundaryPrim)
{
  CH_TIME("EBLevelAdvect::advectToFacesBCG (fluxfab)");
  IntVectSet cfivs; //not used here.  only used in flux interpolation
  m_ebPatchAdvect[a_dit]->setTimeAndDt(a_time, a_dt);
  // EBCellFAB& primState = m_ebPatchAdvect[a_dit]->getPrimState();
  m_ebPatchAdvect[a_dit]->setVelocities(a_normalVel, a_advectionVel);

  //placeholder (no flattening used here)
  EBCellFAB flattening;
  //not reused
  EBCellFAB  slopesPrim[SpaceDim];
  EBCellFAB  slopesSeco[SpaceDim];
  bool verbose = false;

  m_ebPatchAdvect[a_dit]->extrapolateBCG(a_extrapState,
                                         // primState,
                                         slopesPrim,
                                         slopesSeco,
                                         flattening,
                                         a_consState,
                                         a_source,
                                         a_cellBox,
                                         a_dit,
                                         verbose);

  if (a_doBoundaryPrim)
    {
      IntVectSet irregIVS = a_ebisBox.getIrregIVS(a_cellBox);
      m_ebPatchAdvect[a_dit]->computeEBIrregFlux(a_boundaryPrim,
                                                 // primState,
                                                 a_consState,
                                                 slopesPrim,
                                                 irregIVS,
                                                 a_source);
    }
}
/*****************************/
void
EBLevelAdvect::
computeNormalVel(LevelData<EBCellFAB>&                          a_normalVel,
                 const LevelData<EBFluxFAB>&                    a_advectionVel,
                 const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloLo,
                 const LayoutData<Vector<BaseIVFAB<Real> * > >& a_coveredVeloHi,
                 const LayoutData<Vector<Vector<VolIndex> > >&  a_coveredFaceLo,
                 const LayoutData<Vector<Vector<VolIndex> > >&  a_coveredFaceHi) const
{
  CH_TIME("EBLevelAdvect::computeNormalVel");
  for (DataIterator dit = m_thisGrids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& cellBox = m_thisGrids.get(dit());
      Box  grownBox = grow(cellBox, 1);
      grownBox &= m_domain;
      m_ebPatchAdvect[dit]->averageVelToCC(a_normalVel[dit()],
                                           a_advectionVel[dit()],
                                           a_coveredVeloLo[dit()],
                                           a_coveredVeloHi[dit()],
                                           a_coveredFaceLo[dit()],
                                           a_coveredFaceHi[dit()],
                                           grownBox);
    }

}
/*****************************/
#include "NamespaceFooter.H"
