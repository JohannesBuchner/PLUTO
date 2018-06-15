#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBCompositeCCProjector.H"
#include "NeumannPoissonEBBC.H"
#include "NeumannPoissonDomainBC.H"
#include "EBFluxFactory.H"
#include "EBAMRPoissonOp.H"
#include "BaseIFFactory.H"
#include "EBArith.H"
#include "EBAlias.H"
#include "EBAMRDataOps.H"

#include "NamespaceHeader.H"

/*****/
EBCompositeCCProjector::
EBCompositeCCProjector(const Vector<EBLevelGrid>&                      a_eblg,
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
                       EBCompositeMACProjector*                        a_inputMAC)
{
  CH_TIME("EBCompositeCCProjector::EBCompositeCCProjector");
  if (a_numLevels < 0)
    {
      m_numLevels = a_eblg.size();
    }
  else
    {
      CH_assert(a_numLevels <= a_eblg.size());
      m_numLevels = a_numLevels;
    }
  CH_assert(a_refRat.size() >= m_numLevels);

  m_eblg      = a_eblg;
  m_quadCFI   = a_quadCFI;
  m_refRat     = a_refRat;
  m_pwlCFI.resize(m_numLevels);
  int nvarquad = 1;
  int radius =1;
  for (int ilev = 1; ilev < m_numLevels; ilev++)
    {
      m_pwlCFI[ilev]  = RefCountedPtr<EBPWLFillPatch>(new  EBPWLFillPatch(m_eblg[ilev].getDBL(),
                                                                          m_eblg[ilev-1].getDBL(),
                                                                          m_eblg[ilev].getEBISL(),
                                                                          m_eblg[ilev-1].getEBISL(),
                                                                          m_eblg[ilev-1].getDomain(),
                                                                          m_refRat[ilev-1],
                                                                          nvarquad,
                                                                          radius,
                                                                          m_eblg[ilev].getEBISL().getEBIS()));
    }

  m_dx.resize(m_numLevels);
  m_dx[0]     = a_coarsestDx;
  for (int ilev = 1; ilev < m_numLevels; ilev++)
    {
      m_dx[ilev]     = m_dx[ilev-1]/m_refRat[ilev-1];
    }

  if (a_inputMAC != NULL)
    {
      m_macProjector = a_inputMAC;
      m_externMAC    = true;
    }
  else
    {
      m_externMAC    = false;
      m_macProjector = new EBCompositeMACProjector(a_eblg, a_refRat, a_quadCFI,
                                                   a_coarsestDx, a_origin,
                                                   a_baseDomainBCVel,
                                                   a_baseDomainBCPhi,
                                                   a_ebbcPhi,
                                                   a_subtractOffMean,
                                                   a_numLevels,
                                                   a_verbosity,
                                                   a_numPreCondIters,
                                                   a_time,
                                                   a_relaxType,
                                                   a_bottomSolverType);

      //default values
      int numSmooth = 4;
      int itermax = 77;
      int mgCycle = 1;
      Real hang = 1.0e-20;
      Real tolerance = 1.0e-12;
      Real normThresh= 1.0e-12;
      setSolverParams(numSmooth,  itermax, mgCycle, hang, tolerance, a_verbosity,normThresh);
    }

  m_macVelo.resize(m_numLevels, NULL);
  m_macGrad.resize(m_numLevels, NULL);

  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      EBFluxFactory ebfluxfact(m_eblg[ilev].getEBISL());
      m_macVelo[ilev] = new LevelData<EBFluxFAB>(m_eblg[ilev].getDBL(), 1,   IntVect::Unit, ebfluxfact);
      m_macGrad[ilev] = new LevelData<EBFluxFAB>(m_eblg[ilev].getDBL(), 1, 2*IntVect::Unit, ebfluxfact);
    }
}
/*****/
EBCompositeCCProjector::
~EBCompositeCCProjector()
{
  if (!m_externMAC)
    {
      delete m_macProjector;
    }
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      delete m_macVelo[ilev];
      delete m_macGrad[ilev];
    }
}
/*****/
void
EBCompositeCCProjector::
setSolverParams(int  a_numSmooths,
                int  a_iterMax,
                int  a_mgcycle,
                Real a_hang,
                Real a_tolerance,
                int  a_verbosity,
                Real a_normThresh)
{
  CH_TIME("EBCompositeCCProjector::setSolverParams");
  m_macProjector->setSolverParams(a_numSmooths, a_iterMax, a_mgcycle, a_hang, a_tolerance,a_verbosity,a_normThresh);
}
/*****/
void
EBCompositeCCProjector::
setTime(Real a_time)
  {
    m_macProjector->setTime(a_time);
  }
/*****/
int
EBCompositeCCProjector::
project(Vector<LevelData<EBCellFAB>* >&              a_velocity,
        Vector<LevelData<EBCellFAB>* >&              a_gradient,
        const Real&                                  a_gradCoef,
        const Real&                                  a_divuCoef,
        const Vector<LevelData<BaseIVFAB<Real> >* >* a_boundaryVelo)


{
  CH_TIME("EBCompositeCCProjector::project");

  averageVelocityToFaces(m_macVelo, a_velocity);
  bool doDivFreeOutflow = false;
  m_macProjector->enforceVelocityBC(m_macVelo, doDivFreeOutflow);

  const int status = m_macProjector->project(m_macVelo, m_macGrad, a_gradCoef, a_divuCoef, a_boundaryVelo);

  averageFaceToCells(a_gradient, m_macGrad);

  EBAMRDataOps::incr(a_velocity, a_gradient, -1.0);

  return status;
}
/*****/
void
EBCompositeCCProjector::
gradient(Vector<LevelData<EBCellFAB>* >&  a_gradient,
         Vector<LevelData<EBCellFAB>* >&  a_phi)
{
  CH_TIME("EBCompositeCCProjector::gradient");
  m_macProjector->gradient(m_macGrad, a_phi);

  //extrapolate to domain boundaries so that
  //we can have a sufficiently accurate representation
  //of the gradient at domain boundariers
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      m_macGrad[ilev]->exchange();
      ccpExtrapolateToDomainBoundaries(*m_macGrad[ilev],
                                       m_eblg[ilev].getDBL(),
                                       m_eblg[ilev].getEBISL(),
                                       m_eblg[ilev].getDomain(),
                                       m_dx[ilev]);
    }

  averageFaceToCells(a_gradient, m_macGrad);
}
void
EBCompositeCCProjector::
averageVelocityToFaces(Vector<LevelData<EBFluxFAB>* >&  a_macVeloc,
                       Vector<LevelData<EBCellFAB>* >&  a_velocity)
{
  CH_TIME("EBCompositeCCProjector::averageVelocityToFaces");
  //interpolate and then send stuff on through to level function
  // int ncomp = a_velocity[0]->nComp();
  Interval interv(0, SpaceDim-1);
  Vector<LevelData<EBCellFAB> *> amrPhi = m_macProjector->getPhi();
  //  Real time = 0.;
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if (ilev > 0)
        {
          //so it can be reused quadcfi is a one-variable animal
          //when we have ebalias, we can accomplish this without copies.
          //for now, tough luck
          //use phi for scratch space
          LevelData<EBCellFAB>& phiCoar = *amrPhi[ilev-1];
          LevelData<EBCellFAB>& phiFine = *amrPhi[ilev  ];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Interval phiInterv(0, 0);
              Interval velInterv(idir, idir);
              a_velocity[ilev-1]->copyTo(velInterv, phiCoar, phiInterv);
              a_velocity[ilev  ]->copyTo(velInterv, phiFine, phiInterv);
              m_quadCFI[ilev]->interpolate(phiFine,
                                           phiCoar,
                                           phiInterv);
              // m_pwlCFI[ilev]->interpolate(phiFine,
              //                             phiCoar,
              //                             phiCoar,
              //                             time,
              //                             time,
              //                             time,
              //                             phiInterv);

              //on copy back, we need ghost cells, so do the data iterator loop
              for (DataIterator dit = phiFine.dataIterator(); dit.ok(); ++dit)
                {
                  Box region = phiFine[dit()].getRegion(); //includes ghost cells
                  (*a_velocity[ilev])[dit()].copy(region, velInterv, region, phiFine[dit()], phiInterv);
                }

              EBLevelDataOps::setVal(phiFine, 0.0);
              EBLevelDataOps::setVal(phiCoar, 0.0);
            }
        }
      a_velocity[ilev]->exchange(interv);
      ccpAverageVelocityToFaces(*a_macVeloc[ilev],
                                *a_velocity[ilev],
                                m_eblg[ilev].getDBL(),
                                m_eblg[ilev].getEBISL(),
                                m_eblg[ilev].getDomain(),
                                m_dx[ilev],
                                *m_eblg[ilev].getCFIVS());
    }
}
void
EBCompositeCCProjector::
averageVelocityToFaces(Vector<LevelData<EBFluxFAB>* >&  a_macVeloc,
                       Vector<LevelData<EBCellFAB>* >&  a_velocity,
                       const int &                      a_comp)
{
  CH_TIME("EBCompositeCCProjector::averageVelocityToFaces");
  //interpolate and then send stuff on through to level function
  Interval interv(0, SpaceDim-1);
  Vector<LevelData<EBCellFAB> *> amrPhi = m_macProjector->getPhi();
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if (ilev > 0)
        {
          //so it can be reused quadcfi is a one-variable animal
          //when we have ebalias, we can accomplish this without copies.
          //for now, tough luck
          //use phi for scratch space
          LevelData<EBCellFAB>& phiCoar = *amrPhi[ilev-1];
          LevelData<EBCellFAB>& phiFine = *amrPhi[ilev  ];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Interval phiInterv(0, 0);
              Interval velInterv(idir, idir);
              a_velocity[ilev-1]->copyTo(velInterv, phiCoar, phiInterv);
              a_velocity[ilev  ]->copyTo(velInterv, phiFine, phiInterv);
              m_quadCFI[ilev]->interpolate(phiFine,
                                           phiCoar,
                                           phiInterv);

              //on copy back, we need ghost cells, so do the data iterator loop
              for (DataIterator dit = phiFine.dataIterator(); dit.ok(); ++dit)
                {
                  Box region = phiFine[dit()].getRegion(); //includes ghost cells
                  (*a_velocity[ilev])[dit()].copy(region, velInterv, region, phiFine[dit()], phiInterv);
                }

              EBLevelDataOps::setVal(phiFine, 0.0);
              EBLevelDataOps::setVal(phiCoar, 0.0);
            }
        }
      a_velocity[ilev]->exchange(interv);
      ccpAverageVelocityToFaces(*a_macVeloc[ilev],
                                *a_velocity[ilev],
                                m_eblg[ilev].getDBL(),
                                m_eblg[ilev].getEBISL(),
                                m_eblg[ilev].getDomain(),
                                m_dx[ilev],
                                *m_eblg[ilev].getCFIVS(),
                                a_comp);
    }

  EBAMRDataOps::averageDown(a_macVeloc, m_eblg, m_refRat);
}
void
EBCompositeCCProjector::
averageStressToFaces(Vector<LevelData<EBFluxFAB>* >&  a_macVeloc,
                     Vector<LevelData<EBCellFAB>* >&  a_velocity)
{
  CH_TIME("EBCompositeCCProjector::averageVelocityToFaces");
  //interpolate and then send stuff on through to level function
  int ncomp = a_velocity[0]->nComp();
  Interval interv(0, ncomp-1);
  Vector<LevelData<EBCellFAB> *> amrPhi = m_macProjector->getPhi();
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if (ilev > 0)
        {
          //so it can be reused quadcfi is a one-variable animal
          //when we have ebalias, we can accomplish this without copies.
          //for now, tough luck
          //use phi for scratch space
          LevelData<EBCellFAB>& phiCoar = *amrPhi[ilev-1];
          LevelData<EBCellFAB>& phiFine = *amrPhi[ilev  ];
          for (int idir = 0; idir < ncomp; idir++)
            {
              Interval phiInterv(0, 0);
              Interval velInterv(idir, idir);
              a_velocity[ilev-1]->copyTo(velInterv, phiCoar, phiInterv);
              a_velocity[ilev  ]->copyTo(velInterv, phiFine, phiInterv);
              m_quadCFI[ilev]->interpolate(phiFine,
                                           phiCoar,
                                           phiInterv);

              //on copy back, we need ghost cells, so do the data iterator loop
              for (DataIterator dit = phiFine.dataIterator(); dit.ok(); ++dit)
                {
                  Box region = phiFine[dit()].getRegion(); //includes ghost cells
                  (*a_velocity[ilev])[dit()].copy(region, velInterv, region, phiFine[dit()], phiInterv);
                }

              EBLevelDataOps::setVal(phiFine, 0.0);
              EBLevelDataOps::setVal(phiCoar, 0.0);
            }
        }
      a_velocity[ilev]->exchange(interv);
      ccpAverageStressToFaces(*a_macVeloc[ilev],
                              *a_velocity[ilev],
                              m_eblg[ilev].getDBL(),
                              m_eblg[ilev].getEBISL(),
                              m_eblg[ilev].getDomain(),
                              m_dx[ilev],
                              *m_eblg[ilev].getCFIVS());
    }
}
/*****/
void
EBCompositeCCProjector::
averageFaceToCells(Vector<LevelData<EBCellFAB>* >&         a_cellData,
                   const Vector<LevelData<EBFluxFAB>* >&   a_faceData)
{
  CH_TIME("EBCompositeCCProjector::averageGradientToCells");
  //average the face gradients to cell centers
  Interval interv(0,0);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      LevelData<EBFluxFAB>& faceData = (LevelData<EBFluxFAB>&) *a_faceData[ilev];
      faceData.exchange(interv);
      ccpExtrapolateToDomainBoundaries(faceData,
                                       m_eblg[ilev].getDBL(),
                                       m_eblg[ilev].getEBISL(),
                                       m_eblg[ilev].getDomain(),
                                       m_dx[ilev]);
    }
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      ccpAverageFaceToCells(*a_cellData[ilev],
                            *a_faceData[ilev],
                            m_eblg[ilev].getDBL(),
                            m_eblg[ilev].getEBISL(),
                            m_eblg[ilev].getDomain(),
                            m_dx[ilev]);
    }
}
/*****/
void
EBCompositeCCProjector::
kappaDivergence(Vector<LevelData<EBCellFAB>* >&              a_divu,
                Vector<LevelData<EBCellFAB>* >&              a_velo,
                const Vector<LevelData<BaseIVFAB<Real> >* >* a_boundaryVelo)
{
  CH_TIME("EBCompositeCCProjector::kappaDivergence");
  averageVelocityToFaces(m_macVelo, a_velo);
  bool doDivFreeOutflow = false;
  m_macProjector->enforceVelocityBC(m_macVelo, doDivFreeOutflow);
  m_macProjector->kappaDivergence(a_divu, m_macVelo, a_boundaryVelo);
}
/*****/

#include "NamespaceFooter.H"
