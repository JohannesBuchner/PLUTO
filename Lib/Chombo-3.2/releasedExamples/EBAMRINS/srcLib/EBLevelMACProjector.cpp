#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MACProjectorF_F.H"
#include "EBLevelMACProjector.H"
#include "NeumannPoissonEBBC.H"
#include "NeumannPoissonDomainBC.H"
#include "EBFluxFactory.H"
#include "EBAMRPoissonOp.H"
#include "EBArith.H"
#include "FaceIterator.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>

#include "NamespaceHeader.H"

bool EBLevelMACProjector::s_verbose  = false;
int  EBLevelMACProjector::s_curLevel = -1;
string EBLevelMACProjector::s_debugString = string("debug:");

/****/
void
EBLevelMACProjector::
setDebugString(string a_string)
{
  s_debugString = a_string;
}
string
EBLevelMACProjector::
getDebugString()
{
  return s_debugString;
}
/****/
void
EBLevelMACProjector::
setVerbose(bool a_verbose)
{
  s_verbose = a_verbose;
}
/*****/
bool
EBLevelMACProjector::
getVerbose()
{
  return s_verbose;
}
/*****/
int
EBLevelMACProjector::
getCurLevel()
{
  return s_curLevel;
}
/*****/
void
EBLevelMACProjector::
setCurLevel(int a_curLevel)
{
  s_curLevel = a_curLevel;
}
/*****/

/*****/
EBLevelMACProjector::
EBLevelMACProjector(const DisjointBoxLayout&                        a_grids,
                    const EBISLayout&                               a_ebisl,
                    const ProblemDomain&                            a_domain,
                    const RealVect&                                 a_dx,
                    const RealVect&                                 a_origin,
                    const LinearSolver<LevelData<EBCellFAB> >&      a_bottomSolver,
                    const RefCountedPtr<BaseEBBCFactory>&           a_ebbcPhi,
                    const RefCountedPtr<BaseDomainBCFactory>&       a_domainBCFactPhi,
                    const RefCountedPtr<BaseDomainBCFactory>&       a_domainBCFactVel,
                    const int  &                                    a_numSmooths,
                    const int  &                                    a_mgCycle,
                    const int  &                                    a_maxIterations,
                    const Real &                                    a_tolerance,
                    const int  &                                    a_maxDepth,
                    const Real&                                     a_time,
                    const IntVect&                                  a_nghostPhi,
                    const IntVect&                                  a_nghostRhs)
{
  m_grids  = a_grids;
  m_ebisl  = a_ebisl;
  m_domain = a_domain;
  m_dx     = a_dx;
  m_origin = a_origin;

  m_ebbcPhi = a_ebbcPhi;
  m_domainBCFactPhi = a_domainBCFactPhi;
  m_domainBCFactVel = a_domainBCFactVel;

  int numPreCondIters = 4;
  int mgRelaxType = 1;
  Vector<EBLevelGrid> gridsVec(1, EBLevelGrid(m_grids, m_ebisl, m_domain));
  Vector<int> refRat(1,2);
  Real alpha = 0.0;
  Real beta = 1.0;
  m_nghostPhi = a_nghostPhi;
  m_nghostRhs = a_nghostRhs;

  Vector<RefCountedPtr<EBQuadCFInterp> >
    quadcfi(1, RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp()));
  m_levelOpFactory = new EBAMRPoissonOpFactory(gridsVec, refRat, quadcfi, m_dx, m_origin,
                                               numPreCondIters, mgRelaxType, m_domainBCFactPhi, m_ebbcPhi,
                                               alpha, beta, a_time, a_nghostPhi, a_nghostRhs);
  m_solver.define(*m_levelOpFactory, (LinearSolver<LevelData<EBCellFAB> >*)(&a_bottomSolver), m_domain,a_maxDepth);
  m_solver.m_homogeneous = true;

  m_maxIterations = a_maxIterations;
  m_tolerance     = a_tolerance;
  m_solver.m_pre   = a_numSmooths;
  m_solver.m_post  = a_numSmooths;
  m_solver.m_cycle = a_mgCycle;

}
/*****/
void
EBLevelMACProjector::
project(LevelData<EBFluxFAB>&  a_velocity,
        LevelData<EBFluxFAB>&  a_gradPhi,
        const LevelData<BaseIVFAB<Real> >* const a_boundaryVelocity)
{
  CH_assert(a_velocity.nComp() == 1);
  CH_assert(a_gradPhi.nComp() == 1);

  EBFluxFactory fluxfact(m_ebisl);
  EBCellFactory cellfact(m_ebisl);
  LevelData<EBFluxFAB> centroidVelocity(m_grids, 1, IntVect::Unit, fluxfact);

  Interval interv(0, 0);
  LevelData<EBFluxFAB>& velExch = (LevelData<EBFluxFAB>&)(a_velocity);
  velExch.exchange(interv);
  //interpolate velocity to face centroids
  EBArith::interpolateFluxToCentroids(centroidVelocity, a_velocity,
                                      m_grids, m_ebisl, m_domain);

  LevelData<EBCellFAB> phi(m_grids, 1, m_nghostPhi, cellfact);
  LevelData<EBCellFAB> rhs(m_grids, 1, m_nghostRhs, cellfact);

  //set rhs = kappa*div(u)
  macKappaDivergence(rhs, centroidVelocity, m_grids, m_ebisl, m_domain, m_dx, a_boundaryVelocity);

  solve(phi, rhs);

  //compute gradient and subtract off of velocity
  macGradient(a_gradPhi, phi,  m_grids, m_ebisl, m_domain, m_dx);
  //enforce bc's on gradient
  macEnforceGradientBC(a_gradPhi,phi,m_grids,m_ebisl,m_domain,m_dx,m_origin,1.e99,m_domainBCFactPhi);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_velocity[dit()][idir] -= a_gradPhi[dit()][idir];
        }
    }
}
/*****/
void
EBLevelMACProjector::
solve(LevelData<EBCellFAB>&       a_phi,
      const LevelData<EBCellFAB>& a_rhs)
{
  pout() << " EBLevelMACProjector::solve" << endl;
  m_solver.solve(a_phi, a_rhs, m_tolerance, m_maxIterations, 3);
}
/*****/
EBLevelMACProjector::
~EBLevelMACProjector()
{
  if (m_levelOpFactory != NULL)
    {
      delete m_levelOpFactory;
      m_levelOpFactory = NULL;
    }
}
/*****/
void
macGradient(LevelData<EBFluxFAB>&        a_gradPhi,
            const LevelData<EBCellFAB>&  a_phi,
            const DisjointBoxLayout&     a_grids,
            const EBISLayout&            a_ebisl,
            const ProblemDomain&         a_domain,
            const RealVect&              a_dx)
{
  Interval interv(0,0);
  //cast to non-const so we can do exchange.
  //none of the valid data is changed so const-nitude
  //is preserved at least in spirit.
  LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&) a_phi;
  phi.exchange(interv);
  int ibox = 0;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {

      macGradient(a_gradPhi[dit()], a_phi[dit()],
                  a_ebisl[dit()], a_grids.get(dit()),
                  a_domain, a_dx);

      ibox++;
    }
}
/****/
void
macGradient(EBFluxFAB&           a_gradPhi,
            const EBCellFAB&     a_phi,
            const EBISBox&       a_ebisBox,
            const Box&           a_box,
            const ProblemDomain& a_domain,
            const RealVect&      a_dx)
{
  a_gradPhi.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      EBFaceFAB&      gradPhiFAB = a_gradPhi[idir];

      //do interior faces.  domain boundary faces are boundary conditions
      Box interiorBox = a_box;
      interiorBox.grow(idir, 1);
      interiorBox &= a_domain;
      interiorBox.grow(idir, -1);

      Box faceBox = surroundingNodes(interiorBox, idir);
      BaseFab<Real>&   regGradPhi = gradPhiFAB.getSingleValuedFAB();
      const BaseFab<Real>& regPhi = a_phi.getSingleValuedFAB();

      FORT_MACGRADPHI(CHF_FRA1(regGradPhi, 0),
                      CHF_CONST_FRA1(regPhi, 0),
                      CHF_CONST_INT(idir),
                      CHF_CONST_REAL(a_dx[idir]),
                      CHF_BOX(faceBox));

      //boundaries are not included
      FaceStop::WhichFaces stopCrit;
      if (a_domain.isPeriodic(idir))
        {
          stopCrit = FaceStop::SurroundingWithBoundary;
        }
      else
        {
          stopCrit = FaceStop::SurroundingNoBoundary;
        }
      IntVectSet irregIVS = a_ebisBox.getIrregIVS(interiorBox);
      for (FaceIterator faceit(irregIVS, a_ebisBox.getEBGraph(), idir, stopCrit);
          faceit.ok(); ++faceit)
        {
          const FaceIndex& face = faceit();
          Real faceGrad = 0.;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              faceGrad += sign(sit())*a_phi(face.getVoF(sit()), 0);
            }
          gradPhiFAB(face, 0) = faceGrad/a_dx[idir];
        }
    }
}
/*****/
void
macEnforceGradientBC(LevelData<EBFluxFAB>&              a_gradPhi,
                     const LevelData<EBCellFAB>&        a_phi,
                     const DisjointBoxLayout&           a_grids,
                     const EBISLayout&                  a_ebisl,
                     const ProblemDomain&               a_domain,
                     const RealVect&                    a_dx,
                     const RealVect&                    a_origin,
                     const Real&                        a_time,
                     RefCountedPtr<BaseDomainBCFactory> a_domainBCFactPhi)
{
  BaseDomainBC* domainBCPhi = a_domainBCFactPhi->create(a_domain, a_ebisl, a_dx);
  Interval interv(0,0);
  //cast to non-const so we can do exchange.
  //none of the valid data is changed so const-nitude
  //is preserved at least in spirit.
  LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>&) a_phi;
  phi.exchange(interv);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      EBFluxFAB&       gradPhi = a_gradPhi[dit()];
      const EBCellFAB&     phi = a_phi[dit()];
      const Box&        dblBox = a_grids.get(dit());
      const EBISBox&   ebisBox = a_ebisl[dit()];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (!a_domain.isPeriodic(idir))
            {
              EBFaceFAB& gradPhiFAB = gradPhi[idir];
              //boundaries only
              FaceStop::WhichFaces stopCrit = FaceStop::AllBoundaryOnly;
              IntVectSet ivsBox = IntVectSet(dblBox);
              for (FaceIterator faceit(ivsBox, ebisBox.getEBGraph(), idir, stopCrit);
                  faceit.ok(); ++faceit)
                {
                  const FaceIndex& face = faceit();
                  //figure out where the inside of the domain is
                  Side::LoHiSide outSide = Side::Hi;
                  if (!a_domain.contains(face.gridIndex(Side::Lo)))
                    {
                      outSide = Side::Lo;
                    }
                  const int comp = 0;
                  Real flux;
                  domainBCPhi->getFaceGradPhi(flux,face,comp,phi,a_origin,
                                              a_dx,idir,outSide,dit(),a_time,false,RealVect::Zero,false);
                  gradPhiFAB(face, 0) = flux;
                }
            }
        }
    }
  delete domainBCPhi;
}
/****/
void
macKappaDivergence(LevelData<EBCellFAB>&        a_divVel,
                   const LevelData<EBFluxFAB>&  a_velocity,
                   const DisjointBoxLayout&     a_grids,
                   const EBISLayout&            a_ebisl,
                   const ProblemDomain&         a_domain,
                   const RealVect&              a_dx,
                   const LevelData<BaseIVFAB<Real> >*  a_boundaryVel )
{
  CH_TIME("EBLevelMACProjector::macKappaDivergence(level)");
  int ibox = 0;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const BaseIVFAB<Real>* boundaryVelFAB = NULL;
      if (a_boundaryVel != NULL)
        {
          boundaryVelFAB = &((*a_boundaryVel)[dit()]);
        }

      macKappaDivergence(a_divVel[dit()], a_velocity[dit()],
                         a_ebisl[dit()],  a_grids.get(dit()),
                         a_domain, a_dx, boundaryVelFAB);

      ibox++;
    }
}
/*****/
void
macKappaDivergence(EBCellFAB&             a_divVel,
                   const EBFluxFAB&       a_velo,
                   const EBISBox&         a_ebisBox,
                   const Box&             a_box,
                   const ProblemDomain&   a_domain,
                   const RealVect&        a_dx,
                   const BaseIVFAB<Real>* a_boundaryVel)
{
  //set the divergence initially to zero
  //then loop through directions and increment the divergence
  //with each directions flux difference.
  a_divVel.setVal(0.0);
  BaseFab<Real>&       regDivF = a_divVel.getSingleValuedFAB();
  regDivF.setVal(0.);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      /* do the regular vofs */
      const EBFaceFAB& veloDir = a_velo[idir];
      const BaseFab<Real>& regVeloDir = veloDir.getSingleValuedFAB();
      int ncons = 1;
      FORT_MACDIVERGEF( CHF_BOX(a_box),
                        CHF_FRA(regDivF),
                        CHF_CONST_FRA(regVeloDir),
                        CHF_CONST_INT(idir),
                        CHF_CONST_INT(ncons),
                        CHF_CONST_REAL(a_dx[idir]));
    }

  //update the irregular vofs using conservative diff
  IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(a_box);
  for (VoFIterator vofit(ivsIrreg, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Real update = 0.;
      for ( int idir = 0; idir < SpaceDim; idir++)
        {
          const EBFaceFAB& veloDir = a_velo[idir];
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int isign = sign(sit());
              Vector<FaceIndex> faces =
                a_ebisBox.getFaces(vof, idir, sit());
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  const FaceIndex& face = faces[iface];
                  Real areaFrac = a_ebisBox.areaFrac(face);
                  Real faceVelo = veloDir(face, 0);
                  update += isign*areaFrac*faceVelo/a_dx[idir];
                }
            }
        }
      //EB Flux assumed zero if (a_boundaryVel == NULL)
      if (a_boundaryVel != NULL)
        {
          RealVect normal = a_ebisBox.normal(vof);
          Real boundaryArea = a_ebisBox.bndryArea(vof);
          Real boundaryFlux = 0.0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Real boundaryVel = (*a_boundaryVel)(vof, idir);
              boundaryFlux += boundaryArea*normal[idir]*boundaryVel/a_dx[idir];
            }
          update -= boundaryFlux;
        }

      a_divVel(vof, 0) = update;
    }
  a_divVel.setCoveredCellVal(0.0, 0);
}
void
macEnforceVelocityBC(LevelData<EBFluxFAB>&              a_velocity,
                     const DisjointBoxLayout&           a_grids,
                     const EBISLayout&                  a_ebisl,
                     const ProblemDomain&               a_domain,
                     const RealVect&                    a_dx,
                     const RealVect&                    a_origin,
                     RefCountedPtr<BaseDomainBCFactory> a_domainBCFactVel,
                     const Real&                        a_time,
                     const bool&                        a_doDivFreeOutflow,
                     const LevelData<BaseIVFAB<Real> >* const a_boundaryVelocity)

{
  CH_TIME("EBLevelMACProjector::macEnforceVelocityBC");
  BaseDomainBC* domainBC = a_domainBCFactVel->create(ProblemDomain(a_domain), a_ebisl, a_dx);
  domainBC->enforceFaceVel(a_velocity,a_grids,a_ebisl,a_domain,a_dx,a_time,a_origin,a_doDivFreeOutflow);
  delete domainBC;
}
void
macEnforceVelocityBC(LevelData<EBFluxFAB>&              a_velocity,
                     const DisjointBoxLayout&           a_grids,
                     const EBISLayout&                  a_ebisl,
                     const ProblemDomain&               a_domain,
                     const RealVect&                    a_dx,
                     const RealVect&                    a_origin,
                     RefCountedPtr<BaseDomainBCFactory> a_domainBCFactVel,
                     const Real&                        a_time,
                     const bool&                        a_doDivFreeOutflow,
                     const int&                         a_comp,
                     const LevelData<BaseIVFAB<Real> >* const a_boundaryVelocity)

{
  CH_TIME("EBLevelMACProjector::macEnforceVelocityBC");
  BaseDomainBC* domainBC = a_domainBCFactVel->create(ProblemDomain(a_domain), a_ebisl, a_dx);
  domainBC->enforceFaceVel(a_velocity,a_grids,a_ebisl,a_domain,a_dx,a_time,a_origin,a_doDivFreeOutflow,a_comp);
  delete domainBC;
}

#include "NamespaceFooter.H"
