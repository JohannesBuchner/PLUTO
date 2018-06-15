#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBCompositeMACProjector.H"
#include "NeumannPoissonEBBC.H"
#include "NeumannPoissonDomainBC.H"
#include "EBFluxFactory.H"
#include "EBArith.H"
#include "EBAMRDataOps.H"
#include "FaceIterator.H"
#include "MACProjectorF_F.H"

#include "NamespaceHeader.H"

/*********/
EBCompositeMACProjector::
EBCompositeMACProjector(const Vector<EBLevelGrid>&                      a_eblg,
                        const Vector<int>&                              a_refRat,
                        const Vector<RefCountedPtr<EBQuadCFInterp> >&   a_quadCFI,
                        const RealVect&                                 a_coarsestDx,
                        const RealVect&                                 a_origin,
                        const RefCountedPtr<BaseDomainBCFactory>&       a_baseDomainBCVel,
                        const RefCountedPtr<BaseDomainBCFactory>&       a_baseDomainBCPhi,
                        const RefCountedPtr<BaseEBBCFactory>&           a_baseEBBCPhi,
                        const bool&                                     a_subtractOffMean,
                        const int&                                      a_numLevels,
                        const int&                                      a_verbosity,
                        const int&                                      a_numPreCondIters,
                        const Real&                                     a_time,
                        const int&                                      a_relaxType,
                        const int&                                      a_bottomSolverType)
{
  CH_TIME("EBCompositeMACProjector::EBCompositeMACProjector");
  m_eblg = a_eblg;
  m_quadCFI = a_quadCFI;
  if (a_numLevels < 0)
    {
      m_numLevels = a_eblg.size();
    }
  else
    {
      CH_assert(a_numLevels <= a_eblg.size());
      m_numLevels = a_numLevels;
    }
  //a few consistency checks
  CH_assert(a_eblg.size() >= m_numLevels);
  CH_assert(a_refRat.size() >= m_numLevels);

  m_baseDomainBCFactVel  = a_baseDomainBCVel;
  m_baseDomainBCFactPhi  = a_baseDomainBCPhi;
  m_refRat     = a_refRat;

  m_dx.resize(m_numLevels);
  m_dx[0]     = a_coarsestDx;
  m_origin    = a_origin;
  m_time      = a_time;

  for (int ilev = 1; ilev < m_numLevels; ilev++)
    {
      m_dx[ilev]     = m_dx[ilev-1]/m_refRat[ilev-1];
    }

  const EBIndexSpace* ebisPtr = a_eblg[0].getEBIS();
  m_divu.resize(m_numLevels, NULL);
  m_phi.resize(m_numLevels, NULL);
  m_fluxReg.resize(m_numLevels, NULL);
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      EBCellFactory ebcellfact(m_eblg[ilev].getEBISL());
      m_divu[ilev]    = new LevelData<EBCellFAB>(m_eblg[ilev].getDBL(), 1, IntVect::Zero, ebcellfact);
      m_phi[ilev]     = new LevelData<EBCellFAB>(m_eblg[ilev].getDBL(), 1, IntVect::Unit, ebcellfact);
      m_fluxReg[ilev] = new EBFluxRegister();

      if (ilev < (m_numLevels-1))
        {
          //flux register lives with coarser level
          m_fluxReg[ilev]->define(m_eblg[ilev+1].getDBL(),   m_eblg[ilev].getDBL(),
                                  m_eblg[ilev+1].getEBISL(), m_eblg[ilev].getEBISL(),
                                  m_eblg[ilev].getDomain(), m_refRat[ilev], 1, ebisPtr);
        }
    }

  //default values
  int numSmooth = 4;
  int itermax =  100;
  int mgCycle = 1;
  Real hang = 1.0e-20;
  Real tolerance = 1.0e-12;
  Real normThresh= 1.0e-12;
  m_useInitialGuess = false;
  m_subtractOffMean = a_subtractOffMean;
  m_bottomSolverType = a_bottomSolverType;

  ProblemDomain coarsestDomain(m_eblg[0].getDomain());

  Real alpha = 0.0;
  Real beta = 1.0;
  m_opFactory = new EBAMRPoissonOpFactory(m_eblg, m_refRat, m_quadCFI, m_dx[0], m_origin,
                                          a_numPreCondIters,a_relaxType,
                                          a_baseDomainBCPhi, a_baseEBBCPhi, alpha, beta, m_time,
                                          IntVect::Unit, IntVect::Zero, m_numLevels);
  if (m_bottomSolverType==0)
    {//BiCGStab bottom solver
      m_solver.define(coarsestDomain, *m_opFactory, &m_bottomSolverBiCG, m_numLevels);
    }
  else if (m_bottomSolverType==1)
    {//EBSimpleSolver bottom solver
      m_solver.define(coarsestDomain, *m_opFactory, &m_bottomSolverSimp, m_numLevels);
    }
  else if (m_bottomSolverType==2)
    {//GMRES bottom solver
      m_solver.define(coarsestDomain, *m_opFactory, &m_bottomSolverGMRes, m_numLevels);
    }
  else
    {
      MayDay::Error("EBCompositeMACProjector::bad m_bottomSolverType type");
    }

  setSolverParams(numSmooth,  itermax, mgCycle, hang, tolerance, a_verbosity, normThresh);

  m_solver.init(m_phi, m_divu, m_numLevels-1, 0);
}
/*****/
void
EBCompositeMACProjector::
setSolverParams(int  a_numSmooths,
                int  a_iterMax,
                int  a_mgcycle,
                Real a_hang,
                Real a_tolerance,
                int  a_verbosity,
                Real a_normThresh)
{
  CH_TIME("EBCompositeMACProjector::setSolverParams");
  m_solver.setSolverParameters(a_numSmooths,
                               a_numSmooths,
                               a_numSmooths,
                               a_mgcycle,
                               a_iterMax,
                               a_tolerance,
                               a_hang,
                               a_normThresh);

  m_solver.m_verbosity = a_verbosity;

  if (m_bottomSolverType==0)
    {
      m_bottomSolverBiCG.m_verbosity = a_verbosity - 2;
      m_bottomSolverBiCG.m_hang = -0.1;//because BiCGStab can be a bouncy ride
    }
  else if (m_bottomSolverType==1)
    {
      m_bottomSolverSimp.setNumSmooths(128*a_numSmooths);
    }
  else if (m_bottomSolverType==2)
    {
      //GMRES
    }
  else
    {
      MayDay::Error("EBCompositeMACProjector::setSolverParams bad bottomSolverType");
    }
}
/*****/
void
EBCompositeMACProjector::
setTime(Real a_time)
  {
    m_time = a_time;
    EBAMRPoissonOp::setOperatorTime(m_time);
  }
/*****/
void
EBCompositeMACProjector::
setInitialPhi(Vector<LevelData<EBCellFAB>* >&  a_phi)
{
  CH_TIME("EBCompositeMACProjector::setInitialPhi");
  m_useInitialGuess = true;
  EBAMRDataOps::assign(m_phi,a_phi);
}
/*****/
void
EBCompositeMACProjector::
initialize(Vector<LevelData<EBCellFAB>* >&  a_phi)
{
  CH_TIME("EBCompositeMACProjector::initialize");
  if (!m_useInitialGuess)
    {
      EBAMRDataOps::setVal(m_phi,0.0);
    }
  m_useInitialGuess = false;
}
/*****/
int
EBCompositeMACProjector::
project(Vector<LevelData<EBFluxFAB>* >&              a_velocity,
        Vector<LevelData<EBFluxFAB>* >&              a_gradient,
        const Real&                                  a_gradCoef,
        const Real&                                  a_divuCoef,
        const Vector<LevelData<BaseIVFAB<Real> >* >* a_boundaryVelo)

{
  CH_TIME("EBCompositeMACProjector::project");
  //compute divergence everywhere.  This includes refluxing at coarse-fine boundaries
  kappaDivergence(m_divu, a_velocity, a_boundaryVelo);

  EBAMRDataOps::scale(m_divu,a_divuCoef);

  //make an initial guess
  initialize(m_phi);

  //solve laplacian phi = div u
  //zeroPhi=false is inconsequential in solveNoInit because of initialize above
  bool zeroPhi = false;
  m_solver.solveNoInit(m_phi, m_divu, m_numLevels-1, 0, zeroPhi);

  //compute grad phi.  This includes replacing coarse grad phi
  //with average of fine at coarse-fine boundaries
  gradient(a_gradient, m_phi);
  enforceGradientBC(a_gradient,m_phi);

  EBAMRDataOps::scale(a_gradient,a_gradCoef);

  //average the gradient down to coarser levels
  EBAMRDataOps::averageDown(a_gradient,m_eblg,m_refRat);

  if (m_subtractOffMean)
    {
      //Subtract off mean of pressure so that convergence tests can make sense.
      //This is also useful if you are using initial guesses for m_phi with all Neumann bcs
      EBAMRDataOps::subtractOffMean(m_phi, m_eblg, m_refRat);
    }

  //u := u - grad phi
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      for (DataIterator dit = m_eblg[ilev].getDBL().dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              EBFaceFAB&         velFAB = (*a_velocity[ilev])[dit()][idir];
              const EBFaceFAB&  gradFAB = (*a_gradient[ilev])[dit()][idir];
              velFAB -= gradFAB;
            }
        }
    }

  return m_solver.m_exitStatus;
}

/***/
void
EBCompositeMACProjector::
correctVelocityComponent(Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >      &  a_coveredVelLo,
                         Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >      &  a_coveredVelHi,
                         const Vector< LayoutData< Vector< Vector<VolIndex> > >* >&  a_coveredFaceLo,
                         const Vector< LayoutData< Vector< Vector<VolIndex> > >* >&  a_coveredFaceHi,
                         const Vector< LayoutData< Vector< IntVectSet > >* >      &  a_coveredSetsLo,
                         const Vector< LayoutData< Vector< IntVectSet > >* >      &  a_coveredSetsHi,
                         const Vector<LevelData<EBFluxFAB>* >                     &  a_macGradient,
                         int                                                         a_coveredFaceDir,
                         int                                                         a_velComp)
{
  CH_TIME("EBCompositeMACProjector::correctVelocityComponent");
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          LayoutData<Vector<BaseIVFAB<Real>* > >*        velPtr = NULL;
          const LayoutData<Vector<Vector<VolIndex> > >* facePtr = NULL;
          const LayoutData<Vector<IntVectSet>  >*       setsPtr = NULL;
          if (sit() == Side::Lo)
            {
              velPtr  =  a_coveredVelLo[ilev];
              facePtr = a_coveredFaceLo[ilev];
              setsPtr = a_coveredSetsLo[ilev];
            }
          else
            {
              velPtr  =  a_coveredVelHi[ilev];
              facePtr = a_coveredFaceHi[ilev];
              setsPtr = a_coveredSetsHi[ilev];
            }

          const LevelData<EBFluxFAB>& macGradient = *a_macGradient[ilev];
          for (DataIterator dit = m_eblg[ilev].getDBL().dataIterator(); dit.ok(); ++dit)
            {
              const EBFluxFAB       & macGradFAB  = macGradient[dit()];
              const Vector<VolIndex>& coveredFace =  (*facePtr)[dit()][a_coveredFaceDir];
              const IntVectSet      & coveredSets =  (*setsPtr)[dit()][a_coveredFaceDir];
              BaseIVFAB<Real>       & coveredVel  = *((*velPtr)[dit()][a_coveredFaceDir]);
              const EBISBox& ebisBox = m_eblg[ilev].getEBISL()[dit()];

              correctVelocityComponent(coveredVel, coveredFace, coveredSets, macGradFAB, ebisBox,
                                       a_coveredFaceDir, sit(), a_velComp);
            }
        }
    }
}
/***/
Real
getAverageFaceGrad(const VolIndex&   a_vof,
                   const EBFaceFAB&  a_macGradient,
                   const EBISBox&    a_ebisBox,
                   int               a_faceDir)
{
  CH_TIME("EBCompositeMACProjector::getAverageFaceGrad");
  Vector<FaceIndex> hiFaces = a_ebisBox.getFaces(a_vof, a_faceDir, Side::Hi);
  Vector<FaceIndex> loFaces = a_ebisBox.getFaces(a_vof, a_faceDir, Side::Lo);
  bool hasHiFaces = (hiFaces.size() > 0);
  bool hasLoFaces = (loFaces.size() > 0);
  Real hiVal = 0.0;
  Real loVal = 0.0;
  if (hasHiFaces)
    {
      for (int iface = 0; iface < hiFaces.size(); iface++)
        {
          hiVal += a_macGradient(hiFaces[iface], 0);
        }
      hiVal /= hiFaces.size();

    }
  if (hasLoFaces)
    {
      for (int iface = 0; iface < loFaces.size(); iface++)
        {
          loVal += a_macGradient(loFaces[iface], 0);
        }
      loVal /= loFaces.size();
    }
  if (!hasLoFaces && hasHiFaces)
    {
      //only want to do the wacky extrapolation thing if  not multivalued
      if ((hiFaces.size() == 1) && (!hiFaces[0].isBoundary()))
        {
          const VolIndex& hiVoF = hiFaces[0].getVoF(Side::Hi);
          Vector<FaceIndex> farFaces = a_ebisBox.getFaces(hiVoF, a_faceDir, Side::Hi);
          if (farFaces.size() > 0)
            {
              Real farVal = 0.0;
              for (int iface = 0; iface < farFaces.size(); iface++)
                {
                  farVal += a_macGradient(farFaces[iface], 0);
                }
              farVal /= farFaces.size();
              loVal = 2*hiVal - farVal;
            }
          else
            {
              loVal = hiVal;
            }
        }
      else
        {
          loVal = hiVal;
        }
    }


  if (!hasHiFaces && hasLoFaces)
    {
      //only want to do the wacky extrapolation thing if not multivalued
      if ((loFaces.size() == 1) && (!loFaces[0].isBoundary()))
        {
          const VolIndex& hiVoF = loFaces[0].getVoF(Side::Lo);
          Vector<FaceIndex> farFaces = a_ebisBox.getFaces(hiVoF, a_faceDir, Side::Lo);
          if (farFaces.size() > 0)
            {
              Real farVal = 0.0;
              for (int iface = 0; iface < farFaces.size(); iface++)
                {
                  farVal += a_macGradient(farFaces[iface], 0);
                }
              farVal /= farFaces.size();
              hiVal = 2*loVal - farVal;
            }
          else
            {
              hiVal = loVal;
            }
        }
      else
        {
          hiVal = loVal;
        }
    }

  Real retval = 0.5*(hiVal + loVal);
  return retval;
}
/***/
void
EBCompositeMACProjector::
correctVelocityComponent(BaseIVFAB<Real>       &  a_coveredVel,
                         const Vector<VolIndex>&  a_coveredFace,
                         const IntVectSet      &  a_coveredSets,
                         const EBFluxFAB       &  a_macGradient,
                         const EBISBox         &  a_ebisBox,
                         int    a_coveredFaceDir,  Side::LoHiSide a_sd,    int a_faceGradComp)
{
  CH_TIME("EBCompositeMACProjector::correctVelocityComponent");
  //if a_coveredFaceDir  == faceGradComp, just linearly extrapolate gradient and subtract.
  //if a_coveredFaceDir != faceGradComp, need to average to cell centers then extrapolate and subtract
  CH_assert(a_coveredVel.nComp() == 1);
  CH_assert(a_macGradient.nComp() == 1);
  const EBFaceFAB& macGradFAB = a_macGradient[a_faceGradComp];
  for (int ivof = 0; ivof < a_coveredFace.size(); ivof++)
    {
      const VolIndex& vof = a_coveredFace[ivof];
      Vector<FaceIndex> nearFaces, farFaces;
      bool hasNearFaces = false;
      bool hasFarFaces  = false;
      VolIndex  nearVoF, farVoF;
      FaceIndex nearFace, farFace;

      nearFaces = a_ebisBox.getFaces(vof, a_coveredFaceDir, flip(a_sd));
      hasNearFaces = (nearFaces.size()==1) && (!nearFaces[0].isBoundary());
      if (hasNearFaces)
        {
          nearFace = nearFaces[0];
          nearVoF = nearFace.getVoF(flip(a_sd));

          farFaces = a_ebisBox.getFaces(nearVoF, a_coveredFaceDir, flip(a_sd));
          hasFarFaces = (farFaces.size()==1)  && (!farFaces[0].isBoundary());
          if (hasFarFaces)
            {
              farFace = farFaces[0];
              farVoF  = farFace.getVoF(flip(a_sd));
            }
        }
      Real extrapGrad;
      if (a_coveredFaceDir == a_faceGradComp)
        {
          //if a_coveredFaceDir  == faceGradComp, just linearly extrapolate gradient and subtract.
          if (hasNearFaces && hasFarFaces)
            {
              Real nearGrad  = macGradFAB(nearFace, 0);
              Real farGrad   = macGradFAB(farFace,  0);
              extrapGrad = 2.0*nearGrad - farGrad;
            }
          else if (hasNearFaces)
            {
              extrapGrad  = macGradFAB(nearFace, 0);
            }
          else
            {
              extrapGrad = 0.0;
            }

        }
      else
        {
          //if a_coveredFaceDir != faceGradComp, need to average to cell centers then extrapolate and subtract
          Real nearGrad = getAverageFaceGrad(vof,     macGradFAB, a_ebisBox, a_faceGradComp);
          if (hasNearFaces)
            {
              Real farGrad =  getAverageFaceGrad(nearVoF, macGradFAB, a_ebisBox, a_faceGradComp);
              extrapGrad = 1.5*nearGrad - 0.5*farGrad;
            }
          else
            {
              extrapGrad = nearGrad;
            }
        }
      a_coveredVel(vof, 0) -= extrapGrad;
    }
}
/***/
void
EBCompositeMACProjector::
correctVelocityComponent(Vector<LevelData<EBFluxFAB>* >&       a_velocity,
                         Vector<LevelData<EBFluxFAB>* >&       a_gradient,
                         int                                   a_icomp)
{
  CH_TIME("EBCompositeMACProjector::correctVelocityComponent");
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      //this routine makes no sense otherwise
      CH_assert(a_gradient[ilev]->nComp() == 1);
      CH_assert(a_velocity[ilev]->nComp() == 1);

      a_gradient[ilev]->exchange(Interval(0,0));
      for (DataIterator dit = m_eblg[ilev].getDBL().dataIterator(); dit.ok(); ++dit)
        {
          //same comp (a_icomp) of velocity on every face.
          //gradient of phi is ordinary mac gradient
          //so it lives on a_icomp faces
          const EBFaceFAB& gradFAB = (*a_gradient[ilev])[dit()][a_icomp];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              EBFaceFAB& velFAB = (*a_velocity[ilev])[dit()][idir];
              if (idir == a_icomp)
                {
                  velFAB -= gradFAB;
                }
              else
                {
                  correctTangentialVelocity(velFAB, gradFAB,
                                            m_eblg[ilev].getDBL().get(dit()),
                                            m_eblg[ilev].getEBISL()[dit()],
                                            (*m_eblg[ilev].getCFIVS())[dit()]);
                }
            }
        }
      a_velocity[ilev]->exchange();
    }
}

void
EBCompositeMACProjector::
correctTangentialVelocity(EBFaceFAB&        a_velocity,
                          const EBFaceFAB&  a_gradient,
                          const Box&        a_grid,
                          const EBISBox&    a_ebisBox,
                          const IntVectSet& a_cfivs)
{
  CH_TIME("EBCompositeMACProjector::correctTangentialVelocity");
  int velDir = a_velocity.direction();
  int gradDir = a_gradient.direction();
  //veldir is the face on which the velocity lives
  //graddir is the direction of the component
  //and the face on which the mac gradient lives
  CH_assert(velDir != gradDir);
  CH_assert(a_velocity.nComp() == 1);
  CH_assert(a_gradient.nComp() == 1);

  //updating in place in fortran so we have to save a acopy
  EBFaceFAB velSave(a_ebisBox, a_velocity.getCellRegion(),  velDir, 1);
  velSave.copy(a_velocity);

  //interior box is the box where all of the stencil can be reached
  //without going out of the domain
  ProblemDomain domainBox = a_ebisBox.getDomain();
  Box interiorBox = a_grid;
  interiorBox.grow(1);
  interiorBox &= domainBox;
  interiorBox.grow(-1);
  Box interiorFaceBox = surroundingNodes(interiorBox, velDir);

  if (!interiorBox.isEmpty())
    {
      BaseFab<Real>&       regVel  = a_velocity.getSingleValuedFAB();
      const BaseFab<Real>& regGrad = a_gradient.getSingleValuedFAB();
      FORT_REGCORRECTTANVEL(CHF_FRA1(regVel,0),
                            CHF_CONST_FRA1(regGrad,0),
                            CHF_BOX(interiorFaceBox),
                            CHF_INT(velDir),
                            CHF_INT(gradDir));
    }

  //do only irregular and boundary cells pointwise
  IntVectSet ivsIrreg = a_ebisBox.getIrregIVS(a_grid);
  IntVectSet ivsGrid(a_grid);
  ivsGrid -= interiorBox;
  ivsGrid |= ivsIrreg;
  FaceIterator faceit(ivsGrid, a_ebisBox.getEBGraph(), velDir, FaceStop::SurroundingWithBoundary);
  for (faceit.reset(); faceit.ok(); ++faceit)
    {
      //average neighboring grads in grad dir direction
      int numGrads = 0;
      Real gradAve = 0.0;
      const FaceIndex& velFace = faceit();
      for (SideIterator sitVel; sitVel.ok(); ++sitVel)
        {
          const VolIndex& vofSide = velFace.getVoF(sitVel());
          const IntVect&   ivSide = vofSide.gridIndex();
          //do not include stuff over coarse-fine interface and inside the domain
          //cfivs includes cells just outside domain
          if (!a_cfivs.contains(ivSide) && domainBox.contains(ivSide))
            {
              for (SideIterator sitGrad; sitGrad.ok(); ++sitGrad)
                {
                  Vector<FaceIndex> gradFaces = a_ebisBox.getFaces(vofSide, gradDir, sitGrad());
                  for (int iface = 0; iface < gradFaces.size(); iface++)
                    {
                      if (!gradFaces[iface].isBoundary())
                        {
                          numGrads++;
                          gradAve += a_gradient(gradFaces[iface], 0);
                        }
                    }
                }//end loop over gradient sides
            }//end cfivs check
          else
            {
              // inside coarse/fine interface or at domain boundary. extrapolate from neighboring vofs to get the gradient
              const Side::LoHiSide inSide    = flip(sitVel());
              const VolIndex&      inSideVof = velFace.getVoF(inSide);
              const IntVect&       inSideIV  = inSideVof.gridIndex();

              IntVect inSideFarIV = inSideIV;
              inSideFarIV[velDir] += sign(inSide);
              if (domainBox.contains(inSideIV) && domainBox.contains(inSideFarIV))
                {
                  Vector<VolIndex> farVofs = a_ebisBox.getVoFs(inSideFarIV);
                  if (farVofs.size()  == 1)
                    {
                      const VolIndex& inSideFarVof = farVofs[0];
                      for (SideIterator sitGrad; sitGrad.ok(); ++sitGrad)
                        {
                          //get the grad for the face adjoining inSideVof on the sitGrad side, in the gradDir direction
                          Vector<FaceIndex> gradFaces    = a_ebisBox.getFaces(inSideVof   , gradDir, sitGrad());
                          Vector<FaceIndex> gradFarFaces = a_ebisBox.getFaces(inSideFarVof, gradDir, sitGrad());
                          if ( (gradFaces.size() == 1) && (gradFarFaces.size() == 1) )
                            {
                              // if ( (!gradFaces[0].isBoundary()) && (!gradFarFaces[0].isBoundary()) )
                              //   {
                                  const Real& inSideGrad    = a_gradient(gradFaces[0], 0);
                                  const Real& inSideFarGrad = a_gradient(gradFarFaces[0], 0);
                                  Real extrapGrad = 2.0*inSideGrad - inSideFarGrad;
                                  gradAve += extrapGrad;
                                  numGrads++;
                                // }
                            }
                        }
                    }
                }
            }//end cfivs check
        }//end loop over sides of velocity face
      if (numGrads > 1)
        {
          gradAve /= Real(numGrads);
        }

      //remember that the fortran updated the velocity  in place so
      //we have to use the saved velocity
      a_velocity(velFace, 0) = velSave(velFace, 0) - gradAve;
    }
}
/***/
void
EBCompositeMACProjector::
gradient(Vector<LevelData<EBFluxFAB>* >&       a_grad,
         const Vector<LevelData<EBCellFAB>* >& a_phi)
{
  CH_TIME("EBCompositeMACProjector::gradient");
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if (ilev>0)
        {
          m_quadCFI[ilev]->interpolate(*a_phi[ilev],
                                       *a_phi[ilev-1],
                                       Interval(0,0));
        }
      a_phi[ilev]->exchange();
      EBLevelMACProjector::setCurLevel(ilev);
      //compute the gradient ignoring other levels
      macGradient(*a_grad[ilev],  *a_phi[ilev],
                  m_eblg[ilev].getDBL(), m_eblg[ilev].getEBISL(),
                  m_eblg[ilev].getDomain(),   m_dx[ilev]);
    }
}
/***/
void
EBCompositeMACProjector::
enforceGradientBC(Vector<LevelData<EBFluxFAB>* >&       a_grad,
                  const Vector<LevelData<EBCellFAB>* >& a_phi)
{
  CH_TIME("EBCompositeMACProjector::enforceGradientBC");
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      EBLevelMACProjector::setCurLevel(ilev);
      macEnforceGradientBC(*a_grad[ilev],  *a_phi[ilev],
                           m_eblg[ilev].getDBL(), m_eblg[ilev].getEBISL(),
                           m_eblg[ilev].getDomain(),   m_dx[ilev],
                           m_origin,m_time,m_baseDomainBCFactPhi);
    }
}
/***/
void
EBCompositeMACProjector::
enforceVelocityBC(Vector<LevelData<EBFluxFAB>* >&  a_velocity,
                  const bool&                      a_doDivFreeOutflow)
{
  CH_TIME("EBCompositeMACProjector::enforceVelocityBC");
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      macEnforceVelocityBC(*a_velocity[ilev],
                           m_eblg[ilev].getDBL(),
                           m_eblg[ilev].getEBISL(),
                           m_eblg[ilev].getDomain(),
                           m_dx[ilev],
                           m_origin,
                           m_baseDomainBCFactVel,
                           m_time,
                           a_doDivFreeOutflow,
                           NULL);
    }
}

/***/
void
EBCompositeMACProjector::
enforceVelocityBC(Vector<LevelData<EBFluxFAB>* >&  a_velocity,
                  const int&                       a_comp,
                  const bool&                      a_doDivFreeOutflow)
{
  CH_TIME("EBCompositeMACProjector::enforceVelocityBC");
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      macEnforceVelocityBC(*a_velocity[ilev],
                           m_eblg[ilev].getDBL(),
                           m_eblg[ilev].getEBISL(),
                           m_eblg[ilev].getDomain(),
                           m_dx[ilev],
                           m_origin,
                           m_baseDomainBCFactVel,
                           m_time,
                           a_doDivFreeOutflow,
                           a_comp,
                           NULL);
    }
}

/*****/
void
EBCompositeMACProjector::
kappaDivergence(Vector<LevelData<EBCellFAB>* >&              a_divu,
                Vector<LevelData<EBFluxFAB>* >&              a_velo,
                const Vector<LevelData<BaseIVFAB<Real> >* >* a_boundaryVelo)
{
  CH_TIME("EBCompositeMACProjector::kappaDivergence");
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      EBLevelMACProjector::setCurLevel(ilev);

      EBFluxFactory fluxfact(m_eblg[ilev].getEBISL());
      //need one ghost cell so that exchange is meaningful
      LevelData<EBFluxFAB> centroidVelocity(m_eblg[ilev].getDBL(), 1, IntVect::Unit, fluxfact);

      Interval interv(0, 0);
      LevelData<EBFluxFAB>& velExch = (LevelData<EBFluxFAB>&)(*a_velo[ilev]);
      velExch.exchange(interv);

      //interpolate velocity to face centroids
      EBArith::interpolateFluxToCentroids(centroidVelocity, *a_velo[ilev],
                                          m_eblg[ilev].getDBL(), m_eblg[ilev].getEBISL(), m_eblg[ilev].getDomain());

      centroidVelocity.exchange(interv);
      LevelData<BaseIVFAB<Real> >*  levelBoundaryVel = NULL;
      if (a_boundaryVelo != NULL)
        {
          levelBoundaryVel = (*a_boundaryVelo)[ilev];
        }
      //compute the divergence igoring other levels
      macKappaDivergence(*a_divu[ilev], centroidVelocity,
                         m_eblg[ilev].getDBL(), m_eblg[ilev].getEBISL(),
                         m_eblg[ilev].getDomain(),   m_dx[ilev],
                         levelBoundaryVel);

      //use flux registers to correct divergence
      //at coarse-fine interface with velocity from finer
      //level
      if (ilev < (m_numLevels-1))
        {
          Real incrScale = 1.0;
          m_fluxReg[ilev]->setToZero();
          //increment with coarse values
          for (DataIterator dit = m_eblg[ilev].getDBL().dataIterator();dit.ok(); ++dit)
            {
              const EBFluxFAB& veloFlux = (*a_velo[ilev])[dit()];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  // This assumes that embedded boundaries and coarse-fine boundaries do not cross.
                  // To remove this assumption use incrementCoarseBoth.
                  m_fluxReg[ilev]->incrementCoarseRegular(veloFlux[idir], incrScale, dit(), interv, idir);
                }
            }
          //increment with fine velocities
          for (DataIterator dit = m_eblg[ilev+1].getDBL().dataIterator();dit.ok(); ++dit)
            {
              const EBFluxFAB& veloFlux = (*a_velo[ilev+1])[dit()];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      // This assumes that embedded boundaries and coarse-fine boundaries do not cross.
                      // To remove this assumption use incrementFineBoth.
                      m_fluxReg[ilev]->incrementFineRegular(veloFlux[idir], incrScale, dit(), interv, idir, sit());
                    }
                }
            }
          //reflux
          Real reflScale = 1.0/m_dx[ilev][0];
          m_fluxReg[ilev]->reflux(*a_divu[ilev], interv, reflScale);
        }
    }
}
/*****/
EBCompositeMACProjector::
~EBCompositeMACProjector()
{
  delete m_opFactory;
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      delete    m_divu[ilev];
      delete     m_phi[ilev];
      delete m_fluxReg[ilev];
    }
}

#include "NamespaceFooter.H"
