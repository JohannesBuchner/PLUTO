#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "parstream.H"
#include "ParmParse.H"
#include "PolyGeom.H"

#include "DebugOut.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "EBCellFAB.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "NeumannPoissonDomainBC.H"
#include "NeumannPoissonEBBC.H"
#include "SPMD.H"
#include "EBLoadBalance.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "AMRMultiGrid.H"
#include "EBAMRIO.H"
#include "BaseIVFactory.H"
#include "EBViscousTensorOpFactory.H"
#include "EBConductivityOpFactory.H"
#include "EBAMRPoissonOpFactory.H"
#include "KappaSquareNormal.H"

#include "AMRLevel.H"
#include "EBAMRCNS.H"
#include "EBCellFactory.H"
#include "BaseIVFactory.H"
#include "VoFIterator.H"
#include "EBPatchGodunovF_F.H"
#include "EBPlanarShockF_F.H"
#include "EBPatchPolytropicF_F.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "EBAMRDataOps.H"
#include "EBLGIntegrator.H"
#include "EBLevelDataOps.H"
#include "EBBackwardEuler.H"
#include "AMR.H"
#include "GodunovGeom.H"
#include "EBNoFlowIBC.H"
#include "EBExplosionIBCFactory.H"
#include "EBPlanarShockIBCFactory.H"
#include "EBEllipticLoadBalance.H"

#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5
int                                                      EBAMRCNS::s_NewPlotFile  = 0;
#endif
RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > EBAMRCNS::s_tempFactory        =  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > > EBAMRCNS::s_veloFactory        =  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
RefCountedPtr<MomentumBackwardEuler>                     EBAMRCNS::s_veloIntegratorBE   =  RefCountedPtr<MomentumBackwardEuler>();
RefCountedPtr< EBLevelBackwardEuler>                     EBAMRCNS::s_tempIntegratorBE   =  RefCountedPtr< EBLevelBackwardEuler>();

RefCountedPtr<MomentumTGA>                     EBAMRCNS::s_veloIntegratorTGA   =  RefCountedPtr<MomentumTGA>();
RefCountedPtr< EBLevelTGA>                     EBAMRCNS::s_tempIntegratorTGA   =  RefCountedPtr< EBLevelTGA>();
RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > > EBAMRCNS::s_tempSolver         =  RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >();
RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > > EBAMRCNS::s_veloSolver         =  RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >();
//BiCGStabSolver<LevelData<EBCellFAB> >                    EBAMRCNS::s_botSolver          = BiCGStabSolver<LevelData<EBCellFAB> >();
EBSimpleSolver                 EBAMRCNS::s_botSolver;
bool EBAMRCNS::s_solversDefined = false;
bool EBAMRCNS::s_noEBCF = false;
bool dumpStuff = false;


bool 
EBAMRCNS::
convergedToSteadyState()
{
  EBCellFactory      cellFact(m_eblg.getEBISL());
  LevelData<EBCellFAB>  udiff(m_eblg.getDBL(), m_stateNew.nComp(), 4*IntVect::Unit, cellFact);
  EBLevelDataOps::setToZero(udiff);
  EBLevelDataOps::incr(udiff, m_stateOld, -1.0);
  EBLevelDataOps::incr(udiff, m_stateNew,  1.0);
  Real umax, umin, eps;
  Real dmax, dmin; 
  int ivar;
  ParmParse pp;
  pp.get("convergence_metric", eps);
  pp.get("convergence_variable", ivar);
  EBLevelDataOps::getMaxMin(dmax, dmin, udiff, ivar);
  EBLevelDataOps::getMaxMin(umax, umin, m_stateNew, ivar);
  Real denom = 1;
  if(Abs(umax - umin) > eps)
    {
      denom = Abs(umax-umin);
    }
  Real maxdiff = Abs(dmax-dmin)/denom;
  pout() << "max difference in convergence variable = " << maxdiff << ", eps set to " << eps << endl;
  return (maxdiff < eps);
}
//---------------------------------------------------------------------------------------
EBAMRCNS::
EBAMRCNS(const EBAMRCNSParams& a_params,
         const RefCountedPtr<EBPatchPolytropicFactory>& a_godFactory,
         const RefCountedPtr<EBSpaceTimeFunction>& a_initialConditions):
  m_params(a_params),
  m_ebPatchGodunovFactory(a_godFactory),
  m_ICs(a_initialConditions),
  m_originalMass(0.0),
  m_originalMomentum(RealVect::Zero),
  m_originalEnergy(0.0),
  m_isConservative(false)
{
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
EBAMRCNS::
~EBAMRCNS()
{
  if((m_level == 0) && m_params.m_doDiffusion)
    {
      clearSolvers();
    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
clearSolvers()
{
  s_tempFactory  =  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
  s_veloFactory  =  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >();
  s_veloIntegratorBE  =  RefCountedPtr<MomentumBackwardEuler>();
  s_tempIntegratorBE  =  RefCountedPtr< EBLevelBackwardEuler>();
  s_veloIntegratorTGA  =  RefCountedPtr<MomentumTGA>();
  s_tempIntegratorTGA  =  RefCountedPtr< EBLevelTGA>();
  s_tempSolver   =  RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >();
  s_veloSolver   =  RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >();
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
setAirDiffusionCoefficients(const LevelData<EBCellFAB>& a_densCell,
                            const LevelData<EBCellFAB>& a_tempCell,
                            const LevelData<EBFluxFAB>& a_densFace,
                            const LevelData<EBFluxFAB>& a_tempFace)
{
  CH_TIME("set_air_diffusion_coef");
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      Box cellBox = m_eblg.getDBL().get(dit());
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          Box faceBox = surroundingNodes(cellBox, idir);
          FORT_AIRVISCOSITY(CHF_FRA1((*m_eta  )      [dit()][idir].getSingleValuedFAB(), 0),
                            CHF_CONST_FRA1(a_densFace[dit()][idir].getSingleValuedFAB(), 0),
                            CHF_CONST_FRA1(a_tempFace[dit()][idir].getSingleValuedFAB(), 0),
                            CHF_BOX(faceBox));

          FORT_AIRTHERMDIFF(CHF_FRA1((*m_bcoTemp)    [dit()][idir].getSingleValuedFAB(), 0),
                            CHF_CONST_FRA1(a_densFace[dit()][idir].getSingleValuedFAB(), 0),
                            CHF_CONST_FRA1(a_tempFace[dit()][idir].getSingleValuedFAB(), 0),
                            CHF_BOX(faceBox));
          
          Vector<FaceIndex> faces = m_eblg.getEBISL()[dit()].getEBGraph().getMultiValuedFaces(idir, cellBox);
          for(int iface = 0; iface < faces.size(); iface++)
            {
              FORT_POINTAIRVISCOSITY(CHF_REAL((*m_eta  )      [dit()][idir](faces[iface], 0)),
                                     CHF_CONST_REAL(a_densFace[dit()][idir](faces[iface], 0)),
                                     CHF_CONST_REAL(a_tempFace[dit()][idir](faces[iface], 0)));


              FORT_POINTAIRTHERMDIFF(CHF_REAL((*m_bcoTemp)    [dit()][idir](faces[iface], 0)),
                                     CHF_CONST_REAL(a_densFace[dit()][idir](faces[iface], 0)),
                                     CHF_CONST_REAL(a_tempFace[dit()][idir](faces[iface], 0)));
            }
          (*m_lambda)[dit()][idir].copy((*m_eta)[dit()][idir]);
          (*m_lambda)[dit()][idir] *= (-2.0/3.0);

        }
      IntVectSet ivs = m_eblg.getEBISL()[dit()].getIrregIVS(cellBox); //irreg because we are also setting the irregular coeffs here
      for(VoFIterator vofit(ivs, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          FORT_POINTAIRVISCOSITY(CHF_REAL((*m_etaIrreg)    [dit()](vofit(), 0)),
                                 CHF_CONST_REAL(a_densCell [dit()](vofit(), 0)),
                                 CHF_CONST_REAL(a_tempCell [dit()](vofit(), 0)));

                                                                   
          FORT_POINTAIRTHERMDIFF(CHF_REAL((*m_bcoTempIrreg)[dit()](vofit(), 0)),
                                 CHF_CONST_REAL(a_densCell [dit()](vofit(), 0)),
                                 CHF_CONST_REAL(a_tempCell [dit()](vofit(), 0)));

          (*m_lambdaIrreg)[dit()](vofit(), 0) = (-2.0/3.0)*(*m_etaIrreg)[dit()](vofit(), 0);
        }

    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
fillCoefficients(const LevelData<EBCellFAB>& a_state)
{
  CH_TIME("EBAMRCNS::fillcoef");
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      // All of these coefficients are time-independent and can be set here.
      (*m_eta     )    [dit()].setVal(m_params.m_viscosityMu);
      (*m_bcoTemp )    [dit()].setVal(m_params.m_thermalCond);
      (*m_etaIrreg)    [dit()].setVal(m_params.m_viscosityMu);
      (*m_lambda     ) [dit()].setVal(m_params.m_viscosityLa);
      (*m_lambdaIrreg) [dit()].setVal(m_params.m_viscosityLa);
      (*m_bcoTempIrreg)[dit()].setVal(m_params.m_thermalCond);
    }
  if(m_level == 0)
    {
      if((!m_params.m_useAirCoefs))
        {
          pout() << "using constant diffusion coefficients"  ;
          pout() << ": mu = " << m_params.m_viscosityMu; 
          pout() << ", lambda = " << m_params.m_viscosityLa ;
          pout() << ", kappa = " << m_params.m_thermalCond;
          pout() << ", cp = " << m_params.m_specHeatCv  << endl;
        }
    }
  if(m_params.m_useAirCoefs)
    {
      pout() << "specific heat = " << m_params.m_specHeatCv  << endl;
      pout() << "using diffusion coefficients for air"  << endl;
      EBCellFactory cellFact(m_eblg.getEBISL());
      EBFluxFactory fluxFact(m_eblg.getEBISL());
      LevelData<EBCellFAB>  tempCell(m_eblg.getDBL(), 1, 4*IntVect::Unit, cellFact);
      LevelData<EBCellFAB>  densCell(m_eblg.getDBL(), 1, 4*IntVect::Unit, cellFact);
      LevelData<EBFluxFAB>  tempFace(m_eblg.getDBL(), 1,   IntVect::Zero, fluxFact);
      LevelData<EBFluxFAB>  densFace(m_eblg.getDBL(), 1,   IntVect::Zero, fluxFact);

      getTemperature(tempCell, a_state);
      Interval srcComp(CRHO, CRHO);
      Interval dstComp(0, 0);
      a_state.copyTo(srcComp, densCell, dstComp);

      if(m_hasCoarser)
        {
          EBAMRCNS* coarCNS = getCoarserLevel();
          EBPWLFillPatch patcher(m_eblg.getDBL()  , coarCNS->m_eblg.getDBL(),
                                 m_eblg.getEBISL(), coarCNS->m_eblg.getEBISL(),
                                 coarCNS->m_eblg.getDomain(), m_ref_ratio, 1, 4);


          EBCellFactory coarFact(coarCNS->m_eblg.getEBISL());
          LevelData<EBCellFAB>  coarTemp(coarCNS->m_eblg.getDBL(), 1, 4*IntVect::Unit, coarFact);
          LevelData<EBCellFAB>  coarDens(coarCNS->m_eblg.getDBL(), 1, 4*IntVect::Unit, coarFact);

          coarCNS->getTemperature(coarTemp, coarCNS->m_stateNew);
          coarCNS->m_stateNew.copyTo(srcComp, coarDens, dstComp);

          patcher.interpolate(tempCell, coarTemp, coarTemp, m_time, m_time, m_time, Interval(0,0));
          patcher.interpolate(densCell, coarDens, coarDens, m_time, m_time, m_time, Interval(0,0));

        }
      EBLevelDataOps::averageCellToFaces(tempFace, tempCell, m_eblg.getDBL(), m_eblg.getEBISL(), m_eblg.getDomain(), 0);
      EBLevelDataOps::averageCellToFaces(densFace, densCell, m_eblg.getDBL(), m_eblg.getEBISL(), m_eblg.getDomain(), 0);

      setAirDiffusionCoefficients(densCell, tempCell, tempFace, densFace);
    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
defineFactories(bool a_atHalfTime)
{
  CH_TIME("EBAMRCNS::defineFactories");
  if(m_params.m_doDiffusion)
    {
      Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
      int nlevels = hierarchy.size();
      Vector<int>                                           refRat(nlevels);
      Vector<EBLevelGrid>                                   eblgs(nlevels);
      Vector<DisjointBoxLayout>                             grids(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> >        >  eta(nlevels);
      Vector<RefCountedPtr<LevelData<EBCellFAB> >        >  acoVelo(nlevels);
      Vector<RefCountedPtr<LevelData<EBCellFAB> >        >  acoTemp(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> >        >  bcoTemp(nlevels);
      Vector<RefCountedPtr<LevelData<EBFluxFAB> >        >  lambda(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  etaIrreg(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  bcoTempIrreg(nlevels);
      Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >  lambdaIrreg(nlevels);
      Vector<RefCountedPtr<EBQuadCFInterp> >                quadCFI(nlevels);

      EBAMRCNS* coarsestLevel = (EBAMRCNS*)(hierarchy[0]);
      Real           lev0Dx      = (coarsestLevel->m_dx[0]);
      ProblemDomain lev0Dom      =  coarsestLevel->m_eblg.getDomain();


      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          EBAMRCNS* cnsLevel = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);

          EBCellFactory fact(cnsLevel->m_eblg.getEBISL());
          LevelData<EBCellFAB> halfSt(cnsLevel->m_eblg.getDBL(),m_nComp, 4*IntVect::Unit, fact);
          if(a_atHalfTime)
            {
              cnsLevel->getHalfState(halfSt);
            }
          else
            {
              Interval interv(0, m_nComp-1);
              cnsLevel->m_stateOld.copyTo(interv, halfSt, interv);
            }
          cnsLevel->fillCoefficients(halfSt);

          eblgs       [ilev] = cnsLevel->m_eblg;
          grids       [ilev] = cnsLevel->m_eblg.getDBL();
          refRat      [ilev] = cnsLevel->m_ref_ratio;
          acoVelo     [ilev] = cnsLevel->m_acoVelo;
          acoTemp     [ilev] = cnsLevel->m_acoTemp;
          eta         [ilev] = cnsLevel->m_eta;
          etaIrreg    [ilev] = cnsLevel->m_etaIrreg;
          lambda      [ilev] = cnsLevel->m_lambda;
          lambdaIrreg [ilev] = cnsLevel->m_lambdaIrreg;
          quadCFI     [ilev] = cnsLevel->m_quadCFI;
          bcoTemp     [ilev] = cnsLevel->m_bcoTemp;
          bcoTempIrreg[ilev] = cnsLevel->m_bcoTempIrreg;
        }

      //alpha, beta get replaced in tga solves
      Real alpha = 1;
      Real beta  = 1;
      IntVect  giv    = 4*IntVect::Unit;

      // Viscous tensor operator.
      bool noMG = true;
      s_veloFactory =
        RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
        (dynamic_cast<AMRLevelOpFactory<LevelData<EBCellFAB> >*>
         (new EBViscousTensorOpFactory(eblgs, alpha, beta, acoVelo, eta,
                                       lambda, etaIrreg, lambdaIrreg,lev0Dx, refRat,
                                       m_params.m_doBCVelo, m_params.m_ebBCVelo,giv, giv, -1, noMG)));


      EBViscousTensorOp::doLazyRelax(m_params.m_doLazyRelax);

      // Thermal diffusion operator.
      int relaxType = 1;
      s_tempFactory =
        RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
        (dynamic_cast<AMRLevelOpFactory<LevelData<EBCellFAB> >*>
         (new EBConductivityOpFactory(eblgs, quadCFI, alpha, beta, acoTemp, 
                                      bcoTemp, bcoTempIrreg, lev0Dx, 
                                      refRat, m_params.m_doBCTemp, 
                                      m_params.m_ebBCTemp, giv, giv, relaxType)));

    }
  (*m_eta     )    .exchange(Interval(0,0));
  (*m_bcoTemp )    .exchange(Interval(0,0));
  (*m_lambda     ) .exchange(Interval(0,0));
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
defineSolvers()
{
  CH_TIME("EBAMRCNS::defineSolvers");
  ParmParse pp;
  bool tagAllIrregular = false;
  if(pp.contains("tag_all_irregular"))
    {
      pp.get("tag_all_irregular", tagAllIrregular);
    }
  if(tagAllIrregular) s_noEBCF = true;

  EBConductivityOp::setForceNoEBCF( s_noEBCF);
  EBViscousTensorOp::setForceNoEBCF(s_noEBCF);
  defineFactories(true);
  if(m_params.m_doDiffusion)
    {
      Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
      int nlevels = hierarchy.size();
      Vector<int>                                           refRat(nlevels);
      Vector<EBLevelGrid>                                   eblgs(nlevels);
      Vector<DisjointBoxLayout>                             grids(nlevels);

      EBAMRCNS* coarsestLevel = dynamic_cast<EBAMRCNS*>(hierarchy[0]);
      ProblemDomain lev0Dom      =  coarsestLevel->m_eblg.getDomain();
      for(int ilev = 0; ilev < nlevels; ilev++)
        {
          EBAMRCNS* cnsLevel = (EBAMRCNS*)(hierarchy[ilev]);
          eblgs       [ilev] = cnsLevel->m_eblg;
          grids       [ilev] = cnsLevel->m_eblg.getDBL();
          refRat      [ilev] = cnsLevel->m_ref_ratio;
        }

      //alpha, beta get replaced in tga solves
      IntVect  giv    = 4*IntVect::Unit;
      RealVect origin =  RealVect::Zero;
      int numSmooth, numMG, maxIter, mgverb;
      Real tolerance, hang, normThresh;
      ParmParse pp("amrmultigrid");
      pp.get("num_smooth", numSmooth);
      pp.get("num_mg",     numMG);
      pp.get("hang_eps",   hang);
      pp.get("norm_thresh",normThresh);
      pp.get("tolerance",  tolerance);
      pp.get("max_iter",   maxIter);
      pp.get("verbosity",  mgverb);

      s_veloSolver = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >( new AMRMultiGrid< LevelData<EBCellFAB> > ());
      s_tempSolver = RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >( new AMRMultiGrid< LevelData<EBCellFAB> > ());

      s_veloSolver->define(lev0Dom, *s_veloFactory, &s_botSolver, nlevels);
      s_tempSolver->define(lev0Dom, *s_tempFactory, &s_botSolver, nlevels);

      s_veloSolver->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);
      s_tempSolver->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG, maxIter, tolerance, hang, normThresh);

      s_veloSolver->m_verbosity = mgverb;
      s_tempSolver->m_verbosity = mgverb;
      Real bottomCushion = 1.0;
      if(pp.contains("bottom_cushion"))
        {
          pp.get("bottom_cushion", bottomCushion);
        }
      s_veloSolver->m_bottomSolverEpsCushion = bottomCushion;
      s_tempSolver->m_bottomSolverEpsCushion = bottomCushion;
      //      s_botSolver.m_verbosity   = 0;
      //  s_botSolver.m_numRestarts = 0;

      s_tempIntegratorBE  = RefCountedPtr< EBLevelBackwardEuler>( new  EBLevelBackwardEuler(grids, refRat, lev0Dom, s_tempFactory, s_tempSolver));
      s_veloIntegratorBE  = RefCountedPtr<MomentumBackwardEuler>(new  MomentumBackwardEuler(grids, refRat, lev0Dom, s_veloFactory, s_veloSolver));
      s_veloIntegratorBE->setEBLG(eblgs);
      s_tempIntegratorBE->setEBLG(eblgs);

      s_tempIntegratorTGA  = RefCountedPtr< EBLevelTGA>( new  EBLevelTGA(grids, refRat, lev0Dom, s_tempFactory, s_tempSolver));
      s_veloIntegratorTGA  = RefCountedPtr<MomentumTGA>(new  MomentumTGA(grids, refRat, lev0Dom, s_veloFactory, s_veloSolver));
      s_veloIntegratorTGA->setEBLG(eblgs);
      s_tempIntegratorTGA->setEBLG(eblgs);
    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
define(AMRLevel*            a_coarser_level_ptr,
       const ProblemDomain& a_problem_domain,
       int                  a_level,
       int                  a_ref_ratio)
{
  CH_TIME("EBAMRCNS::define");
  if(m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::define, level=" << a_level << endl;
    }

  EBPatchGodunov::useConservativeSource(true);
  m_isDefined = true;
  AMRLevel::define(a_coarser_level_ptr,
                   a_problem_domain,
                   a_level,
                   a_ref_ratio);

  if (a_coarser_level_ptr != NULL)
    {
      EBAMRCNS* amrg_ptr =
        dynamic_cast<EBAMRCNS*>(a_coarser_level_ptr);
      if (amrg_ptr == NULL)
        {
          pout() << "EBAMRG::define:cannot cast  to derived class"
                 << endl;
          MayDay::Error();
        }
      m_params = amrg_ptr->m_params;
    }
  Real dxSize  = m_params.m_domainLength/m_problem_domain.domainBox().size(0);
  m_dx = dxSize*RealVect::Unit;
  m_nGhost = 4;

  m_ebPatchGodunov = RefCountedPtr<EBPatchGodunov>();
  m_ebPatchGodunov = RefCountedPtr<EBPatchGodunov>(m_ebPatchGodunovFactory->create());
  m_ebPatchGodunov->define(m_problem_domain, m_dx);

  m_nComp      = m_ebPatchGodunov->numConserved();
  m_nPrim      = m_ebPatchGodunov->numPrimitives();
  m_stateNames = m_ebPatchGodunov->stateNames();
  m_primNames  = m_ebPatchGodunov->primNames();
  m_ref_ratio  = a_ref_ratio;
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
Real
EBAMRCNS::
advance()
{
  CH_TIME("EBAMRCNS::advance");
  EBPatchGodunov::s_whichLev = m_level;
  EBViscousTensorOp::s_step = AMR::s_step;
  //this is so I can output a rhs
  EBViscousTensorOp::s_whichLev = m_level;
  m_dtOld = m_dt;

  if(m_params.m_variableCoeff && (m_level== 0) && m_params.m_doDiffusion) defineSolvers();

  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS advance for level " << m_level << endl;
    }
  m_stateNew.copyTo(m_stateNew.interval(), m_stateOld, m_stateOld.interval());

  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB>      divergeF(m_eblg.getDBL(),m_nComp          , 4*IntVect::Unit, fact);
  
  pout() << "getting diverence of flux" << endl;

  //PDivU and rhouedele are not multiplied by the volume fraction (that happens in the integrator.
  fluxDivergence( divergeF);

  pout() << "advancing density" << endl;
  if(!m_params.m_doDiffusion)
    {
      explicitAdvance(divergeF);
    }
  else
    {
      CH_TIME("diffusion_dance");
      // step 1--get u*
      pout() << "step 1: getting Ustar (U^n + explicit contribution)" << endl;

      //U* = Un - dt*L^c(U^n)
      LevelData<EBCellFAB>       UStar(m_eblg.getDBL(), m_nComp, 4*IntVect::Unit, fact);
      getUStar(UStar, m_stateOld, divergeF);

      pout() << "step 1.5: doing explicit redistribution into ustar of hyperbolic mass difference" <<endl;
      //redistribute stuff in hyperbolic registers.
      hyperbolicRedistribution(UStar);

      ////debug
      //LevelData<EBCellFAB> consAndPrimStar;
      //fillConsAndPrim(consAndPrimStar, UStar, 1);
      ////end debug

      pout() << "step 2: solving for div(sigma)" <<endl;

      //(\rho^* I - \dt L^m) \ubold^\npo = \dt (\rho \ubold)^*
      //L^m(\ubold^\nph) = \rho^*\left(\frac{\ubold^\npo - \ubold^*}{\dt}\right)
      //(\rho \ubold)^\npo =(\rho \ubold)^* + \dt L^m(U^*)
      LevelData<EBCellFAB>    divSigma(m_eblg.getDBL(),SpaceDim, 4*IntVect::Unit, fact);
      getDivSigma(divSigma, UStar);

      //debug
      //LevelData<EBCellFAB> consAndPrimDS;
      //fillConsAndPrim(consAndPrimDS, UStar, 1);
      //LevelData<EBCellFAB> consAndPrimSN;
      //fillConsAndPrim(consAndPrimSN, m_stateNew, 1);
      //end debug

      pout() << "step 3: getting the components of div Ld = div(sigma u)" << endl;
      //$L^d(U^*) =u \cdot \grad \cdot \sigma) + \sigma \cdot (\grad \ubold)$
      //(\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
      //LevelData<EBCellFAB>  uDotDivSig(m_eblg.getDBL(),       1, 4*IntVect::Unit, fact);
      //LevelData<EBCellFAB>    dissFunc(m_eblg.getDBL(),       1, 4*IntVect::Unit, fact);
      //getSplitLdOfU(dissFunc, uDotDivSig,  UStar, divSigma);
      LevelData<EBCellFAB>   divSigmaU(m_eblg.getDBL(),       1, 4*IntVect::Unit, fact);
      getSingleLdOfU(divSigmaU,  UStar);

      //debug
      //      EBLevelDataOps::setVal(divSigmaU, 0.0);
      //end debug

      //debug
      //LevelData<EBCellFAB> consAndPrimDbst;
      //fillConsAndPrim(consAndPrimDbst, UStar, 1);
      //end debug

      /** (\kappa \rho^\npo C_v I - \frac{\kappa \dt} L^k)T^\npo = \kappa \rho^\npo C_v T^{**}
          (\rho E)^{n+1} = (\rho E)^{**} + \dt L^k(T**)      **/
      pout() << "step 4: solve for conduction term and add to energy" << endl;
      LevelData<EBCellFAB> dtDivKGradT(m_eblg.getDBL(),       1, 4*IntVect::Unit, fact);
      getDivKappaGradT(dtDivKGradT, UStar);

      pout() << "putting state into m_statenew and flooring" << endl;
      finalAdvance(UStar);
    }
  Real new_dt;
  {
    CH_TIME("time_step_calc");
    //deal with time and time step
    Real maxWaveSpeed = m_ebLevelGodunov.getMaxWaveSpeed(m_stateOld);
    new_dt = m_params.m_cfl*m_dx[0]/maxWaveSpeed;
    pout() << "max wave speed = " << maxWaveSpeed << ", dx = " << m_dx[0] <<", new_dt = "  << new_dt << endl;
    m_time += m_dt;
  }

  //save stable timestep to save computational effort
  m_dtNew = new_dt;
  if(m_params.m_checkMaxMin)
    {
      CH_TIME("max_min_check");
      pout() << "after advance max mins for data for level " << m_level << endl;
      for(int icomp = 0; icomp < m_stateNew.nComp(); icomp++)
        {
          Real maxVal, minVal;
          EBLevelDataOps::getMaxMin(maxVal, minVal, m_stateNew, icomp);
          pout() << " "
                 << setprecision(4)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 <<  "comp = " <<  icomp 
                 << ", max = " << maxVal 
                 << ", min = " << minVal << endl;
        }
    }
  return new_dt;
}
/////////
void  
EBAMRCNS::
floorVariable(LevelData<EBCellFAB>& a_data,
              const EBLevelGrid a_eblg,
              Real a_minVal, int a_ivar)

{
  CH_TIME("EBAMRCNS::floorVariable");
  const DisjointBoxLayout& dbl = a_eblg.getDBL();
  for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = dbl[dit()];
      FORT_FLOORSCALAR(CHF_FRA1(a_data[dit()].getSingleValuedFAB(),a_ivar),
                       CHF_CONST_REAL(a_minVal),
                       CHF_BOX(region));
      //only need multivalued cells because this is pointwise
      IntVectSet ivs =a_eblg.getEBISL()[dit()].getMultiCells(region);
      for(VoFIterator vofit(ivs, a_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          Real val = Max(a_data[dit()](vofit(), a_ivar), a_minVal);
          a_data[dit()](vofit(), a_ivar) = val;
        }
    }
}
///
//$L^d(U^*) =u \cdot \grad \cdot \sigma) + \sigma \cdot (\grad \ubold)$
//(\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
void  
EBAMRCNS::
getSingleLdOfU(LevelData<EBCellFAB>      & a_divSigmaU,
               LevelData<EBCellFAB>      & a_UStar)
{

  CH_TIME("EBAMRCNS::getSingleLdOfU");
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> velStar(m_eblg.getDBL(), SpaceDim, m_nGhost*IntVect::Unit, fact);
  LevelData<EBCellFAB> velHalf(m_eblg.getDBL(), SpaceDim, m_nGhost*IntVect::Unit, fact);

  getVelocity(velStar, a_UStar);
  LevelData<EBCellFAB>  kappaConsDivSigmaU(m_eblg.getDBL(), 1,4*IntVect::Unit, fact);
  LevelData<EBCellFAB>    nonConsDivSigmaU(m_eblg.getDBL(), 1,4*IntVect::Unit, fact);
  LevelData<EBCellFAB>* veloCoar = NULL;
  LevelData<EBCellFAB>* tempCoar = NULL;
  getCoarserTemperatureAndVelocity(tempCoar, veloCoar);
  s_veloIntegratorBE->getKappaDivSigmaU(kappaConsDivSigmaU,
                                        velStar,
                                        veloCoar,
                                        m_level);
 
  //now get the non conservative version
  EBLevelDataOps::setToZero(nonConsDivSigmaU);
  EBLevelDataOps::incr     (nonConsDivSigmaU, kappaConsDivSigmaU, 1.0);
  KappaSquareNormal normalizinator(m_eblg);
  normalizinator(nonConsDivSigmaU);

  //(\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
  //     mass diff = kappa(1-kappa)*(kappaConsDissFcn - nonConsDissFcn)
  updateEnergyBySingleLdAndRedistribute(a_divSigmaU, 
                                        kappaConsDivSigmaU,
                                        nonConsDivSigmaU,
                                        a_UStar);

  if(m_hasCoarser)
    {
      delete veloCoar;
      delete tempCoar;
    }
  
}


///
//$L^d(U^*) =u \cdot \grad \cdot \sigma) + \sigma \cdot (\grad \ubold)$
//(\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
void  
EBAMRCNS::
getSplitLdOfU(LevelData<EBCellFAB>      & a_dissFunc, 
              LevelData<EBCellFAB>      & a_uDotDivSig,  
              LevelData<EBCellFAB>      & a_UStar, 
              const LevelData<EBCellFAB>& a_divSigma)
{
  CH_TIME("EBAMRCNS::getSplitLdOfU");
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> velStar(m_eblg.getDBL(), SpaceDim, m_nGhost*IntVect::Unit, fact);
  getVelocity(velStar, a_UStar);

  //get udot divsigma = velstar dot divsigma
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_uDotDivSig[dit()].setVal(0.);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          EBCellFAB incr(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), 1);
          incr.setVal(0.);
          int isrc = idir; int idst = 0; int inco = 1;
          //this makes incr = divsigma_d
          incr.plus(a_divSigma[dit()], isrc, idst, inco);
          //this makes incr = divsigma_d*vel_d
          incr.mult(   velStar[dit()], isrc, idst, inco);
          
          a_uDotDivSig[dit()] += incr;
        }
    }
  LevelData<EBCellFAB>        kappaUDotDivSigma(m_eblg.getDBL(), 1,4*IntVect::Unit, fact);
  EBLevelDataOps::setToZero  (kappaUDotDivSigma);
  EBLevelDataOps::incr       (kappaUDotDivSigma, a_uDotDivSig, 1.0);
  EBLevelDataOps::kappaWeight(kappaUDotDivSigma);
  
  LevelData<EBCellFAB>  kappaDivSigmaU(m_eblg.getDBL(), 1,4*IntVect::Unit, fact);
  LevelData<EBCellFAB>* veloCoar = NULL;
  LevelData<EBCellFAB>* tempCoar = NULL;
  getCoarserTemperatureAndVelocity(tempCoar, veloCoar);
  s_veloIntegratorBE->getKappaDivSigmaU(kappaDivSigmaU,
                                        velStar,
                                        veloCoar,
                                        m_level);
 
  ///kappa diss = kappa div(sigma u) - kappa (udot div sigma)
  LevelData<EBCellFAB>  kappaConsDissFunc(m_eblg.getDBL(), 1,4*IntVect::Unit, fact);
  EBLevelDataOps::setToZero(kappaConsDissFunc);
  EBLevelDataOps::incr     (kappaConsDissFunc, kappaDivSigmaU,     1.0);
  EBLevelDataOps::incr     (kappaConsDissFunc, kappaUDotDivSigma, -1.0);

  //dissipation should only add energy to the flow
  //  floorVariable(kappaConsDissFunc, m_eblg, 0.0, 0);

  //now get the non conservative version
  LevelData<EBCellFAB>   nonConsDissFunc(m_eblg.getDBL(), 1,4*IntVect::Unit, fact);
  EBLevelDataOps::setToZero(nonConsDissFunc);
  EBLevelDataOps::incr     (nonConsDissFunc, kappaConsDissFunc, 1.0);
  KappaSquareNormal normalizinator(m_eblg);
  normalizinator(nonConsDissFunc);  


  //(\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
  //     mass diff = kappa(1-kappa)*(kappaConsDissFcn - nonConsDissFcn)
  updateEnergyBySplitLdAndRedistribute(a_dissFunc,   
                                       a_uDotDivSig, 
                                       kappaConsDissFunc,
                                       nonConsDissFunc,
                                       a_UStar);

  if(m_hasCoarser)
    {
      delete veloCoar;
      delete tempCoar;
    }
  
}
//(\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
//     mass diff = kappa(1-kappa)*(kappaConsDissFcn - nonConsDissFcn)
void  
EBAMRCNS::
updateEnergyBySplitLdAndRedistribute(LevelData<EBCellFAB>      & a_dissFunc, 
                                     LevelData<EBCellFAB>      & a_uDotDivSigma,  
                                     LevelData<EBCellFAB>      & a_kappaConsDissFunc,
                                     LevelData<EBCellFAB>      & a_nonConsDissFunc,
                                     LevelData<EBCellFAB>      & a_UStar)
{
  CH_TIME("EBAMRCNS::updateEnergyBySplitLdAndRedistribute");
  //make the dissipation function = kappa*cons + (1-kappa)*noncons
  EBLevelDataOps::setToZero(a_dissFunc);
  EBLevelDataOps::incr(a_dissFunc, a_kappaConsDissFunc, 1.0);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      IntVectSet ivs = m_eblg.getEBISL()[dit()].getIrregIVS(region);
      for(VoFIterator vofit(ivs, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          Real kappa   =  m_eblg.getEBISL()[dit()].volFrac(vofit());
          Real kapCons =  a_kappaConsDissFunc[dit()](vofit(), 0);
          Real nonCons =    a_nonConsDissFunc[dit()](vofit(), 0);
          a_dissFunc[dit()](vofit(), 0) = kapCons + (1.-kappa)*nonCons;
          m_massDiff[dit()](vofit(), CENG) = (m_dt*(1.-kappa)*kapCons - kappa*(1.-kappa)*nonCons);
        }
    }
  Interval enerint(CENG, CENG);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB incr(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), 1);
      incr.setVal(0.);
      //this makes incr = dissipation function = gradu dot sigma
      incr += a_dissFunc[dit()];
      //this makes incr = dissipation function + udot div sigma
      incr += a_uDotDivSigma[dit()];
      //this makes incr = dt*(dissipation function + udot div sigma)
      incr *= m_dt;
      //this adds effects of dissipation into the energy
      int isrc = 0; int idst = CENG; int inco = 1;
      a_UStar[dit()].plus(incr, isrc, idst, inco);

      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), enerint);
      m_massDiff[dit()].setVal(0.);
    }

  //this smushes the energy difference into the state
  
  a_UStar.exchange(Interval(0, m_nComp-1));
  //  m_ebLevelRedist.resetWeights(a_UStar, CENG);
  m_ebLevelRedist.resetWeights(a_UStar, CENG);

  m_ebLevelRedist.redistribute(m_redisRHS, enerint);
  m_ebLevelRedist.setToZero();
}
//(\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
//     mass diff = kappa(1-kappa)*(kappaConsDissFcn - nonConsDissFcn)
void  
EBAMRCNS::
updateEnergyBySingleLdAndRedistribute(LevelData<EBCellFAB>      & a_divSigmaU,
                                      LevelData<EBCellFAB>      & a_kappaConsDivSigmaU,
                                      LevelData<EBCellFAB>      & a_nonConsDivSigmaU,
                                      LevelData<EBCellFAB>      & a_UStar)
{
  CH_TIME("EBAMRCNS::updateEnergyBySingleLdAndRedistribute");
  //make the dissipation function = kappa*cons + (1-kappa)*noncons
  EBLevelDataOps::setToZero(a_divSigmaU);
  EBLevelDataOps::incr(a_divSigmaU, a_kappaConsDivSigmaU, 1.0);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      IntVectSet ivs = m_eblg.getEBISL()[dit()].getIrregIVS(region);
      for(VoFIterator vofit(ivs, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          Real kappa   =  m_eblg.getEBISL()[dit()].volFrac(vofit());
          Real kapCons =  a_kappaConsDivSigmaU[dit()](vofit(), 0);
          Real nonCons =    a_nonConsDivSigmaU[dit()](vofit(), 0);
          a_divSigmaU[dit()](vofit(), 0) = kapCons + (1.-kappa)*nonCons;
          m_massDiff[dit()](vofit(), CENG) = m_dt*((1.-kappa)*kapCons - kappa*(1.-kappa)*nonCons);
        }
    }
  Interval enerint(CENG, CENG);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB incr(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), 1);
      incr.setVal(0.);
      //this makes incr = dissipation function = divsigma u
      incr += a_divSigmaU[dit()];
      //this makes incr = dt*(divsigmau)
      incr *= m_dt;
      //this adds effects of dissipation into the energy
      int isrc = 0; int idst = CENG; int inco = 1;
      a_UStar[dit()].plus(incr, isrc, idst, inco);

      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), enerint);
      m_massDiff[dit()].setVal(0.);
    }

  //this smushes the energy difference into the state
  
  a_UStar.exchange(Interval(0, m_nComp-1));
  m_ebLevelRedist.resetWeights(a_UStar, CENG);

  m_ebLevelRedist.redistribute(m_redisRHS, enerint);
  m_ebLevelRedist.setToZero();
}
//U* = Un + dt*L^c(U^n) + (dt/2)*L^v(U^n)
void  
EBAMRCNS::
getUStar(LevelData<EBCellFAB>      & a_UStar, 
         const LevelData<EBCellFAB>& a_UN,
         const LevelData<EBCellFAB>& a_divergeF)
{
  CH_TIME("EBAMRCNS::getUStar");

  EBCellFactory fact(m_eblg.getEBISL());
  //advance everything explicitly
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB  dtLcU(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), m_nComp);
      dtLcU.setVal(0.);
      dtLcU += a_divergeF[dit()];
      dtLcU *= m_dt;

      a_UStar[dit()].setVal(0.);
      a_UStar[dit()] += a_UN[dit()];
      a_UStar[dit()] -= dtLcU;
    }

  m_ebLevelGodunov.floorConserved(a_UStar, m_time, m_dt);
}
/////////
void 
EBAMRCNS::
finalAdvance(LevelData<EBCellFAB>& a_Ustar)
{
  CH_TIME("EBAMRCNS::finalAdvance");
  //set stateNew to udbst
  Interval comps(0, m_nComp-1);
  a_Ustar.copyTo(comps, m_stateNew, comps);

  m_ebLevelGodunov.floorConserved(m_stateNew, m_time, m_dt);
}

/** gets divergence of shear stress
     Adds divSigma into the momentum of U*
     (\rho^* I - \dt L^m) \ubold^\npo = \dt (\rho \ubold)^*
     L^m(\ubold^\nph) = \rho^*\left(\frac{\ubold^\npo - \ubold^*}{\dt}\right)
     
  */
//============
void 
EBAMRCNS::
getDivSigma(LevelData<EBCellFAB>& a_divSigma,  
            LevelData<EBCellFAB>& a_UStar)
{
  CH_TIME("EBAMRCNS::getDivSigma");
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> velold(m_eblg.getDBL(), SpaceDim, m_nGhost*IntVect::Unit, fact);
  LevelData<EBCellFAB> velnew(m_eblg.getDBL(), SpaceDim, m_nGhost*IntVect::Unit, fact);
  LevelData<EBCellFAB> rhsZero(m_eblg.getDBL(), SpaceDim, m_nGhost*IntVect::Unit, fact);
  EBLevelDataOps::setToZero(rhsZero);

  getVelocity(   velold,  a_UStar);
  EBFluxRegister*       coarVelFRPtr = NULL;
  EBFluxRegister*       fineVelFRPtr = NULL;
  LevelData<EBCellFAB>* vCoarOldPtr = NULL;
  LevelData<EBCellFAB>* vCoarNewPtr = NULL;
  int ncomp = SpaceDim;
  Real tCoarOld = 0.0;
  Real tCoarNew = 0.0;

  if(m_hasCoarser)
    {
      // Get coarser data if we have it.
      EBAMRCNS* coarCNS = getCoarserLevel();
      coarVelFRPtr = &coarCNS->m_veloFluxRegister;
      tCoarNew = coarCNS->m_time;
      tCoarOld = tCoarNew - coarCNS->m_dt;
      const EBLevelGrid& ceblg = coarCNS->m_eblg;
      EBCellFactory cfact(ceblg.getEBISL());
      vCoarNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), ncomp, 4*IntVect::Unit, cfact);
      vCoarOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), ncomp, 4*IntVect::Unit, cfact);

      coarCNS->getVelocity(*vCoarOldPtr, coarCNS->m_stateOld);
      coarCNS->getVelocity(*vCoarNewPtr, coarCNS->m_stateNew);
    }
  if(m_hasFiner)
    {
      // Get finer flux registers if we have them.
      fineVelFRPtr = &m_veloFluxRegister;
    }

  
  Interval srcInt(CRHO, CRHO);  
  Interval dstInt(0, 0);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      (*m_acoVelo)[dit()].copy(region, dstInt, region, a_UStar[dit()], srcInt);
    }

  if(m_params.m_backwardEuler)
    {
      s_veloIntegratorBE->updateSoln(velnew, velold, rhsZero,  fineVelFRPtr, coarVelFRPtr,
                                     vCoarOldPtr, vCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt,
                                     m_level, true);
    }
  else
    {
      s_veloIntegratorTGA->updateSoln(velnew, velold, rhsZero,  fineVelFRPtr, coarVelFRPtr,
                                      vCoarOldPtr, vCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt,
                                      m_level, true);
    }

  EBLevelDataOps::scale(velold, -1./m_dt);
  EBLevelDataOps::scale(velnew,  1./m_dt);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_divSigma[dit()].setVal(0.);
      //sets divsigma = (velold - velold)/dt
      a_divSigma[dit()] += velnew[dit()];
      a_divSigma[dit()] += velold[dit()]; //see the scale above
      //sets divsigma = rho(velold - velold)/dt
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          int isrc = 0; int idst = idir; int inco = 1;
          a_divSigma[dit()].mult((*m_acoVelo)[dit()], isrc, idst, inco);
        }

      //now add divsigma*dt into ustar
      EBCellFAB dtDivSigma(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), SpaceDim);
      dtDivSigma.setVal(0.);
      dtDivSigma += a_divSigma[dit()];
      dtDivSigma *= m_dt;

      int isrc = 0; int idst = CMOMX; int inco = SpaceDim;
      a_UStar[dit()].plus(dtDivSigma, isrc, idst, inco);
    }

  //the getCoarser foo did a new.
  if(m_hasCoarser)
    {
      delete vCoarOldPtr;
      delete vCoarNewPtr;
    }
}


// (\kappa \rho^\npo C_v I - \frac{\kappa \dt}{2} L^k)T^\npo = \kappa \rho^\npo C_v T^{**}
void 
EBAMRCNS::
getDivKappaGradT(LevelData<EBCellFAB>& a_dtDivKappaGradT,
                 LevelData<EBCellFAB>& a_UStar)
{
  CH_TIME("EBAMRCNS::getDivKappaGradT");

  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB>  Told(m_eblg.getDBL(),        1,4*IntVect::Unit, fact);
  LevelData<EBCellFAB>  Tnew(m_eblg.getDBL(),        1,4*IntVect::Unit, fact);
  getTemperature(Told,  a_UStar);

  EBFluxRegister*       coarTempFRPtr = NULL;
  EBFluxRegister*       fineTempFRPtr = NULL;
  LevelData<EBCellFAB>* TCoarOldPtr = NULL;
  LevelData<EBCellFAB>* TCoarNewPtr = NULL;
  Real tCoarOld= 0;
  Real tCoarNew= 0;
  if(m_hasCoarser)
    {
      // Get coarser data if we have it.
      EBAMRCNS* coarCNS = getCoarserLevel();
      coarTempFRPtr = &coarCNS->m_tempFluxRegister;
      tCoarNew = coarCNS->m_time;
      tCoarOld = tCoarNew - coarCNS->m_dt;

      const EBLevelGrid& ceblg = coarCNS->m_eblg;
      EBCellFactory cfact(ceblg.getEBISL());
      TCoarNewPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), 1, 4*IntVect::Unit, cfact);
      TCoarOldPtr = new LevelData<EBCellFAB>(ceblg.getDBL(), 1, 4*IntVect::Unit, cfact);

      coarCNS->getTemperature(*TCoarOldPtr, coarCNS->m_stateOld);
      coarCNS->getTemperature(*TCoarNewPtr, coarCNS->m_stateNew);
    }
  if(m_hasFiner)
    {

      // Get finer flux registers if we have them.

      fineTempFRPtr = &m_tempFluxRegister;
    }
  // Set the a coefficients for the thermal conduction equation to Cv * rho, 
  // specifying both the old and new values so that the density gets linearly
  // interpolated.
  Interval srcInt(CRHO, CRHO);  
  Interval dstInt(0, 0);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      (*m_acoTemp)[dit()].copy(region, dstInt, region, a_UStar[dit()], srcInt);
      (*m_acoTemp)[dit()] *= m_params.m_specHeatCv;
    }

  LevelData<EBCellFAB> rhsTemp(m_eblg.getDBL(), 1,  m_nGhost*IntVect::Unit, fact);
  EBLevelDataOps::setToZero(rhsTemp);
  Interval srcComp(CENG, CENG);
  Interval dstComp(0,0);

  m_redisRHS.copyTo(srcComp, rhsTemp, dstComp);
  EBLevelDataOps::setToZero(m_redisRHS);

  //defineSolvers();

  if(m_params.m_backwardEuler)
    {
      s_tempIntegratorBE->updateSoln(Tnew, Told, rhsTemp,  fineTempFRPtr, coarTempFRPtr,
                                     TCoarOldPtr, TCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt/2.,
                                     m_level, true);
    }
  else
    {
      s_tempIntegratorTGA->updateSoln(Tnew, Told, rhsTemp,  fineTempFRPtr, coarTempFRPtr,
                                      TCoarOldPtr, TCoarNewPtr, m_time, tCoarOld, tCoarNew, m_dt/2.,
                                      m_level, true);

    }


  EBLevelDataOps::scale(Told, -1.0);
  EBLevelDataOps::scale(Tnew,  1.0);

  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_dtDivKappaGradT[dit()].setVal(0.);
      //sets divkgradt = ((tnew-told)/dt)*dt
      a_dtDivKappaGradT[dit()] += Told[dit()]; //neg sign is in the scale call above
      a_dtDivKappaGradT[dit()] += Tnew[dit()]; 

      //this makes divkgradt = (rho cv (tnew-told)/dt)*dt
      a_dtDivKappaGradT[dit()] *= (*m_acoTemp)[dit()];

      int isrc = 0; int idst = CENG; int inco = 1;
      a_UStar[dit()].plus(a_dtDivKappaGradT[dit()], isrc, idst, inco);
    }

  if(m_hasCoarser)
    {
      delete TCoarOldPtr;
      delete TCoarNewPtr;
    }
}

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
getVelocity(LevelData<EBCellFAB>&       a_velocity,
            const LevelData<EBCellFAB>& a_state) const
{

  CH_TIME("EBAMRCNS::getVelocity");
  Interval velInterv(0, SpaceDim-1);
  Interval momInterv(CMOMX, CMOMX+SpaceDim-1);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      a_velocity[dit()].copy(region, velInterv, region, a_state[dit()], momInterv);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          a_velocity[dit()].divide(a_state[dit()], CRHO, idir, 1);
        }
    }
}
//---------------------------------------------------------------------------------------

void
EBAMRCNS::
getTemperature(LevelData<EBCellFAB>&       a_temperature,
               const LevelData<EBCellFAB>& a_state) const
{
  CH_TIME("EBAMRCNS::getTemperature");
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      BaseFab<Real>& regTemp = a_temperature[dit()].getSingleValuedFAB();
      const BaseFab<Real>& regStat = a_state[dit()].getSingleValuedFAB();
      FORT_GETTEMPERATURE(CHF_FRA1(regTemp,0),
                          CHF_CONST_FRA(regStat),
                          CHF_BOX(region),
                          CHF_CONST_REAL(m_params.m_specHeatCv));
      IntVectSet ivsMulti = m_eblg.getEBISL()[dit()].getMultiCells(region);
      for(VoFIterator vofit(ivsMulti, m_eblg.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          Vector<Real> cvec(CNUM);
          for(int ivar = 0; ivar< CNUM; ivar++)
            {
              cvec[ivar] = a_state[dit()](vofit(), ivar);
            }
          Real temper;
          FORT_POINTGETTEMPERATURE(CHF_REAL(temper),
                                   CHF_VR(cvec),
                                   CHF_CONST_REAL(m_params.m_specHeatCv));

          a_temperature[dit()](vofit(), 0)  = temper;
        }
    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
getHalfState(LevelData<EBCellFAB>& a_stateInt)
{
  //interpolate state to n+1/2
  Real told = 0; Real tnew = 1; Real time = 0.5;
  EBArith::timeInterpolate(a_stateInt, m_stateOld, m_stateNew,
                           m_eblg.getDBL(), time, told, tnew);
}
//---------------------------------------------------------------------------------------
//set  output to div (sigma)
void
EBAMRCNS::
kappaMomentumSource(LevelData<EBCellFAB>& a_kappaDivSigma,
                    const LevelData<EBCellFAB>& a_velocity,
                    const LevelData<EBCellFAB>* a_veloCoar,
                    const LevelData<EBCellFAB>& a_state)//for coefficients
{
  //set the a coefficients to time n
  Interval srcInt(CRHO, CRHO);  
  Interval dstInt(0, 0);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      (*m_acoVelo)[dit()].copy(region, dstInt, region, a_state[dit()], srcInt);
    }

  Real alpha = 0;   Real beta = 1; //want just the div(flux) part of the operator
  // Compute the viscous diffusion term.  coefficient is unity because we want the straight operator
  s_veloIntegratorBE->applyOperator(a_kappaDivSigma,
                                    a_velocity, 
                                    a_veloCoar,
                                    m_level, alpha, beta, true);
}

//---------------------------------------------------------------------------------------
/// set  output to (del dot (kappa grad T))  + div(sigma u)
void
EBAMRCNS::
kappaEnergySource(LevelData<EBCellFAB>& a_kappaEnergySource, 
                  const LevelData<EBCellFAB>& a_velocity,
                  const LevelData<EBCellFAB>* a_veloCoar,
                  const LevelData<EBCellFAB>& a_temperat,
                  const LevelData<EBCellFAB>* a_tempCoar,
                  const LevelData<EBCellFAB>& a_state)//for coefficients
{
  
  CH_TIME("EBAMRCNS::kappaEnergySource");
  //set the a coefficients to time n
  Interval srcInt(CRHO, CRHO);  
  Interval dstInt(0, 0);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& region = m_eblg.getDBL().get(dit());
      (*m_acoTemp)[dit()].copy(region, dstInt, region, a_state[dit()], srcInt);
      (*m_acoTemp)[dit()] *= m_params.m_specHeatCv;
    }


  // Compute the heat conduction term = volfrac(div kappa grad T)
  Real alpha = 0;   Real beta = 1; //want just the div(flux) part of the operator
  s_tempIntegratorBE->applyOperator(a_kappaEnergySource, a_temperat,  a_tempCoar,
                                    m_level, alpha, beta,  true);

  /**/
  /// add volfrac*div( sigma u) to a_kappaEnergySource 
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB>  kappaDivSigmaU(m_eblg.getDBL(), 1,4*IntVect::Unit, fact);
  s_veloIntegratorBE->getKappaDivSigmaU(kappaDivSigmaU,
                                         a_velocity, 
                                         a_veloCoar,
                                         m_level);

  
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      a_kappaEnergySource[dit()] += kappaDivSigmaU[dit()];
    }
  /**/

}
//this includes news so need to call delete
//---------------------------------------------------------------------------------------
void
EBAMRCNS::
getCoarserTemperatureAndVelocity(LevelData<EBCellFAB>* & a_tempCoar,
                                 LevelData<EBCellFAB>* & a_veloCoar)
{
  
  CH_TIME("EBAMRCNS::getCoarserTemperatureAndVelocity");
  LevelData<EBCellFAB>* veloCoar = NULL;
  LevelData<EBCellFAB>* tempCoar = NULL;
  if(m_hasCoarser)
    {
      EBAMRCNS* coarPtr = getCoarserLevel();
      EBCellFactory factCoar(coarPtr->m_eblg.getEBISL());
      veloCoar = new LevelData<EBCellFAB>(coarPtr->m_eblg.getDBL(), SpaceDim,4*IntVect::Unit, factCoar);
      LevelData<EBCellFAB>    veloCoarOld(coarPtr->m_eblg.getDBL(), SpaceDim,4*IntVect::Unit, factCoar);
      LevelData<EBCellFAB>    veloCoarNew(coarPtr->m_eblg.getDBL(), SpaceDim,4*IntVect::Unit, factCoar);
      tempCoar = new LevelData<EBCellFAB>(coarPtr->m_eblg.getDBL(),        1,4*IntVect::Unit, factCoar);
      LevelData<EBCellFAB>    tempCoarOld(coarPtr->m_eblg.getDBL(),        1,4*IntVect::Unit, factCoar);
      LevelData<EBCellFAB>    tempCoarNew(coarPtr->m_eblg.getDBL(),        1,4*IntVect::Unit, factCoar);

      coarPtr->getVelocity(   veloCoarOld, coarPtr->m_stateOld);
      coarPtr->getVelocity(   veloCoarNew, coarPtr->m_stateNew);
      coarPtr->getTemperature(tempCoarOld, coarPtr->m_stateOld);
      coarPtr->getTemperature(tempCoarNew, coarPtr->m_stateNew);

      Real tCoarNew = coarPtr->m_time;
      Real tCoarOld = tCoarNew - coarPtr->m_dt;

      //interpolate coarse solution to fine time
      EBArith::timeInterpolate(*veloCoar, veloCoarOld, veloCoarNew, coarPtr->m_eblg.getDBL(), m_time, tCoarOld, tCoarNew);        
      EBArith::timeInterpolate(*tempCoar, tempCoarOld, tempCoarNew, coarPtr->m_eblg.getDBL(), m_time, tCoarOld, tCoarNew);        
    }

  a_tempCoar = tempCoar;
  a_veloCoar = veloCoar;
}

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
explicitHyperbolicSource(LevelData<EBCellFAB>&       a_momentSource, 
                         LevelData<EBCellFAB>&       a_energySource,
                         const LevelData<EBCellFAB>& a_state,
                         bool a_doNormalization)
{
  CH_TIME("EBAMRCNS::explicit_hyperbolic_source");
  pout() << "computing explicit hyperbolic source" << endl;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB>  velocity(m_eblg.getDBL(), SpaceDim,4*IntVect::Unit, fact);
  LevelData<EBCellFAB>  temperat(m_eblg.getDBL(),        1,4*IntVect::Unit, fact);
  LevelData<EBCellFAB>* veloCoar = NULL;
  LevelData<EBCellFAB>* tempCoar = NULL;
  getVelocity(   velocity,  a_state);
  getTemperature(temperat,  a_state);
  getCoarserTemperatureAndVelocity(tempCoar, veloCoar);

  //velocity is a part of both source terms.  temperature only comes into the energy source term.
  kappaMomentumSource(a_momentSource, velocity, veloCoar,                     a_state);
  kappaEnergySource(  a_energySource, velocity, veloCoar, temperat, tempCoar, a_state);

  //the getcoarser stuff did a new
  if(m_hasCoarser)
    {
      delete veloCoar;
      delete tempCoar;
    }

  if(a_doNormalization)
    {
      //finally all these things are multiplied by kappa so we need to make 
      //them normalized 
      KappaSquareNormal normalizinator(m_eblg);
      normalizinator(a_momentSource);
      normalizinator(a_energySource);
    }
}
//---------------------------------------------------------------------------------------
void
EBAMRCNS::
hyperbolicSource(LevelData<EBCellFAB>&       a_source)
{
  CH_TIME("EBAMRCNS::hyperbolicsource");
  bool zeroHyperbolicSource = false;

  // begin debug 
  zeroHyperbolicSource = true;
  //set source to zero
  if(zeroHyperbolicSource)
    {
      pout() << "zeroing out hyperbolic source"  << endl;
      EBLevelDataOps::setToZero(a_source);
    }
  //end debug
 else
    {
      EBPatchGodunov::useConservativeSource(true);
      //this is important because sometimes there is no diffusion
      EBLevelDataOps::setVal(a_source, 0.0);
      if(m_params.m_doDiffusion)
        {
          //get all the primitives (velocities, temps) 
          //we need from this level and the next coarser level
          EBCellFactory fact(m_eblg.getEBISL());
          LevelData<EBCellFAB>  momeSour(m_eblg.getDBL(), SpaceDim,4*IntVect::Unit, fact);
          LevelData<EBCellFAB>  enerSour(m_eblg.getDBL(),        1,4*IntVect::Unit, fact);
          explicitHyperbolicSource(momeSour, enerSour, m_stateOld, true);

          Interval momeInterv(CMOMX, CMOMX+SpaceDim-1);
          Interval enerInterv(CENG,  CENG);
          Interval vectInterv(0, SpaceDim-1);
          Interval scalInterv(0, 0);
          for(DataIterator dit = a_source.dataIterator(); dit.ok(); ++dit)
            {
              const Box& region = m_eblg.getDBL().get(dit());
              a_source[dit()].copy(region, momeInterv, region, momeSour[dit()], vectInterv);
              a_source[dit()].copy(region, enerInterv, region, enerSour[dit()], scalInterv);

            }
          a_source.exchange();
        }

    }  
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
fluxDivergence( LevelData<EBCellFAB>& a_divergeF)
{
  CH_TIME("EBAMRCNS::fluxdivergence");
  pout() << "taking flux divergence on level " << m_level << " with dt = " << m_dt << endl;
  IntVect ivGhost = 4*IntVect::Unit;
  EBCellFactory fact(m_eblg.getEBISL());
  LevelData<EBCellFAB> source(m_eblg.getDBL(),m_nComp, 4*IntVect::Unit, fact);
  //set up arguments to step
  //undefined lfr in case we need it
  EBFluxRegister lfr;
  //undefined leveldata in case we need it
  const LevelData<EBCellFAB> ld;
  //set arguments to dummy arguments and
  //then fix if they are available
  EBFluxRegister* coarFR = &lfr;
  EBFluxRegister* fineFR = &lfr;
  const LevelData<EBCellFAB>* coarDataOld = &ld;
  const LevelData<EBCellFAB>* coarDataNew = &ld;

  Real told = 0.0;
  Real tnew = 0.0;
  
  
  if(m_hasCoarser)
    {
      EBAMRCNS* coarPtr = getCoarserLevel();
      //recall that my flux register goes between my
      //level and the next finer level
      coarFR = &coarPtr->m_divFFluxRegister;
      coarDataOld = &coarPtr->m_stateOld;
      coarDataNew = &coarPtr->m_stateNew;
      tnew = coarPtr->m_time;
      told = tnew - coarPtr->m_dt;

      // At this level, m_time theoretically can not be smaller than its coarser level time (told); 
      // if it happends due to the precision/floating error, then we resign m_time to be told. 
      // One alternative way is to ensure a tolerance check. However, this current fix is intend to 
      // avoid further error accumulation (based on the choice of tolerance). 
      // Ideally, one would want to have a "generic" fix in either AMR.cpp or postTimeStep().
      if (m_time < told)
        {
          m_time = told;
        } 
      
    }
  if(m_hasFiner)
    {
      //recall that my flux register goes between my
      //level and the next finer level
      fineFR = &m_divFFluxRegister;
    }

  hyperbolicSource(source);
  // flux register manipulation happens in levelgodunov

  //Stuff needed for the energy equation 
  //also come out of here since we need the half time
  //values  of the primitive state.    Not the most intuitive
  //interface, I know.
 
  m_ebLevelGodunov.divergeF(a_divergeF,
                            m_massDiff,
                            *fineFR,
                            *coarFR,
                            m_stateOld,
                            source,
                            *coarDataOld,
                            *coarDataNew,
                            m_time,
                            told,
                            tnew,
                            m_dt);
  coarseFineIncrement();
}
//---------------------------------------------------------------------------------------
void
EBAMRCNS::
coarseFineIncrement()
{
  CH_TIME("EBAMRCNS::coarseFineIncrement");
  Interval interv(0, m_nComp-1);
  // increment redistribution register between this level
  // and next coarser level
  {
    CH_TIME("coarse-fine guano");
    if(m_hasCoarser && (!s_noEBCF))
      {
        for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
          {
            m_ebFineToCoarRedist.increment(m_massDiff[dit()], dit(), interv);
          }
      }

    //initialize redistribution register between this level
    // and next finer level
    //this includes re-redistribution registers
    if(m_hasFiner  && (!s_noEBCF))
      {

        m_ebCoarToFineRedist.setToZero();
        m_ebCoarToCoarRedist.setToZero();
        for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
          {
            BaseIVFAB<Real>& massDiffFAB = m_massDiff[dit()];
            m_ebCoarToFineRedist.increment(massDiffFAB, dit(), interv);
            m_ebCoarToCoarRedist.increment(massDiffFAB, dit(), interv);
          }
      }
  }
}
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
void
EBAMRCNS::
hyperbolicRedistribution(LevelData<EBCellFAB>& a_state)
{
  CH_TIME("EBAMRCNS::hyperbolic redistribution");
  //explicit redistribution here because these are the hyperbolic terms
  m_ebLevelRedist.setToZero();
  if(m_params.m_useMassRedist)
    {
      //if use mass weighting, need to
      //fix weights of redistribution object
      int densityIndex = m_ebPatchGodunov->densityIndex();
      a_state.exchange(Interval(0, m_nComp-1));
      m_ebLevelRedist.resetWeights(a_state, densityIndex);
    }

  Interval consInterv(0, m_nComp-1);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_ebLevelRedist.increment(m_massDiff[dit()], dit(), consInterv);
    }
  //
  m_ebLevelRedist.redistribute(a_state, consInterv);

  m_ebLevelRedist.setToZero();
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
explicitAdvance(const LevelData<EBCellFAB>& a_divergeF)
{
  CH_TIME("EBAMRCNS::explicitAdvance");
  CH_assert(!m_params.m_doDiffusion);
  //advance everything explicitly
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB dtDivergeF(m_eblg.getEBISL()[dit()], m_eblg.getDBL().get(dit()), a_divergeF.nComp());
      dtDivergeF.setVal(0.);
      dtDivergeF += a_divergeF[dit()];
      dtDivergeF *= m_dt;
      m_stateNew[dit()] -= dtDivergeF;
    }
  hyperbolicRedistribution(m_stateNew);

  m_ebLevelGodunov.floorConserved(m_stateNew, m_time, m_dt);
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
int 
EBAMRCNS::
getFinestLevel()
{
  int imax = 0;
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  for(int ilev = 0; ilev < hierarchy.size(); ilev++)
    {
      EBAMRCNS* cnsLevel = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
      int numBoxes = 0;
      if(cnsLevel->m_eblg.getDBL().isClosed()) numBoxes = cnsLevel->m_eblg.getDBL().size();
      if(numBoxes > 0) 
        imax++;
      else
        break;
    }
  return imax-1;
}
//---------------------------------------------------------------------------------------
void
EBAMRCNS::
postTimeStep()
{
  CH_TIME("EBAMRCNS::postTimeStep");
  pout() << " in EBAMRCNS postTimeStep for level " << m_level << endl;

 // Average the finer grid levels to be consistent with this one.
 
  if (m_hasFiner)
    {
      Interval interv(0, m_nComp-1);
      EBAMRCNS* finePtr = getFinerLevel();
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);
    }

  //this does the refluxing and redistribution evil dance
  postTimeStepRefluxRedistDance();
  m_ebLevelGodunov.floorConserved(m_stateNew, m_time, m_dt);

  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
      m_redisRHS[dit()].setVal(0.0);
    }

  if(m_params.m_checkMaxMin)
    {
      CH_TIME("max_min_check_post");
      pout() << "after post time step max mins for data for level " << m_level << endl;
      for(int icomp = 0; icomp < m_stateNew.nComp(); icomp++)
        {
          Real maxVal, minVal;
          EBLevelDataOps::getMaxMin(maxVal, minVal, m_stateNew, icomp);
          pout() << " "
                 << setprecision(4)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 <<  "comp = " <<  icomp 
                 << ", max = " << maxVal 
                 << ", min = " << minVal << endl;
        }
    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
postTimeStepRefluxRedistDance()
{
  CH_TIME("EBAMRCNS::evilRefluxRedistDance");
  resetWeights();
  refluxRedistInteraction();
  coarseFineRedistribution(Interval(0, m_nComp-1));
  explicitReflux(Interval(0, m_nComp-1));
  ParmParse pp;
  bool turn_off_reflux = false;
  pp.query("turn_off_implicit_reflux", turn_off_reflux);
  if(turn_off_reflux)
    {
      pout() << "implicit reflux turned off" << endl;
    }
  if(m_params.m_doDiffusion && !turn_off_reflux)
    {
      refluxRHSEnergyAndMomentum();
      //      defineSolvers();
      implicitReflux();
    }
}
//---------------------------------------------------------------------------------------
void 
EBAMRCNS::
implicitReflux()
{
  CH_TIME("EBAMRCNS::implicitRedistReflux");
  pout() << "EBAMRCNS::implicit redist/reflux for level" << m_level << endl;
  Real crseTime = -1.0;
  Real timeEps = 1.0e-2*m_dt;
  if (m_level > 0) crseTime = m_coarser_level_ptr->time();
  
  if (m_level == 0 || (abs(crseTime - m_time) < timeEps))
    {
       int finestLev = getFinestLevel();
      //if I have not caught up to the next coarser level,
      //then I need to do implicit{reflux,redist} for my level and all
      //finer levels
      Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
      Vector<LevelData<EBCellFAB>* >    deltaVelocity(finestLev+1, NULL);
      Vector<LevelData<EBCellFAB>* >    deltaTemperat(finestLev+1, NULL);
      Vector<LevelData<EBCellFAB>* > dtRefluxDivergeM(finestLev+1, NULL);
      Vector<LevelData<EBCellFAB>* > dtRefluxDivergeE(finestLev+1, NULL);
      for(int ilev = 0; ilev <= finestLev; ilev++)
        {
          EBAMRCNS* cnsLevel = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);

          EBCellFactory fact(cnsLevel->m_eblg.getEBISL());
          deltaVelocity[ilev]    = new LevelData<EBCellFAB>(cnsLevel->m_eblg.getDBL(), SpaceDim, 4*IntVect::Unit, fact);
          deltaTemperat[ilev]    = new LevelData<EBCellFAB>(cnsLevel->m_eblg.getDBL(),        1, 4*IntVect::Unit, fact);
          dtRefluxDivergeE[ilev] = new LevelData<EBCellFAB>(cnsLevel->m_eblg.getDBL(),        1, 4*IntVect::Unit, fact);
          dtRefluxDivergeM[ilev] = new LevelData<EBCellFAB>(cnsLevel->m_eblg.getDBL(), SpaceDim, 4*IntVect::Unit, fact);

          //copy redistRHS to reflux divergence holders
          Interval srcComp, dstComp;
          srcComp = Interval(CMOMX, CMOMX+SpaceDim-1);
          dstComp = Interval(0, SpaceDim-1);
          cnsLevel->m_redisRHS.copyTo(srcComp, *dtRefluxDivergeM[ilev], dstComp);
          srcComp = Interval(CENG, CENG);
          dstComp = Interval(0, 0);
          cnsLevel->m_redisRHS.copyTo(srcComp, *dtRefluxDivergeE[ilev], dstComp);
          EBLevelDataOps::scale(*dtRefluxDivergeM[ilev], cnsLevel->m_dt);
          EBLevelDataOps::scale(*dtRefluxDivergeE[ilev], cnsLevel->m_dt);

          //these MATTER especially with subcycling
          EBLevelDataOps::setVal(*deltaVelocity[ilev], 0.0);
          EBLevelDataOps::setVal(*deltaTemperat[ilev], 0.0);
        }

      pout() << "getting momentum reflux/redist increment"  << endl;
      int baseLev = Max(m_level-1, 0);
      Real baseDt = ((EBAMRCNS*)hierarchy[baseLev])->m_dt;

      //(rho I - dt Lv) delta = dt* Dr(Fm)
      getRefluxDeltaV(  deltaVelocity, dtRefluxDivergeM, baseLev, finestLev, baseDt);

      //mom += dt* Lv(deltav)
      incrMomentumByDeltaV(deltaVelocity, baseLev, finestLev);

      pout() << "getting energy reflux/redist increment"  << endl;
      //(rho I - dt Lt) delta =dt* Dr(Fe)
      getRefluxDeltaT(deltaTemperat, dtRefluxDivergeE, baseLev, finestLev, baseDt);

      //ene += dt* Lt(deltaT)
      incrEnergyByDeltaT(deltaTemperat, baseLev, finestLev);

      for(int ilev = 0; ilev <= finestLev; ilev++)
        {
          delete    deltaVelocity[ilev];
          delete    deltaTemperat[ilev];
          delete dtRefluxDivergeM[ilev];
          delete dtRefluxDivergeE[ilev];
        }

      //reset redistribution RHS to avoid double counting
      for(int ilev = m_level; ilev <= finestLev; ilev++)
        {
          EBAMRCNS* cnsLevel = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
          EBLevelDataOps::setToZero(cnsLevel->m_redisRHS);
        }
    }
}
//---------------------------------------------------------------------------------------
//(rho I - dt Lv) delta = dt*Dr(Fm)
void 
EBAMRCNS::
getRefluxDeltaV(Vector<LevelData<EBCellFAB>* >&  a_deltaVelocity, 
                      Vector<LevelData<EBCellFAB>* >&  a_dtRefluxDivergeM, 
                      int  a_baseLev, int a_finestLev, Real a_dtBase)
{
  CH_TIME("EBAMRCNS::getRefluxDeltaV");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  //set identity coefficients
  for(int ilev = 0; ilev <= a_finestLev; ilev++)
    {
      EBAMRCNS* cnsLevel = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
      LevelData<EBCellFAB>& state =  cnsLevel->m_stateNew;
      LevelData<EBCellFAB>& acoef = *cnsLevel->m_acoVelo;
      Interval srcInt(CRHO, CRHO);  
      Interval dstInt(0, 0);
      for(DataIterator dit = cnsLevel->m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          const Box& region = cnsLevel->m_eblg.getDBL().get(dit());
          acoef[dit()].copy(region, dstInt, region, state[dit()], srcInt);
        }
    }

  //set alpha to 1, beta = -dtbase
  s_veloIntegratorBE->resetSolverAlphaAndBeta(1.0, -a_dtBase);

  //solve equation (rho I - dt Lv) delta = dt*Dr(Fm)
  //rhs already multiplied by dt
  //first true = zero phi
  //second true = force homogeneous bcs (solving for a delta t)
  s_veloSolver->solve(a_deltaVelocity, a_dtRefluxDivergeM, a_finestLev, a_baseLev, true, true);
}
//---------------------------------------------------------------------------------------
//(rho Cv I - dt Lt) delta = dt*Dr(Fe)
void 
EBAMRCNS::
getRefluxDeltaT(Vector<LevelData<EBCellFAB>* >& a_deltaTemperat,
                      Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeE, 
                      int a_baseLev, int a_finestLev, Real a_dtBase)
{
  CH_TIME("EBAMRCNS::getRefluxDeltaT");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  //set identity coefficients
  for(int ilev = 0; ilev <= a_finestLev; ilev++)
    {
      EBAMRCNS* cnsLevel = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
      Interval srcInt(CRHO, CRHO);  
      Interval dstInt(0, 0);
      LevelData<EBCellFAB>& state =  cnsLevel->m_stateNew;
      LevelData<EBCellFAB>& acoef = *cnsLevel->m_acoTemp;
      for(DataIterator dit = cnsLevel->m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          const Box& region = cnsLevel->m_eblg.getDBL().get(dit());
          acoef[dit()].copy(region, dstInt, region, state[dit()], srcInt);
          acoef[dit()] *= m_params.m_specHeatCv;
        }
    }
  //set alpha to 1, beta = -dtbase
  s_tempIntegratorBE->resetSolverAlphaAndBeta(1.0, -a_dtBase);

  //solve equation (rho C_v I - dt Lv) delta = dt*Dr(Fm)
  //rhs already multiplied by dt
  //first true = zero phi
  //second true = force homogeneous bcs (solving for a delta t)
  s_tempSolver->solve(a_deltaTemperat, a_dtRefluxDivergeE, a_finestLev, a_baseLev, true, true);
}
//---------------------------------------------------------------------------------------
//mom += dt* Lv(deltav)
// deltaM gets multiplied by rho 
void 
EBAMRCNS::
incrMomentumByDeltaV(Vector<LevelData<EBCellFAB>* >& a_deltaVelocity, 
                     int a_baseLev, int a_finestLev)
{
  CH_TIME("EBAMRCNS::incrMomentumByDeltaV");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  for(int ilev = a_baseLev; ilev <= a_finestLev; ilev++)
    {
      EBAMRCNS* cnsLevel   = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
      cnsLevel->incrMomentumByDeltaV( *a_deltaVelocity[ilev], cnsLevel->m_stateNew);
    }
}
//---------------------------------------------------------------------------------------
void 
EBAMRCNS::
incrMomentumByDeltaV(LevelData<EBCellFAB>& a_deltaVelocity, 
                     LevelData<EBCellFAB>& a_state)

{
  CH_TIME("EBAMRCNS::incrMomentumByDeltaVLevel");
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& b = m_eblg.getDBL()[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      EBCellFAB deltaMom(ebisBox, b, SpaceDim);
      deltaMom.setVal(0.);
      deltaMom += a_deltaVelocity[dit()];

      EBCellFAB& stateNew  = a_state[dit()];
      //delta mom now = dV 
      //make delta = rho dV
      for(int idir = 0; idir <SpaceDim; idir++)
        {
          int isrc = CRHO; int idst = idir; int inco = 1;
          deltaMom.mult(stateNew, isrc, idst, inco);
        }
      //add change in momentum to the state
      {
        int isrc = 0; int idst = CMOMX; int inco = SpaceDim;
        stateNew.plus(deltaMom, isrc, idst, inco);
      }
    }
}
//---------------------------------------------------------------------------------------
//ene += dt* Lt(deltaT)
// deltaT gets multiplied by rho cv 
void 
EBAMRCNS::
incrEnergyByDeltaT(Vector<LevelData<EBCellFAB>* >& a_deltaTemperat, 
                   int a_baseLev, int a_finestLev)
{
  CH_TIME("EBAMRCNS::incrEnergyByDeltaT");
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  for(int ilev = a_baseLev; ilev <= a_finestLev; ilev++)
    {
      EBAMRCNS* cnsLevel   = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
      cnsLevel->incrEnergyByDeltaT(cnsLevel->m_stateNew, *a_deltaTemperat[ilev]);
    }
} 

//---------------------------------------------------------------------------------------
//ene += dt* Lt(deltaT)
// deltaT gets multiplied by rho cv 
void 
EBAMRCNS::
incrEnergyByDeltaT(LevelData<EBCellFAB>& a_state,
                   LevelData<EBCellFAB>& a_deltaTemperat)
{
  CH_TIME("EBAMRCNS::incrEnergyByDeltaTLevel");
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& b = m_eblg.getDBL()[dit()];
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      EBCellFAB deltaEnergy(ebisBox, b, 1);
      deltaEnergy.setVal(0.);
      deltaEnergy += a_deltaTemperat[dit()];

      //deltaEnergy now = dT 
      //make deltaEnergy = rho  Cv dT
      {
        int isrc = CRHO; int idst = 0; int inco = 1;
        deltaEnergy.mult(a_state[dit()], isrc, idst, inco);
        deltaEnergy *= m_params.m_specHeatCv;
      }
      //add change in energy to the state
      {
        int isrc = 0; int idst = CENG; int inco = 1;
        a_state[dit()].plus(deltaEnergy, isrc, idst, inco);
      }
    }
}
//---------------------------------------------------------------------------------------
void 
EBAMRCNS::
refluxRedistInteraction()
{
  //the flux register must modify the redistribution
  //registers
  CH_TIME("EBAMRCNS::refluxRedistInteraction");
  if(m_hasFiner  && (!s_noEBCF))
    {
      Real scale = -1.0/m_dx[0];
      Interval interv(0, m_nComp-1);
      m_divFFluxRegister.incrementRedistRegister(m_ebCoarToFineRedist,
                                                 interv, scale);

      m_divFFluxRegister.incrementRedistRegister(m_ebCoarToCoarRedist,
                                                 interv, scale);
    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
resetWeights()
{
  CH_TIME("EBAMRCNS::resetWeights");
  //if use mass weighting, need to
  //fix weights of redistribution objects
  if(m_params.m_useMassRedist)
    {
      int densevar = m_ebPatchGodunov->densityIndex();
      m_stateNew.exchange(Interval(0, m_nComp-1));
      m_ebLevelRedist.resetWeights(m_stateNew, densevar);

      if(m_hasCoarser && (!s_noEBCF))
        {
          EBAMRCNS* coarPtr = getCoarserLevel();
          coarPtr->m_stateNew.exchange(Interval(0, m_nComp-1));
          m_ebFineToCoarRedist.resetWeights(coarPtr->m_stateNew, densevar);
        }
      if(m_hasFiner && (!s_noEBCF))
        {
          m_ebCoarToFineRedist.resetWeights(m_stateNew, densevar);
          m_ebCoarToCoarRedist.resetWeights(m_stateNew, densevar);
        }
    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
explicitReflux(const Interval & a_interv)
{
  CH_TIME("EBAMRCNS::explicitReflux");
  if (m_hasFiner)
    {
      Real scale = -1.0/m_dx[0];
      m_divFFluxRegister.reflux(m_stateNew, a_interv, scale);
      m_divFFluxRegister.setToZero();
    }
}
//---------------------------------------------------------------------------------------
void
EBAMRCNS::
coarseFineRedistribution(const Interval & a_interv)
{
  CH_TIME("EBAMRCNS::coarseFineRedistribution");
  if(m_hasCoarser)
    {
      if(m_params.m_doSmushing && (!s_noEBCF))
        {
          //redistibute to coarser level
          EBAMRCNS* coarPtr = getCoarserLevel();
          //put mass directly into state
          m_ebFineToCoarRedist.redistribute(coarPtr->m_stateNew, a_interv);
          m_ebFineToCoarRedist.setToZero();
        }
    }
  if (m_hasFiner)
    {
      EBAMRCNS* finePtr = getFinerLevel();
      if(m_params.m_doSmushing  && (!s_noEBCF))
        {
          //redistibute to finer level
          m_ebCoarToFineRedist.redistribute(finePtr->m_stateNew, a_interv);
          //do the re redistirubtion
          m_ebCoarToCoarRedist.redistribute(         m_stateNew, a_interv);
          m_ebCoarToFineRedist.setToZero();
          m_ebCoarToCoarRedist.setToZero();
        }
    }
}
//---------------------------------------------------------------------------------------

void
EBAMRCNS::
refluxRHSEnergyAndMomentum()
{
  CH_TIME("EBAMRCNS::refluxAndRHSEnergyAndMomentum");
  //this does the refluxing and redistribution evil dance
  Interval interm(CMOMX, CMOMX+SpaceDim-1);
  Interval intere(CENG, CENG);
  if (m_hasFiner)
    {
      // reflux from finer level solution
      CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
      Real scale = -1.0/m_dx[0];
      Interval interfluxm(0, SpaceDim-1);
      Interval interfluxe(0, 0);
      m_veloFluxRegister.reflux(m_redisRHS, interm, interfluxm, scale);
      m_tempFluxRegister.reflux(m_redisRHS, intere, interfluxe, scale);

      m_veloFluxRegister.setToZero();
      m_tempFluxRegister.setToZero();
    }
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
tagCells(IntVectSet& a_tags)
{
  CH_TIME("EBAMRCNS::tagCells");
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::tagCells for level " << m_level << endl;
    }

  // Create tags based on undivided gradient of density
  IntVectSet localTags;
  ParmParse pp;
  bool trackShock = false;
  int imin = -1;
  if(pp.contains("track_shock"))
    {
      pp.get("track_shock", trackShock);
    }
  if(trackShock)
    {
      Real trailBuffer;
      pp.get("trail_buffer", trailBuffer);
      CH_assert((trailBuffer > 0) && (trailBuffer < 1.0));
      Real rtrail = trailBuffer*m_problem_domain.domainBox().size(0);
      int itrail = Max(1, int(rtrail));

      Real shocklocx;
      FORT_GETSHOCKLOCX(CHF_REAL(shocklocx), CHF_CONST_REAL(m_time));
      int iloc = shocklocx/m_dx[0];
      imin = iloc - itrail;
    }

  vector<Real> problo(SpaceDim);
  // If there is a coarser level interpolate undefined ghost cells
  //only interpolate the density

  bool tagOnEnergy = false;
  pp.query("tag_on_energy", tagOnEnergy);
  int densityIndex = CRHO;
  if(tagOnEnergy)
    densityIndex = CENG;

  Interval intervDensity(densityIndex, densityIndex);
  EBCellFactory factory(m_eblg.getEBISL());
  int nCons = m_ebPatchGodunov->numConserved();
  LevelData<EBCellFAB> consTemp(m_eblg.getDBL(), nCons, IntVect::Unit, factory);
  Interval consInterv(0, nCons-1);
  m_stateNew.copyTo(consInterv, consTemp, consInterv);
  if (m_hasCoarser)
    {
      const EBAMRCNS* coarCNS = getCoarserLevel();
      int refRatCrse = coarCNS->refRatio();
      int nghost = 1;
      EBPWLFillPatch patcher(m_eblg.getDBL(),
                             coarCNS->m_eblg.getDBL(),
                             m_eblg.getEBISL(),
                             coarCNS->m_eblg.getEBISL(),
                             coarCNS->m_eblg.getDomain().domainBox(),
                             refRatCrse, m_nComp, nghost);

      Real coarTimeOld = 0.0;
      Real coarTimeNew = 1.0;
      Real fineTime    = 0.0;
      patcher.interpolate(consTemp,
                          coarCNS->m_stateOld,
                          coarCNS->m_stateNew,
                          coarTimeOld,
                          coarTimeNew,
                          fineTime,
                          intervDensity);
    }
  consTemp.exchange(intervDensity);

  // Compute undivided gradient
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const Box& b = m_eblg.getDBL().get(dit());
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      EBCellFAB gradFab(ebisBox, b, SpaceDim);
      const EBCellFAB& stateFab = consTemp[dit()];
      BaseFab<Real>& regGradFab = gradFab.getSingleValuedFAB();
      const BaseFab<Real>& regStateFab = stateFab.getSingleValuedFAB();

      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          const Box bCenter = b & grow(m_eblg.getDomain(),-BASISV(idir));
          const Box bLo     = b & adjCellLo(bCenter,idir);
          const int hasLo = ! bLo.isEmpty();
          const Box bHi     = b & adjCellHi(bCenter,idir);
          const int hasHi = ! bHi.isEmpty();

          FORT_GETRELATIVEGRAD(CHF_FRA1(regGradFab,idir),
                               CHF_CONST_FRA1(regStateFab,densityIndex),
                               CHF_CONST_INT(idir),
                               CHF_BOX(bLo),
                               CHF_CONST_INT(hasLo),
                               CHF_BOX(bHi),
                               CHF_CONST_INT(hasHi),
                               CHF_BOX(bCenter));

          //do one-sided diffs where necessary at irregular cells.
          IntVectSet ivsIrreg = ebisBox.getIrregIVS(b);
          for(VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph());
              vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const IntVect&  iv = vof.gridIndex();
              //one-sided diffs on domain bndry
              bool onLeftDomain = iv[idir] == m_eblg.getDomain().domainBox().smallEnd(idir);
              bool onRighDomain = iv[idir] == m_eblg.getDomain().domainBox().bigEnd(idir);
              bool hasFacesLeft = (ebisBox.numFaces(vof, idir, Side::Lo) > 0) && !onLeftDomain;
              bool hasFacesRigh = (ebisBox.numFaces(vof, idir, Side::Hi) > 0) && !onRighDomain;

              Real valCent = stateFab(vof, densityIndex);
              Real dpl = 0.;
              Real dpr = 0.;
              Real dpc = 0.;

              //compute one-sided diffs where you have them
              if(hasFacesLeft)
                {
                  Vector<FaceIndex> facesLeft =
                    ebisBox.getFaces(vof, idir, Side::Lo);
                  Real valLeft = 0.0;
                  for(int iface = 0; iface <facesLeft.size(); iface++)
                    {
                      VolIndex vofLeft = facesLeft[iface].getVoF(Side::Lo);
                      valLeft += stateFab(vofLeft, densityIndex);
                    }
                  valLeft /= Real(facesLeft.size());
                  dpl = valCent - valLeft;
                }
              if(hasFacesRigh)
                {
                  Vector<FaceIndex> facesRigh =
                    ebisBox.getFaces(vof, idir, Side::Hi);
                  Real valRigh = 0.0;
                  for(int iface = 0; iface <facesRigh.size(); iface++)
                    {
                      VolIndex vofRigh = facesRigh[iface].getVoF(Side::Hi);
                      valRigh += stateFab(vofRigh, densityIndex);
                    }
                  valRigh /= Real(facesRigh.size());
                  dpr = valRigh - valCent;
                }
              if(hasFacesLeft && hasFacesRigh)
                {
                  dpc = 0.5*(dpl+dpr);
                }
              else if(!hasFacesLeft && !hasFacesRigh)
                {
                  dpc = 0.0;
                }
              else if(hasFacesLeft && !hasFacesRigh)
                {
                  dpc = dpl;
                }
              else if(hasFacesRigh && !hasFacesLeft)
                {
                  dpc = dpr;
                }

              gradFab(vof, idir) = dpc/valCent;
            }
        }

      EBCellFAB gradMagFab(ebisBox, b, 1);
      BaseFab<Real>& regGradMagFab = gradMagFab.getSingleValuedFAB();
      FORT_MAGNITUDE(CHF_FRA1(regGradMagFab,0),
                     CHF_CONST_FRA(regGradFab),
                     CHF_BOX(b));

      //pointwise op so just have to iterate over multivalued cells
      IntVectSet ivsMulti = ebisBox.getMultiCells(b);
      for(VoFIterator vofit(ivsMulti, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real mag = 0.0;
          for(int idir = 0; idir < SpaceDim; idir++)
            {
              Real graddir = gradFab(vof, idir);
              mag += graddir*graddir;
            }
          mag = sqrt(mag);
          gradMagFab(vof, 0) = mag;
        }

      // Tag where gradient exceeds threshold

      IntVectSet ivsTot(b);
      for(VoFIterator vofit(ivsTot, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const  IntVect& iv  = vof.gridIndex();

          if ((iv[0] > imin) && (gradMagFab(vof, 0) >= m_params.m_refineThresh))
            {
              localTags |= iv;
            }
        }

      bool tagAllIrregular = false;
      bool tagInflow = false;
      if(pp.contains("tag_all_irregular"))
        {
          pp.get("tag_all_irregular", tagAllIrregular);
        }
      if(pp.contains("tag_inflow"))
        {
          pp.get("tag_inflow", tagInflow);
        }
      if(tagAllIrregular || s_noEBCF)
        {
          pout() << "tagging all irregular cells" << endl;
          //refine all irregular cells
          //this is probably not ideal.
          IntVectSet irregIVS = ebisBox.getIrregIVS(b);
          localTags |= irregIVS;
          s_noEBCF = true;
        }
      if(tagInflow)
        {
          pout() << "tagging cells at inflow" << endl;
          Box lowBox = adjCellLo(m_eblg.getDomain().domainBox(), 0, 1);
          lowBox.shift(0, 1);
          localTags |= lowBox;
        }
    }

  localTags.grow(m_params.m_tagBufferSize);
  bool tagHiLo = false;
  pp.query("tag_hi_and_lo", tagHiLo);
  if(tagHiLo)
    {
      Box loBox = adjCellLo(m_problem_domain.domainBox(), 1, 1);
      Box hiBox = adjCellHi(m_problem_domain.domainBox(), 1, 1);
      loBox.shift(1,  1);
      hiBox.shift(1, -1);
      localTags |= loBox;
      localTags |= hiBox;
    }

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_eblg.getDomain();
  localTags &= localTagsBox;
  a_tags = localTags;
}
//---------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------
void
EBAMRCNS::
tagCellsInit(IntVectSet& a_tags)
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::tagCellsInit for level " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  tagCells(a_tags);
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
updateMomentumAndEnergy(const Vector<LevelData<EBCellFAB>*>& a_newVelocity,
                        const Vector<LevelData<EBCellFAB>*>& a_newTemperature,
                        int a_baseLevel, 
                        int a_maxLevel)
{
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  Interval velocityInterval(0, SpaceDim-1);
  Interval momentumInterval(CMOMX, CMOMX+SpaceDim-1);
  Interval tempInterval(0, 0);
  Interval energyInterval(CENG, CENG);
  for(int ilev = a_baseLevel; ilev <= a_maxLevel; ++ilev)
    {
      EBAMRCNS* cnsLevel = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
      const DisjointBoxLayout& dbl = cnsLevel->m_eblg.getDBL();
      LevelData<EBCellFAB>& state  = cnsLevel->m_stateNew;
      for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          const Box& region = dbl.get(dit());
          EBCellFAB& stateFAB = state[dit()];
          EBCellFAB& velFAB = (*a_newVelocity[ilev])[dit()];
          EBCellFAB& tempFAB = (*a_newTemperature[ilev])[dit()];

          // The momentum density is rho * v.
          stateFAB.copy(region, momentumInterval, region, velFAB, velocityInterval);
          for (int idir = 0; idir < SpaceDim; ++idir)
            stateFAB.mult(stateFAB, CRHO, CMOMX + idir, 1);

          // The total energy density is rho * E = 0.5*rho*v*v + rho*cv*T, where 
          // E is the total specific energy. We thus compute v*v in velFAB and 
          // accumulate everything into the energy density.

          // Add in the specific internal energy.
          tempFAB.mult(m_params.m_specHeatCv);
          stateFAB.copy(region, energyInterval, region, tempFAB, tempInterval);

          // Add in the specific kinetic energy.
          velFAB *= velFAB; // v**2 componentwise.
          for (int idir = 0; idir < SpaceDim; ++idir)
            stateFAB.plus(velFAB, idir, CENG, 1);

          // Multiply by density.
          stateFAB.mult(stateFAB, CRHO, CENG, 1);
        }
    }
}
//-------------------------------------------------------------------------
void
EBAMRCNS::
preRegrid(int a_base_level, const Vector<Vector<Box> >& a_new_grids)
{
  //(\rho v)^{pre} = (\rho I - \dt_b L^v)v^n
  //E^{pre} = (\rho C_v I - \dt_b L^T)T^n + \half \rho |v^{pre}|^2)

}
//-------------------------------------------------------------------------
void
EBAMRCNS::
postRegrid(int a_base_level)
{
  if(m_level == a_base_level) defineSolvers();
}
//-------------------------------------------------------------------------
void
EBAMRCNS::
regrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRCNS::regrid");
  //first save old data
  // save data for later copy
  //not using m_eblg.getEBISL() because it gets wiped later
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  EBISLayout ebislOld;
  ebisPtr->fillEBISLayout(ebislOld, m_eblg.getDBL(), m_eblg.getDomain().domainBox(), m_nGhost);
  EBCellFactory factoryOld(ebislOld);
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  LevelData<EBCellFAB> stateSaved(m_eblg.getDBL(), m_nComp, ivGhost, factoryOld);
  Interval interv(0,m_nComp-1);
  stateSaved.define(m_eblg.getDBL(), m_nComp, ivGhost, factoryOld);
  m_stateNew.copyTo(interv, stateSaved, interv);

  //create grids and ebis layouts
  //the ebisl has to know about the fact that we really have
  //four ghost cells.
  m_level_grids = a_new_grids;
  Vector<int> proc_map;

  EBLoadBalance(proc_map,a_new_grids, m_eblg.getDomain().domainBox());

  DisjointBoxLayout grids(a_new_grids, proc_map);

  m_eblg.define(grids, m_problem_domain, 6, ebisPtr);
  pout() << "regrid: number of grids now in level " << m_level << " = " << grids.size() << endl;
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS regrid for level " << m_level << endl;
    }



  // set up data structures
  levelSetup();

  // interpolate to coarser level
  if (m_hasCoarser)
    {
      EBAMRCNS* coarPtr = getCoarserLevel();
      m_ebFineInterp.interpolate(m_stateNew,
                                 coarPtr->m_stateNew,
                                 interv);
    }

  // copy from old state
  //  pout() << "new grids for level " << m_level << " = " << m_eblg.getDBL() << endl;
  stateSaved.copyTo(interv,m_stateNew, interv);
  m_stateNew.copyTo(interv,m_stateOld, interv);
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
initialGrid(const Vector<Box>& a_new_grids)
{
  CH_TIME("EBAMRCNS::initialGrid");
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS initialGrid for level " << m_level << endl;
    }

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  m_level_grids = a_new_grids;

  // load balance and create boxlayout
  Vector<int> proc_map;
  if(m_params.m_doDiffusion)
    {
      EBEllipticLoadBalance(proc_map,a_new_grids, m_problem_domain.domainBox());
    }
  else
    {
      LoadBalance(proc_map,a_new_grids);
    }
  if(m_params.m_verbosity >= 3)
    {
      pout() << " just loadbalanced " << m_level << endl;
    }

  DisjointBoxLayout grids(a_new_grids,proc_map);
  m_eblg.define(grids, m_problem_domain, m_nGhost, ebisPtr);
  pout() << "initial grid: number of grids for level " << m_level << " = " << grids.size() << endl;
  //  pout() << "initial grids for level " << m_level << " = " << m_eblg.getDBL() << endl;
//  if(m_params.m_verbosity >= 5)
//    {
//      pout() << "EBAMRCNS::initialgrid grids " << endl;
//      DisjointBoxLayout dbl = m_eblg.getDBL();
//      dumpDBL(&dbl);
//    }

  // set up data structures
  levelSetup();
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
initialData()
{

  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS initialData for level " << m_level << endl;
    }

  // If we've been initialized with an IBC adapter function, look to the IBC
  // object for initial conditions.
  if (dynamic_cast<EBSpaceTimeFunctionIBCAdaptor*>(&(*m_ICs)))
    {
      const EBPhysIBC* const ebphysIBCPtr =
        m_ebPatchGodunov->getEBPhysIBC();

      ebphysIBCPtr->initialize(m_stateNew, m_eblg.getEBISL());
      ebphysIBCPtr->initialize(m_stateOld, m_eblg.getEBISL());
    }

  // Otherwise, just initialize both new and old states with our space-time
  // function.
  else
    {
      m_ICs->setDomain(m_problem_domain);
      m_ICs->evaluate(m_stateNew, m_eblg.getEBISL(), m_dx);
      m_ICs->evaluate(m_stateOld, m_eblg.getEBISL(), m_dx);
    }
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
postInitialGrid(const bool a_restart)
{
  if(m_level == 0) defineSolvers();
}
void
EBAMRCNS::
postInitialize()
{
  // Average the finer grid levels to be consistent with this one.
  if (m_hasFiner)
    {
      Interval interv(0, m_nComp-1);
      EBAMRCNS* finePtr = getFinerLevel();
      finePtr->m_ebCoarseAverage.average(m_stateNew,
                                         finePtr->m_stateNew,
                                         interv);
    }

  // Here we record the totals for the conserved quantities on grid level 0.
  if(!m_hasCoarser)
    {
      // Sum the quantities.
      int densityIndex = m_ebPatchGodunov->densityIndex();
      sumConserved(m_originalMass, densityIndex);
      Interval velInt = m_ebPatchGodunov->velocityInterval();
      for (int iv = velInt.begin(); iv < velInt.end(); ++iv)
        sumConserved(m_originalMomentum[iv-velInt.begin()], iv);
      sumConserved(m_originalEnergy, CENG); // Index not otherwise accessible.
    }
  //  if(m_level==0) defineSolvers();
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
setConservative(bool a_conservative)
{
  m_isConservative = a_conservative;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
syncWithFineLevel()
{
  CH_TIME("EBAMRCNS::syncWithFine");
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS syncWithFineLevel for level " << m_level << endl;
    }
  //stuff that needs to be setup from the finer
  //level.  A bunch of objects depend upon the layouts
  //from both levels and the finer level changes more
  //often from regrid so this needs to be called from the finer
  //level
  CH_assert(m_hasFiner);
  if(m_hasFiner)
    {
      EBAMRCNS* finePtr = getFinerLevel();
      int nRefFine = refRatio();
      const EBLevelGrid& finer_eblg = finePtr->m_eblg;
      const DisjointBoxLayout& finer_dbl = finer_eblg.getDBL();
      const EBISLayout& finer_ebisl      = finer_eblg.getEBISL();
      // maintain flux registers
      m_divFFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                m_nComp, Chombo_EBIS::instance(), s_noEBCF);

      m_veloFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                SpaceDim, Chombo_EBIS::instance(), s_noEBCF);

      m_tempFluxRegister.define(finer_dbl,
                                m_eblg.getDBL(),
                                finer_ebisl,
                                m_eblg.getEBISL(),
                                m_eblg.getDomain().domainBox(),
                                nRefFine,
                                1, Chombo_EBIS::instance(), s_noEBCF);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      if (!s_noEBCF)
        {
          m_ebCoarToFineRedist.define(finer_eblg, m_eblg, nRefFine , m_nComp, 1);
          //define coarse to coarse redistribution object
          m_ebCoarToCoarRedist.define(finer_eblg, m_eblg, nRefFine , m_nComp, 1);
        }

      //set all the registers to zero
      if(!s_noEBCF)
        {
          m_ebCoarToFineRedist.setToZero();
          m_ebCoarToCoarRedist.setToZero();
        }
      m_divFFluxRegister.setToZero();
      m_veloFluxRegister.setToZero();
      m_tempFluxRegister.setToZero();
    }

}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
Real
EBAMRCNS::
computeDt()
{
  return m_dtNew;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
Real
EBAMRCNS::
computeInitialDt()
{
  Real maxwavespeed =  m_ebLevelGodunov.getMaxWaveSpeed(m_stateNew);
  CH_assert(Abs(m_dx[0]-m_dx[1])<1.e-9);
  Real newDT = m_params.m_initialDtMultiplier * m_dx[0] /maxwavespeed;
  pout() << "max wave speed = " << maxwavespeed << ", dx = " << m_dx[0] <<", new_dt = "  << newDT << endl;

  return newDT;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
sumConserved(Real& a_sumcons,
             const int& a_ivar) const
{
  CH_TIME("EBAMRCNS::sumConserved");
  Real sumgrid = 0;
  for(DataIterator dit= m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& state = m_stateNew[dit()];
      Box thisBox = m_eblg.getDBL().get(dit());
      IntVectSet uberIVS(thisBox);
      const EBISBox& ebisBox = m_eblg.getEBISL()[dit()];
      for(VoFIterator vofit(uberIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real consVal = state(vof, a_ivar);
          Real volFrac = ebisBox.volFrac(vof);
          Real volume = volFrac;
          sumgrid += consVal*volume;
        }
    }

  Vector<Real> all_sum;
  gather(all_sum,sumgrid,uniqueProc(SerialTask::compute));
  Real sumallgrid = 0.;
  if (procID() == uniqueProc(SerialTask::compute))
    {
      for (int i = 0; i < all_sum.size(); ++i)
        {
          sumallgrid += all_sum[i];
        }
    }
  broadcast(sumallgrid,uniqueProc(SerialTask::compute));
  a_sumcons = sumallgrid;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
EBAMRCNS*
EBAMRCNS::
getCoarserLevel() const
{
  EBAMRCNS* retval = NULL;
  if(m_coarser_level_ptr != NULL)
    {
      retval = dynamic_cast <EBAMRCNS*> (m_coarser_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getCoarserLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
EBAMRCNS*
EBAMRCNS::
getFinerLevel() const
{
  EBAMRCNS* retval = NULL;
  if(m_finer_level_ptr != NULL)
    {
      retval = dynamic_cast<EBAMRCNS*>(m_finer_level_ptr);

      if (retval == NULL)
        {
          pout() << "EBAMRG::getFinerLevel: dynamic cast failed"
                 << endl;
          MayDay::Error();
        }
    }
  return retval;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
levelSetup()
{
  CH_TIME("levelSetup");
  if(m_params.m_verbosity >= 3)
    {
      pout() << " in EBAMRCNS levelSetup for level " << m_level << endl;
    }

  EBCellFactory factoryNew(m_eblg.getEBISL());
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_redisRHS.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  EBLevelDataOps::setVal(m_redisRHS, 0.0);


  m_sets.define(m_eblg.getDBL());
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_sets[dit()] = m_eblg.getEBISL()[dit()].getIrregIVS(m_eblg.getDBL().get(dit()));
    }
  EBCellFactory       cellFact(m_eblg.getEBISL());
  EBFluxFactory       fluxFact(m_eblg.getEBISL());
  BaseIVFactory<Real> bivfFact(m_eblg.getEBISL(), m_sets);

  // Allocate the A coefficients for our viscous and conduction operators.
  m_acoVelo    = RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblg.getDBL(), 1, 4*IntVect::Unit, cellFact));
  m_acoTemp    = RefCountedPtr< LevelData<EBCellFAB> >       (new LevelData<EBCellFAB>       (m_eblg.getDBL(), 1, 4*IntVect::Unit, cellFact));

  // Allocate the viscous tensor coefficients on regular and irregular cells.
  m_eta         = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, 4*IntVect::Unit, fluxFact));
  m_etaIrreg    = RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, 4*IntVect::Unit, bivfFact));
  m_lambda      = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, 4*IntVect::Unit, fluxFact));
  m_lambdaIrreg = RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, 4*IntVect::Unit, bivfFact));

  // Allocate the thermal conductivity on regular and irregular cells.
  m_bcoTemp     = RefCountedPtr< LevelData<EBFluxFAB> >       (new LevelData<EBFluxFAB>       (m_eblg.getDBL(), 1, 4*IntVect::Unit, fluxFact));
  m_bcoTempIrreg= RefCountedPtr< LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >(m_eblg.getDBL(), 1, 4*IntVect::Unit, bivfFact));

  EBAMRCNS* coarPtr = getCoarserLevel();
  EBAMRCNS* finePtr = getFinerLevel();

  m_hasCoarser = (coarPtr != NULL);
  m_hasFiner   = (finePtr != NULL);

  //define redistribution object for this level
  //for now set to volume weighting
  m_ebLevelRedist.define(m_eblg.getDBL(),
                         m_eblg.getEBISL(),
                         m_eblg.getDomain(),
                         m_nComp);

  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      const EBLevelGrid& coEBLG = coarPtr->m_eblg;

      m_ebCoarseAverage.define(m_eblg.getDBL(),
                               coEBLG.getDBL(),
                               m_eblg.getEBISL(),
                               coEBLG.getEBISL(),
                               coEBLG.getDomain().domainBox(),
                               nRefCrse,
                               m_nComp, Chombo_EBIS::instance());
      m_ebFineInterp.define(m_eblg.getDBL(),
                            coEBLG.getDBL(),
                            m_eblg.getEBISL(),
                            coEBLG.getEBISL(),
                            coEBLG.getDomain().domainBox(),
                            nRefCrse,
                            m_nComp);

      // maintain levelgodunov
      m_ebLevelGodunov.define(m_eblg.getDBL(),
                              coEBLG.getDBL(),
                              m_eblg.getEBISL(),
                              coEBLG.getEBISL(),
                              m_eblg.getDomain(),
                              nRefCrse,
                              m_dx,
                              m_params,
                              m_ebPatchGodunov,
                              m_hasCoarser,
                              m_hasFiner);

      //define fine to coarse redistribution object
      //for now set to volume weighting
      if(!s_noEBCF)
        {
          CH_TIME("fineToCoar_defs");
          m_ebFineToCoarRedist.define(m_eblg, coEBLG,
                                      nRefCrse, m_nComp, 
                                      m_params.m_redistRad);
                                   
          m_ebFineToCoarRedist.setToZero();
        }

      int nvarQuad = 1; //temperature

      m_quadCFI = RefCountedPtr<EBQuadCFInterp>
        (new EBQuadCFInterp(m_eblg.getDBL(),
                            coEBLG.getDBL(),
                            m_eblg.getEBISL(),
                            coEBLG.getEBISL(),
                            coEBLG.getDomain(),
                            nRefCrse, nvarQuad,
                            (*m_eblg.getCFIVS()),
                            Chombo_EBIS::instance()));

      coarPtr->syncWithFineLevel();
    }
  else
    {
      m_quadCFI = RefCountedPtr<EBQuadCFInterp>(new EBQuadCFInterp());
      m_ebLevelGodunov.define(m_eblg.getDBL(),
                              DisjointBoxLayout(),
                              m_eblg.getEBISL(),
                              EBISLayout(),
                              m_eblg.getDomain(),
                              m_ref_ratio,
                              m_dx,
                              m_params,
                              m_ebPatchGodunov,
                              m_hasCoarser,
                              m_hasFiner);
    }

  m_sets.define(m_eblg.getDBL());
  for(DataIterator dit = m_eblg.getDBL().dataIterator();
      dit.ok(); ++dit)
    {
      Box thisBox = m_eblg.getDBL().get(dit());
      m_sets[dit()] = m_eblg.getEBISL()[dit()].getIrregIVS(thisBox);
    }
  BaseIVFactory<Real> factory(m_eblg.getEBISL(), m_sets);
  //the ghost cells on the mass redistribution array
  //are tied to the redistribution radius
  m_massDiff.define(m_eblg.getDBL(), m_nComp, ivGhost, factory);
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      m_massDiff[dit()].setVal(0.0);
    }
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
LevelData<EBCellFAB>&
EBAMRCNS::
getNewState()
{
  return m_stateNew;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
LevelData<EBCellFAB>&
EBAMRCNS::
getOldState()
{
  return m_stateOld;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
Vector<LevelData<EBCellFAB>*>
EBAMRCNS::
getNewStates()
{
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  int numLevels = hierarchy.size();
  Vector<LevelData<EBCellFAB>*> states(numLevels);
  for (int ilev = 0; ilev < numLevels; ++ilev)
    {
      EBAMRCNS* solver = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
      CH_assert(solver != 0);
      states[ilev] = &(solver->m_stateNew);
    }
  return states;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
Vector<LevelData<EBCellFAB>*>
EBAMRCNS::
getOldStates()
{
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  int numLevels = hierarchy.size();
  Vector<LevelData<EBCellFAB>*> states(numLevels);
  for (int ilev = 0; ilev < numLevels; ++ilev)
    {
      EBAMRCNS* solver = dynamic_cast<EBAMRCNS*>(hierarchy[ilev]);
      CH_assert(solver != 0);
      states[ilev] = &(solver->m_stateOld);
    }
  return states;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
Real
EBAMRCNS::
getDt() const
{
  return m_dt;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
EBISLayout
EBAMRCNS::
getEBISLayout() const
{
  EBISLayout ebisl = m_eblg.getEBISL();
  return ebisl;
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
fillConsAndPrim(LevelData<EBCellFAB>& a_data, LevelData<EBCellFAB>& a_state, int logflag)
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::fillConsAndPrim" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;

  EBCellFactory ebcellfact(m_eblg.getEBISL());
  a_data.define(m_eblg.getDBL(), consAndPrim, IntVect::Zero, ebcellfact);

  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
      const Box& grid = m_eblg.getDBL().get(dit());
      const EBCellFAB& consfab = a_state[dit()];
      EBCellFAB primfab(ebisbox, grid, nPrim);
      //cfivs, time and timestep fake and not used here
      Real faket = 1.0;
      IntVectSet emptyivs;
      m_ebPatchGodunov->setValidBox(grid, ebisbox, emptyivs, faket, faket);
      m_ebPatchGodunov->consToPrim(primfab, consfab, grid,logflag);

      EBCellFAB& outputfab = a_data[dit()];

      Interval consIntervSrc(0, nCons-1);
      Interval consIntervDst(0, nCons-1);
      Interval primIntervSrc(0, nPrim-1);
      Interval primIntervDst(nCons, consAndPrim-1);

      // copy regular data
      outputfab.copy(grid, consIntervDst,  grid, consfab, consIntervSrc);
      outputfab.copy(grid, primIntervDst,  grid, primfab, primIntervSrc);
      Real coveredVal = -1.2345678e-9;
      for(int ivar = 0; ivar < consAndPrim; ivar++)
        {
          outputfab.setInvalidData(coveredVal, ivar);
        }

    }//end loop over grids
}
//-------------------------------------------------------------------------

#ifdef CH_USE_HDF5
//-------------------------------------------------------------------------
void
EBAMRCNS::
writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  //so i will eliminate the middleman
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
readCheckpointHeader(HDF5Handle& a_handle)
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::readCheckpointHeader" << endl;
    }

  //stuff in non-eb checkpoint header is already
  //set in the define function, such as
  // Setup the number of components
  // Setup the component names
  // So i will eliminate the middleman.
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writeCheckpointLevel" << endl;
    }

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_real["dt"]              = m_dt;


  // Write the header for this level
  header.writeToFile(a_handle);

  // Write the data for this level
  write(a_handle,m_eblg.getDBL());
  write(a_handle,m_stateOld,"dataOld");
  write(a_handle,m_stateNew,"dataNew");
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
readCheckpointLevel(HDF5Handle& a_handle)
{
  CH_assert(m_isDefined);
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  //reasons for deviations from non-eb stuff
  // the tag buffer size is set by the factory.
  // dx is set in the define function
  // the problem domain is set in the define function

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  // Get the grids
  Vector<Box> vboxGrids;
  const int gridStatus = read(a_handle, vboxGrids);
  if (gridStatus != 0)
    {
      MayDay::Error("readCheckpointLevel: file has no grids");
    }

  Vector<int> proc_map;
  EBLoadBalance(proc_map,vboxGrids, m_problem_domain);

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  DisjointBoxLayout grids(vboxGrids,proc_map);
  m_eblg.define(grids, m_problem_domain, m_nGhost, ebisPtr);

  //this keeps the order of the AMRLevel m_level_grids
  //consistent with m_eblg.getDBL()
  LayoutIterator lit = m_eblg.getDBL().layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      Box b = m_eblg.getDBL().get(lit());
      m_level_grids.push_back(b);
    }

  EBCellFactory factoryNew(m_eblg.getEBISL());
  //m_nghost is set in define function
  IntVect ivGhost = m_nGhost*IntVect::Unit;
  m_stateNew.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_stateOld.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  m_redisRHS.define(m_eblg.getDBL(),m_nComp, ivGhost, factoryNew);
  EBLevelDataOps::setVal(m_redisRHS, 0.0);

  // Set up data structures
  levelSetup();
  //  Interval vars(0, m_nComp-1);
  //the false says to not redefine the data
  int dataStatusNew = read<EBCellFAB>(a_handle,
                                      m_stateNew,
                                      "dataNew",
                                      m_eblg.getDBL(),
                                      Interval(),
                                      false);

  int dataStatusOld = read<EBCellFAB>(a_handle,
                                      m_stateOld,
                                      "dataOld",
                                      m_eblg.getDBL(),
                                      Interval(),
                                      false);

  if ((dataStatusNew != 0) || (dataStatusOld != 0))
    {
      MayDay::Error("file does not contain state data");
    }

}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
writePlotHeader(HDF5Handle& a_handle) const
{
  if(s_NewPlotFile == 0)
    {
      writePlotHeaderOld(a_handle);
      return;
    }
  Vector<int> refRatios;
  const EBAMRCNS* current = this;
  int nlevs = 0;
  while(current != NULL){
    refRatios.push_back(current->refRatio());
    nlevs++;
    current = dynamic_cast<const EBAMRCNS*>(current->m_finer_level_ptr);
  }

  headerEBFile(a_handle, nlevs, refRatios,
               m_eblg.getDomain().domainBox(), m_dx, m_stateNew.ghostVect());

  writeCellCenteredNames(a_handle, m_stateNames);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  EBPatchGodunov& patchgod = 
    const_cast<EBPatchGodunov&>(dynamic_cast<const EBPatchGodunov&>(*m_ebPatchGodunov));
  patchgod.expressions(expressions);
  expressions.writeToFile(a_handle);

}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void
EBAMRCNS::
writePlotHeaderOld(HDF5Handle& a_handle) const
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writePlotHeader" << endl;
    }

  HDF5HeaderData header;
  // Setup the number of components
  //have to add in a lot of geometric crap.
  // 3 norms + 6 area fracs + 1 distance + 1 volFrac
  // =  11 extra components
  //forces 3d
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;

  int indexVolFrac = consAndPrim;
  int indexAreaFrac = indexVolFrac+1;
  int indexNormal = indexAreaFrac+ 2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int nCompTotal = indexDist+1;
  CH_assert(nCompTotal == consAndPrim + 3*SpaceDim+2);

  Vector<string> names(nCompTotal);

  for (int i = 0; i < nCons; i++)
    {
      names[i] = m_stateNames[i];
    }
  for (int i = 0; i < nPrim; i++)
    {
      names[nCons + i] = m_primNames[i];
    }

  string volFracName("fraction-0");
  Vector<string> normName(3);
  Vector<string> areaName(6);
  string distName("distance-0");

  normName[0] = "xnormal-0";
  normName[1] = "ynormal-0";
  normName[2] = "znormal-0";

  areaName[0] = "xAreafractionLo-0";
  areaName[1] = "xAreafractionHi-0";
  areaName[2] = "yAreafractionLo-0";
  areaName[3] = "yAreafractionHi-0";
  areaName[4] = "zAreafractionLo-0";
  areaName[5] = "zAreafractionHi-0";

  names[indexVolFrac] = volFracName;

  for (int i=0; i < 2*SpaceDim; i++)
    {
      names[indexAreaFrac+i] = areaName[i];
    }

  for (int i=0; i < SpaceDim; i++)
    {
      names[indexNormal+i] = normName[i];
    }

  names[indexDist] = distName;

  //now output this into the hdf5 handle
  header.m_int["num_components"] = nCompTotal;
  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < nCompTotal; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = names[comp];
    }

  // Write the header to the file
  header.writeToFile(a_handle);

  if (m_params.m_verbosity >= 4)
    {
      pout() << header << endl;
    }
}
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
void EBAMRCNS::writePlotLevel(HDF5Handle& a_handle) const
{
  if(s_NewPlotFile == 0)
    {
      writePlotLevelOld(a_handle);
      return;
    }
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writePlotLevel " << m_level<< endl;
    }
  writeCellCentered(a_handle, m_level, &m_stateNew);
}

//-------------------------------------------------------------------------
void EBAMRCNS::writePlotLevelOld(HDF5Handle& a_handle) const
{

  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRCNS::writePlotLevel" << endl;
    }
  //fill output data including a bunch of geometric crud
  int nCons = m_ebPatchGodunov->numConserved();
  int nPrim = m_ebPatchGodunov->numPrimitives() ;
  int consAndPrim = nCons + nPrim;
  int indexVolFrac = consAndPrim;
  int indexAreaFrac = indexVolFrac+1;
  int indexNormal = indexAreaFrac+ 2*SpaceDim;
  int indexDist = indexNormal+SpaceDim;
  int nCompTotal = indexDist+1;
  CH_assert(nCompTotal == consAndPrim + 3*SpaceDim+2);

  Real coveredVal = -1.2345678e-9;
  Vector<Real> coveredValuesCons(nCons, coveredVal);
  Vector<Real> coveredValuesPrim(nPrim, coveredVal);

  Vector<Real> coveredValues;
  coveredValues.append(coveredValuesCons);
  coveredValues.append(coveredValuesPrim);

  LevelData<FArrayBox> fabData(m_eblg.getDBL(), nCompTotal, IntVect::Zero);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
      const Box& grid = m_eblg.getDBL().get(dit());
      const EBCellFAB& consfab = m_stateNew[dit()];
      EBCellFAB primfab(ebisbox, grid, nPrim);
      //cfivs, time and timestep fake and not used here
      Real faket = 1.0;
      IntVectSet emptyivs;
      ParmParse pp;
      int logflag = 0;
      if(pp.contains("logflag"))
        {
          pp.get("logflag", logflag);
        }
      EBPatchGodunov& patchgod = (EBPatchGodunov&)(*m_ebPatchGodunov);
      patchgod.setValidBox(grid, ebisbox, emptyivs, faket, faket);
      patchgod.consToPrim(primfab, consfab, grid, logflag);

      FArrayBox& currentFab = fabData[dit()];

      // copy regular data
      currentFab.copy(consfab.getSingleValuedFAB(),0,0,nCons);
      currentFab.copy(primfab.getSingleValuedFAB(),0,nCons,nPrim);

      // set default volume fraction
      currentFab.setVal(1.0,indexVolFrac);

      // set default area fractions
      for (int i=0; i < 2*SpaceDim; i++)
        {
          currentFab.setVal(1.0,indexAreaFrac+i);
        }

      // set default normal
      for (int i=0; i < SpaceDim; i++)
        {
          currentFab.setVal(0.0,indexNormal+i);
        }

      // set default distance of EB from corner
      currentFab.setVal(0.0,indexDist);

      // set special values
      // iterate through the current grid
      // NOTE:  this is probably an inefficient way to do this
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          // set special values for covered cells
          if (ebisbox.isCovered(iv))
            {
              for (int icomp = 0; icomp < consAndPrim; icomp++)
                {
                  Real cval = coveredValues[icomp];

                  currentFab(iv,icomp) = cval;
                }
              // volume fraction is zero
              currentFab(iv,indexVolFrac) = 0.0;

              // area fractions are zero
              for (int i=0; i < 2*SpaceDim; i++)
                {
                  currentFab(iv,indexAreaFrac+i) = 0.0;
                }
            }

          // set special values for irregular cells
          if (ebisbox.isIrregular(iv))
            {
              Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
              Real volFrac = ebisbox.volFrac(vofs[0]);
              RealVect normal = ebisbox.normal(vofs[0]);

              // set volume fraction
              currentFab(iv,indexVolFrac) = volFrac;

              // set area fractions--use only the first face you find
              for (int i=0; i < SpaceDim; i++)
                {
                  Vector<FaceIndex> faces;

                  faces = ebisbox.getFaces(vofs[0],i,Side::Lo);
                  if (faces.size() == 0)
                    {
                      currentFab(iv,indexAreaFrac+2*i) = 0.0;
                    }
                  else
                    {
                      currentFab(iv,indexAreaFrac+2*i) =
                        ebisbox.areaFrac(faces[0]);
                    }

                  faces = ebisbox.getFaces(vofs[0],i,Side::Hi);
                  if (faces.size() == 0)
                    {
                      currentFab(iv,indexAreaFrac+2*i+1) = 0.0;
                    }
                  else
                    {
                      currentFab(iv,indexAreaFrac+2*i+1) =
                        ebisbox.areaFrac(faces[0]);
                    }
                }

              // set normal
              for (int i=0; i < SpaceDim; i++)
                {
                  currentFab(iv,indexNormal+i) = normal[i];
                }

              // set distance unless the length of the normal is zero
              Real length = PolyGeom::dot(normal,normal);

              if (length > 0)
                {
                  Real dist = PolyGeom::computeAlpha(volFrac,normal)*m_dx[0];
                  currentFab(iv,indexDist) = -dist;
                }
            } //end if(isIrregular)
        }//end loop over cells
    }//end loop over grids

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx[0];
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_eblg.getDomain().domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (m_params.m_verbosity >= 4)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  write(a_handle,fabData.boxLayout());
  write(a_handle,fabData,"data");
}
//-------------------------------------------------------------------------

#endif
#include "NamespaceFooter.H"
