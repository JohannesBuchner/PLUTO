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
#include <cstdio>
#include <iostream>

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "EBBackwardEuler.H"
#include "AllRegularService.H"
#include "EBLevelRedist.H"
#include "DirichletPoissonDomainBC.H"
#include "DirichletPoissonEBBC.H"
#include "EBAMRPoissonOpFactory.H"
#include "EBLevelTGA.H"
#include "EBAMRDataOps.H"
#include "RedistStencil.H"
#include "BiCGStabSolver.H"
#include "EBAMRIO.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "SphereIF.H"
#include "EBAMRIO.H"
#include "EBAMRPoissonOp.H"

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "UsingNamespace.H"

/************/
RealVect s_center = 0.5*RealVect::Unit;
Real     s_radius = 0.45;
int      s_nCellsFine = 32;
int      s_nStepsFine = 16;
Real     s_diffusivity = 0.1;
IntVect  s_giv = 4*IntVect::Unit;
int s_precond = 4;
int s_reltype = 1;
BiCGStabSolver<LevelData<EBCellFAB> > s_bottomsolv;
/************/
Real
getValue(const VolIndex& a_vof,
         const Real&     a_dx)
{
  RealVect location = EBArith::getVofLocation(a_vof, a_dx*RealVect::Unit, RealVect::Zero);
  location -= s_center;
  Real sumSquare;
  PolyGeom::unifyVector(location, sumSquare);
  Real normrad = sqrt(sumSquare)/s_radius;

  Real value = (normrad)*(1.0-normrad);
  return value;
}
/************/
void
fillData(LevelData<EBCellFAB>& a_data,
         const DisjointBoxLayout a_dbl,
         const EBISLayout&       a_ebisl,
         const Box& a_domain,
         const Real& a_dx)
{
  for (DataIterator dit = a_dbl.dataIterator(); dit.ok(); ++dit)
    {
      IntVectSet ivs(a_dbl.get(dit()));
      for (VoFIterator vofit(ivs, a_ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          for (int ivar = 0; ivar < a_data.nComp(); ivar++)
            {
              Real value = getValue(vofit(), a_dx);
              a_data[dit()](vofit(), ivar) = value;
            }
        }
    }
}
/************/


void
defineHeatSolver(RefCountedPtr< BaseLevelHeatSolver< LevelData<EBCellFAB>,  EBFluxFAB, EBFluxRegister> >& a_heat,
                 RefCountedPtr<AMRMultiGrid<LevelData <EBCellFAB> > >                                   & a_amrmg,
                 const DisjointBoxLayout                                                                & a_grids,
                 const EBISLayout                                                                       & a_ebisl,
                 const Box                                                                              & a_domain,
                 const Real                                                                             & a_dx,
                 int a_whichIntegrator)
{
  ProblemDomain probdom(a_domain);
  Vector<DisjointBoxLayout> grids(1, a_grids);
  Vector<EBISLayout>        ebisl(1, a_ebisl);
  Vector<EBLevelGrid>       eblgs(1, EBLevelGrid(a_grids, a_ebisl, probdom));
  Vector<int>              refrat(1,2);
  Vector<RefCountedPtr<EBQuadCFInterp> >    quadcfi(1);
  Real alpha = 1.0;
  Real time = 0;
  RefCountedPtr<DirichletPoissonDomainBCFactory>  dobc(new DirichletPoissonDomainBCFactory());;
  RefCountedPtr<DirichletPoissonEBBCFactory    >  ebbc(new DirichletPoissonEBBCFactory());
  dobc->setValue(0.);
  ebbc->setValue(0.);
  RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >
    opfact(new EBAMRPoissonOpFactory(eblgs, refrat, quadcfi, a_dx*RealVect::Unit, RealVect::Zero, s_precond,
                                     s_reltype, dobc, ebbc, alpha, s_diffusivity, time, s_giv, s_giv));


  a_amrmg = RefCountedPtr<AMRMultiGrid<LevelData <EBCellFAB> > >
    (new AMRMultiGrid<LevelData <EBCellFAB> > () );
  a_amrmg->m_verbosity = 2;
  s_bottomsolv.m_verbosity = 0;
  a_amrmg->define(probdom, *opfact, &s_bottomsolv, 1);
  a_amrmg->m_verbosity = 2;
  if (a_whichIntegrator == 0)
    {
      pout() << "defining backward euler solver" << endl;
      EBLevelBackwardEuler* heatptr =   new EBLevelBackwardEuler(grids, refrat, probdom, opfact, a_amrmg);
      heatptr->setEBLG(eblgs);
      a_heat = RefCountedPtr< BaseLevelHeatSolver< LevelData<EBCellFAB>,  EBFluxFAB, EBFluxRegister> >(heatptr);
    }
  else if (a_whichIntegrator == 1)
    {
      pout() << "defining crank nicolson solver"  << endl;
      EBLevelCrankNicolson* heatptr =   new EBLevelCrankNicolson(grids, refrat, probdom, opfact, a_amrmg);
      heatptr->setEBLG(eblgs);

      a_heat = RefCountedPtr< BaseLevelHeatSolver< LevelData<EBCellFAB>,  EBFluxFAB, EBFluxRegister> >(heatptr);
    }
  else if (a_whichIntegrator == 2)
    {
      pout() << "defining TGA solver" << endl;
      EBLevelTGA* heatptr =   new EBLevelTGA(grids, refrat, probdom, opfact, a_amrmg);
      heatptr->setEBLG(eblgs);

      a_heat = RefCountedPtr< BaseLevelHeatSolver< LevelData<EBCellFAB>,  EBFluxFAB, EBFluxRegister> >(heatptr);
    }
  else
    {
      pout() << "Derp.  I do not know about integrator type " << a_whichIntegrator << endl;
      MayDay::Error();
    }
}
/************/
void
generateData(LevelData<EBCellFAB>           & a_datum,
             DisjointBoxLayout              & a_grids,
             EBISLayout                     & a_ebisl,
             const Box                      & a_domain,
             const Real                     & a_dx,
             const int                      & a_refToFinestCalc,
             const int                      & a_whichIntegrator)
{

  pout() << "filling data" << endl;
  EBCellFactory fact(a_ebisl);
  LevelData<EBCellFAB> phiOld(a_grids, 1, s_giv, fact);
  LevelData<EBCellFAB> phiNew(a_grids, 1, s_giv, fact);
  LevelData<EBCellFAB> source(a_grids, 1, s_giv, fact);

  fillData(phiOld, a_grids, a_ebisl, a_domain, a_dx);
  fillData(phiNew, a_grids, a_ebisl, a_domain, a_dx);
  fillData(source, a_grids, a_ebisl, a_domain, a_dx);


  RefCountedPtr< BaseLevelHeatSolver< LevelData<EBCellFAB>,  EBFluxFAB, EBFluxRegister> >  heatSolver;
  RefCountedPtr<AMRMultiGrid<LevelData <EBCellFAB> > >  amrmg;
  defineHeatSolver(heatSolver, amrmg, a_grids, a_ebisl, a_domain, a_dx, a_whichIntegrator);

  Real dt = a_dx;
  int nstop = s_nStepsFine;
  nstop /= a_refToFinestCalc;

  Real time = 0;
  LevelData<EBFluxFAB> dummyflux;

  amrmg->init(Vector<LevelData<EBCellFAB>*>(1,&phiOld),Vector<LevelData<EBCellFAB>*>(1, &source), 0, 0);
  for (int istep = 0; istep < nstop; istep++)
    {
      heatSolver->updateSoln(phiNew, phiOld, source,
                             NULL, NULL, NULL, NULL,
                             time, time, time, dt, 0);

      time += dt;
      EBLevelDataOps::assign(phiOld,phiNew);
      pout() << "step = " << istep << ", time = "  << time << endl;
    }

  a_datum.define(a_grids, 1, s_giv, fact);
  Interval interv(0, phiNew.nComp()-1);
  phiNew.copyTo(interv, a_datum, interv);
}
/*****/
void
getLayouts(DisjointBoxLayout& a_grids,
           EBISLayout       & a_ebisl,
           const Box&         a_domain)
{
  Vector<Box> boxes(1, a_domain);
  Vector<int> procs(1, 0);
  a_grids= DisjointBoxLayout(boxes, procs);
  Chombo_EBIS::instance()->fillEBISLayout(a_ebisl, a_grids, a_domain, 4);
}
/***/
void
makeSphereGeometry(Box& a_domainFine, Real& a_dxFine)
{
  a_domainFine = Box(IntVect::Zero, (s_nCellsFine-1)*IntVect::Unit);
  a_dxFine = 1.0/s_nCellsFine;

  RealVect origin = RealVect::Zero;
  pout() << "sphere geometry with fluid inside" << endl;
  SphereIF sphereIF(s_radius, s_center, true);
  GeometryShop workshop(sphereIF, 0,a_dxFine*RealVect::Unit);
  //this generates the new EBIS
  int biggridsize = 2048;
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domainFine, RealVect::Zero, a_dxFine, workshop, biggridsize);
}
/*****/
void
compareError(LevelData<EBCellFAB>& a_errorMedi,
             LevelData<EBCellFAB>& a_errorCoar,
             DisjointBoxLayout   & a_gridsMedi,
             DisjointBoxLayout   & a_gridsCoar,
             EBISLayout          & a_ebislMedi,
             EBISLayout          & a_ebislCoar,
             Box                 & a_domaiCoar,
             int& a_testverbosity)
{

  Vector<LevelData<EBCellFAB>* > errorMedi(1, &a_errorMedi);
  Vector<LevelData<EBCellFAB>* > errorCoar(1, &a_errorCoar);
  Vector<DisjointBoxLayout>      gridsMedi(1, a_gridsMedi);
  Vector<EBISLayout>             ebislMedi(1, a_ebislMedi);
  Vector<DisjointBoxLayout>      gridsCoar(1, a_gridsCoar);
  Vector<EBISLayout>             ebislCoar(1, a_ebislCoar);
  Vector<int> refRat(1, 2);
  Vector<Real> order;
  EBArith::compareError(order,
                        errorMedi, errorCoar,
                        gridsMedi, gridsCoar,
                        ebislMedi, ebislCoar,
                        refRat,  a_domaiCoar,
                        a_testverbosity);
}
void
solutionErrorTest(int testverbosity, int fileout)
{

  if (testverbosity >= 1)
    pout() << "generating geometry" << endl;
  Box domaiFine;
  Real dxFine;
  makeSphereGeometry(domaiFine, dxFine);

  Box domaiMedi = coarsen(domaiFine, 2);
  Box domaiCoar = coarsen(domaiMedi, 2);

  if (testverbosity >= 1)
    pout() << "getting grids" << endl;

  DisjointBoxLayout gridsFine, gridsMedi, gridsCoar;
  EBISLayout        ebislFine, ebislMedi, ebislCoar;

  getLayouts(gridsFine, ebislFine, domaiFine);
  getLayouts(gridsMedi, ebislMedi, domaiMedi);
  getLayouts(gridsCoar, ebislCoar, domaiCoar);

  int refFine = 1;  int refMedi = 2;   int refCoar = 4;

  Real dxMedi = refMedi*dxFine;
  Real dxCoar = refCoar*dxFine;

  //three integrators, three tests
  for (int integ = 0; integ < 3; integ++)
    {
      LevelData<EBCellFAB>       solutFine, solutMedi, solutCoar;
      LevelData<EBCellFAB>                  errorMedi, errorCoar;


      if (integ == 0)
        pout() << "testing backward euler" << endl;
      else if (integ == 1)
        pout() << "testing crank nicolson" << endl;
      else
        pout() << "testing TGA" << endl;

      if (testverbosity >= 1) pout() << "generating fine solution" << endl;
      generateData(solutFine, gridsFine, ebislFine, domaiFine, dxFine, refFine, integ);

      if (testverbosity >= 1) pout() << "generating medi solution" << endl;
      generateData(solutMedi, gridsMedi, ebislMedi, domaiMedi, dxMedi, refMedi, integ);

      if (testverbosity >= 1) pout() << "generating coar solution" << endl;
      generateData(solutCoar, gridsCoar, ebislCoar, domaiCoar, dxCoar, refCoar, integ);


      if (testverbosity >= 1) pout() << "generating medi error from medi and fine solutions" << endl;
      EBAMRDataOps::getErrorFromCoarseAndFine(errorMedi,
                                              solutMedi, gridsMedi, ebislMedi, domaiMedi,
                                              solutFine, gridsFine, ebislFine, domaiFine);


      if (testverbosity >= 1) pout() << "generating coar error from medi and coar solutions" << endl;
      EBAMRDataOps::getErrorFromCoarseAndFine(errorCoar,
                                              solutCoar, gridsCoar, ebislCoar, domaiCoar,
                                              solutMedi, gridsMedi, ebislMedi, domaiMedi);


      if (testverbosity >= 1)
        pout() << "comparing error " << endl;
      compareError(errorMedi, errorCoar,
                   gridsMedi, gridsCoar,
                   ebislMedi, ebislCoar,
                   domaiCoar, testverbosity);


#ifdef CH_USE_HDF5
      if (fileout == 1)
        {
          if (testverbosity > 1)  pout() << "Outputting error to file" << endl;
          char fileCoar[200];
          char fileMedi[200];
          if (integ == 0)
            {
              sprintf(fileCoar, "BackEulerErrorCoar.%dd.hdf5", SpaceDim);
              sprintf(fileMedi, "BackEulerErrorMedi.%dd.hdf5", SpaceDim);
            }
          else if (integ == 1)
            {
              sprintf(fileCoar, "CrankNicErrorCoar.%dd.hdf5", SpaceDim);
              sprintf(fileMedi, "CrankNicErrorMedi.%dd.hdf5", SpaceDim);
            }
          else if (integ == 2)
            {
              sprintf(fileCoar, "TGAErrorCoar.%dd.hdf5", SpaceDim);
              sprintf(fileMedi, "TGAErrorMedi.%dd.hdf5", SpaceDim);
            }
          writeEBLevelname(&errorCoar, fileCoar);
          writeEBLevelname(&errorMedi, fileMedi);
        }
#else
      pout() << "hdf5 not compiled in so I cannot output error to the file" << endl;
#endif

    }
}
/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  {
#endif

    int testverbosity = 1;
    int fileout = 1;
    solutionErrorTest(testverbosity, fileout);

#ifdef CH_MPI
    dumpmemoryatexit();
  }
  MPI_Finalize();
#endif
}

