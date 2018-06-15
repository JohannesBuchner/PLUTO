#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "EBExplosionIBCFactory.H"
#include "EBPatchPolytropicFactory.H"

#include "EBPatchPolytropic.H"

#include "EBAMRCNSFactory.H"
#include "EBAMRCNS.H"
#include "AMRLevel.H"
#include "AMR.H"
#include "BaseIVFactory.H"
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
#include "AllRegularService.H"
#include "EBLevelRedist.H"
#include "RedistStencil.H"
#include "SlabService.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBFABView.H"
#include "memtrack.H"
#include "EBAMRCNSParams.H"
#include "GodunovGeom.H"
#include "EBViscousTensorOpFactory.H"
#include "DirichletPoissonEBBC.H"
#include   "NeumannPoissonEBBC.H"
#include "DirichletPoissonDomainBC.H"
#include   "NeumannPoissonDomainBC.H"
#include "DirichletViscousTensorEBBC.H"
#include   "NeumannViscousTensorEBBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include   "NeumannViscousTensorDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include   "NeumannConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include   "NeumannConductivityEBBC.H"
#include "CH_Attach.H"

#include <iostream>

#include "UsingNamespace.H"

using std::ifstream;
using std::ios;
/***************/
void amrGodunov(const Box&      a_domain,
                const RealVect& a_dx)
{
  EBAMRCNSParams params;
  int iprob = 0; // 0=inflow-outflow.  1=solid walls
  fillAMRParams(params, iprob);
  
  pout() << "explosion initial conditions" << endl;
  pout() << "solid wall boundary conditions" << endl;
  ParmParse ppgodunov;
  bool tagAllIrregular;
  if(ppgodunov.contains("tag_all_irregular"))
    {
      ppgodunov.get("tag_all_irregular", tagAllIrregular);
    }
  EBAMRCNS::s_noEBCF = tagAllIrregular;
  vector<Real> centerpp(SpaceDim,0.5);
  RealVect center;
  ppgodunov.getarr("explosion_center",centerpp,0,SpaceDim);
  for (int i = 0; i < SpaceDim; i++)
    center[i] = centerpp[i];

  Real size = 0.25;
  ppgodunov.get("explosion_size",size);

  Real p0, r0, p1, r1;
  Real gamma = 1.4;
  ppgodunov.get("gamma",gamma);
  ppgodunov.get("explosion_p0", p0);
  ppgodunov.get("explosion_r0", r0);
  ppgodunov.get("explosion_p1", p1);
  ppgodunov.get("explosion_r1", r1);


  int onedprob;
  ppgodunov.get("one_dim_problem", onedprob);
  bool doOneDOnly = (onedprob != 0);
  RealVect oneDNormal = BASISREALV(1);
  RealVect oneDOrigin = RealVect::Zero;
  if(doOneDOnly)
    {
      vector<Real> voned(SpaceDim, 1.);
      ppgodunov.getarr("one_dim_normal",voned, 0, SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
        oneDNormal[i] = voned[i];
      ppgodunov.getarr("one_dim_origin",voned, 0, SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
        oneDOrigin[i] = voned[i];
    }
  bool useLimiting = true;

  EBExplosionIBCFactory bcfactory(gamma, size, p0, r0, p1, r1, center, onedprob, oneDNormal, oneDOrigin);

  //create patch integrator
  int ifourth, iflatten, iartvisc;
  ppgodunov.get("use_fourth_order_slopes", ifourth);
  ppgodunov.get("use_flattening"         , iflatten);
  ppgodunov.get("use_art_visc"           , iartvisc);
  bool useFourthOrderSlopes = (ifourth  == 1);
  bool useFlattening        = (iflatten == 1);
  bool useArtificialVisc    = (iartvisc == 1);

  RefCountedPtr<EBPatchPolytropicFactory> patchGamma = 
    RefCountedPtr<EBPatchPolytropicFactory>
    (new EBPatchPolytropicFactory(&bcfactory,
                                  gamma,
                                  useFourthOrderSlopes,
                                  useFlattening,
                                  useArtificialVisc,
                                  useLimiting,
                                  doOneDOnly));


  EBAMRCNSFactory amrg_fact(params,
                            patchGamma);

  AMR amr;
  setupAMR(amr, params, amrg_fact, a_domain, false, -1.0);

  // run
  int nstop = 0;
  ppgodunov.get("max_step",nstop);

  Real stopTime = 0.0;
  ppgodunov.get("max_time",stopTime);

  amr.run(stopTime,nstop);

  // output last pltfile and statistics
  //cleanup
  amr.conclude();

}


/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  // registerDebugger();
  // setChomboMPIErrorHandler();
#endif
  {
    // Check for an input file
    char* inFile = NULL;

    EBDebugPoint::s_ivd = IntVect(D_DECL(16,45,0));
    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
    // Parse the command line and the input file (if any)
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    Box coarsestDomain;
    RealVect dx;
    // run amrGodunov
    godunovGeometry(coarsestDomain, dx);

    amrGodunov(coarsestDomain, dx);

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  } //scoping trick

#ifdef CH_MPI
  pout() << "dumping timers" << endl;
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  pout() << "Done." << endl;
}
