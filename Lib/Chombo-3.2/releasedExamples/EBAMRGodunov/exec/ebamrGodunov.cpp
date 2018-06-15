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

#include "EBPlanarShockIBCFactory.H"
#include "ModianoIBCFactory.H"
#include "EBPatchPolytropicFactory.H"

#include "EBPatchPolytropic.H"

#include "EBAMRGodunovFactory.H"
#include "EBAMRGodunov.H"
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
#include "GodunovGeom.H"
#include "CH_Attach.H"
#include "EBAMRGLoadBalance.H"

#include <iostream>

#define CHECK_FLOATING_PT 0

#if CHECK_FLOATING_PT==1
#include <fenv.h>
#endif

using std::ifstream;
using std::ios;

#include "UsingNamespace.H"

void amrGodunov(const Box&      a_coarsestDomain,
                const RealVect& a_dx);

/***************/
int
main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  // registerDebugger();
  // setChomboMPIErrorHandler();
#endif
  { //scoping trick
#if CHECK_FLOATING_PT==1
  //  int except =  FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW |  FE_INVALID ;
  int except =  FE_DIVBYZERO |  FE_OVERFLOW |  FE_INVALID ;
  feenableexcept(except);
#endif

 // Check for an input file
  char* inFile = NULL;

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
  //IntVectSet::setMaxDense(10000000);
  godunovGeometry(coarsestDomain, dx);

  amrGodunov(coarsestDomain, dx);
  }
#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
/***************/
void amrGodunov(const Box&      a_domain,
                const RealVect& a_dx)
{
  // read inputs
  ParmParse ppgodunov;

  int verbosity;
  ppgodunov.get("verbosity",verbosity);
 CH_assert(verbosity >= 0);

  Real gamma = 1.4;
  ppgodunov.get("gamma",gamma);


  int nstop = 0;
  ppgodunov.get("max_step",nstop);

  int redistRad = 0;
  ppgodunov.get("redist_radius",redistRad);

  Real stopTime = 0.0;
  ppgodunov.get("max_time",stopTime);

  bool is_periodic[SpaceDim];
  for (int i = 0; i < SpaceDim; ++i)
    {
      is_periodic[i] = false;
    }

  int max_level = 0;
  ppgodunov.get("max_level",max_level);

  int num_read_levels = Max(max_level,1);
  std::vector<int> ref_ratios; // (num_read_levels,1);
  // note this requires a ref_ratio to be defined for the
  // finest level (even though it will never be used)
  ppgodunov.getarr("ref_ratio",ref_ratios,0,num_read_levels+1);

  std::vector<int> regrid_intervals; // (num_read_levels,1);
  ppgodunov.getarr("regrid_interval",regrid_intervals,0,num_read_levels);

  Real refineThresh = 0.3;
  ppgodunov.get ("refine_thresh",refineThresh);

  int block_factor = 1;
  ppgodunov.get("block_factor",block_factor);

  int tagBufferSize = 3;
  ppgodunov.get("tag_buffer_size",tagBufferSize);

  int max_grid_size = 32;
  ppgodunov.get("max_grid_size",max_grid_size);

  Real fill_ratio = 0.75;
  ppgodunov.get("fill_ratio",fill_ratio);

  int checkpoint_interval = 0;
  ppgodunov.get("checkpoint_interval",checkpoint_interval);

  int plot_interval = 0;
  ppgodunov.get("plot_interval",plot_interval);

  Real cfl = 0.8;
  ppgodunov.get("cfl",cfl);

  Real initialCFL = 0.1;
  ppgodunov.get("initial_cfl",initialCFL);

  Real fixed_dt = -1;
  ppgodunov.get("fixed_dt",fixed_dt);

  Real max_dt_growth = 1.1;
  ppgodunov.get("max_dt_growth",max_dt_growth);

  Real dt_tolerance_factor = 1.1;
  ppgodunov.get("dt_tolerance_factor",dt_tolerance_factor);

  vector<Real> length(SpaceDim);
  ppgodunov.getarr("domain_length",length,0,SpaceDim);
  RealVect domainLength;
  for (int idir=0;idir<SpaceDim;idir++)
    {
      domainLength[idir] = length[idir];
    }

  int ifourth, iflatten, iartvisc;
  ppgodunov.get("use_fourth_order_slopes", ifourth);
  ppgodunov.get("use_flattening"         , iflatten);
  ppgodunov.get("use_art_visc"           , iartvisc);
  bool useFourthOrderSlopes = (ifourth  ==1);
  bool useFlattening        = (iflatten ==1);
  bool useArtificialVisc    = (iartvisc ==1);

  ProblemDomain prob_domain(a_domain.smallEnd(),
                            a_domain.bigEnd(),
                            is_periodic);

  // create initial/boundary condition object

  pout() << "planar shock initial conditions" << endl;
  Real ms;
  ppgodunov.get("shock_mach",ms);

  Real center;
  ppgodunov.get("shock_center",center);

  int inormal;
  ppgodunov.get("shock_normal",inormal);

  int ishockback;
  ppgodunov.get("shock_backward",ishockback);
  bool shockbackward = (ishockback == 1);
  //force these for now.
  bool doSmushing = true;
  bool doRZCoords = false;
  bool hasSourceTerm = false;
  bool useLimiting = true;

  EBPlanarShockIBCFactory
    bcfactory(gamma, ms, center, inormal, shockbackward, doRZCoords);

  int iusemassredist;
  ppgodunov.get("use_mass_redist", iusemassredist);
  bool useMassRedist    = (iusemassredist ==1);

  //create patch integrator
  EBPatchPolytropicFactory patchGamma(&bcfactory,
                                      gamma,
                                      useFourthOrderSlopes,
                                      useFlattening,
                                      useArtificialVisc,
                                      useLimiting,
                                      doRZCoords);

  EBAMRGodunovFactory amrg_fact(initialCFL,
                                cfl,
                                redistRad,
                                domainLength,
                                refineThresh,
                                tagBufferSize,
                                verbosity,
                                useMassRedist,
                                doSmushing,
                                doRZCoords,
                                hasSourceTerm,
                                &patchGamma);

  bool useTimedLoadBalance = false;
  ppgodunov.query("use_timed_load_balance", useTimedLoadBalance);
  if (useTimedLoadBalance)
    {
      pout() << "using timed load balance"  << endl;
      EBAMRGodunov::setLoadBalance(EBAMRGLoadBalance);
    }

  AMR amr;
  amr.define(max_level,ref_ratios,prob_domain,&amrg_fact);

  if (fixed_dt > 0)
    {
      amr.fixedDt(fixed_dt);
    }

  // set grid generation parameters
  amr.maxGridSize(max_grid_size);
  amr.blockFactor(block_factor);
  amr.fillRatio(fill_ratio);

  // the hyperbolic codes use a grid buffer of 1
  int gridBufferSize;
  ppgodunov.get("grid_buffer_size",gridBufferSize);
  amr.gridBufferSize(gridBufferSize);

  // set output parameters
  amr.checkpointInterval(checkpoint_interval);
  amr.plotInterval(plot_interval);
  amr.regridIntervals(regrid_intervals);
  amr.maxDtGrow(max_dt_growth);
  amr.dtToleranceFactor(dt_tolerance_factor);

  if (ppgodunov.contains("use_subcycling"))
    {
      bool useSubcycling;
      ppgodunov.get("use_subcycling", useSubcycling);
      if (!useSubcycling)
        {
          pout() << "SUBCYCLING IN TIME TURNED OFF!!!"  << endl;
        }
      amr.useSubcyclingInTime(useSubcycling);
    }
  if (ppgodunov.contains("plot_prefix"))
    {
      std::string prefix;
      ppgodunov.get("plot_prefix",prefix);
      amr.plotPrefix(prefix);
    }

  if (ppgodunov.contains("chk_prefix"))
    {
      std::string prefix;
      ppgodunov.get("chk_prefix",prefix);
      amr.checkpointPrefix(prefix);
    }

  amr.verbosity(verbosity);

  if (!ppgodunov.contains("restart_file"))
    {
      if (!ppgodunov.contains("fixed_hierarchy"))
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
          amr.setupForNewAMRRun();
        }
      else
        {
          std::string grid_file;
          ppgodunov.get("fixed_hierarchy",grid_file);

          // initialize from a list of grids in "grid_file"
          Vector<Vector<Box> > amrGrids(max_level+1);
          godunovFixedGrids(amrGrids,
                            prob_domain,
                            max_level,
                            max_grid_size,
                            block_factor,
                            verbosity,
                            grid_file);
          amr.setupForFixedHierarchyRun(amrGrids,1);
        }
    }
  else
    {
      std::string restart_file;
      ppgodunov.get("restart_file",restart_file);
      pout() << " restarting from file " << restart_file << endl;

#ifdef CH_USE_HDF5
      HDF5Handle handle(restart_file,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      amr.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
    }

  // run
  amr.run(stopTime,nstop);

  // output last pltfile and statistics
  //cleanup
  amr.conclude();

}


