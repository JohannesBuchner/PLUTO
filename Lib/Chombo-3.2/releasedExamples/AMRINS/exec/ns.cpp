#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*****************/
/*****************/

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include "BRMeshRefine.H"
#include "AMRLevel.H"
#include "AMR.H"
#include "AMRIO.H"
#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "AMRNavierStokes.H"
#include "AMRNavierStokesFactory.H"
#include "PhysBCUtil.H"
#include "channelBC.H"
#include "ProblemDomain.H"
#include "AMRPoissonOp.H" // only for s_relaxMode
#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#ifndef CH_NTIMER
#include "memusage.H"
#endif
#include "CH_Timer.H"

#include "UsingNamespace.H"

// define main loop here because that's the way amrgodunov does it!
void mainLoop();

//#define TRAP_FPE  //(should be off by default)
#ifdef TRAP_FPE
// Previous versions of glibc require the following code:
extern "C"
{
#include <fpu_control.h>
}
/* IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask */
static void __attribute__ ((constructor)) trapfpe(void)
{
  pout() << " Turning on floating-point traps! " << endl;
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM);
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_OM | _FPU_MASK_UM);
  fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_UM);
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
  //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_DM | _FPU_MASK_UM);
  //fpu_control_t cw = _FPU_DEFAULT;
   _FPU_SETCW(cw);
   /* On x86, this expands to: */
   /* unsigned int cw = 0x037f & ~(0x01 | 0x04 | 0x08); */
   /* __asm__ ("fldcw %0" : : "m" (*&cw));              */
}
#endif

int main(int argc, char* argv[])
{
#ifdef TRAP_FPE
  trapfpe();
#endif

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#ifdef CH_AIX
  H5dont_atexit();
#endif
  // setChomboMPIErrorHandler();
  MPI_Barrier(Chombo_MPI::comm);  // Barrier #1
#endif

  // read inputs
  char* in_file = argv[1];
  ParmParse pp(argc-2, argv+2,NULL, in_file);
  AMRPoissonOp::s_relaxMode = 1;
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  // run actual algorithm
  mainLoop();

#ifndef CH_NTIMER
  Real end_memory = get_memory_usage_from_OS();

  pout() << "Everything done.  mem= " << end_memory << " MB" << endl;
#endif

#ifdef CH_MPI
  dumpmemoryatexit();
#endif

#ifdef CH_MPI
  MPI_Finalize();
#endif
}

void
mainLoop()
{
  //ReadInput.start();
  // read inputs
  ParmParse ppMain("main");

  int verbosity = 0;
  ppMain.query("verbosity", verbosity);
  CH_assert (verbosity >= 0);

  int nstop = 0;
  ppMain.query("max_step", nstop);

  Real stopTime = 0.;
  ppMain.query("max_time", stopTime);

  Real fixed_dt = -1;
  ppMain.query("fixed_dt", fixed_dt);

  vector<int> tempVectInt;
  Vector<int> num_cells(SpaceDim,0);
  ppMain.queryarr("num_cells", tempVectInt, 0, SpaceDim);
  num_cells = tempVectInt;

  CH_assert ( D_TERM ( ( num_cells[0] > 0 ), &&
                    ( num_cells[1] > 0 ), &&
                    ( num_cells[2] > 0 ) ) );
  CH_assert ( D_TERM ( ( num_cells[0] % 2 == 0 ), &&
                    ( num_cells[1] % 2 == 0 ), &&
                    ( num_cells[2] % 2 == 0 ) ) );

  vector<int> is_periodica(SpaceDim,0);
  bool is_periodic[SpaceDim];

  ppMain.queryarr("is_periodic", is_periodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
    {
      is_periodic[dim] = (is_periodica[dim] == 1);
      if (is_periodic[dim] && verbosity >= 2 && procID() == 0)
        pout() << "Using periodic BCs in direction:" << dim << endl;
    }

  int max_level = 0;
  ppMain.query("max_level", max_level);

  int num_read_levels = Max(max_level,1);
  Vector<int> ref_ratios (num_read_levels+1,1);
  ppMain.queryarr("ref_ratio", tempVectInt,0,num_read_levels+1);
  ref_ratios = tempVectInt;

  Vector<int> regrid_intervals(num_read_levels,1);
  ppMain.queryarr("regrid_interval", tempVectInt,0,num_read_levels);
  regrid_intervals = tempVectInt;

  Real refine_threshold = 0.3;
  ppMain.query("refine_threshold", refine_threshold);

  int block_factor = 2;
  ppMain.query("block_factor", block_factor);

  int max_grid_size = 32;
  ppMain.query("max_grid_size", max_grid_size);

  int max_base_grid_size = max_grid_size;
  ppMain.query("max_base_grid_size", max_base_grid_size);

  Real fill_ratio = 0.75;
  ppMain.query("fill_ratio", fill_ratio);

  int checkpoint_interval = -1;
  ppMain.query("checkpoint_interval", checkpoint_interval);

  int plot_interval = -1;
  ppMain.query("plot_interval", plot_interval);

  Real max_dt_grow = 1.5;
  ppMain.query("max_dt_grow", max_dt_grow);

  bool predefinedGrids = ppMain.contains("gridfile");
  string gridfile;
  if (predefinedGrids)
  {
    ppMain.get("gridfile", gridfile);
  }

  string checkPrefix = "chk";
  ppMain.query("checkPrefix", checkPrefix);

  string plotPrefix = "plt";
  ppMain.query("plotPrefix", plotPrefix);

  Real cfl = 0.5;
  ppMain.query("cfl", cfl);

  // solver parameter which will be passed all the way throughout the code
  int tempBool = 0;
  ppMain.query("limitSolverCoarsening", tempBool);
  bool limitSolverCoarsening = (tempBool == 1);

  if (verbosity >= 2)
  {
    // display inputs
    pout() << "nstop = " << nstop << endl;
    pout() << "stopTime = " << stopTime << endl;
    if (fixed_dt > 0)
      pout() << "fixed dt = " << fixed_dt << endl;

    pout() << "num_cells = " << D_TERM ( num_cells[0] << "  " <<,
                                       num_cells[1] << "  " <<,
                                       num_cells[2] << ) endl;

    pout() << "max level = " << max_level << endl;
    pout() << "refinement ratios = ";
    for (int i = 0; i < ref_ratios.size(); ++i) pout() << ref_ratios[i] << " ";
    pout() << endl;
    pout() << "regrid intervals = ";
    for (int i = 0; i < regrid_intervals.size(); ++i) pout()
      << regrid_intervals[i] << " ";
    pout() << endl;
    pout() << "refinement threshold = " << refine_threshold << endl;

    pout() << "blocking factor = " << block_factor << endl;
    pout() << "max grid size = " << max_grid_size << endl;
    pout() << "max base grid size = " << max_base_grid_size << endl;
    pout() << "fill ratio = " << fill_ratio << endl;

    pout() << "checkpoint interval = " << checkpoint_interval << endl;
    pout() << "plot interval = " << plot_interval << endl;
    pout() << "CFL = " << cfl << endl;
    pout() << "maxDtGrow = " << max_dt_grow << endl;
    pout() << "plotPrefix = " << plotPrefix << endl;
    pout() << "checkPrefix = " << checkPrefix << endl;
    pout() << "verbosity = " << verbosity << endl;
    pout() << "limit solver coarsening = " << limitSolverCoarsening << endl;

    if (predefinedGrids)
    {
      pout() << "gridFile = " << gridfile << endl;
    }
  }

  // create objects
  ProblemDomain prob_domain(IntVect::Zero,
                            IntVect(D_DECL(num_cells[0]-1,
                                           num_cells[1]-1,
                                           num_cells[2]-1)),
                            is_periodic);

  Vector<Vector<Box> > amrGrids;
  if (predefinedGrids)
    {
#ifdef CH_MPI
      // if (procID() ==  uniqueProc(SerialTask::compute))
      MPI_Barrier(Chombo_MPI::comm);
      if (procID() ==  0)
        {
#endif
          // read in predefined grids
          ifstream is(gridfile.c_str(), ios::in);
          if (is.fail())
            {
              pout() << "cannot open grids file " << gridfile << endl;
              MayDay::Error();
            }

#ifdef CH_MPI
      if (verbosity >= 3)
        {
          pout() << "procID: " << procID() << "  opening gridfile \n" << endl;
        }
#endif

      // format of file -- number of levels, then for
      // each level starting with level 1, number of
      // grids on level, list of boxes
      int in_numLevels;
      is >> in_numLevels;
      CH_assert (in_numLevels <= max_level+1);
      if (verbosity >= 3)
        {
          pout() << "numLevels = " << in_numLevels << endl;
        }
      while (is.get() != '\n');
      amrGrids.resize(in_numLevels);
      // check to see if coarsest level needs to be broken up
      domainSplit(prob_domain,amrGrids[0], max_base_grid_size,4);

      if (verbosity >= 3)
        {
          pout() << "level 0: ";
          for (int n=0; n<amrGrids[0].size(); n++)
            {
              pout () << amrGrids[0][n] << endl;
            }
        }
      // now loop over levels, starting with level 1
      int ngrid;
      for (int lev=1; lev<in_numLevels; lev++)
        {
          is >> ngrid;
          if (verbosity >= 3)
            {
              pout() << "level " << lev << " numGrids = " << ngrid << endl;
              pout() << "Grids: ";
            }
          while (is.get() != '\n');
          amrGrids[lev].resize(ngrid);
          for (int i=0; i<ngrid; i++)
            {
              Box bx;
              is >> bx;

              // advance to next box
              while (char ch = is.get())
                {
                  if (ch == '#') break;
                  if (ch == '\n') break;
                }

              // quick check on box size
              Box bxRef(bx);
              // not really sure why i was doing this (holdover from
              // legacy code)
              //bxRef.refine(ref_ratios[lev-1]);
              if (bxRef.longside() > max_grid_size)
                {
                  pout() << "Grid " << bx << " too large" << endl;
                  MayDay::Error();
                }
              if (verbosity >= 3)
                {
                  pout() << bx << endl;
                }
              amrGrids[lev][i] = bx;
            } // end loop over boxes on this level
        } // end loop over levels

#ifdef CH_MPI
        } // end if procID = 0
      //broadcast (amrGrids, uniqueProc(SerialTask::compute));
    broadcast (amrGrids, 0);
#endif

    } // end if grids are predefined
  //ReadInput.stop();

  //pout() << " Beginning AMR Setup.  mem=" << get_memory_usage_from_OS() << "MB" << endl;

  AMRNavierStokesFactory amrns_fact;
  amrns_fact.CFL(cfl);
  amrns_fact.refinementThreshold(refine_threshold);

  channelBC physBC;
  amrns_fact.setPhysBC(physBC);

  amrns_fact.setLimitSolverCoarsening(limitSolverCoarsening);

  AMR thisAMR;

  thisAMR.verbosity(verbosity);
  AMRLevel::verbosity(verbosity);

  thisAMR.define(max_level, ref_ratios, prob_domain,
                 &amrns_fact);

  thisAMR.checkpointInterval(checkpoint_interval);
  thisAMR.plotInterval(plot_interval);
  thisAMR.maxGridSize(max_grid_size);
  thisAMR.maxBaseGridSize(max_base_grid_size);
  thisAMR.fillRatio(fill_ratio);
  thisAMR.blockFactor(block_factor);
  thisAMR.regridIntervals(regrid_intervals);
  thisAMR.maxDtGrow(max_dt_grow);
  if (fixed_dt > 0) thisAMR.fixedDt(fixed_dt);
  thisAMR.checkpointPrefix(checkPrefix);
  thisAMR.plotPrefix(plotPrefix);

  if (!ppMain.contains("restart_file") )
    {
      if (!predefinedGrids)
        {
          // initialize from scratch for an AMR run
          // initialize hierarchy of levels
          thisAMR.setupForNewAMRRun();
        }
      else
        {
          thisAMR.setupForFixedHierarchyRun(amrGrids, 1);
        }
    }
  else
    {
      // initialize from restart file
      std::string restart_file;
      ppMain.query("restart_file", restart_file);
#ifdef CH_USE_HDF5
      HDF5Handle handle(restart_file, HDF5Handle::OPEN_RDONLY);
      thisAMR.setupForRestart(handle);
      handle.close();
#else
      MayDay::Error("amrNavierStokes restart only defined with HDF5");
#endif
      if (verbosity >= 2)
        {
          pout() << "restart_file = " << restart_file << endl;
        }
    } // end if restart

  // run until conclusion
  thisAMR.run(stopTime, nstop);

#ifndef CH_NTIMER
  pout() << "AMR Run completed.  mem= " << get_memory_usage_from_OS()
         << "MB" << endl;

#endif

  // output last plotfile and statistics
  thisAMR.conclude();

}
