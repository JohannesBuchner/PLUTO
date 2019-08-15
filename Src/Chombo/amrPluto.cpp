#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::ifstream;
using std::ios;

// #define HALEM_PROC_SPEED
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "FABView.H"

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"

#include "AMR.H"
#include "AMRLevelPlutoFactory.H"

#include "UsingNamespace.H"

#ifdef CH_LINUX
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

extern "C" {
#include "globals.h"
}

#define NOPT      32

// amrPluto is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrPluto(int, char *av[], char *, Cmd_Line *);

/* ********************************************************************* */
int main(int a_argc, char* a_argv[])
/*
 *
 *
 *
 *
 *
 *
 *
 *********************************************************************** */
{
  Real end_memory, size;
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
#ifdef CH_AIX
  H5dont_atexit();
#endif
//  setChomboMPIErrorHandler();
#endif

  int rank;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  prank = rank;
#else
  rank = 0;
#endif

  char ini_file[128];
  Cmd_Line cmd_line;

  sprintf (ini_file,"pluto.ini");  // default
  ParseCmdLineArgs(a_argc, a_argv, ini_file, &cmd_line);

#ifdef TRAP_FPE
  enableFpExceptions();
#endif

  // Run amrPluto, i.e., do the computation
  amrPluto(a_argc, a_argv, ini_file, &cmd_line);

#ifdef CH_MPI
  // Exit MPI 
  dumpmemoryatexit();
  MPI_Finalize();
#endif

}

/* ********************************************************************* */
void amrPluto(int argc, char *argv[], char *ini_file, Cmd_Line *cmd_line)
/*
 *
 *
 *
 *
 *
 *
 *********************************************************************** */
{
  Real residentSetSize=0.0;
  Real size=0.0;
  
  Runtime runtime;
  Grid    grid;
  RuntimeSetup  (&runtime, cmd_line, ini_file);
  RuntimeSet (&runtime);

  std::string base_name = "pout";
  base_name = runtime.output_dir+string("/")+base_name;
  setPoutBaseName(base_name);

  ShowConfig(argc, argv, ini_file);
 
  pout() << endl << "> Generating grid..." << endl << endl;
  int nghost = GetNghost();
  for (int idim = 0; idim < DIMENSIONS; idim++) {
    grid.nghost[idim]  = nghost;
    grid.np_int[idim]  = grid.np_int_glob[idim] = runtime.npoint[idim];
    grid.np_tot[idim]  = grid.np_tot_glob[idim] = runtime.npoint[idim] + 2*nghost;
    grid.beg[idim]     = grid.gbeg[idim] = grid.lbeg[idim] = nghost;
    grid.end[idim]     = grid.gend[idim] = grid.lend[idim] 
                       = (grid.lbeg[idim] - 1) + grid.np_int[idim];
    grid.lbound[idim]  = runtime.left_bound[idim];
    grid.rbound[idim]  = runtime.right_bound[idim];
  }
  SetGrid (&runtime, &grid);
  
  // VERBOSITY
  
  // This determines the amount of diagnositic output generated
  int verbosity = 1; /* -- change through cmd line args ?? -- */
  CH_assert(verbosity >= 0);
  
  // [Grid]
  
  // Set the physical size of the dimensions of the domain

  Real domainLength = 0.0, xbeg[3], xend[3];
  g_x2stretch = 1.;
  g_x3stretch = 1.;

  for (int dir = 0; dir < DIMENSIONS; dir++){
    xbeg[dir] = runtime.patch_left_node[dir][1];
    xend[dir] = runtime.patch_left_node[dir][2];
    domainLength = MAX(domainLength, xend[dir] - xbeg[dir]);
  }
  domainLength = xend[IDIR] - xbeg[IDIR];

#if (CHOMBO_LOGR == YES)
  domainLength =  log(xend[IDIR]/xbeg[IDIR]);
#endif

#if CH_SPACEDIM > 1
  g_x2stretch = (xend[JDIR] - xbeg[JDIR])/domainLength*(double)runtime.npoint[IDIR]/(double)runtime.npoint[JDIR];
#endif
#if CH_SPACEDIM == 3
  g_x3stretch = (xend[KDIR] - xbeg[KDIR])/domainLength*(double)runtime.npoint[IDIR]/(double)runtime.npoint[KDIR];
#endif

#if GEOMETRY == CARTESIAN
  g_stretch_fact = g_x2stretch*g_x3stretch;
#endif

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);

  D_EXPAND(numCells[0] = runtime.patch_npoint[IDIR][1]; ,
           numCells[1] = runtime.patch_npoint[JDIR][1]; ,
           numCells[2] = runtime.patch_npoint[KDIR][1];)

  CH_assert(D_TERM(   (numCells[0] > 0),
                   && (numCells[1] > 0),
                   && (numCells[2] > 0)));
  CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                   && (numCells[1] % 2 == 0),
                   && (numCells[2] % 2 == 0)));

  // REFINEMENT 
  ParamFileRead (ini_file);

  //Maximum AMR level limit
  int maxLevel = 0;
  maxLevel = atoi(ParamFileGet("Levels",1));
  int numReadLevels = MAX(maxLevel, 1);

  // Refinement ratios between levels
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  std::vector<int> refRatios(numReadLevels + 1);
  int totLevels = 1;
  for (int nlev = 0; nlev <= numReadLevels; nlev++){
    refRatios[nlev] = atoi(ParamFileGet("Ref_ratio",nlev+1));
    if (nlev < maxLevel) totLevels *= refRatios[nlev];
  }
 
  int nprocs = 1;  // number of processors
#ifdef CH_MPI
  MPI_Comm_size(Chombo_MPI::comm, &nprocs);
#endif

  pout() << endl << endl;
  pout() << "> AMR: " << endl << endl;
  pout() << "  Number of levels:      " << maxLevel << endl;
  pout() << "  Equivalent Resolution: " << grid.np_int[IDIR]*totLevels;
  D_EXPAND(                                                  , 
           pout() << " x " << grid.np_int[JDIR]*totLevels;   ,
           pout() << " x " << grid.np_int[KDIR]*totLevels;)
  pout() << endl;
  pout() << "  Number of procs:       " << nprocs << endl << endl;

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals(numReadLevels);
  for (int nlev = 0; nlev < numReadLevels; nlev++)
    regridIntervals[nlev] = atoi(ParamFileGet("Regrid_interval",nlev+1));

  // Threshold that triggers refinement
  Real refineThresh = atof(ParamFileGet("Refine_thresh",1));

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = atoi (ParamFileGet("Tag_buffer_size",1));

  // Minimum dimension of a grid
  int blockFactor = atoi(ParamFileGet("Block_factor",1));

  // Maximum dimension of a grid
  int maxGridSize = atoi(ParamFileGet("Max_grid_size",1));

  NMAX_POINT = maxGridSize + 8;
//  NMAX_POINT = 2*maxGridSize;

  Real fillRatio = atof (ParamFileGet("Fill_ratio",1));

  // [Time] 
//  Real cfl         = runtime.cfl;         // CFL multiplier
//  Real maxDtGrowth = runtime.cfl_max_var; // Limit the time step growth
//  Real stopTime    = runtime.tstop;       // Stop when the simulation time get here
  Real initialDT   = runtime.first_dt;    // Initial time step

  Real dtToleranceFactor = 1.1;  // Let the time step grow by this factor 
                                 // above the "maximum" before reducing it
 
  int nstop = 999999;  // Stop after this number of steps
  if (cmd_line->maxsteps > 0) nstop = cmd_line->maxsteps;

  // [Solver]

  Riemann_Solver *Solver;  // SOLVER
  std::string riem = runtime.solv_type;
  Solver = SetSolver(riem.c_str());
 
  Real fixedDt = -1; // Determine if a fixed or variable time step will be used

  // [Boundary] 
  vector<int> isPeriodica(SpaceDim,0);

  for (int dir = 0; dir < DIMENSIONS; dir++){
    if (runtime.left_bound[dir] == PERIODIC || runtime.right_bound[dir] == PERIODIC){
      isPeriodica[dir] = 1;
    } else {
      isPeriodica[dir] = 0;
    }
  }

  bool isPeriodic[SpaceDim];
  // convert periodic from int->bool
  for (int dim = 0; dim < SpaceDim; dim++)
     {
       isPeriodic[dim] = (isPeriodica[dim] == 1);
       if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
         pout() << "Using Periodic BCs in direction: " << dim << endl;
     }

  // [Output]

  // Set up checkpointing
  Output output;
  GetOutputFrequency (&output, "Checkpoint_interval");

  Real checkpointPeriod   = output.dt;
  Real checkpointClock    = output.dclock;
  int  checkpointInterval = output.dn;

  // Set up plot file writing
  GetOutputFrequency (&output, "Plot_interval");

  Real plotPeriod   = output.dt;
  Real plotClock    = output.dclock;
  int  plotInterval = output.dn;

  if (cmd_line->write == NO) {
    plotPeriod   = checkpointPeriod   = -1.0;
    plotClock    = checkpointClock    = -1.0;
    plotInterval = checkpointInterval = -1;
  }
 
  // [Parameters]

  for (int i = 0; i < USER_DEF_PARAMETERS; ++i) g_inputParam[i] = runtime.aux[i];

  // Initialize PVTE_LAW EOS Table before calling Init()
  
#if EOS == PVTE_LAW && NIONS == 0
  #if TV_ENERGY_TABLE == YES
  MakeInternalEnergyTable();
  #else
  MakeEV_TemperatureTable();
  #endif
  MakePV_TemperatureTable();
#endif
  
  // Call Init once so all processors share global variables
  // assigniments such as g_gamma, g_unitDensity, and so on.
  // This is necessary since not all processors call Init
  // from Startup.

  double u[256];
  Init (u, g_domBeg[0], g_domBeg[1], g_domBeg[2]);

  // Print the parameters
  if ( verbosity >= 2 ) {
    pout() << "maximum step = " << nstop << endl;
    pout() << "maximum time = " << runtime.tstop << endl;

    pout() << "number of cells = " << D_TERM(numCells[0] << "  " <<,
                                             numCells[1] << "  " <<,
                                             numCells[2] << ) endl;

    pout() << "maximum level = " << maxLevel << endl;

    pout() << "refinement ratio = ";
    for (int i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
    pout() << endl;

    pout() << "regrid interval = ";
    for (int i = 0; i < regridIntervals.size(); ++i) pout() << regridIntervals[i] << " ";
    pout() << endl;

    pout() << "refinement threshold = " << refineThresh << endl;

    pout() << "blocking factor = " << blockFactor << endl;
    pout() << "max grid size = " << maxGridSize << endl;
    pout() << "fill ratio = " << fillRatio << endl;

    pout() << "checkpoint interval = " << checkpointInterval << endl;
    pout() << "plot interval = " << plotInterval << endl;
    pout() << "CFL = " << runtime.cfl << endl;
    pout() << "initial dt = " << initialDT << endl;
    if (fixedDt > 0) pout() << "fixed dt = " << fixedDt << endl;
    pout() << "maximum dt growth = " << runtime.cfl_max_var << endl;
    pout() << "dt tolerance factor = " << dtToleranceFactor << endl;
  }

  ShowUnits();

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  initialDT = initialDT*(probDomain.domainBox().bigEnd(0)-probDomain.domainBox().smallEnd(0)+1);

  PatchPluto* patchPluto =  static_cast<PatchPluto*>(new PatchPluto());

  // Set up the patch integrator
  patchPluto->setBoundary(runtime.left_bound, runtime.right_bound);
  patchPluto->setRiemann(Solver);

  // Set up the AMRLevel... factory
  AMRLevelPlutoFactory amrGodFact;

  amrGodFact.define(runtime.cfl,
                    domainLength,
                    verbosity,
                    refineThresh,
                    tagBufferSize,
                    initialDT,
                    patchPluto);

  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrGodFact);

  if (fixedDt > 0)  amr.fixedDt(fixedDt);

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  amr.blockFactor(blockFactor);
  amr.fillRatio(fillRatio);

  // The hyperbolic codes use a grid buffer of 1
  amr.gridBufferSize(1);

  // Set output parameters
  amr.checkpointPeriod(checkpointPeriod);
  amr.checkpointClock (checkpointClock);
  amr.checkpointInterval(checkpointInterval);
  amr.plotPeriod  (plotPeriod);
  amr.plotClock   (plotClock);
  amr.plotInterval(plotInterval);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(runtime.cfl_max_var);
  amr.dtToleranceFactor(dtToleranceFactor);

  // Setup output files
  std::string prefix = "data.";
  prefix = runtime.output_dir+string("/")+prefix;
  amr.plotPrefix(prefix);

  // Setup checkpoint files
  prefix = "chk.";
  prefix = runtime.output_dir+string("/")+prefix;
  amr.checkpointPrefix(prefix);

  amr.verbosity(verbosity);

  // Setup input files
  if (cmd_line->restart == NO && cmd_line->h5restart == NO){
    if (!ParamExist("fixed_hierarchy")) {

      // initialize from scratch for AMR run
      // initialize hierarchy of levels
      amr.setupForNewAMRRun();
    } else       {
      MayDay::Error("Fixed_hierarchy option still disabled");
    }
  } else { 
    char   restartNumber[32];
    sprintf (restartNumber, "%04d.hdf5", cmd_line->nrestart);
    prefix += restartNumber;

    pout() << endl << "> Restarting from file " << prefix << endl;

#ifdef CH_USE_HDF5
    HDF5Handle handle(prefix,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
    // Create dummies to be passed to startup
    // startup should be called at restart too
/*
    Box dummyBox (IntVect::Zero,
                  IntVect(D_DECL(NMAX_POINT-1,
                                 NMAX_POINT-1,
                                 NMAX_POINT-1)));
    FArrayBox dummyData(dummyBox,NVAR);
    dummyData.setVal(0.0);
    patchPluto->startup(dummyData);
*/
#else
    MayDay::Error("amrPluto restart only defined with hdf5");
#endif
  }

  time_t tbeg, tend;
  double exec_time;

  // Run the computation
  pout() << "> Starting Computation" << endl;
  time (&tbeg);

  amr.run(runtime.tstop,nstop);

  time (&tend);
  exec_time = difftime (tend, tbeg);
  pout() << "> Elapsed time: " << exec_time << " s" << endl;
  pout() << "> Done\n";

  // Output the last plot file and statistics
  amr.conclude();
}

/* ******************************************************************* */
void print (const char *fmt, ...)
/* 
 * 
 * General purpose print function used by PLUTO-Chombo.
 * The corresponding static grid version is given in tools.c
 *
 ********************************************************************* */
{
  char buffer[256];
  va_list args;
  va_start (args, fmt);
  vsprintf (buffer,fmt, args);
  pout() << buffer;
}

#undef NOPT
