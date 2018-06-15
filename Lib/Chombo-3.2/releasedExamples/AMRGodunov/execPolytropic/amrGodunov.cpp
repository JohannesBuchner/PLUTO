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
#include "AMRLevel.H"
#include "AMRLevelPolytropicGasFactory.H"
#include "AMRLevelPolytropicGas.H"

#include "PolytropicPhysics.H"

#include "RampIBC.H"
#include "ExplosionIBC.H"
#include "WaveIBC.H"

#include "UsingNamespace.H"



// Possible pressure relationships for the initial condition
#define PRESSURE_ISENTROPIC 0
#define PRESSURE_CONSTANT   1

// amrGodunov is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void amrGodunov();

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
#ifdef CH_AIX
  H5dont_atexit();
#endif
  // setChomboMPIErrorHandler();
#endif

  int rank, number_procs;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
  rank = 0;
  number_procs = 1;
#endif

  if (rank == 0)
    {
      pout() << " number_procs = " << number_procs << endl;
    }

  // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
    {
      inFile = a_argv[1];
    }
  else
    {
      pout() << "Usage:  amrGodunov...ex <inputfile>" << endl;
      pout() << "No input file specified" << endl;
      return -1;
    }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);


  // Run amrGodunov, i.e., do the computation
  amrGodunov();

#ifdef CH_MPI
  // Exit MPI
  CH_TIMER_REPORT();
  pout() << "dumping memory info" << endl;
  dumpmemoryatexit();
  pout() << "mpi finalizing" << endl;
  MPI_Finalize();
#endif
  pout() << "leaving " << endl;
  return 0;
}

void amrGodunov()
{
  // Read inputs that are prefixed with "godonouv."
  ParmParse ppgodunov("godunov");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppgodunov.query("verbosity",verbosity);
  CH_assert(verbosity >= 0);

  // For all gas dynamics
  Real gamma = 1.4;
  ppgodunov.get("gamma",gamma);

  // Stop after this number of steps
  int nstop = 0;
  ppgodunov.get("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppgodunov.get("max_time",stopTime);

  // Set the physical size of the longest dimension of the domain
  Real domainLength = 1.0;
  ppgodunov.get("domain_length",domainLength);

  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i)
    {
      numCells[i] = 0;
    }
  ppgodunov.getarr("num_cells",numCells,0,SpaceDim);

  CH_assert(D_TERM(   (numCells[0] > 0),
                   && (numCells[1] > 0),
                   && (numCells[2] > 0)));
  CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                   && (numCells[1] % 2 == 0),
                   && (numCells[2] % 2 == 0)));

  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppgodunov.getarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim = 0; dim < SpaceDim; dim++)
    {
      isPeriodic[dim] = (isPeriodica[dim] == 1);
      if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
        {
          pout() << "Using Periodic BCs in direction: " << dim << endl;
        }
    }

  // Maximum AMR level limit
  int maxLevel = 0;
  ppgodunov.get("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppgodunov.getarr("ref_ratio",refRatios,0,numReadLevels+1);

  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  ppgodunov.getarr("regrid_interval",regridIntervals,0,numReadLevels);

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppgodunov.get("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineThresh = 0.3;
  ppgodunov.get ("refine_thresh",refineThresh);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppgodunov.get("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppgodunov.get("max_grid_size",maxGridSize);

  Real fillRatio = 0.75;
  ppgodunov.get("fill_ratio",fillRatio);

  // Order of the normal predictor (CTU -> 0, PLM -> 1, PPM -> 2)
  std::string normalPred;
  int normalPredOrder;
  ppgodunov.get("normal_predictor",normalPred);
  if (normalPred == "CTU" || normalPred == "ctu")
    {
      normalPredOrder = 0;
    }
  else if (normalPred == "PLM" || normalPred == "plm")
    {
      normalPredOrder = 1;
    }
  else if (normalPred == "PPM" || normalPred == "ppm")
    {
      normalPredOrder = 2;
    }
  else
    {
      MayDay::Error("Normal predictor must by PLM or PPM");
    }

  // Use fourth order slopes
  int inFourthOrderSlopes = 1;
  bool useFourthOrderSlopes;
  ppgodunov.get("use_fourth_order_slopes",inFourthOrderSlopes);
  useFourthOrderSlopes = (inFourthOrderSlopes == 1);

  // Do slope limiting
  int inPrimLimiting = 1;
  bool usePrimLimiting;
  ppgodunov.get("use_prim_limiting",inPrimLimiting);
  usePrimLimiting = (inPrimLimiting == 1);

  // Do slope limiting using characteristics
  int inCharLimiting = 0;
  bool useCharLimiting;
  ppgodunov.get("use_char_limiting",inCharLimiting);
  useCharLimiting = (inCharLimiting == 1);

  // Do slope flattening
  int inFlattening = 1;
  bool useFlattening;
  ppgodunov.get("use_flattening",inFlattening);
  useFlattening = (inFlattening == 1);

  // Apply artificial viscosity
  int inArtificialViscosity = 1;
  bool useArtificialViscosity;
  ppgodunov.get("use_artificial_viscosity",inArtificialViscosity);
  useArtificialViscosity = (inArtificialViscosity == 1);

  // Artificial viscosity coefficient/multiplier
  Real artificialViscosity = 0.1;
  ppgodunov.get("artificial_viscosity",artificialViscosity);

  // Don't use high-order limiter by default
  int inHighOrderLimiter = 0;
  bool highOrderLimiter;
  ppgodunov.query("high_order_limiter", inHighOrderLimiter);
  highOrderLimiter = (inHighOrderLimiter == 1);

  // Set up checkpointing
  int checkpointInterval = 0;
  ppgodunov.query("checkpoint_interval",checkpointInterval);

  // Set up plot file writing
  int plotInterval = 0;
  ppgodunov.query("plot_interval",plotInterval);

  // CFL multiplier
  Real cfl = 0.8;
  ppgodunov.get("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppgodunov.get("initial_cfl",initialCFL);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppgodunov.query("fixed_dt",fixedDt);

  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  ppgodunov.get("max_dt_growth",maxDtGrowth);

  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  ppgodunov.get("dt_tolerance_factor",dtToleranceFactor);

  // Create and define IBC (initial and boundary condition) object
  PhysIBC* ibc;

  // A minimum pressure needed to construct PolytropicPhysics - used in slope
  // flattening
  Real smallPressure;

  // Don't use source term by default
  bool useSourceTerm = false;

  // Source term multiplier
  Real sourceTermScaling = 0.0;

  // Determine the sample problem specified
  std::string problemString;
  if (ppgodunov.contains("problem"))
  {
    ppgodunov.query("problem",problemString);

    // Print some parameters
    if (verbosity >= 2)
    {
      pout() << "problem = " << problemString << endl;
      pout() << "gamma = " << gamma << endl;
    }

    if (problemString == "ramp")
    {
      if (isPeriodic[0] || isPeriodic[1])
      {
        MayDay::Error("Neither x or y boundaries can be periodic");
      }

      // Ramp problem
      Real alpha = 30.0;
      ppgodunov.get("angle_deg",alpha);

      Real ms = 10.0;
      ppgodunov.get("shock_mach",ms);

      Real xcorner = 0.1;
      ppgodunov.get("xcorner",xcorner);

      if (verbosity >= 2)
      {
        pout() << "alpha = " << alpha << endl;
        pout() << "shock_mach = " << ms << endl;
        pout() << "xcorner = " << xcorner << endl;
      }

      // Define IBC for ramp problem
      RampIBC* rampibc = new RampIBC;
      rampibc->setFortranCommon(smallPressure,
                                gamma,
                                alpha,
                                ms,
                                xcorner,
                                artificialViscosity);
      ibc = rampibc;
    }
    else if (problemString == "explosion")
    {
      // Explosion problem
      Real ms = 10.0;
      ppgodunov.get("shock_mach",ms);

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      Real size = 0.25;
      ppgodunov.get("initial_size",size);

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "shock_mach = " << ms << endl;
        pout() << "initial_center = " << D_TERM(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << ) endl;
        pout() << "initial_size = " << size << endl;
        pout() << "initial_velocity = " << D_TERM(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << ) endl;
      }

      // Define IBC for explosion problem
      ExplosionIBC* explosionibc = new ExplosionIBC;
      explosionibc->setFortranCommon(smallPressure,
                                     gamma,
                                     ms,
                                     center,
                                     size,
                                     velocity,
                                     artificialViscosity);
      ibc = explosionibc;
    }
    else if (problemString == "wave")
    {
      // Plane wave problem
      Real ambientDensity = 1.4;
      ppgodunov.get("ambient_density",ambientDensity);

      Real deltaDensity = 0.014;
      ppgodunov.get("delta_density",deltaDensity);

      int pressure = -1;
      std::string pressureString;
      ppgodunov.get("initial_pressure",pressureString);
      if (pressureString == "isentropic")
      {
        pressure = PRESSURE_ISENTROPIC;
      }
      else if (pressureString == "constant")
      {
        pressure = PRESSURE_CONSTANT;
      }

      if (pressure == -1)
      {
        pout() << "Invalid pressure, \"" << pressureString << "\", specified in input file" << endl << endl;
        return;
      }
      else
      {
        pressure = PRESSURE_ISENTROPIC;
      }

      vector<int> waveNumberpp(SpaceDim,0);
      waveNumberpp[0] = 1;
      IntVect waveNumber;
      ppgodunov.getarr("wave_number",waveNumberpp,0,SpaceDim);
      int norm2 = 0;
      for (int i = 0; i < SpaceDim; i++)
      {
        waveNumber[i] = waveNumberpp[i];
        norm2 += waveNumber[i]*waveNumber[i];
      }

      if (norm2 == 0)
      {
        pout() << "One component of the wave number must be non-zero" << endl << endl;
        return;
      }

      vector<Real> centerpp(SpaceDim,0.5);
      RealVect center;
      ppgodunov.getarr("initial_center",centerpp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        center[i] = centerpp[i];
      }

      vector<Real> velocitypp(SpaceDim,0.0);
      RealVect velocity;
      ppgodunov.getarr("initial_velocity",velocitypp,0,SpaceDim);
      for (int i = 0; i < SpaceDim; i++)
      {
        velocity[i] = velocitypp[i];
      }

      if (verbosity >= 2)
      {
        pout() << "ambient_density = " << ambientDensity << endl;
        pout() << "delta_density = " << deltaDensity << endl;
        pout() << "initial_pressure = " << pressureString << endl;
        pout() << "wave_number = " << D_TERM(waveNumber[0] << "  " <<,
                                             waveNumber[1] << "  " <<,
                                             waveNumber[2] << ) endl;
        pout() << "initial_center = " << D_TERM(center[0] << "  " <<,
                                                center[1] << "  " <<,
                                                center[2] << ) endl;
        pout() << "initial_velocity = " << D_TERM(velocity[0] << "  " <<,
                                                  velocity[1] << "  " <<,
                                                  velocity[2] << ) endl;
      }

      // Define IBC for plane wave problem
      WaveIBC* waveibc = new WaveIBC;
      waveibc->setFortranCommon(smallPressure,
                                gamma,
                                ambientDensity,
                                deltaDensity,
                                pressure,
                                waveNumber,
                                center,
                                velocity,
                                artificialViscosity);
      ibc = waveibc;
    }
    else
    {
      // The sample problem name given isn't valid
      pout() << "Invalid problem, \"" << problemString << "\", specified in input file" << endl << endl;
      return;
    }
  }
  else
  {
    // A sample problem must be specified
    pout() << "\"godunov.problem\" not specified in input file" << endl << endl;
    return;
  }

  if (verbosity >= 2)
  {
    pout() << "verbosity = " << verbosity << endl;

    pout() << "maximum_step = " << nstop << endl;
    pout() << "maximum_time = " << stopTime << endl;
    if (fixedDt > 0)
    {
      pout() << "fixed_dt = " << fixedDt << endl;
    }

    pout() << "number_of_cells = " << D_TERM(numCells[0] << "  " <<,
                                             numCells[1] << "  " <<,
                                             numCells[2] << ) endl;
    pout() << "is_period = " << D_TERM(isPeriodic[0] << "  " <<,
                                       isPeriodic[1] << "  " <<,
                                       isPeriodic[2] << ) endl;

    pout() << "maximum_level = " << maxLevel << endl;
    pout() << "refinement_ratio = ";
    for (int i = 0; i < refRatios.size(); ++i)
    {
      pout() << refRatios[i] << " ";
    }
    pout() << endl;

    pout() << "regrid_interval = ";
    for (int i = 0; i < regridIntervals.size(); ++i)
    {
      pout() << regridIntervals[i] << " ";
    }
    pout() << endl;
    pout() << "tag_buffer_size = " << tagBufferSize << endl;

    pout() << "refinement_threshold = " << refineThresh << endl;

    pout() << "blocking_factor = " << blockFactor << endl;
    pout() << "max_grid_size = " << maxGridSize << endl;
    pout() << "fill_ratio = " << fillRatio << endl;

    pout() << "normal_predictor = ";
    if (normalPredOrder == 0)
    {
      pout() << "CTU" << endl;
    }
    else if (normalPredOrder == 1)
    {
      pout() << "PLM" << endl;
    }
    else if (normalPredOrder == 2)
    {
      pout() << "PPM" << endl;
    }
    else
    {
      pout() << "Unknown (" << normalPredOrder << ")" << endl;
    }

    pout() << "slope_order = "
           << (useFourthOrderSlopes ? "2nd" : "4th") << endl;
    pout() << "use_primitive_slope_limiting = "
           << (usePrimLimiting ? "yes" : "no") << endl;
    pout() << "use_characteristic_slope_limiting = "
           << (useCharLimiting ? "yes" : "no") << endl;
    pout() << "use_slope_flattening = "
           << (useFlattening ? "yes" : "no") << endl;

    pout() << "use_artificial_viscosity = "
           << (useArtificialViscosity ? "yes" : "no") << endl;
    if (useArtificialViscosity)
    {
      pout() << "artificial_viscosity = " << artificialViscosity << endl;
    }

    pout() << "use_source_term = "
           << (useSourceTerm ? "yes" : "no") << endl;
    if (useSourceTerm)
    {
      pout() << "source_term_scaling = " << sourceTermScaling << endl;
    }

    pout() << "checkpoint_interval = " << checkpointInterval << endl;
    pout() << "plot_interval = " << plotInterval << endl;

    pout() << "CFL = " << cfl << endl;
    pout() << "initial_CFL = " << initialCFL << endl;

    pout() << "maximum_dt_growth = " << maxDtGrowth << endl;
    pout() << "dt_tolerance_factor = " << dtToleranceFactor << endl;
  }

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  // Set up the physics for polytropic gas dynamics
  PolytropicPhysics polytropicPhysics(smallPressure);
  polytropicPhysics.setPhysIBC(ibc);

  // Cast to physics base class pointer for technical reasons
  GodunovPhysics* godunovPhysics = static_cast<GodunovPhysics*> (&polytropicPhysics);

  // Set up the AMRLevel... factory
  AMRLevelPolytropicGasFactory amrGodFact;

  amrGodFact.define(cfl,
                    domainLength,
                    verbosity,
                    refineThresh,
                    tagBufferSize,
                    initialCFL,
                    godunovPhysics,
                    normalPredOrder,
                    useFourthOrderSlopes,
                    usePrimLimiting,
                    useCharLimiting,
                    useFlattening,
                    useArtificialViscosity,
                    artificialViscosity,
                    useSourceTerm,
                    sourceTermScaling,
                    highOrderLimiter);

  AMR amr;

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrGodFact);

  if (fixedDt > 0)
  {
    amr.fixedDt(fixedDt);
  }

  // Set grid generation parameters
  amr.maxGridSize(maxGridSize);
  amr.blockFactor(blockFactor);
  amr.fillRatio(fillRatio);

  // The hyperbolic codes use a grid buffer of 1
  amr.gridBufferSize(1);

  // Set output parameters
  amr.checkpointInterval(checkpointInterval);
  amr.plotInterval(plotInterval);
  amr.regridIntervals(regridIntervals);
  amr.maxDtGrow(maxDtGrowth);
  amr.dtToleranceFactor(dtToleranceFactor);

  // Set up output files
  if (ppgodunov.contains("plot_prefix"))
  {
    std::string prefix;
    ppgodunov.query("plot_prefix",prefix);
    amr.plotPrefix(prefix);
  }

  if (ppgodunov.contains("chk_prefix"))
  {
    std::string prefix;
    ppgodunov.query("chk_prefix",prefix);
    amr.checkpointPrefix(prefix);
  }

  amr.verbosity(verbosity);

  // Set up input files
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
      std::string gridFile;
      ppgodunov.query("fixed_hierarchy",gridFile);

      // initialize from a list of grids in "gridFile"
      Vector<Vector<Box> > amrGrids(maxLevel+1);
      setupFixedGrids(amrGrids,
                      probDomain,
                      maxLevel,
                      maxGridSize,
                      blockFactor,
                      verbosity,
                      gridFile);
      amr.setupForFixedHierarchyRun(amrGrids,1);
    }
  }
  else
  {
    std::string restartFile;
    ppgodunov.query("restart_file",restartFile);

#ifdef CH_USE_HDF5
    HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
    // read from checkpoint file
    amr.setupForRestart(handle);
    handle.close();
#else
    MayDay::Error("amrGodunov restart only defined with hdf5");
#endif
  }

  // Run the computation.
  amr.run(stopTime, nstop);

  // Output the last plot file and statistics.
  amr.conclude();
}

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile)
{
  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
  {
    a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));

    // Read in predefined grids
    ifstream is(a_gridFile.c_str(), ios::in);

    if (is.fail())
    {
      MayDay::Error("Cannot open grids file");
    }

    // Format of file:
    //   number of levels, then for each level (starting with level 1):
    //   number of grids on level, list of boxes

    int inNumLevels;
    is >> inNumLevels;

    CH_assert (inNumLevels <= a_maxLevel+1);

    if (a_verbosity >= 3)
    {
      pout() << "numLevels = " << inNumLevels << endl;
    }

    while (is.get() != '\n');

    a_amrGrids.resize(inNumLevels);

    // Check to see if coarsest level needs to be broken up
    domainSplit(a_domain,a_amrGrids[0],a_maxGridSize,a_blockFactor);

    if (a_verbosity >= 3)
    {
      pout() << "level 0: ";
      for (int n = 0; n < a_amrGrids[0].size(); n++)
      {
        pout() << a_amrGrids[0][0] << endl;
      }
    }

    // Now loop over levels, starting with level 1
    int ngrid;
    for (int lev = 1; lev < inNumLevels; lev++)
    {
      is >> ngrid;

      if (a_verbosity >= 3)
      {
        pout() << "level " << lev << " numGrids = " << ngrid << endl;
        pout() << "Grids: ";
      }

      while (is.get() != '\n');

      a_amrGrids[lev].resize(ngrid);

      for (int i = 0; i < ngrid; i++)
      {
        Box bx;
        is >> bx;

        while (is.get() != '\n');

        // Quick check on box size
        Box bxRef(bx);

        if (bxRef.longside() > a_maxGridSize)
        {
          pout() << "Grid " << bx << " too large" << endl;
          MayDay::Error();
        }

        if (a_verbosity >= 3)
        {
          pout() << bx << endl;
        }

        a_amrGrids[lev][i] = bx;
      } // End loop over boxes on this level
    } // End loop over levels
  }

  // Broadcast results to all the processors
  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
}

