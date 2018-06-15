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
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>

#ifndef CH_DISABLE_SIGNALS
#include <signal.h>  // For handling Ctrl-C.
#endif

#include "IntVectSet.H"
#include "parstream.H"
#include "Tuple.H"
#include "BoxIterator.H"
#include "AMR.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_TIMER
using namespace Chombo;
#endif

#define SMALL_TIME 1.0e-8

//-----------------------------------------------------------------------
// RUDIMENTARY SIGNAL HANDLING -- when control-C is pressed, what happens?
// If we're taking timesteps, we don't want to drop everything as soon as
// we hear a Ctrl-C. This signal handler allows us to set a flag to determine
// whether someone has pressed Ctrl-C and then act on it when we please.
// For the record, I'm putting this in here in order to be able to break out
// of a run within a Python interpreter. Python sets up its own SIGINT handler,
// so this enables us to get it back and treat it safely in the context of
// AMR calculations. -Jeff Johnson, 5/19/2010.

// This flag is set to true if an interrupt signal was intercepted, false otherwise.
static bool s_interrupted = false;

#ifndef CH_DISABLE_SIGNALS
//-----------------------------------------------------------------------
static void
handleCtrlC(int signum)
{
  s_interrupted = true;
}
//-----------------------------------------------------------------------
#endif

int AMR::s_step = 0;

//-----------------------------------------------------------------------
void AMR::useSubcyclingInTime(bool a_useSubcycling)
{
  m_useSubcycling = a_useSubcycling;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::setDefaultValues()
{
  m_isDefined = false;
  m_isSetUp   = false;
  m_useSubcycling = true;
  m_max_level = -1;
  m_finest_level = -1;
  m_checkpoint_interval= -1;
  m_plot_interval=-1;
  m_checkForSteadyState = false;
  m_plot_period=-1.0;
  m_next_plot_time=-1.0;
  m_max_grid_size= 0;
  m_max_base_grid_size = m_max_grid_size;
  m_restart_step = -1;
  m_lastcheck_step = -1;
  m_cur_step =0;
  s_step = m_cur_step;
  m_maxDtGrow = 1.1;
  m_time_eps = 1.0e-6;
  m_dt_base = -1;
  m_amrlevels.resize(0);
  m_use_meshrefine=false;
  m_plotfile_prefix = string("pltstate");
  m_checkpointfile_prefix = string("chk");
  m_verbosity = 0;
  m_cur_time = 0;
  m_dt_tolerance_factor = 1.1;
  m_fixedDt = -1;
  m_blockFactor = 4;
#ifdef CH_USE_TIMER
  m_timer = NULL ;
#endif
  m_scheduler = RefCountedPtr<Scheduler>();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AMR::AMR()
{
  setDefaultValues();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::clearMemory()
{
  CH_TIME("AMR::clearMemory");

  // clean up memory
  if (m_isDefined)
    {
      for (int lev = 0; lev < m_amrlevels.size(); ++lev)
        {
          if (m_amrlevels[lev] != NULL)
            {
              delete m_amrlevels[lev];
              m_amrlevels[lev] = NULL;
            }
        }
      m_isDefined = false;
      m_isSetUp = false;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AMR::~AMR()
{
  clearMemory();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::maxGridSize(int a_max_grid_size)
{
  CH_TIME("AMR::maxGridSize");

  CH_assert(a_max_grid_size >= 0);

  if (m_max_base_grid_size == 0)
    {
      m_max_base_grid_size = a_max_grid_size;
    }

  m_max_grid_size = a_max_grid_size;
  m_mesh_refine_ptr->maxSize(a_max_grid_size);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::maxBaseGridSize(int a_max_base_grid_size)
{
  CH_assert(a_max_base_grid_size > 0);

  m_max_base_grid_size = a_max_base_grid_size;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::dtToleranceFactor(Real a_dt_tolerance_factor)
{
  CH_assert(a_dt_tolerance_factor > 0);

  m_dt_tolerance_factor = a_dt_tolerance_factor;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::gridBufferSize(int a_grid_buffer_size)
{
  CH_TIME("AMR::gridBufferSize");

  m_mesh_refine_ptr->bufferSize(a_grid_buffer_size);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int AMR::maxGridSize() const
{
  CH_assert(isDefined());

  return(m_max_grid_size);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int AMR::maxBaseGridSize() const
{
  CH_assert(isDefined());
  return(m_max_base_grid_size);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::plotPrefix(const std::string& a_plotfile_prefix)
{
  CH_assert(isDefined());

  m_plotfile_prefix = a_plotfile_prefix;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::blockFactor(int a_blockFactor)
{
  CH_TIME("AMR::blockFactor");

  CH_assert(a_blockFactor >= 1);

  m_blockFactor = a_blockFactor;
  m_mesh_refine_ptr->blockFactor(a_blockFactor);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::setMeshRefine(RefCountedPtr<MeshRefine> a_mesh_refine_ptr)
{
    CH_assert(!a_mesh_refine_ptr.isNull());
    m_mesh_refine_ptr = a_mesh_refine_ptr;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::fillRatio(Real a_fillRatio)
{
  CH_TIME("AMR::fillRatio");

  CH_assert(isDefined());

  m_fillRatio = a_fillRatio;
  m_mesh_refine_ptr->fillRatio(a_fillRatio);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::checkpointPrefix(const std::string& a_checkpointfile_prefix)
{
  CH_assert(isDefined());

  m_checkpointfile_prefix = a_checkpointfile_prefix;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::regridIntervals(const Vector<int>& a_regridIntervals)
{
  CH_assert(isDefined());
  CH_assert(a_regridIntervals.size() >= m_max_level);

  m_regrid_intervals = a_regridIntervals;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::plotInterval(int a_plot_interval)
{
  m_plot_interval = a_plot_interval;
}

//-----------------------------------------------------------------------
void AMR::checkForSteadyState(bool a_checkForSteadyState)
{
  m_checkForSteadyState = a_checkForSteadyState;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::plotPeriod(Real a_plot_period)
{
  m_plot_period = a_plot_period;
  m_next_plot_time = m_cur_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::checkpointInterval(int a_checkpoint_interval)
{
  m_checkpoint_interval = a_checkpoint_interval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::define(int                          a_max_level,
                 const Vector<int>&           a_ref_ratios,
                 const Box&                   a_prob_domain,
                 const AMRLevelFactory* const a_amrLevelFact)
{
  ProblemDomain physdomain(a_prob_domain);

  define(a_max_level, a_ref_ratios, physdomain, a_amrLevelFact);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::define(int                          a_max_level,
                 const Vector<int>&           a_ref_ratios,
                 const ProblemDomain&         a_prob_domain,
                 const AMRLevelFactory* const a_amrLevelFact)
{
  CH_TIME("AMR::define");

  clearMemory();

  // Function already called in AMR() ... so make sure you set any values
  // _after_ the define() function. (ndk)
  setDefaultValues();

  m_isDefined = true;
  CH_assert(a_ref_ratios.size() > a_max_level);

  if (m_verbosity >= 3)
    {
      pout() << "AMR::define(...)" << endl;
    }

  CH_assert(a_max_level >= 0);

  m_max_level = a_max_level;

  // set to zero to start
  m_finest_level=0;
  m_use_meshrefine = true;

  // import refinement ratios
  m_ref_ratios = a_ref_ratios;
  if (m_ref_ratios.size() == 0)
    {
      m_ref_ratios.resize(1,1);
    }

  // Set the subcycling reduction factors to 1; i.e., no subcycling yet.
  m_reduction_factor.resize(a_ref_ratios.size(),1);
  if (m_reduction_factor.size() == 0)
    {
      m_reduction_factor.resize(1,1);
    }

  // resize vectors
  m_dt_new.resize(m_max_level+1);
  m_dt_cur.resize(m_max_level+1);
  m_steps_since_regrid.resize(m_max_level+1,0);
  m_regrid_intervals.resize(m_max_level+1,0);
  m_cell_updates.resize(m_max_level+1,0);

  // Create hierarchy of levels.
  m_amrlevels.resize(m_max_level+1,NULL);

  // Create base level.
  m_amrlevels[0] = a_amrLevelFact->new_amrlevel();

  m_amrlevels[0]->define(NULL,
                         a_prob_domain,
                         0,
                         m_ref_ratios[0]);

  // Create finer levels.
  int ref_factor = 1;
  for (int level = 0; level < m_max_level; ++level)
    {
      ref_factor *= m_ref_ratios[level];

      const ProblemDomain level_prob_domain = refine(a_prob_domain, ref_factor);

      m_amrlevels[level+1] = a_amrLevelFact->new_amrlevel();
      m_amrlevels[level+1]->define(m_amrlevels[level],
                                   level_prob_domain,
                                   level+1,
                                   m_ref_ratios[level+1]);
      m_amrlevels[level]->finerLevelPtr(m_amrlevels[level+1]);
    }

  // Build MeshRefine object using some default values.
  // These can be reset later, if desired.
  Real defaultFillRatio = 0.75;
  int defaultBlockFactor = 1;
  int defaultBufferSize = 1;

  // If we don't have a MeshRefine instance
  if (m_mesh_refine_ptr.isNull())
    {
      // Default implementation is BRMeshRefine
      m_mesh_refine_ptr = RefCountedPtr<MeshRefine>(new BRMeshRefine());
    }
  m_mesh_refine_ptr->define(a_prob_domain, a_ref_ratios, defaultFillRatio,
                       defaultBlockFactor, defaultBufferSize,m_max_grid_size);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::setupForFixedHierarchyRun(const Vector<Vector<Box> >& a_amr_grids,
                                    int                         a_proper_nest)
{
  CH_TIME("AMR::setupForFixedHierarchyRun");

  CH_assert(isDefined());

  m_isSetUp = true;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::setupForFixedHierarchyRun(...)" << endl;
    }

  m_finest_level = a_amr_grids.size() -1;
  m_finest_level_old = m_finest_level;

  // set to zero to start
  m_finest_level=0;

  m_use_meshrefine = false;

  // set this to -1 initially
  m_dt_base = -1;

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  m_amr_grids = a_amr_grids;
  m_regrid_intervals.resize(m_max_level, -1);
  m_cur_step = 0;
  s_step = m_cur_step;
  m_finest_level = a_amr_grids.size() - 1;

  if (m_finest_level > m_max_level)
    {
      cerr << "AMR::setupForFixedHierarchy: too many levels in input grids" << endl;
      cerr << "max_level from define function = " << m_max_level << endl;;
      cerr << "finest level(top level of grids in setup) = = " << m_finest_level << endl;
      MayDay::Error("AMR::setupForFixedHierarchy: define call and setup call inconsistent in class AMR");
    }

  {
    CH_TIME("Initialize");
    for (int level = 0; level <= m_finest_level; ++level)
      {
        m_amrlevels[level]->initialGrid(m_amr_grids[level]);
      }
    for (int level = m_finest_level; level >= 0; --level)
      {
        m_amrlevels[level]->postInitialGrid(false);
      }
    for (int level = 0; level <= m_finest_level; ++level)
      {
        m_amrlevels[level]->initialData();
      }
  }

  {
    CH_TIME("Regrid");
    for (int level = m_finest_level + 1; level <= m_max_level; ++level)
      {
        m_amrlevels[level]->regrid(Vector<Box>());
      }
  }

  {
    CH_TIME("PostInitialize");
    // call post-initialize once all the levels have been defined
    for (int level = m_finest_level; level >= 0; --level)
      {
        m_amrlevels[level]->postInitialize();
      }
  }

  {
    CH_TIME("InitialDt");
    for (int level = 0; level <= m_finest_level; ++level)
      {
        m_dt_new[level] = m_amrlevels[level]->computeInitialDt();
        m_dt_cur[level] = m_dt_new[level];
      }

    assignDt();
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::setupForNewAMRRun()
{
  CH_TIME("AMR::setupForNewAMRRun");

  CH_assert(m_isDefined);

  m_isSetUp = true;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::setupForNewAMRRun" << endl;
    }

  m_cur_step = 0;
  s_step = m_cur_step;
  initialGrid();

  m_finest_level_old = m_finest_level;
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_dt_new[level] = m_amrlevels[level]->computeInitialDt();
      m_dt_cur[level] = m_dt_new[level];
    }

  assignDt();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
#ifdef CH_USE_HDF5
// read checkpoint file
void AMR::setupForRestart(HDF5Handle& a_handle)
{
  CH_TIME("AMR::setupForRestart");

  CH_assert(m_isDefined);

  m_isSetUp = true;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::restart" << endl;
    }

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (m_verbosity >= 3)
    {
      pout() << "hdf5 header data: " << endl;
      pout() << header << endl;
    }

  // read max level
  if (header.m_int.find("max_level") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain max_level");
    }

  // note that this check should result in a warning rather than an error,
  // because you should be able to restart with a different number of levels
  // (DFM 2/27/02)
  int max_level_check = header.m_int ["max_level"];
  if (max_level_check != m_max_level)
    {
      pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
      pout() << "max level input to define = " << m_max_level << endl;
      pout() << "max level in checkpoint = " << max_level_check << endl;
      MayDay::Warning("AMR::restart: checkpoint file inconsistent with inputs to define ");
    }

  if (m_verbosity >= 2)
    {
      pout() << "read max_level = " << m_max_level << endl;
    }

  // read finest level
  if (header.m_int.find("num_levels") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain num_levels");
    }
  int num_levels = header.m_int ["num_levels"];

  if (m_verbosity >= 2)
    {
      pout() << "read num_levels = " << num_levels << endl;
    }

  m_finest_level = num_levels - 1;
  m_finest_level_old = m_finest_level;

  if (m_finest_level > m_max_level)
    {
      pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
      pout() << "numlevels input to define = " << m_max_level + 1<< endl;
      pout() << "numlevels in checkpoint = " << num_levels << endl;
      MayDay::Error("AMR::restart: checkpoint file inconsistent with inputs to define ");
    }

  if (m_verbosity >= 2)
    {
      pout() << "set finest_level = " << m_finest_level << endl;
    }

  if (m_finest_level > m_max_level)
    {
      MayDay::Error("AMR::restart: finest_level > max_level");
    }

  if (header.m_int.find("iteration") == header.m_int.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain iteration");
    }

  m_cur_step = header.m_int ["iteration"];
  s_step = m_cur_step;
  if (m_verbosity >= 2)
    {
      pout() << "read cur_step = " << m_cur_step << endl;
    }
  m_restart_step = m_cur_step;
  if (m_verbosity >= 2)
    {
      pout() << "set restart_step = " << m_restart_step << endl;
    }

  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMR::restart: checkpoint file does not contain time");
    }
  m_cur_time = header.m_real["time"];
  if (m_verbosity >= 2)
    {
      pout() << "read cur_time = " << m_cur_time << endl;
    }

  // number of steps since regrid: default 0 if not in file
  m_steps_since_regrid.resize(m_max_level+1, 0);
  for (int level = 0; level < m_max_level; ++level)
    {
      char headername[100];
      sprintf(headername, "steps_since_regrid_%d", level);
      if (header.m_int.find(headername) != header.m_int.end())
        {
          m_steps_since_regrid[level] = header.m_int [headername];
        }
    }
  if (m_verbosity >= 2)
    {
      pout() << "read steps_since_regrid =";
      for (int level = 0; level < m_max_level; ++level)
        {
          pout() << " " << m_steps_since_regrid[level];
        }
      pout() << endl;
    }

  // read physics class header data
  for (int level = 0; level <= m_finest_level; ++level)
    {
      // reset to root group to read header info
      a_handle.setGroup("/");
      m_amrlevels[level]->readCheckpointHeader(a_handle);
      m_amrlevels[level]->readCheckpointLevel(a_handle);
    }

  for (int level = 0; level < m_finest_level; ++level)
    {
      int refratio_test = m_amrlevels[level]->refRatio();
      if (refratio_test != m_ref_ratios[level])
        {
          pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
          pout() << "for level " << level << endl;
          pout() << "refratio input to define = " << m_ref_ratios[level] << endl;
          pout() << "refratio in checkpoint = "  << refratio_test << endl;
          MayDay::Error("AMR::restart: checkpoint file inconsistent with inputs to define ");
        }
    }

  // fine to coarse transversal of levels for grid setup
  for (int level = m_finest_level; level >= 0; --level)
    {
      m_amrlevels[level]->postInitialGrid(true);
    }

  // maintain time steps
  m_dt_new.resize(m_max_level+1);
  m_dt_cur.resize(m_max_level+1);

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_dt_new[level] = m_amrlevels[level]->dt();
      m_dt_cur[level] = m_dt_new[level];
    }

  assignDt();

  //restart cell updates(we could also output them to the chk file)
  m_cell_updates.resize(m_max_level+1,0);

  // final thing to do -- call initialGrid and initialData on undefined levels
  // (just in case there are setup things which need to be done there
  // (DFM 2/27/02)
  for (int level = m_finest_level+1; level <= m_max_level; ++level)
    {
      m_amrlevels[level]->initialGrid(Vector<Box>());
      m_amrlevels[level]->initialData();
    }
}
#endif
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::conclude() const
{
  CH_TIME("AMR::conclude");

  CH_assert(isDefined());
  CH_assert(isSetUp());

  if (m_verbosity >= 3)
    {
      pout() << "AMR::conclude" << endl;
    }

  if ((m_plot_interval >= 0) || (m_plot_period >= 0.0))
  {
    writePlotFile();
  }

  if ((m_checkpoint_interval >= 0)     &&
      (m_lastcheck_step != m_cur_step) &&
      (m_restart_step != m_cur_step))
    {
      writeCheckpointFile();
    }

  // Call any scheduled functions. This is placed after plotting and
  // checkpointing so that the plotting functions can congeal plot files.
  if (!m_scheduler.isNull())
  {
    m_scheduler->conclude(m_cur_step, m_cur_time);
  }

  // Call the hooks for the AMRLevel objects.
  for (int level = 0; level <= m_finest_level; ++level)
    m_amrlevels[level]->conclude(m_cur_step);

  if (m_verbosity >= 2)
    {
      long long total_cell_updates = 0;
      for (int ll = 0; ll < m_max_level+1; ll++)
        {
          total_cell_updates += m_cell_updates[ll];

#ifdef CH_OSF1
          pout() << "number of points updated at level " << ll <<" = " <<
              (Real)m_cell_updates[ll] << endl;
#else
          pout() << "number of points updated at level " << ll <<" = " <<
              m_cell_updates[ll] << endl;
#endif
        }
#ifdef CH_OSF1
      pout() << "total number of points updated = " << (Real)total_cell_updates << endl;
#else
      pout() << "total number of points updated = " << total_cell_updates << endl;
#endif
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// go baby go
void AMR::run(Real a_max_time, int a_max_step)
{
  CH_TIME("AMR::run");

  CH_assert(isDefined());
  CH_assert(isSetUp());

#ifdef CH_USE_TIMER
  double last_timestep_time = 0 ;
#endif
  if (m_verbosity >= 3)
    {
      pout() << "AMR::coarseTimeStep:" << endl;
      pout() << "max_time = " << a_max_time << endl;
      pout() << "max_step = " << a_max_step << endl;
    }

#ifndef CH_DISABLE_SIGNALS
  // Replace the current signal handler for Ctrl-C for the duration of
  // this method.
  struct sigaction ctrlC, oldCtrlC;
  ctrlC.sa_handler = handleCtrlC;
  sigemptyset(&ctrlC.sa_mask);
  ctrlC.sa_flags = 0;
  sigaction(SIGINT, &ctrlC, &oldCtrlC);
#endif

  Real old_dt_base = m_dt_base;
  bool steadyStateStop = false;
  for ( ; (m_cur_step < a_max_step) &&
          (a_max_time - m_cur_time > m_time_eps*m_dt_base);
        ++m_cur_step, m_cur_time += old_dt_base)
    {
        s_step = m_cur_step;
#ifdef CH_USE_TIMER
      m_timer->start() ;
#endif
      old_dt_base = m_dt_base;
      for (int level = 0; level <= m_max_level; ++level)
        {
          m_amrlevels[level]->time(m_cur_time);
        }

      // Drop a checkpoint file if it's time.
      if ((m_checkpoint_interval > 0)      &&
          (m_lastcheck_step != m_cur_step) &&
          (m_restart_step != m_cur_step)   &&
          (m_cur_step % m_checkpoint_interval == 0))
        {
          writeCheckpointFile();
          m_lastcheck_step= m_cur_step;
        }

      // Plot if we've gone enough steps.
      if ((m_plot_interval > 0) &&
          (m_cur_step % m_plot_interval == 0))
        {
          writePlotFile();
        }

      // Plot if enough time has passed.
      if ((m_plot_period > 0.0) &&
          (m_cur_time >= m_next_plot_time))
        {
          writePlotFile();
          m_next_plot_time = m_cur_time + m_plot_period;
        }

      // Call any scheduled functions. This is placed here so that
      // the plotter function can assume plot files have already been dumped.
      if (!m_scheduler.isNull())
        m_scheduler->execute(m_cur_step, m_cur_time);

      int level = 0;
      int stepsLeft = 0;
      bool timeBoundary = true;
      (void)timeStep(level,stepsLeft,timeBoundary);

      assignDt();
#ifdef CH_USE_TIMER
      m_timer->stop();
#endif
      if (m_verbosity >= 1)
        {
          std::ios::fmtflags origFlags = pout().flags();
          int origWidth = pout().width();
          int origPrecision = pout().precision();
          pout() << resetiosflags(ios::fixed)
                 << "coarse time step " << setw(3) << m_cur_step
                 << "  old time = " << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << m_cur_time
                 << "  old dt = "   << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << old_dt_base
#ifdef CH_USE_TIMER
                 << "  wallclocktime = "
                 << m_timer->wc_time() - last_timestep_time
#endif
                 << resetiosflags(ios::scientific)
                 << endl;
          pout().flags(origFlags);
          pout().width(origWidth);
          pout().precision(origPrecision);
        }
#ifdef CH_USE_TIMER
      last_timestep_time = m_timer->wc_time() ;
#endif

      // If we have assigned a signal handler for interrupts, check for
      // an interrupt and call the handler.
      if (s_interrupted)
        break;
      if (m_checkForSteadyState)
        {
          steadyStateStop = m_amrlevels[0]->convergedToSteadyState();
          for (int level = 1; level <= m_finest_level; ++level)
            {
              steadyStateStop = ((steadyStateStop)&&(m_amrlevels[level]->convergedToSteadyState()));
            }
        }
      if (steadyStateStop)
        {
          pout() << "AMR: all amr levels said steady state was reached" << endl;
          if (m_plot_interval > 0)
            {
              pout() << "writing plot file" << endl;
              writePlotFile();
            }
          if ((m_checkpoint_interval > 0)      &&
              (m_lastcheck_step != m_cur_step) )
            {
              pout() << "writing checkpoint file" << endl;
              writeCheckpointFile();
            }
          break;
        }
    }

#ifndef CH_DISABLE_SIGNALS
  // Re-instate the old Ctrl-C handler.
  sigaction(SIGINT, &oldCtrlC, NULL);

  // If Ctrl-C was uttered, clear the interrupt flag and notify the old
  // handler.
  if (s_interrupted)
  {
    s_interrupted = false;
    oldCtrlC.sa_handler(SIGINT);
  }
#endif
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// add restriction on growth of new time step
// dtnew = min(dtold*factor, dtwant) default factor = 1.1
void AMR::assignDt()
{
  CH_TIME("AMR::assignDt");

  CH_assert(isDefined());

  if (m_verbosity >= 3)
    {
      pout() << "AMR::assignDt" << endl;
    }

  if (m_useSubcycling)
    {
      if (m_fixedDt > 0)
        {
          m_dt_base = m_fixedDt;
        }
      else
        {
          m_dt_base = m_dt_new[0];

          // Multiply time step for each level by refinement factor, then min
          // across all levels.  Only go to finest_level_old because the higher
          // level dt's have not been set
          for (int level = 1; level <= m_finest_level_old; ++level)
            {
              int ref_factor = 1;

              for (int ilev = 0; ilev < level; ++ilev)
                {
                  ref_factor *= m_ref_ratios[ilev];
                }

              Real dt_base_equiv = m_dt_new[level] * ref_factor;
              m_dt_base = Min(m_dt_base, dt_base_equiv);
            }

          // Check all the actual dt's (scaled by the refinement factor) and only
          // allow a growth of "m_maxDtGrow".
          for (int level = 0; level <= m_finest_level_old; ++level)
            {
              int ref_factor = 1;

              for (int ilev = 0; ilev < level; ++ilev)
                {
                  ref_factor *= m_ref_ratios[ilev];
                }

              Real dt_base_equiv = m_dt_cur[level] * ref_factor * m_maxDtGrow;
              m_dt_base = Min(m_dt_base, dt_base_equiv);
            }
        }

      // refine base time step for all levels
      for (int level = 0; level <= m_max_level; ++level)
        {
          int ref_factor = 1;

          // reset reduction factors as there is no subcycling going on yet
          m_reduction_factor[level] = 1;

          for (int ilev = 0; ilev < level; ++ilev)
            {
              ref_factor *= m_ref_ratios[ilev];
            }

          Real dt_level = m_dt_base / ref_factor;
          m_amrlevels[level]->dt(dt_level);
        }
    }
  else
    {
      ///no subcycling.

      if (m_fixedDt > 0)
        {
          m_dt_base = m_fixedDt;
        }
      else
      {
        m_dt_base = m_dt_new[0];
        for (int level = 0; level <= m_finest_level_old; ++level)
          {
            m_dt_base = Min(m_dt_base, m_dt_new[level]);
            Real dt_base_equiv = m_dt_cur[level]*m_maxDtGrow;
            m_dt_base = Min(m_dt_base, dt_base_equiv);
          }
      }
      //everybody gets the same time step if subcyling is turned off
      for (int level = 0; level <= m_max_level; ++level)
        {
          m_amrlevels[level]->dt(m_dt_base);
        }
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// This function performs a time step of level "a_level" and all finer
// levels.  It does subcycling in time as necessary (see below).  For this
// reason the number of time steps left, "a_stepsLeft", for this level is
// passed in.  If no additional subcycling occurs then this is returned.
// If additional subcycling occurs then this is used to compute the new
// number of time steps remaining (for this level) and this is returned.
int AMR::timeStep(int a_level, int a_stepsLeft, bool a_coarseTimeBoundary)
{
  CH_TIME("AMR::timeStep");
  CH_assert(isDefined());
  CH_assert(isSetUp());
  if (m_verbosity >= 3)
    {
      pout() << "AMR::timeStep(" << a_level << ")" << endl;
    }

  if (m_useSubcycling)
    {
      if (a_level < m_max_level)
        {
          CH_TIME("AMR::timeStep::regrid");
          if (m_verbosity >= 4)
            {
              pout() << "Regrid (level " << a_level << ") needed - ";
            }

          // regrid if necessary
          if (needToRegrid(a_level,a_stepsLeft))
            {
              if (m_verbosity >= 4)
                {
                  pout() << "yes" << endl;
                }

              regrid(a_level);
            }
          else
            {
              if (m_verbosity >= 4)
                {
                  pout() << "no" << endl;
                }
            }
        }

      // decrement a_stepsLeft here
      // (need to do it after checking for regridding)
      a_stepsLeft--;

      // If this wasn't just done by the next coarser level, check to see if
      // it is necessary to do additional subcycling in time.
      if ((!a_coarseTimeBoundary) && (m_fixedDt <= 0))
        {
          CH_TIME("AMR::timeStep::subcycle");

          // The factor by which the current time step at the current level
          // has been divided (so far) for subcycling.
          int maxFactor = m_reduction_factor[a_level];

          // Compute the new subcycling factor for this level and all finer
          // levels and find the maximum
          for (int i = a_level; i <= m_max_level; i++)
            {
              int factor;
              Real dtCur = m_amrlevels[i]->dt();
              Real dtNew = m_dt_new[i];

              // The current factor for level "i"
              factor = m_reduction_factor[i];

              // While the current dt exceeds the new (max) dt by a tolerance
              // double the subcycling factor and half the current dt
              while (dtCur > m_dt_tolerance_factor*dtNew)
                {
                  factor *= 2;
                  dtCur *= 0.5;
                }

              if (factor > maxFactor)
                {
                  maxFactor = factor;
                }
            }

          // More subcycling is necessary
          if (maxFactor > m_reduction_factor[a_level])
            {
              if (m_verbosity >= 3)
                {
                  pout() << "  Subcycling --- maxFactor: " << maxFactor << endl;
                }

              // Adjust the number of time steps left for the current level
              a_stepsLeft = (a_stepsLeft+1)*maxFactor/m_reduction_factor[a_level] - 1;

              // Adjust the dt's on this and all finer levels
              for (int i = a_level; i <= m_max_level; i++)
                {
                  int factor;

                  factor = maxFactor/m_reduction_factor[i];
                  m_amrlevels[i]->dt(m_amrlevels[i]->dt()/factor);

                  if (m_verbosity >= 4)
                    {
                      pout() << "    Level " << i << ": factor: " << factor
                             << " (" << m_reduction_factor[i] << "), dt: "
                             << m_amrlevels[i]->dt() << endl;
                    }

                  m_reduction_factor[i] = maxFactor;
                }
            }
        }


      {
        // advance this level
        CH_TIME("AMR::timeStep::advance");
        m_amrlevels[a_level]->advance();
      }

      Real dt_level;

      {
        // get the new dt
        CH_TIME("AMR::timeStep::newDt");
        dt_level = m_amrlevels[a_level]->computeDt();
      }

      // Save the current dt and the new (max) dt.
      m_dt_cur[a_level] = m_amrlevels[a_level]->dt();
      m_dt_new[a_level] = dt_level;

      // increment counter that gives the number of cells updates.
      long long numPts = 0;
      const Vector<Box >& levelBoxes = m_amrlevels[a_level]->boxes();

      for (int ll = 0;ll < levelBoxes.size();ll++)
        {
          numPts += levelBoxes[ll].numPts();
        }

      m_cell_updates[a_level] += numPts;

      if (a_level < m_max_level)
        {
          ++m_steps_since_regrid[a_level];
        }

      // advance the finer levels by subcycling
      if (a_level < m_finest_level)
        {
          CH_TIME("AMR::timeStep::finerLevels");

          int stepsLeft = m_ref_ratios[a_level];
          bool timeBoundary = true;

          while (stepsLeft > 0)
            {
              // Advance the finer level and take into account possible
              // subcycling by allowing for a change in "stepsLeft".
              //[NOTE: the if () test looks redundant with above, but it is not
              //       because m_finest_level may change during a regrid();
              //       why you would regrid during a subcycle I don't know. <dbs>]
              if (a_level < m_finest_level)
                stepsLeft = timeStep(a_level+1,stepsLeft,timeBoundary);

              // The first time the next finer level time aligns with the current
              // level time.  After that this is not the case.
              //[NOTE: this if () test _is_ redundant. <dbs>]
              if (timeBoundary == true)
                {
                  timeBoundary = false;
                }
            }
        }

      {
        CH_TIME("AMR::timeStep::postTimeStep");
        m_amrlevels[a_level]->postTimeStep();
      }

    }
  else
    {
      //subcycling is off so things are quite a bit simpler
      if (a_level != 0)
        {
          MayDay::Error("AMR::timestep is confused by being called for a_level != 0 with subcycling off");
        }
      if (a_level < m_max_level)
        {
          CH_TIME("AMR::timeStep::regrid");
          if (m_verbosity >= 4)
            {
              pout() << "Regrid (level " << a_level << ") needed - ";
            }

          // regrid if necessary
          //if not, increment number of steps since regrid.
          bool regridThisLevel =
            (m_regrid_intervals[a_level] > 0) &&
            (m_steps_since_regrid[a_level] >= m_regrid_intervals[a_level]);
          if (regridThisLevel)
            {
              regrid(a_level);
              m_steps_since_regrid[a_level] = 0;
            }
          else
            {
              m_steps_since_regrid[a_level]++;
            }
        }
      {
        // advance the levels
        //consistent with subcycling, do in level order
        CH_TIME("AMR::timeStep::advance");
        for (int ilev = 0; ilev <= m_finest_level; ilev++)
          {
            const Vector<Box >& levelBoxes = m_amrlevels[ilev]->boxes();

            long long numPts = 0;
            for (int ll = 0;ll < levelBoxes.size();ll++)
              {
                numPts += levelBoxes[ll].numPts();
              }
            m_cell_updates[ilev] += numPts;

            m_amrlevels[ilev]->advance();
            m_dt_cur[ilev] = m_amrlevels[ilev]->dt();
          }
      }
      {
        // get the new dt and save the old one
        CH_TIME("AMR::timeStep::newDt");
        Real dtMin = m_amrlevels[0]->computeDt();
        for (int ilev = 1; ilev <= m_finest_level; ilev++)
          {
            Real dt_level = m_amrlevels[ilev]->computeDt();
            dtMin = Min(dt_level, dtMin);
          }
        for (int ilev = 0; ilev <= m_finest_level; ilev++)
          {
            m_dt_new[ilev] = dtMin;
          }
      }

      {
        // advance the levels
        CH_TIME("AMR::timeStep::postTimeStep");
        //call postTimeStep.
        //to stay consistent with the way things get done in subcycling land,
        //this gets done from the finest level first
        for (int ilev = m_finest_level; ilev >= 0; ilev--)
          {
            m_amrlevels[ilev]->postTimeStep();
          }
      }

    }

  // Return the (possibly updated) number of time steps left on this level.
  return(a_stepsLeft);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool AMR::needToRegrid(int a_level, int a_stepsLeft) const
{
  CH_TIME("AMR::needToRegrid");

  CH_assert(isDefined());
  CH_assert(isSetUp());

  bool regridThisLevel = false;
  bool regridNextCoarserLevel = false;

  if (m_verbosity >= 3)
    {
      pout() << "AMR::needToRegrid(" << a_level << ")" << endl;
    }

  regridThisLevel =
    (m_regrid_intervals[a_level] > 0) &&
    (m_steps_since_regrid[a_level] >= m_regrid_intervals[a_level]);

  int nextCoarserLevel = a_level-1;

  regridNextCoarserLevel =
    (a_stepsLeft == 0)                         &&
    (nextCoarserLevel >= 0)                    &&
    (m_regrid_intervals[nextCoarserLevel] > 0) &&
    (m_steps_since_regrid[nextCoarserLevel] >= m_regrid_intervals[nextCoarserLevel]);

  return(regridThisLevel && !regridNextCoarserLevel);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// generate new grid hierarchy
void AMR::regrid(int a_base_level)
{
  CH_TIME("AMR::regrid");

  CH_assert(isDefined());
  CH_assert(isSetUp());

  if (m_verbosity >= 2)
    {
      pout() << "AMR::regrid(" << a_base_level << ")" << endl;
    }

  if (m_verbosity >= 2)
    {
      pout() << "AMR::regrid: base level = " << a_base_level << endl;
    }

  for (int level = a_base_level; level <= m_max_level; ++level)
    {
      m_steps_since_regrid[level] = 0;
    }

  m_finest_level_old = m_finest_level;
  int top_level = Min(m_finest_level, m_max_level - 1);

  Vector<Vector<Box> > old_grids(top_level+1);
  Vector<Vector<Box> > new_grids;
  Vector<IntVectSet> tags(top_level+1);
  Vector<ProblemDomain> problem_domains(top_level+1);

  if (m_use_meshrefine)
    {
      for (int level = a_base_level; level <= top_level; ++level)
        {
          m_amrlevels[level]->tagCells(tags[level]);
          old_grids[level] = m_amrlevels[level]->boxes();
          problem_domains[level] = m_amrlevels[level]->problemDomain();

          if (m_verbosity >= 2)
            {
              pout() << "AMR::regrid: problem domain[" << level << "]: " << problem_domains[level] << endl;
            }

          if (m_verbosity >= 4)
            {
              pout() << "AMR::regrid: old_grids[" << level << "]: " << endl;

              for (int i = 0; i< old_grids[level].size(); ++i)
                {
                  pout() << "  " << i << ": " << old_grids[level][i] << endl;
                }

              if (m_verbosity >= 5)
                {
                  pout() << "AMR::regrid: tags[" << level << "]: " << tags[level]
                         << endl;
                }
            }
        }

      int new_finest_level;
      new_finest_level = m_mesh_refine_ptr->regrid(new_grids,
                                                   tags,
                                                   a_base_level,
                                                   top_level,
                                                   old_grids);

      //can only add one level at a time
      new_finest_level = Min(m_finest_level+1, new_finest_level);

      if ((m_finest_level != new_finest_level) && (m_verbosity >= 2))
        {
          pout() << "finest level changes here from "
                 << m_finest_level << " to "
                 << new_finest_level << endl;
        }

      //allow for levels to change
      m_finest_level = Min(new_finest_level, m_max_level);

      //need to assign times if number of levels has grown
      if (m_finest_level > m_finest_level_old)
        {
          Real ratTime = m_amrlevels[m_finest_level_old]->time();

          for (int level = m_finest_level_old+1; level <= m_finest_level; ++level)
            {
              m_amrlevels[level]->time(ratTime);
            }
        }

      if (m_verbosity >= 4)
        {
          if (new_grids.size() == 0)
            {
              pout() << "No new_grids" << endl;
              pout() << endl;
            }
          else
            {
              for (int level = a_base_level; level <= m_finest_level; ++level)
                {
                  pout() << "new_grids[" << level << "]: " << endl;

                  for (int i = 0; i< new_grids[level].size(); ++i)
                    {
                      pout() << "  " << i << ": " << new_grids[level][i] << endl;
                    }

                  pout() << endl;
                }
            }
        }

    }
  else
    {
      // use pre-defined grids
      new_grids = m_amr_grids;
    }

  // before regridding (but after tagging) allow for pre-regridding ops
  // (tjl 8/23/06 - this needed for mapped grid computations and a default
  // implementation is provided in AMRLevel to preserve existing codes)
  for (int level = m_finest_level; level >= a_base_level; --level)
    {
      m_amrlevels[level]->preRegrid(a_base_level, new_grids);
    }

  for (int level = a_base_level + 1; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->regrid(new_grids[level]);
    }

  for (int level = m_finest_level + 1; level <= m_max_level; ++level)
    {
      m_amrlevels[level]->regrid(Vector<Box>());
    }

  // now that the new hierarchy is defined, do post-regridding ops
  // (dfm 8/26/05 -- call postRegrid on base_level as well, to
  // cover the case where all levels finer than base_level are removed)
  for (int level = m_finest_level; level >= a_base_level; --level)
    {
      m_amrlevels[level]->postRegrid(a_base_level);
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<AMRLevel*> AMR::getAMRLevels()
{
  return m_amrlevels;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::initialTime(Real a_initialTime)
{
  // this needs to be called before setup, but after define
  CH_assert(m_isDefined);
  CH_assert(!m_isSetUp);

  m_cur_time = a_initialTime;

  // propagate this time to all of the AMR levels
  for (int lev=0; lev<m_amrlevels.size(); lev++)
    {
      m_amrlevels[lev]->time(a_initialTime);
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real AMR::getCurrentTime() const
{
  return m_cur_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::makeBaseLevelMesh(Vector<Box>& a_grids) const
{
  CH_TIME("AMR::makeBaseLevelMesh");

  CH_assert(m_isDefined);

  if (m_max_base_grid_size == 0)
    {
      // define base level to be single grid
      a_grids.resize(1);
      a_grids[0] = m_amrlevels[0]->problemDomain().domainBox();
    }
  else
    {
      if (m_max_base_grid_size < m_blockFactor)
        {
          MayDay::Abort("Base grid size must be greater than blocking factor");
        }

      // chop base level up into grids of no more than m_max_grid_size on a side
      // first coarsen to enforce the blocking factor
      ProblemDomain problem_domain = coarsen(m_amrlevels[0]->problemDomain(),
                                             m_blockFactor);
      if (refine(problem_domain, m_blockFactor) != m_amrlevels[0]->problemDomain())
        {
          MayDay::Error("level 0 problem domain not coarsenable by blocking factor");
        }
      int max_grid_size = m_max_base_grid_size / m_blockFactor;
      Tuple<Vector<int>,SpaceDim> box_sizes;
      IntVect num_grids;
      IntVect base_size;

      for (int d = 0; d < SpaceDim; ++d)
        {
          int num_div = 1;
          int domain_size = problem_domain.domainBox().size(d);
          while (num_div * max_grid_size < domain_size)
            {
              ++num_div;
            }

          // int(x/y) +(x%y)?1:0 is integer division with rounding upwards, for x,y>0
          base_size[d] = int(domain_size/num_div) +((domain_size%num_div) ? 1 : 0);
          box_sizes[d].resize(num_div, base_size[d]);
          box_sizes[d][num_div-1] = domain_size -(num_div - 1) * base_size[d];
          num_grids[d] = num_div;
        }

      Box b(IntVect::Zero,num_grids - IntVect::Unit);
      const IntVect& domain_hi = problem_domain.domainBox().bigEnd();
      const IntVect& domain_lo = problem_domain.domainBox().smallEnd();

      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          IntVect lo = domain_lo + iv * base_size;
          IntVect hi = min(lo + base_size - IntVect::Unit, domain_hi);
          Box grid(lo,hi);

          a_grids.push_back(refine(grid, m_blockFactor) );
        }
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// generate initial grid hierarchy
void AMR::initialGrid()
{
  CH_TIME("AMR::initialGrid");

  CH_assert(m_isDefined);

  if (m_verbosity >= 3)
    {
      pout() << "AMR::initialGrid" << endl;
    }

  Vector<Vector<Box> > old_grids(1);
  Vector<Vector<Box> > new_grids;

  makeBaseLevelMesh(old_grids[0]);

  // we will keep the old tags around while we generate new ones
  Vector<IntVectSet> old_tags(m_max_level);

  if (m_use_meshrefine)
    {
      for (int top_level = 0;
           top_level < Min(m_finest_level+1, m_max_level);
           ++top_level)
        {
          if (m_verbosity >= 2)
            {
              pout() << "AMR::initialGrid: top level = " << top_level << endl;
            }

          Vector<IntVectSet> tags(top_level+1);
          Vector<ProblemDomain> problem_domains(top_level+1);

          // loop over all levels and initialize grids,
          // then loop over all levels and initialize data
          // then loop over levels and generate tags
          // do this in three separate loops to handle
          // the case where initial data is generated
          // using a multilevel operation(for instance,
          // computing initial velocity from initial
          // vorticity through a multilevel elliptic solve)
          // DFM(11/28/2000)
          for (int level = 0; level <= top_level; ++level)
            {
              m_amrlevels[level]->initialGrid(old_grids[level]);
            }

          for (int level = top_level; level >= 0; --level)
            {
              m_amrlevels[level]->postInitialGrid(false);
            }

          for (int level = 0; level <= top_level; ++level)
            {
              m_amrlevels[level]->initialData();
            }

          for (int level = 0; level <= top_level; ++level)
            {
              m_amrlevels[level]->tagCellsInit(tags[level]);

              // union old tags with current ones.  this prevents
              // you from unrefining a region you previously
              // decided you wanted refined
              tags[level] |= old_tags[level];

              problem_domains[level] = m_amrlevels[level]->problemDomain();

              if (m_verbosity >= 3)
                {
                  pout() << "AMR::initialGrid: problem domain["
                         << level << "]: " << problem_domains[level] << endl;
                }

              if (m_verbosity >= 4)
                {
                  pout() << "AMR::initialGrid: old_grids[" << level << "]: " << endl;

                  for (int i = 0; i< old_grids[level].size(); ++i)
                    {
                      pout() << "  " << i << ": " << old_grids[level][i] << endl;
                    }

                  if (m_verbosity >= 5)
                    {
                      pout() << "AMR::initialGrid: tags[" << level << "]: "
                             << tags[level] << endl;
                    }
                }
            }
          m_finest_level = m_mesh_refine_ptr->regrid(new_grids,
                                                     tags,
                                                     0,
                                                     top_level,
                                                     old_grids);

          // do this only if a new level was generated
          if (m_finest_level > top_level)
            {
              old_grids = new_grids;

              // copy current tags to old_tags
              for (int lev = 0; lev <= top_level; lev++)
                {
                  old_tags[lev] = tags[lev];
                }
            }

          if (m_verbosity >= 4)
            {
              for (int level = 0; level < new_grids.size(); ++level)
                {
                  pout() << "new_grids[" << level << "]: " << endl;

                  for (int i = 0; i< new_grids[level].size(); ++i)
                    {
                      pout() << "  " << i << ": " << new_grids[level][i] << endl;
                    }

                  pout() << endl;
                }
            }
        }

      if (m_finest_level == 0)
        {
          new_grids.resize(1);
          new_grids[0] = old_grids[0];
        }
    }
  else
    {
      // use predefined grids
      CH_assert(m_amr_grids.size() > 0);

      new_grids = m_amr_grids;
      m_finest_level = m_amr_grids.size()-1;
    }

  // separate loops for initialGrid and initialData
  // to ensure that hierarchy is defined before we do initialization
  // (for case where initialization is a multilevel operatation)
  // (dfm 11/7/02)
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->initialGrid(new_grids[level]);
    }

  for (int level = m_finest_level; level >= 0; --level)
    {
      m_amrlevels[level]->postInitialGrid(false);
    }

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->initialData();
    }

  for (int level = m_finest_level + 1; level <= m_max_level; ++level)
    {
      // m_amrlevels[level]->regrid(Vector<Box>());
      m_amrlevels[level]->initialGrid(Vector<Box>());
      m_amrlevels[level]->initialData();
    }

  // call post-initialize once all the levels have been defined
  for (int level = m_finest_level; level >= 0; --level)
    {
      m_amrlevels[level]->postInitialize();
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::writePlotFile() const
{
  CH_TIME("AMR::writePlotFile");

  CH_assert(m_isDefined);

  if (m_verbosity >= 3)
    {
      pout() << "AMR::writePlotFile" << endl;
    }

#ifdef CH_USE_HDF5
  string iter_str = m_plotfile_prefix;

  char suffix[100];
  sprintf(suffix,"%06d.%dd.hdf5",m_cur_step,SpaceDim);

  iter_str += suffix;

  if (m_verbosity >= 2)
    {
      pout() << "plot file name = " << iter_str << endl;
    }

  HDF5Handle handle(iter_str.c_str(), HDF5Handle::CREATE);

  // write amr data
  HDF5HeaderData header;
  header.m_int ["max_level"]  = m_max_level;
  header.m_int ["num_levels"] = m_finest_level + 1;
  header.m_int ["iteration"]  = m_cur_step;
  header.m_real["time"]       = m_cur_time;

  // should steps since regrid be in the checkpoint file?
  header.writeToFile(handle);

  if (m_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // write physics class header data
  m_amrlevels[0]->writePlotHeader(handle);

  // write physics class per-level data
  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->writePlotLevel(handle);
    }

  handle.close();
#endif

  // Here's an extra hook for doing custom plots. :-P -JNJ
  m_amrlevels[0]->writeCustomPlotFile(m_plotfile_prefix, m_cur_step);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::writeCheckpointFile() const
{
  CH_TIME("AMR::writeCheckpointFile");

  CH_assert(m_isDefined);

  if (m_verbosity >= 3)
    {
      pout() << "AMR::writeCheckpointFile" << endl;
    }

#ifdef CH_USE_HDF5
  string iter_str = m_checkpointfile_prefix;

  char suffix[100];
  sprintf(suffix,"%06d.%dd.hdf5",m_cur_step,SpaceDim);

  iter_str += suffix;

  if (m_verbosity >= 2)
    {
      pout() << "checkpoint file name = " << iter_str << endl;
    }

  HDF5Handle handle(iter_str.c_str(), HDF5Handle::CREATE);

  // write amr data
  HDF5HeaderData header;
  header.m_int ["max_level"]  = m_max_level;
  header.m_int ["num_levels"] = m_finest_level + 1;
  header.m_int ["iteration"]  = m_cur_step;
  header.m_real["time"]       = m_cur_time;

  for (int level = 0; level < m_regrid_intervals.size(); ++level)
    {
      char headername[100];
      sprintf(headername, "regrid_interval_%d", level);
      header.m_int[headername] = m_regrid_intervals[level];
    }

  for (int level = 0; level < m_max_level; ++level)
    {
      char headername[100];
      sprintf(headername, "steps_since_regrid_%d", level);
      header.m_int[headername] = m_steps_since_regrid[level];
    }

  // should steps since regrid be in the checkpoint file?
  header.writeToFile(handle);

  if (m_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // write physics class data
  m_amrlevels[0]->writeCheckpointHeader(handle);

  for (int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->writeCheckpointLevel(handle);
    }

  handle.close();
#endif
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::verbosity(int a_verbosity)
{
  m_verbosity = a_verbosity;

  if (m_verbosity >= 2)
    {
      pout() << " AMR::verbosity  m_verbosity=" << m_verbosity << endl;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int AMR::verbosity() const
{
  return(m_verbosity);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::maxDtGrow(Real a_dtGrowFactor)
{
  m_maxDtGrow = a_dtGrowFactor;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real AMR::maxDtGrow() const
{
  return m_maxDtGrow;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::timeEps(Real a_timeEps)
{
  m_time_eps = a_timeEps;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real AMR::timeEps() const
{
  return m_time_eps;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void AMR::fixedDt(Real a_dt)
{
  m_fixedDt = a_dt;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real AMR::fixedDt() const
{
  return m_fixedDt;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool AMR::isDefined() const
{
  return m_isDefined;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool AMR::isSetUp() const
{
  return m_isSetUp;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
#ifdef CH_USE_TIMER
Timer * AMR::timer(Timer *a_timer )
{
  Timer * old_timer = m_timer ;
  if ( a_timer != NULL )
  {
    m_timer = a_timer ;
  }
  return old_timer ;
}
#endif
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AMR::
schedule(RefCountedPtr<Scheduler> a_scheduler)
{
  CH_assert(!a_scheduler.isNull());
  m_scheduler = a_scheduler;

  // Set up the scheduler's periodic functions to interact with
  // this AMR object.
  m_scheduler->setUp(*this);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
