#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMR_H_
#define _AMR_H_

#include <iostream>
#include <string>
#include <ctime>

#include "REAL.H"
#include "Vector.H"
#include "AMRLevel.H"
#include "AMRLevelFactory.H"
#include "BRMeshRefine.H"
#include "RefCountedPtr.H"
#include "ProblemDomain.H"
#include "Box.H"
#include "CH_HDF5.H"
#include "Scheduler.H"
#include "NamespaceHeader.H"

/// Framework for Berger-Oliger timestepping for AMR
/**
   This class is a framework for Berger-Oliger timestepping for
   adaptive mesh refinement of time-dependent problems.  It is
   applicable to both hyperbolic and parabolic problems.  It
   represents a hierarchy of levels of refinement as a collection of
   AMRLevel objects.

   The usage pattern is this:

   <ol>
     <li>Call define() to define the stuff that does not change throughout
         the run (maxlevel, refinement ratios, domain, and operator).
     <li>Modify any parameters you like (blocking factor, max grid size, etc.)
         using parameter-setting functions.
     <li>Call any one of the three setup functions (setupForRestart(),
         setupForNewAMRRun(), or setupForFixedHierarchyRun())
         so AMR can set up all its internal data structures.
     <li>Call run() to run the calculation.
     <li>Call conclude() to output (to stdout) how many cells were updated
         and so forth.
   </ol>
*/
class AMR
{

public:
  //hook to do output of stuff by timestep in objects that do not know AMR
  static int s_step;
  ///
  /**
     Default constructor.  Use must subsequently call a define() function
     and a setup function.
  */
  AMR();

  ///
  /**
     Destructor.
  */
  virtual ~AMR();

  ///
  /**
     Defines this object.  User must call a setup function before running.

     Arguments:

     <ul>
      <li>a_max_level (not modified): the maximum level allowed, where the
          base level is zero.  There will be a total of a_max_level+1
          levels, since level zero and level a_max_level will both
          exist.
      <li>a_ref_ratios (not modified): refinement ratios.  There must be at
          least a_max_level+1 elements,  or an error will result. Element
          zero is the base level.
      <li>a_prob_domain (not modified): problem domain on the base level.
      <li>a_amrlevelFact (not modified): pointer to a physics class factory
          object.  The object it points to is used to construct the
          collection of AMRLevel objects in this AMR as objects of the
          physics class type.  It can be destructed after this define()
          returns.
     </ul>
  */
  void define(int                          a_max_level,
              const Vector<int>&           a_ref_ratios,
              const Box&                   a_prob_domain,
              const AMRLevelFactory* const a_amrLevelFact);

  ///
  /**
     Defines this object.  User must call a setup function before running.

     Arguments:

     <ul>
      <li>a_max_level (not modified): the maximum level allowed, where the
          base level is zero.  There will be a total of a_max_level+1
          levels, since level zero and level a_max_level will both
          exist.
      <li>a_ref_ratios (not modified): refinement ratios.  There must be at
          least a_max_level+1 elements,  or an error will result. Element
          zero is the base level.
      <li>a_prob_domain (not modified): problem domain on the base level.
      <li>a_amrlevelFact (not modified): pointer to a physics class factory
          object.  The object it points to is used to construct the
          collection of AMRLevel objects in this AMR as objects of the
          physics class type.  It can be destructed after this define()
          returns.
     </ul>
  */
  void define(int                          a_max_level,
              const Vector<int>&           a_ref_ratios,
              const ProblemDomain&         a_prob_domain,
              const AMRLevelFactory* const a_amrLevelFact);

#ifdef CH_USE_HDF5
  ///
  /**
     Sets up this object from checkpointed data.  User must have
     previously called define().

     Need to call this function or setupForNewAMRRun() or
     setupforFixedHierarchyRun() before you run.
  */
  void setupForRestart(HDF5Handle& a_handle);
#endif

  ///
  /**
     Sets up this object for cold start.  User must have previously
     called define().

     Need to call this function or setupForRestart() or
     setupforFixedHierarchyRun() before you run.
   */
  void setupForNewAMRRun();

  /**
     This function sets the hierarchy and sets m_regrid_intervals to -1
     (turns off regridding).  If you want to keep regridding
     on, call regridIntervals() after this call.
   */
  void setupForFixedHierarchyRun(const Vector<Vector<Box> >& a_amr_grids,
                                 int                         a_proper_nest = 1);

  ///
  /**
     Runs the calculation.  User must have previously called
     both the define() function and a setup function.
  */
  void run(Real a_max_time, int a_max_step);

  ///
  /**
     You should call this last.  It writes the last
     checkpoint file and tells you how many cells you updated and
     all that.
   */
  void conclude() const;

  /**
     \name Parameter-setting functions
  */

  /**@{*/

  ///
  /**
     Sets the checkpoint file prefix.

     Should be called after define()
     and before setup.
  */
  void checkpointPrefix(const std::string& a_checkpointfile_prefix);

  //! Tells AMR to write plot files after every \a a_plot_interval steps.
  void plotInterval(int a_plot_interval);


  //! Tells AMR to check each level for steady state and stop if we get there
  void checkForSteadyState(bool a_steadyState);

  //! Tells AMR to write plot files after every \a a_plot_period time units.
  void plotPeriod(Real a_plot_period);

  //! Sets up a schedule for periodically-called functions.
  void schedule(RefCountedPtr<Scheduler> a_scheduler);

  ///
  /**
     Sets the interval to write checkpoint files, in terms of the base level
     time step.
  */
  void checkpointInterval(int a_checkpoint_interval);

  ///
  /**
     Set the maximum grid size.  Should be called after define()
     and before setup.
   */
  void maxGridSize(int a_max_grid_size);

  ///
  /**
     Set the maximum grid size for level 0 grids.

     Defaults to m_max_grid_size, so should be called after maxGridSize()
     and before setup.
   */
  void maxBaseGridSize(int a_max_base_grid_size);

  ///
  /**
     Set the factor by which the current dt must exceed the new (max) dt
     for time subcycling to occur (i.e., reduction of the current dt by
     powers of 2).
   */
  void dtToleranceFactor(Real a_dt_tolerance_factor);

  ///
  /**
     Set the MeshRefine instance.  Should be called before define()
     if you want a_mesh_refine_ptr defined.
   */
  void setMeshRefine(RefCountedPtr<MeshRefine> a_mesh_refine_ptr);

  ///
  /**
     Set the fill ratio for MeshRefine.  Should be called after define()
     and before setup.
   */
  void fillRatio(Real a_fillRat);

  ///
  /**
     Set the blocking factor for MeshRefine.  Should be called after define()
     and before setup.
   */
  void blockFactor(int a_blockFactor);

  ///
  /**
     Set the buffering for MeshRefine. Should be called after define()
     and before setup.
  */
  void gridBufferSize(int a_grid_buffer_size);

  ///
  /**
     Sets verbosity level to a_verbosity.

     -inf to 0: print nothing.<br>
     1: prints a message every coarse time step.<br>
     2: also prints a message when regridding and writing files.<br>
     3: also prints function trace, and some more inf during file I/O.<br>
     4 to inf: also prints list of boxes and processor maps during regrid.<br>

     This should be OK to call any time after define and before run.
  */
  void verbosity (int a_verbosity);

  ///
  /**
     Sets the regridding intervals.  This should be OK to call any time
     after define() and before run.
   */
  void regridIntervals(const Vector<int>& a_regridIntervals);

  ///
  /** Set maximum factor by which a timestep can grow. */
  void maxDtGrow(Real a_dtGrowFactor);

  ///
  /** Set amount by which two times may differ and still be considered
      equal (used to determine whether we've reached stop time yet)
      In practice, this is multiplied by the current base-level timestep,
      so we consider ourselves finished if (stopTime - time) < timeEps*dt0
  */
  void timeEps(Real a_timeEps);

  ///
  /** Set a fixed timestep. This must be called before calling
      setupForFixedHierarchyRun(), setupForNewAMRRun(), or setupForRestart().
      Note also that fixedDt is not currently saved in checkpoint files
  */
  void fixedDt(Real a_dt);

  ///
  /**
      Sets the initial time before starting compututation; this must be
      called before calling setupForFixedHierarchyRun() or setupForNewAMRRun().
      Note that this will be over-ridden by m_checkpointfile if
      restarting using setupForRestart().
  */
  void initialTime(Real a_initialTime);

  ///
  /**
     Turn subcycling in time off or on.   Default is true (use cycling).
  */
  void useSubcyclingInTime(bool a_useSubcycling);

  ///
  /**
     Sets the plot file prefix.

     Should be called after define()
     and before setup.
  */
  void plotPrefix(const std::string& a_plotfile_prefix);

  /**@}*/

  /**
     \name Access functions
  */

  /**@{*/

  ///
  /**
     Has a define() function been called?  Lots of these
     functions will assert fail if not.
   */
  bool isDefined() const;

  ///
  /**
     Has a setup function been called?  If not,
     you can't call run().
   */
  bool isSetUp() const;

  ///
  /**
     Returns the maximum grid size.
   */
  int maxGridSize() const;

  ///
  /**
     Returns the maximum grid size for level 0 grids.
   */
  int maxBaseGridSize() const;

  ///
  /**
     Returns current verbosity level.

     -inf to 0: print nothing.<br>
     1: prints a message every coarse time step.<br>
     2: also prints a message when regridding and writing files.<br>
     3: also prints function trace, and some more inf during file I/O.<br>
     4 to inf: also prints list of boxes and processor maps during regrid.<br>
  */
  int verbosity () const;

  ///
  /** Returns current maximum factor by which a timestep can grow. */
  Real maxDtGrow() const;

  ///
  /** Returns amount by which two times may differ and still be considered
      equal (used to determine whether we've reached stop time yet)
  */
  Real timeEps() const;

  ///
  /** Returns the fixed timestep. */
  Real fixedDt() const;

  ///
  /**
     Returns a vector of all the AMR levels
   */
  Vector<AMRLevel*> getAMRLevels();

  ///
  /**
     Returns the current time
   */
  Real getCurrentTime() const;

#ifdef CH_USE_TIMER
  ///
  /**
     Returns the Chombo timer, and resets it if the arg is not NULL.
  */
  Chombo::Timer * timer(Chombo::Timer *a_timer = NULL );
#endif

  /**@}*/

protected:

  bool m_useSubcycling;
  // advance by dt on this level and all finer levels and return the
  // number of steps left - given the number of steps left on entry.
  int timeStep(int a_level, int a_stepsLeft, bool a_coarseTimeBoundary);

  // internal use only
  void clearMemory();

  // make new grids.
  void regrid(int a_base_level);

  // make initial grids.
  // not an accurate name, since there are no grids to re.
  void initialGrid();

  void writePlotFile() const;

  void writeCheckpointFile() const;

  // computes maximum stable time step given the maximum stable time
  // step on the individual levels.
  void assignDt();

  // whether regridding should be done now.
  bool needToRegrid(int a_level, int a_numStepsLeft) const;

  void makeBaseLevelMesh (Vector<Box>& a_grids) const;

  void setDefaultValues();

  int  m_blockFactor;
  Real m_fillRatio;
  int  m_max_level;
  int  m_finest_level_old;
  int  m_finest_level;
  int  m_checkpoint_interval;
  int  m_plot_interval;
  Real m_plot_period;
  Real m_next_plot_time;
  int  m_max_grid_size;
  int  m_max_base_grid_size;
  Real m_dt_tolerance_factor;
  Real m_fixedDt;
  bool              m_checkForSteadyState;
  bool              m_isDefined;
  bool              m_isSetUp;
  Vector<AMRLevel*> m_amrlevels;
  Vector<int>       m_ref_ratios;
  Vector<int>       m_reduction_factor;
  Vector<int>       m_regrid_intervals;

  Real         m_dt_base;
  // New (maximum) dt
  Vector<Real> m_dt_new;
  // Current dt
  Vector<Real> m_dt_cur;
  Real         m_maxDtGrow;
  Real         m_time_eps;

  RefCountedPtr<MeshRefine> m_mesh_refine_ptr;
  bool         m_use_meshrefine;

  Vector<Vector<Box> > m_amr_grids;

  int m_cur_step;
  int m_restart_step;
  int m_lastcheck_step;

  Real         m_cur_time;
  Vector<int>  m_steps_since_regrid;
  Vector<long long> m_cell_updates;

  std::string m_plotfile_prefix;
  std::string m_checkpointfile_prefix;

  int m_verbosity;

  RefCountedPtr<Scheduler> m_scheduler;
#ifdef CH_USE_TIMER
  Chombo::Timer *m_timer;  //assumes the application manages the memory
#endif
};

#include "NamespaceFooter.H"
#endif
