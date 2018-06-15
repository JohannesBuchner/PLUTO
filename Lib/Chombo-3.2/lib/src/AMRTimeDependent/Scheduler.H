#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SCHEDULER_H_
#define _SCHEDULER_H_

#include <map>
#include "CH_assert.H"
#include "REAL.H"
#include "RefCountedPtr.H"
#include "NamespaceHeader.H"

// Forward declaration of AMR class.
class AMR;

//! \class Scheduler
//! This class executes functions at periodic intervals during an
//! AMR simulation.
class Scheduler
{
  public:

  //! \class PeriodicFunction
  //! This function is called periodically by the Scheduler. Its
  //! call operator must be overridden to define its behavior.
  class PeriodicFunction
  {
    public:

    //! Default construction.
    PeriodicFunction();

    //! Destructor.
    virtual ~PeriodicFunction();

    //! Override this method to prepare the periodic function to interact
    //! with the given AMR object when set to an interval in time steps.
    //! By default this does nothing.
    //! \param a_AMR The AMR object with which this function will interact
    //!              during periodic function calls.
    //! \param a_interval The interval (in steps) at which the periodic function
    //!                   is to be called.
    virtual void setUp(AMR& a_AMR, int a_interval);

    //! Override this method to prepare the periodic function to interact
    //! with the given AMR object when set to an interval in time units.
    //! By default this does nothing.
    //! \param a_AMR The AMR object with which this function will interact
    //!              during periodic function calls.
    //! \param a_interval The interval (in simulation time units) at which
    //! the periodic function is to be called.
    virtual void setUp(AMR& a_AMR, Real a_interval);

    //! Override this operator to define the behavior of the periodic function.
    //! \param a_step The step at which the function is called.
    //! \param a_time The simulation time at which the function is called.
    virtual void operator()(int a_step, Real a_time) = 0;

    //! Override this operator to perform a task at the conclusion of a
    //! simulation. This is called by the conclude() method of AMR and
    //! is fed the final step and simulation time. By default, this does
    //! nothing.
    //! \param a_step The step at which the function is called.
    //! \param a_time The simulation time at which the function is called.
    virtual void conclude(int a_step, Real a_time);

    private:

    PeriodicFunction(const PeriodicFunction&);
    PeriodicFunction& operator=(const PeriodicFunction&);
  };

  //! Ordering operator for pointers to periodic functions.
  struct PeriodicFunctionLessThan
  {
    bool operator()(const RefCountedPtr<PeriodicFunction>& a_lhs,
                    const RefCountedPtr<PeriodicFunction>& a_rhs) const
    {
      CH_assert(!a_lhs.isNull() && !a_rhs.isNull());
      // Just use the pointer's address to order the two.
      return (&*a_lhs < &*a_rhs);
    }
  };

  //! Default constructor. Makes an empty schedule.
  Scheduler();

  //! Destructor.
  virtual ~Scheduler();

  //! Add a periodic function that is called every \a a_interval steps.
  //! \param a_function The function to be called periodically.
  //! \param a_interval The number of steps that elapse between calls to \a a_function.
  void schedule(RefCountedPtr<PeriodicFunction> a_function,
                int a_interval);

  //! Add a periodic function that is called every \a a_interval time units.
  //! \param a_function The function to be called periodically.
  //! \param a_interval The number of simulation time units that elapse between calls to \a a_function.
  void schedule(RefCountedPtr<PeriodicFunction> a_function,
                Real a_interval);

  //! Execute the schedule on the given step and at the given time.
  //! \param a_step The step in the simulation at which the schedule is executed.
  //! \param a_time The simulation time at which the schedule is executed.
  void execute(int a_step, Real a_time) const;

  //! This function is called by the associated AMR object to set up
  //! interactions between itself and the periodic functions within the scheduler.
  void setUp(AMR& a_AMR);

  //! This function is called by the associated AMR object upon the conclusion
  //! of a simulation.
  //! \param a_step The step in the simulation at which the schedule is executed.
  //! \param a_time The simulation time at which the schedule is executed.
  void conclude(int a_step, Real a_time) const;

  private:

  // Step-triggered functions.
  mutable std::map<RefCountedPtr<PeriodicFunction>, int, PeriodicFunctionLessThan> m_stepTriggeredFunctions;

  // Time-triggered functions.
  mutable std::map<RefCountedPtr<PeriodicFunction>, std::pair<Real, Real>, PeriodicFunctionLessThan> m_timeTriggeredFunctions;

  // Forbidden.
  Scheduler(const Scheduler&);
  Scheduler& operator=(const Scheduler&);
};

//! \class PlotterPeriodicFunction
//! This placebo allows one to enable periodic plots using the Scheduler
//! mechanism.
class PlotterPeriodicFunction: public Scheduler::PeriodicFunction
{
  public:

  explicit PlotterPeriodicFunction(const std::string& a_prefix);
  void operator()(int a_step, Real a_time);

  void setUp(AMR& a_AMR, int a_interval);
  void setUp(AMR& a_AMR, Real a_interval);
  void conclude(int a_step, Real a_time);

  private:

  std::string m_prefix;
};

//! \class CheckpointPeriodicFunction
//! This placebo allows one to enable periodic checkpoints using the Scheduler
//! mechanism.
class CheckpointPeriodicFunction: public Scheduler::PeriodicFunction
{
  public:

  explicit CheckpointPeriodicFunction(const std::string& a_prefix);

  void setUp(AMR& a_AMR, int a_interval);
  void operator()(int a_step, Real a_time);

  private:

  std::string m_prefix;
};

#include "NamespaceFooter.H"
#endif
