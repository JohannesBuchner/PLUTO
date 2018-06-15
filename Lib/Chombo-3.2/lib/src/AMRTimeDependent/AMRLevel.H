#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVEL_H_
#define _AMRLEVEL_H_

#include "Vector.H"
#include "REAL.H"
#include "Box.H"
#include "ProblemDomain.H"
#include "DataIterator.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "CH_HDF5.H"
#include "NamespaceHeader.H"

//class HDF5Handle;
//class IntVectSet;
//template <class TTYPE> class Vector;

/// Abstract base class for time-dependent data at a level of refinement.
/**
   This is an abstract base class for data at the same level of
   refinement within a hierarchy of levels.  The concrete class
   derived from AMRLevel is called a  physics class.  The domain
   of a level is a disjoint union of rectangles in a logically
   rectangular index space.  Data is defined within this domain.
   There is also a problem domain, which may be larger, within which
   data can, in theory, be interpolated from some coarser level.

   AMRLevel is intimately entangled with the AMR class.
   The AMR contains a collection of AMRLevel objects which
   represent a hierarchy of levels of refinement.  The class  AMR
   is responsible for calling the correct member functions of
   AMRLevel.  The user is responsible for implementing the physics
   class, and might never call any AMRLevel functions described
   here.
*/

class AMRLevel
{
public:
  ///
  /**
     Default constructor.
  */
  AMRLevel();

  ///
  /**
     Destructor.
  */
  virtual
    ~AMRLevel();

  ///
  /**
     Defines this AMRLevel.

     <ul>
       <li>a_coarser_level_ptr (not modified): pointer to next coarser level
           object.
       <li>a_problem_domain (not modified): problem domain of this level.
       <li>a_level (not modified): index of this level.  The base level is
           zero.
       <li>a_ref_ratio (not modified): the refinement ratio between this level
           and the next finer level.
     </ul>
  */
  virtual
    void define(AMRLevel*  a_coarser_level_ptr,
                const Box& a_problem_domain,
                int        a_level,
                int        a_ref_ratio);

  ///
  /**
     Defines this AMRLevel.

     <ul>
       <li>a_coarser_level_ptr (not modified): pointer to next coarser level
           object.
       <li>a_problem_domain (not modified): problem domain of this level.
       <li>a_level (not modified): index of this level.  The base level is
           zero.
       <li>a_ref_ratio (not modified): the refinement ratio between this level
           and the next finer level.
     </ul>
  */
  virtual
    void define(AMRLevel*            a_coarser_level_ptr,
                const ProblemDomain& a_problem_domain,
                int                  a_level,
                int                  a_ref_ratio);

  ///
  /**
     Advances this level by one time step.  Returns an estimate of the
     new time step.

     This is a pure virtual function and MUST be defined in the derived
     class.

  */
  virtual
    Real advance() = 0;

  ///
  /**
     Return true if the solution has converged to  a steady state for this level.
     Default version always returns false.
   */
  virtual
  bool convergedToSteadyState()
  {
    return false;
  }
    
  ///
  /**
     Things to do after advancing this level by one time step.

     This is a pure virtual function and  MUST be defined in the derived
     class.

  */
  virtual
    void postTimeStep() = 0;

  ///
  /**
     Creates tagged cells for dynamic mesh refinement.

     This is a pure virtual function and  MUST be defined in the derived
     class.

  */
  virtual
    void tagCells(IntVectSet& a_tags) = 0;

  ///
  /**
     Creates tagged cells for mesh refinement at initialization.

     This is a pure virtual function and MUST be defined in the derived
     class.

  */
  virtual
    void tagCellsInit(IntVectSet& a_tags) = 0;

  ///
  /**
      Performs any pre-regridding operations which are necessary.

      This is not a pure virtual function to preserve compatibility
      with earlier versions of AMRLevel.  The AMRLevel::preRegrid()
      instantiation is a no-op.
  */
  virtual
    void preRegrid(int a_base_level, const Vector<Vector<Box> >& a_new_grids);

  ///
  /**
     Redefines this level to have the specified domain a_new_grids.

     This is a pure virtual function and MUST be defined in the derived
     class.

  */
  virtual
    void regrid(const Vector<Box>& a_new_grids) = 0;

  ///
  /**
      Performs any post-regridding operations which are necessary.

      This is not a pure virtual function to preserve compatibility
      with earlier versions of AMRLevel.  The AMRLevel::postRegrid()
      instantiation is a no-op.
  */
  virtual
    void postRegrid(int a_base_level);

  ///
  /**
     Initializes this level to have the specified domain a_new_grids.

     This is a pure virtual function and MUST be defined in the derived
     class.

  */
  virtual
    void initialGrid(const Vector<Box>& a_new_grids) = 0;

  ///
  /**
     Performs operations required after the grid has been defined but
     before data initialization.  This will also be called after
     readCheckpointLevel during a restart procedure with argument
     a_restart set to true.

     Levels are accessed from finest to coarsest.  The
     AMRLevel::postInitialGrid() instantiation is a no-op.

  */
  virtual
    void postInitialGrid(const bool a_restart);

  ///
  /**
     Initializes data.

     This is a pure virtual function and  MUST be defined in the derived
     class.

  */
  virtual
    void initialData() = 0;

  ///
  /**
     Things to do after initialization.

     This is a pure virtual function and  MUST be defined in the derived
     class.

  */
  virtual
    void postInitialize() = 0;

  //! Override this method to have an AMRLevel subclass perform some
  //! operation upon the conclusion of a simulation. This is called
  //! when AMR::conclude() is called. The final step is passed to the
  //! method.
  //! \param a_step The last step in the simulation.
  virtual void conclude(int a_step) const;

  /**
     \name I/O functions
  */
  /**@{*/

#ifdef CH_USE_HDF5
  ///
  /**
     Writes checkpoint header.

     This is a pure virtual function and  MUST be defined in the derived
     class.

  */
  virtual
    void writeCheckpointHeader (HDF5Handle& a_handle) const = 0;

  ///
  /**
     Write checkpoint data for this level.

     This is a pure virtual function and  MUST be defined in the derived
     class.

  */
  virtual
    void writeCheckpointLevel (HDF5Handle& a_handle) const = 0;

  ///
  /**
     Reads checkpoint header.

     This is a pure virtual function and  MUST be defined in the derived
     class.

  */
  virtual
    void readCheckpointHeader (HDF5Handle& a_handle) = 0;

  ///
  /**
     Reads checkpoint data for this level.

     This is a pure virtual function and MUST be defined in the derived
     class.

  */
  virtual
    void readCheckpointLevel (HDF5Handle& a_handle) = 0;

  ///
  /**
     Writes plot header.

     This is a pure virtual function and MUST be defined in the derived
     class.

  */
  virtual
    void writePlotHeader (HDF5Handle& a_handle) const = 0;

  ///
  /**
     Write plot file for this level.

     This is a pure virtual function and MUST be defined in the derived
     class.

  */
  virtual
    void writePlotLevel (HDF5Handle& a_handle) const = 0;
#endif

  //! This allows one to write a plot file in a non-HDF5 format. It is called only at
  //! refinement level 0, so AMR data will have to be handled by the implementer.
  //! \param a_prefix A prefix for the custom plot file name.
  //! \param a_step The current time step.
  virtual void writeCustomPlotFile(const std::string& a_prefix,
                                   int a_step) const;

  /**@}*/

  /**
     \name Parameter-setting functions
  */
  /**@{*/

  ///
  /**
     Sets the pointer-to-finer-level member to a_finer_level_ptr.
  */
  virtual
    void finerLevelPtr(AMRLevel* a_finer_level_ptr);

  ///
  /**
     Sets the time step to a_dt.
  */
  virtual
    void dt(Real a_dt);

  ///
  /**
     Sets the time to a_time.

  */
  virtual
    void time(Real a_time);

  ///
  /**
     Sets the initial dt multiplier to a_initial_dt_multiplier.
  */
  virtual
    void initialDtMultiplier(Real a_initial_dt_multiplier);

  /**@}*/

  /**
     \name Access functions
  */
  /**@{*/

  ///
  /**
     Returns the current value of the time step.

  */
  virtual
    Real dt() const;

  ///
  /**
     Returns the current value of the time on this level.
  */
  virtual
    Real time() const;

  ///
  /**
     Returns the initial dt multiplier.
  */
  virtual
    Real initialDtMultiplier() const;

  ///
  /**
     Returns the problem domain of this level.

  */
  virtual
    const ProblemDomain& problemDomain() const;

  ///
  /**
     Returns the domain of this level.

  */
  virtual
    Vector<Box> boxes() const;

  ///
  /**
     Returns true if any AMRLevel::define function has been called,
     false otherwise.
  */
  bool isDefined() const;

  ///
  /**
     Returns true if a coarser level exists, is defined, and has a grid.
  */
  bool hasCoarserLevel() const;

  ///
  /**
     Returns true if a finer level exists, is defined, and has a grid.
  */
  bool hasFinerLevel() const;

  ///
  /**
     Returns the index of this level
  */
  virtual
    int level() const;

  ///
  /**
     Returns the refinement ratio between this level and the next finer level.

  */
  virtual
    int refRatio() const;

  ///
  /**
     Returns maximum stable time step for this level.

     This is a pure virtual function and  MUST be defined in the derived
     class.

  */
  virtual
    Real computeDt() = 0;

  ///
  /**
      Returns maximum stable time step for this level with initial data.

      This is a pure virtual function and  MUST be defined in the derived
      class.

  */
  virtual
    Real computeInitialDt() = 0;

  //! Retrieve an array of all of the AMRLevel objects in the entire hierarchy.
  Vector<AMRLevel*> getAMRLevelHierarchy();

  /**@}*/

  ///
  /**
     Returns current verbosity level.  Minimum verbosity is 0, for
     which nothing is printed.

  */
  static
    int verbosity();

  ///
  /**
     Sets verbosity level to a_verbosity.  Minimum verbosity is 0, for
     which nothing is printed.
  */
  static
    void verbosity(int a_verbosity);

protected:

  // verbosity level
  static int s_verbosity;

  // the problem domain
  ProblemDomain m_problem_domain;

  //
  Vector<Box> m_level_grids;

  // the level
  int m_level;

  // refinement ratio between this level and the next finer
  int m_ref_ratio;

  // initial time step multipier
  Real m_initial_dt_multiplier;

  // time step
  Real m_dt;

  // time
  Real m_time;

  // pointer to next coarser level
  AMRLevel* m_coarser_level_ptr;

  // pointer to next finer level
  AMRLevel* m_finer_level_ptr;

  bool m_isDefined;
};

#include "NamespaceFooter.H"
#endif
