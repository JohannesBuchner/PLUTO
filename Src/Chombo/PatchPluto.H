#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCH_GODUNOV_H_
#define _PATCH_GODUNOV_H_

#include "FArrayBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "RealVect.H"
#include "AMRIO.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "FluxBox.H"

#include "PatchGrid.H"

#if (COOLING != NO)
 #define SKIP_SPLIT_CELLS
#endif
#define SKIP_SPLIT_CELLS

#include <string>
using std::string;

#include "NamespaceHeader.H"

///
/**
   The base class PatchPluto provides an implementation of a second-order,
   unsplit Pluto method acting on a single grid/patch.  PatchPluto
   provides an interface to the level integrator, LevelPluto, which manages
   the entire level and flux corrections (via flux registers) to the coarser
   and finer levels.  In addition, the physics dependent code is not provided
   in PatchPluto but is supplied by the user by subclassing PatchPluto and
   implementing the pure virtual functions.  Some parameters can also be
   adjusted to modify the algorithm.  All functions are virtual so any of them
   can be reimplemented by the user.

   There are three types of grid variables that appear in the unsplit Pluto
   method: conserved variables, primitive variables, and fluxes, denoted below
   by U, W, F, respectively.  It is often convenient to have the number of
   primitive variables and fluxes exceed the number of conserved variables.
   In the case of primitive variables, redundant quantities are carried that
   parameterize the equation of state in order to avoid multiple calls to that
   the equation of state function.  In the case of fluxes, it is often
   convenient to split the flux for some variables into multiple components,
   e.g., dividing the momentum flux into advective and pressure terms.  The
   API given here provides the flexibility to support these various
   possibilities.
 */
class PatchPluto
{
public:
  /// Constructor
  /**
   */
  PatchPluto();

  /// Destructor
  /**
   */
  virtual ~PatchPluto();

  /// Define the object
  /**
   */
  virtual void define(ProblemDomain& a_domain,
                      const Real&    a_dx,
                      const int&     a_level,
                      const int&     a_numGhost);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to new PatchPluto object with its boundary
     conditions and Riemann solver information defined.
   */
  virtual PatchPluto* new_patchPluto() const;

  /// Set the current time before calling updateState()
  /**
   */
  virtual void setCurrentTime(const Real& a_currentTime);

  /// Set the current box before calling updateState()
  /**
   */
  virtual void setCurrentBox(const Box& a_currentBox);

  /// Set the value of conserved variables on ghost cells out of
  /// the computational domain.

  virtual void setBoundary(const int *leftbound, 
                           const int *rightbound);

  virtual void setRiemann(Riemann_Solver riem);

  virtual void setGrid(const Box&    a_box,
                       struct GRID  *grid,
                       FArrayBox& a_dV);

  virtual void initialize(LevelData<FArrayBox>& a_U);

  virtual void convertFArrayBox(FArrayBox& U);

  virtual void saveFluxes (const State_1D *, int, int, Grid *);
  virtual void showPatch (Grid *);
  virtual void startup(FArrayBox& a_U);

  /// Update the conserved variables and return the final fluxes used for this
  /**
     Update the conserved variables and return the final fluxes that were used
     for this.  Compute the fluxes using a second-order, unsplit Pluto method
     based on the input conserved variables, a_U. Also return the maximum wave 
     speed. 
   */

  virtual void advanceStep(FArrayBox&       a_U,
                           FArrayBox&       a_Utmp,
                           const FArrayBox& a_dV,
                           FArrayBox&       a_split_tags,
                           BaseFab<unsigned char>&   a_Flags,
                           FluxBox&         a_F,
                           Time_Step        *Dts,
                           const Box&       a_box,
                           Grid *grid);

  virtual void getPrimitiveVars (Data_Arr U, Data *d, Grid *grid);

  virtual void computeRefGradient(FArrayBox&, FArrayBox&, const FArrayBox &, const Box&);

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
   */
  virtual int numConserved();

  /// Names of the conserved variables
  /**
     Return the names of the conserved and primitive variables.  
     A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numConserved()-1.
   */
  virtual Vector<string> ConsStateNames();
  virtual Vector<string> PrimStateNames();

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
   */
  virtual int numFluxes();

  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
   */
  virtual bool isDefined() const;

protected:

  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  virtual int numPrimitives();

  /// Component index within the primitive variables of the pressure
  /**
     Return the component index withn the primitive variables for the
     pressure.  Used for slope flattening (slope computation).
   */
  virtual int pressureIndex();

  // Has define() been called
  bool m_isDefined;

  // Problem domain, grid spacing, level number and number of ghost cells
  ProblemDomain m_domain;
  Real m_dx;
  int  m_level;
  int  m_numGhost;

  // Current time and has it been set
  Real m_currentTime;
  bool m_isCurrentTimeSet;

  // Current box and has it been set
  Box m_currentBox;
  bool m_isCurrentBoxSet;

  // Boundary type
  int left_bound_side[3];
  int right_bound_side[3];
  bool m_isBoundarySet;

  // Riemann Solver
  Riemann_Solver *rsolver;
  bool m_isRiemannSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const PatchPluto& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PatchPluto(const PatchPluto& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
