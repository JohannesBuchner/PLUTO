#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _OLDPATCHGODUNOV_H_
#define _OLDPATCHGODUNOV_H_

#include "OldPhysIBC.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"

#include <string>

#include "UsingNamespace.H"

using std::string;

///
/**
   The base class OldPatchGodunov provides an implementation of a second-order,
   unsplit Godunov method acting on a single grid/patch.  OldPatchGodunov
   provides an interface to the level integrator, OldLevelGodunov, which manages
   the entire level and flux corrections (via flux registers) to the coarser
   and finer levels.  In addition, the physics dependent code is not provided
   in OldPatchGodunov but is supplied by the user by subclassing OldPatchGodunov and
   implementing the pure virtual functions.  Some parameters can also be
   adjusted to modify the algorithm.  All functions are virtual so any of them
   can be reimplemented by the user.

   There are three types of grid variables that appear in the unsplit Godunov
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
class OldPatchGodunov
{
public:
  /// Constructor
  /**
   */
  OldPatchGodunov();

  /// Destructor
  /**
   */
  virtual ~OldPatchGodunov();

  /// Define the object
  /**
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&    a_dx);

  /// Set the initial and boundary condition object
  /**
   */
  virtual void setPhysIBC(OldPhysIBC* a_bc);

  /// Get the initial and boundary condition object
  /**
   */
  virtual OldPhysIBC* getPhysIBC() const;

  /// Set parameters for slope computations
  /**
   */
  virtual void setSlopeParameters(bool a_fourthOrderSlopes,
                                  bool a_flattening,
                                  bool a_limitSlopes = true);

  /// Fourth-order slope query
  /**
     Return true if you are using fourth-order slopes and return false if
     you are using second-order slopes.
   */
  virtual bool useFourthOrderSlopes();

  /// Slope flattening query
  /**
     Return true if the application is using slope flattening.
   */
  virtual bool useFlattening();

  /// Slope limiter query
  /**
     Return true if the application is using slope limiting
   */
  virtual bool limitSlopes();

  /// Set parameters for artificial viscosity
  /**
   */
  virtual void setArtificialViscosity(bool a_useArtificialViscosity,
                                      Real a_artificialViscosity);

  /// Artificial viscosity query
  /**
     Return true if the application is using artificial viscosity.
   */
  virtual bool useArtificialViscosity();

  /// Artificial viscosity coefficient.
  /**
   */
  virtual Real artificialViscosityCoefficient();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to new OldPatchGodunov object with its initial and boundary
     condtions, slope parameters, and artificial viscosity information defined.
   */
  virtual OldPatchGodunov* new_patchGodunov() const = 0;

  /// Set the current time before calling updateState()
  /**
   */
  virtual void setCurrentTime(const Real& a_currentTime);

  /// Set the current box before calling updateState()
  /**
   */
  virtual void setCurrentBox(const Box& a_currentBox);

  /// Update the conserved variables and return the final fluxes used for this
  /**
     Update the conserved variables and return the final fluxes that were used
     for this.  Defines the container for the fluxes (a_F) and calls
     computeFluxes to compute fluxes for the update. If there are no source
     terms then a_S should be null constructed.  Also return the maximum
     wave speed.
   */
  virtual void updateState(FArrayBox&       a_U,
                           FArrayBox        a_F[CH_SPACEDIM],
                           Real&            a_maxWaveSpeed,
                           const FArrayBox& a_S,
                           const Real&      a_dt,
                           const Box&       a_box);

  /// Compute the fluxes which will be used for the update
  /**
     Compute the fluxes using a second-order, unsplit Godunov method
     based on the input conserved variables, a_U, and source terms, a_S.  If
     there are no source terms then a_S should be null constructed.
     Assumes a_F has been defined already.
   */
  virtual void computeFluxes(FArrayBox&       a_U,
                             FArrayBox        a_F[CH_SPACEDIM],
                             const FArrayBox& a_S,
                             const Real&      a_dt,
                             const Box&       a_box);

  /// Compute the maximum wave speed
  /**
   */
  virtual Real getMaxWaveSpeed(const FArrayBox& a_U,
                               const Box&       a_box) = 0;

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
   */
  virtual int numConserved() = 0;

  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numConserved()-1.
   */
  virtual Vector<string> stateNames();

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
   */
  virtual int numFluxes() = 0;


protected:
  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
   */
  virtual bool isDefined() const;

  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  virtual int numPrimitives() = 0;

  /// Number of primitive variables for which slopes are computed
  /**
     Return the number of primitive variables for which slopes are computed.
     Only slopes corresponding to primitive variables in the interval 0 to
     numSlopes() - 1 are computed and only primitive variables in that interval
     are updated using the slopes.
   */
  virtual int numSlopes() = 0;

  /// Compute the primitive variables from the conserved variables within a_box
  /**
   */
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box) = 0;

  /// Compute the slope flattening coefficients
  /**
     Compute the slope flattening coefficients, a_flattening, using the
     primitive variables, a_W, within a_box.
   */
  virtual void computeFlattening(FArrayBox&       a_flattening,
                                 const FArrayBox& a_W,
                                 const Box&       a_box);

  /// Compute the 2nd or 4th order slopes of the primitive variables
  /**
     Compute the limited slope, a_dW, of the primitive variables, a_W,
     over the range of indices, 0 to numSlopes()-1.  This only used
     a_flattening if 4th order slopes are being computed and slope
     flattening is turned on.  This also calls applyLimiter() to do slope
     limiting.
   */
  virtual void slope(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const FArrayBox& a_flattening,
                     const int&       a_dir,
                     const Box&       a_box);

  /// Extrapolate the primitive variables to the cell faces
  /**
     Extrapolate the primitive variables, a_W, in the minus and plus direction,
     a_dir, using the slopes a_dW within a_box.  See document for
     details.
   */
  virtual void normalPred(FArrayBox&       a_WMinus,
                          FArrayBox&       a_WPlus,
                          const FArrayBox& a_W,
                          const FArrayBox& a_dW,
                          const Real&      a_scale,
                          const int&       a_dir,
                          const Box&       a_box) = 0;

  /// Increment the primitive variables by a source term
  /**
     The default implementation does nothing to the primitive variables.
   */
  virtual void incrementWithSource(FArrayBox&       a_W,
                                   const FArrayBox& a_S,
                                   const Real&      a_scale,
                                   const Box&       a_box);

  /// Compute a Riemann problem and generate fluxes at the faces
  /**
     Given input left and right states in a direction, a_dir, compute a
     Riemann problem and generate fluxes at the faces within a_box.
   */
  virtual void riemann(FArrayBox&       a_F,
                       const FArrayBox& a_WLeft,
                       const FArrayBox& a_WRight,
                       const int&       a_dir,
                       const Box&       a_box) = 0;

  /// Update the primitive variables using fluxes (of the conserved variables)
  /**
     Given the fluxes, a_F, in a direction, a_dir, and a scaling factor,
     a_scale, update the primitive variables, a_WMinus and a_WPlus, within
     a_box.  Note:  The fluxes are fluxes of conserved variables.
   */
  virtual void updatePrim(FArrayBox&       a_WMinus,
                          FArrayBox&       a_WPlus,
                          const FArrayBox& a_F,
                          const Real&      a_scale,
                          const int&       a_dir,
                          const Box&       a_box) = 0;

  /// Apply artificial viscosity to the fluxes
  /**
     Return fluxes, a_F, updated using artificial viscosity based on the
     conserved variables, a_U, and the face centered divergence of the
     velocity, a_divVel, in a given direction, a_dir, within a_box.
   */
  virtual void artificialViscosity(FArrayBox&       a_F,
                                   const FArrayBox& a_U,
                                   const FArrayBox& a_divVel,
                                   const int&       a_dir,
                                   const Box&       a_box);

  /// Update the conserved variable using fluxes and a scaling factor
  /**
     Given the fluxes, a_F, in a direction, a_dir, and a scaling factor,
     a_scale, update the conserved variables, a_U, within a_box:

     a_U = a_U - a_scale * (a_F(a_dir,HiSide) - a_F(a_dir,LoSide))
   */
  virtual void updateCons(FArrayBox&       a_U,
                          const FArrayBox& a_F,
                          const Real&      a_scale,
                          const int&       a_dir,
                          const Box&       a_box) = 0;

  /// Perform final update of conserved variable using fluxes
  /**
     Given the fluxes, a_F, in a direction, a_dir, and a scaling factor,
     a_scale, update the conserved variables, a_U, within a_box:

     a_U = a_U - a_scale * (a_F(a_dir,HiSide) - a_F(a_dir,LoSide))
   */
  virtual void finalUpdate(FArrayBox&       a_U,
                           const FArrayBox& a_F,
                           const Real&      a_scale,
                           const int&       a_dir,
                           const Box&       a_box) = 0;

  virtual void postUpdateCons(FArrayBox&       a_U,
                          const FArrayBox& a_Uold,
                          const Real&      a_dt,
                          const Real&      a_dx,
                          const Box&       a_box);

  /// Interval within the primitive variables corresponding to the velocities
  /**
     Return the interval of component indices within the primitive variable
     of the velocities.  Used for slope flattening (slope computation) and
     computing the divergence of the velocity (artificial viscosity).
   */
  virtual Interval velocityInterval() = 0;

  /// Component index within the primitive variables of the pressure
  /**
     Return the component index withn the primitive variables for the
     pressure.  Used for slope flattening (slope computation).
   */
  virtual int pressureIndex() = 0;

  /// Component index within the primitive variables of the bulk modulus
  /**
     Return the component index withn the primitive variables for the
     bulk modulus.  Used for slope flattening (slope computation) used
     as a normalization to measure shock strength.
   */
  virtual int bulkModulusIndex() = 0;

  /// Apply a slope limiter to computed slopes
  /**
     Given the center difference, a_dW, and the left and right differences,
     a_dWLeft and a_dWRight, in direction, a_dir, apply a slope limiter
     within a_box to generate final slopes, a_dW.  This is called by slope().

     A default implementation is provided which implements a van Leer limiter,
     see the documentation for details.
   */
  virtual void applyLimiter(FArrayBox&       a_dW,
                            const FArrayBox& a_dWLeft,
                            const FArrayBox& a_dWRight,
                            const int&       a_dir,
                            const Box&       a_box);

  /// Compute the face centered divergence of the velocity
  /**
     Compute the face centered divergence of the velocity, a_divVel, from
     the primitive variables, a_W, in a direction, a_dir, within a_box.
   */
  virtual void divVel(FArrayBox&       a_divVel,
                      const FArrayBox& a_W,
                      const int        a_dir,
                      const Box&       a_box);

  // Has define() been called
  bool m_isDefined;

  // Problem domain and grid spacing
  ProblemDomain m_domain;
  Real m_dx;

  // Slope computation flags and have they been set
  bool m_useFourthOrderSlopes;
  bool m_useFlattening;
  bool m_limitSlopes;

  // Have they been set
  bool m_isSlopeSet;

  // Artificial viscosity flag, coefficient, and have they been set
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;
  bool m_isArtViscSet;

  // Initial and boundary condition object and has it been set
  OldPhysIBC* m_bc;
  bool m_isBCSet;

  // Current time and has it been set
  Real m_currentTime;
  bool m_isCurrentTimeSet;

  // Current box and has it been set
  Box m_currentBox;
  bool m_isCurrentBoxSet;


private:
  // Disallowed for all the usual reasons
  void operator=(const OldPatchGodunov& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  OldPatchGodunov(const OldPatchGodunov& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
