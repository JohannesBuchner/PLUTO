#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GODUNOVPHYSICS_H_
#define _GODUNOVPHYSICS_H_

#include <string>
using std::string;

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "PhysIBC.H"
#include "FluxBox.H"

#include "GodunovUtilities.H"
#include "NamespaceHeader.H"

class HDF5HeaderData;

///
/**
   The base class GodunovPhysics provides the physics-dependent components
   for a higher-order Godunov method for a single patch: characteristic
   analysis, Riemann solver, quasilinear update, conservative update,
   and transformations between conserved, primitive, and flux variables.
   This class is essentially pure virtual; i.e., all of its member functions
   are virtual; and the ones that have default implementations are ones
   that are optionally defined; i.e., the default definition is to send
   an error message. Physics-dependent versions of this class that are
   required in real applications are derived from this class by inheritance.
*/
class GodunovPhysics
{
public:
  /// Constructor
  /**
   */
  GodunovPhysics();

  /// Get the initial and boundary condition object
  /**
   */
  PhysIBC* getPhysIBC() const;

  /// Set the initial and boundary condition object
  /**
   */
  void setPhysIBC(PhysIBC* a_bc);

  /// Destructor
  /**
   */
  virtual ~GodunovPhysics();

  /// Define the object
  /**
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Set the current box (default implementation - do nothing)
  /**
   */
  virtual void setCurrentBox(const Box& a_currentBox);

  /// Compute the maximum wave speed.
  /**
   */
  virtual Real getMaxWaveSpeed(const FArrayBox& a_U,
                               const Box&       a_box) = 0;

  /// Compute the maximum wave speed
  /**
   */
  virtual void soundSpeed(FArrayBox& a_speed,
                          const FArrayBox& a_U,
                          const Box&       a_box);

  /// Object factory for this class
  /**
   */
  virtual GodunovPhysics* new_godunovPhysics() const = 0;

  /// Compute the increment in the conserved variables from face variables.
  /**
     Compute dU = dt*dUdt, the change in the conserved variables over
     the time step. The fluxes are returned are suitable for use in refluxing.
     This has a default implementation but can be redefined as needed.
   */
  virtual void computeUpdate(FArrayBox&       a_dU,
                             FluxBox&         a_F,
                             const FArrayBox& a_U,
                             const FluxBox&   a_WHalf,
                             const bool&      a_useArtificialViscosity,
                             const Real&      a_artificialViscosity,
                             const Real&      a_currentTime,
                             const Real&      a_dx,
                             const Real&      a_dt,
                             const Box&       a_box);

  /// Compute the fluxes from primitive variable values on a face
  /**
     This has a default implementation which throws an error.  The method is
     here so that the default implementation of "computeUpdate" can use it
     and the user can supply it.  It has an implementation so if the user
     redefines "computeUpdate" they aren't force to implement "getFlux" -
     which is only used by the default implementation of "computeUpdate".
   */
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_WHalf,
                       const int&       a_dir,
                       const Box&       a_box);

  /// Transform a_dW from primitive to characteristic variables
  /**
     On input, a_dW contains the increments of the primitive variables. On
     output, it contains the increments in the characteristic variables.

     IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
     smallest eigenvalue first, the largest eigenvalue last, and orders the
     characteristic variables accordingly.
  */
  virtual void charAnalysis(FArrayBox&       a_dW,
                            const FArrayBox& a_W,
                            const int&       a_dir,
                            const Box&       a_box) = 0;

  /// Transform a_dW from characteristic to primitive variables
  /**
     On input, a_dW contains the increments of the characteristic variables.
     On output, it contains the increments in the primitive variables.

     IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
     smallest eigenvalue first, the largest eigenvalue last, and orders the
     characteristic variables accordingly.
  */
  virtual void charSynthesis(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Box&       a_box) = 0;

  /// Compute the characteristic values (eigenvalues)
  /**
     Compute the characteristic values (eigenvalues).

     IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
     smallest eigenvalue first, the largest eigenvalue last, and orders the
     characteristic variables accordingly.
   */
  virtual void charValues(FArrayBox&       a_lambda,
                          const FArrayBox& a_W,
                          const int&       a_dir,
                          const Box&       a_box) = 0;

  /// Add to (increment) the source terms given the current state
  /**
     On input, a_S contains the current source terms.  On output, a_S has
     had any additional source terms (based on the current state, a_W)
     added to it.  This should all be done on the region defined by a_box.
  */
  virtual void incrementSource(FArrayBox&       a_S,
                               const FArrayBox& a_W,
                               const Box&       a_box) = 0;

  /// Compute the solution to the Riemann problem.
  /**
     Given input left and right states in a direction, a_dir, compute a
     Riemann problem and generate fluxes at the faces within a_box.
  */
  virtual void riemann(/// face-centered solution to Riemann problem
                       FArrayBox&       a_WStar,
                       /// left state, on cells to left of each face
                       const FArrayBox& a_WLeft,
                       /// right state, on cells to right of each face
                       const FArrayBox& a_WRight,
                       /// state on cells, used to set boundary conditions
                       const FArrayBox& a_W,
                       /// current time
                       const Real&      a_time,
                       /// direction of faces
                       const int&       a_dir,
                       /// face-centered box on which to set a_WStar
                       const Box&       a_box) = 0;

  /// Post-normal predictor calculation.
  /**
     Add increment to normal predictor, e.g. to account for source terms due to
     spatially-varying coefficients, to bound primitive variable ranges.
  */
  virtual void postNormalPred(FArrayBox&       a_dWMinus,
                              FArrayBox&       a_dWPlus,
                              const FArrayBox& a_W,
                              const Real&      a_dt,
                              const Real&      a_dx,
                              const int&       a_dir,
                              const Box&       a_box) = 0;

  /// Compute the quasilinear update A*dW/dx.
  /**
   */
  virtual void quasilinearUpdate(FArrayBox&       a_AdWdx,
                                 const FArrayBox& a_wHalf,
                                 const FArrayBox& a_W,
                                 const Real&      a_scale,
                                 const int&       a_dir,
                                 const Box&       a_box) = 0;

  /// Compute primitive variables from conserved variables.
  /**
   */
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box) = 0;

  /// Compute the artificial viscosity contribution to the flux
  /**
     Compute the artificial viscosity contribution to the flux.  This has
     a default implementation but this can be overridded as needed.
   */
  virtual void artVisc(FArrayBox&       a_F,
                       const FArrayBox& a_U,
                       const Real&      a_artificialViscosity,
                       const Real&      a_currentTime,
                       const int&       a_dir,
                       const Box&       a_box);

#ifdef CH_USE_HDF5
  virtual void expressions(HDF5HeaderData& a_holder) const
  {
  }
#endif

  /**
     \name Access functions
  */
  /**@{*/

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
  */
  virtual int numConserved() = 0;

  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.
  */
  virtual Vector<string> stateNames() = 0;

  /// Returns true if 4th-order artificial viscosity is defined.
  /**
   */
  virtual bool fourthOrderArtificialViscosityIsDefined() const;

  /// Defines fourth-order artifical viscosity strong shock threshold.
  /**
   */
  virtual void setFourthOrderArtificialViscosityParameter(const Real& M0sq);
 /// Returns fourth-order artifical viscosity strong shock threshold.
  /**
   */
  virtual Real getFourthOrderArtificialViscosityParameter() const;

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
  */
  virtual int numFluxes() = 0;

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

  /// Used to limit the absolute value of a "pressure" difference
  /**
     Return a value that is used by slope flattening to limit (away from
     zero) the absolute value of a slope in the pressureIndex() component
     (slope computation).
   */
  virtual Real smallPressure() = 0;

  /// Component index within the primitive variables of the bulk modulus
  /**
     Return the component index withn the primitive variables for the
     bulk modulus.  Used for slope flattening (slope computation) used
     as a normalization to measure shock strength.
   */
  virtual int bulkModulusIndex() = 0;

 /// Component index within the primitive variables of the density.
  /**
     Return the component index within the primitive variables for the
     density.  Used for fourth-order accurate artificial viscosity.
   */
  virtual int densityIndex();

  /**@}*/

protected:
  // Has this object been defined
  bool m_isDefined;

  // Problem domain and grid spacing
  ProblemDomain m_domain;
  Real          m_dx;

  // Object containing various methods for the Godunov computation
  GodunovUtilities m_util;

  // Flag to use fourth-order accurate artificial viscosity, strong shock
  // threshhold M0sq (looks like the shock Mach number squared).
  bool m_useFourthOrderArtificialViscosity;
  Real m_M0sq;

  // Initial and boundary condition object and has it been set
  PhysIBC* m_bc;
  bool     m_isBCSet;


private:

  // Disallowed for all the usual reasons
  void operator=(const GodunovPhysics&);
  GodunovPhysics(const GodunovPhysics&);
};

#include "NamespaceFooter.H"
#endif
