#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POLYTROPICPHYSICS_H_
#define _POLYTROPICPHYSICS_H_

#include <string>
using std::string;

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "CH_HDF5.H"

#include "GodunovPhysics.H"

#include "NamespaceHeader.H"

///
/**
   The base class PolytropicPhysics provides the physics-dependent components
   for a higher-order method for a single patch: characteristic
   analysis, Riemann solver, quasilinear update, conservative update,
   and transformations between conserved, primitive, and flux variables.
   This class is essentially pure; i.e., all of its member functions
   are; and the ones that have default implementations are ones
   that are optionally defined; i.e., the default definition is to send
   an error message. Physics-dependent versions of this class that are
   required in real applications are derived from this class by inheritance.
*/
class PolytropicPhysics: public GodunovPhysics
{
public:
  /// Constructor
  /**
   */
  PolytropicPhysics(const Real& a_smallPressure);

  /// Destructor
  /**
   */
  ~PolytropicPhysics();

  /// Compute the maximum wave speed
  /**
   */
  Real getMaxWaveSpeed(const FArrayBox& a_U,
                       const Box&       a_box);

  /// Compute the speed of sound
  /**
   */
  void soundSpeed(FArrayBox& a_speed,
                  const FArrayBox& a_U,
                  const Box&       a_box);

  /// Object factory for this class
  /**
   */
  virtual GodunovPhysics* new_godunovPhysics() const;

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
  */
  int numConserved();

  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numConserved()-1.
  */
  Vector<string> stateNames();

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
  */
  int numFluxes();

  /// Component index within the primitive variables of the density.
  /**
     Return the component index within the primitive variables for the
     density.  Used for fourth-order accurate artificial viscosity.
   */
  int densityIndex();

  /// Compute a flux from primitive variable values on a face
  /**
   */
  void getFlux(FArrayBox&       a_flux,
               const FArrayBox& a_whalf,
               const int&       a_dir,
               const Box&       a_box);

  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
  */
  int numPrimitives();

  /// Transform a_dW from primitive to characteristic variables
  /**
     On input, a_dW contains the increments of the primitive variables. On
     output, it contains the increments in the characteristic variables.

     IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
     smallest eigenvalue first, the largest eigenvalue last, and orders the
     characteristic variables accordingly.
  */
  void charAnalysis(FArrayBox&       a_dW,
                    const FArrayBox& a_W,
                    const int&       a_dir,
                    const Box&       a_box);

  /// Transform a_dW from characteristic to primitive variables
  /**
     On input, a_dW contains the increments of the characteristic variables.
     On output, it contains the increments in the primitive variables.

     IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
     smallest eigenvalue first, the largest eigenvalue last, and orders the
     characteristic variables accordingly.
  */
  void charSynthesis(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Box&       a_box);

  /// Compute the characteristic values (eigenvalues)
  /**
     Compute the characteristic values (eigenvalues)

     IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
     smallest eigenvalue first, the largest eigenvalue last, and orders the
     characteristic variables accordingly.
   */
  void charValues(FArrayBox&       a_lambda,
                  const FArrayBox& a_W,
                  const int&       a_dir,
                  const Box&       a_box);

  /// Add to (increment) the source terms given the current state
  /**
     On input, a_S contains the current source terms.  On output, a_S has
     had any additional source terms (based on the current state, a_W)
     added to it.  This should all be done on the region defined by a_box.
  */
  void incrementSource(FArrayBox&       a_S,
                       const FArrayBox& a_W,
                       const Box&       a_box);

  /// Compute the solution to the Riemann problem.
  /**
     Given input left and right states in a direction, a_dir, compute a
     Riemann problem and generate fluxes at the faces within a_box.
  */
  void riemann(FArrayBox&       a_WStar,
               const FArrayBox& a_WLeft,
               const FArrayBox& a_WRight,
               const FArrayBox& a_W,
               const Real&      a_time,
               const int&       a_dir,
               const Box&       a_box);

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
                              const Box&       a_box);

  /// Compute the quasilinear update A*dW/dx.
  /**
   */
  void quasilinearUpdate(FArrayBox&       a_dWdx,
                         const FArrayBox& a_WHalf,
                         const FArrayBox& a_W,
                         const Real&      a_scale,
                         const int&       a_dir,
                         const Box&       a_box);

  /// Compute primitive variables from conserved variables.
  /**
   */
  void consToPrim(FArrayBox&       a_W,
                  const FArrayBox& a_U,
                  const Box&       a_box);

  /// Interval within the primitive variables corresponding to the velocities
  /**
     Return the interval of component indices within the primitive variable
     of the velocities.  Used for slope flattening (slope computation) and
     computing the divergence of the velocity (artificial viscosity).
   */
  virtual Interval velocityInterval();

  /// Interval within the flux variables corresponding to vector flux
  virtual Interval vectorFluxInterval();

  /// Component index within the primitive variables of the pressure
  /**
     Return the component index withn the primitive variables for the
     pressure.  Used for slope flattening (slope computation).
   */
  virtual int pressureIndex();

  /// Used to limit the absolute value of a "pressure" difference (away from zero)
  /**
     Return a value that is used by slope flattening to limit (away from
     zero) the absolute value of a slope in the pressureIndex() component
     (slope computation).
   */
  virtual Real smallPressure();

  /// Component index within the primitive variables of the bulk modulus
  /**
     Return the component index withn the primitive variables for the
     bulk modulus.  Used for slope flattening (slope computation) used
     as a normalization to measure shock strength.
   */
  virtual int bulkModulusIndex();

#ifdef CH_USE_HDF5
  virtual void expressions(HDF5HeaderData& a_holder) const;
#endif

protected:
  // Used to limit the absolute value of a "pressure" difference (away from zero)
  Real m_smallPressure;


private:

  // Disallowed for all the usual reasons
  void operator=(const PolytropicPhysics&);
  PolytropicPhysics(const PolytropicPhysics&);
};

#include "NamespaceFooter.H"

#endif
