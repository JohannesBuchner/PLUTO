#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTPHYSICS_H_
#define _ADVECTPHYSICS_H_

#include "GodunovPhysics.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

/// A class derived from GodunovPhysics for simple advection-diffusion problems
/**
 */
class AdvectPhysics : public GodunovPhysics
{
public:
  /// Constructor
  /**
   */
  AdvectPhysics();

  /// Destructor
  /**
   */
  ~AdvectPhysics();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to new AdvectPhysics object with the same definition
     as this object.
  */
  virtual GodunovPhysics* new_godunovPhysics() const;

  /// Compute the maximum wave speed
  /**
   */
  virtual Real getMaxWaveSpeed(const FArrayBox& a_U,
                               const Box&       a_box);

  /// COMPUTE fluxes from primitive values on a face ( advVel*wHalf)
  /** Fluxes are computed as advVel*wHalf
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
                            const Box&       a_box)
  {

    // This is the identity mapping for advection - do nothing
  }

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
                             const Box&       a_box)
  {

    // This is the identity mapping for advection - do nothing
  }

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
                          const Box&       a_box);

  /// Add to (increment) the source terms given the current state
  /**
     On input, a_S contains the current source terms.  On output, a_S has
     had any additional source terms (based on the current state, a_W)
     added to it.  This should all be done on the region defined by a_box.
  */
  virtual void incrementSource(FArrayBox&       a_S,
                               const FArrayBox& a_W,
                               const Box&       a_box)
  {

    // The source term does not explicitly depend on the current primitive state
  }

  /// Compute the solution to the Riemann problem.
  /**
     Given input left and right states in a direction, a_dir, compute a
     Riemann problem and generate fluxes at the faces within a_box.
  */
  virtual void riemann(FArrayBox&       a_WGdnv,
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
  virtual void quasilinearUpdate(FArrayBox&       a_AdWdx,
                                 const FArrayBox& a_WHalf,
                                 const FArrayBox& a_W,
                                 const Real&      a_scale,
                                 const int&       a_dir,
                                 const Box&       a_box);

  /// Compute primitive variables from conserved variables.
  /**
   */
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box);

  /// Set cell-centered and face centered advection velocity
  void setVelocities(FArrayBox* a_celVelPtr,
                     FluxBox*   a_advVelPtr)
  {
    m_isVelSet = true;
    m_advVelPtr  = a_advVelPtr;
    m_cellVelPtr = a_celVelPtr;
  }

  /**
     \name Access functions
  */
  /**@{*/

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
   */
  virtual int numConserved()
  {
    return 1;
  }

  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.
   */
  virtual Vector<string> stateNames()
  {
    Vector<string> retval(1, string("scalar"));

    return retval;
  }

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the
     number of conserved variables if addition fluxes/face-centered
     quantities are computed.
  */
  virtual int numFluxes()
  {
    return 1;
  }

  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
  */
  virtual bool isDefined() const
  {
    return m_isDefined;
  }

  /// Number of primitve variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  virtual int numPrimitives()
  {
    return 1;
  }

  /// Interval within the primitive variables corresponding to the velocities
  /**
     Return the interval of component indices within the primitive variable
     of the velocities.  Used for slope flattening (slope computation) and
     computing the divergence of the velocity (artificial viscosity).
   */
  virtual Interval velocityInterval()
  {
    MayDay::Error("AdvectPhysics::velocityInterval - not defined");

    Interval retval(-1,-1);
    return retval;
  }

  /// Component index within the primitive variables of the pressure
  /**
     Return the component index withn the primitive variables for the
     pressure.  Used for slope flattening (slope computation).
     // Component index within the primitive variables of the pressure
     // since this doesn't apply to this set of equations, return a bogus value
   */
  virtual int pressureIndex()
  {
    MayDay::Error("AdvectPhysics::pressureIndex - not defined");

    return -1;
  }

  /// Used to limit the absolute value of a "pressure" difference
  /**
     Return a value that is used by slope flattening to limit (away from
     zero) the absolute value of a slope in the pressureIndex() component
     (slope computation).
     // Used to limit the absolute value of a "pressure" difference
     // since this doesn't apply to this set of equations, return a bogus value
   */
  virtual Real smallPressure()
  {
    MayDay::Error("AdvectPhysics::smallPressure - not defined");

    return -1.0;
  }

  /// Component index within the primitive variables of the bulk modulus
  /**
     Return the component index withn the primitive variables for the
     bulk modulus.  Used for slope flattening (slope computation) used
     as a normalization to measure shock strength.
   */
  virtual int bulkModulusIndex()
  {
    MayDay::Error("AdvectPhysics::bulkModulusIndex - not defined");

    return -1;
  }

  /// Component index within the primitive variables of the density
  /**
     Return the component index withn the primitive variables for the
     density. Not defined for AdvectPhysics.
   */
  virtual int densityIndex()
  {
    MayDay::Error("AdvectPhysics::bulkModulusIndex - not defined");

    return -1;
  }

  /**@}*/

protected:
  bool m_isVelSet;
  /// face-centered advection velocity
  FluxBox* m_advVelPtr;

  /// cell-centered advection velocity (centered at old time)
  FArrayBox* m_cellVelPtr;

};

#include "NamespaceFooter.H"
#endif
