#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTIONPHYSICS_H_
#define _ADVECTIONPHYSICS_H_

#include "GodunovPhysics.H"
#include "FluxBox.H"

#include "UsingNamespace.H"

/// An GodunovPhysics-derived class for simple advection-diffusion problems
/**
 */
class AdvectionPhysics : public GodunovPhysics
{
public:
  /// Constructor
  /**
   */
  AdvectionPhysics();

  /// Destructor
  /**
   */
  ~AdvectionPhysics();

  /// Set the number of components
  /**
   */
  virtual void setNComp(int a_nComp);

  /// Compute the maximum wave speed
  /**
   */
  virtual Real getMaxWaveSpeed(const FArrayBox& a_U,
                               const Box&       a_box);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to new AdvectionPhysics object with the same definition
     as this object.
  */
  virtual GodunovPhysics* new_godunovPhysics() const;

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
   */
  virtual int numConserved();

  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.
   */
  virtual Vector<string> stateNames();

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the
     number of conserved variables if addition fluxes/face-centered
     quantities are computed.
  */
  virtual int numFluxes();

  /// compute fluxes from primitive values on a face ( advVel*wHalf)
  /** Fluxes are computed as advVel*wHalf
   */
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_WHalf,
                       const int&       a_dir,
                       const Box&       a_box);

  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
  */
  virtual bool isDefined() const;

  /// Number of primitve variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  virtual int numPrimitives();

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
                            const Box&       a_box);

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
                             const Box&       a_box);

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
                               const Box&       a_box);

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

  /// Interval within the primitive variables corresponding to the velocities
  /**
     Return the interval of component indices within the primitive variable
     of the velocities.  Used for slope flattening (slope computation) and
     computing the divergence of the velocity (artificial viscosity).
   */
  virtual Interval velocityInterval();

  /// Component index within the primitive variables of the pressure
  /**
     Return the component index withn the primitive variables for the
     pressure.  Used for slope flattening (slope computation).
   */
  virtual int pressureIndex();

  /// Used to limit the absolute value of a "pressure" difference
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

  /// Set face-centered advection velocity
  void setAdvVelPtr(FluxBox* a_advVelPtr);

  /// Get face-centered advection velocity
  const FluxBox* advectionVelPtr() const;

  /// Set cell-centered advection velocity
  void setCellVelPtr(FArrayBox* a_cellVelPtr);

  /// Get cell-centered advection velocity
  const FArrayBox* cellVelPtr() const;


protected:
  /// face-centered advection velocity
  FluxBox* m_advVelPtr;

  /// cell-centered advection velocity (centered at old time)
  FArrayBox* m_cellVelPtr;

  /// Have the numbers below been defined
  bool m_nCompDefined;

  /// number of conserved variables
  int m_numCons;

  /// number of primitive variables
  int m_numPrim;

  /// number of fluxes
  int m_numFlux;
};

#endif
