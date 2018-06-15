#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCHADVECTION_H_
#define _PATCHADVECTION_H_

#include "OldPatchGodunov.H"
#include "FluxBox.H"
#include "UsingNamespace.H"

/// An OldPatchGodunov-derived class for simple advection-diffusion problems
/**
 */
class PatchAdvection : public OldPatchGodunov
{
public:
  ///
  /**
   */
  PatchAdvection();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to new OldPatchGodunov object with its initial and boundary
     condtions, slope parameters, and artificial viscosity information defined.
  */
  virtual OldPatchGodunov* new_patchGodunov() const;


  /// Compute the maximum wave speed
  /**
   */
  virtual Real getMaxWaveSpeed(const FArrayBox& a_U,
                               const Box&       a_box);

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
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
   */
  virtual int numFluxes();

  /// set face-centered advection velocity for this patch
  void setAdvVelPtr(FluxBox* a_advVelPtr);

  /// access face-centered advection velocities
  const FluxBox* advectionVelPtr() const;

  /// set cell-centered velocity field for this patch
  void setCellVelPtr(FArrayBox* a_cellVelPtr);

  /// access cell-centered velicity field pointer
  const FArrayBox* cellVelPtr() const;

  /// set the number of conserved (and primitive, and flux) variables
  /**
      Note that for this simple tracing, conservative == primitive
   */
  virtual void setNumVar(int a_numVar);


protected:
  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  virtual int numPrimitives();

  ///  Number of primitive variables for which slopes are computed
  /**
     Return the number of primitive variables for which slopes are computed.
     Only slopes corresponding to primitive variables in the interval 0 to
     numSlopes() - 1 are computed and only primitive variables in that interval
     are updated using the slopes.
   */
  virtual int numSlopes();


  /// Compute the primitive variables from the conserved variables within a_box
  /**
   */
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box);

  /// Compute the conserved variables from the primitive variables within a_box
  /**
   */
  virtual void primToCons(FArrayBox&       a_U,
                          const FArrayBox& a_W,
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
                          const Box&       a_box);

  /// increment the primitive variables by a source term
  /**
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
                       const Box&       a_box);

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
                          const Box&       a_box);


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

  /// Component index within the primitive variables of the bulk modulus
  /**
     Return the component index withn the primitive variables for the
     bulk modulus.  Used for slope flattening (slope computation) used
     as a normalization to measure shock strength.
   */
  virtual int bulkModulusIndex();

  /// face-centered advection velocity, set through the advVel() intervace
  FluxBox* m_advVelPtr;

  /// cell-centered velocity field (centered at old time)
  FArrayBox* m_cellVelPtr;

  /// number of conserved variables
  /** For this simple adection tracing, we assume that conservative and
      primitive are the same
  */
  int m_numCons;

  /// number of primitive variables
  int m_numPrim;

  /// number of fluxes
  int m_numFlux;

  /// number of variabls in slope computation
  int m_numSlope;

};

#endif









