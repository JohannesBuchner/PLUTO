#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PATCHGODUNOV_H_
#define _PATCHGODUNOV_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "FluxBox.H"

#include "GodunovPhysics.H"
#include "GodunovUtilities.H"
#include "NamespaceHeader.H"

class HDF5HeaderData;

///
/**
   The base class PatchGodunov provides an implementation of a second-order,
   unsplit Godunov method acting on a single grid/patch.  PatchGodunov
   provides an interface to the level integrator, LevelGodunov, which manages
   the entire level and flux corrections (via flux registers) to the coarser
   and finer levels.
 */
class PatchGodunov
{
public:
  /// Constructor
  /**
   */
  PatchGodunov();

  /// Destructor
  /**
   */
  virtual ~PatchGodunov();

  /// Define the object
  /**
   */
  virtual void define(/// problem domain
                      const ProblemDomain&        a_domain,
                      /// grid spacing
                      const Real&                 a_dx,
                      /// physics class
                      const GodunovPhysics* const a_gdnvPhysicsPtr,
                      /// order of the normal predictorr: 0 for CTU, 1 for PLM, 2 for PPM
                      const int&                  a_normalPredOrder,
                      /// whether to use 4th-order slopes
                      const bool&                 a_useFourthOrderSlopes,
                      /// whether to apply slope limiting to primitive variables
                      const bool&                 a_usePrimLimiting,
                      /// whether to apply slope limiting to characteristic variables
                      const bool&                 a_useCharLimiting,
                      /// whether to apply slope flattening
                      const bool&                 a_useFlattening,
                      /// whether to apply artificial viscosity
                      const bool&                 a_useArtificialViscosity,
                      /// artificial viscosity coefficient
                      const Real&                 a_artificialViscosity);

  /// Set the current time before calling updateState(), computeWHalf(), computeUpdate().
  /**
   */
  virtual void setCurrentTime(const Real& a_currentTime);

  /// Set the current box before calling updateState(), computeWHalf(), computeUpdate().
  /**
   */
  virtual void setCurrentBox(const Box& a_currentBox);

  /// Update the conserved variables and return the fluxes used to do this.
  /**
   */
  virtual void updateState(/// conserved variables, updated in this routine
                           FArrayBox&       a_U,
                           /// fluxes from 2nd-order unsplit Godunov method, returned
                           FluxBox&         a_F,
                           /// maximum wave speed, returned
                           Real&            a_maxWaveSpeed,
                           /// source terms, or null constructed if none
                           const FArrayBox& a_S,
                           /// time step
                           const Real&      a_dt,
                           /// Box of a_U
                           const Box&       a_box);

  virtual void updateState(/// conserved variables, updated in this routine
                           FArrayBox&       a_U,
                           /// fluxes from 2nd-order unsplit Godunov method, returned
                           FluxBox&         a_F,
                           /// primitive variables extrapolated to cell faces and a half time step
                           FluxBox&         a_wHalf,
                           /// maximum wave speed, returned
                           Real&            a_maxWaveSpeed,
                           /// source terms, or null constructed if none
                           const FArrayBox& a_S,
                           /// time step
                           const Real&      a_dt,
                           /// Box of a_U
                           const Box&       a_box);

  /// Compute the time-centered values of the primitive variables on cell faces.
  /**
     This API is used in cases where some operation over the whole
     level must be performed on the face-centered variables prior to the
     final difference update. Examples include incompressible flow and MHD,
     in which it is necessary to compute the projection of a face-centered
     vector field on its divergence-free part. To complete the differencing,
     it is necessary to call the member function computeUpdate().
   */
  virtual void computeWHalf(/// primitive variables extrapolated to cell faces and a half time step
                            FluxBox&         a_WHalf,
                            /// conserved variables at beginning of time step
                            const FArrayBox& a_U,
                            /// source terms, or null constructed if none
                            const FArrayBox& a_S,
                            /// time step
                            const Real&      a_dt,
                            /// Box of a_U
                            const Box&       a_box);

  /// Compute the increment in the conserved variables from face variables.
  /**
   Compute dU = dt*dU/dt, the change in the conserved variables over
   the time step. The fluxes are returned are suitable for use in refluxing.
   */
  virtual void computeUpdate(/// change in conserved variables over the time step
                             FArrayBox&       a_dU,
                             /// fluxes for refluxing
                             FluxBox&         a_F,
                             /// conserved variables at beginning of time step
                             const FArrayBox& a_U,
                             /// primitive variables extrapolated to cell faces and a half time step
                             const FluxBox&   a_WHalf,
                             /// time step
                             const Real&      a_dt,
                             /// Box of a_U
                             const Box&       a_box);

  /// Compute fluxes from primitive variables on faces
  /**
   */
  void computeFluxes(FluxBox&        a_F,
                     const FluxBox&  a_WHalf,
                     const Box&      a_box);

  /// CTU normal predictor
  /**
     Compute the increments in the characteristic amplitudes using CTU
     (for CTU, increments are zero)
   */
  void CTUNormalPred(FArrayBox&       a_WMinus,
                     FArrayBox&       a_WPlus,
                     const Real&      a_dt,
                     const Real&      a_dx,
                     const FArrayBox& a_W,
                     const FArrayBox& a_flat,
                     const int&       a_dir,
                     const Box&       a_box);

  /// PLM normal predictor
  /**
     Compute the increments in the characteristic amplitudes using PLM
   */
  void PLMNormalPred(FArrayBox&       a_WMinus,
                     FArrayBox&       a_WPlus,
                     const Real&      a_dt,
                     const Real&      a_dx,
                     const FArrayBox& a_W,
                     const FArrayBox& a_flat,
                     const int&       a_dir,
                     const Box&       a_box);

  /// PPM normal predictor
  /**
     Compute the increments in the characteristic amplitudes using PPM
   */
  void PPMNormalPred(FArrayBox&       a_WMinus,
                     FArrayBox&       a_WPlus,
                     const Real&      a_dt,
                     const Real&      a_dx,
                     const FArrayBox& a_W,
                     const FArrayBox& a_flat,
                     const int&       a_dir,
                     const Box&       a_box);

  /// Accessor to get at the GodunovPhysics object
  /**
   */
  GodunovPhysics* getGodunovPhysicsPtr();

  //! Get the problem domain for this object.
  const ProblemDomain& problemDomain() const;

  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
   */
  virtual bool isDefined() const;

  /// Set whether to use high-order limiter.
  void highOrderLimiter(bool a_highOrderLimiter);

  //! Returns the grid spacing at this level.
  Real dx() const;

#ifdef CH_USE_HDF5
  virtual void expressions(HDF5HeaderData& a_holder) const
  {
  }
#endif

protected:
  // Problem domain and grid spacing
  ProblemDomain m_domain;
  Real          m_dx;

  // Object providing physics information
  GodunovPhysics*  m_gdnvPhysics;

  // Object containing various methods for the Godunov computation
  GodunovUtilities m_util;

  // Order of the normal predictor (1 -> PLM, 2-> PPM)
  int m_normalPredOrder;

  // Use 4th order slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  // Do slope limiting in the primitive or characteristic variables, respect.
  bool m_usePrimLimiting;
  bool m_useCharLimiting;

  // Do slope flattening - MUST BE USING 4th order slopes
  bool m_useFlattening;

  // Apply artificial viscosity of a set value
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;

  // Current time and has it been set
  Real m_currentTime;
  bool m_isCurrentTimeSet;

  // Current box and has it been set
  Box  m_currentBox;
  bool m_isCurrentBoxSet;

  // Has this object been defined
  bool m_isDefined;

  // Use a high-order limiter?  (default false)
  bool m_highOrderLimiter;

private:

  // Disallowed for all the usual reasons
  void operator=(const PatchGodunov&);
  PatchGodunov(const PatchGodunov&);
};

#include "NamespaceFooter.H"
#endif
