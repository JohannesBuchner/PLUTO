#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GODUNOVUTILITIES_H_
#define _GODUNOVUTILITIES_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "NamespaceHeader.H"

class GodunovPhysics;

///
/**
   Utility class for higher-order Godunov methods: slopes, parabolic
   interpolants, limiters. Contains no physics-dependent methods, but
   one of the member functions (PPMFaceValues()) may require a pointer
   to a GodunovPhysics class in order to set boundary slopes.
*/
class GodunovUtilities
{
public:
  /// Constructor
  /**
   */
  GodunovUtilities();

  /// Destructor
  /**
   */
  ~GodunovUtilities();

  /// Define the object
  /**
   */
  void define(const ProblemDomain& a_domain,
              const Real&          a_dx);

  /**
     \name Slopes
  */

  /**@{*/

  /// Compute componentwise van Leer slopes.
  /**
     Given cell averages W, compute van Leer slopes dW on a
     component-by-component basis.
   */
  void vanLeerSlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_numSlopes,
                     const bool&      a_useLimiting,
                     const int&       a_dir,
                     const Box&       a_box);

  /// Compute fourth-order slopes.
  /**
     Given cell averages W and van Leer slopes dWvL, compute fourth-order
     slopes dW4. Limiting is performed in a separate pass.
  */
  void fourthOrderSlopes(FArrayBox&       a_dW4,
                         const FArrayBox& a_W,
                         const FArrayBox& a_dWvL,
                         const int&       a_numSlopes,
                         const int&       a_dir,
                         const Box&       a_box);

  /// Compute slopes (dW- and dW+) using one sided differences
  /**
   */
  void oneSidedDifferences(FArrayBox&       a_dWMinus,
                           FArrayBox&       a_dWPlus,
                           const FArrayBox& a_W,
                           const int&       a_dir,
                           const Box&       a_box);

  /// Compute slopes (dW (center), dW-, and dW+)
  /**
     a_dwCent, a_dwMinus, a_dWPlus, and a_dW all live on
     cell-centered a_entireBox.

     For i in a_centerBox:
     a_dwCent[i] = (a_W[i+1] - a_W[i-1])/2;
     a_dWMinus[i] = a_W[i] - a_W[i-1];
     a_dWPlus[i] = a_W[i+1] - a_W[i].

     For i in a_loBox, set only a_dWPlus[i] = a_W[i+1] - a_W[i]
     and copy it to a_dWCent[i] and a_dWMinus[i].

     For i in a_hiBox, set only a_dWMinus[i] = a_W[i] - a_W[i-1]
     and copy it to a_dWCent[i] and a_dWPlus[i].
   */
  void slopes(FArrayBox&       a_dWCent,
              FArrayBox&       a_dWMinus,
              FArrayBox&       a_dWPlus,
              const FArrayBox& a_W,
              const int&       a_numSlopes,
              const int&       a_dir,
              const Box&       a_loBox,
              const Box&       a_hiBox,
              const Box&       a_centerBox,
              const Box&       a_entireBox,
              const int&       a_hasLo,
              const int&       a_hasHi);

  /**@}*/

  /**
     \name Parabolic interpolants
  */

  /**@{*/

  /// PLM normal predictor.
  /**
     Compute the increments in the characteristic amplitudes.
  */
  void PLMNormalPred(FArrayBox&       a_dWCharMinus,
                     FArrayBox&       a_dWCharPlus,
                     const FArrayBox& a_dWChar,
                     const FArrayBox& a_Lambda,
                     const Real&      a_dtbydx,
                     const Box&       a_box);

  /// PPM face-centered interpolant.
  /**
     Given the cell average a_W, compute fourth-order accurate
     face-centered values WFace on a_box by differentiating the indefinite
     integral. Limiting is performed in a separate pass.
  */
  void PPMFaceValues(FArrayBox&            a_WFace,
                     const FArrayBox&      a_W,
                     const int&            a_numSlopes,
                     const bool&           a_useLimiting,
                     const int&            a_dir,
                     const Box&            a_box,
                     const Real&           a_time,
                     const GodunovPhysics* a_physPtr = NULL);


  /// PPM normal predictor.
  /**
     On input, dW(Minus,Plus), contain the characteristic
     expansions of the differences between the (minus, plus) face values
     and the cell average. On output, dW(Minus,Plus) contain the
     characteristic amplitudes of the corrections required to compute
     the normal predictor.
  */
  void PPMNormalPred(FArrayBox&       a_dWMinus,
                     FArrayBox&       a_dWPlus,
                     const FArrayBox& a_Lambda,
                     const Real&      a_dtbydx,
                     const int&       a_numSlopes,
                     const Box&       a_box);

  /**@}*/

  /**
     \name Limiters
  */

  /**@{*/

  /// van Leer slope limiter.
  /**
     On input, dW contains the centered, unlimited slopes, and
     dW(Minus,Plus) contain the one-sided slopes from the minus, plus sides.
     On output, dW contains the limited slopes.
     slopes dW4. Limiting is performed in a separate pass.
  */
  void slopeLimiter(FArrayBox&       a_dW,
                    const FArrayBox& a_dWMinus,
                    const FArrayBox& a_dWPlus,
                    const int&       a_numSlopes,
                    const Box&       a_box);

  /// extremum-preserving van Leer slope limiter.
  /**
     On input, dW contains the centered, unlimited slopes, and
     dW(Minus,Plus) contain the one-sided slopes from the minus, plus sides.
     On output, dW contains the limited slopes.
     slopes dW4. Limiting is performed in a separate pass.
  */
  void slopeLimiterExtPreserving(FArrayBox&       a_dW,
                                 const FArrayBox& a_dWMinus,
                                 const FArrayBox& a_dWPlus,
                                 const int&       a_numSlopes,
                                 const Box&       a_box,
                                 const int&       a_dir);

  /// PPM Limiter.
  /**
     a_dWMinus and a_dWPlus are the differences between the face values
     on the minus and plus sides of cells and the average in the cell.
     That is,
     a_dWMinus[i] = WFace[i - e/2] - a_W[i]
     a_dWPlus[i] = WFace[i + e/2] - a_W[i]
     where e is the unit vector in dimension a_dir.
     The PPM limiter is applied to these values to obtain a monotone
     interpolant in the cell.
     The function returns the limited a_dWMinus and a_dWPlus on a_box.
     petermc, 4 Sep 2008:  included a_W in argument list

     If m_highOrderLimiter,
     then need a_dWMinus and a_dWPlus on a_box,
     and need a_W on on a_box grown by 3 in dimension a_dir.
     Returns limited a_dWMinus and a_dWPlus on a_box.
  */
  void PPMLimiter(FArrayBox& a_dWMinus,
                  FArrayBox& a_dWPlus,
                  const FArrayBox& a_W,
                  const int& a_numSlopes,
                  const int& a_dir,
                  const Box& a_box);

  /// Set whether to use high-order limiter.
  void highOrderLimiter(bool a_highOrderLimiter);

  /**@}*/

  /**
     \name Dissipation mechanisms
  */

  /**@{*/

  /// Compute the slope flattening coefficients
  /**
     Compute the slope flattening coefficients, a_flattening, using the
     primitive variables, a_W, within a_box.
  */
  void computeFlattening(/// flattening coeffs, 1 component on a_box
                         FArrayBox&       a_flattening,
                         /// primitive variables, on a_box grown by 3 within m_domain
                         const FArrayBox& a_W,
                         /// interval of a_W with velocity; length SpaceDim
                         const Interval&  a_velInt,
                         /// index of a_W with pressure
                         const int&       a_presInd,
                         /// minimum pressure to use in ratio
                         const Real&      a_smallPres,
                         /// index of a_W with bulk modulus
                         const int&       a_bulkModulusInd,
                         /// cell-centered box on which a_flattening lives
                         const Box&       a_box);

  /// Apply the flattening to slopes.
  /**
     Multiply every component of a_dW by a_flat, on a_box.
   */
  void applyFlattening(/// slopes to be flattened, on at least a_box
                       FArrayBox&       a_dW,
                       /// cell-centered flattening coefficients, on at least a_box
                       const FArrayBox& a_flat,
                       /// cell-centered box on which to flatten
                       const Box&       a_box);

  /// Compute face-centered artificial viscosity flux.
  /**
     Increments face-centered flux in the a_dir direction with quadratic
     artificial viscosity.
  */
  void artificialViscosity(FArrayBox&       a_F,
                           const FArrayBox& a_U,
                           const FArrayBox& a_divVel,
                           const Real&      a_scale,
                           const int&       a_dir,
                           const Box&       a_box);

  /**@}*/


  /**
     \name Divergence
  */

  /**@{*/

  /// Compute face-centered velocity divergence.
  /**
     Returns face-centered velocity divergence on faces in the
     direction a_dir. The velocities are the components a_velInterval
     of a_W.
  */
  void divVel(FArrayBox&       a_divVel,
              const FArrayBox& a_W,
              const Interval&  a_velInt,
              const int&       a_dir,
              const Box&       a_box);

  /// Compute high-order face-centered velocity divergence for artificial viscosity
  /**
     Computes a face-centered nonlinear function of the divergence
     suitable for use as an artificial viscosity coefficient for
     a fourth-order method.

     F_visc = -dx*K*dU/dx
     K = max(-dx*div(u)*min(1,h^2*|div(u)/(c^2)/a_M0sq|),0)
  */
  void divVelHO(FArrayBox&       a_divVel,
                const FArrayBox& a_W,
                const int&       a_dir,
                const Box&       a_box,
                GodunovPhysics*  a_physPtr);

  /**@}*/

  /**
     \name Deconvolution and convolution
  */

  /**@{*/

  /// Deconvolve from cell-centered to cell-averaged data (or vice versa)
  /**
     Deconvolves (or convolves, if a_sign = -1) by the formula

     a_avgFab += a_sign * laplacian(a_cellFab) * m_dx^2 / 24, on a_cellFab.box()
  */
  void deconvolve(FArrayBox&        a_avgFab,
                  const FArrayBox&  a_cellFab,
                  int               a_sign = 1);

  /// Deconvolve from cell-centered to cell-averaged data (or vice versa) on a specified box
  /**
     Deconvolves (or convolves, if a_sign = -1) by the formula

     a_avgFab += a_sign * laplacian(a_cellFab) * m_dx^2 / 24, on a_box
  */
  void deconvolve(FArrayBox&        a_avgFab,
                  const FArrayBox&  a_cellFab,
                  const Box&        a_box,
                  int               a_sign = 1);

  /// Deconvolve from face-centered to face-averaged data (or vice versa)
  /**
     Deconvolves (or convolves, if a_sign = -1) by the formula

     a_avgFlux += a_sign * laplacian(a_cellFlux) * m_dx^2 / 24, on faces of a_cellFlux.box()
  */
  void deconvolveFace(FluxBox&        a_avgFlux,
                      const FluxBox&  a_cellFlux,
                      int             a_sign = 1);

  /// Deconvolve from face-centered to face-averaged data (or vice versa) on a specified box
  /**
     Deconvolves (or convolves, if a_sign = -1) by the formula

     a_avgFlux += a_sign * laplacian(a_cellFlux) * m_dx^2 / 24, on faces of a_box
  */
  void deconvolveFace(FluxBox&        a_avgFlux,
                      const FluxBox&  a_cellFlux,
                      const Box&      a_box,
                      int             a_sign = 1);

  /**@}*/

protected:
  // Problem domain and grid spacing
  ProblemDomain m_domain;
  Real          m_dx;

  // Has this object been defined
  bool          m_isDefined;

  // Use a high-order limiter?  (default false)
  bool          m_highOrderLimiter;

private:
  // We may allow copying and assignment later.
  // Disallowed for all the usual reasons
  void operator=(const GodunovUtilities&);
  GodunovUtilities(const GodunovUtilities&);
};

#include "NamespaceFooter.H"
#endif
