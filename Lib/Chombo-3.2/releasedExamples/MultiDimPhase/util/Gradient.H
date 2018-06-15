#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GRADIENT_H_
#define _GRADIENT_H_

// uses the "new" gradient, which handles corner cells better for
// transverse gradients, at a cost of more complication, along with
// (most likely) worsened scaling in parallel

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "QuadCFInterp.H"

#include "NamespaceHeader.H"

/// Class to encapsulate Gradient functions (both CC and face-centered)
class Gradient
{

public:

  /// computes cell-centered level-operator gradient of cell-centered phi
  /** if phiCrse != NULL, does coarse-fine boundary
      conditions for phi (quadratic interpolation) */
  static void levelGradientCC(
                              ///
                              LevelData<FArrayBox>& a_grad,
                              ///
                              LevelData<FArrayBox>& a_phi,
                              ///
                              const LevelData<FArrayBox>* a_phiCrsePtr,
                              ///
                              const Real a_dx,
                              ///
                              const int a_nRefCrse,
                              ///
                              const ProblemDomain& a_dProblem);

  ///computes cell-centered level-operator gradient of cell-centered phi
  /** if phiCrse != NULL, does coarse-fine boundary
      conditions for phi (quadratic interpolation).  This (deprecated)
      interface uses a Box instead of a ProblemDomain */
  static void levelGradientCC(
                              ///
                              LevelData<FArrayBox>& a_grad,
                              ///
                              LevelData<FArrayBox>& a_phi,
                              ///
                              const LevelData<FArrayBox>* a_phiCrsePtr,
                              ///
                              const Real a_dx,
                              ///
                              const int a_nRefCrse,
                              ///
                              const Box& a_dProblem);

  /// computes cell-centered level-operator gradient of cell-centered phi
  /** if phiCrse != NULL, does coarse-fine boundary
      conditions for phi (quadratic C/F interpolation); predefined
      QuadCFInterp object is passed in for efficiency */
  static void levelGradientCC(
                              ///
                              LevelData<FArrayBox>& a_grad,
                              ///
                              LevelData<FArrayBox>& a_phi,
                              ///
                              const LevelData<FArrayBox>* a_phiCrsePtr,
                              ///
                              const Real a_dx,
                              ///
                              const LayoutData<IntVectSet>& a_gridIVS,
                              ///
                              QuadCFInterp& a_cfInterp);

  /// computes cell-centered, level-operator gradient of cell-centered phi
  /** in this case, assume that all relevant
      BC's (coarse-fine and physical) have already been set, so
      phi can be a const variable */
  static void levelGradientCC(
                              ///
                              LevelData<FArrayBox>& a_grad,
                              ///
                              const LevelData<FArrayBox>& a_phi,
                              ///
                              const Real a_dx);

  /// computes cell-centered composite gradient of cell centered phi
  /**
      Uses same coarse-level C/F BC's as LevelGradientCC, if phiFinePtr !=
      NULL, then also uses one-sided differencing to compute gradient
      on coarse side of coarse-fine interface. */
  static void compGradientCC(
                             ///
                             LevelData<FArrayBox>& a_grad,
                             ///
                             LevelData<FArrayBox>& a_phi,
                             ///
                             const LevelData<FArrayBox>* a_phiCrsePtr,
                             ///
                             const LevelData<FArrayBox>* a_phiFinePtr,
                             ///
                             const Real a_dx,
                             ///
                             const int a_nRefCrse,
                             ///
                             const int a_nRefFine,
                             ///
                             const ProblemDomain& a_dProblem);

  /// computes cell-centered composite gradient of cell centered phi
  /**
      uses same coarse-level C/F BC's as LevelGradientCC, if phiFinePtr !=
      NULL, then also uses one-sided differencing to compute gradient
      on coarse side of coarse-fine interface. This (deprecated)
      interface uses a Box instead of a ProblemDomain */
  static void compGradientCC(
                             ///
                             LevelData<FArrayBox>& a_grad,
                             ///
                             LevelData<FArrayBox>& a_phi,
                             ///
                             const LevelData<FArrayBox>* a_phiCrsePtr,
                             ///
                             const LevelData<FArrayBox>* a_phiFinePtr,
                             ///
                             const Real a_dx,
                             ///
                             const int a_nRefCrse,
                             ///
                             const int a_nRefFine,
                             ///
                             const Box& a_dProblem);

  /// computes cell-centered composite gradient of cell-centered phi
  /**
      uses same coarse-level C/F BC's as LevelGradientCC, if phiFinePtr !=
      NULL, then also uses one-sided differencing to compute gradient
      on coarse side of coarse-fine interface. A predefined QuadCFInterp
      is also passed in for efficiency.  Note that no fine QuadCFInterp
      is necessary because we use one-sided differencing for coarse
      cells adjacent to finer-level coarse-fine interfaces. Note also
      that gradient is only really defined in valid regions of grids. */
  static void compGradientCC(
                             ///
                             LevelData<FArrayBox>& a_grad,
                             ///
                             LevelData<FArrayBox>& a_phi,
                             ///
                             const LevelData<FArrayBox>* a_phiCrsePtr,
                             ///
                             const LevelData<FArrayBox>* a_phiFinePtr,
                             ///
                             const Real a_dx,
                             ///
                             const int a_nRefFine,
                             ///
                             const ProblemDomain& a_dProblem,
                             ///
                             const LayoutData<IntVectSet>& a_gridIVS,
                             ///
                             QuadCFInterp& a_cfInterpCrse);

  /// computes cell-centered composite gradient of cell-centered phi
  /**
      uses same coarse-level C/F BC's as LevelGradientCC, if phiFinePtr !=
      NULL, then also uses one-sided differencing to compute gradient
      on coarse side of coarse-fine interface. A predefined QuadCFInterp
      is also passed in for efficiency.  Note that no fine QuadCFInterp
      is necessary because we use one-sided differencing for coarse
      cells adjacent to finer-level coarse-fine interfaces. Note also
      that gradient is only really defined in valid regions of grids.
      This (deprecated) interface uses a Box instead of a ProblemDomain */
  static void compGradientCC(
                             ///
                             LevelData<FArrayBox>& a_grad,
                             ///
                             LevelData<FArrayBox>& a_phi,
                             ///
                             const LevelData<FArrayBox>* a_phiCrsePtr,
                             ///
                             const LevelData<FArrayBox>* a_phiFinePtr,
                             ///
                             const Real a_dx,
                             ///
                             const int a_nRefFine,
                             ///
                             const Box& a_dProblem,
                             ///
                             const LayoutData<IntVectSet>& a_gridIVS,
                             ///
                             QuadCFInterp& a_cfInterpCrse);

  /// computes cell-centered composite gradient of cell-centered phi
  /**
      this one assumes that all ghost-cell values have already been set;
      if phiFinePtr != NULL, then also uses one-sided differencing to compute
      gradient on coarse side of corarse-fine interface.  note that gradient
      is only really defined in valid regions of grids. */
  static void compGradientCC(
                             ///
                             LevelData<FArrayBox>& a_Grad,
                             ///
                             const LevelData<FArrayBox>& a_phi,
                             ///
                             const LevelData<FArrayBox>* a_phiFinePtr,
                             ///
                             const Real a_dx,
                             ///
                             const int a_nRefFine,
                             ///
                             const ProblemDomain& a_dProblem);

  /// computes cell-centered composite gradient of cell-centered phi
  /**  this one assumes that all ghost-cell values have already been set;
       if phiFinePtr != NULL, then also uses one-sided differencing to
       compute gradient on coarse side of corarse-fine interface.
       note that gradient is only really defined in valid regions of
       grids. This (deprecated) interface uses a Box instead of a
       ProblemDomain */
  static void compGradientCC(
                             ///
                             LevelData<FArrayBox>& a_Grad,
                             ///
                             const LevelData<FArrayBox>& a_phi,
                             ///
                             const LevelData<FArrayBox>* a_phiFinePtr,
                             ///
                             const Real a_dx,
                             ///
                             const int a_nRefFine,
                             ///
                             const Box& a_dProblem);

  /// computes edge-centered level-operator gradient of cell-centered phi
  /** if phiCrsePtr != NULL, does quadratic interpolation
      to compute coarse-fine boundary conditions for phi */
  static void levelGradientMAC(
                               ///
                               LevelData<FluxBox>& a_edgeGrad,
                               ///
                               LevelData<FArrayBox>& a_phi,
                               ///
                               const LevelData<FArrayBox>* a_phiCrsePtr,
                               ///
                               const Real a_dx,
                               ///
                               const int a_nRefCrse,
                               ///
                               const ProblemDomain& a_dProblem);

  /// computes edge-centered level-operator gradient of cell-centered phi
  /** if phiCrsePtr != NULL, does quadratic interpolation
      to compute coarse-fine boundary conditions for phi. This (deprecated)
      interface uses a Box instead of a ProblemDomain */
  static void levelGradientMAC(
                               ///
                               LevelData<FluxBox>& a_edgeGrad,
                               ///
                               LevelData<FArrayBox>& a_phi,
                               ///
                               const LevelData<FArrayBox>* a_phiCrsePtr,
                               ///
                               const Real a_dx,
                               ///
                               const int a_nRefCrse,
                               ///
                               const Box& a_dProblem);

  /// computes edge-centered level-operator gradient of cell-centered phi
  /**
     if phiCrsePtr != NULL, does quadratic interpolation
     to compute coarse-fine boundary conditions for phi.  A predefined
     QuadCFInterp object is passed in for efficiency */
  static void levelGradientMAC(
                               ///
                               LevelData<FluxBox>& a_edgeGrad,
                               ///
                               LevelData<FArrayBox>& a_phi,
                               ///
                               const LevelData<FArrayBox>* a_phiCrsePtr,
                               ///
                               const Real a_dx,
                               ///
                               const LayoutData<IntVectSet>& a_gridIVS,
                               ///
                               QuadCFInterp& a_cfInterpCrse);

  /// computes edge-centered level-operator gradient of cell-centered phi
  /**
      assumes _ALL_ ghost cell values have been preset (so phi can be const)
  */
  static void levelGradientMAC(
                               ///
                               LevelData<FluxBox>& a_edgeGrad,
                               ///
                               const LevelData<FArrayBox>& a_phi,
                               ///
                               const Real a_dx);

  /// computes edge-centered composite gradient of cell-centered phi
  /** if phiCrsePtr != NULL, does quadratic interpolation to compute
      coarse-fine boundary conditions for phi; since the edge between
      this level and finer levels is not considered to be a part of
      this level, there is no fine-level coarse-fine BC.  because of this,
      this function produces exactly the same results as LevelGradientMAC
      -- it's just included for completeness... */
  static void compGradientMAC(
                              ///
                              LevelData<FluxBox>& a_edgeGrad,
                              ///
                              LevelData<FArrayBox>& a_phi,
                              ///
                              const LevelData<FArrayBox>* a_phiCrse,
                              ///
                              const LevelData<FArrayBox>* a_phiFine,
                              ///
                              const Real a_dx,
                              ///
                              const int a_nRefCrse,
                              ///
                              const int a_nRefFine,
                              ///
                              const Box& a_dProblem);

  /// computes edge-centered composite gradient of cell-centered phi
  /**  if phiCrsePtr != NULL, does quadratic interpolation to compute
       coarse-fine boundary conditions for phi; since the edge between
       this level and finer levels is not considered to be a part of
       this level, there is no fine-level coarse-fine BC.  because of
       this, this function produces exactly the same results as
       LevelGradientMAC -- it's just included for completeness... */
  static void compGradientMAC(
                              ///
                              LevelData<FluxBox>& a_edgeGrad,
                              ///
                              LevelData<FArrayBox>& a_phi,
                              ///
                              const LevelData<FArrayBox>* a_phiCrse,
                              ///
                              const LevelData<FArrayBox>* a_phiFine,
                              ///
                              const Real a_dx,
                              ///
                              const int a_nRefCrse,
                              ///
                              const int a_nRefFine,
                              ///
                              const ProblemDomain& a_dProblem);

  /// computes edge-centered composite gradient of cell-centered phi
  /** if phiCrsePtr != NULL, does quadratic interpolation to compute
      coarse-fine boundary conditions for phi; since the edge between
      this level and finer levels is not considered to be a part of
      this level, there is no fine-level coarse-fine BC.  because of
      this, this function produces exactly the same results as
      LevelGradientMAC -- it's just included for completeness... In
      this one, a predefined QuadCFInterp object is passed in.  This
      (deprecated) interface uses a Box instead of a ProblemDomain object */
  static void compGradientMAC(
                              ///
                              LevelData<FluxBox>& a_edgeGrad,
                              ///
                              LevelData<FArrayBox>& a_phi,
                              ///
                              const LevelData<FArrayBox>* a_phiCrse,
                              ///
                              const LevelData<FArrayBox>* a_phiFine,
                              ///
                              const Real a_dx,
                              ///
                              const int a_nRefFine,
                              ///
                              const LayoutData<IntVectSet>& a_gridIVS,
                              ///
                              QuadCFInterp& a_cfInterpCrse);

  /// utility function for internal use
  static void singleBoxMacGrad(FArrayBox& a_gradFab,
                               const FArrayBox& a_phiFab,
                               int a_gradComp,
                               int a_phiComp,
                               int a_numComp,
                               const Box& a_edgeBox,
                               Real a_dx,
                               int a_dir,
                               int a_edgeDir,
                               const IntVectSet& a_gridIVS);

  /// create grid IVS for all grids
  /** creates an IVS representation of the DisjointBoxLayout
      in the neighborhood around each box in the DBL.
  */
  static void createGridIVS(LayoutData<IntVectSet>& a_gridsIVS,
                            const DisjointBoxLayout& a_grids,
                            const int a_nGrow);

  /// internal utility function
  /** creates IVS representation of the DisjointBoxLayout
      in the neighborhood of localBox
  */
  static void createGridIVS(IntVectSet& a_gridIVS,
                            const DisjointBoxLayout& a_grids,
                            const Box& a_localBox,
                            const int a_nGrow = 1);


};

#include "NamespaceFooter.H"
#endif
