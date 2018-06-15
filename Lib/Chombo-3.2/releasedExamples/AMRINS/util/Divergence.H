#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIVERGENCE_H_
#define _DIVERGENCE_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelFluxRegister.H"
#include "QuadCFInterp.H"

#include "UsingNamespace.H"

/// Class to encapsulate Divergence functions

class Divergence
{

public:

  /** computes cell-centered level-operator divergence of
      cell-centered vector field u; if uCrse != NULL, does coarse-fine
      boundary conditions for u -- if  quadInterp == true, uses quadratic
      coarse-fine boundary conditions, otherwise, use extrap BC's */
  static void levelDivergenceCC(
                                ///
                                LevelData<FArrayBox>& a_div,
                                ///
                                LevelData<FArrayBox>& a_u,
                                ///
                                LevelData<FArrayBox>* a_uCrsePtr,
                                ///
                                const Real a_dx,
                                ///
                                const int a_nRefCrse,
                                ///
                                const ProblemDomain& a_dProblem,
                                ///
                                const bool a_quadInterp);

  /** computes cell-centered level-operator divergence of
      cell-centered vector field u; if uCrse != NULL, does coarse-fine
      boundary conditions for u -- if  quadInterp == true, uses quadratic
      coarse-fine boundary conditions, otherwise, use extrap BC's.
      This (deprecated) interface uses a Box instead of a ProblemDomain. */
  static void levelDivergenceCC(
                                ///
                                LevelData<FArrayBox>& a_div,
                                ///
                                LevelData<FArrayBox>& a_u,
                                ///
                                LevelData<FArrayBox>* a_uCrsePtr,
                                ///
                                const Real a_dx,
                                ///
                                const int a_nRefCrse,
                                ///
                                const Box& a_dProblem,
                                ///
                                const bool a_quadInterp);

  /** computes cell-centered level-operator divergence of
      cell-centered vector field u; if uCrse != NULL, does coarse-fine
      boundary conditions for u -- if  quadInterp == true, uses quadratic
      coarse-fine boundary conditions, otherwise, use extrap BC's.  This
      one takes a pre-allocated QuadCFInterp object as an argument */
  static void levelDivergenceCC(
                                ///
                                LevelData<FArrayBox>& a_div,
                                ///
                                LevelData<FArrayBox>& a_u,
                                ///
                                LevelData<FArrayBox>* a_uCrsePtr,
                                ///
                                const Real a_dx,
                                ///
                                const bool a_quadInterp,
                                ///
                                QuadCFInterp& a_cfInterp
                                );

  /** computes cell-centered composite-operator divergence
      of cell-centered vector field u; uses same coarse-fine BC's as
      levelDivergenceCC -- if uFine != NULL, also does flux-matching
      BC with finer level uFine, using same coarse-fine BC's as used
      for coarse level BC. */
  static void compDivergenceCC(
                               ///
                               LevelData<FArrayBox>& a_div,
                               ///
                               LevelData<FArrayBox>& a_u,
                               ///
                               LevelData<FArrayBox>* a_uCrsePtr,
                               ///
                               LevelData<FArrayBox>* a_uFinePtr,
                               ///
                               const Real a_dx,
                               ///
                               const int a_nRefCrse,
                               ///
                               const int a_nRefFine,
                               ///
                               const ProblemDomain& a_dProblem,
                               ///
                               const bool a_quadInterp);

  /** computes cell-centered composite-operator divergence
      of cell-centered vector field u; uses same coarse-fine BC's as
      levelDivergenceCC -- if uFine != NULL, also does flux-matching
      BC with finer level uFine, using same coarse-fine BC's as used
      for coarse level BC. This (deprecated) interface uses a Box instead
      of a ProblemDomain */
  static void compDivergenceCC(
                               ///
                               LevelData<FArrayBox>& a_div,
                               ///
                               LevelData<FArrayBox>& a_u,
                               ///
                               LevelData<FArrayBox>* a_uCrsePtr,
                               ///
                               LevelData<FArrayBox>* a_uFinePtr,
                               ///
                               const Real a_dx,
                               ///
                               const int a_nRefCrse,
                               ///
                               const int a_nRefFine,
                               ///
                               const Box& a_dProblem,
                               ///
                               const bool a_quadInterp);

  /** computes cell-centered composite-operator divergence
      of cell-centered vector field u; uses same coarse-fine BC's as
      levelDivergenceCC -- if uFine != NULL, also does flux-matching
      BC with finer level uFine, using same coarse-fine BC's as used
      for coarse level BC. In this one, pre-allocated LevelFluxRegister
      and QuadCFInterp objects are passed in as an argument (saves
      allocation inside fn.)*/
  static void compDivergenceCC(
                               ///
                               LevelData<FArrayBox>& a_div,
                               ///
                               LevelData<FArrayBox>& a_u,
                               ///
                               LevelData<FArrayBox>* a_uCrsePtr,
                               ///
                               LevelData<FArrayBox>* a_uFinePtr,
                               ///
                               const Real a_dx,
                               ///
                               const int a_nRefFine,
                               ///
                               const int a_nRefCrse,
                               ///
                               const ProblemDomain& a_dProblem,
                               ///
                               const bool a_quadInterp,
                               ///
                               LevelFluxRegister* a_fluxRegFinePtr,
                               ///
                               QuadCFInterp& a_cfInterpCrse,
                               ///
                               QuadCFInterp& a_cfInterpFine);

  /** computes cell-centered composite-operator divergence
      of cell-centered vector field u; uses same coarse-fine BC's as
      levelDivergenceCC -- if uFine != NULL, also does flux-matching
      BC with finer level uFine, using same coarse-fine BC's as used
      for coarse level BC. In this one, pre-allocated LevelFluxRegister
      and QuadCFInterp objects are passed in as an argument (saves
      allocation inside fn.)  This (deprecated) interface uses a
      Box instead of a ProblemDomain. */
  static void compDivergenceCC(
                               ///
                               LevelData<FArrayBox>& a_div,
                               ///
                               LevelData<FArrayBox>& a_u,
                               ///
                               LevelData<FArrayBox>* a_uCrsePtr,
                               ///
                               LevelData<FArrayBox>* a_uFinePtr,
                               ///
                               const Real a_dx,
                               ///
                               const int a_nRefFine,
                               ///
                               const int a_nRefCrse,
                               ///
                               const Box& a_dProblem,
                               ///
                               const bool a_quadInterp,
                               ///
                               LevelFluxRegister* a_fluxRegFinePtr,
                               ///
                               QuadCFInterp& a_cfInterpCrse,
                               ///
                               QuadCFInterp& a_cfInterpFine);

  /** computes cell-centered level-operator divergence of
      edge-centered vector field uEdge; assumes all coarse-fine BC's
      have already been set. */
  static void levelDivergenceMAC(
                                 ///
                                 LevelData<FArrayBox>& a_div,
                                 ///
                                 const LevelData<FluxBox>& a_uEdge,
                                 ///
                                 const Real a_dx);

  /// really basic single-fab MAC divergence
  static void simpleDivergenceMAC(///
                                  FArrayBox& a_div,
                                  ///
                                  const FluxBox& a_uEdge,
                                  ///
                                  const Real a_dx);

  /** computes composite cell-centered divergence of edge-centered
      vector field uEdge; if finer level data uEdgeFine exists, use
      flux-matching condition to compute divergence along coarse-fine
      interface. */
  static void compDivergenceMAC(
                                ///
                                LevelData<FArrayBox>& a_div,
                                ///
                                LevelData<FluxBox>& a_uEdge,
                                ///
                                LevelData<FluxBox>* a_uEdgeFinePtr,
                                ///
                                const Real a_dx,
                                ///
                                const Real* a_dxFine,
                                ///
                                const int a_nRefFine,
                                ///
                                const ProblemDomain& a_dProblem);

  /** computes composite cell-centered divergence of edge-centered
      vector field uEdge; if finer level data uEdgeFine exists, use
      flux-matching condition to compute divergence along coarse-fine
      interface. This (deprecated) interface uses a Box instead of a
      ProblemDomain */
  static void compDivergenceMAC(
                                ///
                                LevelData<FArrayBox>& a_div,
                                ///
                                LevelData<FluxBox>& a_uEdge,
                                ///
                                LevelData<FluxBox>* a_uEdgeFinePtr,
                                ///
                                const Real a_dx,
                                ///
                                const Real* a_dxFine,
                                ///
                                const int a_nRefFine,
                                ///
                                const Box& a_dProblem);

  /** just like normal compDivergenceMAC, but pass in a predefined
      flux register (more efficient) */
  static void compDivergenceMAC(
                                ///
                                LevelData<FArrayBox>& a_div,
                                ///
                                LevelData<FluxBox>& a_uEdge,
                                ///
                                LevelData<FluxBox>* a_uEdgeFinePtr,
                                ///
                                LevelFluxRegister* a_fluxRegPtr,
                                ///
                                const Real a_dx,
                                ///
                                const Real* a_dxFine,
                                ///
                                const int a_nRefine,
                                ///
                                const ProblemDomain& a_dProblem);

  /** just like normal compDivergenceMAC, but pass in a predefined
      flux register (more efficient). This (deprecated) interface uses a
      Box instead of a ProblemDomain. */
  static void compDivergenceMAC(
                                ///
                                LevelData<FArrayBox>& a_div,
                                ///
                                LevelData<FluxBox>& a_uEdge,
                                ///
                                LevelData<FluxBox>* a_uEdgeFinePtr,
                                ///
                                LevelFluxRegister* a_fluxRegPtr,
                                ///
                                const Real a_dx,
                                ///
                                const Real* a_dxFine,
                                ///
                                const int a_nRefine,
                                ///
                                const Box& a_dProblem);
};

#endif
