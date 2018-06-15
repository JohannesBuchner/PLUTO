#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTUTIL_H_
#define _ADVECTUTIL_H_

#include "LevelData.H"
#include "FluxBox.H"
#include "FArrayBox.H"
#include "PatchGodunov.H"

#include "UsingNamespace.H"

/// utility function to do tracing part of advection velocity computation
void
traceAdvectionVel(LevelData<FluxBox>& a_advVel,
                  LevelData<FArrayBox>& a_old_vel,
                  const LevelData<FArrayBox>& a_viscousSource,
                  PatchGodunov& a_patchGodVelocity,
                  Real a_old_time, Real a_dt);

void
computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal,
                           LevelData<FArrayBox>& a_old_scal,
                           LevelData<FluxBox>& a_adv_vel,
                           LevelData<FArrayBox>& a_old_vel,
                           LevelData<FArrayBox>& a_diffusiveSrc,
                           PatchGodunov& a_patchGodScalar,
                           Real a_old_time, Real a_dt);

void
computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal,
                           LevelData<FluxBox>& a_edgeScalHalf,
                           LevelData<FArrayBox>& a_old_scal,
                           LevelData<FluxBox>& a_adv_vel,
                           LevelData<FArrayBox>& a_old_vel,
                           LevelData<FArrayBox>& a_diffusiveSrc,
                           PatchGodunov& a_patchGodScalar,
                           Real a_old_time, Real a_dt);

void
computePredictedVelocities(LevelData<FluxBox>& a_uHalf,
                           LevelData<FArrayBox>& a_traceVel,
                           LevelData<FluxBox>& a_advVel,
                           LevelData<FArrayBox>& a_viscousSource,
                           PatchGodunov& a_patchGodVelocity,
                           LevelData<FluxBox>& a_grad_eLambda,
                           LevelData<FluxBox>& a_gradPhi,
                           bool a_applyFreestreamCorrection,
                           Real a_old_time, Real a_dt);


#endif
