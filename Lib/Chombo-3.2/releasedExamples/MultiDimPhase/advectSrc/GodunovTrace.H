#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GODUNOVTRACE_H_
#define _GODUNOVTRACE_H_

#include "FluxBox.H"
#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "UsingNamespace.H"

// Set of Functions which perform upwinded tracing


/*@ManMemo: trace state to edges in direction dir at time t + dt/2
  -- this version uses edge-centered advection velocities which
  are passed in */
void TraceState(/// state at time t+dt/2 on edges in direction dir
                FArrayBox& a_stateHalf,
                /// cell-centered state at time t
                const FArrayBox& a_state,
                /// cell-centered velocity at time t
                const FArrayBox& a_cellVel,
                /// edge-centered advection velocity at time t+dt/2
                const FluxBox& a_advectionVel,
                /// cell-centered source
                const FArrayBox& a_source,
                /// Physical domain
                const ProblemDomain& a_dProblem,
                /// interior of grid patch
                const Box& a_gridBox,
                /// timeStep
                const Real a_dt,
                /// cell-spacing
                const Real a_dx,
                /// direction in which to perform tracing
                const int dir,
                /// which components to trace
                const Interval& a_srcComps,
                ///  where to put traced components in a_stateHalf,
                const Interval& a_destComps);


/*@ManMemo: trace state to edges in direction dir at time t + dt/2
  -- this version uses edge-centered advection velocities which
  are passed in; deprecated interface using a Box for the domain
  instead of a ProblemDomain (will default to non-periodic case) */
void TraceState(/// state at time t+dt/2 on edges in direction dir
                FArrayBox& a_stateHalf,
                /// cell-centered state at time t
                const FArrayBox& a_state,
                /// cell-centered velocity at time t
                const FArrayBox& a_cellVel,
                /// edge-centered advection velocity at time t+dt/2
                const FluxBox& a_advectionVel,
                /// cell-centered source
                const FArrayBox& a_source,
                /// Physical domain
                const Box& a_dProblem,
                /// interior of grid patch
                const Box& a_gridBox,
                /// timeStep
                const Real a_dt,
                /// cell-spacing
                const Real a_dx,
                /// direction in which to perform tracing
                const int dir,
                /// which components to trace
                const Interval& a_srcComps,
                ///  where to put traced components in a_stateHalf,
                const Interval& a_destComps);


/*@ManMemo: trace state to edges in direction dir at time t+dt/2
  -- this version does not have edge-centered advection velocities
  passed in; uses cell-centered velocity field to approximate them */
void TraceState(/// state at time t+dt/2 on edges in direction dir
                FArrayBox& a_stateHalf,
                /// cell-centered state at time t
                const FArrayBox& a_state,
                /// cell-centered velocity at time t
                const FArrayBox& a_cellVel,
                /// cell-centered source
                const FArrayBox& a_source,
                /// Physical domain
                const ProblemDomain& a_dProblem,
                /// interior of grid patch
                const Box& a_gridBox,
                /// timeStep
                const Real a_dt,
                /// cell-spacing
                const Real a_dx,
                /// direction in which to perform tracing
                const int dir,
                /// which components to trace
                const Interval& a_srcComps,
                ///  where to put traced components in a_stateHalf,
                const Interval& a_destComps);


/*@ManMemo: trace state to edges in direction dir at time t+dt/2
  -- this version does not have edge-centered advection velocities
  passed in; uses cell-centered velocity field to approximate them.
  deprecated interface using a Box for the domain
  instead of a ProblemDomain (will default to non-periodic case) */

void TraceState(/// state at time t+dt/2 on edges in direction dir
                FArrayBox& a_stateHalf,
                /// cell-centered state at time t
                const FArrayBox& a_state,
                /// cell-centered velocity at time t
                const FArrayBox& a_cellVel,
                /// cell-centered source
                const FArrayBox& a_source,
                /// Physical domain
                const Box& a_dProblem,
                /// interior of grid patch
                const Box& a_gridBox,
                /// timeStep
                const Real a_dt,
                /// cell-spacing
                const Real a_dx,
                /// direction in which to perform tracing
                const int dir,
                /// which components to trace
                const Interval& a_srcComps,
                ///  where to put traced components in a_stateHalf,
                const Interval& a_destComps);







#endif


