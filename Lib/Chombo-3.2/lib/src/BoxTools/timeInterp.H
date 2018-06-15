#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TIMEINTERP_H_
#define _TIMEINTERP_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "Interval.H"
#include "NamespaceHeader.H"

/// interpolate data in time -- assumes one-to-one correspondence in comps
void
timeInterp(LevelData<FArrayBox>& a_dest, Real a_time,
           const LevelData<FArrayBox>& a_old_phi, Real a_old_time,
           const LevelData<FArrayBox>& a_new_phi, Real a_new_time,
           const Interval& a_comps);

/// interpolate data in time -- dest comps can be misaligned with src comps
void
timeInterp(LevelData<FArrayBox>& a_dest, Real a_time,
           const LevelData<FArrayBox>& a_old_phi, Real a_old_time,
           const LevelData<FArrayBox>& a_new_phi, Real a_new_time,
           const Interval& a_src_comps, const Interval& a_dest_comps);

#include "NamespaceFooter.H"

#endif
