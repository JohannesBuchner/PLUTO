#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPUTENORM_H_
#define _COMPUTENORM_H_

#include "REAL.H"
#include "Vector.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Interval.H"
#include "FluxBox.H"
#include "RealVect.H"
#include "NamespaceHeader.H"

/// Returns the volume-weighted p-norm of phi over all valid regions
Real computeNorm(const Vector<LevelData<FArrayBox>* >& a_phi,
                 const Vector<int>&                    a_nRefFine,
                 const Real                            a_dxCrse,
                 const Interval                        a_comps = Interval(0,0),
                 const int                             a_p = 2,
                 const int                             a_lBase = 0);

/// Returns the volume-weighted p-norm of phi over the valid region
Real computeNorm(const LevelData<FArrayBox>& a_phi,
                 const DisjointBoxLayout*    a_finerGrids,
                 const int                   a_nRefFine,
                 const Real                  a_dx,
                 const Interval              a_comps = Interval(0,0),
                 const int                   a_p = 2);

/// Returns the volume-weighted p-norm of phi over all valid regions: anisotropic version
Real computeNorm(const Vector<LevelData<FArrayBox>* >& a_phi,
                 const Vector<IntVect>&                a_nRefFine,
                 const RealVect&                       a_dxCrse,
                 const Interval                        a_comps = Interval(0,0),
                 const int                             a_p = 2,
                 const int                             a_lBase = 0);

/// Returns the volume-weighted p-norm of phi over the valid region: anisotropic version
Real computeNorm(const LevelData<FArrayBox>& a_phi,
                 const DisjointBoxLayout*    a_finerGrids,
                 const IntVect&              a_nRefFine,
                 const RealVect&             a_dx,
                 const Interval              a_comps = Interval(0,0),
                 const int                   a_p = 2);

/// Returns the max value of phi over all valid regions
Real computeMax(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Interval                        a_comps = Interval(0,0),
                const int                             a_lBase = 0);

/// Returns the maximum value of phi over the valid region
Real computeMax(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGrids,
                const int                   a_nRefFine,
                const Interval              a_comps = Interval(0,0));

/// Returns the min value of phil over all valid regions
Real computeMin(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<int>&                    a_nRefFine,
                const Interval                        a_comps = Interval(0,0),
                const int                             a_lBase = 0);

/// Returns the minimum value of phi over the valid region
Real computeMin(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGrids,
                const int                   a_nRefFine,
                const Interval              a_comps = Interval(0,0));

/// Extension of the norm(...) function for FluxBoxes
Real norm(const BoxLayoutData<FluxBox>& a_layout,
          const Interval&               a_interval,
          const int&                    a_p);
#include "NamespaceFooter.H"
#endif
