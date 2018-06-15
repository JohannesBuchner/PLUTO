#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _1DPOTENTIALSOLVE_H_
#define _1DPOTENTIALSOLVE_H_


//#define SAVE_SPACEDIM CH_SPACEDIM
#undef CH_SPACEDIM
#define CH_SPACEDIM  2
#include "LevelData.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "Vector.H"
//#undef CH_SPACEDIM
//#define CH_SPACEDIM SAVE_SPACEDIM
//#define CH_SPACEDIM 2
//#undef SAVE_SPACEDIM

using namespace Chombo;


int doOneDimensionalPotentialSolve(Vector<D2::LevelData<D2::FArrayBox>* >& a_rhs,
                                   Vector<D2::LevelData<D2::FArrayBox>* >& a_gradient,
                                   Vector<D2::ProblemDomain>& a_amrDomains,
                                   Vector<Real>& a_amrDx,
                                   Vector<int>& a_amrRefRatios,
                                   int a_finest_level,
                                   int a_lBase,
                                   int a_numComps,
                                   int a_verbosity);


#endif
