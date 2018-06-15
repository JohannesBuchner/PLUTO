#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SETVALLEVEL_H_
#define _SETVALLEVEL_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "UsingNamespace.H"

///
/**
   Set a_levels[a_min:a_max] to zero.
 */
void setValLevels(Vector< LevelData<FArrayBox>* >&   a_levels,
                  int                                a_min,
                  int                                a_max,
                  Real                               a_val);

///
/**
   Set a_level to zero.
 */
void setValLevel(LevelData<FArrayBox>&   a_level,
                 Real                    a_val);


///
/**
   Set a_level to zero.
 */
void setValLevel(LevelData<FluxBox>&   a_level,
                 Real                  a_val);

#endif
