#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_LANG_CC
// the entire reason for this file is to provide a consistent switch between
// old and new gradient implementations

// if NEW_GRADIENT is defined, we use the "new" gradient, which
// handles corner cells better for
// transverse gradients, at a cost of more complication, along with
// (most likely) worsened scaling in parallel
// this switch is implemented in the "GradientType.H" file in order to
// ensure that c++ and fortran files are consistent
#endif

#ifndef _GRADIENTTYPE_H_
#define _GRADIENTTYPE_H_

#define NEW_GRADIENT 1

#endif
