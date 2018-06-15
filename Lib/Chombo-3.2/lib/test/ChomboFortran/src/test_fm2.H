#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TEST_FM2_H_
#define _TEST_FM2_H_

#include "FORT_PROTO.H"
#include "REAL.H"

extern "C"
{

// Prototype for Fortran procedure test_FM1 ...
//
#define FORT_TEST_FM2 FORTRAN_NAME( TEST_FM2 ,test_fm2 )
void
FORT_TEST_FM2(
      CHFp_INT(a_i) );

}

#endif
