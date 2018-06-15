#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REAL_H_
#define _REAL_H_

#ifdef CH_USE_FLOAT /* precision */
#  undef  CH_USE_DOUBLE
#  undef  CH_USE_FLOAT
#  define CH_USE_FLOAT  1
#else   /* precision */
#  undef  CH_USE_FLOAT
#  undef  CH_USE_DOUBLE
#  define CH_USE_DOUBLE 1
#endif  /* precision */

#if !defined(CH_LANG_FORT) /* language */
#  ifdef CH_USE_FLOAT /* precision */
     typedef float Real;
#    define MPI_CH_REAL         MPI_FLOAT
#    define H5T_NATIVE_REAL     H5T_NATIVE_FLOAT
#    define BASEFAB_REAL_SETVAL 1.23456789e+30
#    define CH_BADVAL           (Real)1.0e30
#  else   /* precision */
     typedef double Real;
#    define MPI_CH_REAL         MPI_DOUBLE
#    define H5T_NATIVE_REAL     H5T_NATIVE_DOUBLE
#    define BASEFAB_REAL_SETVAL 1.23456789e+300
#    define CH_BADVAL           (Real)1.0e99
#  endif  /* precision */

#else /* language */

#  ifdef CH_USE_FLOAT /* precision */
#    define REAL_T REAL
#    define COMPLEX_T COMPLEX
#  else /* precision */
#    if defined(CH_CRAY)
#      define REAL_T REAL
#      define COMPLEX_T COMPLEX
#    else /*cray */
#      define REAL_T REAL*8
#      define COMPLEX_T COMPLEX*16
#    endif /*cray*/
#  endif /* precision */

#endif /* language */

#if defined(CH_LANG_FORT)
#  define Real_T REAL_T
#  define Real_t REAL_T
#  define real_t REAL_T
#  define Complex_T COMPLEX_T
#  define complex_T COMPLEX_T
#  define complex_t COMPLEX_T
#endif

#if !defined(CH_USE_FLOAT) && !defined(CH_USE_DOUBLE)
#error Exactly one of CH_USE_FLOAT or CH_USE_DOUBLE must be defined
#endif

#endif
