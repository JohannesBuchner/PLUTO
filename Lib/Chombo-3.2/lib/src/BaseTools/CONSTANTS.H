#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_
#include "BaseNamespaceHeader.H"

#if defined(CH_USE_FLOAT) && !defined(CH_LANG_CC)
#define zero (0.0e0)
#endif
#if defined(CH_USE_FLOAT) || defined(CH_LANG_CC)

#define one      (1.0e0)
#define two      (2.0e0)
#define three    (3.0e0)
#define four     (4.0e0)
#define five     (5.0e0)
#define six      (6.0e0)
#define seven    (7.0e0)
#define eight    (8.0e0)
#define nine     (9.0e0)
#define ten     (10.0e0)
#define twelve  (12.0e0)
#define fifteen (15.0e0)
#define sixteen (16.0e0)
#define twenty  (20.0e0)
#define thirty  (30.0e0)
#define tenth   (0.100e0)
#define eighth  (0.125e0)
#define sixth   (1.000e0 /  6.000e0)
#define twelfth (1.000e0 / 12.000e0)
#define thirtieth (1.000e0 / 30.000e0)
#define fourth  (0.250e0)
#define third   (1.000e0 /  3.000e0)
#define half    (0.500e0)
#define two3rd  (2.000e0 /  3.000e0)
#define Pi      (3.14159265358979323846264338327950288e0)

#define ZERO     (0.0E0)
#define ONE      (1.0E0)
#define TWO      (2.0E0)
#define THREE    (3.0E0)
#define FOUR     (4.0E0)
#define FIVE     (5.0E0)
#define SIX      (6.0E0)
#define SEVEN    (7.0E0)
#define EIGHT    (8.0E0)
#define NINE     (9.0E0)
#define TEN     (10.0E0)
#define TWELVE  (12.0E0)
#define FIFTEEN (15.0E0)
#define SIXTEEN (16.0E0)
#define TWENTY  (20.0E0)
#define THIRTY  (30.0E0)
#define TENTH   (0.100E0)
#define EIGHTH  (0.125E0)
#define SIXTH   (1.000E0 / 6.000E0)
#define TWELFTH (1.000E0 / 12.000E0)
#define THIRTIETH (1.000E0 / 30.000E0)
#define FOURTH  (0.250E0)
#define THIRD   (1.000E0 / 3.000E0)
#define HALF    (0.500E0)
#define TWO3RD  (2.000E0 / 3.000E0)
#define PI      (3.14159265358979323846264338327950288E0)

#elif defined(CH_USE_DOUBLE) && !defined(CH_LANG_CC)

#define zero     (0.0d0)
#define one      (1.0d0)
#define two      (2.0d0)
#define three    (3.0d0)
#define four     (4.0d0)
#define five     (5.0d0)
#define six      (6.0d0)
#define seven    (7.0d0)
#define eight    (8.0d0)
#define nine     (9.0d0)
#define ten     (10.0d0)
#define twelve  (12.0d0)
#define fifteen (15.0d0)
#define sixteen (16.0d0)
#define twenty  (20.0d0)
#define thirty  (30.0d0)
#define tenth   (0.100d0)
#define eighth  (0.125d0)
#define sixth   (1.000d0 / 6.000d0)
#define twelfth (1.000d0 / 12.000d0)
#define thirtieth (1.000d0 / 30.000d0)
#define fourth  (0.250d0)
#define third   (1.000d0 / 3.000d0)
#define half    (0.500d0)
#define two3rd  (2.000d0 / 3.000d0)
#define Pi      (3.14159265358979323846264338327950288d0)

#define ZERO    (0.0D0)
#define ONE     (1.0D0)
#define TWO     (2.0D0)
#define THREE   (3.0D0)
#define FOUR    (4.0D0)
#define FIVE    (5.0D0)
#define SIX     (6.0D0)
#define SEVEN   (7.0D0)
#define EIGHT   (8.0D0)
#define NINE    (9.0D0)
#define TEN     (10.0D0)
#define TWELVE  (12.0D0)
#define FIFTEEN (15.0D0)
#define SIXTEEN (16.0D0)
#define TWENTY  (20.0D0)
#define THIRTY  (30.0D0)
#define TENTH   (0.100D0)
#define EIGHTH  (0.125D0)
#define SIXTH   (1.000D0 / 6.000D0)
#define TWELFTH (1.000D0 / 12.000D0)
#define THIRTIETH (1.000D0 / 30.000D0)
#define FOURTH  (0.250D0)
#define THIRD   (1.000D0 / 3.000D0)
#define HALF    (0.500D0)
#define TWO3RD  (2.000D0 / 3.000D0)
#define PI      (3.14159265358979323846264338327950288D0)

#else

#error One of CH_USE_(FLOAT,DOUBLE) must be defined

#endif

#include "BaseNamespaceFooter.H"
#endif /*CH_CONSTANTS_H*/
