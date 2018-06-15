#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_COMPLEX_H_
#define _CH_COMPLEX_H_

#ifdef CH_USE_COMPLEX

#include <complex>
using std::complex;

#include "REAL.H"

#include "BaseNamespaceHeader.H"

/// Complex numbers suitable for use in FABs
/** Chombo has its own complex number class so it can:
 *  1) provide the same class name regardless of the current precision
 *     (just like the Chombo "Real" type)
 *  2) provide write access to the real and imaginary parts of the
 *     complex number without going through a constructor (this is not
 *     allowed by the C++ 1999 standard).
 */

class Complex : public complex<Real>
{
public:
  ///
  inline Complex()
  {
  }

  ///
  inline Complex(const Real& re, const Real& im):complex<Real>(re, im)
  {
  }

  /// return a non-const reference to the real part of the complex number
  inline Real & re();

  /// return a copy of the real part of the complex number
  inline Real   re() const
  {
    return real();
  };

  /// return a non-const reference to the imaginary part of the complex number
  inline Real & im();

  /// return a copy of the imaginary part of the complex number
  inline Real   im() const
  {
    return imag();
  };
};

//XXX dont check the compiler -- just assume this will work
//XXX#if defined(__GNUC__) || defined(__INTEL_COMPILER__) || defined(__IBMCPP__) || defined(__HP_aCC) || defined(__DECCXX)

// The compilers we've seen implement complex as a 2-array of
// numbers of the appropriate type, with the real part first.
// Run the test program "Chombo/lib/test/BoxTools/test_complex.cpp" to check.
// This is the most likely implementation, and the only real problems are
// whether the compiler puts the data at the beginning of the class's
// memory, and whether it will allow the cast.
inline Real & Complex::re()
{
  return (((Real*)this)[0]);
}

inline Real & Complex::im()
{
  return (((Real*)this)[1]);
}

//XXX#else
//XXX#error Complex type is not implemented for this compiler; contact the Chombo developers.
//XXX#endif

#include "BaseNamespaceFooter.H"
#endif //CH_USE_COMPLEX
#endif //_CH_COMPLEX_H_
