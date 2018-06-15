#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RealVect.H"

#include "functionsF_F.H"
#include "RZPolyValue.H"

RZPolyValue::RZPolyValue()
{
}

RZPolyValue::~RZPolyValue()
{
}

Real RZPolyValue::value(const RealVect& a_point,
                        const RealVect& a_normal,
                        const Real&     a_time,
                        const int&      a_comp) const
{
  Real value;

  FORT_GETPHIRZPOLY(CHF_REAL(value),
                    CHF_CONST_REALVECT(a_point));

  return value;
}
