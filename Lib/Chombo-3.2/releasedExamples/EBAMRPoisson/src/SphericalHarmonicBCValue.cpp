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
#include "SphericalHarmonicBCValue.H"

SphericalHarmonicBCValue::SphericalHarmonicBCValue()
{
  m_isDefined = false;

  m_lmphase = RealVect::Unit;
}

SphericalHarmonicBCValue::~SphericalHarmonicBCValue()
{
}

void SphericalHarmonicBCValue::define(const RealVect& a_lmphase)
{
  m_isDefined = true;

  m_lmphase = a_lmphase;
}
#include "EBAMRPoissonOp.H"
Real SphericalHarmonicBCValue::value(const RealVect& a_point,
                                     const RealVect& a_normal,
                                     const Real&     a_time,
                                     const int&      a_comp) const
{
 CH_assert(m_isDefined);

  Real value;

  FORT_GETSHPHIPOINT(CHF_REAL(value),
                     CHF_CONST_REALVECT(m_lmphase),
                     CHF_CONST_REALVECT(a_point),
                     CHF_CONST_REAL(a_time));

  return value;
}
