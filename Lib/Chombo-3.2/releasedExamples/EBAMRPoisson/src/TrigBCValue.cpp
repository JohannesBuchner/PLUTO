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
#include "TrigBCValue.H"

TrigBCValue::TrigBCValue()
{
  m_isDefined = false;

  m_trig = RealVect::Unit;
}

TrigBCValue::~TrigBCValue()
{
}

void TrigBCValue::define(const RealVect& a_trig)
{
  m_isDefined = true;

  m_trig = a_trig;
}
#include "EBAMRPoissonOp.H"
Real TrigBCValue::value(const RealVect& a_point,
                        const RealVect& a_normal,
                        const Real&     a_time,
                        const int&      a_comp) const
{
 CH_assert(m_isDefined);

  Real value;

  FORT_GETPHIPOINT(CHF_REAL(value),
                   CHF_CONST_INTVECT(m_trig),
                   CHF_CONST_REALVECT(a_point),
                   CHF_CONST_REAL(a_time));

  return value;
}
