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
#include "TrigBCBetaFlux.H"

TrigBCBetaFlux::TrigBCBetaFlux()
{
  m_isDefined = false;

  m_trig = RealVect::Unit;
}

TrigBCBetaFlux::~TrigBCBetaFlux()
{
}

void TrigBCBetaFlux::define(const RealVect& a_trig)
{
  m_isDefined = true;

  m_trig = a_trig;
}
#include "EBAMRPoissonOp.H"
Real TrigBCBetaFlux::value(const RealVect& a_point,
                           const RealVect& a_normal,
                           const Real&     a_time,
                           const int&      a_comp) const
{
 CH_assert(m_isDefined);

  RealVect gradient;
  FORT_GETBETAGRADPHIPOINT(CHF_REALVECT(gradient),
                           CHF_CONST_REALVECT(m_trig),
                           CHF_CONST_REALVECT(a_point),
                           CHF_CONST_REAL(a_time));

  Real flux = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    flux += gradient[idir] * a_normal[idir];
  }

  return flux;
}
