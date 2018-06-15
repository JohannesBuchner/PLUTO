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
#include "SphericalHarmonicBCBetaFlux.H"

SphericalHarmonicBCBetaFlux::SphericalHarmonicBCBetaFlux()
{
  m_isDefined = false;

  m_lmphase = RealVect::Unit;
}

SphericalHarmonicBCBetaFlux::~SphericalHarmonicBCBetaFlux()
{
}

void SphericalHarmonicBCBetaFlux::define(const RealVect& a_lmphase)
{
  m_isDefined = true;

  m_lmphase = a_lmphase;
}
#include "EBAMRPoissonOp.H"
Real SphericalHarmonicBCBetaFlux::value(const RealVect& a_point,
                                        const RealVect& a_normal,
                                        const Real&     a_time,
                                        const int&      a_comp) const
{
 CH_assert(m_isDefined);

  RealVect gradient;
  FORT_GETBETAGRADSHPHIPOINT(CHF_REALVECT(gradient),
                             CHF_CONST_REALVECT(m_lmphase),
                             CHF_CONST_REALVECT(a_point),
                             CHF_CONST_REAL(a_time));

  Real flux = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    flux += gradient[idir] * a_normal[idir];
  }

  return flux;
}
