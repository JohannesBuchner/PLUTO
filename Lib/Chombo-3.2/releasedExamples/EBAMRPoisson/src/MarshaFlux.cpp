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
#include "MarshaFlux.H"

MarshaFlux::MarshaFlux()
{
}

MarshaFlux::~MarshaFlux()
{
}

Real MarshaFlux::value(const RealVect& a_point,
                       const RealVect& a_normal,
                       const Real&     a_time,
                       const int&      a_comp)  const
{
  RealVect gradient;
  FORT_GETMARSHAGRADPHIPOINT(CHF_REALVECT(gradient),
                             CHF_CONST_REALVECT(a_point));

  Real flux = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    flux += gradient[idir] * a_normal[idir];
  }

  return flux;
}
