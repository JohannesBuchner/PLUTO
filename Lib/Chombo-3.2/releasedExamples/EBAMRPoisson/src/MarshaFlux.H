#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MARSHAFLUX_H_
#define _MARSHAFLUX_H_

#include "RealVect.H"

#include "BaseBCValue.H"
#include "UsingNamespace.H"

class MarshaFlux: public BaseBCValue
{
public:
  MarshaFlux();

  virtual ~MarshaFlux();

  virtual Real value(const RealVect& a_point,
                     const RealVect& a_normal,
                     const Real&     a_time,
                     const int&      a_comp) const;

private:
};

#endif
