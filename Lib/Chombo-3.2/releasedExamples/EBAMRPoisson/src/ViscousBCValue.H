#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VISCOUSBCVALUE_H_
#define _VISCOUSBCVALUE_H_

#include "RealVect.H"

#include "BaseBCFuncEval.H"
#include "MayDay.H"
#include "UsingNamespace.H"

class ViscousTrigBCValue: public BaseBCFuncEval
{
public:
  ViscousTrigBCValue();

  virtual ~ViscousTrigBCValue();

  virtual Real value(const RealVect&       a_point,
                     const int&            a_comp) const;

  virtual Real derivative(const RealVect&       a_point,
                          const int&            a_comp,
                          const int&            a_derivDir
                          ) const;

  virtual void define(const RealVect& a_trig);

private:
  bool m_isDefined;

  RealVect m_trig;
};

#endif
