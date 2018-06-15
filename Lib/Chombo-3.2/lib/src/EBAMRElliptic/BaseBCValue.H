#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEBCVALUE_H_
#define _BASEBCVALUE_H_

#include "RealVect.H"
#include "EBArith.H"
#include "NamespaceHeader.H"

class BaseBCValue
{
public:
  BaseBCValue()
  {
  }

  virtual ~BaseBCValue()
  {
  }

  virtual Real value(const RealVect& a_point,
                     const RealVect& a_normal,
                     const Real&     a_time,
                     const int&      a_comp) const = 0;

  virtual Real value(const FaceIndex&      a_face,
                     const Side::LoHiSide& a_side,
                     const DataIndex&      a_dit,
                     const RealVect&       a_point,
                     const RealVect&       a_normal,
                     const Real&           a_time,
                     const int&            a_comp) const
  {
    return value(a_point,a_normal,a_time,a_comp);
  }

  virtual Real value(const VolIndex&  a_vof,
                     const RealVect&  a_centroid,
                     const RealVect&  a_normal,
                     const RealVect&  a_dx,
                     const RealVect&  a_probLo,
                     const DataIndex& a_dit,
                     const Real&      a_time,
                     const int&       a_comp) const
  {
    RealVect point = a_centroid;
    point *= a_dx;
    point += EBArith::getVofLocation(a_vof,a_dx,a_probLo);
    return value(point,a_normal,a_time,a_comp);
  }

  virtual Real value(const IntVect&   a_iv,
                     const DataIndex& a_dit,
                     const RealVect&  a_point,
                     const RealVect&  a_normal,
                     const Real&      a_time,
                     const int&       a_comp) const
  {
    return value(a_point,a_normal,a_time,a_comp);
  }

  virtual Real beta(const RealVect& a_point,
                    const Real&     a_time) const
  {
    MayDay::Error("BaseBCValue::beta virtual function not implemented");
    return 1.2345678e20;
  }

};

#include "NamespaceFooter.H"
#endif
