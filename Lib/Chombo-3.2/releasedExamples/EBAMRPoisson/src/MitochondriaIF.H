#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MITOCHONDRIAIF_H_
#define _MITOCHONDRIAIF_H_

#include "BaseIF.H"
#include "TransformIF.H"

#include "UsingNamespace.H"

class SphereDimpleIF: public BaseIF
{
public:
  SphereDimpleIF();

  virtual ~SphereDimpleIF();

  virtual Real value(const RealVect & a_point) const;

  virtual BaseIF* newImplicitFunction() const;
};

class CylinderDimpleIF: public BaseIF
{
public:
  CylinderDimpleIF();

  virtual ~CylinderDimpleIF();

  virtual Real value(const RealVect & a_point) const;

  virtual BaseIF* newImplicitFunction() const;
};

class MitochondriaIF: public BaseIF
{
public:
  MitochondriaIF();

  virtual ~MitochondriaIF();

  virtual Real value(const RealVect & a_point) const;

  virtual BaseIF* newImplicitFunction() const;

private:
  BaseIF* m_dimple;
#if 0
  TransformIF* m_dimple2;
  TransformIF* m_dimple3;
#endif
};

#endif
