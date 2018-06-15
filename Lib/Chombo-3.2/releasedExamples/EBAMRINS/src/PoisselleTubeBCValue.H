#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#ifndef _POISSELLETUBEBCVALUE_H_
#define _POISSELLETUBEBCVALUE_H_

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "EBIBC.H"
#include "EBIBCFactory.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "BaseEBBC.H"
#include "PolyGeom.H"

#include "UsingNamespace.H"

///instantiation of basebcvalue that always returns a parabola around center of the domain
/**
    vel  = inflowvel at r=0
    dv/dr = 0 at r=0
    vel = 0 at r=a_poisslleRadius
    V = (1-r^2/r^2_0)*V_0
*/
class PoisselleTubeBCValue: public BaseBCValue
{
public:

  PoisselleTubeBCValue(const RealVect& a_centerPt,
                       const RealVect& a_tubeAxis,
                       const Real&     a_tubeRadius,
                       const Real&     a_maxVel,
                       const int&      a_velComp)
  {
    m_centerPt   = a_centerPt;
    m_tubeAxis   = a_tubeAxis;
    m_tubeRadius = a_tubeRadius;
    m_maxVel     = a_maxVel;
    m_velComp    = a_velComp;
    Real sum;
    PolyGeom::unifyVector(m_tubeAxis, sum);
  }

  ~PoisselleTubeBCValue()
  {
  }

  Real value(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const;

  Real     getRadius(const RealVect& a_point) const;
  RealVect    getVel(const Real&     a_radius) const;

  PoisselleTubeBCValue()
  {
    m_tubeRadius = -1;
  }

private:
  RealVect      m_centerPt;
  RealVect      m_tubeAxis;
  Real          m_tubeRadius;
  Real          m_maxVel;
  int           m_velComp;
};

#endif
