#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PoisselleTubeBCValue.H"
#include "PolyGeom.H"

#include "UsingNamespace.H"

/***/
RealVect
PoisselleTubeBCValue::
getVel(const Real& a_radius) const
{
  Real radsq = a_radius*a_radius;
  Real tubeRadSq = m_tubeRadius*m_tubeRadius;
  //    V = (1-r^2/r^2_0)*V_0
  Real vectSize = (1.0- (radsq/tubeRadSq))*m_maxVel;
  RealVect retval = m_tubeAxis*vectSize;
  return retval;
}
/***/
Real
PoisselleTubeBCValue::
getRadius(const RealVect& a_pt) const
{
  RealVect trueNorm, closestPt;
  PolyGeom::pointToLine(closestPt,  trueNorm, a_pt,
                        m_centerPt, m_tubeAxis);
  RealVect ebRadVec = closestPt;
  ebRadVec -= a_pt;

  Real ebRadius = PolyGeom::dot(ebRadVec, ebRadVec);
  ebRadius = sqrt(ebRadius);

  return ebRadius;
}
/***/
Real
PoisselleTubeBCValue::
value(const RealVect& a_point, const RealVect& a_normal, const Real& a_time, const int& a_comp) const
{
  Real radius = this->getRadius(a_point);
  RealVect velocity = getVel(radius);
  return velocity[m_velComp];
}

