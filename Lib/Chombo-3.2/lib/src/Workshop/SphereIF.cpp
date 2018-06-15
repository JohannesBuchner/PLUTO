#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphereIF.H"

#include "NamespaceHeader.H"

SphereIF::SphereIF(const Real&     a_radius,
                   const RealVect& a_center,
                   const bool&     a_inside)
{
  // Remember the parameters
  m_radius = a_radius;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}

SphereIF::SphereIF(const SphereIF& a_inputIF)
{
  // Remember the parameters
  m_radius  = a_inputIF.m_radius;
  m_center  = a_inputIF.m_center;
  m_inside  = a_inputIF.m_inside;

  m_radius2 = a_inputIF.m_radius2;
}

SphereIF::~SphereIF()
{
}

GeometryService::InOut SphereIF::InsideOutside(const RealVect& lo,
                                               const RealVect& hi) const
{

  // Fast Sphere-Box intersection from "Graphics Gems" pp 335-338, 1990

  Real dmin = 0; //distance from sphere center to nearest  point in box (squared)
  Real dmax = 0; //distance from sphere center to farthest point in box (squared)
  Real ai, bi, a, b;
  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      ai = m_center[i] - lo[i];
      bi = m_center[i] - hi[i];
      a = ai * ai;
      b = bi * bi;
      dmax = dmax + Max(a,b);
      if (m_center[i] < lo[i] || m_center[i] > hi[i]) dmin = dmin + Min(a,b);
    }

  if (m_inside)
    {
      if (dmin >= m_radius2) return GeometryService::Covered;
      if (dmax <  m_radius2) return GeometryService::Regular;
    }
  else
    {
      if (dmin >  m_radius2) return GeometryService::Regular;
      if (dmax <= m_radius2) return GeometryService::Covered;
    }
  return GeometryService::Irregular;
}


void SphereIF::GetParams(Real&     a_radius,
                         RealVect& a_center,
                         bool&     a_inside) const
{
  // Copy parameter information over
  a_radius = m_radius;
  a_center = m_center;
  a_inside = m_inside;
}

void SphereIF::SetParams(const Real&     a_radius,
                         const RealVect& a_center,
                         const bool&     a_inside)
{
  // Set parameter information
  m_radius = a_radius;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}

Real SphereIF::value(const RealVect& a_point) const
{
  Real retval;

  // The distance squared for m_center to a_point
  Real distance2;

  // Compute the distance squared
  distance2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    Real cur;
    cur = a_point[idir] - m_center[idir];

    distance2 += cur*cur;
  }

  // Return the difference between the sqaures (zero on the sphere)
  retval = distance2 - m_radius2;

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* SphereIF::newImplicitFunction() const
{
  SphereIF* spherePtr = new SphereIF(m_radius,
                                     m_center,
                                     m_inside);

  return static_cast<BaseIF*>(spherePtr);
}

#include "NamespaceFooter.H"
