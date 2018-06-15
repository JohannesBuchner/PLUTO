#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SphereArrayIF.H"

#include "NamespaceHeader.H"

SphereArrayIF::SphereArrayIF(const Real&     a_radius,
                             const RealVect& a_firstCenter,
                             const RealVect& a_spacing)
{
  // Remember the parameters
  m_radius = a_radius;
  m_firstCenter = a_firstCenter;
  m_spacing = a_spacing;

  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}

SphereArrayIF::SphereArrayIF(const SphereArrayIF& a_inputIF)
{
  // Remember the parameters
  m_radius      = a_inputIF.m_radius;
  CH_assert(m_radius > 0);
  m_firstCenter = a_inputIF.m_firstCenter;
  m_spacing     = a_inputIF.m_spacing;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_spacing[idir] < 2.5*m_radius)
        {
          MayDay::Abort("SphereArrayIF: spheres too close");
        }
    }
  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}

SphereArrayIF::~SphereArrayIF()
{
}

Real SphereArrayIF::value(const RealVect& a_point) const
{
  //compute which center we are closest to
  RealVect dist = a_point;
  dist -= m_firstCenter;
  dist /= m_spacing;
  //make distance integer numbers in all directions
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      dist[idir] = rint(dist[idir]);
    }
  //now multiply the distance by the spacing and add it to the firstCenter;
  dist *= m_spacing;
  RealVect center = m_firstCenter;
  center += dist;

  // The distance squared for m_center to a_point
  Real distance2;

  // Compute the distance squared
  distance2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    Real cur;
    cur = a_point[idir] - center[idir];

    distance2 += cur*cur;
  }

  // Return the difference between the sqaures (zero on the sphere)
  Real retval = m_radius2 - distance2;

  return retval;
}
                                                       //
BaseIF* SphereArrayIF::newImplicitFunction() const
{
  SphereArrayIF* spherePtr = new SphereArrayIF(m_radius,
                                               m_firstCenter,
                                               m_spacing);

  return static_cast<BaseIF*>(spherePtr);
}

GeometryService::InOut SphereArrayIF::InsideOutside(const RealVect& a_low,
                                                    const RealVect& a_high) const
{
  RealVect lo(a_low);
  RealVect hi(a_high);


  lo-=m_firstCenter;
  hi-=m_firstCenter;

  // first, find nearst spheres to hi and lo
  IntVect nearestLo, nearestHi;
  bool differByTwo = true;
  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      nearestLo[i] = (int)(lo[i]/m_spacing[i]);
      nearestHi[i] = (int)(hi[i]/m_spacing[i]);
      if ( nearestHi[i] - nearestLo[i] <= 1 ) differByTwo = false;
    }
  if (differByTwo)
    {
      return GeometryService::Irregular;
    }

  nearestLo-=IntVect::Unit;
  nearestHi+=IntVect::Unit;
  Box near(nearestLo, nearestHi);
  BoxIterator n(near);

  for (n.reset(); n.ok(); ++n)
  {
    Real dmin = 0;
    Real dmax = 0;
    Real ai, bi, a, b, c;
    for (int i=0; i<CH_SPACEDIM; ++i)
      {
        c = n()[i]*m_spacing[i];
        ai= c-lo[i];
        bi= c-hi[i];
        a = ai*ai;
        b = bi*bi;
        dmax = dmax + Max(a,b);
        if (c<lo[i] || c> hi[i]) dmin = dmin + Min(a,b);
      }
    if (dmax < m_radius2)
      {
        return GeometryService::Covered;
      }
    if (dmin <= m_radius2 && dmax >= m_radius2) return GeometryService::Irregular;
  }

  return GeometryService::Regular;

}

#include "NamespaceFooter.H"
