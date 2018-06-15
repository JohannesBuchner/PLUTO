#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GaussianSphereIF.H"

#include "NamespaceHeader.H"

GaussianSphereIF::GaussianSphereIF(const Real&     a_sigma,
                                   const RealVect& a_center,
                                   const bool&     a_inside)
{
  // Remember the parameters
  m_sigma  = a_sigma;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the sigma squared and the normalizing factor
  m_sigma2 = m_sigma * m_sigma;
  m_normal = 1.0 / (m_sigma * sqrt(2.0*M_PI));
}

GaussianSphereIF::GaussianSphereIF(const GaussianSphereIF& a_inputIF)
{
  // Remember the parameters
  m_sigma  = a_inputIF.m_sigma;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;

  m_sigma2 = a_inputIF.m_sigma2;
  m_normal = a_inputIF.m_normal;
}

GaussianSphereIF::~GaussianSphereIF()
{
}

void GaussianSphereIF::GetParams(Real&     a_sigma,
                                 RealVect& a_center,
                                 bool&     a_inside) const
{
  // Copy parameter information over
  a_sigma  = m_sigma;
  a_center = m_center;
  a_inside = m_inside;
}

void GaussianSphereIF::SetParams(const Real&     a_sigma,
                                 const RealVect& a_center,
                                 const bool&     a_inside)
{
  // Set parameter information
  m_sigma  = a_sigma;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the sigma squared and the normalizing factor
  m_sigma2 = m_sigma * m_sigma;
  m_normal = 1.0 / (m_sigma * sqrt(2.0*M_PI));
}

Real GaussianSphereIF::value(const RealVect& a_point) const
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

  retval = m_normal * exp(-distance2/(2*m_sigma2));

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

Real GaussianSphereIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  RealVect point;
  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      point[idir] = a_point[idir];
    }

  return value(point);
}

BaseIF* GaussianSphereIF::newImplicitFunction() const
{
  GaussianSphereIF* spherePtr = new GaussianSphereIF(m_sigma,
                                                     m_center,
                                                     m_inside);

  return static_cast<BaseIF*>(spherePtr);
}

#include "NamespaceFooter.H"
