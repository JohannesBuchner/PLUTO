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
#include "IntersectionIF.H"
#include "MultiSphereIF.H"

#include "NamespaceHeader.H"

MultiSphereIF::MultiSphereIF(const Vector<Real>&     a_radii,
                             const Vector<RealVect>& a_centers,
                             const bool&             a_inside)
{
  // Remember the number of spheres
  m_numSpheres = a_radii.size();

  // Remember the parameters
  m_radii   = a_radii;
  m_centers = a_centers;
  m_inside  = a_inside;

  // Create all the spheres
  Vector<BaseIF*> spheres(m_numSpheres);
  for (int i = 0; i < m_numSpheres; i++)
  {
    spheres[i] = new SphereIF(m_radii[i],m_centers[i],false);
  }

  // Intersect their exteriors and complement if needed
  IntersectionIF multiSphere(spheres);
  m_multiSphere = new ComplementIF(multiSphere,m_inside);

  // Delete all the spheres (they've been copied into m_multiSphere)
  for (int i = 0; i < m_numSpheres; i++)
  {
    delete spheres[i];
  }
}

MultiSphereIF::MultiSphereIF(const MultiSphereIF& a_inputIF)
{
  // Remember the parameters
  m_numSpheres = a_inputIF.m_numSpheres;
  m_radii      = a_inputIF.m_radii;
  m_centers    = a_inputIF.m_centers;
  m_inside     = a_inputIF.m_inside;

  if (a_inputIF.m_multiSphere == NULL)
  {
    m_multiSphere = NULL;
  }
  else
  {
    m_multiSphere = (ComplementIF *)a_inputIF.m_multiSphere->newImplicitFunction();
  }
}

MultiSphereIF::~MultiSphereIF()
{
  // Delete the IF object
  if (m_multiSphere != NULL)
  {
    delete m_multiSphere;
  }
}

void MultiSphereIF::GetParams(Vector<Real>&     a_radii,
                              Vector<RealVect>& a_centers,
                              bool&             a_inside) const
{
  // Copy parameter information over
  a_radii   = m_radii;
  a_centers = m_centers;
  a_inside  = m_inside;
}

void MultiSphereIF::SetParams(const Vector<Real>&     a_radii,
                              const Vector<RealVect>& a_centers,
                              const bool&             a_inside)
{
  // Delete the IF object
  if (m_multiSphere != NULL)
  {
    delete m_multiSphere;
  }

  // Set the number of spheres
  m_numSpheres = a_radii.size();

  // Set parameter information
  m_radii   = a_radii;
  m_centers = a_centers;
  m_inside  = a_inside;

  // Create all the spheres
  Vector<BaseIF*> spheres(m_numSpheres);
  for (int i = 0; i < m_numSpheres; i++)
  {
    spheres[i] = new SphereIF(m_radii[i],m_centers[i],false);
  }

  // Intersect their exteriors and complement if needed
  IntersectionIF multiSphere(spheres);
  m_multiSphere = new ComplementIF(multiSphere,m_inside);

  // Delete all the spheres (they've been copied into m_multiSphere)
  for (int i = 0; i < m_numSpheres; i++)
  {
    delete spheres[i];
  }
}

Real MultiSphereIF::value(const RealVect& a_point) const
{
  Real retval;

  retval = 0.0;

  if (m_multiSphere != NULL)
  {
    retval = m_multiSphere->value(a_point);
  }

  return retval;
}

BaseIF* MultiSphereIF::newImplicitFunction() const
{
  MultiSphereIF* spherePtr = new MultiSphereIF(m_radii,
                                               m_centers,
                                               m_inside);

  return static_cast<BaseIF*>(spherePtr);
}


GeometryService::InOut MultiSphereIF::InsideOutside(const RealVect& a_low, const RealVect& a_high) const
{
  return m_multiSphere->InsideOutside(a_low, a_high);
}

#include "NamespaceFooter.H"
