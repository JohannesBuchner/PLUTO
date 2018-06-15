#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TorusIF.H"

#include "NamespaceHeader.H"

TorusIF::TorusIF(const Real&     a_majorRadius,
                 const Real&     a_minorRadius,
                 const RealVect& a_center,
                 const bool&     a_inside)
{
  // Remember parameters
  m_majorRadius = a_majorRadius;
  m_minorRadius = a_minorRadius;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the minor radius squared
  m_minorRadius2 = m_minorRadius * m_minorRadius;
}

TorusIF::TorusIF(const TorusIF& a_inputIF)
{
  // Remember parameters
  m_majorRadius = a_inputIF.m_majorRadius;
  m_minorRadius = a_inputIF.m_minorRadius;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;

  // Precompute the minor radius squared
  m_minorRadius2 = m_minorRadius * m_minorRadius;
}

TorusIF::~TorusIF()
{
}

void TorusIF::GetParams(Real&     a_majorRadius,
                        Real&     a_minorRadius,
                        RealVect& a_center,
                        bool&     a_inside) const
{
  // Copy parameter information over
  a_majorRadius = m_majorRadius;
  a_minorRadius = m_minorRadius;
  a_center = m_center;
  a_inside = m_inside;
}

void TorusIF::SetParams(const Real&     a_majorRadius,
                        const Real&     a_minorRadius,
                        const RealVect& a_center,
                        const bool&     a_inside)
{
  // Set parameter information
  m_majorRadius = a_majorRadius;
  m_minorRadius = a_minorRadius;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the minor radius squared
  m_minorRadius2 = m_minorRadius * m_minorRadius;
}

Real TorusIF::value(const RealVect& a_point) const
{
  Real retval;

  // Get the radius of a_point (excluding the z coordinate in 3D)
  Real radius1;

  radius1 = 0.0;
  for (int idir = 0; idir < 2; idir++)
  {
    Real cur;
    cur = a_point[idir] - m_center[idir];

    radius1 += cur*cur;
  }

  // Move to the center of the minor circle
  radius1 = sqrt(radius1) - m_majorRadius;

  // Get the radius squared in this coordinate system
  Real radius2;

  radius2 = radius1*radius1;
  for (int idir = 2; idir < SpaceDim; idir++)
  {
    Real cur;
    cur = a_point[idir] - m_center[idir];

    radius2 += cur*cur;
  }

  // Return the difference between the squares (zero on the torus)
  retval = radius2 - m_minorRadius2;

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* TorusIF::newImplicitFunction() const
{
  TorusIF* torusPtr = new TorusIF(m_majorRadius,
                                  m_minorRadius,
                                  m_center,
                                  m_inside);

  return static_cast<BaseIF*>(torusPtr);
}

#include "NamespaceFooter.H"
