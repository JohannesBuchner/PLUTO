#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RhodoneaIF.H"

#include "NamespaceHeader.H"

RhodoneaIF::RhodoneaIF(const Real&     a_innerRadius,
                       const Real&     a_outerRadius,
                       const int&      a_frequency,
                       const RealVect& a_center,
                       const bool&     a_inside)
{
  CH_assert(SpaceDim == 2);
  CH_assert(a_innerRadius > 0.);
  CH_assert(a_outerRadius > 0.);
  CH_assert(a_frequency != 0);

  // Remember the parameters
  m_innerRadius = a_innerRadius;
  m_outerRadius = a_outerRadius;
  m_frequency   = a_frequency;
  m_center = a_center;
  m_inside = a_inside;
}

RhodoneaIF::RhodoneaIF(const RhodoneaIF& a_inputIF)
{
  // Remember the parameters
  m_innerRadius = a_inputIF.m_innerRadius;
  m_outerRadius = a_inputIF.m_outerRadius;
  m_frequency   = a_inputIF.m_frequency;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;
}

RhodoneaIF::~RhodoneaIF()
{
  CH_assert(SpaceDim == 2);
}

void RhodoneaIF::GetParams(Real&     a_innerRadius,
                           Real&     a_outerRadius,
                           int&      a_frequency,
                           RealVect& a_center,
                           bool&     a_inside) const
{
  // Copy parameter information over
  a_innerRadius = m_innerRadius;
  a_outerRadius = m_outerRadius;
  a_frequency   = m_frequency;
  a_center = m_center;
  a_inside = m_inside;
}

void RhodoneaIF::SetParams(const Real&     a_innerRadius,
                           const Real&     a_outerRadius,
                           const int&      a_frequency,
                           const RealVect& a_center,
                           const bool&     a_inside)
{
  CH_assert(a_innerRadius > 0.);
  CH_assert(a_outerRadius > 0.);
  CH_assert(a_frequency != 0);

  // Set parameter information
  m_innerRadius = a_innerRadius;
  m_outerRadius = a_outerRadius;
  m_frequency   = a_frequency;
  m_center = a_center;
  m_inside = a_inside;
}

Real RhodoneaIF::value(const RealVect& a_point) const
{
  Real retval;

  // Compute polar coordinates of this point, relative to center
  Real rad = 0.;
  for (int idir=0; idir<SpaceDim; idir++)
    {
      Real cur = a_point[idir]-m_center[idir];
      rad += cur*cur;
    }
  rad = sqrt(rad);
  Real theta = atan2(a_point[1]-m_center[1],a_point[0]-m_center[0]);

  // Compute the radius of the rhodonea at this theta
  Real rRhod = m_innerRadius + m_outerRadius*sin(m_frequency*theta);

  // The ratio should be 1.0 on the surface of the rhodonea
  Real ratio = rad/rRhod;

  retval = ratio - 1.0;

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* RhodoneaIF::newImplicitFunction() const
{
  RhodoneaIF* ellipsoidPtr = new RhodoneaIF(m_innerRadius,
                                            m_outerRadius,
                                            m_frequency,
                                            m_center,
                                            m_inside);

  return static_cast<BaseIF*>(ellipsoidPtr);
}

#include "NamespaceFooter.H"
