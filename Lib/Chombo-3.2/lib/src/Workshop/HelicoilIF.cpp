#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EllipsoidIF.H"
#include "HelicoilIF.H"

#include "NamespaceHeader.H"

HelicoilIF::HelicoilIF(const Real& a_helixR,
                       const Real& a_helixPitch,
                       const Real& a_circleR,
                       const bool& a_inside)
{
  // Remember the parameters
  m_helixR     = a_helixR;
  m_helixPitch = a_helixPitch;
  m_circleR    = a_circleR;

  // Save inside flag
  m_inside = a_inside;

  Real slope = m_helixPitch / (2*M_PI*m_helixR);

  if (slope > 1.0)
  {
    m_ellipseR = m_circleR / sin(atan(slope));
    m_vertical = false;
  }
  else
  {
    m_ellipseR = m_circleR / cos(atan(slope));
    m_vertical = true;
  }

  RealVect radii(D_DECL(m_circleR,m_ellipseR,1.0));
  RealVect center(D_DECL(m_helixR,0.0,0.0));
  bool inside = true;

  EllipsoidIF crossSection(radii,center,inside);
  Real rate = 2*M_PI / m_helixPitch;

  m_helixIF = new HelixIF(crossSection,rate,m_inside,m_vertical);
}

HelicoilIF::HelicoilIF(const HelicoilIF& a_inputIF)
{
  // Remember the parameters
  m_helixR     = a_inputIF.m_helixR;
  m_helixPitch = a_inputIF.m_helixPitch;
  m_circleR    = a_inputIF.m_circleR;

  // Save inside flag
  m_inside = a_inputIF.m_inside;

  m_ellipseR = a_inputIF.m_ellipseR;
  m_vertical = a_inputIF.m_vertical;

  m_helixIF = a_inputIF.m_helixIF->newImplicitFunction();
}

HelicoilIF::~HelicoilIF()
{
  delete m_helixIF;
}

Real HelicoilIF::value(const RealVect& a_point) const
{
  Real retval = m_helixIF->value(a_point);

  return retval;
}

BaseIF* HelicoilIF::newImplicitFunction() const
{
  HelicoilIF* helicoilPtr;

  helicoilPtr = new HelicoilIF(m_helixR,m_helixPitch,m_circleR,m_inside);

  return static_cast<BaseIF*>(helicoilPtr);
}

#include "NamespaceFooter.H"
