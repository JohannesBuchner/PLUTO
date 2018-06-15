#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TylerChannelIF.H"

#include "NamespaceHeader.H"

TylerChannelIF::TylerChannelIF(const Real& a_x1,
                 const Real& a_x2,
                 const Real& a_y1,
                 const Real& a_y2,
                 const Real& a_yDomainLength)
{
  m_x1 = a_x1;
  m_x2 = a_x2;

  m_y1 = a_y1;
  m_y2 = a_y2;

  m_yDomainLength = a_yDomainLength;
}

TylerChannelIF::TylerChannelIF(const TylerChannelIF& a_inputIF)
{
  m_x1 = a_inputIF.m_x1;
  m_x2 = a_inputIF.m_x2;

  m_y1 = a_inputIF.m_y1;
  m_y2 = a_inputIF.m_y2;

  m_yDomainLength = a_inputIF.m_yDomainLength;
}

TylerChannelIF::~TylerChannelIF()
{
}

Real TylerChannelIF::value(const RealVect& a_point) const
{
  Real retval;

  Real x = a_point[0];
  Real y = a_point[1];

  Real xmid = (m_x1 + m_x2) / 2.0;
  Real ymid = m_yDomainLength / 2.0;

  Real theta = 4.0 * (x - xmid) / (m_x2 - m_x1);
  Real frac = (tanh(theta) + 1.0) / 2.0;

  Real ydist = m_y1*(1.0-frac) + m_y2*frac;

  retval = Abs(y - ymid) - ydist;

  return retval;
}

BaseIF* TylerChannelIF::newImplicitFunction() const
{
  TylerChannelIF* TylerChannelPtr = new TylerChannelIF(m_x1,
                                                       m_x2,
                                                       m_y1,
                                                       m_y2,
                                                       m_yDomainLength);

  return static_cast<BaseIF*>(TylerChannelPtr);
}

#include "NamespaceFooter.H"
