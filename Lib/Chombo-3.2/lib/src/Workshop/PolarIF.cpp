#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PolarIF.H"

#include "NamespaceHeader.H"

PolarIF::PolarIF(const Real& a_primaryRadius,
                 const Real& a_perturbation,
                 const int & a_frequency,
                 const bool& a_inside)
{
  // Remember the parameters
  m_primaryRadius = a_primaryRadius;
  m_perturbation  = a_perturbation;
  m_frequency     = a_frequency;
  m_inside        = a_inside;
}

PolarIF::PolarIF(const PolarIF& a_inputIF)
{
  // Remember the parameters
  m_primaryRadius = a_inputIF.m_primaryRadius;
  m_perturbation  = a_inputIF.m_perturbation;
  m_frequency     = a_inputIF.m_frequency;
  m_inside        = a_inputIF.m_inside;
}

PolarIF::~PolarIF()
{
}

void PolarIF::GetParams(Real& a_primaryRadius,
               Real& a_perturbation,
               int & a_frequency,
               bool& a_inside)const
{
  // Copy parameter information over
  a_primaryRadius = m_primaryRadius;
  a_perturbation  = m_perturbation;
  a_frequency     = m_frequency;
  a_inside        = m_inside;
 }

void PolarIF::SetParams(const Real& a_primaryRadius,
               const Real& a_perturbation,
               const int & a_frequency,
               const bool& a_inside)
{
  // Set the parameters
  m_primaryRadius = a_primaryRadius;
  m_perturbation  = a_perturbation;
  m_frequency     = a_frequency;
  m_inside        = a_inside;
}

Real PolarIF::value(const RealVect& a_point) const
{
  Real retval;

  Real radius;
  Real theta;

  radius = sqrt(a_point[0]*a_point[0] + a_point[1]*a_point[1]);
  theta = atan2(a_point[1],a_point[0]);

  retval = radius - (m_primaryRadius + m_perturbation*cos(m_frequency*theta));

  if (!m_inside)
    {
      retval = -retval;
    }
  return retval;
}

BaseIF* PolarIF::newImplicitFunction() const
{
  PolarIF* polarPtr = new PolarIF( m_primaryRadius,
                                   m_perturbation,
                                   m_frequency,
                                   m_inside);

  return static_cast<BaseIF*>(polarPtr);
}

#include "NamespaceFooter.H"
