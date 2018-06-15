#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ArteryIF.H"

#include "NamespaceHeader.H"

ArteryIF::ArteryIF(const int&      a_type,
                   const Real&     a_radius,
                   const Real&     a_amplitude,
                   const Real&     a_start,
                   const Real&     a_end,
                   const RealVect& a_center,
                   const bool&     a_inside)
{
  if (a_type != 1 && a_type != 2)
  {
    MayDay::Error("ArteryIF::ArteryIF: Artery type must be 1 or 2");
  }

  // Remember the parameters
  m_type      = a_type;
  m_radius    = a_radius;
  m_amplitude = a_amplitude;
  m_start     = a_start;
  m_end       = a_end;
  m_center    = a_center;
  m_inside    = a_inside;
}

ArteryIF::ArteryIF(const ArteryIF& a_inputIF)
{
  // Remember the parameters
  m_type      = a_inputIF.m_type;
  m_radius    = a_inputIF.m_radius;
  m_amplitude = a_inputIF.m_amplitude;
  m_start     = a_inputIF.m_start;
  m_end       = a_inputIF.m_end;
  m_center    = a_inputIF.m_center;
  m_inside    = a_inputIF.m_inside;
}

ArteryIF::~ArteryIF()
{
}

Real ArteryIF::value(const RealVect& a_point) const
{
  Real retval;

  Real x = a_point[0];
  Real xbar;

  // Adjust x so that it is a constant 0.0 if x < m_start, 1.0 if x >= x_end,
  // and varies linearly between 0.0 and 1.0 if m_start <= x < x_end
  if (x < m_start)
  {
    xbar = 0.0;
  }
  else if (x < m_end)
  {
    xbar = (x - m_start) / (m_end - m_start);
  }
  else
  {
    xbar = 1.0;
  }

  Real func = 0.5 * (cos(2*M_PI*    xbar)
                   - cos(2*M_PI*3.5*xbar));

  Real perturb = m_amplitude * func;

  if (m_type == 1)
  {
    Real radius = 0.0;

    for (int idir = 1; idir < SpaceDim; idir++)
    {
      radius += (a_point[idir]-0.5)*(a_point[idir]-0.5);
    }

    radius = sqrt(radius);

    retval = radius - (m_radius + perturb);
  }

  if (m_type == 2)
  {
    Real radius = 0.0;

    if (SpaceDim == 2)
    {
      radius += ((a_point[1]-m_center[1]) - perturb*sin(2*M_PI*xbar))
              * ((a_point[1]-m_center[1]) - perturb*sin(2*M_PI*xbar));
    }

    if (SpaceDim == 3)
    {
      radius += ((a_point[1]-m_center[1]) - perturb*cos(2*M_PI*xbar))
              * ((a_point[1]-m_center[1]) - perturb*cos(2*M_PI*xbar));

      radius += ((a_point[2]-m_center[2]) - perturb*sin(2*M_PI*xbar))
              * ((a_point[2]-m_center[2]) - perturb*sin(2*M_PI*xbar));
    }

    radius = sqrt(radius);

    retval = radius - (m_radius + perturb);
  }

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* ArteryIF::newImplicitFunction() const
{
  ArteryIF* arteryPtr = new ArteryIF(m_type,
                                     m_radius,
                                     m_amplitude,
                                     m_start,
                                     m_end,
                                     m_center,
                                     m_inside);

  return static_cast<BaseIF*>(arteryPtr);
}

#include "NamespaceFooter.H"
