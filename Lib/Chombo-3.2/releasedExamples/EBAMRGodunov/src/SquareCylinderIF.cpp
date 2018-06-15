#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "RealVect.H"
#include "PolyGeom.H"
#include "SquareCylinderIF.H"

#include "NamespaceHeader.H"

SquareCylinderIF::SquareCylinderIF(const Real&     a_radius,
                                   const bool&     a_inside)
{
  // Remember the parameters
  m_radius = a_radius;
  m_inside = a_inside;
}

SquareCylinderIF::SquareCylinderIF(const SquareCylinderIF& a_inputIF)
{
  // Remember the parameters
  m_radius = a_inputIF.m_radius;
  m_inside = a_inputIF.m_inside;
}

SquareCylinderIF::~SquareCylinderIF()
{
}

void SquareCylinderIF::GetParams(Real&     a_radius,
                                 bool&     a_inside) const
{
  // Copy parameter information over
  a_radius = m_radius;
  a_inside = m_inside;
}

void SquareCylinderIF::SetParams(const Real&     a_radius,
                                 const bool&     a_inside)
{
  // Set parameter information
  m_radius = a_radius;
  m_inside = a_inside;
}

Real SquareCylinderIF::value(const RealVect& a_point) const
{
  Real dist = 0.0;
  //dist is max dist along coord line from x axis (y = 0, z = 0 line)
  for (int idir = 1; idir < SpaceDim; idir++)
    {
      Real tandist  =  Abs(a_point[idir]);
      dist = Max(dist, tandist);
    }

 CH_assert(m_radius > 0);
  Real retval = dist - m_radius;

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* SquareCylinderIF::newImplicitFunction() const
{
  SquareCylinderIF* squareCylinderPtr = new SquareCylinderIF(m_radius,
                                                             m_inside);

  return static_cast<BaseIF*>(squareCylinderPtr);
}

#include "NamespaceFooter.H"
