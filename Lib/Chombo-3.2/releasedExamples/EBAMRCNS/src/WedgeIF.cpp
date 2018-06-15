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
#include "WedgeIF.H"

#include "NamespaceHeader.H"

WedgeIF::WedgeIF(const RealVect& a_center,
                 const Real&     a_slope,
                 const int&      a_depVar,
                 const int&      a_indVar)
{
  // Remember the parameters
  m_center= a_center;
  m_slope = a_slope ;
  m_depVar = a_depVar ;
  m_indVar = a_indVar ;
}

WedgeIF::WedgeIF(const WedgeIF& a_inputIF)
{
  // Remember the parameters
  m_center= a_inputIF.m_center;
  m_slope = a_inputIF.m_slope ;
  m_depVar = a_inputIF.m_depVar ;
  m_indVar = a_inputIF.m_indVar ;
}

WedgeIF::~WedgeIF()
{
}


Real WedgeIF::value(const RealVect& a_point) const
{
 CH_assert(m_depVar >= 0);
 CH_assert(m_depVar < SpaceDim);
 CH_assert(m_indVar >= 0);
 CH_assert(m_indVar < SpaceDim);
 CH_assert(m_indVar != m_depVar);

  Real indDist = a_point[m_indVar] - m_center[m_indVar];
  Real funcVal = m_center[m_depVar] - m_slope*Abs(indDist);
  Real retval  = funcVal - a_point[m_depVar];

  return retval;
}

BaseIF* WedgeIF::newImplicitFunction() const
{
  WedgeIF* wedgePtr = new WedgeIF(m_center,
                                  m_slope,
                                  m_depVar,
                                  m_indVar);

  return static_cast<BaseIF*>(wedgePtr);
}

#include "NamespaceFooter.H"
