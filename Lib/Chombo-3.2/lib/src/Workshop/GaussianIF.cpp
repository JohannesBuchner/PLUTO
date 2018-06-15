#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GaussianIF.H"

#include "NamespaceHeader.H"

GaussianIF::GaussianIF(const RealVect& a_origin,
                       const Real&     a_height,
                       const RealVect& a_sig2,
                       const int&      a_upDir,
                       const bool&     a_inside)
{
  // Remember the parameters
  m_origin = a_origin;
  m_height = a_height;
  m_sig2   = a_sig2;
  m_upDir  = a_upDir;
  m_inside = a_inside;
}

GaussianIF::GaussianIF(const GaussianIF& a_inputIF)
{
  // Remember the parameters
  m_origin = a_inputIF.m_origin;
  m_height = a_inputIF.m_height;
  m_sig2   = a_inputIF.m_sig2;
  m_upDir  = a_inputIF.m_upDir;
  m_inside = a_inputIF.m_inside;
}

GaussianIF::~GaussianIF()
{
}

void GaussianIF::GetParams(RealVect& a_origin,
                           Real&     a_height,
                           RealVect& a_sig2,
                           int&      a_upDir,
                           bool&     a_inside) const
{
  // Copy parameter information over
  a_origin = m_origin;
  a_height = m_height;
  a_sig2   = m_sig2;
  a_upDir  = m_upDir;
  a_inside = m_inside;
}

void GaussianIF::SetParams(const RealVect& a_origin,
                           const Real&     a_height,
                           const RealVect& a_sig2,
                           const int&      a_upDir,
                           const bool&     a_inside)
{
  // Set parameter information
  m_origin = a_origin;
  m_height = a_height;
  m_sig2   = a_sig2;
  m_upDir  = a_upDir;
  m_inside = a_inside;
}

Real GaussianIF::value(const RealVect& a_point) const
{
  Real retval = 0.0;

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != m_upDir)
        {
          retval += pow(a_point[idir]-m_origin[idir],2)/(-2.0*m_sig2[idir]);
        }
    }

  retval  = exp(retval);
  retval *= m_height;
  retval -= (a_point[m_upDir]-m_origin[m_upDir]);

  // Change the sign to change inside to outside
  if (!m_inside)
    {
      retval = -retval;
    }

  CH_assert(!isnan(retval) && !isinf(retval));

  return retval;
}

BaseIF* GaussianIF::newImplicitFunction() const
{
  GaussianIF* gaussianPtr = new GaussianIF(m_origin,
                                           m_height,
                                           m_sig2,
                                           m_upDir,
                                           m_inside);

  return static_cast<BaseIF*>(gaussianPtr);
}

#include "NamespaceFooter.H"
