#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MollifyIF.H"

#include "NamespaceHeader.H"

MollifyIF::MollifyIF(const BaseIF& a_impFunc,
                     const BaseIF& a_mollifier,
                     const Real&   a_min,
                     const Real&   a_max,
                     const int&    a_numPts)
{
  // Copy the implicit functions
  m_impFunc   = a_impFunc.newImplicitFunction();
  m_mollifier = a_mollifier.newImplicitFunction();

  // Save the sampled neighborhood specifications
  m_min    = a_min;
  m_max    = a_max;
  m_numPts = a_numPts;

  // Spacing of the sample points
  if (m_numPts > 1)
  {
    m_dx = (m_max - m_min) / (m_numPts - 1);
  }
  else
  {
    m_dx = 0.0;
  }

  // Box for the sampled mollifier
  m_sampleBox.define(IntVect::Zero, (m_numPts-1) * IntVect::Unit);

  // Create the FArrayBox for the sampled mollifier
  int nComp = 1;

  m_sampledMollifier.define(m_sampleBox,nComp);

  // Sum the mollifier values so they can be normalized
  m_mollifierSum = 0.0;

  // Sample the mollifier on the necessary grid
  for (BoxIterator bit(m_sampleBox); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect point;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      point[idir] = m_max - iv[idir] * m_dx;
    }

    Real value = m_mollifier->value(point);

    m_mollifierSum += value;
    m_sampledMollifier(iv,0) = value;
  }

  // Normalize the sampled mollifier to be a partition of unity (i.e. to sum
  // to 1.0)
  m_sampledMollifier *= 1.0 / m_mollifierSum;
}

MollifyIF::MollifyIF(const MollifyIF& a_inputIF)
{
  // Copy things over
  if (a_inputIF.m_impFunc == NULL)
  {
    m_impFunc = NULL;
  }
  else
  {
    m_impFunc = a_inputIF.m_impFunc->newImplicitFunction();
  }

  if (a_inputIF.m_mollifier == NULL)
  {
    m_mollifier = NULL;
  }
  else
  {
    m_mollifier = a_inputIF.m_mollifier->newImplicitFunction();
  }

  m_min    = a_inputIF.m_min;
  m_max    = a_inputIF.m_max;
  m_numPts = a_inputIF.m_numPts;

  m_dx     = a_inputIF.m_dx;

  m_sampleBox = a_inputIF.m_sampleBox;

  m_sampledMollifier.define(a_inputIF.m_sampledMollifier.box(),
                            a_inputIF.m_sampledMollifier.nComp());
  m_sampledMollifier.copy  (a_inputIF.m_sampledMollifier);
}

MollifyIF::~MollifyIF()
{
  // Delete the IFs (if they exist)
  if (m_impFunc != NULL)
  {
    delete m_impFunc;
  }

  if (m_mollifier != NULL)
  {
    delete m_mollifier;
  }
}

Real MollifyIF::value(const RealVect& a_point) const
{
  Real retval = 0.0;

  // Apply the mollifier to the function at the current point (a_point)
  for (BoxIterator bit(m_sampleBox); bit.ok(); ++bit)
  {
    IntVect iv = bit();
    RealVect point(a_point);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      point[idir] += m_min + iv[idir] * m_dx;
    }

    retval += m_sampledMollifier(iv,0) * m_impFunc->value(point);
  }

  return retval;
}

Real MollifyIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  RealVect point;
  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      point[idir] = a_point[idir];
    }

  return value(point);
}

BaseIF* MollifyIF::newImplicitFunction() const
{
  MollifyIF* mollifyPtr = new MollifyIF(*m_impFunc, *m_mollifier,
                                        m_min, m_max, m_numPts);

  return static_cast<BaseIF*>(mollifyPtr);
}

#include "NamespaceFooter.H"
