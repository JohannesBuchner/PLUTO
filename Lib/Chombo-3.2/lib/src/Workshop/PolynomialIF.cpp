#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PolynomialIF.H"
#include "NamespaceHeader.H"

void PolyTerm::partialDerivative(PolyTerm       & a_partial,
                                 const IntVect  & a_whichPartialOp,
                                 const PolyTerm & a_polyterm) const
{
  a_partial = a_polyterm;

  // The direction and  order of differentation is summarzized by an IntVect
  for (int idir = 0; idir < SpaceDim; ++idir)
  {
    // Differentiate in the idir direction a_whichPartialOp number of times.
    int curPartialOp = a_whichPartialOp[idir];
    for (int ipartial = curPartialOp; ipartial > 0; --ipartial)
    {
      firstOrderPartial(a_partial,idir);
    }
  }
}

void PolyTerm::firstOrderPartial(PolyTerm       & a_partial,
                                 const int      & a_whichPartialDir) const
{
  // Take one derivative in the a_whichPartialDir direction
  a_partial.coef *= a_partial.powers[a_whichPartialDir];
  if (a_partial.powers[a_whichPartialDir] > 0)
  {
    a_partial.powers[a_whichPartialDir] -= 1;
  }
}

PolynomialIF::PolynomialIF(const Vector<PolyTerm>& a_polynomial,
                           const bool&             a_inside)
{
  // Copy polynomial
  int size = a_polynomial.size();

  m_polynomial.resize(size);

  for (int iterm = 0; iterm < size; iterm++)
    {
      m_polynomial[iterm] = a_polynomial[iterm];
    }

  // Save inside flag
  m_inside = a_inside;
}

PolynomialIF::PolynomialIF(const PolynomialIF& a_inputIF)
{
  // Copy polynomial
  int size = a_inputIF.m_polynomial.size();

  m_polynomial.resize(size);

  for (int iterm = 0; iterm < size; iterm++)
    {
      m_polynomial[iterm] = a_inputIF.m_polynomial[iterm];
    }

  // Save inside flag
  m_inside = a_inputIF.m_inside;
}

PolynomialIF::~PolynomialIF()
{
}

void PolynomialIF::GetParams(Vector<PolyTerm>& a_polynomial,
                             bool&             a_inside) const
{
  // Copy polynomial
  int size = m_polynomial.size();

  a_polynomial.resize(size);

  for (int iterm = 0; iterm < size; iterm++)
    {
      a_polynomial[iterm] = m_polynomial[iterm];
    }

  // Save inside flag
  a_inside = m_inside;
}

void PolynomialIF::SetParams(const Vector<PolyTerm>& a_polynomial,
                             const bool&             a_inside)
{
  // Copy polynomial
  int size = a_polynomial.size();

  m_polynomial.resize(size);

  for (int iterm = 0; iterm < size; iterm++)
    {
      m_polynomial[iterm] = a_polynomial[iterm];
    }

  // Save inside flag
  m_inside = a_inside;
}

Real PolynomialIF::value(const RealVect         & a_point,
                         const Vector<PolyTerm> & a_polynomial) const
{
  Real retval;

  int size = a_polynomial.size();

  // Evaluate the polynomial
  retval = 0.0;
  for (int iterm = 0; iterm < size; iterm++)
    {
      Real cur;

      cur = a_polynomial[iterm].coef;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          cur *= pow(a_point[idir],a_polynomial[iterm].powers[idir]);
        }

      retval += cur;
    }

  // Change the sign to change inside to outside
  if (!m_inside)
    {
      retval = -retval;
    }

  return retval;
}

Real PolynomialIF::value(const RealVect& a_point) const
{
  return value(a_point,m_polynomial);
}

Real PolynomialIF::value(const IndexTM<int,GLOBALDIM>  & a_partialDerivativeOp,
                         const IndexTM<Real,GLOBALDIM> & a_point) const
{
#ifdef CH_REALM
  // 0.0 is correct if a_partialDerivativeOp[GLOBALDIM] != 0
  Real retval = 0.0;

  Vector<PolyTerm> derivedPoly(m_polynomial);
  RealVect point;
  IntVect partialDerivativeOp;

  if (a_partialDerivativeOp[GLOBALDIM - 1] == 0)
    {
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          // Project on the SpaceDim dimensions, preliminary to calling the
          // generic value function. The reference surface polynomial has no z-terms.
          point[idir] = a_point[idir];

          // Project on the SpaceDim dimensions preliminary to calling the
          // generic partial derivative function.
          partialDerivativeOp[idir] = a_partialDerivativeOp[idir];
        }

      // Calculate the derivative of the polynomial
      partialDerivative(derivedPoly,partialDerivativeOp,m_polynomial);

      // Evaluate the derivative at the point in space
      retval = value(point,derivedPoly);
    }

  return retval;
#else
  MayDay::Abort("This function uses GLOBALDIM, which doesn't make sense outside of realm");

  return 0.0;
#endif
}

void PolynomialIF::partialDerivative(Vector<PolyTerm>       & a_partial,
                                     const IntVect          & a_whichPartialOp,
                                     const Vector<PolyTerm> & a_polynomial)const
{
  //The direction and order of differentiation is summarized by an IntVect
  int size = m_polynomial.size();

  for (int iterm = 0; iterm < size; iterm++)
    {
      a_partial[iterm].partialDerivative(a_partial[iterm],a_whichPartialOp,a_polynomial[iterm]);
    }
}

BaseIF* PolynomialIF::newImplicitFunction() const
{
  PolynomialIF* polynomialPtr = new PolynomialIF(m_polynomial,
                                                 m_inside);

  return static_cast<BaseIF*>(polynomialPtr);
}



#include "NamespaceFooter.H"
