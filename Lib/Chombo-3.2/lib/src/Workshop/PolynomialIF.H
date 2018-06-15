#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POLYNOMIALIF_H_
#define _POLYNOMIALIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"
#include "IndexTM.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This represents one term in a general polynomial
 */
class PolyTerm
{
public:
  ///
  /**
     Partial derivative
  */
  void partialDerivative(PolyTerm       & a_partial,
                         const IntVect  & a_whichPartial,
                         const PolyTerm & a_poly) const;

  ///
  /**
     First order partial
  */
  void firstOrderPartial(PolyTerm       & a_partial,
                         const int      & a_whichPartial) const;


  ///
  /**
     Coefficient of this polynomial term
  */
  Real coef;

  ///
  /**
     Powers of this polynomial term
  */
  IntVect powers;
};

///
/**
   This implicit function specifies a general polynomial function.
*/
class PolynomialIF: public BaseIF
{
public:
  ///
  /**
     Constructor specifying the term composing the polynomial (a_polynomial)
     and whether the domain is on the inside (a_inside), i.e. where the
     function is negative.
  */
  PolynomialIF(const Vector<PolyTerm>& a_polynomial,
               const bool&             a_inside);

  /// Copy constructor
  PolynomialIF(const PolynomialIF& a_inputIF);

  /// Destructor
  virtual ~PolynomialIF();

  ///
  /**
     Return the parameter information
  */
  virtual void GetParams(Vector<PolyTerm>& a_polynomial,
                         bool&             a_inside) const;

  ///
  /**
     Set the parameter information
  */
  virtual void SetParams(const Vector<PolyTerm>& a_polynomial,
                         const bool&             a_inside);

  ///
  /**
     Return the value of the function at a_point using a_polynomial
  */
  virtual Real value(const RealVect         & a_point,
                     const Vector<PolyTerm> & a_polynomial) const;

 ///
  /**
     Return the value of the function at a_point.
  */
  virtual Real value(const RealVect& a_point) const;

  ///
  /**
     Return the value of the derivative at a_point.
  */
  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;
  ///
  /**
     Partial derivative
  */
  virtual void partialDerivative(Vector<PolyTerm>       & a_partial,
                                 const IntVect          & a_whichPartialOp,
                                 const Vector<PolyTerm> & a_polynomial)const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  Vector<PolyTerm> m_polynomial;  // polynomial (term by term)
  bool             m_inside;      // inside flag

private:
  PolynomialIF()
  {
    MayDay::Abort("PolynomialIF uses strong construction");
  }

  void operator=(const PolynomialIF& a_inputIF)
  {
    MayDay::Abort("PolynomialIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
