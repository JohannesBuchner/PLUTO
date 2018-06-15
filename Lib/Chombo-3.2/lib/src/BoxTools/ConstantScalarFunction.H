#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CONSTANTSCALARFUNCTION_H_
#define _CONSTANTSCALARFUNCTION_H_

#include "ScalarFunction.H"
#include "NamespaceHeader.H"

//! \class ConstantScalarFunction
//! This ScalarFunction subclass represents a function that returns a
//! constant value everywhere in space and time.
class ConstantScalarFunction: public ScalarFunction
{
  public:

  //! Creates a constant scalar function returning the given value.
  //! \param a_value The value returned by the function.
  explicit ConstantScalarFunction(Real a_value);

  //! Destructor.
  ~ConstantScalarFunction();

  // Overridden methods.
  Real operator()(const RealVect& a_x, Real a_t) const;

  Real derivative(const IntVect& a_order,
                  const RealVect& a_x,
                  Real a_t) const;

  private:

  Real m_value;
};

#include "NamespaceFooter.H"
#endif
