#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CONSTANTVECTORFUNCTION_H_
#define _CONSTANTVECTORFUNCTION_H_

#include "VectorFunction.H"
#include "NamespaceHeader.H"

//! \class ConstantVectorFunction
//! This VectorFunction subclass represents a function that returns a
//! constant value everywhere in space and time.
class ConstantVectorFunction: public VectorFunction
{
  public:

  //! Creates a constant scalar function returning the given value.
  //! \param a_value The value returned by the function.
  explicit ConstantVectorFunction(const RealVect& a_value);

  //! Destructor.
  ~ConstantVectorFunction();

  // Overridden methods.
  RealVect operator()(const RealVect& a_x, Real a_t) const;

  RealVect derivative(const IntVect& a_order,
                      const RealVect& a_x,
                      Real a_t) const;

  private:

  RealVect m_value;
};

#include "NamespaceFooter.H"
#endif
