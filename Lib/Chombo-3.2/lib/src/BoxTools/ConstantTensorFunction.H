#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CONSTANTTENSORFUNCTION_H_
#define _CONSTANTTENSORFUNCTION_H_

#include "TensorFunction.H"
#include "NamespaceHeader.H"

//! \class ConstantVectorFunction
//! This VectorFunction subclass represents a function that returns a
//! constant value everywhere in space and time.
class ConstantTensorFunction: public TensorFunction
{
  public:

  //! Creates a constant function returning the given value.
  //! \param a_value The value returned by the function.
  explicit ConstantTensorFunction(const RealTensor& a_value);

  //! Destructor.
  ~ConstantTensorFunction();

  // Overridden methods.
  RealTensor operator()(const RealVect& a_x, Real a_t) const;

  private:

  RealTensor m_value;
};

#include "NamespaceFooter.H"
#endif
