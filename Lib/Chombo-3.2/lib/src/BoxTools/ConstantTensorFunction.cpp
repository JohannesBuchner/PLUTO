#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ConstantTensorFunction.H"
#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
ConstantTensorFunction::
ConstantTensorFunction(const RealTensor& a_value):
  TensorFunction(true, true),
  m_value(a_value)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
ConstantTensorFunction::
~ConstantTensorFunction()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RealTensor
ConstantTensorFunction::
operator()(const RealVect& a_x, Real a_t) const
{
  return m_value;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
