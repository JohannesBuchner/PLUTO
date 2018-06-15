#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ConstantVectorFunction.H"
#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
ConstantVectorFunction::
ConstantVectorFunction(const RealVect& a_value):
  VectorFunction(true, true),
  m_value(a_value)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
ConstantVectorFunction::
~ConstantVectorFunction()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RealVect
ConstantVectorFunction::
operator()(const RealVect& a_x, Real a_t) const
{
  return m_value;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RealVect
ConstantVectorFunction::
derivative(const IntVect& a_order,
           const RealVect& a_x,
           Real a_t) const
{
  // All derivatives are zero.
  if (a_order == IntVect::Zero)
    return m_value;
  else
    return RealVect::Zero;

}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
