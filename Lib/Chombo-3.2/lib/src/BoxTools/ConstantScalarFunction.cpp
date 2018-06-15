#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ConstantScalarFunction.H"
#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
ConstantScalarFunction::
ConstantScalarFunction(Real a_value):
  ScalarFunction(true, true),
  m_value(a_value)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
ConstantScalarFunction::
~ConstantScalarFunction()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
ConstantScalarFunction::
operator()(const RealVect& a_x, Real a_t) const
{
  return m_value;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
ConstantScalarFunction::
derivative(const IntVect& a_order,
           const RealVect& a_x,
           Real a_t) const
{
  // All derivatives are zero.
  if (a_order == IntVect::Zero)
    return m_value;
  else
    return 0.0;

}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
