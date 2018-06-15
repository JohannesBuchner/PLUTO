#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ScalarFunction.H"
#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
ScalarFunction::
ScalarFunction(bool a_homogeneous,
               bool a_constant):
  m_isHomogeneous(a_homogeneous),
  m_isConstant(a_constant)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
ScalarFunction::
~ScalarFunction()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real
ScalarFunction::
derivative(const IntVect& a_order,
           const RealVect& a_x,
           Real a_t) const
{
  return 0.0;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
ScalarFunction::
hasDerivative(const IntVect& a_order) const
{
  return false;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
