#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "VectorFunction.H"
#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
VectorFunction::
VectorFunction(bool a_homogeneous,
               bool a_constant):
  m_isHomogeneous(a_homogeneous),
  m_isConstant(a_constant)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
VectorFunction::
~VectorFunction()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
RealVect
VectorFunction::
derivative(const IntVect& a_order,
           const RealVect& a_x,
           Real a_t) const
{
  return RealVect();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool
VectorFunction::
hasDerivative(const IntVect& a_order) const
{
  return false;
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
