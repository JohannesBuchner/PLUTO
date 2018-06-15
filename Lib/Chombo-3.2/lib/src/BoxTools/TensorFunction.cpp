#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TensorFunction.H"
#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
TensorFunction::
TensorFunction(bool a_homogeneous,
               bool a_constant):
  m_isHomogeneous(a_homogeneous),
  m_isConstant(a_constant)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
TensorFunction::
~TensorFunction()
{
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
