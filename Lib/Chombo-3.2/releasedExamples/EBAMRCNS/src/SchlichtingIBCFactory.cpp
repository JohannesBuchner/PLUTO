#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SchlichtingIBCFactory.H"
#include "SchlichtingIBC.H"

#include "NamespaceHeader.H"

/******************/
SchlichtingIBCFactory::
SchlichtingIBCFactory(const SchlichtingParams& a_params)
  :EBPhysIBCFactory()
{
  m_params = a_params;
}
/******************/
SchlichtingIBCFactory::
~SchlichtingIBCFactory()
{;}
/******************/
EBPhysIBC*
SchlichtingIBCFactory::
create() const
{
  SchlichtingIBC* retval =
    new SchlichtingIBC(m_params);

  return static_cast<EBPhysIBC*>(retval);
}
/******************/

#include "NamespaceFooter.H"
