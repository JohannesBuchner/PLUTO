#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "AMRLevel.H"
#include "EBAMRCNS.H"
#include "EBAMRCNSFactory.H"

#include "NamespaceHeader.H"

/****************/
AMRLevel*
EBAMRCNSFactory::
new_amrlevel() const
{
  EBAMRCNS* amrg_ptr = new EBAMRCNS(m_params, m_factory, m_ICs);
  //  amrgptr->setConservative(false);
  return (static_cast <AMRLevel*> (amrg_ptr));
}
/****************/

#include "NamespaceFooter.H"
