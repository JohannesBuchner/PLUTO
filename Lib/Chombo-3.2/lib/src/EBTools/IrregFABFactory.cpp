#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL
#include "IrregFABFactory.H"
#include "NamespaceHeader.H"
/***************/
/***************/
IrregFABFactory::~IrregFABFactory()
{
}
/***************/
/***************/
IrregFABFactory::IrregFABFactory(const EBISLayout& a_ebisl)
{
  m_ebisl = a_ebisl;
}
/***************/
/***************/
IrregFAB*
IrregFABFactory::create(const Box& a_box, int a_ncomps,
                      const DataIndex& a_dit) const
{
  return new IrregFAB(a_box, m_ebisl[a_dit].getEBGraph(),  a_ncomps);
}
/***************/
/***************/
#include "NamespaceFooter.H"
