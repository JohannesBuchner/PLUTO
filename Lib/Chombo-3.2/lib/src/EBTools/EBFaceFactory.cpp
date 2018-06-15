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
#include "EBFaceFactory.H"
#include "NamespaceHeader.H"
/***************/
/***************/
EBFaceFactory::~EBFaceFactory()
{
}
/***************/
/***************/
EBFaceFactory::EBFaceFactory(const EBISLayout& a_ebisl,
                             const int& a_idir)
{
  m_ebisl = a_ebisl;
  m_idir = a_idir;
}
/***************/
/***************/
EBFaceFAB*
EBFaceFactory::create(const Box& a_box, int a_ncomps,
                      const DataIndex& a_dit) const
{
  return new EBFaceFAB(m_ebisl[a_dit], a_box, m_idir, a_ncomps);
}
/***************/
/***************/
#include "NamespaceFooter.H"
