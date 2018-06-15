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

#ifndef _BASEIFFACTORYI_H_
#define _BASEIFFACTORYI_H_
#include "NamespaceHeader.H"

/***************/
/***************/
template <class T> inline
BaseIFFactory<T>::~BaseIFFactory()
{
}
/***************/
/***************/
template <class T> inline
BaseIFFactory<T>::BaseIFFactory(const EBISLayout& a_ebisl,
                                const LayoutData<IntVectSet>& a_sets,
                                const int& a_facedir)
{
  m_ebisl = a_ebisl;
  m_facedir = a_facedir;
  m_sets.define(a_sets.boxLayout());
  for (DataIterator dit = m_sets.dataIterator(); dit.ok(); ++dit)
    {
      m_sets[dit()] = a_sets[dit()];
    }
}
/***************/
/***************/
template <class T> inline
BaseIFFAB<T>*
BaseIFFactory<T>::create(const Box& a_box, int a_ncomps,
                         const DataIndex& a_dit) const
{
  //intersects with input box because leveldata will expect
  //that anyway when it tries to do a copy or a linearization.
  IntVectSet ivsInter = m_sets[a_dit];
  ivsInter &= a_box;
  return new BaseIFFAB<T>(ivsInter, m_ebisl[a_dit].getEBGraph(), m_facedir, a_ncomps);
}
/***************/
/***************/

#include "NamespaceFooter.H"
#endif
