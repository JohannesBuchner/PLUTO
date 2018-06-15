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

#ifndef _BASEIVFACTORYI_H_
#define _BASEIVFACTORYI_H_
#include "NamespaceHeader.H"

/***************/
/***************/
template <class T> inline
BaseIVFactory<T>::~BaseIVFactory()
{
}
/***************/
/***************/
template <class T> inline
BaseIVFactory<T>::BaseIVFactory(const EBISLayout& a_ebisl,
                                const LayoutData<IntVectSet>& a_sets)
{
  m_ebisl = a_ebisl;
  m_sets.define(a_sets.boxLayout());
  for (DataIterator dit = m_sets.dataIterator(); dit.ok(); ++dit)
    {
      m_sets[dit()] = a_sets[dit()];
    }
}

template <class T> inline
BaseIVFactory<T>::BaseIVFactory(const EBISLayout& a_ebisl)
{
  m_ebisl = a_ebisl;
  m_sets.define(a_ebisl.getDisjointLayout());
  DisjointBoxLayout layout = a_ebisl.getDisjointLayout();
  for (DataIterator dit = m_sets.dataIterator(); dit.ok(); ++dit)
    {
      Box region =  m_ebisl[dit()].getRegion();
      m_sets[dit()] = m_ebisl[dit()].boundaryIVS(region);
    }
}
/***************/
/***************/
template <class T> inline
BaseIVFAB<T>*
BaseIVFactory<T>::create(const Box& a_box, int a_ncomps,
                         const DataIndex& a_dit) const
{
  //intersects with input box because leveldata will expect
  //that anyway when it tries to do a copy or a linearization.
  IntVectSet ivsInter = m_sets[a_dit];
  ivsInter &= a_box;
  return new BaseIVFAB<T>(ivsInter, m_ebisl[a_dit].getEBGraph(), a_ncomps);
}
/***************/
/***************/

#include "NamespaceFooter.H"
#endif
