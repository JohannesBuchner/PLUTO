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

#ifndef _BASEEBCELLFACTORYI_H_
#define _BASEEBCELLFACTORYI_H_
#include "NamespaceHeader.H"

/***************/
/***************/
template <class T>
BaseEBCellFactory<T>::~BaseEBCellFactory()
{
}
/***************/
/***************/
template <class T>
BaseEBCellFactory<T>::BaseEBCellFactory(const EBISLayout& a_ebisl)
{
  m_ebisl = a_ebisl;
}
/***************/
/***************/
template <class T>
BaseEBCellFAB<T>*
BaseEBCellFactory<T>::create(const Box& a_box, int a_ncomps,
                             const DataIndex& a_dit) const
{
  return new BaseEBCellFAB<T>(m_ebisl[a_dit], a_box, a_ncomps);
}
/***************/
/***************/

#include "NamespaceFooter.H"
#endif
