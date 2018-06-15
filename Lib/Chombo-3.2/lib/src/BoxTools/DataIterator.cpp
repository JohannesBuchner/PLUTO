#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DataIterator.H"
#include "NamespaceHeader.H"

#ifdef CH_MPI
DataIterator::DataIterator(const BoxLayout& plan,
                           const int* layoutID)
  :m_layout(plan), m_indices(plan.m_dataIndex), m_current(0)
{

}

void DataIterator::end()
{
  m_current = m_indices->size();
}

#endif
#include "NamespaceFooter.H"
