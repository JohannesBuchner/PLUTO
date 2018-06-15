#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NeighborIterator.H"
#include "NamespaceHeader.H"

NeighborIterator::NeighborIterator(const DisjointBoxLayout& dbl)
  :m_dblPtr(&dbl), m_lindex(-1, -1, dbl.m_layout)
{
}

void NeighborIterator::begin(const DataIndex& a_dataIndex)
{
  const Vector<Vector<std::pair<int, LayoutIndex> > >&  vv = *(m_dblPtr->m_neighbors);
  const std::vector<std::pair<int,   LayoutIndex> >& v = vv[a_dataIndex.intCode()].constStdVector();
  m_current = v.begin();
  m_end = v.end();
  if (ok())
    {
       m_lindex = m_current->second;
    }
}

Box NeighborIterator::box() const
{
  CH_assert(m_dblPtr != NULL);
  Box rtn=m_dblPtr->get(m_lindex);
  if (m_current->first == -1) return rtn;
  m_dblPtr->physDomain().shiftIt(rtn, m_current->first);
  return rtn;
}

Box NeighborIterator::unshift(const Box& a_box) const
{
  if (m_current->first == -1) return a_box;
  Box rtn(a_box);
  m_dblPtr->physDomain().unshiftIt(rtn, m_current->first);
  return rtn;
}

#include "NamespaceFooter.H"
