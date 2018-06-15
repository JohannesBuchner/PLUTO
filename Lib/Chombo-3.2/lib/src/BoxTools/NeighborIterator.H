#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEIGHBORITERATOR_H_
#define _NEIGHBORITERATOR_H_

#include "DataIndex.H"
#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "NamespaceHeader.H"

///An Iterator based on a DisjointBoxLayout object for neighboring boxes.
/**

*/
class NeighborIterator
{
public:

  NeighborIterator():m_dblPtr(NULL)
  {
  }

  ///
  NeighborIterator(const DisjointBoxLayout& dbl);

  /// return the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  inline const LayoutIndex& operator()() const ;

  Box box() const ;

  /// move the iterator to the next index in the neighbor list
  inline void operator++();

  /// return true if this iterator is still in the neighbor list
  inline bool ok() const;

  /// initialize this iterator to the first index in the neighbor list
  void begin(const DataIndex& a_dataIndex);

  ///  OK, this one requires a little explanation
  /**
     if the *current* box is a periodic image box, then this operation returns the unmapping
     of a_box.
  */
  Box unshift(const Box& a_box) const ;

private:

  const DisjointBoxLayout* m_dblPtr;
  std::vector<std::pair<int, LayoutIndex> >::const_iterator m_current, m_end;
  LayoutIndex m_lindex;
};

inline bool NeighborIterator::ok() const
{
  CH_assert(m_dblPtr!=NULL);
  return m_current != m_end;
}

inline void NeighborIterator::operator++()
{
  ++m_current;
  if (m_current != m_end)
    m_lindex = m_current->second;
}

inline const LayoutIndex& NeighborIterator::operator()() const
{
  CH_assert(m_dblPtr!= NULL);
  CH_assert(this->ok());
  return m_lindex;
}

#include "NamespaceFooter.H"
#endif
