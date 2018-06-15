#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DATAINDEX_H_
#define _DATAINDEX_H_

#include <cstdlib>
#include <iostream>
#include "MayDay.H"
#include "parstream.H"
#include "NamespaceHeader.H"

#ifndef WRAPPER
class BoxLayout;

class DataIterator;
class LayoutIterator;

#endif /* WRAPPER */

using std::endl;
///An index for LayoutIterator
class LayoutIndex
{
public:
  ///
  LayoutIndex();

  virtual ~LayoutIndex()
  {}

  // default copy and assign should be fine.

  ///
  bool operator == (const LayoutIndex& a_rhs) const;

  bool eq(const LayoutIndex& a_rhs) const
  {
    return *this == a_rhs;
  }

  ///
  bool operator != (const LayoutIndex& a_rhs) const
  {
    return !(*this == a_rhs);
  }

  /// returns 'true' if this DataIndex has been null constructed.
  bool isNull() const;

  // not user function. breaks design encapsulation. method not intended for
  // users. The alternative is to make this class a friend of BoxLayoutData<T>,
  // which implicitly makes this class templated, which makes no bloody sense.
  // The behaviour of this method is undefined for users.  bvs
  int intCode() const;

  int datInd() const
  {
#ifdef CH_MPI
    return m_datInd;
#else
    //no MPI means dataindex is just a layout index
    return m_index;
#endif
  }

  LayoutIndex(const LayoutIndex& a_input)
  {
    m_index        = a_input.m_index;
    m_datInd       = a_input.m_datInd;
    m_layoutIntPtr = a_input.m_layoutIntPtr;
  }
  ///
  /**
     Print to the given output stream in ASCII.
  */
  virtual void output() const
  {
    pout() << "(" << m_index << "," << m_datInd << ")" << endl;
  }

private:
  friend class DataIterator;
  friend class LayoutIterator;
  friend class TimedDataIterator;
  friend class NeighborIterator;
  friend class BoxLayout;
  friend class DataIndex;
  friend class DisjointBoxLayout;

  int        m_index;
  int        m_datInd;
  const int* m_layoutIntPtr;

  LayoutIndex(int        a_indexIntoBoxes,
              int        a_indexIntoData,
              const int* a_layoutID)
    :
    m_index (a_indexIntoBoxes),
    m_datInd(a_indexIntoData),
    m_layoutIntPtr(a_layoutID)
  {}
};

class DataIndex : public LayoutIndex
{
public:
  explicit DataIndex(const LayoutIndex& a_promotion)
    :
    LayoutIndex(a_promotion)
  {}
  DataIndex()
    :
    LayoutIndex()
  {}

  virtual ~DataIndex()
  {}

private:
  friend class LayoutIterator;
  friend class DataIterator;
  friend class TimedDataIterator;
  friend class NeighborIterator;
  friend class BoxLayout;
  friend class Copier;

  DataIndex(int        a_indexIntoBox,
            int        a_indexIntoData,
            const int* a_layoutID)
    :
    LayoutIndex(a_indexIntoBox, a_indexIntoData, a_layoutID)
  {
  }
};

#ifndef WRAPPER
inline LayoutIndex::LayoutIndex()
  :
  m_index(0),
  m_layoutIntPtr(NULL)
{
}

inline bool LayoutIndex::operator == (const LayoutIndex& a_rhs) const
{
  return (m_index == a_rhs.m_index);
}

inline int LayoutIndex::intCode() const
{
  return m_index;
}

inline bool LayoutIndex::isNull() const
{
  return m_layoutIntPtr == NULL;
}
#endif

#include "NamespaceFooter.H"
#endif
