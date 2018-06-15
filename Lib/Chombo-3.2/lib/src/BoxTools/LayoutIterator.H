#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LAYOUTITERATOR_H_
#define _LAYOUTITERATOR_H_

#include "DataIndex.H"
#include "BoxLayout.H"

#include "NamespaceHeader.H"

///An Iterator based on a BoxLayout object.
/**
   An Iterator based on a BoxLayout object.  It does not
   support a dereferencing operation(1), since it is intended
   to work with all of BoxLayouts, DisjointBoxLayouts, BoxLayoutDatas
   LevelData's, and any object that is built on top
   of a BoxLayout object.  LayoutIterator accesses the data in a BoxLayout-based
   object in a NON-data-parallel manner (i.e. every processor iterates through
   \b all the Boxes in the BoxLayout).  This differs from the DataIterator class.

   BoxLayout-based objects can act as the Factory for the LayoutIterator.

 (1) STL-speak.  not critical for comprehension, but can help people
 familiar with STL iterators and expecting similar behaviour.

*/
class LayoutIterator
{
public:
  /// a null constructed LayoutIterator will return false on ok()
  LayoutIterator()
  {}

  // Default copy and null constructor should be fine.
  // Not useful to someone using the iterator.  Used by
  // other classes in Chombo.  Could make it private and
  // let BoxLayout have access....

  virtual ~LayoutIterator()
  {}

  /// return the index that this iterator is at
  const LayoutIndex& operator () () const;

  /// return a copy of the index that this iterator is at
  LayoutIndex i() const
  {
    return this->operator()();
  }

  /// move the iterator to the next Box in the layout
  virtual void operator ++ ();

  /// move the iterator to the next Box in the layout
  void incr()
  {
    ++(*this);
  }

  /// return true if this iterator is still in its Layout
  virtual bool ok() const;

  /// initialize this iterator to the first Box in its Layout
  void begin();

  /// same as begin()
  void reset();

  /// move this iterator to after the last Box in the layout
  /** The iterator will be !ok() afterwards. */
  void end();



  const LayoutIndex& operator[](int ivec) const
  {
    return (*this->m_indicies)[ivec];
  }

  LayoutIndex& operator[](int ivec)
  {
    return (*this->m_indicies)[ivec];
  }

protected:
  friend class BoxLayout;
  friend class DisjointBoxLayout;
  friend class TimedDataIterator;
  LayoutIterator(const BoxLayout& a_boxlayout,
                 const int*       a_layoutID);

  BoxLayout    m_layout;
  RefCountedPtr<Vector<LayoutIndex> > m_indicies;

  int m_current;
};

inline const LayoutIndex& LayoutIterator::operator () () const
{
  CH_assert(ok());
  return (*m_indicies)[m_current];
}

inline void LayoutIterator::operator ++ ()
{
  ++m_current;
}

inline bool LayoutIterator::ok() const
{
  return m_current < m_layout.size();
}

inline void LayoutIterator::begin()
{
  m_current = 0;
}

#include "NamespaceFooter.H"
#endif
