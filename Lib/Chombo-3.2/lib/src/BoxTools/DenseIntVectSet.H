#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DENSEINTVECTSET_H_
#define _DENSEINTVECTSET_H_

#include "BitSet.H"
#include "IntVect.H"
#include "Box.H"
#include "Vector.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

class ProblemDomain;
class DenseIntVectSetIterator;
/// Dense representation implementation of IntVectSet class

/**  Performance-oriented class that reproduces most of the
  functionality of IntVectSet.  This object stores a dense
  representation of an IntVectSet over the box of definition

  Storage is performed with a BitSet.  This stores 0 or 1 on
  a per-bit basis for each location in the IntVectSet.

  For an explanation of undocumented functions look at IntVectSet

@see IntVectSet

*/
class DenseIntVectSet
{
public:
  ///
  DenseIntVectSet()
  {}

  ///
  /** you can either have the domain begin all true, or all false */
  DenseIntVectSet(const Box& a_domain, bool init = true);
  // copy, and operator= should be fine

  ///
  /**
     make this set the complement of itself and the input
  */
  DenseIntVectSet& operator-=(const Box& box);

  ///
  /**
     make this set the complement of itself and the input
  */
  DenseIntVectSet& operator-=(const IntVect& intvect);

  ///
  /**
     make this set the complement of itself and the input
  */
  DenseIntVectSet& operator-=(const DenseIntVectSet& ivs);

  ///
  /**
     make this set the intersection of itself and the input
  */
  DenseIntVectSet& operator&=(const Box& box);

  ///
  /**
     make this set the intersection of itself and the input
  */
  DenseIntVectSet& operator&=(const ProblemDomain& domain);

  ///
  /**
     make this set the union of itself and the input.  intvect MUST
         be within the domain of of the DenseIntVectSet or this will throw
         an error
  */
  DenseIntVectSet& operator|=(const IntVect& intvect);

  ///
  /**
     make this set the union of itself and the input.  b MUST
         be within the domain of of the DenseIntVectSet or this will throw
         an error
  */
  DenseIntVectSet& operator|=(const Box& b);

  ///
  /**
     make this set the union of itself and the input. resulting
         DenseIntVectSet has domain of the minBox holding both oeprator
         and operand
  */
  DenseIntVectSet& operator|=(const DenseIntVectSet& b);

  ///
  /**
     make this set the intersection of itself and the input
  */
  DenseIntVectSet& operator&=(const DenseIntVectSet& ivs);

  ///Shift every cell in the set by the input.  Fast
  /**
     Shift every cell in the set by the input.
  */
  void shift(const IntVect& iv);

  /// O(1) time inquiry of containment
  /* slower, constant time pointwise access.
     You should not build a BoxIterator, then access this
     method if you intend to iterate over the objects true
     or false terms.  you should build a DenseIntVectSetIterator */
  bool operator[](const IntVect& index) const;

  /// O(1) time inquiry of containment
  /** returns 'true' if the entire set of points specified
      by the box 'b' is a member of the IntvectSet */
  bool contains(const Box& box) const;

  ///
  /**
     coarsen the set by the input
  */
  void coarsen(int iref);

  ///
  /**
     refine the set by the input
  */
  void refine(int iref);

  ///
  /**
     grow the set by the input
  */
  void grow(int igrow);

  ///
  /**
     grow the set by the input
  */
  void grow(int idir, int igrow);

  ///
  /**
     @see IntVectSet::growHi
  */
  void growHi();

  ///
  /**
     @see IntVectSet::growHi(int)
  */
  void growHi(const int a_dir);

  ///
  /**
      Make empty by setting all bits to false
  */
  void makeEmptyBits();

  ///
  /** Chop the DenseIntVectSet  at the chop_pnt in the dir direction
      Returns one DenseIntVectSet and  modifies the object DenseIntVectSet.
      The union of the two is the original IntVectSet.
      The modified DenseIntVectSet is the low end, the returned DenseIntvectSet
      is the high end.  @see IntVectSet::chop */
  DenseIntVectSet chop(int dir, int chop_pnt);

  inline const Box& box() const;

  ///
  /** @see IntVectSet::nestingRegion
   */
  void nestingRegion(int a_radius, const Box& a_domain);

  ///

  void nestingRegion(int a_radius, const ProblemDomain& a_domain);

  ///cheaper than numPts()==0 by wide margin
  /**
     return true if the set has no points
  */
  bool isEmpty() const; //cheaper than numPts by wide margin

  ///cheaper than numPts()==box().numPts(), by significant margin
  /**
     return true if the set consists of a set of a full box
  */
  bool isFull() const;  // cheaper than numPts, bu significant margin

  ///
  /**
     Return true if all bits on the high hyperfaces are zero
  */
  bool isHighEmpty() const;

  ///
  /**
     Return true if all bits on the high hyperface are zero
  */
  bool isHighEmpty(const int a_dir) const;

  ///
  /**
     return the number of points in the set
  */
  int  numPts() const;

  ///
  bool operator==(const DenseIntVectSet& a_lhs) const;

  /**
      Primary sorting criterion: Box::operator<() applied to m_domain.
      Secondary sorting criterion: BitSet::operator<() applied to m_bits.
  */
  bool operator<(const DenseIntVectSet& a_ivs) const;

  ///
  /**
         turn DenseIntVectSet into a Vector of Boxes
  */
  Vector<Box> createBoxes() const;
  ///
  /**
     set object to its minimum representation.  changes data, but
     is a logically const operation
  */
  void compact() const ;

  void recalcMinBox() const;

  const Box& mBox() const
  {
    return m_minBox;
  }

  /** \name Linearization routines */
  /*@{*/
  int linearSize() const;

  void linearIn(const void* const inBuf);

  void linearOut(void* const a_outBuf) const;

  /*@}*/
private:

  void grow(const IntVect& iv);
  // version without checking intersect or empty
  DenseIntVectSet& intersect(const DenseIntVectSet& rhs);
  bool isHighEmpty(const int a_dirS, const int a_dirE) const;

  Box m_domain;
  BitSet m_bits;
  friend class DenseIntVectSetIterator;
  Box m_minBox;
};

/// Iterate over all the 'true' members of a DenseIntVectSet set
/** This class is used by IVSIterator to implement its iterator when IntVectSet
 *  is stored as a DenseIntVectSet.
 */
class DenseIntVectSetIterator
{
public:

  ///
  DenseIntVectSetIterator();

  ///
  DenseIntVectSetIterator(const DenseIntVectSet& ivs);

  //default null, assignment, copy and destructor should work fine.

  ///
  void define(const DenseIntVectSet& ivs);

  ///
  const IntVect& operator()() const;

  ///
  bool ok() const;

  ///
  void operator++();

  ///
  void begin();

  ///
  void end();

  static DenseIntVectSet emptyDenseIntVectSet;

private:

  ///
  void thisIntVect(const int a_linearPos);

  ///
  void setStride();

  BitSetTrueIterator     m_iterator;
  const DenseIntVectSet* m_ivsPtr;
  IntVect                m_current;
  IntVect                m_stride;
  int                    m_prevLinearPos;
};

//===================================================================

// inline functions

inline
const Box& DenseIntVectSet::box() const
{
  return m_domain;
}

inline
DenseIntVectSet& DenseIntVectSet::operator-=(const IntVect& intvect)
{
  if (m_domain.contains(intvect))
    m_bits.setFalse(m_domain.index(intvect));
  return *this;
}

inline void DenseIntVectSet::shift(const IntVect& iv)
{
  m_domain.shift(iv);
  m_minBox.shift(iv);
}

/** After default construction, ok() is false so it looks like an empty bitset
 *  but if operator++() is called, it will fail.
 */
inline
DenseIntVectSetIterator::DenseIntVectSetIterator()
  :
  m_ivsPtr(0),
//   m_current(IntVect::Zero),  // Leave uninitialized if empty bitset
  m_stride(IntVect::Unit),      // Avoid /0 if op++ called on empty bitset (if
  m_prevLinearPos(-2)           //   defined later)
{
}

/** Never really used.  Normal construction by IntVectSet is default then
 *  define().
 */
inline
DenseIntVectSetIterator::DenseIntVectSetIterator(const DenseIntVectSet& ivs)
  :
  m_iterator(ivs.m_bits),       // Finds first linear position
  m_ivsPtr(&ivs),
//   m_current(IntVect::Zero),  // Leave uninitialized if empty bitset
  m_stride(IntVect::Unit),      // Avoid /0 if op++ called on empty bitset
  m_prevLinearPos(-2)
{
  if (ok())
    {
      setStride();
      thisIntVect(m_iterator());
    }
}

inline void
DenseIntVectSetIterator::define(const DenseIntVectSet& ivs)
{
  m_iterator.define(ivs.m_bits);  // Finds first linear position
  m_ivsPtr = &ivs;
  m_prevLinearPos = -2;
  if (ok())
    {
      setStride();
      thisIntVect(m_iterator());
    }
}

inline const IntVect&
DenseIntVectSetIterator::operator()() const
{
  return m_current;
}

inline bool
DenseIntVectSetIterator::ok() const
{
  return m_iterator.ok();
}

inline void
DenseIntVectSetIterator::operator++()
{
  CH_assert(m_ivsPtr != 0);  // Catch the case where op++ is used on a
  ++m_iterator;              //   default construction
  thisIntVect(m_iterator());
}

inline void
DenseIntVectSetIterator::begin()
{
  m_iterator.begin();
  m_prevLinearPos = -2;
  if (ok())
    {
      thisIntVect(m_iterator());
    }
}

inline void
DenseIntVectSetIterator::end()
{
  m_iterator.end();
}

inline void
DenseIntVectSetIterator::setStride()
{
  D_TERM6(
    m_stride[0] = 1;,
    m_stride[1] = std::max(1, m_stride[0]*m_ivsPtr->m_domain.size(0));,
    m_stride[2] = std::max(1, m_stride[1]*m_ivsPtr->m_domain.size(1));,
    m_stride[3] = std::max(1, m_stride[2]*m_ivsPtr->m_domain.size(2));,
    m_stride[4] = std::max(1, m_stride[3]*m_ivsPtr->m_domain.size(3));,
    m_stride[5] = std::max(1, m_stride[4]*m_ivsPtr->m_domain.size(4));)
}

#include "NamespaceFooter.H"
#endif
