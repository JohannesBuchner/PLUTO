#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INTVECTSET_H_
#define _INTVECTSET_H_

#include "SPACE.H"

#ifndef WRAPPER
#include <iostream>
#include "Box.H"
#include "IntVect.H"
#include "TreeIntVectSet.H"
#include "DenseIntVectSet.H"
#include "parstream.H"
#include "NamespaceHeader.H"

#endif

/// An irregular domain on an integer lattice
/**
   \nosubgrouping
   IntVectSet represents an irregular region in an integer lattice
   SpaceDim-dimensional index space as an arbitrary collection of
   IntVects.  A full set calculus is defined.  Any IntVect or
   cell-centered Box can be fully represented as an IntVectSet.  There
   is an iterator that provides access to the contents of an
   IntVectSet.  IntVectSets are implicitly cell-centered.
   Intersection, union, and complement operations are defined for
   pairs of IntVectSets (and, by extension, an IntVectSet and an
   IntVect or an IntVectSet and a Box).  The minimum Box of an
   IntVectSet is defined as the smallest Box that contains every
   IntVect in the IntVectSet.

   The IntVects in an IntVectSet do not have a canonical ordering.
*/
class IntVectSet
{
public:
  friend
  class IVSIterator;

  /**
        \name Constructors, Destructor, and defines
  */

  /*@{*/

  ~IntVectSet();

  /// default constructor: defines an empty IntVectSet.
  IntVectSet();
  /// default define: modifies this IntVectSet to be empty.
  void
  define();

  /// copy constructor
  IntVectSet(const IntVectSet& ivs);
  /// copy define: modifies this IntVectSet into a copy of \a ivs
  void
  define(const IntVectSet& ivs);
  /// copy define: same as \c define(IntVectSet)
  void
  define_intvectset(const IntVectSet& ivs)
  {
    define(ivs);
  }
  /// return a copy of this IntVectSet
  IntVectSet
  copy() const
  {
    return *this;
  }

  /// conversion constructor
  explicit
  IntVectSet(const DenseIntVectSet& a_dense);
  /// conversion define
  void
  define (const DenseIntVectSet& a_dense);
  /// conversion constructor
  explicit
  IntVectSet(const TreeIntVectSet& a_tree);
  /// conversion define
  void
  define (const TreeIntVectSet& a_tree);

  /// IntVect constructor
  /** construct this to be an IntVectSet with just one IntVect. */
  explicit
  IntVectSet(const IntVect& iv);
  /// IntVect define: modifies this IntVectSet to have just one IntVect
  void
  define(const IntVect& iv);
  /// IntVect define: same as \c define(IntVect)
  void
  define_intvect(const IntVect& iv)
  {
    define(iv);
  }

  /// Box constructor
  /** construct this to be an IntVectSet with all the IntVects in the Box \b. */
  explicit
  IntVectSet(const Box& b);
  /// Box define
  /** modifies this IntVect to have all the IntVects in the Box \b. */
  void
  define(const Box& b);
  /// Box define: same as \c define(Box)
  void
  define_box(const Box& b)
  {
    define(b);
  }
  /// Define this IntVectSet to have all the corners of the input box.
  void
  define_boxCorners(const Box& b);

  /*@}*/

  /**
    \name Union operators
  */
  /*@{*/

  /// unions another IntVectSet \a ivs into this IntVectSet
  IntVectSet&
  operator|=(const IntVectSet& ivs);

  /// same as operator|= except it doesn't return *this
  void
  or_intvectset(const IntVectSet& ivs)
  {
    *this |= ivs;
  }

  /// unions a single IntVect \a iv into this IntVectSet
  IntVectSet&
  operator|=(const IntVect& iv);
  /// same as operator|= except it doesn't return *this
  void
  or_intvect(const IntVect& ivs)
  {
    *this |= ivs;
  }

  /// unions the IntVects in the Box \b into this IntVectSet
  IntVectSet&
  operator|=(const Box& b);
  /// same as operator|= except it doesn't return *this
  void
  or_box(const Box& b)
  {
    *this |= b;
  }

  /// Returns a new IntVectSet that is the union of two IntVectSets
  friend
  IntVectSet
  operator|(const IntVectSet& ivs1, const IntVectSet& ivs2);

  /// Returns a new IntVectSet that is the union of an IntVectSet and an IntVect
  friend
  IntVectSet
  operator|(const IntVectSet& ivs, const IntVect& iv);

  /// Returns a new IntVectSet that is the union of an IntVectSet and an IntVect
  friend
  IntVectSet
  operator|(const IntVect& iv, const IntVectSet& ivs);

  /// Returns a new IntVectSet that is the union of an IntVectSet and a Box
  friend
  IntVectSet
  operator|(const IntVectSet& ivs, const Box& b);

  /// Returns a new IntVectSet that is the union of an IntVectSet and a Box
  friend
  IntVectSet
  operator|(const Box& b, const IntVectSet& ivs);

  /*@}*/

  /**
         \name Complement operators
  */
 /*@{*/

  /// Returns the complement of the IntVectSet \a ivs within this IntVectSet
  IntVectSet
  operator-(const IntVectSet& ivs) const;

  /// Returns the complement of the Box \a b within this IntVectSet
  IntVectSet
  operator-(const Box& b) const;

  /// Returns the complement of the IntVect \a iv within this IntVectSet
  IntVectSet
  operator-(const IntVect& iv) const;

  /// Remove the IntVects in the IntVectSet \a ivs from this IntVectSet
  /**
   * Modifies this IntVectSet to be the complement of the IntVectSet \a ivs
   * within this IntVectSet.
   */
  IntVectSet&
  operator-=(const IntVectSet& ivs);
  /// same as operator-= except it doesn't return *this
  void
  minus(const IntVectSet& ivs)
  {
    *this -= ivs;
  }

  /// Remove the IntVects in the Box \a b from this IntVectSet
  /**
   * Modifies this IntVectSet to be the complement of the argument
   *  Box within this IntVectSet.
   */
  IntVectSet&
  operator-=(const Box& b);
  /// same as operator-= except it doesn't return *this
  void
  minus_box(const Box& b)
  {
    *this -=b;
  }

  /// Remove the IntVect \a iv from this IntVectSet
  /**
     Modifies this IntVectSet to be the complement of the argument
     IntVect within this IntVectSet.
  */
  IntVectSet&
  operator-=(const IntVect& iv);
  /// same as operator-= except it doesn't return *this
  void
  minus_intvect(const IntVect& iv)
  {
    *this -=iv;
  }

  /*@}*/

  /**
     \name Intersection operators */
  /*@{*/

  /// Returns a new IntVectSet that is the intersection of two IntVectSets
  /**
   *  The result may be empty.
   */
  friend
  IntVectSet
  operator&(const IntVectSet& ivs1, const IntVectSet& ivs2);

  /// Returns a new IntVectSet that is the intersection of an IntVectSet and a Box
  /**
   *  The result may be empty.
   */
  friend
  IntVectSet
  operator&(const IntVectSet& ivs, const Box& b);

  /// Returns a new IntVectSet that is the intersection of an IntVectSet and a Box
  /**
   *  The result may be empty.
   */
  friend
  IntVectSet
  operator&(const Box& b, const IntVectSet& ivs);

  /// Modifies this IntVectSet to its intersection with another IntVectSet
  IntVectSet&
  operator&=(const IntVectSet& ivs);
  /// same as operator&= except it doesn't return *this
  void
  and_intvectset(const IntVectSet& ivs)
  {
    *this &= ivs;
  }

  /// Modifies this IntVectSet to be its intersection with a Box
  IntVectSet&
  operator&=(const Box& b);
  /// same as operator&= except it doesn't return *this
  void
  and_box(const Box& b)
  {
    *this &= b;
  }

  /// Modifies this IntVectSet to be its intersection with the Box in a ProblemDomain
  IntVectSet&
  operator&=(const ProblemDomain& domain);
  /// same as operator&= except it doesn't return *this
  void
  and_domain(const ProblemDomain& d)
  {
    *this &= d;
  }

  /*@}*/

  /**
     \name Modification functions */
  /*@{*/

  /// Add IntVects to this IntVectSet in all directions
  /**
     Modifies this IntVectSet by growing it by the specified number
     of IntVects in all directions, including diagonal.  It is an error
     if \a igrow < 0.
  */
  void
  grow(int igrow);

  /// Add IntVects to an IntVectSet in all directions
  /**
     Creates a new IntVectSet that is a copy of the argument IntVectSet \a ivs
     grown by the specified number of IntVects in all directions \a igrow, including
     diagonal.  It is an error if \a igrow < 0.
  */
  friend
  IntVectSet
  grow(const IntVectSet& ivs, int igrow);

  /// Add IntVects to this IntVectSet in one direction
  /**
     Modifies this IntVectSet by growing it by the specified number of
     IntVects \a igrow in the specified coordinate direction \a idir.
     Directions are zero-based.
  */
  IntVectSet&
  grow(int idir, int igrow);
  /// same as \c grow(idir,igrow) except doesn't return *this
  void
  grow_dir(int idir, int igrow)
  {
    grow(idir, igrow);
  }

  /// Analogous to surroundingNodes() for a Box.
  /**
     Modifies this IntVectSet by adding IntVects to the high-side of each
     direction.
  */
  void
  growHi();

  /// Analogous to surroundingNodes(dir) for a Box.
  /**
     Modifies this IntVectSet by adding IntVects to the high-side of a specified
     direction.
  */
  void
  growHi(const int a_dir);

  /// Refine all the IntVects in this IntVecSet
  /**
     Modifies this IntVectSet by refining it by the factor \a iref.  It is
     an error if \a iref <= 0.  Definition of refinement: for each IntVect \p iv
     in the original IntVectSet, the refined IntVectSet will contain the Box
     defined by <b>refine( Box(iv,iv), iref )</b>.
   */
  IntVectSet&
  refine(int iref = 2);

  /// Refine all the IntVects in an IntVectSet
  /**
     Creates a new IntVectSet that is a copy of the argument IntVectSet \a ivs
     refined by the factor \a iref.  It is an error if \a iref <= 0.
     Definition of refinement: for each IntVect \p iv in the original
     IntVectSet \a ivs, the refined IntVectSet will contain the Box defined by
     <b>refine( Box(iv,iv), iref )<b>.
   */
  friend
  IntVectSet
  refine(const IntVectSet& ivs, int iref);

  /// Coarsen all the IntVects in this IntVectSet
  /**
     Modifies this IntVectSet by coarsening it by the factor \a iref.  It is
     an error if \a iref <= 0.  Definition of coarsening: for each IntVect
     \p iv in the original IntVectSet, the refined IntVectSet will contain the
     IntVect defined by <b>coarsen( iv, iref )</b>.
   */
  IntVectSet&
  coarsen(int iref = 2);

  /// Coarsen all the IntVects in an IntVectSet
  /**
     Creates a new IntVectSet that is a copy of the argument IntVectSet \a ivs coarsened
     by the factor \a iref.  It is an error if \a iref <= 0.
     Definition of coarsening: for each IntVect \p iv in the original
     IntVectSet, the refined IntVectSet will contain the IntVect defined by
     <b>coarsen( iv, iref )</b>.
   */
  friend
  IntVectSet
  coarsen(const IntVectSet& ivs, int iref);

  /// Increment all the IntVects in this IntVectSet by the IntVect \a iv
  void
  shift(const IntVect& iv);

  /// Make this IntVectSet be properly nested
  /**
     Modifies this IntVectSet to remove all IntVects that are not at least
     \a radius from the the edge of the IntVectSet in any direction.  \a
     radius must be positive.  IntVects that border the Box \a domain are
     spared from this trimming.
  */
  void
  nestingRegion(int radius, const Box& domain, int granularity  = 1);

  /// Make this IntVectSet be properly nested
  /**
     Modifies this IntVectSet to remove all IntVects that are not at least
     \a radius from the the edge of the IntVectSet in any direction.  \a
     radius must be positive.  IntVects that border non-periodic boundaries
     are spared from this trimming.  Radius extends across periodic images.
  */
  void
  nestingRegion(int radius, const ProblemDomain& probdomain, int granularity = 1);
  /// same as \c nestingRegion() except it doesn't return *this
  void
  nestingRegion_prob(int radius, const ProblemDomain& probdomain)
  {
    nestingRegion(radius, probdomain);
  }

  /// Modifies this IntVectSet to be empty
  void
  makeEmpty();

  /// Modifies this IntVectSet to empty but leaves dense domain box unchanged.
  /**
     This leaves unchanged the domain box for a dense IntVectSet - it only sets
     all bits to zero.  So you can still do unions later (in the domain box) and
     not be converted to a tree representation.  There is no difference from
     makeEmpty() for a tree IntVectSet
  */
  void
  makeEmptyBits();

  /// Chop the IntVectSet at the \a chop_pnt index in the \a dir direction.
  /**
     Returns one IntVectSet and modifies this IntVectSet.  The union of the
     two is the original IntVectSet.  This IntVectSet gets the IntVects with
     \a dir index >= \a chop_pnt, and the returned IntVectSet gets the
     IntVects with \a dir index < \a chop_pnt.  It is an error if \a dir is
     invalid.
  */
  IntVectSet chop(int dir, int chop_pnt);

  void chop(int dir, int chop_pnt, IntVectSet& a_hi);
  ///
  /**
     Set the max box size for IntVectSet::define(Box) which will
     make the IntVectSet dense. Default is 6.4^6.
   */
  static void setMaxDense(const int& a_maxDense);

  /*@}*/

  /**
     \name Accessor and Inquiry functions */
  /*@{*/

  /// Returns the number of IntVects in this IntVectSet
  int
  numPts() const;

  /// Returns the minimum enclosing box of this IntVectSet
  const Box&
  minBox() const;

  /// Forces recalculation of the minimum enclosing box of this IntVectSet
  void
  recalcMinBox() const;

  /// Returns true if no IntVects are in this IntVectSet
  bool
  isEmpty() const;

  /// Returns true if this IntVectSet is currently being represented in a dense fashion
  bool
  isDense() const;

  /// Returns true if this IntVectSet contains \a iv
  bool
  contains(const IntVect& iv) const;

  /// Returns true if this IntVectSet contains all the IntVects in \a ivs
  bool
  contains(const IntVectSet& ivs) const;
  /// same as \c contains(ivs)
  bool
  contains_intvectset(const IntVectSet& ivs) const
  {
    return contains(ivs);
  }

  /// Returns true if this IntVectSet contains all the IntVects in \a box
  /**
   * Note: this is done using an algorithm that is much faster than looping
   * through all the IntVects.
   */
  bool
  contains(const Box& box) const;
  /// same as \c contains(b)
  bool
  contains_box(const Box& b) const
  {
    return contains(b);
  }

  /// Returns a Vector<Box> representation of this IntVectSet.
  /**
   * The returned Boxes contain only the IntVects in this IntVectSet, so
   * the union of the Boxes in the Vector == this IntVectSet.
   */
  Vector<Box>
  boxes() const;

  /// Optimize memory usage of this IntVectSet
  /**
     Attempts to optimize storage for IntVectSet.  This will also
     optimize iteration through the IntVectSet later.  Best to call
     it when you are no longer modifying the IntVectSet
  */
  void
  compact() const;

  /*@}*/

  /**
     \name Linearization functions */
  /*@{*/

  /// Returns the number of bytes in a linear representation of *this
  int
  linearSize() const;

  /// Modify *this using the data in the linear representation in \a a_inBuf
  void
  linearIn(const void* const a_inBuf);

  /// Write a linear representation of *this to \a a_outBuf
  /**
     Assumes a_outBuf is at least linearSize() bytes long.
   */
  void
  linearOut(void* const a_outBuf) const;

  /*@}*/

  /// Returns true if this IntVectSet has the same IntVects as \a a_ivs
  bool
  operator==(const IntVectSet& a_ivs) const;

  /** Primary sorting criterion: being dense; if *this is dense and \a a_ivs is
      not, then *this is smaller, and vice versa.
      Secondary sorting criterion: operator< as defined on DenseIntVectSet or
      TreeIntVectSet.
      In a total tie, returns false.

      These criteria might not seem natural, but that doesn't matter as the only
      reason for this operator is to support using IntVectSet as the key to an
      std::map.
  */
  bool
  operator<(const IntVectSet& a_ivs) const;

  /// Writes a Vector<Box> representation to an output stream.
  void
  printBoxes(std::ostream& a_ostream) const;

  /// Writes a Vector<Box> representation to the parallel stdout.
  void
  p() const
  {
    printBoxes(pout());
  }

  /// Writes a text representation of an IntVectSet to an output stream
  friend
  std::ostream&
  operator<<(std::ostream& os, const IntVectSet& ivs);

  void convert() const; // turn dense rep into Tree rep.  very costly.
                        // it is 'logically' const, but does modify data structures;

  // not for public consumption.  used in memory tracking.
  static long int count;
  static long int peakcount;

  //set to 6400000 as default.  resettable.
  static int s_maxDense;

private:

  bool m_isdense;
  TreeIntVectSet m_ivs;
  DenseIntVectSet m_dense;
  // not a user function.  called by memory tracking system on
  // exit to clean up static allocation pools used for the optimization
  // of these routines.
  friend void dumpmemoryatexit();
  static void clearStaticMemory(); // after this functon is called, you cannot
                              // do any other operations with IntVectSet.
};

/// Iterator for an IntVectSet
/**
   IVSIterator iterates over every point (IntVect)
   in an IntVectSet.  It has exactly the same
   syntax and sematic as BoxIterator.  Typical usage:

\code
   IntVectSet ivs;
   ...
   IVSIterator ivsit (ivs);
   for (ivsit.begin(); ivsit.ok(); ++ivsit)
   {
     IntVect iv = ivsit();
     (do operations involving iv)
   }
\endcode
*/
class IVSIterator
{
public:

  /**
   * A default constructed iterator iterates over an empty IntVectSet.
   * It starts in the \c begin() state, and is never \c ok().
   */
  IVSIterator():m_isdense(true)
  {}

  /**
   * Iterates over IntVectSet \a ivs.  There is no canonical ordering.
   * Starts in the \c begin() state.
   */
  IVSIterator(const IntVectSet& ivs);

  ~IVSIterator()
  {}

  /**
   * Modifies this IVSIterator to iterate over the IntVectSet \a ivs.
   * There is no canonical ordering.  Starts in the \c begin() state.
   */
  void define(const IntVectSet& ivs);

  /// returns the IntVect that this iterator is at
  const IntVect& operator()() const ;

  /// same as operator()
  const IntVect& iv() const
  {
    return this->operator()();
  }

  /// returns true if this iterator is still in its IntVectSet
  bool ok() const;

  /// move iterator to the next IntVect in its IntVectSet
  void operator++();

  /// same as operator++
  void incr()
  {
    ++(*this);
  }

  /// initialize this iterator to the first IntVect in its IntVectSet
  void begin();

  /// same as begin()
  void reset();

  /// move this iterator to after the last IntVect in the set
  /** The iterator will be !ok() afterwards. */
  void end();

private:
  bool m_isdense;
  DenseIntVectSetIterator m_dense;
  TreeIntVectSetIterator  m_tree;
};

#ifndef WRAPPER

inline const IntVect& IVSIterator::operator()() const
{
  if (m_isdense) return m_dense();
  return m_tree();
}

inline bool IVSIterator::ok() const
{
  if (m_isdense) return m_dense.ok();
  return m_tree.ok();
}

inline void  IVSIterator::operator++()
{
  if (m_isdense) ++m_dense;
  else          ++m_tree;
}
inline void IVSIterator::reset()
{
  begin();
}

inline void IVSIterator::begin()
{
  if (m_isdense) m_dense.begin();
  else          m_tree.begin();
}

inline void IVSIterator::end()
{
  if (m_isdense) m_dense.end();
  else          m_tree.end();
}

inline IntVectSet::IntVectSet(): m_isdense(true)
{
  count++;
  if (count > peakcount) peakcount = count;
}

inline void IntVectSet::define(const IntVectSet& ige_in)
{
  *this = ige_in;
}

inline IntVectSet::IntVectSet(const IntVectSet& ige_in)
{
  count++;
  if (count > peakcount) peakcount = count;

  *this = ige_in;
}

inline   bool
IntVectSet::isDense() const
{
  return m_isdense;
}

#endif /* WRAPPER */

#include "NamespaceFooter.H"
#endif
