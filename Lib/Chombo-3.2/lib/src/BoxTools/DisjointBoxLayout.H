#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DISJOINTBOXLAYOUT_H_
#define _DISJOINTBOXLAYOUT_H_

#include "Vector.H"
#include "BoxLayout.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

struct SliceSpec;

///A BoxLayout that has a concept of disjointedness
/**
   DisjointBoxLayout is-a BoxLayout that has a concept of disjointedness
   of the boxes.  DisjointBoxLayout is only different from BoxLayout in
   that it has a method isDisjoint().  DisjointBoxLayouts also have
   a ProblemDomain member to ensure disjointness in the case of periodic
   boundary conditions, among other things.
*/

class DisjointBoxLayout: public BoxLayout
{

public:

  friend class Copier;

  /**
     \name Constructors, destructors, defines
  */
  /*@{*/

  ///
  /**
     Construct BoxLayout with no boxes.
   */
  DisjointBoxLayout();

  ///
  /** Construct from a Vector of Boxes and optionally a Vector of
      processor assignments.  If the processor assignment Vector
      is present, it must be either zero-length or the same length
      as the Box Vector.  On exit, the DisjointBoxLayout will be closed.
  */
  DisjointBoxLayout(const Vector<Box>& a_boxes,
                    const Vector<int>& a_procIDs);

  ///
  /** Construct from a Vector of Boxes and optionally a Vector of
      processor assignments.  If the processor assignment Vector
      is present, it must be either zero-length or the same length
      as the Box Vector.  On exit, the DisjointBoxLayout will be closed.
  */
  DisjointBoxLayout(const Vector<Box>& a_boxes,
                    const Vector<int>& a_procIDs,
                    const ProblemDomain& a_physDomain);

  ///
  /**
     Ref-counted destruction.  Once the last reference to the
     implementation of this class is destroyed, the data members
     are cleaned up.
  */
  virtual
  ~DisjointBoxLayout()
  {}

  ///
  /** Define this DisjointBoxLayout from a Vector of Boxes.  This is a
      wrapper for a call to BoxLayout::define, required because
      DisjointBoxLayout::define is overloaded and the compiler will not
      look into the base class for name resolution.
  */
  virtual void
  define(const Vector<Box>& a_boxes,
         const Vector<int>& a_procIDs);

  ///
  /**
    Like the define() that takes the same two arguments, except this one
    also initializes a_procIDs if it's not NULL.  In that case, the passed-in
    a_procIDs should point to a Vector of size 0.
  */
  void
  defineAndLoadBalance(const Vector<Box>& a_boxes,
                       Vector<int>* a_procIDs);

  ///
  /**
    Like the define() that takes the same two arguments, except this one
    also initializes a_procIDs if it's not NULL.  In that case, the passed-in
    a_procIDs should point to a Vector of size 0.
  */
  void
  defineAndLoadBalance(const Vector<Box>& a_boxes,
                       Vector<int>* a_procIDs,
                       const ProblemDomain& a_physDomain);

  ///
  /** Shallow define. Only way to promote a BoxLayout.  If BoxLayout
      has been closed, then this method checks isDisjoint and throws an
      error if not disjoint.  If <i>a_layout</i> is disjoint, then this object
      becomes a closed DisjointBoxLayout object.
  */
  void
  define(const BoxLayout& a_layout);

  ///
  /** Shallow define. Only way to promote a BoxLayout.  If BoxLayout
      has been closed, then this method checks isDisjoint and throws an
      error if not disjoint.  If <i>a_layout</i> is disjoint, then this object
      becomes a closed DisjointBoxLayout object.  If BoxLayout has been closed,
      and is sorted, the neighbour iterator is rebuilt (now with knowledge of
      the problem domain).
  */
  void
  define(const BoxLayout& a_layout, const ProblemDomain& a_physDomain);

  ///
  /** Define this DisjointBoxLayout from a Vector of Boxes.  This is a
      wrapper for a call to BoxLayout::define, required because
      DisjointBoxLayout::define is overloaded and the compiler will not
      look into the base class for name resolution.
  */
  virtual void
  define(const Vector<Box>& a_boxes,
         const Vector<int>& a_procIDs,
         const ProblemDomain& a_physDomain);

  ///
  /** Define this DisjointBoxLayout from a Vector of Boxes.
      Same as define.
  */
  virtual void
  define_pd(const Vector<Box>& a_boxes,
            const Vector<int>& a_procIDs,
            const ProblemDomain& a_physDomain)
  {define(a_boxes, a_procIDs,  a_physDomain);}

  /*@}*/

  /**
     \name Checks
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object contains a disjoint
      union of Boxes.

      The Disjoint testing algorithm assumes that the
      boxes are CELL-Centered boxes.  We can expand the algorithm
      when someone needs it.  bvs
  */
  bool
  isDisjoint() const;

  ///
  /** Checks to see that problem domains are compatible.
      To be compatible:
      - <i>a_domain</i> and the ProblemDomain in
      this DisjointBoxLayout must have the same periodicity.
      - The ProblemDomains must be periodic in the same dimensions.
      - In all periodic directions, the period must be the same.

      There are no checks in non-periodic directions, since
      DisjointBoxLayouts don't care about physical domains
      except in the periodic case.
   */
  bool checkPeriodic(const ProblemDomain& a_domain) const;

  /** Performs problem domain check when we don't have direct
      access to a problemDomain.

      Returns <tt>true</tt> if this and <i>a_dbl</i>
      have compatible ProblemDomains.
  */
  bool checkDomains(const DisjointBoxLayout& a_dbl) const;

  /*@}*/

  /**
     \name Modification functions
  */
  /*@{*/

  ///
  /** Mark this DisjointBoxLayout as complete and unchangeable.
   */
  virtual void
  close();

  ///
  /** Set ProblemDomain.  Useful when using addBox() and trying to use closeNoSort()
   */
  void setDomain(const ProblemDomain& a_domain)
  {
    if (*m_closed)
      {
        MayDay::Error("attempt to setDomain on closed DisjointBoxLayout");
      }
    m_physDomain = a_domain;
  }

  ///
  /** Create new open DisjointBoxLayout with copied data.

      Doesn't check for disjointness, since <i>a_source</i> is already
      guaranteed to be disjoint.  Also copies problem domain
      info. */
  virtual void
  deepCopy(const DisjointBoxLayout& a_source);

  ///
  /** Create new open DisjointBoxLayout with copied data.

      Checks that the source BoxLayout isDisjoint, throws an error
      if it is not. Otherwise, same as BoxLayout::deepCopy. */
  virtual void
  deepCopy(const BoxLayout& a_source);

  ///
  /** Create new open DisjointBoxLayout with copied data.

      Checks that the source BoxLayout isDisjoint, throws an error
      if it is not. Otherwise, same as BoxLayout::deepCopy. */
  virtual void
  deepCopy(const BoxLayout& a_source, const ProblemDomain& a_physDomain);

  /**
    Makes a_to a one-cell-thick version of *this -- sliced at a_ss.position
    along the a_ss.direction'th axis.
    If a_ss.position is outside the range of *this, that's a fatal error.
  */
  void degenerate( DisjointBoxLayout& a_to, const SliceSpec& a_ss ) const;

  ///
  /**

     Coarsen a DisjointBoxLayout:
     - <i>output</i> must be open
     - <i>input</i> must be closed
     - <i>refinement</i> must be a positive non-zero integer
     - <i>output</i> and <i>input</i> do not share an
     implementation.

     <i>output</i> is first deepCopy'ed from <i>input</i>,
     then coarsen(refinement) is called on each box of <i>output</i>.

     <i>output</i> returns from this function closed.

     LayoutIterators and DataIterators from <i>input</i> and <i>output</i>
     can be used interchangeably.
  */
  friend void coarsen(DisjointBoxLayout& output,
                      const DisjointBoxLayout& input,
                      int refinement);

  ///
  /**
     Refine a DisjointBoxLayout:
     - <i>output</i> must be open
     - <i>input</i> must be closed
     - <i>refinement</i> must be a positive non-zero integer
     - <i>output</i> and <i>input</i> do not share an
     implementation.

     <i>output</i> is first deepCopy'ed from <i>input</i>,
     then refine(refinement) is called on each box of <i>output</i>.

     <i>output</i> returns from this function closed.

     LayoutIterators and DataIterators from <i>input</i> and <i>output</i>
     can be used interchangeably.
  */
  friend void refine(DisjointBoxLayout& output,
                     const DisjointBoxLayout& input,
                     int refinement);

  /*@}*/

  /**
     \name Neighbor functions
  */
  /*@{*/

  ///
  /** Returns the set of boxes which you get from calling Box::adjCellLo
      on each box in input. Output must be open, is closed upon exit from
      this function. A negative length results in the set of boxes
      of width a_len on the low side of, but _inside_ the original boxes.
  */
  friend void adjCellLo(DisjointBoxLayout& a_output,
                        const DisjointBoxLayout& a_input,
                        int a_dir, int a_len);

  friend void adjCellLo(DisjointBoxLayout& a_output,
                        const DisjointBoxLayout& a_input,
                        int a_dir) { adjCellLo(a_output,a_input,a_dir,1);}
  ///
  /** Returns the set of boxes which result from calling Box::adjCellHi
      on each box in input.  Output must be open, will be closed upon
      exiting this function.  A negative length results in the set of boxes
      of width a_len on the low side of, but _inside_ the original boxes.
  */
  friend void adjCellHi(DisjointBoxLayout& a_output,
                        const DisjointBoxLayout& a_input,
                        int a_dir, int a_len);

  friend void adjCellHi(DisjointBoxLayout& a_output,
                        const DisjointBoxLayout& a_input,
                        int a_dir) { adjCellHi(a_output,a_input,a_dir,1);}

  /*@}*/

  const ProblemDomain& physDomain() const;

  virtual void closeNoSort(); // close without sorting; used by AdjCellHi, etc.

protected:
  ///
  friend class NeighborIterator;
  ProblemDomain m_physDomain;

  // The "pair" here : the first int is the index of the periodic shift iterator, the second unsigned int is
  // the box index.  an unshifted box is stored with a -1 shift index.
  RefCountedPtr<Vector<Vector<std::pair<int, LayoutIndex> > > > m_neighbors; //what a mouthful

  void computeNeighbors();
  virtual void closeN( RefCountedPtr<Vector<Vector<std::pair<int, LayoutIndex> > > > neighbors);
  virtual void closeNO(); // close, but don't bother building neighbors

};
void adjCellLo_dbl(DisjointBoxLayout& a_output,
                   const DisjointBoxLayout& a_input,
                   int a_dir, int a_len);
void adjCellHi_dbl(DisjointBoxLayout& a_output,
                   const DisjointBoxLayout& a_input,
                   int a_dir, int a_len);
void coarsen_dbl(DisjointBoxLayout& output,
                 const DisjointBoxLayout& input,
                 int refinement);
void refine_dbl(DisjointBoxLayout& output,
                const DisjointBoxLayout& input,
                int refinement);

inline
void adjCellLo_dbl(DisjointBoxLayout& a_output,
                   const DisjointBoxLayout& a_input,
                   int a_dir, int a_len)
{
  adjCellLo(a_output,a_input,a_dir,a_len);
}
inline
void adjCellHi_dbl(DisjointBoxLayout& a_output,
                   const DisjointBoxLayout& a_input,
                   int a_dir, int a_len)
{
  adjCellHi(a_output,a_input,a_dir,a_len);
}
inline
void coarsen_dbl(DisjointBoxLayout& output,
                 const DisjointBoxLayout& input,
                 int refinement)
{
  coarsen(output, input, refinement);
}
inline
void refine_dbl(DisjointBoxLayout& output,
                const DisjointBoxLayout& input,
                int refinement)
{
  refine(output, input, refinement);
}

#include "NamespaceFooter.H"
#endif
