#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeFArrayBox.H
// adapted from dMartin/Chombo-IAMR/util/FluxBox.H
// by Dan Martin, Fri, Jan 14, 2000
// petermc, 1 Nov 2000
// 8 Jan 2001 added

#ifndef _NODEFARRAYBOX_H_
#define _NODEFARRAYBOX_H_

#include "FArrayBox.H"
#include "NamespaceHeader.H"

/// A wrapper for an FArrayBox to contain NODE-centered data.
/** This is a class to contain NODE-centered data on a box.<br>
    <b>Question</b>:  Why not just use NODE-centered FArrayBox?<br>
    <b>Answer</b>:  We want to use LevelData<NodeFArrayBox>,
    which requires that the underlying BoxLayout be a DisjointBoxLayout.
    Adjacent NODE-centered boxes intersect, but adjacent CELL-centered
    boxes do not.<br>
    The underlying DisjointBoxLayout will be CELL-centered.

    In a LevelData<NodeFArrayBox>, some of the grids may share nodes.
    The data on shared nodes must be the same.

    Example:  A four-grid layout is shown on the left.  As shown on
    the right, the nodes marked "2" are shared by two grids, and
    those marked "3" are shared by three grids.
    <pre>
    +-----+                     YYYY--+
    |     |                     |     |
    |     |                     |     |
    |     |                     |     |
    +---+--+--+---+             +---2223222---+
    |      |      |             |      2      |
    |      |      |             |      2      |
    |      |      |             |      2      |
    +---+------+------+         X---22222223------X
    |          |                |          |
    |          |                |          |
    |          |                |          |
    +----------+                +-------YYYY
    </pre>
    Additional nodes may be shared if the boxes extend to the faces
    of a domain that is periodic in one or more directions.
    In the example above, these nodes are marked "X" and "Y".<br>
    Horizontally periodic domain:  node "X" shared by two grids.<br>
    Vertically periodic domain:  nodes "Y" shared by two grids.

    One of the most useful features of LevelData<NodeFArrayBox> is the
    copyTo() function.  THE copyTo() FUNCTION MUST BE USED WITH CAUTION
    if the source and destination have different layouts.

    Consider an example where A and B are on layouts of one box
    each, and these two boxes abut:
    <pre>
    +-------+
    |       |
    |   A   |
    |       |
    =========
    |       |
    |   B   |
    |       |
    +-------+
    </pre>
    If we do A.copyTo(B), then the data on the nodes of the interface
    (marked "=") are NOT copied, because the underlying CELL-centered
    DisjointBoxLayouts of the LevelDatas do not intersect.

    So whenever we want to use src.copyTo(dest), unless we are sure
    that src and dest have the same layout, we do the following:
    - define BoxLayoutData<NodeFArrayBox> destExpanded on the BoxLayout of
    dest with each box expanded by one in each dimension;
    - do src.copyTo(destExpanded);
    - box by box, copy destExpanded to dest.
*/
class NodeFArrayBox
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.  User must subsequently call define().
   */
  NodeFArrayBox();

  ///
  /** Constructs NodeFArrayBox on CELL-centered box <i>a_bx</i>
      with <i>a_nComp</i> components using the aliased data space.
  */
  NodeFArrayBox(const Box& a_bx, int a_nComp, Real* a_alias = NULL);

  ///
  /** Constructs an aliased NodeFArrayBox.
   */
  NodeFArrayBox(const Interval& a_comps, NodeFArrayBox& a_original);

  ///
  /** Destructor.
   */
  ~NodeFArrayBox();

  ///
  /** Defines NodeFArrayBox on CELL-centered box <i>a_bx</i>
      with <i>a_nComp</i> components.  If called more than once
      on the same instance, the box and FAB will be resize()d.
  */
  void define(const Box& a_bx, int a_nComp=1);

  /**
     Change this NodeFArrayBox so it covers the Box a_bx with a_nComps
     components.  If a_alias is not NULL, it is used as the data memory
     (and is assumed to be large enough).
  */
  void resize(const Box& a_bx, int a_nComp=1, Real* a_alias=NULL);

  /*@}*/

  /**
     \name Accessors
  */
  /*@{*/

  ///
  /** Returns the CELL-centered domain where the array is defined.
   */
  const Box& box() const;

  ///
  /** Returns a modifiable reference to the NODE-centered FArrayBox
      containing the data.
  */
  FArrayBox& getFab();

  ///
  /** Returns a constant reference to the NODE-centered FArrayBox
      containing the data.
  */
  const FArrayBox& getFab() const;

  Real dotProduct(const NodeFArrayBox& a_fab2) const
  {
    return m_fab.dotProduct(a_fab2.m_fab);
  }

  Real dotProduct(const NodeFArrayBox& a_fab2, const Box& box) const
  {
    return m_fab.dotProduct(a_fab2.m_fab, box);
  }

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Modifies the data in this NodeFArrayBox by copying data
      from <i>a_src</i> into it, over the nodes that they have in common.<br>
      If this NodeFArrayBox and <i>a_src</i> have abutting grids,
      as shown:
      <pre>
      +--------+
      |        |
      |  this  |
      |        |
      ==========   <- data copied here from a_src to this
      |        |
      |  a_src |
      |        |
      +--------+
      </pre>
      then the data in this NodeFArrayBox WILL be modified along the
      edge nodes where they abut -- even though this.box() and
      <i>a_src</i>.box(), being CELL-centered, do NOT intersect.

      All components are copied.
  */
  void copy(const NodeFArrayBox& a_src);

  ///
  /** Modifies the data in this NodeFArrayBox by copying the data
      from <i>a_src</i> into it, over the intersection of:
      - the nodes surrounding the CELL-centered box <i>a_regionFrom</i>;
      - the nodes surrounding the CELL-centered box <i>a_regionTo</i>;
      - the nodes containing data in this NodeFArrayBox;
      - the nodes containing data in <i>a_src</i>.

      The components in the interval <i>a_Csrc</i> in <i>a_src</i>
      are copied to the components in the interval <i>a_Cdest</i>
      in this NodeFArrayBox.

      This function is required in order to have BoxLayoutData<NodeFArrayBox>.
  */
  void copy(const Box& a_regionFrom,
            const Interval& a_Cdest,
            const Box& a_regionTo,
            const NodeFArrayBox& a_src,
            const Interval& a_Csrc);

  /*@}*/

  /**
     \name Linearization functions
  */
  /*@{*/

  ///
  /** Returns size, in number of bytes, of a flat linear
      representation of the data in this object in the area
      defined by the nodes in CELL-centered box <i>a_R</i> and the
      component Interval <i>a_comps</i>.

      This function is required in order to have BoxLayoutData<NodeFArrayBox>.
  */
  int size(const Box& a_R, const Interval& a_comps) const;

  ///
  /** Writes into <i>a_buf</i> a linear representation of the internal data for
      the nodes surrounding CELL-centered box <i>a_R</i>, over the component
      Interval <i>a_comps</i>.

      Assumes that sufficient
      memory for the buffer has already been allocated by the caller.

      This function is required in order to have BoxLayoutData<NodeFArrayBox>.
  */
  void linearOut(void* a_buf, const Box& a_R, const Interval& a_comps) const;

  ///
  void linearIn(void* a_buf, const Box& a_R, const Interval& a_comps);

  /// These functions are required for broadcast & gather.
  void linearOut(void* a_buf) const;
  ///
  void linearIn(const void* const a_buf);
  ///
  int linearSize(void) const;

  ///
  static int preAllocatable()
  {
    return 0;
  }

  /*@}*/

  ///
  /**
     Modifies this NodeFArrayBox so that all values of all components are set to
     the given value, a_x.  (FArrayBox inherits such a function from BaseFab.)
  */
  void setVal(Real a_x);

  /**
     The setVal functions set subregions in the `BaseFab' to a constant value.
     This most general form specifies the subbox, the starting component
     number, and the number of components to be set.
  */
  void setVal(Real          a_x,
              const Box& a_bx,
              int        a_nstart,
              int        a_numcomp);

  ///
  /** HERE
      Returns the Lp-norm of this FAB using components
      (a_comp : a_comp+a_numcomp-1) and within the a_subbox.
      a_p < 0 -> ERROR
      a_p = 0 -> infinity norm (max norm)
      a_p = 1 -> sum of ABS(FAB)
      a_p > 1 -> Lp-norm
  */
  Real norm(const Box& a_subbox,
            int        a_p = 2,
            int        a_comp = 0,
            int        a_numcomp = 1) const
  {
    return m_fab.norm(a_subbox, a_p, a_comp, a_numcomp);
  }

  ///
  /**
     Returns the Lp-norm of this FAB using components
     (a_comp : a_comp+a_numcomp-1).
     a_p < 0  -> ERROR
     a_p = 0  -> infinity norm (max norm)
     a_p = 1  -> sum of ABS(FAB)
     a_p > 1  -> Lp-norm
  */
  Real norm(int a_p = 2,
            int a_comp = 0,
            int a_numcomp = 1) const
  {
    return m_fab.norm(a_p, a_comp, a_numcomp);
  }

  ///
  /**
     Returns sum of pow(fab[i,c],p): i in a_subbox, a_comp <= c <
     a_comp+a_numcomp, a_p >= 2 only
  */
  Real sumPow(const Box& a_subbox,
              int        a_p = 2,
              int        a_comp = 0,
              int        a_numcomp = 1) const

  {
    return m_fab.sumPow(a_subbox, a_p, a_comp, a_numcomp);
  }

  ///
  /**
     Returns the minimum value of given component of this NodeFArrayBox.
  */
  Real min(int a_comp = 0) const
  {
    return m_fab.min(a_comp);
  }

  ///
  /**
     Returns the minimum value of given component of this NodeFArrayBox in
     given a_subbox.

  */
  Real min(const Box& a_subbox,
           int        a_comp = 0) const
  {
    return m_fab.min(a_subbox, a_comp);
  }

  ///
  /**
     Returns the maximum value of given component of this NodeFArrayBox.
  */
  Real max(int a_comp = 0) const
  {
    return m_fab.max(a_comp);
  }

  ///
  /**
     Returns the maximum value of given component of this NodeFArrayBox in
     given a_subbox.

  */
  Real max(const Box& a_subbox,
           int        a_comp = 0) const

  {
    return m_fab.max(a_subbox, a_comp);
  }

  ///
  /**
     Finds location of minimum value in given component of this NodeFArrayBox.
  */
  IntVect minIndex(int a_comp = 0) const
  {
    return m_fab.minIndex(a_comp);
  }

  ///
  /**
     Returns location of minimum value in given component  of this NodeFArrayBox
     in given a_subbox.
  */
  IntVect minIndex(const Box& a_subbox,
                   int        a_comp = 0) const
  {
    return m_fab.minIndex(a_subbox, a_comp);
  }

  ///
  /**
     Returns location of maximum value in given component  of this NodeFArrayBox.
  */
  IntVect maxIndex(int a_comp = 0) const
  {
    return m_fab.maxIndex(a_comp);
  }

  ///
  /**
     Returns location of maximum value in given component of this NodeFArrayBox
     in given a_subbox.
  */
  IntVect maxIndex(const Box& a_subbox,
                   int        a_comp = 0) const
  {
    return m_fab.maxIndex(a_subbox, a_comp);
  }

  ///
  /**
     Computes a_mask array with value of 1 in cells where this NodeFArrayBox
     has value less than a_val, 0 otherwise.  a_mask is resized by this
     function.  The number of cells marked with 1 returned.
  */
  int maskLT(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const
  {
    return m_fab.maskLT(a_mask, a_val, a_comp);
  }

  ///
  /**
     Computes a_mask array with value of 1 in cells where this NodeFArrayBox
     has value less than or equal to a_val, 0 otherwise.  a_mask is
     resized by this function.  The number of cells marked with 1
     returned.
  */
  int maskLE(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const
  {
    return m_fab.maskLE(a_mask, a_val, a_comp);
  }

  ///
  /**
     Computes a_mask array with value of 1 in cells where this NodeFArrayBox
     has value equal to a_val, 0 otherwise.  a_mask is resized by this
     function.  The number of cells marked with 1 returned.

  */
  int maskEQ(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const
  {
    return m_fab.maskEQ(a_mask, a_val, a_comp);
  }

  ///
  /**
     Computes a_mask array with value of 1 in cells where this NodeFArrayBox
     has value greater than a_val, 0 otherwise.  a_mask is resized by this
     function.  The number of cells marked with 1 returned.
  */
  int maskGT(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const
  {
    return m_fab.maskGT(a_mask, a_val, a_comp);
  }

  ///
  /**
     Computes a_mask array with value of 1 in cells where this NodeFArrayBox
     has value greater than or equal to a_val, 0 otherwise.  a_mask is
     resized by this function.  The number of cells marked with 1 returned.
  */
  int maskGE(BaseFab<int>& a_mask,
             Real          a_val,
             int           a_comp = 0) const
  {
    return m_fab.maskGE(a_mask, a_val, a_comp);
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing each value with its absolute value.
  */
  void abs()
  {
    m_fab.abs();
  }

  ///
  /**
   */
  int nComp() const
  {
    return m_fab.nComp();
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing each value with its absolute value,
     for components (a_comp : a_comp+a_numcomp-1).
  */
  void abs(int a_comp,
           int a_numcomp = 1)
  {
    m_fab.abs(a_comp, a_numcomp);
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing eahc value with its absolute value,
     for components (a_comp : a_comp+a_numcomp-1) and within the a_subbox.
  */
  void abs (const Box& a_subbox,
            int        a_comp = 0,
            int        a_numcomp = 1)
  {
    m_fab.abs(a_subbox, a_comp, a_numcomp);
  }

  ///
  /**
     Returns sum of given component of NodeFArrayBox.
  */
  Real sum(int a_comp,
           int a_numcomp = 1) const
  {
    return m_fab.sum(a_comp, a_numcomp);
  }

  ///
  /**
     Returns sum of component of this NodeFArrayBox in given a_subbox.
  */
  Real sum(const Box& a_subbox,
           int        a_comp,
           int        a_numcomp = 1) const
  {
    return m_fab.sum(a_subbox, a_comp, a_numcomp);
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing each value x with a_r/x.
  */
  NodeFArrayBox& invert(Real a_r)
  {
    m_fab.invert(a_r);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing each value x with a_r/x.  For
     given range of components.
  */
  NodeFArrayBox& invert(Real a_r,
                        int  a_comp,
                        int  a_numcomp = 1)
  {
    m_fab.invert(a_r, a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing each value x with a_r/x.  For
     given range of components and within given a_subbox.
  */
  NodeFArrayBox& invert(Real       a_r,
                        const Box& a_subbox,
                        int        a_comp = 0,
                        int        a_numcomp = 1)
  {
    m_fab.invert(a_r, a_subbox, a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing each value with its additive
     inverse.  For given range of components and within given a_subbox.
  */
  NodeFArrayBox& negate(const Box& a_subbox,
                        int        a_comp = 0,
                        int        a_numcomp = 1)
  {
    m_fab.negate(a_subbox, a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing each value with its additive
     inverse.  For given range of components.
  */
  NodeFArrayBox& negate(int a_comp,
                        int a_numcomp = 1)
  {
    m_fab.negate(a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by replacing each value with its additive
     inverse.
  */
  NodeFArrayBox& negate()
  {
    m_fab.negate();
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by adding the scalar Real a_r to all values.  For
     given range of components and within given a_subbox.
  */
  NodeFArrayBox& plus(Real       a_r,
                      const Box& a_subbox,
                      int        a_comp = 0,
                      int        a_numcomp = 1)
  {
    m_fab.plus(a_r, a_subbox, a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by adding the scalar Real a_r to all values.  For
     given range of components.
  */
  NodeFArrayBox& plus(Real a_r,
                      int  a_comp,
                      int  a_numcomp = 1)
  {
    m_fab.plus(a_r,  a_comp, a_numcomp);
    return *this;
  }

  ///
  Real& operator() (const IntVect& a_p,
                    int N = 0)
  {
    return m_fab(a_p,N);
  }

  ///
  const Real& operator() (const IntVect& a_p,
                          int N = 0) const
  {
    return m_fab(a_p,N);
  }
  ///
  /**
     Modifies this NodeFArrayBox by adding the scalar Real a_r to all values.
  */
  NodeFArrayBox& operator += (Real a_r)
  {
    m_fab += a_r;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise addition of the values of the
     argument NodeFArrayBox.  You might come to grief if the domains of the
     NodeFArrayBoxes don't match.
  */
  NodeFArrayBox& operator += (const NodeFArrayBox& a_x)
  {
    m_fab += a_x.m_fab;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by adding the scalar Real a_r to all values.
  */
  NodeFArrayBox& plus(Real a_r)
  {
    m_fab += a_r;
    return *this;
  }

  NodeFArrayBox& plus_real(Real a_r)
  {
    return this->plus(a_r);
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise addition of the values of the
     argument NodeFArrayBox.  You might come to grief if the domains of the
     NodeFArrayBoxes don't match.  The same as += operator.

  */
  NodeFArrayBox& plus(const NodeFArrayBox& a_x)
  {
    m_fab += a_x.m_fab;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise scaled addition of the
     argument NodeFArrayBox (a[i] <- a[i] + a_scale * a_src[i]).  Uses domain
     of the intersection of these two NodeFArrayBoxes.
  */
  NodeFArrayBox& plus(const NodeFArrayBox& a_src,
                      const Real&          a_scale)
  {
    m_fab.plus(a_src.m_fab, a_scale);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise scaled addition of the
     argument NodeFArrayBox (a[i] <- a[i] + a_scale * a_src[i]).  Uses domain
     of the intersection of these two NodeFArrayBoxes.
  */
  NodeFArrayBox& plus(const NodeFArrayBox& a_src,
                      const Real&      a_scale,
                      int              a_srccomp,
                      int              a_destcomp,
                      int              a_numcomp = 1)
  {
    m_fab.plus(a_src.m_fab, a_scale, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise addition of values in the argument
     NodeFArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     to this NodeFArrayBox's components (a_destcomp : a_destcomp+a_numcomp-1)
     where the domains of the two NodeFArrayBoxes intersect.
  */
  NodeFArrayBox& plus(const NodeFArrayBox& a_src,
                      int              a_srccomp,
                      int              a_destcomp,
                      int              a_numcomp = 1)
  {
    m_fab.plus(a_src.m_fab,  a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise addition of values in the argument
     NodeFArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     to this NodeFArrayBox's components (a_destcomp : a_destcomp+numcomp-1)
     where the domain of this NodeFArrayBox intersects the a_subbox.  NOTE:
     a_subbox must be contained in this FAB.
  */
  NodeFArrayBox& plus(const NodeFArrayBox& a_src,
                      const Box&       a_subbox,
                      int              a_srccomp,
                      int              a_destcomp,
                      int              a_numcomp = 1)
  {
    m_fab.plus(a_src.m_fab, a_subbox,  a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise addition of values in the argument
     NodeFArrayBox.  Adds src's components (a_srccomp : a_srccomp+a_numcomp-1)
     in the Box a_srcbox to this NodeFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) in the Box a_destbox.  Corresponding locations
     within the two NodeFArrayBoxes are indexed relative to a_srcbox and a_destbox,
     and will in general not be the same.  The a_srcbox and a_destbox must be
     same size.  The results are UNDEFINED if the a_src and dest NodeFArrayBoxes
     are the same and the a_srcbox and a_destbox overlap.

  */
  NodeFArrayBox& plus(const NodeFArrayBox& a_src,
                      const Box&       a_srcbox,
                      const Box&       a_destbox,
                      int              a_srccomp,
                      int              a_destcomp,
                      int              a_numcomp = 1)
  {
    m_fab.plus(a_src.m_fab, a_srcbox, a_destbox,  a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  NodeFArrayBox& plus(const NodeFArrayBox& a_src,
                      const Box&       a_srcbox,
                      const Box&       a_destbox,
                      const Real&      a_scale,
                      int              a_srccomp,
                      int              a_destcomp,
                      int              a_numcomp = 1)
  {
    m_fab.plus(a_src.m_fab, a_srcbox, a_destbox,  a_scale, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by subtracting the scalar Real a_r to all values.
     Note: use plus(-a_r) for more general operations.
  */
  NodeFArrayBox& operator -= (Real a_r)
  {
    m_fab -= a_r;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise subtraction of the values of the
     argument NodeFArrayBox.  You might come to grief if the domains of the
     NodeFArrayBoxes don't match.
  */
  NodeFArrayBox& operator -= (const NodeFArrayBox& a_x)
  {
    m_fab -= a_x.m_fab;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise subtraction of the values of the
     argument NodeFArrayBox.  You might come to grief if the domains of the
     NodeFArrayBoxes don't match.  The same as -= operator.
  */
  NodeFArrayBox& minus(const NodeFArrayBox& a_x)
  {
    m_fab -= a_x.m_fab;
    return *this;
  }

  /**
     Modifies this NodeFArrayBox by pointwise subtraction of values in the
     argument NodeFArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) from this NodeFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domains of the two NodeFArrayBoxes
     intersect.
  */
  NodeFArrayBox& minus(const NodeFArrayBox& a_src,
                       int              a_srccomp,
                       int              a_destcomp,
                       int              a_numcomp = 1)
  {
    m_fab.minus(a_src.m_fab, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  /**
     Modifies this NodeFArrayBox by pointwise subtraction of values in the
     argument NodeFArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) from this NodeFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domain of this NodeFArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  NodeFArrayBox& minus(const NodeFArrayBox& a_src,
                       const Box&       a_subbox,
                       int              a_srccomp,
                       int              a_destcomp,
                       int              a_numcomp = 1)
  {
    m_fab.minus(a_src.m_fab, a_subbox, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise subtraction of values in the
     argument NodeFArrayBox.  Subtracts a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) in the Box a_srcbox from this NodeFArrayBox's
     components (a_destcomp : a_destcomp+a_numcomp-1) in the Box a_destbox.
     Corresponding locations within the two NodeFArrayBoxes are indexed relative
     to a_srcbox and a_destbox, and will in general not be the same.  The
     a_srcbox and a_destbox must be same size.  The results are UNDEFINED
     if the a_src and dest NodeFArrayBoxes are the same and the a_srcbox and
     a_destbox overlap.
  */
  NodeFArrayBox& minus(const NodeFArrayBox& a_src,
                       const Box&       a_srcbox,
                       const Box&       a_destbox,
                       int              a_srccomp,
                       int              a_destcomp,
                       int              a_numcomp = 1)
  {
    m_fab.minus(a_src.m_fab, a_srcbox, a_destbox, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by multiplying all values by the scalar Real a_r.
  */
  NodeFArrayBox& operator *= (Real a_r)
  {
    m_fab *= a_r;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by multiplying all values by the scalar Real a_r.
  */
  NodeFArrayBox& mult(Real a_r)
  {
    m_fab *= a_r;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by multiplying all values by the scalar
     Real a_r.  For given range of components.
  */
  NodeFArrayBox& mult(Real a_r,
                      int  a_comp,
                      int  a_numcomp = 1)
  {
    m_fab.mult(a_r, a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by multiplying all values by the scalar
     Real a_r.  For given range of components and within given a_subbox.
  */
  NodeFArrayBox& mult(Real       a_r,
                      const Box& a_subbox,
                      int        a_comp = 0,
                      int        a_numcomp = 1)
  {
    m_fab.mult(a_r, a_subbox, a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise multiplication of the values by the
     argument NodeFArrayBox.  You might come to grief if the domains of the
     NodeFArrayBoxes don't match.
  */
  NodeFArrayBox& operator *= (const NodeFArrayBox& a_x)
  {
    m_fab *= a_x.m_fab;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise multiplication by the values in the
     argument NodeFArrayBox.  You might come to grief if the domains of the
     NodeFArrayBoxes don't match.  The same as the *= operator.
  */
  NodeFArrayBox& mult(const NodeFArrayBox& a_x)
  {
    m_fab *= a_x.m_fab;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise multiplication by values in the
     argument NodeFArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) by this NodeFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domains of the two NodeFArrayBoxes
     intersect.
  */
  NodeFArrayBox& mult(const NodeFArrayBox& a_src,
                      int              a_srccomp,
                      int              a_destcomp,
                      int              a_numcomp = 1)
  {
    m_fab.mult(a_src.m_fab, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise multiplication by values in the
     argument NodeFArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) by this NodeFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) where the domain of this NodeFArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  NodeFArrayBox& mult(const NodeFArrayBox& a_src,
                      const Box&       a_subbox,
                      int              a_srccomp,
                      int              a_destcomp,
                      int              a_numcomp = 1)
  {
    m_fab.mult(a_src.m_fab, a_subbox, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise multiplication by values in the
     argument NodeFArrayBox.  Multiplies a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) in the Box a_srcbox by this NodeFArrayBox's
     components (a_destcomp : a_destcomp+a_numcomp-1) in the Box a_destbox.
     Corresponding locations within the two NodeFArrayBoxes are indexed relative
     to a_srcbox and a_destbox, and will in general not be the same.  The
     a_srcbox and a_destbox must be same size.  The results are UNDEFINED if
     the a_src and dest NodeFArrayBoxes are the same and the a_srcbox and a_destbox
     overlap.
  */
  NodeFArrayBox& mult(const NodeFArrayBox& a_src,
                      const Box&       a_srcbox,
                      const Box&       a_destbox,
                      int              a_srccomp,
                      int              a_destcomp,
                      int              a_numcomp = 1)
  {
    m_fab.mult(a_src.m_fab, a_srcbox, a_destbox, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by dividing all values by the scalar Real a_r.
  */
  NodeFArrayBox& operator /= (Real a_r)
  {
    m_fab /= a_r;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by dividing all values by the scalar Real a_r.
  */
  NodeFArrayBox& divide(Real a_r)
  {
    m_fab /= a_r;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by dividing all values by the scalar Real a_r.
     For given range of components.
  */
  NodeFArrayBox& divide(Real a_r,
                        int  a_comp,
                        int  a_numcomp = 1)
  {
    m_fab.divide(a_r, a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by dividing all values by the scalar Real
     a_r.  For given range of components and within given a_subbox.
  */
  NodeFArrayBox& divide(Real       a_r,
                        const Box& a_subbox,
                        int        a_comp = 0,
                        int        a_numcomp = 1)
  {
    m_fab.divide(a_r, a_subbox, a_comp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise division of the values by the
     argument NodeFArrayBox.  You might come to grief if the domains of the
     NodeFArrayBoxes don't match.
  */
  NodeFArrayBox& operator /= (const NodeFArrayBox& a_x)
  {
    m_fab /= a_x.m_fab;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise division by the values in the
     argument NodeFArrayBox.  You might come to grief if the domains of the
     NodeFArrayBoxes don't match.  The same as the /= operator.
  */
  NodeFArrayBox& divide(const NodeFArrayBox& a_x)
  {
    m_fab /= a_x.m_fab;
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise division by values in the argument
     NodeFArrayBox.  Divides this NodeFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) by a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) where the domains of the two NodeFArrayBoxes intersect.
  */
  NodeFArrayBox& divide(const NodeFArrayBox& a_src,
                        int              a_srccomp,
                        int              a_destcomp,
                        int              a_numcomp = 1)
  {
    m_fab.divide(a_src.m_fab, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise division by values in the argument
     NodeFArrayBox.  Divides this NodeFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) by a_src's components (a_srccomp :
     a_srccomp+a_numcomp-1) where the domain of this NodeFArrayBox intersects
     the a_subbox.  NOTE: a_subbox must be contained in this FAB.
  */
  NodeFArrayBox& divide(const NodeFArrayBox& a_src,
                        const Box&       a_subbox,
                        int              a_srccomp,
                        int              a_destcomp,
                        int              a_numcomp = 1)
  {
    m_fab.divide(a_src.m_fab, a_subbox, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

  ///
  /**
     Modifies this NodeFArrayBox by pointwise division by values in the argument
     NodeFArrayBox.  Divides this NodeFArrayBox's components (a_destcomp :
     a_destcomp+a_numcomp-1) in the Box a_destbox by a_src's components
     (a_srccomp : a_srccomp+a_numcomp-1) in the Box a_srcbox.  Corresponding
     locations within the two NodeFArrayBoxes are indexed relative to a_srcbox and
     a_destbox, and will in general not be the same.  The a_srcbox and
     a_destbox must be same size.  The results are UNDEFINED if the a_src and
     dest NodeFArrayBoxes are the same and the a_srcbox and a_destbox overlap.
  */
  NodeFArrayBox& divide(const NodeFArrayBox& a_src,
                        const Box&       a_srcbox,
                        const Box&       a_destbox,
                        int              a_srccomp,
                        int              a_destcomp,
                        int              a_numcomp = 1)
  {
    m_fab.divide(a_src.m_fab, a_srcbox, a_destbox, a_srccomp, a_destcomp, a_numcomp);
    return *this;
  }

protected:
  // CELL-centered box. data are defined on the surrounding nodes.
  Box m_box;

  // NODE-centered data
  FArrayBox m_fab;

private:
  /// These are disallowed for performance reasons

  NodeFArrayBox (const NodeFArrayBox&);
  NodeFArrayBox& operator = (const NodeFArrayBox&);

};

#include "NamespaceFooter.H"
#endif
