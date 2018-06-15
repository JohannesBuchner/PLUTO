#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeSetOperations.H
// petermc, 12 June 2003

#ifndef _NODESETOPERATIONS_H_
#define _NODESETOPERATIONS_H_

#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

/// Finds interior nodes of a DisjointBoxLayout in a particular box.
/** Returns nodes inside and surrounding CELL-centered box <i>a_box</i>
    that are NOT on the interface of (CELL-centered) <i>a_boxes</i>
    with the next coarser level.

    Example:  for two particular boxes in the four-grid layout below,
    the interior nodes are those that are marked "O".  For the
    second box, the nodes marked "Y" are interior nodes if and
    only if the domain is periodic in the vertical direction and
    the boxes extend to the top and bottom faces of the domain.
<pre>
                  +-----+                     +YY---+
                  |     |                     |OOOOO|
                  |     |                     |OOOOO|
                  |     |                     |OOOOO|
              +---+OOO--+---+             +---+OOOOO+---+
              |OOOOOOO      |             |      |      |
              |OOOOOOO      |             |      |      |
              |OOOOOOO      |             |      |      |
          +---+OOOOOO+------+         +---+------+------+
          |          |                |          |
          |          |                |          |
          |          |                |          |
          +----------+                +----------+
</pre>
    The set of interior nodes of a box is the disjoint union of
    the nodes surrounded by cells of the box, and the interior
    boundary nodes of the box.
 */
void interiorNodes(/// interior nodes of <i>a_box</i> in the layout <i>a_boxes</i>
                   IntVectSet& a_ivs,
                   /// CELL-centered physical domain containing <i>a_boxes</i>
                   const ProblemDomain& a_base_domain,
                   /// layout of CELL-centered boxes
                   const DisjointBoxLayout& a_boxes,
                   /// the box in <i>a_boxes</i> for which we find interior nodes
                   const Box& a_box);


/// Finds interior nodes of a DisjointBoxLayout in a particular box.
/** Returns nodes inside and surrounding CELL-centered box <i>a_box</i>
    that are NOT on the interface of (CELL-centered) <i>a_boxes</i>
    with the next coarser level.

    Example:  for two particular boxes in the four-grid layout below,
    the interior nodes are those that are marked "O".  For the
    second box, the nodes marked "Y" are interior nodes if and
    only if the domain is periodic in the vertical direction and
    the boxes extend to the top and bottom faces of the domain.
<pre>
                  +-----+                     +YY---+
                  |     |                     |OOOOO|
                  |     |                     |OOOOO|
                  |     |                     |OOOOO|
              +---+OOO--+---+             +---+OOOOO+---+
              |OOOOOOO      |             |      |      |
              |OOOOOOO      |             |      |      |
              |OOOOOOO      |             |      |      |
          +---+OOOOOO+------+         +---+------+------+
          |          |                |          |
          |          |                |          |
          |          |                |          |
          +----------+                +----------+
</pre>
    The set of interior nodes of a box is the disjoint union of
    the nodes surrounded by cells of the box, and the interior
    boundary nodes of the box.
 */
void interiorNodes(/// interior nodes of <i>a_box</i> in the layout <i>a_boxes</i>
                   IntVectSet& a_ivs,
                   /// CELL-centered physical domain containing <i>a_boxes</i>
                   const Box& a_base_domain,
                   /// layout of CELL-centered boxes
                   const DisjointBoxLayout& a_boxes,
                   /// the box in <i>a_boxes</i> for which we find interior nodes
                   const Box& a_box);


/// Finds interior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each box in <i>a_boxes</i>,
    a Vector of IntVectSets that contain the indices of the
    nodes that lie on faces of the box but are not on the
    interface of <i>a_boxes</i> with the next coarser level.

    Example:  for two particular boxes, marked "me", in the
    four-grid layout below, the interior boundary nodes are those
    that are marked with "O".  For the second box, the nodes marked "Y"
    are also interior boundary nodes if and only if the domain is
    periodic in the vertical direction and the grids extend to the
    top and bottom faces of the domain.
<pre>
                  +-----+                     +YY---+
                  |     |                     |     |
                  |     |                     | me  |
                  |     |                     |     |
              +---+OOO--+---+             +---+OOOOO+---+
              |      O      |             |      |      |
              |  me  O      |             |      |      |
              |      O      |             |      |      |
          +---+OOOOOO+------+         +---+------+------+
          |          |                |          |
          |          |                |          |
          |          |                |          |
          +----------+                +----------+
</pre>
    The result <i>a_IVSV</i> is needed by copyInteriorNodes().

    <b>Question</b>:  Why not simply store the indices of all interior
    boundary nodes of <i>a_boxes</i> in a single IntVectSet, rather than
    a LayoutData< Vector<IntVectSet> >? <br>
    <b>Answer</b>:  Because of the way IntVectSets are stored, this would be
    very slow.

    The different components of the Vector<IntVectSet> contain
    indices of nodes on different faces of the box.
 */
void interiorBoundaryNodes(/// object containing interior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<IntVectSet> >& a_IVSV,
                           /// layout of CELL-centered boxes
                           const DisjointBoxLayout& a_boxes,
                           /// CELL-centered physical domain containing <i>a_boxes</i>
                           const ProblemDomain& a_domain);


/// Finds interior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each box in <i>a_boxes</i>,
    a Vector of IntVectSets that contain the indices of the
    nodes that lie on faces of the box but are not on the
    interface of <i>a_boxes</i> with the next coarser level.

    Example:  for two particular boxes, marked "me", in the
    four-grid layout below, the interior boundary nodes are those
    that are marked with "O".  For the second box, the nodes marked "Y"
    are also interior boundary nodes if and only if the domain is
    periodic in the vertical direction and the grids extend to the
    top and bottom faces of the domain.
<pre>
                  +-----+                     +YY---+
                  |     |                     |     |
                  |     |                     | me  |
                  |     |                     |     |
              +---+OOO--+---+             +---+OOOOO+---+
              |      O      |             |      |      |
              |  me  O      |             |      |      |
              |      O      |             |      |      |
          +---+OOOOOO+------+         +---+------+------+
          |          |                |          |
          |          |                |          |
          |          |                |          |
          +----------+                +----------+
</pre>
    The result <i>a_IVSV</i> is needed by copyInteriorNodes().

    <b>Question</b>:  Why not simply store the indices of all interior
    boundary nodes of <i>a_boxes</i> in a single IntVectSet, rather than
    a LayoutData< Vector<IntVectSet> >? <br>
    <b>Answer</b>:  Because of the way IntVectSets are stored, this would be
    very slow.

    The different components of the Vector<IntVectSet> contain
    indices of nodes on different faces of the box.
 */
void interiorBoundaryNodes(/// object containing interior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<IntVectSet> >& a_IVSV,
                           /// layout of CELL-centered boxes
                           const DisjointBoxLayout& a_boxes,
                           /// CELL-centered physical domain containing <i>a_boxes</i>
                           const Box& a_domain);


/// Finds interior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each destination box,
    a Vector of IntVectSets that together represent the nodes
    of the destination box that are also interior boundary nodes
    of the source boxes.

    Specifically, if IBN(<i>a_src</i>) is the set of interior boundary nodes
    of <i>a_src</i>, then the Vector<IntVectSet> <i>a_IVSV</i>[dit()] holds
    the intersection of IBN(<i>a_src</i>) and the nodes
    surrounding <i>a_dest</i>[dit()].

    The result <i>a_IVSV</i> is used in copyInteriorNodes().
 */
void interiorBoundaryNodes(/// object containing interior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<IntVectSet> >& a_IVSV,
                           /// layout of CELL-centered destination boxes
                           const DisjointBoxLayout& a_dest,
                           /// layout of CELL-centered source boxes
                           const DisjointBoxLayout& a_src,
                           /// CELL-centered physical domain containing <i>a_dest</i> and <i>a_src</i>
                           const ProblemDomain& a_domain);


/// Finds interior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each destination box,
    a Vector of IntVectSets that together represent the nodes
    of the destination box that are also interior boundary nodes
    of the source boxes.

    Specifically, if IBN(<i>a_src</i>) is the set of interior boundary nodes
    of <i>a_src</i>, then the Vector<IntVectSet> <i>a_IVSV</i>[dit()] holds
    the intersection of IBN(<i>a_src</i>) and the nodes
    surrounding <i>a_dest</i>[dit()].

    The result <i>a_IVSV</i> is used in copyInteriorNodes().
 */
void interiorBoundaryNodes(/// object containing interior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<IntVectSet> >& a_IVSV,
                           /// layout of CELL-centered destination boxes
                           const DisjointBoxLayout& a_dest,
                           /// layout of CELL-centered source boxes
                           const DisjointBoxLayout& a_src,
                           /// CELL-centered physical domain containing <i>a_dest</i> and <i>a_src</i>
                           const Box& a_domain);


/// Finds exterior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each box, the exterior boundary
    nodes of the box.

    The exterior boundary nodes are the nodes on a level that are
    not interior nodes:  they lie on the physical boundary and on
    the interface with the coarser level.

    Example:  for two particular boxes, marked "me", in the
    four-grid layout below, the exterior boundary nodes are those
    that are marked "O".  For the second box, the nodes marked "Y"
    are also exterior boundary nodes, unless the domain is periodic
    in the vertical direction and the grids extend to the
    top and bottom faces of the domain.
<pre>
                  +-----+                     OYYOOOO
                  |     |                     O     O
                  |     |                     O me  O
                  |     |                     O     O
              OOOOO--+--+---+             +---O--+--O---+
              O      |      |             |      |      |
              O  me  |      |             |      |      |
              O      |      |             |      |      |
          +---O------O------+         +---+------+------+
          |          |                |          |
          |          |                |          |
          |          |                |          |
          +----------+                +----------+
</pre>
    The argument <i>a_interior</i> is obtained from<br>
    <tt>
    interiorBoundaryNodes(a_interior, a_boxes, a_domain)
    </tt>
    where <i>a_domain</i> is the physical domain.

    You don't need the physical domain <i>a_domain</i> as an argument
    to this function, even to get periodicities, because everything needed
    from <i>a_domain</i> is put into <i>a_interior</i>.
 */
void exteriorBoundaryNodes(/// object containing exterior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<IntVectSet> >& a_exterior,
                           /// object containing interior boundary nodes of <i>a_boxes</i>
                           const LayoutData< Vector<IntVectSet> >& a_interior,
                           /// layout of CELL-centered boxes
                           const DisjointBoxLayout& a_boxes);


/// Sets data to zero on specified boundary nodes.
/** Sets <i>a_dest</i> to zero on the nodes in <i>a_IVSV</i>.

    Specifically, <i>a_dest</i>[dit()] is set to zero on
    the nodes in <i>a_IVSV</i>[dit()].

    We usually call this function with <i>a_IVSV</i> representing the
    exterior boundary nodes of <i>a_dest</i>.  These are obtained with
    the sequence of calls:<br>
    <tt>
    interiorBoundaryNodes(a_IVSVint, a_dest.boxLayout(), a_domain)
    <tt><br>
    where <i>a_domain</i> is the physical domain, and<br>
    <tt>
    exteriorBoundaryNodes(a_IVSV, a_IVSVint, a_dest.boxLayout())
    <tt>.
*/
void zeroBoundaryNodes(/// the data
                       BoxLayoutData<NodeFArrayBox>& a_dest,
                       /// object containing boundary nodes
                       const LayoutData< Vector<IntVectSet> >& a_IVSV);

/// Copies data at interior nodes from one LevelData<NodeFArrayBox> to another.
/** Copies data from <i>a_src</i> to <i>a_dest</i>,
    over the interior nodes of <i>a_src</i>.
    They need not have the same layout.

    The argument <i>a_IVSV</i> is obtained from<br>
    <tt>
    interiorBoundaryNodes(a_IVSV,
    a_dest.boxLayout(), a_src.boxLayout(), a_domain)
    </tt><br>
    where <i>a_domain</i> is the physical domain.

    If <i>a_dest</i> and <i>a_src</i> have the same layout,
    then we may obtain <i>a_IVSV</i> more efficiently from<br>
    <tt>
    interiorBoundaryNodes(a_IVSV, a_src.boxLayout(), a_domain)
    </tt>.<br>
    Do not use this shorter form if <i>a_dest</i> and <i>a_src</i>
    have different layouts.

    Example:  if <i>a_src</i> has the four-grid layout on the left below, then
    the nodes from which data are copied to <i>a_dest</i> are marked with
    "X" and "O" on the right below.  (The distinction between "X" and "O"
    is that nodes marked "O" are on interior boundaries of the layout.)
    Data are copied from the nodes marked "Y" if and only if the domain
    is periodic in the vertical direction and the grids extend to the
    top and bottom faces of the domain.
<pre>
                  +-----+                     +YY---+
                  |     |                     |OOOOO|
                  |     |                     |OOOOO|
                  |     |                     |OOOOO|
              +---+--+--+---+             +---+XXXXX+---+
              |      |      |             |OOOOOOXOOOOOO|
              |      |      |             |OOOOOOXOOOOOO|
              |      |      |             |OOOOOOXOOOOOO|
          +---+------+------+         +---+XXXXXX+------+
          |          |                |OOOOOOOOOO|
          |          |                |OOOOOOOOOO|
          |          |                |OOOOOOOOOO|
          +----------+                +--------YY+
</pre>
    Both pictures represent the layout of <i>a_src</i>.
    The layout of <i>a_dest</i> may be completely different from that
    of <i>a_src</i>.

    Data at nodes in <i>a_dest</i> that are not interior boundary nodes
    of <i>a_src</i> remain unchanged by this function.
 */
void copyInteriorNodes(/// destination data
                       LevelData<NodeFArrayBox>& a_dest,
                       /// source data
                       const LevelData<NodeFArrayBox>& a_src,
                       /// object containing interior boundary nodes of source layout that are also nodes of destination layout
                       const LayoutData< Vector<IntVectSet> >& a_IVSV);

/// Determines which IntVectSets are rectangular.
/** Sets the bit <i>a_IVSVfull</i>[dit()][comp] to be true if and only if
    <i>a_IVSV</i>[dit()][comp] is a full rectangle.
 */
void
fullIntVectSets(/// object containing bit settings of whether each component is rectangular
                LayoutData< BitSet >& a_IVSVfull,
                /// object containing interior boundary nodes of a layout
                const LayoutData< Vector<IntVectSet> >& a_IVSV);

/*
  These routines are the same as above, but use LayoutData< Vector<Box> >
  instead of LayoutData< Vector<IntVectSet> >.
*/

/// Finds interior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each box in <i>a_boxes</i>,
    a Vector of IntVectSets that contain the indices of the
    nodes that lie on faces of the box but are not on the
    interface of <i>a_boxes</i> with the next coarser level.

    Example:  for two particular boxes, marked "me", in the
    four-grid layout below, the interior boundary nodes are those
    that are marked with "O".  For the second box, the nodes marked "Y"
    are also interior boundary nodes if and only if the domain is
    periodic in the vertical direction and the grids extend to the
    top and bottom faces of the domain.
<pre>
                  +-----+                     +YY---+
                  |     |                     |     |
                  |     |                     | me  |
                  |     |                     |     |
              +---+OOO--+---+             +---+OOOOO+---+
              |      O      |             |      |      |
              |  me  O      |             |      |      |
              |      O      |             |      |      |
          +---+OOOOOO+------+         +---+------+------+
          |          |                |          |
          |          |                |          |
          |          |                |          |
          +----------+                +----------+
</pre>
    The result <i>a_IVSV</i> is needed by copyInteriorNodes().

    <b>Question</b>:  Why not simply store the indices of all interior
    boundary nodes of <i>a_boxes</i> in a single IntVectSet, rather than
    a LayoutData< Vector<Box> >? <br>
    <b>Answer</b>:  Because of the way IntVectSets are stored, this would be
    very slow.

    The different components of the Vector<Box> contain
    indices of nodes on different faces of the box.
 */
void interiorBoundaryNodes(/// object containing interior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<Box> >& a_IVSV,
                           /// layout of CELL-centered boxes
                           const DisjointBoxLayout& a_boxes,
                           /// CELL-centered physical domain containing <i>a_boxes</i>
                           const Box& a_domain);


/// Finds interior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each box in <i>a_boxes</i>,
    a Vector of IntVectSets that contain the indices of the
    nodes that lie on faces of the box but are not on the
    interface of <i>a_boxes</i> with the next coarser level.

    Example:  for two particular boxes, marked "me", in the
    four-grid layout below, the interior boundary nodes are those
    that are marked with "O".  For the second box, the nodes marked "Y"
    are also interior boundary nodes if and only if the domain is
    periodic in the vertical direction and the grids extend to the
    top and bottom faces of the domain.
<pre>
                  +-----+                     +YY---+
                  |     |                     |     |
                  |     |                     | me  |
                  |     |                     |     |
              +---+OOO--+---+             +---+OOOOO+---+
              |      O      |             |      |      |
              |  me  O      |             |      |      |
              |      O      |             |      |      |
          +---+OOOOOO+------+         +---+------+------+
          |          |                |          |
          |          |                |          |
          |          |                |          |
          +----------+                +----------+
</pre>
    The result <i>a_IVSV</i> is needed by copyInteriorNodes().

    <b>Question</b>:  Why not simply store the indices of all interior
    boundary nodes of <i>a_boxes</i> in a single IntVectSet, rather than
    a LayoutData< Vector<Box> >? <br>
    <b>Answer</b>:  Because of the way IntVectSets are stored, this would be
    very slow.

    The different components of the Vector<Box> contain
    indices of nodes on different faces of the box.
 */
void interiorBoundaryNodes(/// object containing interior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<Box> >& a_IVSV,
                           /// layout of CELL-centered boxes
                           const DisjointBoxLayout& a_boxes,
                           /// CELL-centered physical domain containing <i>a_boxes</i>
                           const ProblemDomain& a_domain);


/// Finds interior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each destination box,
    a Vector of IntVectSets that together represent the nodes
    of the destination box that are also interior boundary nodes
    of the source boxes.

    Specifically, if IBN(<i>a_src</i>) is the set of interior boundary nodes
    of <i>a_src</i>, then the Vector<Box> <i>a_IVSV</i>[dit()] holds
    the intersection of IBN(<i>a_src</i>) and the nodes
    surrounding <i>a_dest</i>[dit()].

    The result <i>a_IVSV</i> is used in copyInteriorNodes().
 */
void interiorBoundaryNodes(/// object containing interior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<Box> >& a_IVSV,
                           /// layout of CELL-centered destination boxes
                           const DisjointBoxLayout& a_dest,
                           /// layout of CELL-centered source boxes
                           const DisjointBoxLayout& a_src,
                           /// CELL-centered physical domain containing <i>a_dest</i> and <i>a_src</i>
                           const Box& a_domain);


/// Finds interior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each destination box,
    a Vector of IntVectSets that together represent the nodes
    of the destination box that are also interior boundary nodes
    of the source boxes.

    Specifically, if IBN(<i>a_src</i>) is the set of interior boundary nodes
    of <i>a_src</i>, then the Vector<Box> <i>a_IVSV</i>[dit()] holds
    the intersection of IBN(<i>a_src</i>) and the nodes
    surrounding <i>a_dest</i>[dit()].

    The result <i>a_IVSV</i> is used in copyInteriorNodes().
 */
void interiorBoundaryNodes(/// object containing interior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<Box> >& a_IVSV,
                           /// layout of CELL-centered destination boxes
                           const DisjointBoxLayout& a_dest,
                           /// layout of CELL-centered source boxes
                           const DisjointBoxLayout& a_src,
                           /// CELL-centered physical domain containing <i>a_dest</i> and <i>a_src</i>
                           const ProblemDomain& a_domain);


/// Finds exterior boundary nodes of a DisjointBoxLayout.
/** Returns a LayoutData containing, for each box, the exterior boundary
    nodes of the box.

    The exterior boundary nodes are the nodes on a level that are
    not interior nodes:  they lie on the physical boundary and on
    the interface with the coarser level.

    Example:  for two particular boxes, marked "me", in the
    four-grid layout below, the exterior boundary nodes are those
    that are marked "O".  For the second box, the nodes marked "Y"
    are also exterior boundary nodes, unless the domain is periodic
    in the vertical direction and the grids extend to the
    top and bottom faces of the domain.
<pre>
                  +-----+                     OYYOOOO
                  |     |                     O     O
                  |     |                     O me  O
                  |     |                     O     O
              OOOOO--+--+---+             +---O--+--O---+
              O      |      |             |      |      |
              O  me  |      |             |      |      |
              O      |      |             |      |      |
          +---O------O------+         +---+------+------+
          |          |                |          |
          |          |                |          |
          |          |                |          |
          +----------+                +----------+
</pre>
    The argument <i>a_interior</i> is obtained from<br>
    <tt>
    interiorBoundaryNodes(a_interior, a_boxes, a_domain)
    </tt>
    where <i>a_domain</i> is the physical domain.

    You don't need the physical domain <i>a_domain</i> as an argument
    to this function, even to get periodicities, because everything needed
    from <i>a_domain</i> is put into <i>a_interior</i>.
 */
void exteriorBoundaryNodes(/// object containing exterior boundary nodes of <i>a_boxes</i>
                           LayoutData< Vector<Box> >& a_exterior,
                           /// object containing interior boundary nodes of <i>a_boxes</i>
                           const LayoutData< Vector<Box> >& a_interior,
                           /// layout of CELL-centered boxes
                           const DisjointBoxLayout& a_boxes);

/// Copies data at interior nodes from one LevelData<NodeFArrayBox> to another.
/** Copies data from <i>a_src</i> to <i>a_dest</i>,
    over the interior nodes of <i>a_src</i>.
    They need not have the same layout.

    The argument <i>a_IVSV</i> is obtained from<br>
    <tt>
    interiorBoundaryNodes(a_IVSV,
    a_dest.boxLayout(), a_src.boxLayout(), a_domain)
    </tt><br>
    where <i>a_domain</i> is the physical domain.

    If <i>a_dest</i> and <i>a_src</i> have the same layout,
    then we may obtain <i>a_IVSV</i> more efficiently from<br>
    <tt>
    interiorBoundaryNodes(a_IVSV, a_src.boxLayout(), a_domain)
    </tt>.<br>
    Do not use this shorter form if <i>a_dest</i> and <i>a_src</i>
    have different layouts.

    Example:  if <i>a_src</i> has the four-grid layout on the left below, then
    the nodes from which data are copied to <i>a_dest</i> are marked with
    "X" and "O" on the right below.  (The distinction between "X" and "O"
    is that nodes marked "O" are on interior boundaries of the layout.)
    Data are copied from the nodes marked "Y" if and only if the domain
    is periodic in the vertical direction and the grids extend to the
    top and bottom faces of the domain.
<pre>
                  +-----+                     +YY---+
                  |     |                     |OOOOO|
                  |     |                     |OOOOO|
                  |     |                     |OOOOO|
              +---+--+--+---+             +---+XXXXX+---+
              |      |      |             |OOOOOOXOOOOOO|
              |      |      |             |OOOOOOXOOOOOO|
              |      |      |             |OOOOOOXOOOOOO|
          +---+------+------+         +---+XXXXXX+------+
          |          |                |OOOOOOOOOO|
          |          |                |OOOOOOOOOO|
          |          |                |OOOOOOOOOO|
          +----------+                +--------YY+
</pre>
    Both pictures represent the layout of <i>a_src</i>.
    The layout of <i>a_dest</i> may be completely different from that
    of <i>a_src</i>.

    Data at nodes in <i>a_dest</i> that are not interior boundary nodes
    of <i>a_src</i> remain unchanged by this function.
 */
void copyInteriorNodes(/// destination data
                       LevelData<NodeFArrayBox>& a_dest,
                       /// source data
                       const LevelData<NodeFArrayBox>& a_src,
                       /// object containing interior boundary nodes of source layout that are also nodes of destination layout
                       const LayoutData< Vector<Box> >& a_IVSV);


/// Sets data to zero on specified boundary nodes.
/** Sets <i>a_dest</i> to zero on the nodes in <i>a_IVSV</i>.

    Specifically, <i>a_dest</i>[dit()] is set to zero on
    the nodes in <i>a_IVSV</i>[dit()].

    We usually call this function with <i>a_IVSV</i> representing the
    exterior boundary nodes of <i>a_dest</i>.  These are obtained with
    the sequence of calls:<br>
    <tt>
    interiorBoundaryNodes(a_IVSVint, a_dest.boxLayout(), a_domain)
    <tt><br>
    where <i>a_domain</i> is the physical domain, and<br>
    <tt>
    exteriorBoundaryNodes(a_IVSV, a_IVSVint, a_dest.boxLayout())
    <tt>.
*/
void zeroBoundaryNodes(/// the data
                       BoxLayoutData<NodeFArrayBox>& a_dest,
                       /// object containing boundary nodes
                       const LayoutData< Vector<Box> >& a_IVSV);

/*
  Following functions are for Vector<Box> calculus.
*/

///
/*
  To a_boxes, append the new points of the boxes in a_new.
*/
void appendBoxes(Vector<Box>&         a_boxes,
                 const Vector<Box>&   a_new);

///
/*
  From a_boxes, remove all points of a_remove.
*/
void removeBoxFromBoxes(Vector<Box>&   a_boxes,
                        const Box&     a_remove);

///
/*
  Return in a_boxes the result of removing the boxes in a_remove from a_base.
*/
void removeBoxesFromBox(Vector<Box>&         a_boxes,
                        const Vector<Box>&   a_remove,
                        const Box&           a_base);

///
/*
  Return in a_boxes the result of removing the points of a_remove from a_base.
*/
void removeBoxFromBox(Vector<Box>&   a_boxes,
                      const Box&     a_remove,
                      const Box&     a_base);

///
/*
  From a_boxes, remove all points in the boxes in a_remove.
*/
void removeBoxesFromBoxes(Vector<Box>&         a_boxes,
                          const Vector<Box>&   a_remove);

///
/*
  Replace each box in a_boxes with its intersection with the union of
  the boxes in a_new.
*/
void intersectBoxes(Vector<Box>&         a_boxes,
                    const Vector<Box>&   a_new);

///
/*
  Replace each box in a_boxes with its intersection with the box a_new.
*/
void intersectBoxes(Vector<Box>&   a_boxes,
                    const Box&     a_new);

///
/*
  Shift all boxes in a_boxes by a_offset.
*/
void shiftBoxes(Vector<Box>&     a_boxes,
                const IntVect&   a_offset);

///
/*
  Convert box centerings from CELL to NODE.
*/
void cellsToNodes(Vector<Box>&   a_boxes);

///
/*
  Convert box centerings from NODE to CELL.
*/
void nodesToCells(Vector<Box>&   a_boxes);

///
/*
  At each level, return exterior boundary nodes.
  At each level except for the finest, return the interior boundary nodes
  of the coarsened finer level.
*/
void exteriorAndInteriorNodes(
                              /// object containing exterior boundary nodes of <i>a_boxes</i> at each level
                              Vector< LayoutData< Vector<Box> >* >& a_exterior,
                              /// object containing interior boundary nodes of coarsened <i>a_boxes</i> at each level except finest
                              Vector< LayoutData< Vector<Box> >* >& a_intFinerCoarsened,
                              /// layouts of CELL-centered boxes
                              const Vector<DisjointBoxLayout>& a_layouts,
                              /// domain of each level
                              const Vector<ProblemDomain>& a_domain,
                              /// refinement ratio between levels
                              const Vector<int>& a_nRefFine);


///
/*
  At each level, return exterior boundary nodes.
  At each level except for the finest, return the interior boundary nodes
  of the coarsened finer level.
*/
void exteriorAndInteriorNodes(
                              /// object containing exterior boundary nodes of <i>a_boxes</i> at each level
                              Vector< LayoutData< Vector<IntVectSet> >* >& a_exterior,
                              /// object containing interior boundary nodes of coarsened <i>a_boxes</i> at each level except finest
                              Vector< LayoutData< Vector<IntVectSet> >* >& a_intFinerCoarsened,
                              /// layouts of CELL-centered boxes
                              const Vector<DisjointBoxLayout>& a_layouts,
                              /// domain of each level
                              const Vector<ProblemDomain>& a_domain,
                              /// refinement ratio between levels
                              const Vector<int>& a_nRefFine);


/// Returns a mask on interior nodes of a DisjointBoxLayout.
/** Returns a LevelData on nodes of a_dest that is
    a_onoff on interior nodes of *a_srcPtr (or of a_dest if a_srcPtr is null);
    1 - a_onoff on all other nodes.
 */
void getMaskInteriorNodes(/// mask built on a_dest, representing interior nodes of a_src.
                          LevelData<NodeFArrayBox>& a_mask,
                          /// layout of CELL-centered destination boxes
                          const DisjointBoxLayout& a_dest,
                          /// layout of CELL-centered source boxes
                          const DisjointBoxLayout* a_srcPtr,
                          /// CELL-centered physical domain containing <i>a_dest</i> and <i>a_src</i>
                          const Box& a_domain,
                          /// 1 or 0, for interior nodes
                          int a_onoff = 1);


/// Returns a mask on interior nodes of a DisjointBoxLayout.
/** Returns a LevelData on nodes of a_dest that is
    a_onoff on interior nodes of *a_srcPtr (or of a_dest if a_srcPtr is null);
    1 - a_onoff on all other nodes.
 */
void getMaskInteriorNodes(/// mask built on a_dest, representing interior boundary nodes of a_src.
                          LevelData<NodeFArrayBox>& a_mask,
                          /// layout of CELL-centered destination boxes
                          const DisjointBoxLayout& a_dest,
                          /// layout of CELL-centered source boxes
                          const DisjointBoxLayout* a_srcPtr,
                          /// CELL-centered physical domain containing <i>a_dest</i> and <i>a_src</i>
                          const ProblemDomain& a_domain,
                          /// 1 or 0, for interior nodes
                          int a_onoff = 1);


/// Returns a mask on interior nodes of a DisjointBoxLayout.
/** Returns a LevelData on nodes of a_dest that is
    1. on interior nodes of a_dest,
    0. on all other nodes.
 */
void getMaskInteriorNodes(/// mask built on a_dest, representing interior nodes of a_dest.
                          LevelData<NodeFArrayBox>& a_mask,
                          /// layout of CELL-centered destination boxes
                          const DisjointBoxLayout& a_dest,
                          /// CELL-centered physical domain containing <i>a_dest</i>
                          const Box& a_domain);


/// Returns a mask on valid nodes on a level.
/** Returns a LevelData on a_layout that is
    1. on valid nodes of a_layout,
    0. on all other nodes.
 */
void getMaskValidNodes(/// mask built on a_dest, representing valid nodes of a_layout.
                       LevelData<NodeFArrayBox>& a_mask,
                       /// layout of CELL-centered boxes on this level
                       const DisjointBoxLayout& a_layout,
                       /// layout of CELL-centered coarsened boxes on finer level
                       const DisjointBoxLayout* a_finerCoarsenedPtr,
                       /// CELL-centered physical domain containing <i>a_layout</i> and <i>a_finerCoarsened</i>
                       const ProblemDomain& a_domain);


/// Returns a mask on valid nodes on a level.
/** Returns a LevelData on a_layout that is
    1. on valid nodes of a_layout,
    0. on all other nodes.
 */
void getMaskValidNodes(/// mask built on a_dest, representing valid nodes of a_layout.
                       LevelData<NodeFArrayBox>& a_mask,
                       /// layout of CELL-centered boxes on this level
                       const DisjointBoxLayout& a_layout,
                       /// layout of CELL-centered coarsened boxes on finer level
                       const DisjointBoxLayout* a_finerCoarsenedPtr,
                       /// CELL-centered physical domain containing <i>a_layout</i> and <i>a_finerCoarsened</i>
                       const Box& a_domain);


///
/*
  At each level, return mask for valid nodes.
*/
void getMaskValidNodes(
                       /// at each level, a mask for valid nodes
                       Vector< LevelData<NodeFArrayBox>* >& a_masks,
                       /// layouts of CELL-centered boxes
                       const Vector<DisjointBoxLayout>& a_layouts,
                       /// domain of each level
                       const Vector<ProblemDomain>& a_domain,
                       /// refinement ratio between levels
                       const Vector<int>& a_nRefFine);


///
/** Return the NODE-centered layout of NODEs surrounding a
    CELL-centered layout.
*/
void getSurroundingNodes(BoxLayout&         a_gridsNodes,
                         const BoxLayout&   a_gridsCells);


///
/** Return the CELL-centered layout of CELLs enclosed by a
    NODE-centered layout.
*/
void getEnclosedCells(BoxLayout&         a_gridsCells,
                      const BoxLayout&   a_gridsNodes);

///
/** Divide a_intFab by 2 on faces of a_bx (and by 4 on edges, 8 on corners)
 */
void halveIntFaces(BaseFab<int>& a_intFab,
                   const Box& a_bx);

#include "NamespaceFooter.H"
#endif
