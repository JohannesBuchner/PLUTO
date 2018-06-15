#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeDotProduct.H
// petermc, 11 June 2003

#ifndef _NODEDOTPRODUCT_H_
#define _NODEDOTPRODUCT_H_

#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

/// Computes dot product of two sets of data on the same set of grids.
/**
    The dot product computed by this function is the sum over
    all components and all points of the product
    of <i>a_dataOne</i> and <i>a_dataTwo</i>,
    weighted by trapezoidal integration rule weights.
 */
Real DotProductNodes(/// first set of data
                     const BoxLayoutData<NodeFArrayBox>& a_dataOne,
                     /// second set of data
                     const BoxLayoutData<NodeFArrayBox>& a_dataTwo,
                     /// CELL-centered layout, on the nodes of which we calculate the dot product
                     const BoxLayout& a_dblIn);


/// Computes dot product of two sets of data on the same set of grids.
/**
    The dot product computed by this function is the sum over
    the components in <i>a_comps</i> and over all points of the product
    of <i>a_dataOne</i> and <i>a_dataTwo</i>,
    weighted by trapezoidal integration rule weights.
 */
Real DotProductNodes(/// first set of data
                     const BoxLayoutData<NodeFArrayBox>& a_dataOne,
                     /// second set of data
                     const BoxLayoutData<NodeFArrayBox>& a_dataTwo,
                     /// CELL-centered layout, on the nodes of which we calculate the dot product
                     const BoxLayout& a_dblIn,
                     /// components over which to take sum
                     const Interval& a_comps);


/// Computes dot product of two sets of data on the same set of grids.
/**
    The dot product computed by this function is the sum over
    the components in a_comps and over all VALID nodes of the product
    of <i>a_dataOne</i> and <i>a_dataTwo</i>.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_dataOne.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_dataOne.getBoxes();
    </tt>
 */
Real DotProductNodes(/// first set of data
                     const LevelData<NodeFArrayBox>& a_dataOne,
                     /// second set of data
                     const LevelData<NodeFArrayBox>& a_dataTwo,
                     /// physical domain
                     const ProblemDomain& a_domain,
                     /// external boundary nodes
                     const LayoutData< Vector<IntVectSet> >& a_IVSVext,
                     /// components over which to take sum
                     const Interval& a_comps);


/// Computes dot product of two sets of data on the same set of grids.
/**
    The dot product computed by this function is the sum over
    the components in a_comps and over all VALID nodes of the product
    of <i>a_dataOne</i> and <i>a_dataTwo</i>.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_dataOne.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_dataOne.getBoxes();
    </tt>
 */
Real DotProductNodes(/// first set of data
                     const LevelData<NodeFArrayBox>& a_dataOne,
                     /// second set of data
                     const LevelData<NodeFArrayBox>& a_dataTwo,
                     /// physical domain
                     const Box& a_domain,
                     /// external boundary nodes
                     const LayoutData< Vector<IntVectSet> >& a_IVSVext,
                     /// components over which to take sum
                     const Interval& a_comps);


/// Computes dot product of two sets of data on the same set of grids.
/**
    The dot product computed by this function is the sum over
    all components and all points of the product
    of <i>a_dataOne</i> and <i>a_dataTwo</i>,
    weighted by trapezoidal integration rule weights.

    This function computes exterior boundary nodes.  It is more efficient
    to precompute these boundary-node objects and call the function that
    takes <i>a_IVSVext</i> as an argument.
 */
Real DotProductNodes(/// first set of data
                     const LevelData<NodeFArrayBox>& a_dataOne,
                     /// second set of data
                     const LevelData<NodeFArrayBox>& a_dataTwo,
                     /// physical domain
                     const ProblemDomain& a_domain,
                     /// components over which to take sum
                     const Interval& a_comps);

/// Computes dot product of two sets of data on the same set of grids.
/**
    The dot product computed by this function is the sum over
    all components and all points of the product
    of <i>a_dataOne</i> and <i>a_dataTwo</i>,
    weighted by trapezoidal integration rule weights.

    This function computes exterior boundary nodes.  It is more efficient
    to precompute these boundary-node objects and call the function that
    takes <i>a_IVSVext</i> as an argument.
 */
Real DotProductNodes(/// first set of data
                     const LevelData<NodeFArrayBox>& a_dataOne,
                     /// second set of data
                     const LevelData<NodeFArrayBox>& a_dataTwo,
                     /// physical domain
                     const Box& a_domain,
                     /// components over which to take sum
                     const Interval& a_comps);

#include "NamespaceFooter.H"
#endif
