#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeIntegrals.H
// petermc, 17 Oct 2005
// used to be in Norms.H

#ifndef _NODEINTEGRALS_H_
#define _NODEINTEGRALS_H_

#include "LevelData.H"
#include "NodeFArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

/// Computes integral of a subbox of a NodeFArrayBox.
/** Computes integral of the data in a subbox of <i>a_nfab</i>.
 */
Real integral(/// the data on which to take the integral
              const NodeFArrayBox& a_nfab,
              /// mesh spacing at this level
              const Real a_dx,
              /// CELL-centered subbox over which to take integral
              const Box& a_subbox,
              /// starting component
              const int a_startComp = 0,
              /// number of components
              const int a_numComp = 1);

/// Computes integral of a NodeFArrayBox.
/** Computes integral of the data in <i>a_nfab</i>.
 */
Real integral(/// the data on which to take the integral
              const NodeFArrayBox& a_nfab,
              /// mesh spacing at this level
              const Real a_dx,
              /// starting component
              const int a_startComp = 0,
              /// number of components
              const int a_numComp = 1);

/// Computes integral of a BoxLayoutData<NodeFArrayBox>.
/** Computes integral of <i>a_layout</i>.  Ignores ghost cells for LevelData.

    This is an internal function that should NOT normally
    be called by the user.
 */
Real integral(/// the data on which to take the norm
              const BoxLayoutData<NodeFArrayBox>& a_layout,
              /// mesh spacing at this level
              const Real a_dx,
              /// interval of components to use
              const Interval& a_interval,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);

/// Computes integral of multilevel array, counting only valid nodes at each level.
/**
    Computes integral of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.
 */
Real integral(/// multilevel data on which to take the norm
              const Vector<LevelData<NodeFArrayBox>* >& a_phi,
              /// CELL-centered physical domain on each level
              const Vector<ProblemDomain>& a_domain,
              /// refinement ratios between successive levels
              const Vector<int>& a_nRefFine,
              /// mesh spacing at base level <i>a_lBase</i>
              const Real a_dxCrse,
              /// components to use in computing norm
              const Interval a_comps,
              /// index of base level to use in Vectors
              const int a_lBase,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);

/// Computes integral of multilevel array, counting only valid nodes at each level.
/**
    Computes integral of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.
 */
Real integral(/// multilevel data on which to take the norm
              const Vector<LevelData<NodeFArrayBox>* >& a_phi,
              /// CELL-centered physical domain on each level
              const Vector<Box>& a_domain,
              /// refinement ratios between successive levels
              const Vector<int>& a_nRefFine,
              /// mesh spacing at base level <i>a_lBase</i>
              const Real a_dxCrse,
              /// components to use in computing norm
              const Interval a_comps,
              /// index of base level to use in Vectors
              const int a_lBase,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);


/// Computes integral of single-level array, counting only valid nodes.
/**
    Computes integral of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    This function computes exterior boundary nodes, and interior boundary
    nodes of the coarsened finer level, if any.  It is more efficient
    to precompute these boundary-node objects and call one of the other
    integral functions below.
*/
Real integral(/// the data on which to take the norm
              const LevelData<NodeFArrayBox>& a_phi,
              /// CELL-centered physical domain on each level
              const ProblemDomain& a_domain,
              /// pointer to CELL-centered grids on the next finer level, or NULL if there is none.
              const DisjointBoxLayout* a_finerGridsPtr,
              /// refinement ratio to next finer level
              const int a_nRefFine,
              /// mesh spacing at this level
              const Real a_dx,
              /// interval of components to use
              const Interval a_comps,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);

/// Computes integral of single-level array, counting only valid nodes.
/**
    Computes integral of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    This function computes exterior boundary nodes, and interior boundary
    nodes of the coarsened finer level, if any.  It is more efficient
    to precompute these boundary-node objects and call one of the other
    integral functions below.
*/
Real integral(/// the data on which to take the norm
              const LevelData<NodeFArrayBox>& a_phi,
              /// CELL-centered physical domain on each level
              const Box& a_domain,
              /// pointer to CELL-centered grids on the next finer level, or NULL if there is none.
              const DisjointBoxLayout* a_finerGridsPtr,
              /// refinement ratio to next finer level
              const int a_nRefFine,
              /// mesh spacing at this level
              const Real a_dx,
              /// interval of components to use
              const Interval a_comps,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);

/// Computes integral of single-level array, counting only valid nodes.
/**
    Computes integral of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>

    The <i>a_IVSVintFinerCoarsened</i> argument is obtained from the call:<br>
    <tt>
    interiorBoundaryNodes(a_IVSVintFinerCoarsened,
    a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real integral(/// the data on which to take the norm
              const LevelData<NodeFArrayBox>& a_phi,
              /// CELL-centered physical domain on each level
              const ProblemDomain& a_domain,
              /// CELL-centered grids on the next finer level
              const DisjointBoxLayout& a_finerGridsCoarsened,
              /// exterior boundary nodes of this level
              const LayoutData< Vector<IntVectSet> >& a_IVSVext,
              /// interior boundary nodes of coarsened finer-level grids
              const LayoutData< Vector<IntVectSet> >& a_IVSVintFinerCoarsened,
              /// refinement ratio to next finer level
              const int a_nRefFine,
              /// mesh spacing at this level
              const Real a_dx,
              /// interval of components to use
              const Interval a_comps,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);

/// Computes integral of single-level array, counting only valid nodes.
/**
    Computes integral of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>

    The <i>a_IVSVintFinerCoarsened</i> argument is obtained from the call:<br>
    <tt>
    interiorBoundaryNodes(a_IVSVintFinerCoarsened,
    a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real integral(/// the data on which to take the norm
              const LevelData<NodeFArrayBox>& a_phi,
              /// CELL-centered physical domain on each level
              const Box& a_domain,
              /// CELL-centered grids on the next finer level
              const DisjointBoxLayout& a_finerGridsCoarsened,
              /// exterior boundary nodes of this level
              const LayoutData< Vector<IntVectSet> >& a_IVSVext,
              /// interior boundary nodes of coarsened finer-level grids
              const LayoutData< Vector<IntVectSet> >& a_IVSVintFinerCoarsened,
              /// refinement ratio to next finer level
              const int a_nRefFine,
              /// mesh spacing at this level
              const Real a_dx,
              /// interval of components to use
              const Interval a_comps,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);


/// Computes integral of single-level array, counting only valid nodes, with no finer level.
/**
    Computes integral of single-level array <i>a_phi</i>, where there is no
    finer level.  Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>
*/
Real integral(/// the data on which to take the norm
              const LevelData<NodeFArrayBox>& a_phi,
              /// CELL-centered physical domain on each level
              const ProblemDomain& a_domain,
              /// exterior boundary nodes of this level
              const LayoutData< Vector<IntVectSet> >& a_IVSVext,
              /// mesh spacing at this level
              const Real a_dx,
              /// interval of components to use
              const Interval a_comps,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);


/// Computes integral of single-level array, counting only valid nodes, with no finer level.
/**
    Computes integral of single-level array <i>a_phi</i>, where there is no
    finer level.  Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>
*/
Real integral(/// the data on which to take the norm
              const LevelData<NodeFArrayBox>& a_phi,
              /// CELL-centered physical domain on each level
              const Box& a_domain,
              /// exterior boundary nodes of this level
              const LayoutData< Vector<IntVectSet> >& a_IVSVext,
              /// mesh spacing at this level
              const Real a_dx,
              /// interval of components to use
              const Interval a_comps,
              /// verbosity, whether to print norms of each box
              bool a_verbose = true);

#include "NamespaceFooter.H"
#endif
