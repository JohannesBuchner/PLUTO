#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Norms.H
// petermc, 12 Aug 2002
// petermc, 17 Oct 2005, put integral functions into NodeIntegrals.H

#ifndef _NODENORMS_H_
#define _NODENORMS_H_

#include "LevelData.H"
#include "NodeFArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

/// Computes L^p norm of a subbox of a NodeFArrayBox.
/** Computes L^p norm of the data in a subbox of <i>a_nfab</i>.
 */
Real norm(/// the data on which to take the norm
          const NodeFArrayBox& a_nfab,
          /// mesh spacing at this level
          const Real a_dx,
          /// CELL-centered subbox over which to take norm
          const Box& a_subbox,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// starting component
          const int a_startComp = 0,
          /// number of components
          const int a_numComp = 1);

/// Computes max norm of a subbox of a NodeFArrayBox.
/** Computes max norm of the data in a subbox of <i>a_nfab</i>.
 */
Real maxnorm(/// the data on which to take the norm
             const NodeFArrayBox& a_nfab,
             /// CELL-centered subbox over which to take norm
             const Box& a_subbox,
             /// starting component
             const int a_startComp = 0,
             /// number of components
             const int a_numComp = 1);

/// Computes L^p norm of a BoxLayoutData<NodeFArrayBox>.
/** Computes L^p norm of <i>a_layout</i>.
    Ignores ghost cells for LevelData.

    This is an internal function that should NOT normally
    be called by the user.
 */
Real norm(/// the data on which to take the norm
          const BoxLayoutData<NodeFArrayBox>& a_layout,
          /// mesh spacing at this level
          const Real a_dx,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// interval of components to use
          const Interval& a_interval,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes L^p norm of a BoxLayoutData<NodeFArrayBox>.
/** Computes L^p norm of <i>a_layout</i>.
    Ignores ghost cells for LevelData.

    This function uses the masks.
 */
Real norm(/// the data on which to take the norm
          const BoxLayoutData<NodeFArrayBox>& a_layout,
          /// mask for the data
          const LevelData<NodeFArrayBox>& a_mask,
          /// domain of this level
          const ProblemDomain& a_domain,
          /// mesh spacing at this level
          const Real a_dx,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// interval of components to use
          const Interval& a_interval,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes max norm of a BoxLayoutData<NodeFArrayBox>.
/** Computes max norm of <i>a_layout</i>.
    Ignores ghost cells for LevelData.

    This function uses the masks.
 */
Real maxnorm(/// the data on which to take the norm
             const BoxLayoutData<NodeFArrayBox>& a_layout,
             /// mask for the data
             const LevelData<NodeFArrayBox>& a_mask,
             /// domain of the data
             const ProblemDomain& a_domain,
             /// interval of components to use
             const Interval& a_interval,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);

/// Computes max norm of a BoxLayoutData<NodeFArrayBox>.
/** Computes max norm of <i>a_layout</i>.
    Ignores ghost cells for LevelData.

    This is an internal function that should NOT normally
    be called by the user.
 */
Real maxnorm(/// the data on which to take the norm
             const BoxLayoutData<NodeFArrayBox>& a_layout,
             /// components to use in computing norm
             const Interval& a_interval,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);

/// Computes norm of multilevel array, counting only valid nodes at each level.
/**
    Computes norm of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.
 */
Real norm(/// multilevel data on which to take the norm
          const Vector<LevelData<NodeFArrayBox>* >& a_phi,
          /// CELL-centered physical domain on each level
          const Vector<ProblemDomain>& a_domain,
          /// refinement ratios between successive levels
          const Vector<int>& a_nRefFine,
          /// mesh spacing at base level <i>a_lBase</i>
          const Real a_dxCrse,
          /// components to use in computing norm
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// index of base level to use in Vectors
          const int a_lBase,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes norm of multilevel array, counting only valid nodes at each level.
/**
    Computes norm of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.

    a_IVSVext and a_IVSVintFinerCoarsened can be found by:
    exteriorAndInteriorNodes(a_IVSVext, a_IVSVintFinerCoarsened,
                             layouts, a_domain, a_nRefFine);
    where layouts is Vector<DisjointBoxLayout> on which the data in
    a_phi live.
 */
Real norm(/// multilevel data on which to take the norm
          const Vector<LevelData<NodeFArrayBox>* >& a_phi,
          /// CELL-centered physical domain on each level
          const Vector<ProblemDomain>& a_domain,
          /// exterior boundary nodes at all level
          const Vector<LayoutData< Vector<Box> >* >& a_IVSVext,
          /// interior boundary nodes of coarsened finer-level grids at all levels except finest
          const Vector<LayoutData< Vector<Box> >* >& a_IVSVintFinerCoarsened,
          /// coarsened finer-level grids at all levels except finest
          const Vector<DisjointBoxLayout>& a_layoutsFinerCoarsened,
          /// refinement ratios between successive levels
          const Vector<int>& a_nRefFine,
          /// mesh spacing at base level <i>a_lBase</i>
          const Real a_dxCrse,
          /// components to use in computing norm
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// index of base level to use in Vectors
          const int a_lBase,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes norm of multilevel array, counting only valid nodes at each level.
/**
    Computes norm of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.

    a_IVSVext and a_IVSVintFinerCoarsened can be found by:
    exteriorAndInteriorNodes(a_IVSVext, a_IVSVintFinerCoarsened,
                             layouts, a_domain, a_nRefFine);
    where layouts is Vector<DisjointBoxLayout> on which the data in
    a_phi live.
 */
Real norm(/// multilevel data on which to take the norm
          const Vector<LevelData<NodeFArrayBox>* >& a_phi,
          /// CELL-centered physical domain on each level
          const Vector<ProblemDomain>& a_domain,
          /// exterior boundary nodes at all level
          const Vector<LayoutData< Vector<IntVectSet> >* >& a_IVSVext,
          /// interior boundary nodes of coarsened finer-level grids at all levels except finest
          const Vector<LayoutData< Vector<IntVectSet> >* >& a_IVSVintFinerCoarsened,
          /// coarsened finer-level grids at all levels except finest
          const Vector<DisjointBoxLayout>& a_layoutsFinerCoarsened,
          /// refinement ratios between successive levels
          const Vector<int>& a_nRefFine,
          /// mesh spacing at base level <i>a_lBase</i>
          const Real a_dxCrse,
          /// components to use in computing norm
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// index of base level to use in Vectors
          const int a_lBase,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes norm of multilevel array, counting only valid nodes at each level.
/**
    Computes norm of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted, with the mask vector.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.
 */
Real norm(/// multilevel data on which to take the norm
          const Vector<LevelData<NodeFArrayBox>* >& a_phi,
          /// mask at each level
          const Vector<LevelData<NodeFArrayBox>* >& a_mask,
          /// domain at each level
          const Vector<ProblemDomain>& a_vectPD,
          /// refinement ratios between successive levels
          const Vector<int>& a_nRefFine,
          /// mesh spacing at base level <i>a_lBase</i>
          const Real a_dxCrse,
          /// components to use in computing norm
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// index of base level to use in Vectors
          const int a_lBase,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes norm of multilevel array, counting only valid nodes at each level.
/**
    Computes norm of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.
 */
Real norm(/// multilevel data on which to take the norm
          const Vector<LevelData<NodeFArrayBox>* >& a_phi,
          /// CELL-centered physical domain on each level
          const Vector<Box>& a_domain,
          /// refinement ratios between successive levels
          const Vector<int>& a_nRefFine,
          /// mesh spacing at base level <i>a_lBase</i>
          const Real a_dxCrse,
          /// components to use in computing norm
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// index of base level to use in Vectors
          const int a_lBase,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes norm of single-level array, counting only valid nodes.
/**
    Computes norm of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    This function computes exterior boundary nodes, and interior boundary
    nodes of the coarsened finer level, if any.  It is more efficient
    to precompute these boundary-node objects and call one of the other
    norm functions below.
 */
Real norm(/// the data on which to take the norm
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
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes norm of single-level array, counting only valid nodes.
/**
    Computes norm of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    This function computes exterior boundary nodes, and interior boundary
    nodes of the coarsened finer level, if any.  It is more efficient
    to precompute these boundary-node objects and call one of the other
    norm functions below.
 */
Real norm(/// the data on which to take the norm
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
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);

/// Computes norm of single-level array, counting only valid nodes.
/**
    Computes norm of single-level array <i>a_phi</i>.
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
    interiorBoundaryNodes(a_IVSVintFinerCoarsened, a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real norm(/// the data on which to take the norm
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
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);


/// Computes norm of single-level array, counting only valid nodes.
/**
    Computes norm of single-level array <i>a_phi</i>.
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
    interiorBoundaryNodes(a_IVSVintFinerCoarsened, a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real norm(/// the data on which to take the norm
          const LevelData<NodeFArrayBox>& a_phi,
          /// CELL-centered physical domain on each level
          const ProblemDomain& a_domain,
          /// CELL-centered grids on the next finer level
          const DisjointBoxLayout& a_finerGridsCoarsened,
          /// exterior boundary nodes of this level
          const LayoutData< Vector<Box> >& a_IVSVext,
          /// interior boundary nodes of coarsened finer-level grids
          const LayoutData< Vector<Box> >& a_IVSVintFinerCoarsened,
          /// refinement ratio to next finer level
          const int a_nRefFine,
          /// mesh spacing at this level
          const Real a_dx,
          /// interval of components to use
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);


/// Computes norm of single-level array, counting only valid nodes.
/**
    Computes norm of single-level array <i>a_phi</i>.
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
    interiorBoundaryNodes(a_IVSVintFinerCoarsened, a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real norm(/// the data on which to take the norm
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
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);


/// Computes norm of single-level array, counting only valid nodes.
/**
    Computes norm of single-level array <i>a_phi</i>.
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
    interiorBoundaryNodes(a_IVSVintFinerCoarsened, a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real norm(/// the data on which to take the norm
          const LevelData<NodeFArrayBox>& a_phi,
          /// CELL-centered physical domain on each level
          const Box& a_domain,
          /// CELL-centered grids on the next finer level
          const DisjointBoxLayout& a_finerGridsCoarsened,
          /// exterior boundary nodes of this level
          const LayoutData< Vector<Box> >& a_IVSVext,
          /// interior boundary nodes of coarsened finer-level grids
          const LayoutData< Vector<Box> >& a_IVSVintFinerCoarsened,
          /// refinement ratio to next finer level
          const int a_nRefFine,
          /// mesh spacing at this level
          const Real a_dx,
          /// interval of components to use
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);


/// Computes norm of single-level array, counting only valid nodes, with no finer level.
/**
    Computes norm of single-level array <i>a_phi</i>, where there is no
    finer level.  Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>
*/
Real norm(/// the data on which to take the norm
          const LevelData<NodeFArrayBox>& a_phi,
          /// CELL-centered physical domain on each level
          const ProblemDomain& a_domain,
          /// exterior boundary nodes of this level
          const LayoutData< Vector<IntVectSet> >& a_IVSVext,
          /// mesh spacing at this level
          const Real a_dx,
          /// interval of components to use
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);


/// Computes norm of single-level array, counting only valid nodes, with no finer level.
/**
    Computes norm of single-level array <i>a_phi</i>, where there is no
    finer level.  Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>
*/
Real norm(/// the data on which to take the norm
          const LevelData<NodeFArrayBox>& a_phi,
          /// CELL-centered physical domain on each level
          const ProblemDomain& a_domain,
          /// exterior boundary nodes of this level
          const LayoutData< Vector<Box> >& a_IVSVext,
          /// mesh spacing at this level
          const Real a_dx,
          /// interval of components to use
          const Interval& a_comps,
          /// norm to take, or 0 for max norm
          const int a_p,
          /// verbosity, whether to print norms of each box
          bool a_verbose = false);


/// Computes max norm of multilevel array, counting only valid nodes at each level.
/**
    Computes max norm of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.
 */
Real maxnorm(/// multilevel data on which to take the norm
             const Vector<LevelData<NodeFArrayBox>* >& a_phi,
             /// CELL-centered physical domain on each level
             const Vector<ProblemDomain>& a_domain,
             /// refinement ratios between successive levels
             const Vector<int>& a_nRefFine,
             /// components to use in computing norm
             const Interval& a_comps,
             /// index of base level to use in Vectors
             const int a_lBase,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);

/// Computes max norm of multilevel array, counting only valid nodes at each level.
/**
    Computes max norm of multilevel array <i>a_phi</i>.
    Only the valid nodes at each level are counted.
    The valid nodes at a level are those that are not covered
    by the interior of any finer level.
 */
Real maxnorm(/// multilevel data on which to take the norm
             const Vector<LevelData<NodeFArrayBox>* >& a_phi,
             /// CELL-centered physical domain on each level
             const Vector<Box>& a_domain,
             /// refinement ratios between successive levels
             const Vector<int>& a_nRefFine,
             /// components to use in computing norm
             const Interval& a_comps,
             /// index of base level to use in Vectors
             const int a_lBase,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);

/// Computes max norm of single-level array, counting only valid nodes.
/**
    Computes max norm of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    This function computes exterior boundary nodes, and interior boundary
    nodes of the coarsened finer level, if any.  It is more efficient
    to precompute these boundary-node objects and call one of the other
    maxnorm functions below.
 */
Real maxnorm(/// the data on which to take the norm
             const LevelData<NodeFArrayBox>& a_phi,
             /// CELL-centered physical domain on each level
             const ProblemDomain& a_domain,
             /// pointer to CELL-centered grids on the next finer level, or NULL if there is none.
             const DisjointBoxLayout* a_finerGridsPtr,
             /// refinement ratio to next finer level
             const int a_nRefFine,
             /// interval of components to use
             const Interval& a_comps,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);

/// Computes max norm of single-level array, counting only valid nodes.
/**
    Computes max norm of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    This function computes exterior boundary nodes, and interior boundary
    nodes of the coarsened finer level, if any.  It is more efficient
    to precompute these boundary-node objects and call one of the other
    maxnorm functions below.
 */
Real maxnorm(/// the data on which to take the norm
             const LevelData<NodeFArrayBox>& a_phi,
             /// CELL-centered physical domain on each level
             const Box& a_domain,
             /// pointer to CELL-centered grids on the next finer level, or NULL if there is none.
             const DisjointBoxLayout* a_finerGridsPtr,
             /// refinement ratio to next finer level
             const int a_nRefFine,
             /// interval of components to use
             const Interval& a_comps,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);


/// Computes max norm of single-level array, counting only valid nodes.
/**
    Computes max norm of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>

    The <i>a_IVSVintFinerCoarsened</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(a_IVSVintFinerCoarsened,
    a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real maxnorm(/// the data on which to take the norm
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
             /// interval of components to use
             const Interval& a_comps,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);


/// Computes max norm of single-level array, counting only valid nodes.
/**
    Computes max norm of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>

    The <i>a_IVSVintFinerCoarsened</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(a_IVSVintFinerCoarsened,
    a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real maxnorm(/// the data on which to take the norm
             const LevelData<NodeFArrayBox>& a_phi,
             /// CELL-centered physical domain on each level
             const ProblemDomain& a_domain,
             /// CELL-centered grids on the next finer level
             const DisjointBoxLayout& a_finerGridsCoarsened,
             /// exterior boundary nodes of this level
             const LayoutData< Vector<Box> >& a_IVSVext,
             /// interior boundary nodes of coarsened finer-level grids
             const LayoutData< Vector<Box> >& a_IVSVintFinerCoarsened,
             /// refinement ratio to next finer level
             const int a_nRefFine,
             /// interval of components to use
             const Interval& a_comps,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);


/// Computes max norm of single-level array, counting only valid nodes.
/**
    Computes max norm of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>

    The <i>a_IVSVintFinerCoarsened</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(a_IVSVintFinerCoarsened,
    a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real maxnorm(/// the data on which to take the norm
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
             /// interval of components to use
             const Interval& a_comps,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);


/// Computes max norm of single-level array, counting only valid nodes.
/**
    Computes max norm of single-level array <i>a_phi</i>.
    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level that are not
    projections of interior nodes of the finer level, if any.

    The <i>a_IVSVext</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>

    The <i>a_IVSVintFinerCoarsened</i> argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(a_IVSVintFinerCoarsened,
    a_phi.getBoxes(), a_finerGridsCoarsened, a_domain);
    </tt>
*/
Real maxnorm(/// the data on which to take the norm
             const LevelData<NodeFArrayBox>& a_phi,
             /// CELL-centered physical domain on each level
             const Box& a_domain,
             /// CELL-centered grids on the next finer level
             const DisjointBoxLayout& a_finerGridsCoarsened,
             /// exterior boundary nodes of this level
             const LayoutData< Vector<Box> >& a_IVSVext,
             /// interior boundary nodes of coarsened finer-level grids
             const LayoutData< Vector<Box> >& a_IVSVintFinerCoarsened,
             /// refinement ratio to next finer level
             const int a_nRefFine,
             /// interval of components to use
             const Interval& a_comps,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);


/// Computes max norm of single-level array, counting only valid nodes, with no finer level.
/**
    Computes max norm of single-level array <i>a_phi</i>, where there is no
    finer level.

    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level.

    The a_IVSVext argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>
*/
Real maxnorm(/// the data on which to take the norm
             const LevelData<NodeFArrayBox>& a_phi,
             /// CELL-centered physical domain on each level
             const ProblemDomain& a_domain,
             /// exterior boundary nodes of this level
             const LayoutData< Vector<IntVectSet> >& a_IVSVext,
             /// interval of components to use
             const Interval& a_comps,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);

/// Computes max norm of single-level array, counting only valid nodes, with no finer level.
/**
    Computes max norm of single-level array <i>a_phi</i>, where there is no
    finer level.

    Only the valid nodes at this level are counted.
    The valid nodes are the interior nodes of this level.

    The a_IVSVext argument is obtained from the calls:<br>
    <tt>
    interiorBoundaryNodes(IVSVint, a_phi.getBoxes(), a_domain);<br>
    exteriorBoundaryNodes(a_IVSVext, IVSVint, a_phi.getBoxes();
    </tt>
*/
Real maxnorm(/// the data on which to take the norm
             const LevelData<NodeFArrayBox>& a_phi,
             /// CELL-centered physical domain on each level
             const Box& a_domain,
             /// exterior boundary nodes of this level
             const LayoutData< Vector<Box> >& a_IVSVext,
             /// interval of components to use
             const Interval& a_comps,
             /// verbosity, whether to print norms of each box
             bool a_verbose = false);

#include "NamespaceFooter.H"
#endif
