#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeQuadCFInterp2.H
// original NodeQuadCFInterp by petermc, Mon, Apr 23, 2001
// petermc, 31 Oct 2001
#ifndef _NODEQUADCFINTERP2_H_
#define _NODEQUADCFINTERP2_H_

#include <cmath>
#include <cstdlib>
#include "BaseFab.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "NodeCFIVS.H"
#include "NodeBC.H"
#include "NamespaceHeader.H"


///  Class to interpolate quadratically at coarse/fine interface when refinement ratio is 2
class NodeQuadCFInterp2
/** Class to interpolate quadratically at interface between this level
    and next coarser level, when the refinement ratio is 2.
    This class should be considered internal to NodeQuadCFInterp.

    <b>Long Description:</b>

    The interface has codimension one.<br>
    If a fine node coincides with a coarse node, then we merely project
    from the coarse node to the fine node.<br>
    Otherwise, we take the mean of the neighboring coarse nodes
    and subtract  1/8 * dx^2  times the sum of the second derivatives
    in the appropriate directions.

    The interpolation is performed in function coarseFineInterp().
    The constructor computes <i>m_loCFIVS</i> and <i>m_hiCFIVS</i> to
    determine the fine nodes at the interface with the coarse level.

    The constructor also takes <i>m_loCFIVS</i> and <i>m_hiCFIVS</i> to
    determine the fine nodes at the interface with the coarse level.
    Calling getFineIVS() on <i>m_loCFIVS</i>[idir][dit()] gives us the
    IntVectSet of nodes of <i>m_grids</i>[dit()] on the face in the low
    direction in dimension idir, that lie on the interface with the
    coarser level.  Similarly with <i>m_hiCFIVS</i>[idir][dit()] for the
    opposite face, in the high direction in dimension idir.

    <b>2-D Description:</b>

    In the 2-D problem, the interface is 1-D.  Between coarse nodes
    at 0 and 1, we approximate the value at the fine node by<br>
    f(1/2) ~ (f(0) + f(1))/2 - 1/8 * f''(1/2)<br>
    where we estimate the second derivative f''(1/2) from the
    coarse values f(0), f(1), and either f(-1) or f(2) or both.

    If the points -1 and 2 are both on the grid:
<pre>
      o---o-x-o---o
     -1   0   1   2
</pre>
    then we use<br>
    f''(1/2) ~ (f(-1) - f(0) - f(1) + f(2))/2.

    If the point -1 is on the grid but 2 is not:
<pre>
      o---o-x-o   o
     -1   0   1   2
</pre>
    then we approximate f''(1/2) by f''(0) and use<br>
    f''(0) ~ (f(-1) - 2 * f(0) + f(1)).

    If the point 2 is on the grid but -1 is not:
<pre>
      o   o-x-o---o
     -1   0   1   2
</pre>
    then we approximate f''(1/2) by f''(1) and use<br>
    f''(1) ~ (f(0) - 2 * f(1) + f(2)).

    <b>3-D Description:</b>

    In the 3-D problem, the interface is 2-D.  For any given fine node
    along the interface, look at its three coordinates.  At least one
    of the coordinates is divisible by 2, meaning that it lies on a
    plane of coarse nodes.  If all of the coordinates are even, then
    the fine node coincides with a coarse node, and we project the
    value.  If only one coordinate is odd, then this reduces to the
    problem of a 1-D interface described above in the 2-D case.

    We are left with the problem of interpolating f(1/2, 1/2).
<pre>
    (0,1)   (1,1)
      o-------o
      |       |
      |   x   |
      |       |
      o-------o
    (0,0)   (1,0)
</pre>
    We use<br>
    f(1/2,1/2) ~ (f(0,0) + f(0,1) + f(1,0) + f(1,1))/4 -
    1/8 * ( d^2 f/dx^2 (1/2,1/2) + d^2 f/dy^2 (1/2,1/2))<br>
    where the second derivatives are estimated from the four
    neighboring coarse nodes and their neighbors.

    In particular, d^2 f/dx^2 (1/2,1/2) is approximated by the
    mean of d^2 f/dx^2 (1/2, 0) and d^2 f/dx^2 (1/2, 1).<br>
    These second derivatives are estimated in the same way as
    described above for the 1-D interface in the 2-D case.
*/
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.  User must subsequently call define().
   */
  NodeQuadCFInterp2();

  ///
  /** Constructor calls setDefaultValues() and then calls
      define() with the same arguments.
  */
  NodeQuadCFInterp2(const DisjointBoxLayout& a_grids,
                    const Box& a_domain,
                    const LayoutData<NodeCFIVS>* const a_loCFIVS,
                    const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                    bool a_interfaceOnly,
                    int a_interpolationDegree,
                    int a_ncomp = 1,
                    bool a_verbose = false);

  ///
  /** Constructor calls setDefaultValues() and then calls
      define() with the same arguments.
  */
  NodeQuadCFInterp2(const DisjointBoxLayout& a_grids,
                    const ProblemDomain& a_domain,
                    const LayoutData<NodeCFIVS>* const a_loCFIVS,
                    const LayoutData<NodeCFIVS>* const a_hiCFIVS,
                    bool a_interfaceOnly,
                    int a_interpolationDegree,
                    int a_ncomp = 1,
                    bool a_verbose = false);

  ///
  /** Destructor.
   */
  ~NodeQuadCFInterp2();

  ///
  /** Full define function.  Makes all coarse-fine
      information and sets internal variables.  The current level
      is taken to be the fine level.
  */
  void define(/// CELL-centered grids at this level
              const DisjointBoxLayout& a_grids,
              /// CELL-centered physical domain at this level
              const ProblemDomain& a_domain,
              /// pointer to object storing coarse/fine interface nodes
              const LayoutData<NodeCFIVS>* const a_loCFIVS,
              /// pointer to object storing coarse/fine interface nodes
              const LayoutData<NodeCFIVS>* const a_hiCFIVS,
              /// whether interpolation is from interface only, meaning that off-interface data should not be used
              bool a_interfaceOnly,
              /// degree of interpolation; 1 for (bi)linear, 2 for (bi)quadratic
              int a_interpolationDegree,
              /// number of components of data
              int a_ncomp = 1,
              /// verbose output flag
              bool a_verbose = false);

  ///
  /** Full define function.  Makes all coarse-fine
      information and sets internal variables.  The current level
      is taken to be the fine level.
  */
  void define(/// CELL-centered grids at this level
              const DisjointBoxLayout& a_grids,
              /// CELL-centered physical domain at this level
              const Box& a_domain,
              /// pointer to object storing coarse/fine interface nodes
              const LayoutData<NodeCFIVS>* const a_loCFIVS,
              /// pointer to object storing coarse/fine interface nodes
              const LayoutData<NodeCFIVS>* const a_hiCFIVS,
              /// whether interpolation is from interface only, meaning that off-interface data should not be used
              bool a_interfaceOnly,
              /// degree of interpolation; 1 for (bi)linear, 2 for (bi)quadratic
              int a_interpolationDegree,
              /// number of components of data
              int a_ncomp = 1,
              /// verbose output flag
              bool a_verbose = false);


  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object was created with the defining
      constructor or if define() has been called.
  */
  bool isDefined() const;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Set whether to give output.  Default is <tt>false</tt>.
   */
  void setVerbose( bool a_verbose );

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Coarse / Fine (inhomogeneous) interpolation operator.
      Fill the nodes of <i>a_phi</i> on the coarse/fine interface with
      interpolated data from <i>a_phiCoarse</i>.
  */
  void coarseFineInterp(/// data at this level
                        LevelData<NodeFArrayBox>& a_phi,
                        /// data at the next coarser level
                        const LevelData<NodeFArrayBox>& a_phiCoarse);

  /*@}*/

protected:

  /** CELL-centered grids at the current level (the finer level)
   */
  DisjointBoxLayout m_grids;

  /** CELL-centered physical domain of this level
   */
  ProblemDomain m_domain;

  /** NODE-centered box of nodes at physical domain at the coarser level
   */
  Box m_domainCoarseNodes;

  /** CELL-centered <i>m_grids</i> coarsened by 2
   */
  DisjointBoxLayout m_coarsenedGrids;

  /** copy of coarse phi, used in CFInterp().
   */
  LevelData<NodeFArrayBox> m_coarseCopy;

  /** interpolating from interface only?
   */
  bool m_interfaceOnly;

  /** degree of interpolation:  1 for (bi)linear, 2 for (bi)quadratic
   */
  int m_interpolationDegree;

  /** number of components of data, needed for setting size of work array
   */
  int m_ncomp;

  /** has full define function been called?
   */
  bool m_isDefined;

  /** pointer to object storing coarse/fine interface nodes
      between this level and next coarser level
  */
  const LayoutData<NodeCFIVS>* m_loCFIVS;

  /** pointer to object storing coarse/fine interface nodes
      between this level and next coarser level
  */
  const LayoutData<NodeCFIVS>* m_hiCFIVS;
  // LayoutData<NodeCFIVS> m_loCFIVS[SpaceDim];
  // LayoutData<NodeCFIVS> m_hiCFIVS[SpaceDim];

  // weights for coarse/fine interpolation
  // Tuple<FArrayBox*, SpaceDim> m_wtLo, m_wtHi, m_wtC;

  /** if <tt>true</tt>, print out extra information
   */
  bool m_verbose;

private:
  //internally useful functions

  void clearMemory();

  void setDefaultValues();

  /** Interpolate from <i>a_psiFab</i> to <i>a_fineFab</i>
      at <i>a_iv</i> along a line.
   */
  void interpLine(/// fine data on NODE-centered FAB
                  FArrayBox& a_fineFab,
                  /// coarse data on NODE-centered FAB
                  const FArrayBox& a_psiFab,
                  /// nodes of <i>a_psiFab</i> from which data can be used
                  const IntVectSet& a_psiFabNodes,
                  /// point at which to interpolate
                  const IntVect& a_iv,
                  /// direction (0 to SpaceDim-1) along which to interpolate
                  int a_idirOther);

};

#include "NamespaceFooter.H"
#endif
