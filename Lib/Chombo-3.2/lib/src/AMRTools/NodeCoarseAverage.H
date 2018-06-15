#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NODECOARSEAVERAGE_H_
#define _NODECOARSEAVERAGE_H_

#include "REAL.H"
#include "BaseFab.H"
#include "IntVectSet.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "BitSet.H"
#include "NamespaceHeader.H"

/// Replaces coarse-level data with an average of fine-level data

class NodeCoarseAverage
/**
   This class replaces data at a coarse level of refinement with an
   average of data at a finer level of refinement, in areas where fine
   data is present.  Coarse-level data is not modified where fine-level
   data is not present.
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
  NodeCoarseAverage();

  ///
  /** Defining constructor.  Calls define() function with
      identical arguments.
  */
  NodeCoarseAverage(const DisjointBoxLayout& a_gridsFine,
                    const DisjointBoxLayout& a_gridsCoarse,
                    int                      a_numcomps,
                    int                      a_refRatio,
                    const ProblemDomain&     a_domainFine);

  ///
  /** Defining constructor.  Calls define() function with
      identical arguments.
  */
  NodeCoarseAverage(const DisjointBoxLayout& a_gridsFine,
                    const DisjointBoxLayout& a_gridsCoarse,
                    int                      a_numcomps,
                    int                      a_refRatio,
                    const Box&               a_domainFine);

  ///
  /** Defining constructor.  Calls define() function with
      identical arguments.
  */
  NodeCoarseAverage(const DisjointBoxLayout& a_gridsCoarse,
                    int                      a_numcomps,
                    int                      a_refRatio,
                    const ProblemDomain&     a_domainFine);

  ///
  /** Defining constructor.  Calls define() function with
      identical arguments.
  */
  NodeCoarseAverage(const DisjointBoxLayout& a_gridsCoarse,
                    int                      a_numcomps,
                    int                      a_refRatio,
                    const Box&               a_domainFine);

  ///
  /** Destructor.
   */
  ~NodeCoarseAverage();

  ///
  /** Defines this object.  Existing information is overriden.
  */
  void define(/// CELL-centered grids at the fine level
              const DisjointBoxLayout& a_gridsFine,
              /// CELL-centered grids at the coarse level
              const DisjointBoxLayout& a_gridsCoarse,
              /// number of components of data
              int                      a_numcomps,
              /// refinement ratio between the two levels
              int                      a_refRatio,
              /// CELL-centered physical domain at the fine level
              const ProblemDomain&     a_domainFine);

  ///
  /** Defines this object.  Existing information is overriden.
  */
  void define(/// CELL-centered grids at the fine level
              const DisjointBoxLayout& a_gridsFine,
              /// CELL-centered grids at the coarse level
              const DisjointBoxLayout& a_gridsCoarse,
              /// number of components of data
              int                      a_numcomps,
              /// refinement ratio between the two levels
              int                      a_refRatio,
              /// CELL-centered physical domain at the fine level
              const Box&               a_domainFine);


  ///
  /** Defines this object.  Existing information is overriden.
      Fine grids are refinements of coarse grids.
  */
  void define(/// CELL-centered grids at the coarse level
              const DisjointBoxLayout& a_gridsCoarse,
              /// number of components of data
              int                      a_numcomps,
              /// refinement ratio between the two levels
              int                      a_refRatio,
              /// CELL-centered physical domain at the fine level
              const ProblemDomain&     a_domainFine);

  ///
  /** Defines this object.  Existing information is overriden.
      Fine grids are refinements of coarse grids.
  */
  void define(/// CELL-centered grids at the coarse level
              const DisjointBoxLayout& a_gridsCoarse,
              /// number of components of data
              int                      a_numcomps,
              /// refinement ratio between the two levels
              int                      a_refRatio,
              /// CELL-centered physical domain at the fine level
              const Box&               a_domainFine);

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

  ///
  /** Set whether the solver gives output while solving.
      Default is <tt>true</tt>.
   */
  void setVerbose( bool a_verbose );


  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Replaces coarse data with the average of nearby cells' fine data
      on the interior nodes of the coarsened fine-level domain.
      Elsewhere, the coarse data is unchanged.

      It is an error to call if not this->isDefined().

      The domain of <i>a_fine</i> should be the same as <i>a_gridsFine</i>
      specified in the most recent call to define().

      It is expected that the coarse and fine levels'
      domains are properly nested.  Both <i>a_coarse</i> and <i>a_fine</i>
      should have the same number of components specified in
      the most recent call to define().

      The fine data is modified only in exchange().
      A work array is also modified.
  */
  //     this method would like to be const, but the work array is changed.
  //     this suggests that the work array should not be persistent.
  void averageToCoarse(/// coarse data
                       LevelData<NodeFArrayBox>& a_coarse,
                       /// fine data, modified in exchange()
                       LevelData<NodeFArrayBox>& a_fine);

  /*@}*/

protected:

  /** has this NodeCoarseAverage object been defined yet?
   */
  bool is_defined;

  /** flag for verbose output
   */
  bool m_verbose;

  /** the refinement ratio between the two levels
   */
  int m_refRatio;

  /** the number of components of data, used in work array
   */
  int m_numcomps;

  /** CELL-centered physical domain at the coarser level
   */
  ProblemDomain m_domainCoarse;

  /** box with indices from -<i>m_refRatio</i>/2 to +<i>m_refRatio</i>/2
      in all directions
  */
  Box m_refbox;

  /** weights for coarsening
   */
  FArrayBox m_weights;

  /** are the fine grids merely refined versions of the coarse grids?
   */
  bool m_sameGrids;

  /** CELL-centered grids at the finer level, coarsened by <i>m_refRatio</i>
   */
  DisjointBoxLayout m_coarsenedGrids;

  /** work array for the coarsening of the fine data, based on the layout
      <i>m_coarsenedGrids</i> */
  LevelData<NodeFArrayBox> m_coarsenedFine;

  /** interior boundary nodes of coarsened finer-level grids,
      distributed among the coarser-level grids  */
  LayoutData< Vector<IntVectSet> > m_IVSV;

  /** whether each section of <i>m_IVSVsame</i> is a complete box
   */
  LayoutData< BitSet > m_IVSVfull;

  /** interior boundary nodes of coarsened finer-level grids,
      distributed among the same grids */
  LayoutData< Vector<IntVectSet> > m_IVSVsame;
};

#include "NamespaceFooter.H"
#endif
