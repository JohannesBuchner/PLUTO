#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeMGInterp.H
// petermc, 20 Jun 2001
#ifndef _NODEMGINTERP_H_
#define _NODEMGINTERP_H_

#include "REAL.H"
#include "LevelData.H"
#include "NodeFArrayBox.H"
#include "NamespaceHeader.H"

/// Class to interpolate coarse-level data to a fine level
/**
   NodeMGInterp manages the interpolation of coarse-level data
   to fine-level data for multigrid.
   This class should be considered internal
   to AMRNodeSolver and not a part of the Chombo API.
*/
class NodeMGInterp
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.  User must subsequently call define().
   */
  NodeMGInterp();

  ///
  /** Destructor.
   */
  ~NodeMGInterp();

  ///
  /** Defining constructor.  Calls define() function with
      identical arguments.
  */
  NodeMGInterp(const DisjointBoxLayout& a_grids,
               int a_numcomps,
               int a_refRatio,
               const Box& a_domain);

  ///
  /** Defining constructor.  Calls define() function with
      identical arguments.
  */
  NodeMGInterp(const DisjointBoxLayout& a_grids,
               int a_numcomps,
               int a_refRatio,
               const ProblemDomain& a_domain);

  ///
  /** Defines this object.  Existing information is overridden.
  */
  void define(/// CELL-centered grids at the finer level
              const DisjointBoxLayout& a_grids,
              /// number of components of data
              int a_numcomps,
              /// refinement ratio between the two levels
              int a_refRatio,
              /// CELL-centered physical domain at the finer level
              const ProblemDomain& a_domain);

  ///
  /** Defines this object.  Existing information is overridden.
  */
  void define(/// CELL-centered grids at the finer level
              const DisjointBoxLayout& a_grids,
              /// number of components of data
              int a_numcomps,
              /// refinement ratio between the two levels
              int a_refRatio,
              /// CELL-centered physical domain at the finer level
              const Box& a_domain);

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
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Adds interpolated coarse-level data to fine-level data.
  */
  void interpToFine(/// data at the finer level
                    LevelData<NodeFArrayBox>& a_fine,
                    /// data at the coarser level
                    const LevelData<NodeFArrayBox>& a_coarse,
                    /// whether the finer-level grids are refinements of the coarser-level grids
                    bool a_sameGrids = false);

  /*@}*/

protected:

  /** has this NodeMGInterp object been defined yet?
   */
  bool is_defined;

  /** the refinement ratio between the two levels
   */
  int m_refRatio;

  /** work array for the coarse-level data in a domain that is the
      outline of the fine level domain on the coarse level
  */
  BoxLayoutData<NodeFArrayBox> m_coarsenedFine;

  /** CELL-centered grids at the finer level
   */
  DisjointBoxLayout m_grids;

  /** CELL-centered physical domain at the finer level
   */
  ProblemDomain m_domain;

  /** reference box used in interpolation
   */
  Box m_boxRef;

  /** weights of corners in interpolation
   */
  FArrayBox m_weights;
};

#include "NamespaceFooter.H"
#endif
