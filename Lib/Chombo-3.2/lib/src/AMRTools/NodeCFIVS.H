#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeCFIVS.H
// petermc, 1 May 2001
#ifndef _NODECFIVS_H_
#define _NODECFIVS_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>
#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "LoHiSide.H"
#include "DisjointBoxLayout.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

/// Class to get fine IntVectSet at coarse/fine interface
/**
    Class to get IntVectSet of nodes (of a particular face of a
    particular box) that lie on interface with next coarser level.

    This class should be considered internal to AMRNodeSolver and should
    not be considered part of the Chombo API.
*/
class NodeCFIVS
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.  User must subsequently call define().
   */
  NodeCFIVS();

  ///
  /** Destructor.
   */
  ~NodeCFIVS();

  ///
  /** Constructor calls setDefaultValues() and then
      calls define() with the same arguments.
  */
  NodeCFIVS(const Box& a_domain,
            const Box& a_box,
            const DisjointBoxLayout& a_levelBoxes,
            int a_idir,
            Side::LoHiSide a_hiorlo);

  ///
  /** Constructor calls setDefaultValues() and then
      calls define() with the same arguments.
  */
  NodeCFIVS(const ProblemDomain& a_domain,
            const Box& a_box,
            const DisjointBoxLayout& a_levelBoxes,
            int a_idir,
            Side::LoHiSide a_hiorlo);

  ///
  /** Full define function.  The current level is taken to be the fine level.
  */
  void define(/// CELL-centered physical domain at this level
              const ProblemDomain& a_domain,
              /// the CELL-centered box at this level
              const Box& a_box,
              /// the set of all CELL-centered boxes at this level
              const DisjointBoxLayout& a_levelBoxes,
              /// direction (0 to SpaceDim-1) of face of box on which we find interface nodes
              int a_idir,
              /// high or low side of box
              Side::LoHiSide a_hiorlo);

  ///
  /** Full define function.  The current level is taken to be the fine level.
  */
  void define(/// CELL-centered physical domain at this level
              const Box& a_domain,
              /// the CELL-centered box at this level
              const Box& a_box,
              /// the set of all CELL-centered boxes at this level
              const DisjointBoxLayout& a_levelBoxes,
              /// direction (0 to SpaceDim-1) of face of box on which we find interface nodes
              int a_idir,
              /// high or low side of box
              Side::LoHiSide a_hiorlo);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  bool isDefined() const;

  ///
  /** Returns <tt>true</tt> if this coarse-fine IntVectSet is empty.
  */
  bool isEmpty() const;

  ///
  /** Returns <tt>true</tt> if this coarse-fine IntVectSet can be represented
      as just a Box, i.e., a user can perform a dense computation instead of
      a pointwise calculation using IVSIterator.

      If isPacked() then you can use packedBox() instead of getFineIVS().
  */
  bool isPacked() const;

  ///
  /** If isPacked() returns <tt>true</tt>,
      then the box returned by this function
      is a suitable substitute for this IntVectSet returned by getFineIVS().
  */
  const Box& packedBox() const;

  ///
  /** Returns indices of fine nodes, on the face of the box,
      that lie on the interface with the coarser level and where
      data need to be interpolated.

      This will be empty if isEmpty() returns <tt>true</tt>.
      If isPacked() then you can use packedBox() instead of getFineIVS().
  */
  const IntVectSet& getFineIVS() const ;

  /*@}*/

  static long long m_packCount, m_sparseCount;

protected:

  //for internal use
  void setDefaultValues();

  /** indices of fine nodes that need to be interpolated
   */
  IntVectSet m_fineInterpIVS;

  bool m_packed;

  bool m_empty;

  Box m_packedBox;

  bool m_isdefined;

private:
  void operator= (const NodeCFIVS& levcfs_in)
  {
  }

  NodeCFIVS(const NodeCFIVS& levcfs_in)
  {
  }
};


inline
bool NodeCFIVS::isEmpty() const
{
  return m_empty;
}


inline
bool NodeCFIVS::isPacked() const
{
  return m_packed;
}


inline
const Box& NodeCFIVS::packedBox() const
{
  return m_packedBox;
}

#include "NamespaceFooter.H"
#endif
