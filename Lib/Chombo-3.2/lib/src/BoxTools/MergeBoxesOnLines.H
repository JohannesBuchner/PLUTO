#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MERGEBOXESONLINES_H_
#define _MERGEBOXESONLINES_H_

#include <climits>

#include "Vector.H"
#include "Box.H"
#include "IntVectSet.H"
#include "REAL.H"
#include "Interval.H"
#include <list>
#include "NamespaceHeader.H"

/// Class which merges boxes in a direction without changing the overall box coverage
/**

*/
class MergeBoxesOnLines
{
public:
  MergeBoxesOnLines();

  /// Destructor
  ~MergeBoxesOnLines();

  /// Merges the boxes along lines in the set direction for a single level
  /**
  */
  void
  mergeBoxes(Vector<Box>&      a_mesh,
             const int&        a_dir);

  /// Merges the boxes along lines in the set direction for a vector of levels
  /**
  */
  void
  mergeBoxes(Vector<Vector<Box> >& a_vvbox,
             const int&            a_dir);

};

#include "NamespaceFooter.H"
#endif
