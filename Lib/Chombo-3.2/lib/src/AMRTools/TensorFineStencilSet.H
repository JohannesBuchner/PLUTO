#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TENSORFINESTENCILSET_H_
#define _TENSORFINESTENCILSET_H_

#include "SPACE.H"
#include <stdlib.h>
#include "IntVectSet.H"
#include "ProblemDomain.H"
#include "QuadCFStencil.H"
#include "NamespaceHeader.H"

/// Class to encapsulate fine-level tangential gradient stencil computation
/** This is a helper class which essentially wraps three IntVectSets which
    define which type of stencil to use in computing fine-level tangential
    gradients.  This is only relevant in the presence of a physical boundary
    (since the only place where the fine-level stencil is shifted is at a
    physical boundary; in a periodic world, all differences are centered.)
*/
class TensorFineStencilSet
{
public:
  /// constructors, destructors, defines

  /// default constructor -- leaves class in undefined state
  TensorFineStencilSet();

  /// full constructor -- defines stencil IntVectSets
  TensorFineStencilSet(const IntVectSet& a_fineIVS,
                       const ProblemDomain& a_domain,
                       int a_normalDir);

  /// destructor
  ~TensorFineStencilSet();

  /// define function
  void define(const IntVectSet& a_fineIVS,
              const ProblemDomain& a_domain,
              int a_dir);

  /// accessor
  /** accessor -- returns cells where fine-level tangential gradient
      in a_tanDir direction will use a centered-difference stencil.  If
      a_tanDir is normal to the face (i.e, isn't a tangential direction
      after all), will return an empty IntVectSet.
  */
  const IntVectSet& getCenteredStencilSet(int a_tanDir) const;

  /// accessor
  /** accessor -- returns cells where fine-level tangential gradient
      in a_tanDir will use a forward-difference stencil.  If a_tanDir
      is normal to the face (i.e. it isn't a tangential direction
      after all), will return an empty IntVectSet.
  */
  const IntVectSet& getForwardStencilSet(int a_tanDir) const;

  /// accessor
  /** accessor -- returns cells where fine-level tangential gradient
      in a_tanDir will use a backward-difference stencil.  If a_tanDir
      is normal to the face, will return an empty IntVectSet.
  */
  const IntVectSet& getBackwardStencilSet(int a_tanDir) const;

  /// has this object been defined?
  bool isDefined() const;

  /// does this contain any non-empty IntVectSets?
  bool isEmpty() const;

protected:
  ///
  bool m_isDefined;

  ///
  bool m_isEmpty;

  ///
  int m_normalDir;

#if CH_SPACEDIM == 1
  ///
  int m_tangentialDirections[SpaceDim];
#else
  ///
  int m_tangentialDirections[SpaceDim-1];
#endif

  /// cells which will use centered diffs
#if CH_SPACEDIM == 1
  IntVectSet m_centeredSet[SpaceDim];
#else
  IntVectSet m_centeredSet[SpaceDim-1];
#endif

#if CH_SPACEDIM == 1
  /// cells which will use forward diffs
  IntVectSet m_forwardSet[SpaceDim];
#else
  IntVectSet m_forwardSet[SpaceDim-1];
#endif

  /// cells which will use backward diffs
#if CH_SPACEDIM == 1
  IntVectSet m_backwardSet[SpaceDim];
#else
  IntVectSet m_backwardSet[SpaceDim-1];
#endif

};

#include "NamespaceFooter.H"
#endif

