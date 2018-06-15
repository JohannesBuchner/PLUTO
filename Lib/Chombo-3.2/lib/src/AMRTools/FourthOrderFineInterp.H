#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERFINEINTERP_H_
#define _FOURTHORDERFINEINTERP_H_

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "FourthOrderPatchInterp.H"
#include "UsingNamespace.H"

/// Fourth-order interpolation in space

/**
 */
class FourthOrderFineInterp
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  FourthOrderFineInterp();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~FourthOrderFineInterp();

  /// Actual constructor.
  /**
     Set up object.
   */
  void define(/// layout at this level
              const DisjointBoxLayout&  a_layout,
              /// number of variables
              const int&                a_numStates,
              /// refinement ratio between this level and next coarser level
              const int&                a_refineCoarse,
              /// problem domain on this level
              const ProblemDomain&      a_domain,
              /// number of coarse ghost cells
              const int&                a_coarseGhostsFill = 0,
              /// dimensions that are fixed, not interpolated
              Interval                  a_fixedDims = Interval() );

  /// return stencil indices
  /**
     Get stencil indices for patch a_dit.
   */
  const BaseFab<IntVect>& getStencil(const DataIterator&    a_dit) const;

  /// return stencil indices
  /**
     Get stencil indices for patch a_dind.
   */
  const BaseFab<IntVect>& getStencil(const DataIndex&    a_dind) const;

  /// return reference to coarse data
  /**
     Return reference to m_coarseData.
   */
  LevelData<FArrayBox>& coarsenedFineData();

  /// return patch interpolator
  /**
     Get object for interpolating on a single patch.
   */
  FourthOrderPatchInterp& patchInterp();

  /// Interpolate in space.
  /**
     Interpolate in space to a_fine from a_coarse.
   */
  void interpToFine(/// interpolated solution on this level
                    LevelData<FArrayBox>&         a_fine,
                    /// coarse solution
                    const LevelData<FArrayBox>&   a_coarse);

  /// Interpolate in space on a single patch.
  /**
     Interpolate in space to a_fine from a_coarse, where both
     are related to index a_dit.
   */
  void interpOnPatch(
                     /// interpolated solution on this level
                     FArrayBox&            a_fine,
                     /// coarse solution
                     const FArrayBox&      a_coarse,
                     /// index
                     const DataIterator&   a_dit);

  /// Interpolate in space on a single patch.
  /**
     Interpolate in space to a_fine from a_coarse, where both
     are related to index a_dind.
   */
  void interpOnPatch(
                     /// interpolated solution on this level
                     FArrayBox&            a_fine,
                     /// coarse solution
                     const FArrayBox&      a_coarse,
                     /// index
                     const DataIndex&      a_dind);


  /// Interpolate in space on a single patch, to a set of cells.
  /**
     Interpolate in space to a_fine from a_coarse, where both
     are related to index a_dit, to the fine cells within a_ivs.
   */
  void interpOnPatch(
                     /// interpolated solution on this level
                     FArrayBox&            a_fine,
                     /// coarse solution
                     const FArrayBox&      a_coarse,
                     /// index
                     const DataIterator&   a_dit,
                     /// we fill in fine cells within these coarse cells
                     const IntVectSet&     a_ivs);

  /// Interpolate in space on a single patch.
  /**
     Interpolate in space to a_fine from a_coarse, where both
     are related to index a_dind, to the fine cells within a_ivs.
  */
  void interpOnPatch(
                     /// interpolated solution on this level
                     FArrayBox&            a_fine,
                     /// coarse solution
                     const FArrayBox&      a_coarse,
                     /// index
                     const DataIndex&      a_dind,
                     /// we fill in fine cells within these coarse cells
                     const IntVectSet&     a_ivs);

protected:

  /// Box layout for this level
  DisjointBoxLayout m_layout;

  /// Box layout for the coarse level:  m_layout coarsened by m_refineCoarse
  DisjointBoxLayout m_layoutCoarsened;

  /// object for patch-based operations
  FourthOrderPatchInterp m_patchInterp;

  /// on m_layoutCoarsened; which stencil to use at each coarse point
  LevelData< BaseFab<IntVect> > m_stencilHere;

  /// ghost vector on m_coarseData
  IntVect m_ghostVect;

  /// copied coarse data in interpToFine, lives on m_layoutCoarsened
  LevelData<FArrayBox> m_coarseData;

  /// Problem domain - index space for this level
  ProblemDomain m_domain;

  /// Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  /// dimensions that are fixed, not interpolated
  Interval m_fixedDims;

  /// 1 in m_fixedDims, m_refineCoarse in other dimensions
  IntVect m_refineVect;

  /// Number of variables
  int m_numStates;

  /// maximum distance of stencil from domain boundary
  int m_maxStencilDist;

  /// number of layers of ghost cells at the coarser level whose fine cells are to be filled in
  int m_coarseGhostsFill;

  /// define() has been called
  bool m_defined;

private:

  // Disallowed for all the usual reasons
  void operator=(const FourthOrderFineInterp&);
  FourthOrderFineInterp(const FourthOrderFineInterp&);
};

#endif
