#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FOURTHORDERPATCHINTERP_H_
#define _FOURTHORDERPATCHINTERP_H_

#include "FArrayBox.H"
#include "ProblemDomain.H"
#include "FourthOrderInterpStencil.H"
#include "IntVectSet.H"
#include "UsingNamespace.H"

/// Fourth-order interpolation in space on a single patch

/**
 */
class FourthOrderPatchInterp
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  FourthOrderPatchInterp();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  virtual ~FourthOrderPatchInterp();

  /// Actual constructor.
  /**
     Set up object.
   */
  virtual void define(
                      /// problem domain on this level
                      const ProblemDomain&      a_domain,
                      /// refinement ratio between this level and next coarser level
                      const int&                a_refineCoarse,
                      /// maximum distance of stencil from domain boundary
                      const int&        a_maxStencilDist,
                      /// dimensions that are fixed, not interpolated
                      Interval               a_fixedDims = Interval() );

  /// Set the coarse box before calling setStencil or interpToFine.
  /**
   */
  virtual void setCoarseBox(const Box& a_coarseBox);

  /// Set stencil for current coarse box and domain
  /**
      Set a_stencil for m_coarseBox.
   */
  virtual void setStencil(BaseFab<IntVect>&  a_stencil);

  /// Interpolate in space.
  /**
     Interpolate in space to a_fine from a_coarse.
   */
  virtual void interpToFine(/// interpolated solution on this level
                            FArrayBox&                a_fine,
                            /// coarse solution
                            const FArrayBox&          a_coarse,
                            /// stencils
                            const BaseFab<IntVect>&   a_stencils);

  /// Interpolate in space.
  /**
     Interpolate in space to a_fine from a_coarse at
     fine cells within the coarse cells listed in a_ivs.
   */
  virtual void interpToFine(/// interpolated solution on this level
                            FArrayBox&                a_fine,
                            /// coarse solution
                            const FArrayBox&          a_coarse,
                            /// stencils
                            const BaseFab<IntVect>&   a_stencils,
                            /// we fill in fine cells within these coarse cells
                            const IntVectSet&         a_ivs);

protected:

  /// all possible stencils, on (-m_maxStencilDist:+m_maxStencilDist)^SpaceDim
  BaseFab<FourthOrderInterpStencil*> m_stencils;

  /// Problem domain - index space for this level
  ProblemDomain m_domain;

  // m_domain coarsened by m_refineCoarse
  ProblemDomain m_coarseDomain;

  /// Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  /// maximum distance of stencil from domain boundary
  int m_maxStencilDist;

  /// max degree of approximations
  int m_degree;

  /// dimensions that are fixed, not interpolated
  Interval m_fixedDims;

  /// 1 in m_fixedDims, m_refineCoarse in other dimensions
  IntVect m_refineVect;

  /// define() has been called
  bool m_defined;

  /// current coarse box
  Box m_coarseBox;

  /// m_coarseBox is set
  bool m_isCoarseBoxSet;

private:

  // Disallowed for all the usual reasons
  void operator=(const FourthOrderPatchInterp& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  FourthOrderPatchInterp(const FourthOrderPatchInterp& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
