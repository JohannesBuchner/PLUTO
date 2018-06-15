#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _QUADCFSTENCIL_H_
#define _QUADCFSTENCIL_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "IntVectSet.H"
#include "CFStencil.H"
#include "DerivStencil.H"
#include "NamespaceHeader.H"

/// class to encapsulate CF info for quadratic interpolation
/**
    QuadCFStencil is a wrapper around a bunch of things
    that allow LevelOp to do coarse-fine interpolation
    on a particular face of a grid.
  */
class QuadCFStencil
{

public:
  /// {\bf Derivative functions}

  ///
  /**
     compute second derivative in devdir_a direction
     at coarse point a_ivin
     Uses centered finite diff approximation if
     ivs_standard_m.contains(a_ivin ).
     Otherwise, it uses the stencil from second_m
     Asserts that a_derivdir != direction of face
  */
  Real computeSecondDerivative(
                               const BaseFab<Real> & a_phic,
                               int a_derivdir,
                               int a_ivar,
                               const IntVect& a_ivin,
                               Real a_dx) const;

  ///
  /**
     compute first derivative in devdir_a direction
     at coarse point a_ivin
     Uses centered finite diff approximation if
     ivs_standard_m.contains(a_ivin ).
     Otherwise, it uses the stencil from firstd_m
     Asserts that a_derivdir != direction of face
  */
  Real computeFirstDerivative(
                              const BaseFab<Real> & a_phic,
                              int a_derivdir,
                              int a_ivar,
                              const IntVect& a_ivin,
                              Real a_dx) const;

  ///
  /**
     compute mixed derivative (direction is unambiguous.
     x and y are the two directions tangential to face)
     at coarse point a_ivin. \\
     In two dimensions, always returns 0. \\
     Uses centered finite diff approximation if
     ivs_standard_m.contains(a_ivin ).
     It uses the stencil from mixed_sten_m.\\
     Returns 0 if SpaceDim != 3
  */
  Real computeMixedDerivative(
                              const BaseFab<Real> & a_phic,
                              int a_ivar,
                              const IntVect& a_ivin,
                              Real a_dx) const;

  ///
  QuadCFStencil();

  ///
  ~QuadCFStencil();

  ///
  QuadCFStencil(
                const Box& a_fine_domain,
                const Box& a_grid,
                const DisjointBoxLayout& a_fineBoxes,
                const DisjointBoxLayout& a_coarBoxes,
                int a_refRatio,
                int a_direction,
                Side::LoHiSide a_hiorlo);

  ///
  QuadCFStencil(
                const ProblemDomain& a_fine_domain,
                const Box& a_grid,
                const DisjointBoxLayout& a_fineBoxes,
                const DisjointBoxLayout& a_coarBoxes,
                int a_refRatio,
                int a_direction,
                Side::LoHiSide a_hiorlo);

  ///
  void   define(
                const Box& a_fine_domain,
                const Box& a_grid,
                const DisjointBoxLayout& a_fineBoxes,
                const DisjointBoxLayout& a_coarBoxes,
                int a_refRatio,
                int a_direction,
                Side::LoHiSide a_hiorlo);

  ///
  void   define(
                const ProblemDomain& a_fine_domain,
                const Box& a_grid,
                const DisjointBoxLayout& a_fineBoxes,
                const DisjointBoxLayout& a_coarBoxes,
                int a_refRatio,
                int a_direction,
                Side::LoHiSide a_hiorlo);


  ///  faster constructor that requires the use of CFStencil::buildPeriodicVector
  void   define(
                const ProblemDomain& a_fine_domain,
                const Box& a_grid,
                const Vector<Box>& a_periodicBoxes,
                const Vector<Box>& a_coarsenedPeriodicBoxes,
                const DisjointBoxLayout& a_coarBoxes,
                int a_refRatio,
                int a_direction,
                Side::LoHiSide a_hiorlo);

  /// {\bf Access functions}

  ///
  /**
     has full define function been called?  return true if so
  */
  bool isDefined() const;

  ///
  /**
     are there any interpolation points?  \\
     returns false if so.*/
  bool isEmpty() const
  {
    CH_assert(m_baseCFS.isDefined());
    return m_baseCFS.isEmpty();
  }

  ///
  /** get fine intvects which need to be interpolated  \\
      This will be empty if isEmpty() returns true*/
  const IntVectSet& getFineIVS() const
  {
    CH_assert(m_baseCFS.isDefined());
    return m_baseCFS.getFineIVS();
  }

  ///
  /** get coarse intvects that underly fiinterpivs.
      This will be empty if isEmpty() returns true*/
  const IntVectSet& getCoarIVS() const
  {
    CH_assert(m_baseCFS.isDefined());
    return m_baseCFS.getCoarIVS();
  }

  bool finePacked() const
  {
    return m_baseCFS.isPacked();
  }

  const Box& packedBox() const
  {
    return m_baseCFS.packedBox();
  }


protected:

  //has define fcn been called?
  bool m_isDefined;

  //direction of face
  int m_direction;

  ///
  /**
     set of all points which have ALL their
     derivatives taken normally.  This means
     that no stencil either leaves the domain or
     intersects the next finer level.
  */
  IntVectSet m_ivsStandard;

  ///
  /**
     set of all points which have at least
     one derivative shifted
  */
  IntVectSet m_ivsQuadd;

  /**
     Derivative stencils.  Only defined at intvects where there
     is at least one non-standard derivative to take
     (ie. at all points of ivs_Quadd)
  */

  ///
  /**
     stencil information for mixed derivatives
     undefined if SpaceDim != 3
  */
  BaseFab<DerivStencil> m_mixedSten;

  ///
  /**
     stencil information for unmixed derivatives
     (second derivatives)
  */
  BaseFab<DerivStencil> m_second[SpaceDim];

  ///
  /**
     stencil information for unmixed derivatives
     (first derivatives)
  */
  BaseFab<DerivStencil> m_firstD[SpaceDim];

  ///
  /**
     whether to drop order at any
     one-sided point
  */
  BaseFab<bool> m_dropOrd;

  ///
  CFStencil m_baseCFS;


private:
  // internally useful but not for public consumption

  //
  void  setDefaultValues();

  /*
    add fab to the stencil.  For internal use only.
  */
  void addFabToSten(const BaseFab<Real> & fabin_a,
                    DerivStencil&  sten_a);

  // there is no assignment operator for this class
  void operator= (const QuadCFStencil& cfs_in)
  {
  }

  // there is no copy constructor for this class
  QuadCFStencil(const QuadCFStencil& cfs_in)
  {
  }

  void buildStencils(const IntVectSet& ivsAllGood);

};

#include "NamespaceFooter.H"
#endif
