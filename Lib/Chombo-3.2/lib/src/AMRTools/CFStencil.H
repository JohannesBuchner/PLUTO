#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CFSTENCIL_H_
#define _CFSTENCIL_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "DisjointBoxLayout.H"
#include "Box.H"
#include "LoHiSide.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

/// Class to encapsulate coarse-fine information for AMR operations.
/**
    CFStencil is a wrapper around a bunch of things
    that allow one to do coarse-fine interpolation.
    This class should be considered internal to QuadCFInterp
    and should not be considered part of the Chombo API.
  */
class CFStencil
{

private:

protected:
  //for internal use
  void  setDefaultValues();

  //the following variables are mainly for self-identification

  // direction normal to interface
  int m_direction;

  // interface is on high or low side of grid
  Side::LoHiSide m_hiorlo;

  //
  DataIndex m_dataIndex;

  //the following variables are to be used by others

  /* fine intvects which need to be interpolated */
  IntVectSet m_fineIVS;

  /* coar intvects where slopes can be computed
      == coarsened version of fiinterp ivs*/
  IntVectSet m_coarIVS;

  bool m_isDefined;

  bool m_isPacked;

  Box  m_packedBox;

public:
  /// {\bf Access functions}

  ///
  /** returns true if CFStencil is defined */
  bool isDefined() const;

  ///
  /**
     are there any interpolation points?  \\
     returns false if so.*/
  bool isEmpty() const;

  ///
  bool isPacked() const
  {
    return  m_isPacked;
  }
  ///
  const Box& packedBox() const

  {
    return m_packedBox;
  }

  ///
  /** get fine intvects which need to be interpolated  \\
      This will be empty if isEmpty() returns true*/
  const IntVectSet& getFineIVS() const;

  ///
  /** get coarse intvects that underly fiinterpivs.
      This will be empty if isEmpty() returns true
      this is the same as getCoarInterpIVS
  */
  const IntVectSet& getCoarIVS() const ;

  ///
  CFStencil& operator= (const CFStencil& cfs_in);

  ///
  CFStencil(const CFStencil& cfs_in);

  /*dummy constructor so i can use this with
     LocalBoxLayoutData
   */
  CFStencil(const Box& a_b, int a_nComp);

  void define(const Box& a_b, int a_nComp);

  ///
  CFStencil();

  ///
  ~CFStencil();

  ///
  CFStencil(
            const Box& FineDomain,
            const Box& grid,
            const DisjointBoxLayout& Levboxes,
            const DisjointBoxLayout& LevCBoxes,
            int RefRatio,
            int Direction,
            Side::LoHiSide hiorlo);

  ///
  CFStencil(
            const ProblemDomain& FineDomain,
            const Box& grid,
            const DisjointBoxLayout& Levboxes,
            const DisjointBoxLayout& LevCBoxes,
            int RefRatio,
            int Direction,
            Side::LoHiSide hiorlo);

  ///
  void   define(
                const Box& FineDomain,
                const Box& grid,
                const DisjointBoxLayout& fine_boxes,
                const DisjointBoxLayout& coar_boxes,
                int Refratio,
                int Direction,
                Side::LoHiSide hiorlo);

  ///
  void   define(
                const ProblemDomain& FineDomain,
                const Box& grid,
                const DisjointBoxLayout& fine_boxes,
                const DisjointBoxLayout& coar_boxes,
                int Refratio,
                int Direction,
                Side::LoHiSide hiorlo);

  ///  faster define function that uses a sorted periodic vector
  /** user calls the helper function "buildPeriodicVector", then
      uses that returned Vector<Box> in the subsequent calls to this define.
  */
  void   define(
                const ProblemDomain& a_fineDomain,
                const Box& a_grid,
                const Vector<Box>& a_periodicVector,
                int a_refRatio,
                int a_direction,
                Side::LoHiSide a_hiorlo);
  ///
  void clear();


  static void buildPeriodicVector(Vector<Box>& a_periodicVector,
                                  const ProblemDomain& a_fineDomain,
                                  const DisjointBoxLayout& a_fineBoxes);

};
#include "NamespaceFooter.H"
#endif
