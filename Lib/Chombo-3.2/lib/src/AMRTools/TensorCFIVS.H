#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TENSORCFIVS_H_
#define _TENSORCFIVS_H_

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
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

///Class to get fine intvectset at CF interface
/**
Class to get fine intvectset at cf interface.
This class should be considered internal to AMRSolver and should
not be considered part of the Chombo API.
*/
class TensorCFIVS
{

private:

protected:
  //for internal use
  void  setDefaultValues();

  //the following variables are mainly for self-identification

  /* fine intvects which need to be interpolated */
  IntVectSet fiinterp_ivs;

  bool isdefined;

public:

  ///
  bool isDefined() const;

  ///
  /** Get fine intvects which need to be interpolated.
      This will be empty if isEmpty() returns true*/
  const IntVectSet& getFineIVS() const ;

  ///
  TensorCFIVS();

  ///
  ~TensorCFIVS();

  ///
  TensorCFIVS(
        const Box& Domain,
        const Box& box_in,
        const DisjointBoxLayout& Levboxes,
        int Direction,
        Side::LoHiSide hiorlo);


  ///
  TensorCFIVS(
        const ProblemDomain& Domain,
        const Box& box_in,
        const DisjointBoxLayout& Levboxes,
        int Direction,
        Side::LoHiSide hiorlo);

  ///
  void define(
              const Box& Domain,
              const Box& box_in,
              const DisjointBoxLayout& fine_boxes,
              int Direction,
              Side::LoHiSide hiorlo);

  ///
  void define(
              const ProblemDomain& Domain,
              const Box& box_in,
              const DisjointBoxLayout& fine_boxes,
              int Direction,
              Side::LoHiSide hiorlo);

private:
  void operator= (const TensorCFIVS& levcfs_in)
  {
  }

  TensorCFIVS(const TensorCFIVS& levcfs_in)
  {
  }
};
#include "NamespaceFooter.H"
#endif
