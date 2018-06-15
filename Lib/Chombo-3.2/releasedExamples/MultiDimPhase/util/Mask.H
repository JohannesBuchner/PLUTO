#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MASK_H_
#define _MASK_H_

#include "BaseFab.H"
#include "BoxLayout.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

// class to determine coarse-fine validity info
class Mask
{

public:
  enum cellType
  {
    maskCovered  = -2,
    maskPhysical = -1,
    maskCopy     =  0,
    maskCoarse   =  1
  };

  void buildMask(BaseFab<int>& a_mask,
                 const ProblemDomain& a_dProblem,
                 const BoxLayout& a_grids,
                 const BoxLayout* a_fineGridsPtr=NULL,
                 int a_nRefFine = -1);

  void buildMasks(LevelData<BaseFab <int> >& a_masks,
                  const ProblemDomain& a_dProblem,
                  const BoxLayout& a_grids,
                  const BoxLayout* a_fineGridsPtr=NULL,
                  int a_nRefFine = -1);

  // deprecated interface
  void buildMask(BaseFab<int>& a_mask, const Box& a_dProblem,
                 const BoxLayout& a_grids,
                 const BoxLayout* a_fineGridsPtr=NULL,
                 int a_nRefFine = -1);

  // deprecated interface
  void buildMasks(LevelData<BaseFab <int> >& a_masks,
                  const Box& a_dProblem,
                  const BoxLayout& a_grids,
                  const BoxLayout* a_fineGridsPtr=NULL,
                  int a_nRefFine = -1);

};



#include "NamespaceFooter.H"

#endif


