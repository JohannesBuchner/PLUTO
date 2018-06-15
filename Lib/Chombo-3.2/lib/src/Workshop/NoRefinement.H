#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NOREFINEMENT_H_
#define _NOREFINEMENT_H_

#include "RefinementCriterion.H"

#include "NamespaceHeader.H"

///
/**
   This refinement criterion always false.  Used to stop subdivision.
 */
template <int dim> class NoRefinement: public RefinementCriterion<dim>
{
public:
  /// Null constructor
  NoRefinement();

  /// Destructor
  virtual ~NoRefinement();

  /// Never refine - return false and "a_refineInDir" zero in all components.
  /**
     Never refine - return false and "a_refineInDir" zero in all components.
   */
  virtual bool doRefine(IndexTM<int,dim>          & a_refineInDir,
                        const CutCellMoments<dim> & a_ccm,
                        const int                 & a_numberOfRefinements);

protected:
};

#include "NamespaceFooter.H"

#include "NoRefinementImplem.H"

#endif
