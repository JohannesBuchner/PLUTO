#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FIXEDREFINEMENT_H_
#define _FIXEDREFINEMENT_H_

#include "RefinementCriterion.H"

#include "NamespaceHeader.H"

///
/**
   This refinement criterion refines a fixed number of levels.
 */
template <int dim> class FixedRefinement: public RefinementCriterion<dim>
{
public:
  /// Null constructor
  FixedRefinement();

  /// Constructor which sets number of refinements
  FixedRefinement(const int & a_maxNumberOfRefinements);

  /// Destructor
  virtual ~FixedRefinement();

  /// Refine a fixed number of times
  /**
     Refine if the number of refinements already finished is smaller than the  max number of refinements.
   */
  virtual bool doRefine(IndexTM<int,dim>          & a_refineInDir,
                        const CutCellMoments<dim> & a_ccm,
                        const int                 & a_numberOfRefinements);

  /// Set the number of refinements
  virtual void setMaxNumberOfRefinements(const int & a_maxNumberOfRefinements);

  /// Get the number of refinements
  virtual int getMaxNumberOfRefinements();

protected:
  int m_maxNumberOfRefinements;
};

#include "NamespaceFooter.H"

#include "FixedRefinementImplem.H"

#endif
