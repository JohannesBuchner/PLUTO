#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _FIXEDREFINEMENTIMPLEM_H_
#define _FIXEDREFINEMENTIMPLEM_H_

#include "MayDay.H"

#include "NamespaceHeader.H"

template <int dim> FixedRefinement<dim>::FixedRefinement()
{
}

template <int dim> FixedRefinement<dim>::FixedRefinement(const int & a_maxNumberOfRefinements)
{
  setMaxNumberOfRefinements(a_maxNumberOfRefinements);
}

template <int dim> FixedRefinement<dim>::~FixedRefinement()
{
}

template <int dim> bool FixedRefinement<dim>::doRefine(IndexTM<int,dim>          & a_refineInDir,
                                                       const CutCellMoments<dim> & a_ccm,
                                                       const int                 & a_numberOfRefinements)
{
  bool retval = false;

  if (a_numberOfRefinements < m_maxNumberOfRefinements)
  {
    retval = true;

    for (int idir = 0; idir < dim; idir++)
    {
      a_refineInDir[idir] = 1;
    }
  }

  return retval;
}

template <int dim> void FixedRefinement<dim>::setMaxNumberOfRefinements(const int & a_maxNumberOfRefinements)
{
  if (a_maxNumberOfRefinements < 0)
  {
    MayDay::Abort("FixedRefinement<dim>::setNumberOfRefinements - maxNumberOfRefinements must be >= 0");
  }

  m_maxNumberOfRefinements = a_maxNumberOfRefinements;
}

template <int dim> int FixedRefinement<dim>::getMaxNumberOfRefinements()
{
  return m_maxNumberOfRefinements;
}

#include "NamespaceFooter.H"

#endif
