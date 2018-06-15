#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NOREFINEMENTIMPLEM_H_
#define _NOREFINEMENTIMPLEM_H_

#include "NamespaceHeader.H"

template <int dim> NoRefinement<dim>::NoRefinement()
{
}

template <int dim> NoRefinement<dim>::~NoRefinement()
{
}

template <int dim> bool NoRefinement<dim>::doRefine(IndexTM<int,dim>          & a_refineInDir,
                                                    const CutCellMoments<dim> & a_ccm,
                                                    const int                 & a_numberOfRefinements)
{
  for (int idir = 0; idir < dim; idir++)
  {
    a_refineInDir[idir] = 0;
  }

  return false;
}

#include "NamespaceFooter.H"

#endif
