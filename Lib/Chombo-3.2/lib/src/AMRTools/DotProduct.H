#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DOTPRODUCT_H_
#define _DOTPRODUCT_H_

#include "BoxLayout.H"
#include "BoxLayoutData.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

// this is the all-component dot-product
extern Real DotProduct(const BoxLayoutData<FArrayBox>& a_dataOne,
                       const BoxLayoutData<FArrayBox>& a_dataTwo,
                       const BoxLayout&                a_dblIn);

extern Real DotProduct(const BoxLayoutData<FArrayBox>& a_dataOne,
                       const BoxLayoutData<FArrayBox>& a_dataTwo,
                       const BoxLayout& a_dblIn,
                       const Interval& a_comps);

#include "NamespaceFooter.H"
#endif
