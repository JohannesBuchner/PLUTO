#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "MayDay.H"

#include "NoRefinement.H"
#include "LSProblem.H"
#include "ComputeCutCellMoments.H"

#include "NamespaceHeader.H"

// Null constructor
ComputeCutCellMoments<1>::ComputeCutCellMoments()
{
}

// Copy constructor
ComputeCutCellMoments<1>::ComputeCutCellMoments(const ComputeCutCellMoments<1>& a_thisComputeCutCellMoments)

  :m_cutCellMoments(a_thisComputeCutCellMoments.m_cutCellMoments)
{
}

// This constructor is used in the recursion
ComputeCutCellMoments<1>::ComputeCutCellMoments(const IFData<1> & a_info)
  :m_cutCellMoments(a_info)
{
}

// Destructor
ComputeCutCellMoments<1>::~ComputeCutCellMoments()
{
}

// Integrate along line segments aligned in a coordinate direction
#if RECURSIVE_GEOMETRY_GENERATION == 0
void ComputeCutCellMoments<1>::computeMoments(const int              & a_order,
                                              const int              & a_degree,
#else
void ComputeCutCellMoments<1>::computeMoments(const int              & a_orderPmax,
                                              const int              & a_degreePmax,
#endif
                                              const bool             & a_useConstraints,
                                              RefinementCriterion<1> & a_refinementCriterion)
{
  int lo = 0;
  int hi = 1;
  int loSign = m_cutCellMoments.m_IFData.m_cornerSigns[lo];
  int hiSign = m_cutCellMoments.m_IFData.m_cornerSigns[hi];

  // If entire edge out of the fluid, then moments = 0.0
  if (loSign <= ON && hiSign<= ON)
  {
#if RECURSIVE_GEOMETRY_GENERATION == 0
    for (int iDegree = 0; iDegree <= a_degree; ++iDegree)
#else
    for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
#endif
    {
      // Definition of m_cutCellMoments.m_moments typedef requires that we define degree as
      // a oneTuple:
      IndexTM<int,1> degree;
      degree[0] = iDegree;
      m_cutCellMoments.m_moments[degree] = 0.0;
    }
  }
  else
  {
    // Assign loPt and hiPt in m_cutCellMoments.m_IFData.m_parentCoord system

    Real loPt = LARGEREALVAL;
    Real hiPt = LARGEREALVAL;

    // m_origin is an IndexTM<Real,1>, which implies we need to use [0] everywhere
    // m_intersection is undefined if hiSign >= ON  && loSign >= ON
    if (loSign >= ON)
    {
      loPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
                 -0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0],
                 m_cutCellMoments.m_IFData.m_cellCenterCoord,
                 0);
    }
    else
    {
      CH_assert(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]
                <= m_cutCellMoments.m_IFData.m_intersection
                &&
                m_cutCellMoments.m_IFData.m_intersection
                <= 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]);

      loPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
                 m_cutCellMoments.m_IFData.m_intersection,
                 m_cutCellMoments.m_IFData.m_cellCenterCoord,
                 0);
    }

    if (hiSign >= ON)
    {
      hiPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
                 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0],
                 m_cutCellMoments.m_IFData.m_cellCenterCoord,
                 0);
    }
    else
    {
      CH_assert(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]
                <= m_cutCellMoments.m_IFData.m_intersection
                &&
                m_cutCellMoments.m_IFData.m_intersection
                <= 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]);

      hiPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
                 m_cutCellMoments.m_IFData.m_intersection,
                 m_cutCellMoments.m_IFData.m_cellCenterCoord,
                 0);
    }

    // Integrate x^degree over the line segment[loPt,hiPt]
#if RECURSIVE_GEOMETRY_GENERATION == 0
    computeMomentsUsingBinomial(loPt,hiPt,loSign,hiSign,a_degree);
#else
    computeMomentsUsingBinomial(loPt,hiPt,loSign,hiSign,a_degreePmax);
#endif
  }
}

void ComputeCutCellMoments<1>::simpleComputeMoments(const Real & a_loPt,
                                                    const Real & a_hiPt,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                                                    const int  & a_degree)
#else
                                                    const int  & a_degreePmax)
#endif
{
#if RECURSIVE_GEOMETRY_GENERATION == 0
  for (int iDegree = 0; iDegree <= a_degree; ++iDegree)
#else
  for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
#endif
  {
    //definition of m_cutCellMoments.m_moments typedef requires that we define degree thus:
    IndexTM<int,1>degree;
    degree[0] = iDegree;
    m_cutCellMoments.m_moments[degree] = pow(a_hiPt,iDegree + 1) - pow(a_loPt,iDegree +1);
    //    Real dxFactor = pow(m_cutCellMoments.m_IFData.m_globalCoord.m_dx[0],iDegree + 1);
    //m_cutCellMoments.m_moments[degree] *= dxFactor;
    m_cutCellMoments.m_moments[degree] /= (iDegree + 1);
  }
}

void ComputeCutCellMoments<1>::computeMomentsUsingBinomial(const Real & a_loPt,
                                                           const Real & a_hiPt,
                                                           const int  & a_loSign,
                                                           const int  & a_hiSign,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                                                           const int  & a_degree)
#else
                                                           const int  & a_degreePmax)
#endif
{
#if RECURSIVE_GEOMETRY_GENERATION == 0
  for (int iDegree = 0; iDegree <= a_degree; ++iDegree)
#else
  for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
#endif
  {
    Real epsilon = a_hiPt - a_loPt;

    // Definition of m_cutCellMoments.m_moments typedef requires that we define degree thus:
    IndexTM<int,1> degree;
    degree[0] = iDegree;
    m_cutCellMoments.m_moments[degree] = 0.0;

    // New method without substracting higher order terms
    for (int j = 1; j <= iDegree + 1; j++)
    {
      int bigger = j;
      int smaller = j;
      if (iDegree + 1 - j > j)
      {
        bigger = iDegree + 1 - j;
      }
      else
      {
        smaller = iDegree + 1 - j;
      }

      int numerator = 1;
      for (int i = bigger + 1; i <= iDegree + 1; ++i)
      {
        numerator *= i;
      }

      int denominator = 1;
      for (int i = 1; i <= smaller; ++i)
      {
        denominator *= i;
      }

      Real factor = numerator / denominator;
      if (a_loSign >= ON)
      {
        m_cutCellMoments.m_moments[degree] += factor * pow(a_loPt,iDegree + 1 - j) * pow(epsilon,j);
      }
      else if (a_hiSign >= ON)
      {
        m_cutCellMoments.m_moments[degree] -= factor * pow(a_hiPt,iDegree + 1 - j) * pow(epsilon,j) * pow(-1.0,j);
      }
    }

    //Real dxFactor = pow(m_cutCellMoments.m_IFData.m_globalCoord.m_dx[0],iDegree + 1);
    //m_cutCellMoments.m_moments[degree] *= dxFactor;
    m_cutCellMoments.m_moments[degree] /= iDegree + 1;
  }
}

void ComputeCutCellMoments<1>::print(ostream & a_out) const
{
  m_cutCellMoments.print(a_out);
}

void ComputeCutCellMoments<1>::operator=(const ComputeCutCellMoments<1> & a_computeCutCellMoments)
{
  // Only copy of the two objects are distinct
  if (this != &a_computeCutCellMoments)
  {
    m_cutCellMoments = a_computeCutCellMoments.m_cutCellMoments;
  }
}

#include "NamespaceFooter.H"
