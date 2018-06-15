#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPUTECUTCELLMOMENTSIMPLEM_H_
#define _COMPUTECUTCELLMOMENTSIMPLEM_H_

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

#include "NamespaceHeader.H"

// Null constructor
template <int dim> ComputeCutCellMoments<dim>::ComputeCutCellMoments()
{
}

// Copy constructor
template <int dim> ComputeCutCellMoments<dim>::ComputeCutCellMoments(const ComputeCutCellMoments<dim>& a_computeCutCellMoments)
  :m_cutCellMoments(a_computeCutCellMoments.m_cutCellMoments)
{
}

// Use this for initializing
template <int dim> ComputeCutCellMoments<dim>::ComputeCutCellMoments(const IFData<dim>& a_info)
  :m_cutCellMoments(a_info)
{
  m_cutCellMoments.m_bdCCOn = false;

  for (int hilo = 0; hilo < 2; ++hilo)
  {
    for (int idir = 0; idir < dim; ++idir)
    {
      // Identifier for which boundary cutCellMoment
      Iv2 bdId;
      bdId[BDID_DIR] = idir;
      bdId[BDID_HILO] = hilo;

#if RECURSIVE_GEOMETRY_GENERATION == 0
      IFData<dim-1> reducedInfo(a_info,idir,hilo);
#else
      IFData<dim-1> reducedInfo(a_info,a_info.m_maxOrder+1,idir,hilo);
#endif
      CutCellMoments<dim-1>bdCutCellMoments(reducedInfo);

      m_cutCellMoments.m_bdCutCellMoments[bdId] = bdCutCellMoments;


      // Notice whether at least one lower dimensional cutCell is on the
      // interface
      if (reducedInfo.m_allVerticesOn)
      {
        m_cutCellMoments.m_bdCCOn = true;
      }
    }
  }
}

// Destructor
template <int dim> ComputeCutCellMoments<dim>::~ComputeCutCellMoments()
{
}

#if RECURSIVE_GEOMETRY_GENERATION == 0
template <int dim> void ComputeCutCellMoments<dim>::computeMoments(const int                & a_order,
                                                                   const int                & a_degreeP,
                                                                   const bool               & a_useConstraints,
                                                                   RefinementCriterion<dim> & a_refinementCriterion,
                                                                   const int                & a_numberOfRefinements)
{
  CH_assert(m_cutCellMoments.m_IFData.m_maxOrder >= a_order);

  Vector<Real> RNorm(3);
  for (int i = 0; i < 3; i++)
  {
    RNorm[i] = LARGEREALVAL;
  }

  for (int i = 0; i < a_degreeP + 1; i++)
  {
    m_cutCellMoments.m_residual.push_back(RNorm);
  }

  if (m_cutCellMoments.m_IFData.m_allVerticesOut)
  {
    for (int iDegree = a_degreeP; iDegree >= 0; --iDegree)
    {
      LSProblem<dim> lsp(iDegree,a_useConstraints);

      // Fill moments
      const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
      for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
      {
        m_cutCellMoments.m_EBmoments[it->first] = 0.0;
      }

      const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
      for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin(); it != monDegreePLess1.end(); ++it)
      {
        m_cutCellMoments.m_moments[it->first] = 0.0;
      }
    }
  }
  else if (m_cutCellMoments.m_IFData.m_allVerticesIn && !m_cutCellMoments.m_bdCCOn)
  {
    for (int iDegree = a_degreeP; iDegree >= 0; --iDegree)
    {
      LSProblem<dim> lsp(iDegree,a_useConstraints);

      // Fill moments of degree P and P-1
      const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
      for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
      {
        m_cutCellMoments.m_EBmoments[it->first] = 0.0;
      }

      const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
      for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin(); it != monDegreePLess1.end(); ++it)
      {
        m_cutCellMoments.m_moments[it->first] = m_cutCellMoments.fullCellQuadrature(it->first,m_cutCellMoments.m_IFData.m_parentCoord);
      }
    }
  }
  else
  {
    for (typename BdCutCellMoments::iterator it = m_cutCellMoments.m_bdCutCellMoments.begin(); it != m_cutCellMoments.m_bdCutCellMoments.end(); ++it)
    {
      NoRefinement<dim-1> noRefinement;

      ComputeCutCellMoments<dim-1> subProblem(it->second.m_IFData);
      subProblem.computeMoments(a_order,a_degreeP+1,a_useConstraints,noRefinement);
      it->second = subProblem.m_cutCellMoments;
      //pout()<< "Subproblem constraints: " << it->second.m_numActiveBounds << endl;
    }

    for (int iDegree = a_degreeP; iDegree >= 0; --iDegree)
    {
      // Make a LSProb
      IvDim zeroDerivative = IndexTM<int,dim>::Zero;
      LSProblem<dim> lsp(a_order,iDegree,a_useConstraints,m_cutCellMoments.m_IFData.m_normalDerivatives[zeroDerivative]);

      Vector<Real> rhs = computeRhs(lsp,a_order);

      if (a_useConstraints)
      {
        lsp.computeBounds(m_cutCellMoments.m_IFData.m_globalCoord.m_dx,m_cutCellMoments);
      }


      // Solve the problem and return residual
      int lsCode = lsp.invertNormalEq(rhs,m_cutCellMoments.m_residual[iDegree]);
      if (lsCode != 0)
        {
          pout() << "Geometry generation least squares problem failed with residual: "
                 << m_cutCellMoments.m_residual[iDegree]<< endl;
          lsp.print(pout());

          pout () << endl << "Problem occurred generating geometry for these cut cell moments:\n"
                  << m_cutCellMoments << endl;

          //MayDay::Error("Geometry generation error.");
          bool status = false;
          a_refinementCriterion.setConstrantSuccessStatus(status);
        }

      // Record number of constraint violations
      m_cutCellMoments.m_numActiveBounds += lsp.numActiveBounds();

      //if (m_cutCellMoments.m_numActiveBounds > 0)
      //  {
      //    pout() <<" Greater than zero. Dim= " << dim << " num="<<m_cutCellMoments.m_numActiveBounds << endl;
      //      }

      // Fill moments
      const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
      for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin();
           it != monDegreePLess1.end(); ++it)
      {
        m_cutCellMoments.m_moments[it->first] = lsp.getUnknown(it->second + lsp.getNumberDegP());
      }

      const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
      for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
      {
        m_cutCellMoments.m_EBmoments[it->first] = lsp.getUnknown(it->second);
      }
    }

    // Move moments from local coord to parent coord
    IndexTM<Real,dim> delta = m_cutCellMoments.m_IFData.m_parentCoord.m_origin;
    delta                  -= m_cutCellMoments.m_IFData.m_localCoord .m_origin;

    PthMoment copyMoments = m_cutCellMoments.m_moments;
    for (typename PthMoment::const_iterator it = copyMoments.begin();
         it != copyMoments.end(); ++it)
      {
        IvDim mono = it->first;
        m_cutCellMoments.m_moments[mono] = m_cutCellMoments.changeMomentCoordinates(copyMoments, mono, delta);
      }

    PthMoment copyEBMoments = m_cutCellMoments.m_EBmoments;
    for (typename PthMoment::const_iterator it = copyEBMoments.begin(); it != copyEBMoments.end(); ++it)
      {
        IvDim mono = it->first;
        m_cutCellMoments.m_EBmoments[mono] = m_cutCellMoments.changeMomentCoordinates(copyEBMoments, mono, delta);
      }

    // Move bdCutCell moments from parent coord to cell center coord. From here on bdCutCell moments will
    // not be used in any least squares problem

    for (typename BdCutCellMoments::iterator it = m_cutCellMoments.m_bdCutCellMoments.begin(); it != m_cutCellMoments.m_bdCutCellMoments.end(); ++it)
    {
      it->second.changeMomentCoordinatesToCellCenter();
    }
  }

  IndexTM<int,dim> refineInDir;

  if (!m_cutCellMoments.isCovered() && !m_cutCellMoments.isRegular() && a_refinementCriterion.baseDoRefine(refineInDir,m_cutCellMoments,a_numberOfRefinements))
  {
    // Refine creates a vector of 2^dim cutcellmoments corresponding to
    // the refined cells inside the current cell
    Vector<CutCellMoments<dim> > refinedCCMoms = refine(a_order,a_degreeP,a_useConstraints,a_refinementCriterion,refineInDir,a_numberOfRefinements);

    // addMomentMaps add up the refined moments in the current
    // CutCellMoments maps
    addMomentMaps(refinedCCMoms,a_degreeP,a_useConstraints);

    // computeResiduals(a_order,a_degreeP,a_useConstraints);
    computeResiduals(refinedCCMoms,a_degreeP);


    // Move the origin of moments from cell center to parent center.
    IndexTM<Real,dim> delta = m_cutCellMoments.m_IFData.m_parentCoord    .m_origin;
    delta                  -= m_cutCellMoments.m_IFData.m_cellCenterCoord.m_origin;

    PthMoment copyMoments = m_cutCellMoments.m_moments;
    for (typename PthMoment::const_iterator it = copyMoments.begin();
         it != copyMoments.end(); ++it)
      {
        IvDim mono = it->first;
        m_cutCellMoments.m_moments[mono] = m_cutCellMoments.changeMomentCoordinates(copyMoments, mono, delta);
      }

    PthMoment copyEBMoments = m_cutCellMoments.m_EBmoments;
    for (typename PthMoment::const_iterator it = copyEBMoments.begin(); it != copyEBMoments.end(); ++it)
      {
        IvDim mono = it->first;
        m_cutCellMoments.m_EBmoments[mono] = m_cutCellMoments.changeMomentCoordinates(copyEBMoments, mono, delta);
      }
  }
}
#else
// Solve
template <int dim> void ComputeCutCellMoments<dim>::computeMoments(const int                & a_orderPmax,
                                                                   const int                & a_degreePmax,
                                                                   const bool               & a_useConstraints,
                                                                   RefinementCriterion<dim> & a_refinementCriterion,
                                                                   const int                & a_numberOfRefinements)
{
  int zerothOrderOfAccuracy = 0;
  int highestDegree         = a_orderPmax + a_degreePmax;

  Vector<Real> RNorm(3);
  for (int i = 0; i < 3; i++)
    {
      RNorm[i] = LARGEREALVAL;
    }

  for (int i = 0; i <= highestDegree; i++)
    {
      m_cutCellMoments.m_residual.push_back(RNorm);
    }

  m_boundaryMomentsComputed = false;

  computeMomentsRecursively(zerothOrderOfAccuracy,
                            highestDegree,
                            a_useConstraints,
                            a_refinementCriterion,
                            a_numberOfRefinements);

  IndexTM<int,dim> refineInDir;

  if (!m_cutCellMoments.isCovered() && !m_cutCellMoments.isRegular() && a_refinementCriterion.doRefine(refineInDir,m_cutCellMoments,a_numberOfRefinements))
    {
      // Refine creates a vector of 2^dim cutcellmoments corresponding to
      // the refined cells inside the current cell
      Vector<CutCellMoments<dim> > refinedCCMoms = refine(a_orderPmax,a_degreePmax,a_useConstraints,a_refinementCriterion,refineInDir,a_numberOfRefinements);

      // addMomentMaps add up the refined moments in the current
      // CutCellMoments maps
      addMomentMaps(refinedCCMoms,a_degreePmax,a_useConstraints);

      // computeResiduals(a_orderPmax,a_degreePmax,a_useConstraints);
      computeResiduals(refinedCCMoms,a_degreePmax);

      // Move the origin of moments from cell center to parent center.
      IndexTM<Real,dim> delta = m_cutCellMoments.m_IFData.m_parentCoord.m_origin;
      delta                  -= m_cutCellMoments.m_IFData.m_cellCenterCoord .m_origin;

      PthMoment copyMoments = m_cutCellMoments.m_moments;
      for (typename PthMoment::const_iterator it = copyMoments.begin();
           it != copyMoments.end(); ++it)
        {
          IvDim mono = it->first;
          m_cutCellMoments.m_moments[mono] = m_cutCellMoments.changeMomentCoordinates(copyMoments, mono, delta);
        }

      PthMoment copyEBMoments = m_cutCellMoments.m_EBmoments;
      for (typename PthMoment::const_iterator it = copyEBMoments.begin(); it != copyEBMoments.end(); ++it)
        {
          IvDim mono = it->first;
          m_cutCellMoments.m_EBmoments[mono] = m_cutCellMoments.changeMomentCoordinates(copyEBMoments, mono, delta);
        }

      // Move bdCutCell moments from parent coord to cell center coord. From here on bdCutCell moments will
      // not be used in any least squares problem

      for (typename BdCutCellMoments::iterator it = m_cutCellMoments.m_bdCutCellMoments.begin(); it != m_cutCellMoments.m_bdCutCellMoments.end(); ++it)
        {
          it->second.changeMomentCoordinatesToCellCenter();
        }
    }

  // Move moments from local coord to parent coord
  IndexTM<Real,dim> delta = m_cutCellMoments.m_IFData.m_parentCoord.m_origin;
  delta                  -= m_cutCellMoments.m_IFData.m_localCoord .m_origin;

  PthMoment copyMoments = m_cutCellMoments.m_moments;
  for (typename PthMoment::const_iterator it = copyMoments.begin();
       it != copyMoments.end(); ++it)
    {
      IvDim mono = it->first;
      m_cutCellMoments.m_moments[mono] = m_cutCellMoments.changeMomentCoordinates(copyMoments, mono, delta);
    }

  PthMoment copyEBMoments = m_cutCellMoments.m_EBmoments;
  for (typename PthMoment::const_iterator it = copyEBMoments.begin(); it != copyEBMoments.end(); ++it)
    {
      IvDim mono = it->first;
      m_cutCellMoments.m_EBmoments[mono] = m_cutCellMoments.changeMomentCoordinates(copyEBMoments, mono, delta);
    }
}

template <int dim> void ComputeCutCellMoments<dim>::computeMomentsRecursively(const int                & a_orderPmax,
                                                                              const int                & a_degreePmax,
                                                                              const bool               & a_useConstraints,
                                                                              RefinementCriterion<dim> & a_refinementCriterion,
                                                                              const int                & a_numberOfRefinements)
{
  CH_assert(m_cutCellMoments.m_IFData.m_maxOrder >= a_orderPmax);

  if (m_cutCellMoments.m_IFData.m_allVerticesOut)
    {
      for (int iOrder = 0; iOrder <= a_orderPmax; ++iOrder)
        {
          LSProblem<dim> lsp(iOrder + a_degreePmax, a_useConstraints);

          // Fill moments
          const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
          for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
            {
              m_cutCellMoments.m_EBmoments[it->first] = 0.0;
            }

          const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
          for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin(); it != monDegreePLess1.end(); ++it)
            {
              m_cutCellMoments.m_moments[it->first] = 0.0;
            }
        }
    }
  else if (m_cutCellMoments.m_IFData.m_allVerticesIn && !m_cutCellMoments.m_bdCCOn)
    {
      for (int iOrder = 0; iOrder <= a_orderPmax; ++iOrder)
        {
          LSProblem<dim> lsp(iOrder + a_degreePmax, a_useConstraints);

          // Fill moments of degree P and P-1
          const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
          for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
            {
              m_cutCellMoments.m_EBmoments[it->first] = 0.0;
            }

          const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
          for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin(); it != monDegreePLess1.end(); ++it)
            {
              m_cutCellMoments.m_moments[it->first] = m_cutCellMoments.fullCellQuadrature(it->first,m_cutCellMoments.m_IFData.m_parentCoord);
            }
        }
    }
  else
    {
      // Only compute the boundary moments if they haven't already been
      // computed (earlier in the recursion)
      if (!m_boundaryMomentsComputed)
      {
        for (typename BdCutCellMoments::iterator it = m_cutCellMoments.m_bdCutCellMoments.begin(); it != m_cutCellMoments.m_bdCutCellMoments.end(); ++it)
          {
            NoRefinement<dim-1> noRefinement;

            ComputeCutCellMoments<dim-1> subProblem(it->second.m_IFData);
            subProblem.computeMoments(a_orderPmax,a_degreePmax+1,a_useConstraints,noRefinement);
            it->second = subProblem.m_cutCellMoments;
          }

        m_boundaryMomentsComputed = true;
      }

      // Make a LSProb
      IvDim zeroDerivative = IndexTM<int,dim>::Zero;
      LSProblem<dim> lsp(a_orderPmax,a_degreePmax,a_useConstraints,m_cutCellMoments.m_IFData.m_normalDerivatives[zeroDerivative]);

      Vector<Real> rhs = computeRhs(lsp,a_orderPmax);

      if (a_useConstraints)
        {
          lsp.computeBounds(m_cutCellMoments.m_IFData.m_globalCoord.m_dx,m_cutCellMoments);
        }

      // Solve the problem and return residual
      int lsCode=lsp.invertNormalEq(rhs,m_cutCellMoments.m_residual[a_degreePmax]);
      if (lsCode != 0)
        {
          pout() << "Geometry generation least squares problem failed with residual:"
                 << m_cutCellMoments.m_residual[a_degreePmax]<< endl;
          lsp.print(pout());
          pout () << "Problem occurred generating geometry for these cut cell moments: " << m_cutCellMoments << endl;
          MayDay::Error("Geometry generation error.[2]");
        }


      // Fill moments
      const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
      for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin();
           it != monDegreePLess1.end(); ++it)
        {
          m_cutCellMoments.m_moments[it->first] = lsp.getUnknown(it->second + lsp.getNumberDegP());
        }

      const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
      for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
        {
          m_cutCellMoments.m_EBmoments[it->first] = lsp.getUnknown(it->second);
        }

    }

  if (a_degreePmax > 0)
    {
      computeMomentsRecursively(a_orderPmax + 1,
                                a_degreePmax - 1,
                                a_useConstraints,
                                a_refinementCriterion,
                                a_numberOfRefinements);

    }
}
#endif

template <int dim> Vector<Real> ComputeCutCellMoments<dim>::computeRhs(LSProblem<dim> & a_lsp,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                                                                       const int      & a_order)
#else
                                                                       const int      & a_orderPmax)
#endif
{
  // Resize rhs
  int numEq = dim*a_lsp.getNumberDegP();
  Vector<Real> rhs(numEq);

  // For each moment iterate thru bd CutCellMoments incrementing same comp
  // of rhs
  const LocPthMoment& locMap = a_lsp.getLocMonomialMapDegreeP();
  for (typename LocPthMoment::const_iterator it = locMap.begin(); it != locMap.end(); ++it)
  {
    int jth = it->first;
    IvDim mono = it->second;

    int hiSide = 1;
    int loSide = 0;
    Iv2 bdId;

    for (int idir = 0; idir < dim; ++idir)
    {
      // Which lower dimensional monomial corresponds (mono,j)
      IndexTM<int,dim-1> mono1Less;
      for (int jdir = 0; jdir < dim; ++jdir)
      {
        if (jdir < idir)
        {
          mono1Less[jdir] = mono[jdir];
        }
        else if (jdir > idir)
        {
          mono1Less[jdir-1] = mono[jdir];
        }
      }

      bdId[0] = idir;
      bdId[1] = hiSide;

      Real hiMom = m_cutCellMoments.m_bdCutCellMoments[bdId].m_moments[mono1Less];

      bdId[1] = loSide;

      Real loMom = m_cutCellMoments.m_bdCutCellMoments[bdId].m_moments[mono1Less];
      int exponent = it->second[idir];

      Real loSideValue;
      Real hiSideValue;

      loSideValue = m_cutCellMoments.m_IFData.m_localCoord.convertDir(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[idir],
                                                     m_cutCellMoments.m_IFData.m_cellCenterCoord,
                                                     idir);

      hiSideValue = m_cutCellMoments.m_IFData.m_localCoord.convertDir( 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[idir],
                                                     m_cutCellMoments.m_IFData.m_cellCenterCoord,
                                                     idir);

      Real loFactor = pow(loSideValue,exponent);
      Real hiFactor = pow(hiSideValue,exponent);

      rhs[(dim*jth) + idir] = hiMom*hiFactor - loMom*loFactor;

      // Add the Taylor series terms
#if RECURSIVE_GEOMETRY_GENERATION == 0
      for (int order = 1; order <= a_order; order++)
#else
      for (int order = 1; order <= a_orderPmax; order++)
#endif
      {
        Vector<IvDim> taylorMonomials;

        generateMultiIndices(taylorMonomials,order);

        for (int i = 0; i < taylorMonomials.size(); i++)
        {
          const IvDim & taylorMonomial = taylorMonomials[i];

          IvDim totalMonomial = mono + taylorMonomial;

          if (m_cutCellMoments.m_EBmoments.find(totalMonomial) !=
              m_cutCellMoments.m_EBmoments.end())
          {
            Real normalDerivative = m_cutCellMoments.m_IFData.m_normalDerivatives[taylorMonomial][idir];
            Real fact = factorial(taylorMonomial);

            Real moment = m_cutCellMoments.m_EBmoments[totalMonomial];

            rhs[(dim*jth) + idir] += normalDerivative * moment / fact;
          }
#if RECURSIVE_GEOMETRY_GENERATION != 0
          else
          {
            MayDay::Error("Unable to find needed monomial for Taylor series");
          }
#endif
        }
      }
    }
  }

  return rhs;
}

#if RECURSIVE_GEOMETRY_GENERATION == 0
template <int dim> Vector< CutCellMoments<dim> > ComputeCutCellMoments<dim>::refine(const int                & a_order,
                                                                                    const int                & a_degreeP,
#else
template <int dim> Vector< CutCellMoments<dim> > ComputeCutCellMoments<dim>::refine(const int                & a_orderPmax,
                                                                                    const int                & a_degreePmax,
#endif
                                                                                    const bool               & a_useConstraints,
                                                                                    RefinementCriterion<dim> & a_refinementCriterion,
                                                                                    const IndexTM<int,dim>   & a_refineInDir,
                                                                                    const int                & a_numberOfRefinements)
{
  Vector<CutCellMoments<dim> > refinedCutCellMoments;

  // This function refines the current cell by a factor of two in all
  // directions where a_refineInDir is non-zero
  IndexTM<Real,dim> refinedDx = m_cutCellMoments.m_IFData.m_globalCoord.m_dx;

  // Current cell iterator and maximum iterator
  IndexTM<int,dim> curIter;
  IndexTM<int,dim> maxIter;

  for (int idir = 0; idir < dim; idir++)
  {
    if (a_refineInDir[idir] != 0)
    {
      // This direction will be divided in two
      refinedDx[idir] /= 2.0;
      maxIter[idir]    = 1;
    }
    else
    {
      // This direction is undivided
      maxIter[idir]    = 0;
    }

    // Initial the start to zero
    curIter[idir] = 0;
  }

  // Iterate through all the subcells
  while (1)
  {
    // Computes the cell center for the iCell refined cell
    IndexTM<Real,dim> refinedCellCenter = m_cutCellMoments.m_IFData.m_globalCoord.m_origin;

    for (int idir = 0; idir < dim; idir++)
    {
      if (a_refineInDir[idir] != 0)
      {
        int sign = 2*curIter[idir] - 1;
        refinedCellCenter[idir] += sign * 0.5*refinedDx[idir];
      }
    }

    // Creates the IF Data for the refined cell, it uses the global
    // constructor of IFData and recomputes all the intersection points
    // and normal/grad normal from the implicit function
#if RECURSIVE_GEOMETRY_GENERATION == 0
    IFData<dim> refinedIFData(m_cutCellMoments.m_IFData.m_function,refinedDx,refinedCellCenter,a_order);
#else
    IFData<dim> refinedIFData(m_cutCellMoments.m_IFData.m_function,refinedDx,refinedCellCenter,a_orderPmax);
#endif
#if 0
    pout() << "refinedIFData" << "\n";
    pout() << refinedIFData<<endl;
#endif
    // Creates a cutcellmoments for the refined cell
    ComputeCutCellMoments<dim> refinedComputeCCMom(refinedIFData);

    // Track current refinement level
    int numberOfRefinements = a_numberOfRefinements + 1;

    // Compute the moments for the refined cell
#if RECURSIVE_GEOMETRY_GENERATION == 0
    refinedComputeCCMom.computeMoments(a_order,a_degreeP,a_useConstraints,a_refinementCriterion,numberOfRefinements);
#else
    refinedComputeCCMom.computeMoments(a_orderPmax,a_degreePmax,a_useConstraints,a_refinementCriterion,numberOfRefinements);
#endif
    refinedCutCellMoments.push_back(refinedComputeCCMom.m_cutCellMoments);

    // Move to the next subcell (if there is one)
    int idir;
    for (idir = 0; idir < dim; idir++)
    {
      if (curIter[idir] < maxIter[idir])
      {
        // If the current index is in range, increment it and continue
        curIter[idir]++;
        break;
      }
      else
      {
        // If the current index at the maximum, reset it and move to the next
        // index
        curIter[idir] = 0;
      }
    }

    // All the subcells have been done
    if (idir == dim)
    {
      break;
    }
  }

  return refinedCutCellMoments;
}

template <int dim> void ComputeCutCellMoments<dim>::addMomentMaps(const Vector<CutCellMoments<dim> > & a_refinedCutCellVector,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                                                                  const int                          & a_degreeP,
#else
                                                                  const int                          & a_degreePmax,
#endif
                                                                  const bool                         & a_useConstraints)
{
  int numberOfRefinedCC = a_refinedCutCellVector.size();
  CH_assert(numberOfRefinedCC>0);

  IndexTM<int,2> bdId;

  // Pick one refined CutCellMoments and use its maps to initialize the maps
  // of the current CutCellMoment to zero the chosen cutcellmoment need to
  // have initialized boundary moments, i.e., be irregular
  CutCellMoments<dim> cutCell;
  int k = 0;

  while (k < numberOfRefinedCC && (a_refinedCutCellVector[k].m_IFData.m_allVerticesOut ||
                                  (a_refinedCutCellVector[k].m_IFData.m_allVerticesIn && !m_cutCellMoments.m_bdCCOn)))
  {
    k++;
  }

  if (k >= numberOfRefinedCC)
  {
    MayDay::Abort("Exceeded the number of refined ComputeCutCellMoments in vector");
  }

  else
  {
    cutCell = a_refinedCutCellVector[k];
  }

  m_cutCellMoments.initialize(cutCell);

  // Add up the moments on the refined cells to obtain the moment in the
  // coarse cell. Move moments to the m_parentCoord, then add.
  for (int i = 0; i < numberOfRefinedCC; i++)
  {
    // The difference between the refined cell center and the current cell
    // center is needed to change the coordinate system of the refined moments
    // to the coordinate system where they are needed.
    CutCellMoments<dim> refinedCutCell = a_refinedCutCellVector[i];
    IndexTM<int,dim> hilo = LARGEINTVAL*IndexTM<int,dim>::Unit;
    IndexTM<Real,dim> refinedCenterDelta = IndexTM<Real,dim>::Zero;

    int ii = i;
    for (int j = 0; j < dim; ++j)
    {
      hilo[j] = (ii & 1);
      ii = ii >> 1;
    }

    IndexTM<int,dim> sign = 2*hilo-1;

    for (int idir = 0; idir < dim; idir++)
    {
      refinedCenterDelta[idir] = 0.25*sign[idir]*m_cutCellMoments.m_IFData.m_globalCoord.m_dx[idir];
    }

    // Add up the moments in the volume
    addMoments(m_cutCellMoments.m_moments,refinedCutCell.m_moments,refinedCenterDelta);
    addMoments(m_cutCellMoments.m_EBmoments,refinedCutCell.m_EBmoments,refinedCenterDelta);

    // Add up the moments on the boundary
    for (int idir = 0; idir < dim; idir++)
    {
      bdId[0] = idir;

      // The difference in cell centers in dim-1 is needed to change the
      // coordinate system of the refined moments to the coordinate
      // system where they are needed.
      IndexTM<Real,dim-1> localRefinedCenterDelta;
      IndexTM<int,dim-1> localHilo;

      for (int jdir = 0; jdir < dim; jdir++)
      {
        if (jdir < idir)
        {
          localRefinedCenterDelta[jdir] = refinedCenterDelta[jdir];
          localHilo[jdir] = hilo[jdir];
        }
        else if (jdir > idir)
        {
          localRefinedCenterDelta[jdir-1] = refinedCenterDelta[jdir];
          localHilo[jdir-1] = hilo[jdir];
        }
      }

      if (hilo[idir] != LARGEINTVAL)
      {
        bdId[1] = hilo[idir];

        // Drop one dimension to add up the boundary moments dropping
        // one dimension allows an easier implementation as 1D is a
        // particular case
        (refinedCutCell.m_bdCutCellMoments[bdId]).addBdMoments(
            m_cutCellMoments.m_bdCutCellMoments[bdId],
            refinedCutCell.m_IFData,
#if RECURSIVE_GEOMETRY_GENERATION == 0
            a_degreeP,
#else
            a_degreePmax,
#endif
            a_useConstraints,
            localRefinedCenterDelta,
            localHilo);
      }
      else
      {
        for (int side = 0; side < 2; side++)
        {
          bdId[1] = side;
          (refinedCutCell.m_bdCutCellMoments[bdId]).addBdMoments(
              m_cutCellMoments.m_bdCutCellMoments[bdId],
              refinedCutCell.m_IFData,
#if RECURSIVE_GEOMETRY_GENERATION == 0
              a_degreeP+1,
#else
              a_degreePmax+1,
#endif
              a_useConstraints,
              localRefinedCenterDelta,
              localHilo);
        }
      }
    }
  }
}

template <int dim> void ComputeCutCellMoments<dim>::addMoments(PthMoment               & a_momentMap,
                                                               PthMoment               & a_refinedMomentMap,
                                                               const IndexTM<Real,dim> & a_refinedCenterDelta)
{
  // Iterate through the REFINEDmap and add the refined moment in the right
  // coordinate system to the coarse map
  for (typename PthMoment::const_iterator it = a_refinedMomentMap.begin(); it != a_refinedMomentMap.end(); ++it)
  {
    IndexTM<int,dim> monomial = it->first;
    Real moment = m_cutCellMoments.changeMomentCoordinates(a_refinedMomentMap,monomial,a_refinedCenterDelta);

    a_momentMap[it->first] += moment;
  }
}

#if RECURSIVE_GEOMETRY_GENERATION == 0
template <int dim> void ComputeCutCellMoments<dim>::computeResiduals(const int  & a_order,
                                                                     const int  & a_degreeP,
#else
template <int dim> void ComputeCutCellMoments<dim>::computeResiduals(const int  & a_orderPmax,
                                                                     const int  & a_degreePmax,
#endif
                                                                     const bool & a_useConstraints)
{
#if RECURSIVE_GEOMETRY_GENERATION == 0
  for (int iDegree = a_degreeP; iDegree >= 0; iDegree--)
#else
  for (int iDegree = a_degreePmax; iDegree >= 0; iDegree--)
#endif
  {
    // Initialize residuals to 0
    int nNorm = 3;
    m_cutCellMoments.m_residual[iDegree].resize(nNorm);

    for (int i = 0; i < nNorm; i++)
    {
      m_cutCellMoments.setResidual(0.0,iDegree,i);
    }

    // Create a lsp problem
    IvDim zeroDerivative = IndexTM<int,dim>::Zero;
#if RECURSIVE_GEOMETRY_GENERATION == 0
    LSProblem<dim> lsp(a_order,iDegree,a_useConstraints,m_cutCellMoments.m_IFData.m_normalDerivatives[zeroDerivative]);
#else
    LSProblem<dim> lsp(a_orderPmax,iDegree,a_useConstraints,m_cutCellMoments.m_IFData.m_normalDerivatives[zeroDerivative]);
#endif


    // Compute the right hand side
#if RECURSIVE_GEOMETRY_GENERATION == 0
    Vector<Real> rhs = computeRhs(lsp,a_order);
#else
    Vector<Real> rhs = computeRhs(lsp,a_orderPmax);
#endif

    // Get the unknowns from the moment maps
    Vector<Real> unknowns(lsp.m_numP+lsp.m_numPLess1);

    for (typename PthMomentLoc::const_iterator it = lsp.m_monoLocP.begin(); it != lsp.m_monoLocP.end(); ++it)
    {
      unknowns[it->second] = m_cutCellMoments.m_EBmoments[it->first];
    }

    for (typename PthMomentLoc::const_iterator it = lsp.m_monoLocPLess1.begin(); it != lsp.m_monoLocPLess1.end(); ++it)
    {
      unknowns[it->second + lsp.m_numP] = m_cutCellMoments.m_moments[it->first];
    }

    // Compute the residuals
    Real maxRi = 0.0;

    for (int i = 0; i < lsp.m_numP*dim; i++)
    {
      Real AtimeX = 0.0;

      for (int j = 0; j < lsp.m_numP + lsp.m_numPLess1; j++)
      {
        AtimeX += lsp.m_matrix[i][j] * unknowns[j];
      }

      Real ri = AtimeX - rhs[i];

      if (Abs(ri) > maxRi)
      {
        m_cutCellMoments.setResidual(Abs(ri),iDegree,0);
        maxRi = Abs(ri);
      }

      m_cutCellMoments.setResidual(m_cutCellMoments.getResidual(iDegree,1) + Abs(ri),iDegree,1);
      m_cutCellMoments.setResidual(m_cutCellMoments.getResidual(iDegree,2) + ri * ri,iDegree,2);
    }

    m_cutCellMoments.setResidual(sqrt(m_cutCellMoments.getResidual(iDegree,2)),iDegree,2);
  }
}

template <int dim> void ComputeCutCellMoments<dim>::computeResiduals(const Vector< CutCellMoments<dim> > & a_refinedCCMoms,
#if RECURSIVE_GEOMETRY_GENERATION == 0
                                                                     const int                           & a_degreeP)
#else
                                                                     const int                           & a_degreePmax)
#endif
{
#if RECURSIVE_GEOMETRY_GENERATION == 0
  for (int iDegree = 0; iDegree <= a_degreeP; iDegree++)
#else
  for (int iDegree = 0; iDegree <= a_degreePmax; iDegree++)
#endif
  {
    Real maxRes = 0.0;
    Real tempRes2 = 0.0;
    Real tempRes1 = 0.0;

    for (int i = 0; i < a_refinedCCMoms.size(); i++)
    {
      Real res0 = a_refinedCCMoms[i].getResidual(iDegree,0);
      Real res1 = a_refinedCCMoms[i].getResidual(iDegree,1);
      Real res2 = a_refinedCCMoms[i].getResidual(iDegree,2);

      if (res0 > maxRes && res0 != LARGEREALVAL)
      {
        maxRes = res0;
        m_cutCellMoments.setResidual(res0,iDegree,0);
      }

      if (res1 != LARGEREALVAL)
      {
        tempRes1 += res1;
      }

      if (res2 != LARGEREALVAL)
      {
        tempRes2 += res2*res2;
      }
    }

    m_cutCellMoments.setResidual(tempRes1,iDegree,1);

    if (tempRes2 != LARGEREALVAL)
    {
      m_cutCellMoments.setResidual(sqrt(tempRes2),iDegree,2);
    }
  }
}

template <int dim> void ComputeCutCellMoments<dim>::print(ostream& a_out) const
{
  m_cutCellMoments.print(a_out);
}

template <int dim> void ComputeCutCellMoments<dim>::dump() const
{
  print(pout());
}

// Operators
template <int dim> void ComputeCutCellMoments<dim>::operator=(const ComputeCutCellMoments<dim> & a_computeCutCellMoments)
{
  // Only copy if the objects are distinct
  if (this != &a_computeCutCellMoments)
  {
    m_cutCellMoments = a_computeCutCellMoments.m_cutCellMoments;
  }
}

template <int dim> ostream& operator<<(ostream                   & a_out,
                                       const ComputeCutCellMoments<dim> & a_computeCutCellMoments)
{
  a_computeCutCellMoments.print(a_out);
  return a_out;
}

template <int dim> Real ComputeCutCellMoments<dim>::factorial(const IvDim & a_multiIndex) const
{
  Real fact = 1.0;

  for (int i = 0; i < dim; i++)
  {
    for (int j = 2; j <= a_multiIndex[i]; j++)
    {
      fact *= j;
    }
  }

  return fact;
}

#include "NamespaceFooter.H"

#endif
