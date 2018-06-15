#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeIntegrals.cpp
// petermc, 17 Oct 2005
// Previously, these functions were in Norms.cpp.

#include <cmath>
#include <iostream>
#include "NodeIntegrals.H"
#include "NodeSetOperations.H"
#include "NodeDotProdF_F.H"
#include "NamespaceHeader.H"

using std::cout;
using std::cerr;
using std::endl;

// ---------------------------------------------------------
// added 11 June 2003, removed from NodeFArrayBox
Real
integral(const NodeFArrayBox& a_nfab,
         const Real a_dx,
         const Box& a_subbox, // CELL-centered
         const int a_startComp,
         const int a_numComp)
{
  const FArrayBox& fab = a_nfab.getFab();
  Box bx(a_nfab.box());

  int nFabComp = fab.nComp();
  CH_assert(a_startComp + a_numComp <= nFabComp);
  CH_assert(bx.contains(a_subbox));

  Box nodeBox(surroundingNodes(a_subbox));

  // Set dv = dx^SpaceDim
  Real dv = 1.;
  for (int idir = 0; idir < SpaceDim; idir++) dv *= a_dx;

  // Weights for trapezoidal rule of integration.
  FArrayBox weight(nodeBox, 1);
  FORT_TRAPWEIGHTS(CHF_FRA1(weight, 0),
                   CHF_BOX(nodeBox),
                   CHF_CONST_REAL(dv));

  Real integralTotal;
  FORT_NODEINTEGRAL(CHF_REAL(integralTotal),
                    CHF_CONST_FRA(fab),
                    CHF_CONST_FRA1(weight, 0),
                    CHF_BOX(nodeBox),
                    CHF_CONST_INT(a_startComp),
                    CHF_CONST_INT(a_startComp + a_numComp-1));

  return integralTotal;
}


// ---------------------------------------------------------
// added 11 June 2003, removed from NodeFArrayBox
Real
integral(const NodeFArrayBox& a_nfab,
         const Real a_dx,
         const int a_startComp,
         const int a_numComp)
{
  Box bx(a_nfab.box());
  return integral(a_nfab, a_dx, bx, a_startComp, a_numComp);
}


// ---------------------------------------------------------
// 28 March 2003:
// This is called by other functions, and should not be called directly.
Real
integral(const BoxLayoutData<NodeFArrayBox>& a_layout,
         const Real a_dx,
         const Interval& a_interval,
         bool a_verbose)
{
  Real integralTotal = 0.;

  for (DataIterator it = a_layout.dataIterator(); it.ok(); ++it)
    {
      const NodeFArrayBox& thisNfab = a_layout[it()];
      const Box& thisBox(a_layout.box(it())); // CELL-centered
      Real thisNfabIntegral =
        integral(thisNfab, a_dx, thisBox,
                 a_interval.begin(), a_interval.size());
      integralTotal += thisNfabIntegral;
    }
# ifdef CH_MPI
  Real recv;
  // add up
  int result = MPI_Allreduce(&integralTotal, &recv, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on integral");
    }
  integralTotal = recv;
# endif
  return integralTotal;
}


// --------------------------------------------------------------
// version of 28 March 2003
Real
integral(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
         const Vector<ProblemDomain>& a_domain,
         const Vector<int>& a_nRefFine,
         const Real a_dxCrse,
         const Interval a_comps,
         const int a_lBase,
         bool a_verbose)
{
  Real dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real integralTotal = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
    {
      // in case there are extra levels which are not defined
      if (a_phi[lev] != NULL)
        {
          CH_assert(a_phi[lev]->isDefined());

          LevelData<NodeFArrayBox>& thisPhi = *(a_phi[lev]);
          const DisjointBoxLayout* finerGridsPtr = NULL;

          if (lev < numLevels-1)
            finerGridsPtr = &(a_phi[lev+1]->getBoxes());

          const ProblemDomain levDomain = a_domain[lev];
          Real integralLevel = integral(thisPhi, levDomain, finerGridsPtr,
                                        a_nRefFine[lev], dxLevel, a_comps, a_verbose);

          integralTotal += integralLevel;
        }

      // update  dxLevel
      dxLevel = dxLevel / Real(a_nRefFine[lev]);
    }

  // shouldn't need to do broadcast/gather thing

  return integralTotal;
}


// --------------------------------------------------------------
Real
integral(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
         const Vector<Box>& a_domain,
         const Vector<int>& a_nRefFine,
         const Real a_dxCrse,
         const Interval a_comps,
         const int a_lBase,
         bool a_verbose)
{
  int numlevels = a_domain.size();
  Vector<ProblemDomain> probdomain(numlevels);
  for (int ilev = 0; ilev < numlevels; ilev++)
    probdomain[ilev] = ProblemDomain(a_domain[ilev]);

  return integral(a_phi, probdomain, a_nRefFine,
                  a_dxCrse, a_comps, a_lBase, a_verbose);
}


// ------------------------------------------------------------
// version of 28 March 2003
Real
integral(const LevelData<NodeFArrayBox>& a_phi,
         const Box& a_domain,
         const DisjointBoxLayout* a_finerGridsPtr,
         const int a_nRefFine,
         const Real a_dx,
         const Interval a_comps,
         bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return integral(a_phi, probdomain, a_finerGridsPtr,
                  a_nRefFine, a_dx, a_comps, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
integral(const LevelData<NodeFArrayBox>& a_phi,
         const ProblemDomain& a_domain,
         const DisjointBoxLayout* a_finerGridsPtr,
         const int a_nRefFine,
         const Real a_dx,
         const Interval a_comps,
         bool a_verbose)
{
  const DisjointBoxLayout& grids = a_phi.getBoxes();

  LayoutData< Vector<IntVectSet> > IVSVint, IVSVext;
  interiorBoundaryNodes(IVSVint, grids, a_domain);
  exteriorBoundaryNodes(IVSVext, IVSVint, grids);

  Real integralLevel = 0.;
  if (a_finerGridsPtr == NULL)
    {
      integralLevel = integral(a_phi, a_domain,
                       IVSVext,
                       a_dx, a_comps, a_verbose);
    }
  else
    {
      DisjointBoxLayout coarsenedFinerGrids;
      coarsen(coarsenedFinerGrids, *a_finerGridsPtr, a_nRefFine);
      LayoutData< Vector<IntVectSet> > IVSVintFinerCoarsened;
      interiorBoundaryNodes(IVSVintFinerCoarsened, grids,
                            coarsenedFinerGrids, a_domain);

      integralLevel = integral(a_phi, a_domain, coarsenedFinerGrids,
                               IVSVext, IVSVintFinerCoarsened,
                               a_nRefFine, a_dx, a_comps, a_verbose);
    }

  return integralLevel;
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
integral(const LevelData<NodeFArrayBox>& a_phi,
         const Box& a_domain,
         const DisjointBoxLayout& a_finerGridsCoarsened,
         const LayoutData< Vector<IntVectSet> >& a_IVSVext,
         const LayoutData< Vector<IntVectSet> >& a_IVSVintFinerCoarsened,
         const int a_nRefFine,
         const Real a_dx,
         const Interval a_comps,
         bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return integral(a_phi, probdomain, a_finerGridsCoarsened,
                  a_IVSVext, a_IVSVintFinerCoarsened,
                  a_nRefFine, a_dx, a_comps, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
integral(const LevelData<NodeFArrayBox>& a_phi,
         const ProblemDomain& a_domain,
         const DisjointBoxLayout& a_finerGridsCoarsened,
         const LayoutData< Vector<IntVectSet> >& a_IVSVext,
         const LayoutData< Vector<IntVectSet> >& a_IVSVintFinerCoarsened,
         const int a_nRefFine,
         const Real a_dx,
         const Interval a_comps,
         bool a_verbose)
{
  // Idea:  copy a_phi to temp, then zero out temp on:
  // - exterior nodes of grids at this level;
  // - projections of interior nodes of the finer grids.
  int ncomps = a_comps.size();
  const DisjointBoxLayout& grids = a_phi.getBoxes();
  LevelData<NodeFArrayBox> temp(grids, ncomps);

  // Copy a_phi to temp.
  Interval newcomps(0, ncomps-1);
  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit)
    {
      const NodeFArrayBox& nfab = a_phi[dit()];
      const Box& bx = grids.get(dit());
      temp[dit()].copy(bx, newcomps, bx, nfab, a_comps);
    }

  // Zero out temp on exterior nodes.
  zeroBoundaryNodes(temp, a_IVSVext);

  // Define zeroCoarsened to be all zero on the coarsened finer grids.
  LevelData<NodeFArrayBox>
    zeroCoarsened(a_finerGridsCoarsened, ncomps, IntVect::Zero);
  for (DataIterator dit(a_finerGridsCoarsened.dataIterator()); dit.ok(); ++dit)
    zeroCoarsened[dit()].getFab().setVal(0.);

  // Set temp to zero on interior nodes of coarsened finer grids.
  copyInteriorNodes(temp, zeroCoarsened, a_IVSVintFinerCoarsened);

  Real integralLevel = integral(temp, a_dx, newcomps, a_verbose);

  return integralLevel;
}


// ------------------------------------------------------------
// version of 27 March 2003:  no finer level
Real
integral(const LevelData<NodeFArrayBox>& a_phi,
         const Box& a_domain,
         const LayoutData< Vector<IntVectSet> >& a_IVSVext,
         const Real a_dx,
         const Interval a_comps,
         bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return integral(a_phi, probdomain,
                  a_IVSVext,
                  a_dx, a_comps, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003:  no finer level
Real
integral(const LevelData<NodeFArrayBox>& a_phi,
         const ProblemDomain& a_domain,
         const LayoutData< Vector<IntVectSet> >& a_IVSVext,
         const Real a_dx,
         const Interval a_comps,
         bool a_verbose)
{
  // Idea:  copy a_phi to temp, then zero out temp on
  // exterior nodes of grids at this level.
  int ncomps = a_comps.size();
  const DisjointBoxLayout& grids = a_phi.getBoxes();
  LevelData<NodeFArrayBox> temp(grids, ncomps);

  // Copy a_phi to temp.
  Interval newcomps(0, ncomps-1);
  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit)
    {
      const NodeFArrayBox& nfab = a_phi[dit()];
      const Box& bx = grids.get(dit());
      temp[dit()].copy(bx, newcomps, bx, nfab, a_comps);
    }

  // Zero out temp on exterior nodes.
  zeroBoundaryNodes(temp, a_IVSVext);

  Real integralLevel = integral(temp, a_dx, newcomps, a_verbose);

  return integralLevel;
}

#include "NamespaceFooter.H"
