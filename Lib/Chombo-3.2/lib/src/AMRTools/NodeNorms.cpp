#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Norms.cpp
// petermc, 12 Aug 2002
// petermc, 9 Sep 2002, revised again:
// For L_p norms where p is finite, use trapezoidal rule.
// petermc, 28 Mar 2003, revised again:
// Use only valid nodes.
// For L_p norms where p is finite, weight by dx^SpaceDim.
// petermc, 17 Oct 2005, put the integral functions into NodeIntegrals.cpp

#include <cmath>
#include <iostream>
#include "NodeNorms.H"
#include "NodeSetOperations.H"
#include "NodeDotProdF_F.H"
#include "NamespaceHeader.H"

using std::cout;
using std::cerr;
using std::endl;

// ---------------------------------------------------------
// added 11 June 2003, removed from NodeFArrayBox
Real
norm(const NodeFArrayBox& a_nfab,
     const Real a_dx,
     const Box& a_subbox, // CELL-centered
     const int a_p,
     const int a_startComp,
     const int a_numComp)
{
  if (a_p == 0)
    return maxnorm(a_nfab, a_subbox, a_startComp, a_numComp);

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

  Real normTotal;
  FORT_NODENORM(CHF_REAL(normTotal),
                CHF_CONST_INT(a_p),
                CHF_CONST_FRA(fab),
                CHF_CONST_FRA1(weight, 0),
                CHF_BOX(nodeBox),
                CHF_CONST_INT(a_startComp),
                CHF_CONST_INT(a_startComp + a_numComp-1));

  return normTotal;
}


// ---------------------------------------------------------
// added 11 June 2003, removed from NodeFArrayBox
Real
maxnorm(const NodeFArrayBox& a_nfab,
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

  Real normTotal;
  FORT_NODEMAXNORM(CHF_REAL(normTotal),
                   CHF_CONST_FRA(fab),
                   CHF_BOX(nodeBox),
                   CHF_CONST_INT(a_startComp),
                   CHF_CONST_INT(a_startComp + a_numComp-1));

  return normTotal;
}


// ---------------------------------------------------------
// 27 March 2003:
// This is called by other functions, and should not be called directly.
Real
norm(const BoxLayoutData<NodeFArrayBox>& a_layout,
     const Real a_dx,
     const int a_p,
     const Interval& a_interval,
     bool a_verbose)
{
  if (a_p == 0)
    return maxnorm(a_layout, a_interval, a_verbose);

  Real normTotal = 0.;

  for (DataIterator it = a_layout.dataIterator(); it.ok(); ++it)
    {
      const NodeFArrayBox& thisNfab = a_layout[it()];
      const Box& thisBox(a_layout.box(it())); // CELL-centered
      Real thisNfabNorm = norm(thisNfab, a_dx, thisBox, a_p,
                               a_interval.begin(), a_interval.size());
      if (a_verbose)
        cout << a_p << "norm(" << thisBox << ") = " << thisNfabNorm << endl;

      if (a_p == 1)
        {
          normTotal += thisNfabNorm;
        }
      else if (a_p == 2)
        {
          normTotal += thisNfabNorm * thisNfabNorm;
        }
      else
        {
          normTotal += pow(thisNfabNorm, Real(a_p));
        }
    }
# ifdef CH_MPI
  Real recv;
  // add up (a_p is not 0)
  int result = MPI_Allreduce(&normTotal, &recv, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on norm");
    }
  normTotal = recv;
# endif
  // now do sqrt, etc
  if (a_p == 2)
    normTotal = sqrt(normTotal);
  else
    if ((a_p != 0) && (a_p != 1))
      normTotal = pow(normTotal, (Real)1.0/Real(a_p));

  return normTotal;
}


// ---------------------------------------------------------
// 7 Dec 2005
Real
norm(const BoxLayoutData<NodeFArrayBox>& a_layout,
     const LevelData<NodeFArrayBox>& a_mask,
     const ProblemDomain& a_domain,
     const Real a_dx,
     const int a_p,
     const Interval& a_interval,
     bool a_verbose)
{
  if (a_p == 0)
    return maxnorm(a_layout, a_mask, a_domain, a_interval, a_verbose);

  Real normTotal = 0.;
  int ncomp = a_interval.size();
  Box domBox = a_domain.domainBox();
  for (DataIterator it = a_layout.dataIterator(); it.ok(); ++it)
    {
      const NodeFArrayBox& thisNfab = a_layout[it()];
      const FArrayBox& dataFab = thisNfab.getFab();
      const FArrayBox& maskFab = a_mask[it()].getFab();
      const Box& thisBox(a_layout.box(it())); // CELL-centered
      NodeFArrayBox dataMasked(thisBox, ncomp);
      FArrayBox& dataMaskedFab = dataMasked.getFab();
      dataMaskedFab.copy(dataFab);
      // dataMaskedFab *= maskFab;
      for (int comp = a_interval.begin(); comp <= a_interval.end(); comp++)
        {
          // Set dataMaskedFab[comp] *= maskFab[0].
          dataMaskedFab.mult(maskFab, 0, comp);
        }

      Real thisNfabNorm = 0.;
      if (thisBox.intersects(domBox))
        {
          Box thisBoxInDomain = thisBox & domBox;
          thisNfabNorm = norm(dataMasked, a_dx, thisBoxInDomain, a_p,
                              a_interval.begin(), a_interval.size());
        }

      if (a_verbose)
        cout << a_p << "norm(" << thisBox << ") = " << thisNfabNorm << endl;

      if (a_p == 1)
        {
          normTotal += thisNfabNorm;
        }
      else if (a_p == 2)
        {
          normTotal += thisNfabNorm * thisNfabNorm;
        }
      else
        {
          normTotal += pow(thisNfabNorm, Real(a_p));
        }
    }
# ifdef CH_MPI
  Real recv;
  // add up (a_p is not 0)
  int result = MPI_Allreduce(&normTotal, &recv, 1, MPI_CH_REAL,
                             MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on norm");
    }
  normTotal = recv;
# endif
  // now do sqrt, etc
  if (a_p == 2)
    normTotal = sqrt(normTotal);
  else
    if ((a_p != 0) && (a_p != 1))
      normTotal = pow(normTotal, (Real)1.0/Real(a_p));

  return normTotal;
}


// ---------------------------------------------------------
// 7 Dec 2005
Real
maxnorm(const BoxLayoutData<NodeFArrayBox>& a_layout,
        const LevelData<NodeFArrayBox>& a_mask,
        const ProblemDomain& a_domain,
        const Interval& a_interval,
        bool a_verbose)
{
  Real normTotal = 0.;
  // a_p == 0:  max norm
  int ncomp = a_interval.size();
  for (DataIterator it = a_layout.dataIterator(); it.ok(); ++it)
    {
      const NodeFArrayBox& thisNfab = a_layout[it()];
      const FArrayBox& dataFab = thisNfab.getFab();
      const FArrayBox& maskFab = a_mask[it()].getFab();
      const Box& thisBox(a_layout.box(it())); // CELL-centered
      NodeFArrayBox dataMasked(thisBox, ncomp);
      FArrayBox& dataMaskedFab = dataMasked.getFab();
      dataMaskedFab.copy(dataFab);
      // dataMaskedFab *= maskFab;
      for (int comp = a_interval.begin(); comp <= a_interval.end(); comp++)
        {
          // Set dataMaskedFab[comp] *= maskFab[0].
          dataMaskedFab.mult(maskFab, 0, comp);
        }
      Real thisNfabNorm =
        maxnorm(dataMasked, thisBox, a_interval.begin(), a_interval.size());
      if (a_verbose)
        cout << "maxnorm(" << thisBox << ") = " << thisNfabNorm << endl;
       normTotal = Max(normTotal, thisNfabNorm);
    }
# ifdef CH_MPI
  Real recv;
  // add up (a_p is not 0)
  int result = MPI_Allreduce(&normTotal, &recv, 1, MPI_CH_REAL,
                             MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on norm");
    }
  normTotal = recv;
# endif
  return normTotal;
}


// ---------------------------------------------------------
// 27 March 2003:
// This is called by other functions, and should not be called directly.
Real
maxnorm(const BoxLayoutData<NodeFArrayBox>& a_layout,
        const Interval& a_interval,
        bool a_verbose)
{
  Real normTotal = 0.;
  // a_p == 0:  max norm
  for (DataIterator it = a_layout.dataIterator(); it.ok(); ++it)
    {
      const Box& thisBox(a_layout.box(it())); // CELL-centered
      const NodeFArrayBox& thisNfab = a_layout[it()];
      Real thisNfabNorm =
        maxnorm(thisNfab, thisBox, a_interval.begin(), a_interval.size());
      if (a_verbose)
        cout << "maxnorm(" << thisBox << ") = " << thisNfabNorm << endl;
      normTotal = Max(normTotal, thisNfabNorm);
    }
# ifdef CH_MPI
  Real recv;
  int result = MPI_Allreduce(&normTotal, &recv, 1, MPI_CH_REAL,
                             MPI_MAX, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on maxnorm");
    }
  normTotal = recv;
# endif
  return normTotal;
}


// --------------------------------------------------------------
Real
norm(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
     const Vector<Box>& a_domain,
     const Vector<int>& a_nRefFine,
     const Real a_dxCrse,
     const Interval& a_comps,
     const int a_p,
     const int a_lBase,
     bool a_verbose)
{
  int numlevels = a_domain.size();
  Vector<ProblemDomain> probdomain(numlevels);
  for (int ilev = 0; ilev < numlevels; ilev++)
    probdomain[ilev] = ProblemDomain(a_domain[ilev]);

  return norm(a_phi, probdomain, a_nRefFine,
              a_dxCrse, a_comps, a_p, a_lBase, a_verbose);
}


// --------------------------------------------------------------
// version of 27 March 2003
Real
norm(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
     const Vector<ProblemDomain>& a_domain,
     const Vector<int>& a_nRefFine,
     const Real a_dxCrse,
     const Interval& a_comps,
     const int a_p,
     const int a_lBase,
     bool a_verbose)
{
  // deleted 27 March 2003
  // if (a_p == 0)
  // return maxnorm(a_phi, a_domain, a_nRefFine, a_comps, a_lBase, a_verbose);

  Real dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real normTotal = 0.0;

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
          Real normLevel = norm(thisPhi, levDomain, finerGridsPtr,
                                a_nRefFine[lev], dxLevel, a_comps, a_p, a_verbose);

          if (a_p == 0)
            {
              // max norm:  reset if new max found at this level
              if (normLevel > normTotal)
                normTotal = normLevel;
            }
          else if (a_p == 1)
            normTotal += normLevel;
          else if (a_p == 2)
            normTotal += normLevel * normLevel;
          else
            normTotal += pow(normLevel, Real(a_p));
        }

      // update  dxLevel
      dxLevel = dxLevel / Real(a_nRefFine[lev]);
    }

  // shouldn't need to do broadcast/gather thing

  // now do sqrt, etc
  if (a_p == 2)
    normTotal = sqrt(normTotal);
  else
    if ((a_p != 0) && (a_p != 1))
      normTotal = pow(normTotal, (Real)1.0/Real(a_p));

  return normTotal;
}


// --------------------------------------------------------------
// New, 17 Nov 2005
Real
norm(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
     const Vector<ProblemDomain>& a_domain,
     const Vector<LayoutData< Vector<Box> >* >& a_IVSVext,
     const Vector<LayoutData< Vector<Box> >* >& a_IVSVintFinerCoarsened,
     const Vector<DisjointBoxLayout>& a_layoutsFinerCoarsened,
     const Vector<int>& a_nRefFine,
     const Real a_dxCrse,
     const Interval& a_comps,
     const int a_p,
     const int a_lBase,
     bool a_verbose)
{
  Real dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real normTotal = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
    {
      // in case there are extra levels which are not defined
      if (a_phi[lev] != NULL)
        {
          CH_assert(a_phi[lev]->isDefined());

          LevelData<NodeFArrayBox>& thisPhi = *(a_phi[lev]);

          const LayoutData< Vector<Box> >& exterior = *(a_IVSVext[lev]);

          const ProblemDomain levDomain = a_domain[lev];
          Real normLevel;
          if (lev < numLevels-1)
            {
              normLevel = norm(thisPhi, levDomain, a_layoutsFinerCoarsened[lev],
                               exterior, *(a_IVSVintFinerCoarsened[lev]),
                               a_nRefFine[lev], dxLevel, a_comps, a_p, a_verbose);
            }
          else // max level
            {
              normLevel = norm(thisPhi, levDomain,
                               exterior,
                               dxLevel, a_comps, a_p, a_verbose);
            }

          if (a_p == 0)
            {
              // max norm:  reset if new max found at this level
              if (normLevel > normTotal)
                normTotal = normLevel;
            }
          else if (a_p == 1)
            normTotal += normLevel;
          else if (a_p == 2)
            normTotal += normLevel * normLevel;
          else
            normTotal += pow(normLevel, Real(a_p));
        }

      // update  dxLevel
      dxLevel = dxLevel / Real(a_nRefFine[lev]);
    }

  // shouldn't need to do broadcast/gather thing

  // now do sqrt, etc
  if (a_p == 2)
    normTotal = sqrt(normTotal);
  else
    if ((a_p != 0) && (a_p != 1))
      normTotal = pow(normTotal, (Real)1.0/Real(a_p));

  return normTotal;
}


// --------------------------------------------------------------
// New, 18 Nov 2005
Real
norm(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
     const Vector<ProblemDomain>& a_domain,
     const Vector<LayoutData< Vector<IntVectSet> >* >& a_IVSVext,
     const Vector<LayoutData< Vector<IntVectSet> >* >& a_IVSVintFinerCoarsened,
     const Vector<DisjointBoxLayout>& a_layoutsFinerCoarsened,
     const Vector<int>& a_nRefFine,
     const Real a_dxCrse,
     const Interval& a_comps,
     const int a_p,
     const int a_lBase,
     bool a_verbose)
{
  Real dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real normTotal = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
    {
      // in case there are extra levels which are not defined
      if (a_phi[lev] != NULL)
        {
          CH_assert(a_phi[lev]->isDefined());

          LevelData<NodeFArrayBox>& thisPhi = *(a_phi[lev]);

          const LayoutData< Vector<IntVectSet> >& exterior = *(a_IVSVext[lev]);

          const ProblemDomain levDomain = a_domain[lev];
          Real normLevel;
          if (lev < numLevels-1)
            {
              normLevel = norm(thisPhi, levDomain, a_layoutsFinerCoarsened[lev],
                               exterior, *(a_IVSVintFinerCoarsened[lev]),
                               a_nRefFine[lev], dxLevel, a_comps, a_p, a_verbose);
            }
          else // max level
            {
              normLevel = norm(thisPhi, levDomain,
                               exterior,
                               dxLevel, a_comps, a_p, a_verbose);
            }

          if (a_p == 0)
            {
              // max norm:  reset if new max found at this level
              if (normLevel > normTotal)
                normTotal = normLevel;
            }
          else if (a_p == 1)
            normTotal += normLevel;
          else if (a_p == 2)
            normTotal += normLevel * normLevel;
          else
            normTotal += pow(normLevel, Real(a_p));
        }

      // update  dxLevel
      dxLevel = dxLevel / Real(a_nRefFine[lev]);
    }

  // shouldn't need to do broadcast/gather thing

  // now do sqrt, etc
  if (a_p == 2)
    normTotal = sqrt(normTotal);
  else
    if ((a_p != 0) && (a_p != 1))
      normTotal = pow(normTotal, (Real)1.0/Real(a_p));

  return normTotal;
}


// --------------------------------------------------------------
// New, 7 Dec 2005
Real
norm(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
     const Vector<LevelData<NodeFArrayBox>* >& a_mask,
     const Vector<ProblemDomain>& a_vectPD,
     const Vector<int>& a_nRefFine,
     const Real a_dxCrse,
     const Interval& a_comps,
     const int a_p,
     const int a_lBase,
     bool a_verbose)
{
  Real dxLevel = a_dxCrse;
  int numLevels = a_phi.size();
  Real normTotal = 0.0;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++)
    {
      // in case there are extra levels which are not defined
      if (a_phi[lev] != NULL)
        {
          CH_assert(a_phi[lev]->isDefined());
          CH_assert(a_mask[lev]->isDefined());

          Real normLevel = norm(*(a_phi[lev]), *(a_mask[lev]), a_vectPD[lev],
                                dxLevel, a_p, a_comps, a_verbose);

          if (a_p == 0)
            {
              // max norm:  reset if new max found at this level
              if (normLevel > normTotal)
                normTotal = normLevel;
            }
          else if (a_p == 1)
            normTotal += normLevel;
          else if (a_p == 2)
            normTotal += normLevel * normLevel;
          else
            normTotal += pow(normLevel, Real(a_p));
        }

      // update dxLevel if there is a next finer level
      if (lev < numLevels - 1) dxLevel = dxLevel / Real(a_nRefFine[lev]);
    }

  // shouldn't need to do broadcast/gather thing

  // now do sqrt, etc
  if (a_p == 2)
    normTotal = sqrt(normTotal);
  else
    if ((a_p != 0) && (a_p != 1))
      normTotal = pow(normTotal, (Real)1.0/Real(a_p));

  return normTotal;
}


// ------------------------------------------------------------
// version of 28 March 2003
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const Box& a_domain,
     const DisjointBoxLayout* a_finerGridsPtr,
     const int a_nRefFine,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
     bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return norm(a_phi, probdomain, a_finerGridsPtr,
              a_nRefFine, a_dx, a_comps, a_p, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const ProblemDomain& a_domain,
     const DisjointBoxLayout* a_finerGridsPtr,
     const int a_nRefFine,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
     bool a_verbose)
{
  const DisjointBoxLayout& grids = a_phi.getBoxes();

  LayoutData< Vector<IntVectSet> > IVSVint, IVSVext;
  // LayoutData< Vector<Box> > IVSVint, IVSVext;
  interiorBoundaryNodes(IVSVint, grids, a_domain);
  exteriorBoundaryNodes(IVSVext, IVSVint, grids);

  Real normLevel = 0.;
  if (a_finerGridsPtr == NULL)
    {
      normLevel = norm(a_phi, a_domain,
                       IVSVext,
                       a_dx, a_comps, a_p, a_verbose);
    }
  else
    {
      DisjointBoxLayout coarsenedFinerGrids;
      coarsen(coarsenedFinerGrids, *a_finerGridsPtr, a_nRefFine);
      LayoutData< Vector<IntVectSet> > IVSVintFinerCoarsened;
      // LayoutData< Vector<Box> > IVSVintFinerCoarsened;
      interiorBoundaryNodes(IVSVintFinerCoarsened, grids,
                            coarsenedFinerGrids, a_domain);

      normLevel = norm(a_phi, a_domain, coarsenedFinerGrids,
                       IVSVext, IVSVintFinerCoarsened,
                       a_nRefFine, a_dx, a_comps, a_p, a_verbose);
    }

  return normLevel;
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const Box& a_domain,
     const DisjointBoxLayout& a_finerGridsCoarsened,
     const LayoutData< Vector<IntVectSet> >& a_IVSVext,
     const LayoutData< Vector<IntVectSet> >& a_IVSVintFinerCoarsened,
     const int a_nRefFine,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
     bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return norm(a_phi, probdomain, a_finerGridsCoarsened,
              a_IVSVext, a_IVSVintFinerCoarsened,
              a_nRefFine, a_dx, a_comps, a_p, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const Box& a_domain,
     const DisjointBoxLayout& a_finerGridsCoarsened,
     const LayoutData< Vector<Box> >& a_IVSVext,
     const LayoutData< Vector<Box> >& a_IVSVintFinerCoarsened,
     const int a_nRefFine,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
     bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return norm(a_phi, probdomain, a_finerGridsCoarsened,
              a_IVSVext, a_IVSVintFinerCoarsened,
              a_nRefFine, a_dx, a_comps, a_p, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const ProblemDomain& a_domain,
     const DisjointBoxLayout& a_finerGridsCoarsened,
     const LayoutData< Vector<IntVectSet> >& a_IVSVext,
     const LayoutData< Vector<IntVectSet> >& a_IVSVintFinerCoarsened,
     const int a_nRefFine,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
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

  Real normLevel = norm(temp, a_dx, a_p, newcomps, a_verbose);

  return normLevel;
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const ProblemDomain& a_domain,
     const DisjointBoxLayout& a_finerGridsCoarsened,
     const LayoutData< Vector<Box> >& a_IVSVext,
     const LayoutData< Vector<Box> >& a_IVSVintFinerCoarsened,
     const int a_nRefFine,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
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

  Real normLevel = norm(temp, a_dx, a_p, newcomps, a_verbose);

  return normLevel;
}


// ------------------------------------------------------------
// version of 27 March 2003:  no finer level
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const Box& a_domain,
     const LayoutData< Vector<IntVectSet> >& a_IVSVext,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
     bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return norm(a_phi, probdomain,
              a_IVSVext,
              a_dx, a_comps, a_p, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003:  no finer level
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const Box& a_domain,
     const LayoutData< Vector<Box> >& a_IVSVext,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
     bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return norm(a_phi, probdomain,
              a_IVSVext,
              a_dx, a_comps, a_p, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003:  no finer level
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const ProblemDomain& a_domain,
     const LayoutData< Vector<IntVectSet> >& a_IVSVext,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
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

  Real normLevel = norm(temp, a_dx, a_p, newcomps, a_verbose);

  return normLevel;
}


// ------------------------------------------------------------
// version of 27 March 2003:  no finer level
Real
norm(const LevelData<NodeFArrayBox>& a_phi,
     const ProblemDomain& a_domain,
     const LayoutData< Vector<Box> >& a_IVSVext,
     const Real a_dx,
     const Interval& a_comps,
     const int a_p,
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

  Real normLevel = norm(temp, a_dx, a_p, newcomps, a_verbose);

  return normLevel;
}


// --------------------------------------------------------------
Real
maxnorm(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
        const Vector<Box>& a_domain,
        const Vector<int>& a_nRefFine,
        const Interval& a_comps,
        const int a_lBase,
        bool a_verbose)
{
  int numlevels = a_domain.size();
  Vector<ProblemDomain> probdomain(numlevels);
  for (int ilev = 0; ilev < numlevels; ilev++)
    probdomain[ilev] = ProblemDomain(a_domain[ilev]);

  Real normAll = maxnorm(a_phi, probdomain, a_nRefFine,
                         a_comps, a_lBase, a_verbose);

  return normAll;
}


// --------------------------------------------------------------
Real
maxnorm(const Vector<LevelData<NodeFArrayBox>* >& a_phi,
        const Vector<ProblemDomain>& a_domain,
        const Vector<int>& a_nRefFine,
        const Interval& a_comps,
        const int a_lBase,
        bool a_verbose)
{
  Real dx = 1.;  // arbitrary dx isn't used anyway
  Real normAll = norm(a_phi, a_domain, a_nRefFine,
                      dx, a_comps, 0, a_lBase, a_verbose);
  return normAll;
}


// ------------------------------------------------------------
// version of 28 March 2003
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const Box& a_domain,
        const DisjointBoxLayout* a_finerGridsPtr,
        const int a_nRefFine,
        const Interval& a_comps,
        bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return maxnorm(a_phi, probdomain, a_finerGridsPtr,
                 a_nRefFine, a_comps, a_verbose);
}


// ------------------------------------------------------------
// version of 27 March 2003
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const ProblemDomain& a_domain,
        const DisjointBoxLayout* a_finerGridsPtr,
        const int a_nRefFine,
        const Interval& a_comps,
        bool a_verbose)
{
  Real dx = 1.;  // arbitrary dx isn't used anyway
  Real normLevel = norm(a_phi, a_domain, a_finerGridsPtr,
                        a_nRefFine, dx, a_comps, 0, a_verbose);
  return normLevel;
}


// ------------------------------------------------------------
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const ProblemDomain& a_domain,
        const DisjointBoxLayout& a_finerGridsCoarsened,
        const LayoutData< Vector<IntVectSet> >& a_IVSVext,
        const LayoutData< Vector<IntVectSet> >& a_IVSVintFinerCoarsened,
        const int a_nRefFine,
        const Interval& a_comps,
        bool a_verbose)
{
  Real dx = 1.;  // arbitrary dx isn't used anyway
  Real normLevel = norm(a_phi, a_domain, a_finerGridsCoarsened,
                        a_IVSVext, a_IVSVintFinerCoarsened,
                        a_nRefFine, dx, a_comps, 0, a_verbose);
  return normLevel;
}


// ------------------------------------------------------------
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const ProblemDomain& a_domain,
        const DisjointBoxLayout& a_finerGridsCoarsened,
        const LayoutData< Vector<Box> >& a_IVSVext,
        const LayoutData< Vector<Box> >& a_IVSVintFinerCoarsened,
        const int a_nRefFine,
        const Interval& a_comps,
        bool a_verbose)
{
  Real dx = 1.;  // arbitrary dx isn't used anyway
  Real normLevel = norm(a_phi, a_domain, a_finerGridsCoarsened,
                        a_IVSVext, a_IVSVintFinerCoarsened,
                        a_nRefFine, dx, a_comps, 0, a_verbose);
  return normLevel;
}


// ------------------------------------------------------------
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const Box& a_domain,
        const DisjointBoxLayout& a_finerGridsCoarsened,
        const LayoutData< Vector<Box> >& a_IVSVext,
        const LayoutData< Vector<Box> >& a_IVSVintFinerCoarsened,
        const int a_nRefFine,
        const Interval& a_comps,
        bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return maxnorm(a_phi, probdomain, a_finerGridsCoarsened,
                 a_IVSVext, a_IVSVintFinerCoarsened,
                 a_nRefFine, a_comps, a_verbose);
}


// ------------------------------------------------------------
// no finer level
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const ProblemDomain& a_domain,
        const LayoutData< Vector<IntVectSet> >& a_IVSVext,
        const Interval& a_comps,
        bool a_verbose)
{
  Real dx = 1.;  // arbitrary dx isn't used anyway
  Real normLevel = norm(a_phi, a_domain,
                        a_IVSVext,
                        dx, a_comps, 0, a_verbose);
  return normLevel;
}


// ------------------------------------------------------------
// no finer level
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const ProblemDomain& a_domain,
        const LayoutData< Vector<Box> >& a_IVSVext,
        const Interval& a_comps,
        bool a_verbose)
{
  Real dx = 1.;  // arbitrary dx isn't used anyway
  Real normLevel = norm(a_phi, a_domain,
                        a_IVSVext,
                        dx, a_comps, 0, a_verbose);
  return normLevel;
}


// ------------------------------------------------------------
// no finer level
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const Box& a_domain,
        const LayoutData< Vector<IntVectSet> >& a_IVSVext,
        const Interval& a_comps,
        bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return maxnorm(a_phi, probdomain,
                 a_IVSVext,
                 a_comps, a_verbose);
}


// ------------------------------------------------------------
// no finer level
Real
maxnorm(const LevelData<NodeFArrayBox>& a_phi,
        const Box& a_domain,
        const LayoutData< Vector<Box> >& a_IVSVext,
        const Interval& a_comps,
        bool a_verbose)
{
  ProblemDomain probdomain(a_domain);
  return maxnorm(a_phi, probdomain,
                 a_IVSVext,
                 a_comps, a_verbose);
}

#include "NamespaceFooter.H"
