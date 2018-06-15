#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeDotProduct.cpp
// petermc, 11 June 2003

#include "NodeDotProduct.H"
#include "NodeSetOperations.H"
#include "NodeDotProdF_F.H"
#include "NamespaceHeader.H"

// ---------------------------------------------------------
Real
DotProductNodes(const BoxLayoutData<NodeFArrayBox>& a_dataOne,
                const BoxLayoutData<NodeFArrayBox>& a_dataTwo,
                const BoxLayout& a_dblIn)
{
  Interval comps = a_dataOne.interval();
  return DotProductNodes(a_dataOne, a_dataTwo, a_dblIn, comps);
}


// ---------------------------------------------------------
Real
DotProductNodes(const BoxLayoutData<NodeFArrayBox>& a_dataOne,
                const BoxLayoutData<NodeFArrayBox>& a_dataTwo,
                const BoxLayout& a_dblIn,
                const Interval& a_comps)
{
  const int startComp = a_comps.begin();
  const int endComp = a_comps.end();
  CH_assert(a_dataOne.nComp() == a_dataTwo.nComp());
  Real rhodot = 0.0;
  //calculate the single-processor dot product
  DataIterator dit = a_dataOne.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      Box fabbox(surroundingNodes(a_dblIn.get(dit())));
      const FArrayBox& onefab = a_dataOne[dit()].getFab();
      const FArrayBox& twofab = a_dataTwo[dit()].getFab();
      CH_assert(onefab.box().contains(fabbox));
      CH_assert(twofab.box().contains(fabbox));

      // Weights for trapezoidal rule.
      FArrayBox weight(fabbox, 1);
      Real dv = 1.;
      FORT_TRAPWEIGHTS(CHF_FRA1(weight, 0),
                       CHF_BOX(fabbox),
                       CHF_CONST_REAL(dv));

      Real dotgrid = 0;
      // petermc, 5 Apr 2002, changed DOTPROD to DOTPRODUCT
      // petermc, 15 Apr 2002 changed parameters
      // petermc, 22 Oct 2002, corrected bug that was looping over
      // all components
      // petermc, 22 Oct 2002, changed to NODEDOTPRODUCT
      FORT_NODEDOTPRODUCT(CHF_REAL(dotgrid),
                          CHF_CONST_FRA(onefab),
                          CHF_CONST_FRA(twofab),
                          CHF_CONST_FRA1(weight, 0),
                          CHF_BOX(fabbox),
                          CHF_CONST_INT(startComp),
                          CHF_CONST_INT(endComp));

      rhodot += dotgrid;
    }

  // now for the multi-processor fandango

  //gather all the rhodots onto a vector and add them up
  int baseProc = 0;
  Vector<Real> dotVec;
  gather(dotVec, rhodot, baseProc);

  Real rhodotTot = 0.0;
  if (procID() == baseProc)
    {
      CH_assert(dotVec.size() == numProc());
      for (int ivec = 0; ivec < dotVec.size(); ivec++)
        {
          rhodotTot += dotVec[ivec];
        }
    }
  //broadcast the sum to all processors.
  broadcast(rhodotTot, baseProc);

  //return the total
  return rhodotTot;
}


// ---------------------------------------------------------
Real
DotProductNodes(const LevelData<NodeFArrayBox>& a_dataOne,
                const LevelData<NodeFArrayBox>& a_dataTwo,
                const Box& a_domain,
                const LayoutData< Vector<IntVectSet> >& a_IVSVext,
                const Interval& a_comps)
{
  ProblemDomain probdomain(a_domain);
  return DotProductNodes(a_dataOne, a_dataTwo, probdomain, a_IVSVext, a_comps);
}


// ---------------------------------------------------------
Real
DotProductNodes(const LevelData<NodeFArrayBox>& a_dataOne,
                const LevelData<NodeFArrayBox>& a_dataTwo,
                const ProblemDomain& a_domain,
                const LayoutData< Vector<IntVectSet> >& a_IVSVext,
                const Interval& a_comps)
{
  // a_dataOne and a_dataTwo must have the same layout.

  // Idea:  copy a_dataOne to tempOne, then zero out tempOne on
  // exterior nodes of grids at this level.
  // Do the same with a_dataTwo to tempTwo.
  int ncomps = a_comps.size();
  const DisjointBoxLayout& grids = a_dataOne.getBoxes();
  LevelData<NodeFArrayBox> tempOne(grids, ncomps);
  LevelData<NodeFArrayBox> tempTwo(grids, ncomps);

  // Copy a_dataOne to tempOne, a_dataTwo to tempTwo.
  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit)
    {
      tempOne[dit()].copy(a_dataOne[dit()]);
      tempTwo[dit()].copy(a_dataTwo[dit()]);
    }

  // Zero out tempOne and tempTwo on exterior nodes.
  zeroBoundaryNodes(tempOne, a_IVSVext);
  zeroBoundaryNodes(tempTwo, a_IVSVext);

  Real thisdot = DotProductNodes(tempOne, tempTwo, grids, a_comps);

  return thisdot;
}


// ---------------------------------------------------------
Real
DotProductNodes(const LevelData<NodeFArrayBox>& a_dataOne,
                const LevelData<NodeFArrayBox>& a_dataTwo,
                const Box& a_domain,
                const Interval& a_comps)
{
  ProblemDomain probdomain(a_domain);
  return DotProductNodes(a_dataOne, a_dataTwo, probdomain, a_comps);
}


// ---------------------------------------------------------
Real
DotProductNodes(const LevelData<NodeFArrayBox>& a_dataOne,
                const LevelData<NodeFArrayBox>& a_dataTwo,
                const ProblemDomain& a_domain,
                const Interval& a_comps)
{
  // a_dataOne and a_dataTwo must have the same layout.

  // Idea:  copy a_dataOne to tempOne, then zero out tempOne on
  // exterior nodes of grids at this level.
  // Do the same with a_dataTwo to tempTwo.
  const DisjointBoxLayout& grids = a_dataOne.getBoxes();

  LayoutData< Vector<IntVectSet> > IVSVint, IVSVext;
  interiorBoundaryNodes(IVSVint, grids, a_domain);
  exteriorBoundaryNodes(IVSVext, IVSVint, grids);

  Real thisdot = DotProductNodes(a_dataOne, a_dataTwo, a_domain,
                                 IVSVext, a_comps);

  return thisdot;
}

#include "NamespaceFooter.H"
