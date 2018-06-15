#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Weds, July 21, 1999

#include "DotProdF_F.H"
#include "DotProduct.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

Real DotProduct(const BoxLayoutData<FArrayBox>& a_dataOne,
                const BoxLayoutData<FArrayBox>& a_dataTwo,
                const BoxLayout&                a_dblIn)
{
  Interval comps = a_dataOne.interval();
  return DotProduct(a_dataOne, a_dataTwo, a_dblIn, comps);
}

Real DotProduct(const BoxLayoutData<FArrayBox>& a_dataOne,
                const BoxLayoutData<FArrayBox>& a_dataTwo,
                const BoxLayout&                a_dblIn,
                const Interval&                 a_comps)
{
  Vector<Real> rhodot;
  rhodot.reserve(10);
  const int startcomp = a_comps.begin();
  const int endcomp = a_comps.end();
  //calculate the single-processor dot product
  DataIterator dit = a_dataOne.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    Box fabbox = a_dblIn.get(dit());
    const FArrayBox& onefab = a_dataOne[dit()];
    const FArrayBox& twofab = a_dataTwo[dit()];
    CH_assert(onefab.box().contains(fabbox));
    CH_assert(twofab.box().contains(fabbox));

    Real dotgrid = 0;
    FORT_DOTPRODUCT(CHF_REAL(dotgrid),
                    CHF_CONST_FRA(onefab),
                    CHF_CONST_FRA(twofab),
                    CHF_BOX(fabbox),
                    CHF_CONST_INT(startcomp),
                    CHF_CONST_INT(endcomp));

    rhodot.push_back(dotgrid);
  }

  // now for the multi-processor fandango

  //gather all the rhodots onto a vector and add them up
  int baseProc = 0;
  Vector<Vector<Real> >dotVec;
  gather(dotVec, rhodot, baseProc);

  Real rhodotTot = 0.0;
  if (procID() == baseProc)
  {
    CH_assert(dotVec.size() == numProc());

    rhodot.resize(a_dblIn.size());

    int index = 0;
    for (int p = 0; p < dotVec.size(); p++)
    {
      Vector<Real>& v = dotVec[p];
      for (int i = 0; i < v.size(); i++, index++)
      {
        rhodot[index] = v[i];
      }
    }

    rhodot.sort();

    for (int ivec = 0; ivec < rhodot.size(); ivec++)
    {
      rhodotTot += rhodot[ivec];
    }
  }

  //broadcast the sum to all processors.
  broadcast(rhodotTot, baseProc);

  //return the total
  return rhodotTot;
}
#include "NamespaceFooter.H"
