#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>

#include "DebugOut.H"
#include "LayoutIterator.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

static const char* indent2 = "      " ;

void dumpNodeLDFPar(const LevelData<NodeFArrayBox>* a_ldfabPtr)
{
  const LevelData<NodeFArrayBox>& ldfab = *a_ldfabPtr;
  Vector<Box> boxes;

  const DisjointBoxLayout& memDBL = ldfab.getBoxes();
  LayoutIterator lit  = memDBL.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
  {
    boxes.push_back(memDBL.get(lit()));
  }

  Vector<int> assign(boxes.size(), 0);

  DisjointBoxLayout locDBL(boxes, assign);
  locDBL.close();

  LevelData<NodeFArrayBox> locLDF(locDBL,
                                  ldfab.nComp(),
                                  ldfab.ghostVect());

  const Interval& interv = ldfab.interval();
  ldfab.copyTo(interv,
               locLDF,
               interv);

  DataIterator dit = locLDF.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDF[dit()].getFab();
    Box nodeGrid = surroundingNodes(locDBL.get(dit()));
    BoxIterator bit(nodeGrid);
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit();
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
  }
}

void dumpNodeLDFLoc(const LevelData<NodeFArrayBox>* a_ldfabPtr)
{
  const LevelData<NodeFArrayBox>& locLDF = *a_ldfabPtr;

  DataIterator dit = locLDF.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDF[dit()].getFab();

    DisjointBoxLayout dbl =  locLDF.disjointBoxLayout();
    Box nodeGrid = surroundingNodes(dbl.get(dit()));
    BoxIterator bit(nodeGrid);
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
  }
}

void dumpNodeFAB(const NodeFArrayBox* a_fabPtr)
{
  const FArrayBox& fab = a_fabPtr->getFab();
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
}
void dumpFAB(const FArrayBox* a_fabPtr)
{
  const FArrayBox& fab = *a_fabPtr;
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
}

void dumpBFR(const BaseFab<Real>* a_fabPtr)
{
  const BaseFab<Real>& fab = *a_fabPtr;
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
}
void dumpBFI(const BaseFab<int>* a_fabPtr)
{
  const BaseFab<int>& fab = *a_fabPtr;
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
}

void dumpFAB2DSlicePretty(const FArrayBox *const a_fabPtr,
                          const int              a_comp,
                          IntVect                a_ivSml,
                          IntVect                a_ivBig,
                          const int              a_prec,
                          std::ostream&          a_out)
{
  const int oldPrec = a_out.precision();
  const int width = a_prec+7;
  a_out.setf(std::ios_base::scientific, std::ios_base::floatfield);
  a_out.precision(a_prec);
  if (a_ivSml == IntVect::Zero && a_ivBig == IntVect::Zero)
    {
      a_ivSml = a_fabPtr->smallEnd();
      a_ivBig = a_fabPtr->bigEnd();
    }
  int sliceDir[2];
  int numSliceDir = 0;
  for (int i = 0; i != SpaceDim; ++i)
    {
      if (a_ivSml[i] != a_ivBig[i])
        {
          sliceDir[std::min(1, numSliceDir++)] = i;
        }
    }
  if (SpaceDim == 2 && numSliceDir == 1)  // Allows 1-D slice in 2D
    {
      sliceDir[1] = (sliceDir[0] == 0) ? 1 : 0;
      numSliceDir = 2;
    }
  CH_assert(numSliceDir == 2);
  const int s0 = sliceDir[0];
  const int s1 = sliceDir[1];
  char lbl[64];
  {  // Use C to find where to put big IntVect printout
    int k = 1;
    for (int i = 0; i != SpaceDim; ++i)
      {
        k += std::sprintf(lbl, "%d,", a_ivBig[i]);
      }
    for (int n = (width+1)*(a_ivBig[s0] - a_ivSml[s0] + 1) - 1 - k; n--;)
      a_out << ' ';
  }
  a_out << a_ivBig << std::endl;
  IntVect iv = a_ivSml;
  for (int j = a_ivBig[s1]; j >= a_ivSml[s1]; --j)
    {
      iv[s1] = j;
      iv[s0] = a_ivSml[s0];
      a_out << std::setw(width) << (*a_fabPtr)(iv, a_comp);
      for (int i = a_ivSml[s0] + 1; i <= a_ivBig[s0]; ++i)
        {
          iv[s0] = i;
          a_out << ' ' << std::setw(width) << (*a_fabPtr)(iv, a_comp);
        }
      a_out << std::endl;
    }
  a_out << a_ivSml << std::endl;
  a_out.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
  a_out.precision(oldPrec);
}

void dumpMaxMin(const LevelData<FArrayBox>* a_ldfabPtr)
{
  for (int icomp = 0; icomp < a_ldfabPtr->nComp(); icomp++)
    {
      Real glomax = -1.0e20;
      Real glomin = 1.0e20;
      for (DataIterator dit = a_ldfabPtr->dataIterator(); dit.ok(); ++dit)
        {
          Real locmax = (*a_ldfabPtr)[dit()].max(icomp);
          Real locmin = (*a_ldfabPtr)[dit()].min(icomp);
          glomax = Max(glomax, locmax);
          glomin = Min(glomin, locmin);
        }
      pout() << setprecision(4) 
             << setiosflags(ios::showpoint) 
             << setiosflags(ios::scientific)
             << "for component " << icomp 
             << ", max = " << glomax 
             << ", min = " << glomin << endl;
    }
}

void dumpLDFPar(const LevelData<FArrayBox>* a_ldfabPtr)
{
  const LevelData<FArrayBox>& ldfab = *a_ldfabPtr;
  Vector<Box> boxes;

  const DisjointBoxLayout& memDBL = ldfab.getBoxes();
  LayoutIterator lit  = memDBL.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
  {
    boxes.push_back(memDBL.get(lit()));
  }

  Vector<int> assign(boxes.size(), 0);

  DisjointBoxLayout locDBL(boxes, assign);
  locDBL.close();

  LevelData<FArrayBox> locLDF(locDBL,
                              ldfab.nComp(),
                              ldfab.ghostVect());

  const Interval& interv = ldfab.interval();
  ldfab.copyTo(interv,
               locLDF,
               interv);

  DataIterator dit = locLDF.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDF[dit()];
    BoxIterator bit(locDBL.get(dit()));
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
  }
}

void dumpLDFLoc(const LevelData<FArrayBox>* a_ldfabPtr)
{
  const LevelData<FArrayBox>& locLDF = *a_ldfabPtr;

  DataIterator dit = locLDF.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDF[dit()];
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          pout() << "\t" << fab(bit(),ivar);
        }
      pout() << endl;
    }
  }
}

void dumpIVSFAB(const IVSFAB<Real>* a_ivsfabPtr)
{
  const IVSFAB<Real>& ivsfab = *a_ivsfabPtr;
  const IntVectSet& ivs = ivsfab.getIVS();
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
  {
    const IntVect& iv = ivsit();

    for (int icomp = 0; icomp < ivsfab.nComp(); icomp++)
    {
      pout() << ivsfab(iv, icomp) << "  " << endl;
    }
  }
}

void dumpDBL(const DisjointBoxLayout* a_dblInPtr)
{
  pout() << indent2 <<"DisjointBoxLayout ";

  if (a_dblInPtr->isClosed())
  {
    pout() << "closed\n";
  }
  else
  {
    pout() << "open\n";
  }

  pout() << *a_dblInPtr << endl;
}

void dumpBL(const BoxLayout* a_dblInPtr)
{
  pout() << indent2 <<"BoxLayout ";

  if (a_dblInPtr->isClosed())
  {
    pout() << "closed\n";
  }
  else
  {
    pout() << "open\n";
  }

  pout() << *a_dblInPtr << endl;
}

void dumpIVS(const IntVectSet* a_ivsPtr)
{
  const IntVectSet& ivs = *a_ivsPtr;
  IVSIterator it(ivs);

  pout() << indent2 << ": IntVects in the IVS are:" << endl;
  pout() << indent2;

  for (it.begin(); it.ok(); ++it)
  {
    pout() << it() << "  ";
  }

  pout() << endl;
}

void dumpBox(const Box* a_boxPtr)
{
  pout() << indent2 << *a_boxPtr << endl;
}

void dumpVBox(const Vector<Box>* a_vectPtr)
{
  Vector<Box> vect = *a_vectPtr;
  for (int ivec = 0; ivec < vect.size(); ivec++)
  {
    pout() << indent2 << vect[ivec] << endl;
  }
}

void dumpVVBox(const Vector<Vector<Box> >* a_vectPtr)
{
  Vector<Vector<Box> > vect = *a_vectPtr;
  for (int iveco = 0; iveco < vect.size(); iveco++)
  {
    pout() << indent2;

    for (int iveci = 0; iveci < vect[iveco].size(); iveci++)
    {
      pout() <<  vect[iveco][iveci] << "   ";
    }

    pout() << endl;
  }
}

void dumpLDDBL(const LevelData<FArrayBox>* a_ldfabPtr)
{
  const DisjointBoxLayout& dbl = a_ldfabPtr->disjointBoxLayout();
  dumpDBL(&dbl);
}
void dumpNodeLDDBL(const LevelData<NodeFArrayBox>* a_ldfabPtr)
{
  const DisjointBoxLayout& dbl = a_ldfabPtr->disjointBoxLayout();
  dumpDBL(&dbl);
}
#include "NamespaceFooter.H"
