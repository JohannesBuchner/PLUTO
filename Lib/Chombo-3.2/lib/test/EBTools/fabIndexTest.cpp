#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <string>
#include "DataIterator.H"
#include "AllRegularService.H"
#include "FaceIterator.H"
#include "SPMD.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "ParmParse.H"
#include "parstream.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "SlabService.H"
#include "BoxIterator.H"
#include "BaseIVFAB.H"
#include "BaseIFFAB.H"
#include "BaseEBCellFAB.H"
#include "BaseEBFaceFAB.H"
#include "VoFIterator.H"
#include "DebugDump.H"
#include "UsingNamespace.H"

/******************/
/// define grids by splitting up domain
/******************/
int
makeLayout(DisjointBoxLayout& dbl1,
           const Box& domain);

/***************/
// define a slab EBIS.
/***************/
int makeGeometry(Box& domain,
                 Real& dx);

/***************/
//make the corresponding layout
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& nghost);
/***************/
/***************/
int checkFabIndex(const EBISLayout& a_ebisl,
                  const DisjointBoxLayout& a_grids,
                  const Box& a_domain,
                  const int& nghost);
/***************/
/***************/
int
getFabVal(const IntVect& a_iv, const int& comp);

/***************/
/***************/

int
main(int argc, char** argv)
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //begin scoping trick
  {
    const char* in_file = "fabindex.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    int eekflag = 0;
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    Box domain;
    Real dx;
    //make the gometry.  this makes the first Chombo_EBIS
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domain,  dx);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    //make grids
    DisjointBoxLayout grids;
    eekflag = makeLayout(grids, domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeLayouts");
      }

    ///create ebislayout
    int nghost = 2;
    EBISLayout ebisl;
    eekflag = makeEBISL(ebisl, grids, domain, nghost);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeEBISL");
      }

    eekflag = checkFabIndex(ebisl, grids, domain, nghost);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in checkFabIndex");
      }

    pout() << "fab index test passed" << endl;
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return 0;
}
/***************/
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& a_nghost)
{
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
 CH_assert(ebisPtr->isDefined());
  ebisPtr->fillEBISLayout(a_ebisl, a_grids, a_domain, a_nghost);
  return 0;
}
/***************/
/***************/
int
makeLayout(DisjointBoxLayout& a_dbl,
           const Box& a_domain)
{
  ParmParse pp;
  int eekflag= 0;
  int ipieces;
  ipieces = Max(ipieces, 1);
  int maxsize;
  pp.get("maxboxsize",maxsize);
  Vector<Box> vbox(1, a_domain);
  domainSplit(a_domain, vbox,  maxsize);
  if (eekflag != 0)
    {
      pout() << "problem in domainsplit" << endl;
      return eekflag;
    }
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vbox);
  if (eekflag != 0)
    {
      pout() << "problem in loadbalance" << endl;
      return eekflag;
    }
  a_dbl.define(vbox, procAssign);
  return eekflag;
}
/**********/
/**********/
int makeGeometry(Box& a_domain, Real& a_dx)
{
  //parse input file
  ParmParse pp;

  vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);

 CH_assert(n_cell.size() == SpaceDim);
  IntVect lo = IntVect::Zero;
  IntVect hi;
  for (int ivec = 0; ivec < SpaceDim; ivec++)
    {
      if (n_cell[ivec] <= 0)
        {
          pout() << " bogus number of cells input = " << n_cell[ivec];
          return(-1);
        }
      hi[ivec] = n_cell[ivec] - 1;
    }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  vector<Real> prob_lo(SpaceDim, 0.0);
  vector<Real> prob_hi(SpaceDim, 1.0);
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.getarr("prob_hi",prob_hi,0,SpaceDim);
  a_dx = (prob_hi[0]-prob_lo[0])/n_cell[0];

  int whichgeom;
  pp.get("whichgeom", whichgeom);
  if (whichgeom == 1)
    {
      vector<int> slab_lo(SpaceDim);
      pp.getarr("slab_lo",slab_lo,0,SpaceDim);
      vector<int> slab_hi(SpaceDim);
      pp.getarr("slab_hi",slab_hi,0,SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          lo[idir] = slab_lo[idir];
          hi[idir] = slab_hi[idir];
        }
      Box coveredBox(lo,hi);
      SlabService slab(coveredBox);
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      RealVect origin = RealVect::Zero;
      ebisPtr->define(a_domain, origin, a_dx, slab);
    }
  else if (whichgeom == 0)
    {
      AllRegularService service;
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      RealVect origin = RealVect::Zero;
      ebisPtr->define(a_domain, origin, a_dx, service);
    }

  return 0;
}
#if 0
int checkFabIndex(const EBISLayout& a_ebisl,
                  const DisjointBoxLayout& a_grids,
                  const Box& a_domain,
                  const int& a_nghost)
{
  int nvar = 1;
  for (DataIterator dit=a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = a_grids.get(dit());
      IntVectSet localIVS(localBox);
      const EBISBox& ebisBox = a_ebisl[dit()];

      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          BaseEBFaceFAB<int> baseebfacefab(ebisBox, localBox, faceDir, nvar);
        }
    }
  return 0;
}
#endif
#if 1
/***************/
/***************/
int checkFabIndex(const EBISLayout& a_ebisl,
                  const DisjointBoxLayout& a_grids,
                  const Box& a_domain,
                  const int& a_nghost)
{
  int eekflag = 0;
  int nvar = SpaceDim;
  for (DataIterator dit=a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = a_grids.get(dit());
      IntVectSet localIVS(localBox);
      const EBISBox& ebisBox = a_ebisl[dit()];

      BaseFab<int> basefab(grow(localBox,1), nvar);
      for (BoxIterator bit(grow(localBox, 1)); bit.ok(); ++bit)
        {
          for (int ivar = 0; ivar < nvar; ivar++)
            {
              basefab(bit(), ivar) = getFabVal(bit(), ivar);
            }
        }

      BaseIVFAB<int>        baseivfab(localIVS, ebisBox.getEBGraph(), nvar);
      BaseEBCellFAB<int> baseebcellfab(ebisBox, localBox, nvar);
      baseivfab.setVal(-100);
      baseebcellfab.setVal(-100);
      VoFIterator vofit(localIVS, ebisBox.getEBGraph());
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const IntVect& iv = vofit().gridIndex();
          for (int ivar = 0; ivar < nvar; ivar++)
            {
              int rightans = getFabVal(iv, ivar);
              baseivfab(vofit(), ivar) = rightans;
              baseebcellfab(vofit(), ivar) = rightans;
              basefab(iv, ivar) = rightans;
            }
        }
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const IntVect& iv = vofit().gridIndex();
          for (int ivar = 0; ivar < nvar; ivar++)
            {
              int rightans = getFabVal(iv, ivar);
              int ivfabans = baseivfab(vofit(), ivar);
              int basefabans = basefab(iv, ivar);
              int baseebcellfabans = baseebcellfab(vofit(), ivar);
              if (basefabans != rightans)
                {
                  eekflag = 1;
                  return eekflag;
                }
              if (ivfabans != rightans)
                {
                  eekflag = 2;
                  return eekflag;
                }
              if (baseebcellfabans != rightans)
                {
                  eekflag = 3;
                  return eekflag;
                }
            }
        }

      // now check setvals
      for (int ivar=0; ivar<nvar; ivar++)
        {
          basefab.setVal(100+ivar, ivar);
          baseivfab.setVal(ivar, 100+ivar);
          baseebcellfab.setVal(ivar, 100+ivar);
        }

      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const IntVect& iv = vofit().gridIndex();
          for (int ivar = 0; ivar < nvar; ivar++)
            {
              int rightans = 100+ivar;
              int ivfabans = baseivfab(vofit(), ivar);
              int basefabans = basefab(iv, ivar);
              int baseebcellfabans = baseebcellfab(vofit(), ivar);
              if (basefabans != rightans)
                {
                  eekflag = 11;
                  return eekflag;
                }
              if (ivfabans != rightans)
                {
                  eekflag = 12;
                  return eekflag;
                }
              if (baseebcellfabans != rightans)
                {
                  eekflag = 13;
                  return eekflag;
                }
            }
        }


    }
  //now do the face fabs.
  for (DataIterator dit=a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = a_grids.get(dit());
      IntVectSet localIVS(localBox);
      const EBISBox& ebisBox = a_ebisl[dit()];

      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          BaseEBFaceFAB<int> baseebfacefab(ebisBox, localBox, faceDir, nvar);
          BaseIFFAB<int> baseiffab(localIVS, ebisBox.getEBGraph(), faceDir, nvar);
          FaceIterator faceit(localIVS, ebisBox.getEBGraph(), faceDir,
                              FaceStop::SurroundingWithBoundary);
          for (faceit.reset(); faceit.ok(); ++faceit)
            {
              const FaceIndex& face = faceit();
              IntVect iv = face.gridIndex(Side::Lo);
              for (int ivar = 0; ivar < nvar; ivar++)
                {
                  int rightans = getFabVal(iv, ivar);
                  baseebfacefab(face, ivar) = rightans;
                  baseiffab(face, ivar) = rightans;
                }

            }

          for (faceit.reset(); faceit.ok(); ++faceit)
            {
              const FaceIndex& face = faceit();
              IntVect iv = face.gridIndex(Side::Lo);
              for (int ivar = 0; ivar < nvar; ivar++)
                {
                  int rightans = getFabVal(iv, ivar);
                  int baseebfacefabans = baseebfacefab(face, ivar);
                  int baseiffabans = baseiffab(face, ivar);
                  baseebfacefab(face, ivar) = rightans;
                  /**/
                  if (baseebfacefabans != rightans)
                    {
                      eekflag = 4;
                      return eekflag;
                    }
                  if (baseiffabans != rightans)
                    {
                      eekflag = 5;
                      return eekflag;
                    }
                  /**/
                }
            }
        }
    }
  return eekflag;
}
#endif

/***************/
/***************/
int
getFabVal(const IntVect& a_iv, const int& comp)
{
  int rval = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      rval += a_iv[idir]*a_iv[idir];
    }
  rval *= (comp+1);
  return rval;
}
/***************/
/***************/
