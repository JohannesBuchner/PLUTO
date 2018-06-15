#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves Tues Nov 20 2001

#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "SlabService.H"
#include "RedistStencil.H"
#include "EBLevelRedist.H"
#include "EBCoarToFineRedist.H"
#include "EBFineToCoarRedist.H"
#include "EBCoarToCoarRedist.H"
#include "BRMeshRefine.H"
#include "DebugDump.H"

#include "UsingNamespace.H"

/******************/
/******************/
int
makeLayouts(Vector<DisjointBoxLayout>& vecDbl,
            const Box& domainCoarsest,
            const Vector<int>& vecRefRat);
/***************/
/***************/
int makeGeometry(Box& domainCoarsest,
                 Vector<int>& vecRefRat);
/***************/
/***************/
int testConservation(const Vector<DisjointBoxLayout>& a_grids,
                     const Box& a_coarsestDomain,
                     const Vector<int>& vecRefRat);

/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/
int
main(int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int eekflag = 0;
  //begin forever present scoping trick
  {

    const char* in_file = "fullredist.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    Box domainCoarsest, domainFinest;
    Vector<DisjointBoxLayout> vecGrids;
    Vector<int> vecRefRat;
    //make the gometry.  this makes the first Chombo_EBIS
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domainFinest,  vecRefRat);
    if (eekflag != 0)
      {
        pout() << "fullredisttest: problem in makegeom" << endl;
        return eekflag;
      }
    int numlevels = vecRefRat.size();
    CH_assert( numlevels == 3 );
    domainCoarsest = domainFinest;
    for (int ilev = 0; ilev < numlevels-1; ilev++)
      {
        domainCoarsest.coarsen(vecRefRat[ilev]);
      }
    eekflag = makeLayouts(vecGrids, domainCoarsest,vecRefRat);
    if (eekflag != 0)
      {
        pout() << "fullredisttest: problem in makelayouts" << endl;
        return eekflag;
      }

    eekflag = testConservation(vecGrids, domainCoarsest, vecRefRat);
    if (eekflag != 0)
      {
        pout() << "fullredisttest: problem in testConservation" << endl;
        return eekflag;
      }
    if (eekflag != 0) return eekflag;
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  pout() << "full redistribution test passed" << endl;
  return eekflag;
}
/***************/
/***************/
Real massFunc(const VolIndex& a_vof, const Box& a_domain)
{
  const IntVect& iv = a_vof.gridIndex();
  Real retval;
  Real problen = a_domain.size(0);
  Real x = iv[0];
  Real y = iv[1];
  x /= problen;
  y /= problen;

  retval = x*x*x + y*y*y;

//  IntVect ivdebug(6,7);
//  if (iv == ivdebug)
//    retval = 4.5;
//  else
//    retval = 0.0;
  return retval;
}
/***************/
/***************/
int testConservation(const Vector<DisjointBoxLayout>& a_vecGrids,
                     const Box& a_coarsestDomain,
                     const Vector<int>& a_vecRefRat)
{
  int eekflag = 0;
  CH_assert(a_vecRefRat.size() == 3);
  CH_assert(a_vecGrids.size() == 3);
  //We start with a known amount of mass to redistribute.
  //on the middle level of refinement.  We do the level redistribution
  //the the coarse-fine redistribution, then the fine-coarse
  //redistribution and finally  the coarse-coarse redistribution.
  //Finally we add up the mass on all valid regions and see
  //if it is the same as the original sum.

  Box fineDomain, mediDomain, coarDomain;
  coarDomain = a_coarsestDomain;
  mediDomain = refine(coarDomain, a_vecRefRat[0]);
  fineDomain = refine(mediDomain, a_vecRefRat[1]);
  const DisjointBoxLayout& coarGrids =a_vecGrids[0];
  const DisjointBoxLayout& mediGrids =a_vecGrids[1];
  const DisjointBoxLayout& fineGrids =a_vecGrids[2];
  EBISLayout coarEBISL,mediEBISL,fineEBISL;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  int redistRad;
  ParmParse pp;
  pp.get("redist_radius",redistRad);
  int nghost = 3*redistRad;
  ebisPtr->fillEBISLayout(coarEBISL, coarGrids, coarDomain, nghost);
  ebisPtr->fillEBISLayout(mediEBISL, mediGrids, mediDomain, nghost);
  ebisPtr->fillEBISLayout(fineEBISL, fineGrids, fineDomain, nghost);
  int iverbose;
  pp.get("verbose", iverbose);
  if (iverbose == 1)
    {
      for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& fineBox = fineGrids.get(dit());
          const EBISBox& ebisBox = fineEBISL[dit()];
          IntVectSet ivs = ebisBox.getIrregIVS(fineBox);
          dumpIVS(&ivs);
        }
      for (DataIterator dit = mediGrids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& mediBox = mediGrids.get(dit());
          const EBISBox& ebisBox = mediEBISL[dit()];
          IntVectSet ivs = ebisBox.getIrregIVS(mediBox);
          dumpIVS(&ivs);
        }
      for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
        {
          const Box& coarBox = coarGrids.get(dit());
          const EBISBox& ebisBox = coarEBISL[dit()];
          IntVectSet ivs = ebisBox.getIrregIVS(coarBox);
          dumpIVS(&ivs);
        }
    }

  EBCellFactory fineFact(fineEBISL);
  EBCellFactory mediFact(mediEBISL);
  EBCellFactory coarFact(coarEBISL);
  IntVect ivghost = IntVect::Zero;
  int nvar = 1;
  LevelData<EBCellFAB> fineSoln(fineGrids, nvar, ivghost, fineFact);
  LevelData<EBCellFAB> mediSoln(mediGrids, nvar, ivghost, mediFact);
  LevelData<EBCellFAB> coarSoln(coarGrids, nvar, ivghost, coarFact);
  //set initial solution to zero
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    fineSoln[dit()].setVal(0.0);
  for (DataIterator dit = mediGrids.dataIterator(); dit.ok(); ++dit)
    mediSoln[dit()].setVal(0.0);
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
    coarSoln[dit()].setVal(0.0);

  IntVectSet fineInvalidIVS;
  IntVectSet mediInvalidIVS;
  IntVectSet coarInvalidIVS;
  //medi grids are valid where there is no fine grid
  for (LayoutIterator lit=fineGrids.layoutIterator();lit.ok();++lit)
    {
      Box cedBox = coarsen(fineGrids.get(lit()), a_vecRefRat[1]);
      mediInvalidIVS |= cedBox;
    }
  //coar grids are valid where there is no medi grid
  for (LayoutIterator lit=mediGrids.layoutIterator();lit.ok();++lit)
    {
      Box cedBox = coarsen(mediGrids.get(lit()), a_vecRefRat[0]);
      coarInvalidIVS |= cedBox;
    }
  //all mass starts from the medium grids.
  Real localSumOld = 0.0;
  LayoutData<BaseIVFAB<Real> > mass(mediGrids);
  for (DataIterator dit = mediGrids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& mediBox = mediGrids.get(dit());
      const EBISBox& ebisBox = mediEBISL[dit()];
      IntVectSet irregIVS = ebisBox.getIrregIVS(mediBox);
      BaseIVFAB<Real>& massFAB= mass[dit()];
      massFAB.define(irregIVS, ebisBox.getEBGraph(), 1);

      for (VoFIterator vofit(irregIVS, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          Real massVoF = massFunc(vofit(), mediDomain);
          massFAB(vofit(), 0) = massVoF;
          if (!mediInvalidIVS.contains(vofit().gridIndex()))
            {
              localSumOld += massVoF;
            }
        }
    }
  //do all the redistribution fandango.
  //Remember Alice? This is a song about Alice.
  RedistStencil mediRDSten(mediGrids, mediEBISL,
                           mediDomain, redistRad);

  pp.get("verbose", iverbose);
  if (iverbose == 1)
    {
      for (DataIterator dit = mediGrids.dataIterator(); dit.ok(); ++dit)
        {
          const BaseIVFAB<VoFStencil>& stenFAB = mediRDSten[dit()];
          VoFIterator vofit(stenFAB.getIVS(), mediEBISL[dit()].getEBGraph());
          for (vofit.reset(); vofit.ok(); ++vofit)
            {
              const VoFStencil& sten = stenFAB(vofit(), 0);

              ////debug
              //if (vofit().gridIndex() == IntVect(6,7))
              {
                pout() << "srcvof=" << vofit() << endl;
                for (int isten = 0; isten < sten.size(); isten++)
                  {
                    const VolIndex& vof  = sten.vof(isten);
                    const Real& weight   = sten.weight(isten);
                    pout() << isten
                           << "  dstvof=" << vof
                           << ", weight=" << weight << endl;
                  }
                pout() << endl;
              }
            }
        }
    }
  EBLevelRedist mediRedist;
  EBFineToCoarRedist fineToCoar;
  EBCoarToFineRedist coarToFine;
  EBCoarToCoarRedist coarToCoar;

  mediRedist.define(mediGrids, mediEBISL,
                    mediDomain,nvar, redistRad);

  fineToCoar.define(mediGrids, coarGrids,
                    mediEBISL, coarEBISL,
                    coarDomain, a_vecRefRat[0],
                    nvar, redistRad);
  coarToFine.define(fineGrids, mediGrids,
                    fineEBISL, mediEBISL,
                    mediDomain, a_vecRefRat[1],
                    nvar, redistRad, Chombo_EBIS::instance());
  coarToCoar.define(fineGrids, mediGrids,
                    fineEBISL, mediEBISL,
                    mediDomain, a_vecRefRat[1],
                    nvar, redistRad);

  //set all the buffers to zero
  mediRedist.setToZero();
  fineToCoar.setToZero();
  coarToFine.setToZero();
  coarToCoar.setToZero();

  //store the mass in each buffer.
  Interval interv(0,0);
  for (DataIterator dit = mediGrids.dataIterator(); dit.ok(); ++dit)
    {
      const BaseIVFAB<Real>& massFAB = mass[dit()];
      mediRedist.increment(massFAB, dit(), interv);
      fineToCoar.increment(massFAB, dit(), interv);
      coarToFine.increment(massFAB, dit(), interv);
      coarToCoar.increment(massFAB, dit(), interv);
    }

  //redistribute the mass on this (the medium) level
  mediRedist.redistribute(mediSoln, interv);
  //redistribute the mass to the next coarser level that went
  //off the level
  fineToCoar.redistribute(coarSoln, interv);
  //redistribute the mass to the next finer level that went
  //under the finer level
  coarToFine.redistribute(fineSoln, interv);
  //unredistribute mass that came from invalid regions
  coarToCoar.redistribute(mediSoln, interv);

  //add up mass in VALID areas of each solution
  //all fine grids are valid
  Real localSumNew = 0.0;
  Real localSumFine=0.0;
  Real localSumMedi=0.0;
  Real localSumCoar=0.0;

  Real fineCellVol = 1;
  Real mediCellVol = 1;
  Real coarCellVol = 1;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      coarCellVol *= a_vecRefRat[0];
      fineCellVol /= a_vecRefRat[1];
    }
  for (DataIterator dit = fineGrids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = fineGrids.get(dit());
      const EBISBox& ebisBox = fineEBISL[dit()];
      IntVectSet ivs(localBox);
      const EBCellFAB& soln = fineSoln[dit()];
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex vof = vofit();
          //since invalid ivs is empty, this should always pass
          if (!fineInvalidIVS.contains(vof.gridIndex()))
            {
              Real density = soln(vof, 0);
              ////debug
              //if (density > 1.0e-8)
              {
                Real volfrac = ebisBox.volFrac(vof);
                localSumNew += density*volfrac*fineCellVol;
                localSumFine +=density*volfrac*fineCellVol;
              }
            }
        }
    }
  for (DataIterator dit = mediGrids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = mediGrids.get(dit());
      const EBISBox& ebisBox = mediEBISL[dit()];
      IntVectSet ivs(localBox);
      const EBCellFAB& soln = mediSoln[dit()];
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex vof = vofit();
          if (!mediInvalidIVS.contains(vof.gridIndex()))
            {
              Real density = soln(vof, 0);
              Real volfrac = ebisBox.volFrac(vof);
              localSumNew +=  density*volfrac*mediCellVol;
              localSumMedi += density*volfrac*mediCellVol;
            }
        }
    }
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& localBox = coarGrids.get(dit());
      const EBISBox& ebisBox = coarEBISL[dit()];
      IntVectSet ivs(localBox);
      const EBCellFAB& soln = coarSoln[dit()];
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex vof = vofit();
          if (!coarInvalidIVS.contains(vof.gridIndex()))
            {
              Real density = soln(vof, 0);
              Real volfrac = ebisBox.volFrac(vof);
              localSumNew +=  density*volfrac*coarCellVol;
              localSumCoar += density*volfrac*coarCellVol;
            }
        }
    }

  //gather what each processor thinks is the sum old and new
  Vector<Real> sumVecOld, sumVecNew;
  int baseproc = 0;
  gather(sumVecOld, localSumOld, baseproc);
  gather(sumVecNew, localSumNew, baseproc);
  Real massTotOld = 0;
  Real massTotNew = 0;
  if (procID() == baseproc)
    {
      CH_assert(sumVecOld.size() == numProc());
      CH_assert(sumVecNew.size() == numProc());
      for (int ivec = 0; ivec < numProc(); ivec++)
        {
          massTotOld += sumVecOld[ivec];
          massTotNew += sumVecNew[ivec];
        }
    }
  //broadcast the sum to all processors.
  broadcast(massTotOld, baseproc);
  broadcast(massTotNew, baseproc);

  if (Abs(massTotOld - massTotNew) > 1.0e-6)
    eekflag = 42;
  return eekflag;
}
/***************/
/***************/
int
makeTags(IntVectSet& a_tags,
         const Box& a_domainCoar)
{
#if (CH_SPACEDIM ==2)
  int nc = a_domainCoar.size(0);
  if (nc < 16)
    {
      pout() << "coarsest level must be at least 16x16" << endl;
      return -22;
    }
  int nmi = nc/2;//16
  int nqu = nc/4;//8
  int ntf = (nc*3)/4;  //24
  int nte = (nc*3)/8; //12
  int nfe = (nc*5)/8; //20
  Box boxf2(IntVect(nmi,nte), IntVect(ntf-1,nfe-1));
  Box boxf4(IntVect(nfe,nqu), IntVect(nc -1,nte-1));
#else
  int nc = a_domainCoar.size(0);
  int nmi = nc/2;//16
  int nqu = nc/4;//8
  int ntf = (nc*3)/4;  //24
  Box boxf2(IntVect(0, nqu,nqu), IntVect(nmi-1,ntf-1,ntf-1));
  Box boxf4(IntVect(0, 0, 0), IntVect(nqu-1,nqu-1,nqu-1));
#endif
  a_tags |= boxf2;
  a_tags |= boxf4;
  return 0;
}
int
makeLayouts(Vector<DisjointBoxLayout>& a_vecDBL,
            const Box& a_domainCoarsest,
            const Vector<int>& a_vecRefRat)
{
  Vector<Vector<Box> > newMeshes;
  ParmParse pp;
  int eekflag= 0;
  int topLevel = a_vecRefRat.size() - 2;
  CH_assert(topLevel >= 0);
  int nlevels = topLevel + 2;
  int baseLevel = 0;
  int isimple;
  pp.get("simple_layout", isimple);
  bool simpleLayout = (isimple==1);
  if (!simpleLayout)
    {
      int blockFactor, bufferSize, maxSize;
      Real fillRat;
      pp.get("block_factor", blockFactor);
      pp.get("buffer_size", bufferSize);
      pp.get("maxboxsize", maxSize);
      pp.get("fill_ratio", fillRat);

      BRMeshRefine mesher(a_domainCoarsest, a_vecRefRat,
                          fillRat, blockFactor, bufferSize, maxSize);

      //tags at base level
      IntVectSet tags;
      eekflag = makeTags(tags, a_domainCoarsest);
      if (eekflag < 0) return eekflag;

      Vector<Vector<Box> > oldMeshes(nlevels);
      oldMeshes[0] = Vector<Box>(1, a_domainCoarsest);
      Box finerDomain = a_domainCoarsest;
      for (int ilev = 1; ilev < nlevels; ilev++)
        {
          finerDomain.refine(a_vecRefRat[ilev]);
          oldMeshes[ilev] = Vector<Box>(1, finerDomain);
        }
      mesher.regrid(newMeshes, tags, baseLevel, topLevel, oldMeshes);
    }
  else
    {
      newMeshes.resize(3);
      IntVect lo,hi;
      lo = IntVect::Zero;
      hi = 7*IntVect::Unit;
      Box coarBox(lo,hi);
      lo = 4*IntVect::Unit;
      hi =11*IntVect::Unit;
      Box mediBox(lo,hi);
      lo =12*IntVect::Unit;
      hi =19*IntVect::Unit;
      Box fineBox(lo,hi);
      newMeshes[0] = Vector<Box>(1, coarBox);
      newMeshes[1] = Vector<Box>(1, mediBox);
      newMeshes[2] = Vector<Box>(1, fineBox);
    }

  a_vecDBL.resize(nlevels);
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      Vector<int> procAssign;
      eekflag = LoadBalance(procAssign, newMeshes[ilev]);
      if (eekflag != 0) return eekflag;
      a_vecDBL[ilev].define(newMeshes[ilev], procAssign);
    }
  int iverbose;
  pp.get("verbose", iverbose);
  if (iverbose == 1)
    {
      for (int ilev = 0; ilev < nlevels; ilev++)
        {
          pout() << "the grids at level " << ilev
                 << " = " << a_vecDBL[ilev] << endl;
        }

    }
  return eekflag;
}
/**********/
int makeGeometry(Box& a_domainFinest,
                 Vector<int>& a_vecRefRat)
{
  int eekflag =  0;
  //parse input file
  ParmParse pp;
  RealVect origin = RealVect::Zero;

#if (CH_SPACEDIM==2)
  Vector<int> n_cell(SpaceDim, 64);
#else
  Vector<int> n_cell(SpaceDim, 16);
#endif

  int numLevels = 3;
  pp.getarr("ref_ratio",a_vecRefRat,0, numLevels);

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

  a_domainFinest.setSmall(lo);
  a_domainFinest.setBig(hi);

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  Real dx = (prob_hi-prob_lo[0])/n_cell[0];
  if (dx < 1.0e-10)
    {
      pout() << " probhi lower than problo " << endl;
      return(-2);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = prob_lo[idir];
    }

  int whichgeom;
  pp.get("which_geom",whichgeom);

  if (whichgeom == 0)
    {
      //allregular
      pout() << "all regular geometry" << endl;
      AllRegularService regserv;
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domainFinest, origin, dx, regserv);
    }
  else if (whichgeom == 1)
    {
      pout() << "ramp geometry" << endl;
      int upDir;
      int indepVar;
      Real startPt;
      Real slope;
      pp.get("up_dir",upDir);
      pp.get("indep_var",indepVar);
      pp.get("start_pt", startPt);
      pp.get("ramp_slope", slope);

      RealVect normal = RealVect::Zero;
      normal[upDir] = 1.0;
      normal[indepVar] = -slope;

      RealVect point = RealVect::Zero;
      point[upDir] = -slope*startPt;

      bool normalInside = true;

      PlaneIF ramp(normal,point,normalInside);

      RealVect vectDx = RealVect::Unit;
      vectDx *= dx;

      GeometryShop workshop(ramp,0,vectDx);
      //this generates the new EBIS
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domainFinest, origin, dx, workshop);
    }
  else if (whichgeom == 2)
    {
      pout() << "slab geometry" << endl;
      vector<int> slab_lo(SpaceDim);
      pp.getarr("slab_lo",slab_lo,0,SpaceDim);
      vector<int> slab_hi(SpaceDim);
      pp.getarr("slab_hi",slab_hi,0,SpaceDim);
      IntVect lo, hi;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          lo[idir] = slab_lo[idir];
          hi[idir] = slab_hi[idir];
        }
      Box coveredBox(lo,hi);
      SlabService slab(coveredBox);
      //this generates the new EBIS
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      RealVect origin = RealVect::Zero;
      ebisPtr->define(a_domainFinest, origin, dx, slab);
    }
  else
    {
      //bogus which_geom
      pout() << " bogus which_geom input = " << whichgeom << endl;
      eekflag = 33;
    }

  return eekflag;
}
/***************/
/***************/
