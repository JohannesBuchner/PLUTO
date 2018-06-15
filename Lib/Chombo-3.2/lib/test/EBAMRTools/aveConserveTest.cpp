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
#include "SlabService.H"
#include "PlaneIF.H"
#include "RedistStencil.H"
#include "EBLevelRedist.H"
#include "EBCoarToFineRedist.H"
#include "EBFineToCoarRedist.H"
#include "EBCoarToCoarRedist.H"
#include "BRMeshRefine.H"
#include "DebugDump.H"
#include "EBLevelDataOps.H"
#include "EBCoarseAverage.H"
#include "EBFABView.H"
#include "EBDebugDump.H"
#include "UsingNamespace.H"

int f_ilev;
IntVect  f_debugIVCoar(D_DECL(26,11,0));
/***************/
Real denseFunc(const VolIndex& a_vof,  const Real& a_dx)
{
  const IntVect& iv = a_vof.gridIndex();
  Real x = a_dx*(iv[0]+ 0.5);
  Real y = a_dx*(iv[1]+ 0.5);

  Real retval = x*x*x + y*y*y;

//  IntVect ivdebug(6,7);
//  if (iv == ivdebug)
//    retval = 4.5;
//  else
//    retval = 0.0;
//  retval = 4.0;
  return retval;
}
/***************/
void setDenseData(Vector<LevelData<EBCellFAB> *>&      a_data,
                 const Vector<DisjointBoxLayout>&      a_grids,
                 const Vector<EBISLayout>&             a_ebisl,
                 const Box&                            a_coarsestDomain,
                 const Vector<int>&                    a_refRat,
                 const Real&                           a_dxCoar,
                 bool                                  a_doingRZ)
{
  Box curDomain = a_coarsestDomain;
  Real curDx = a_dxCoar;
  for (int ilev = 0; ilev < a_data.size(); ilev++)
    {
      for (DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          const Box& box = a_grids[ilev].get(dit());
          IntVectSet ivsIter(box);
          const EBISBox ebisBox = a_ebisl[ilev][dit()];
          for (VoFIterator vofit(ivsIter, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              Real massVal = denseFunc(vofit(), curDx);
              (*a_data[ilev])[dit()](vofit(), 0) = massVal;
            }
        }
      curDomain.refine(a_refRat[ilev]);
      curDx /= a_refRat[ilev];
    }
}


Real sumMassLevel(LevelData<EBCellFAB>&         a_data,
                  const DisjointBoxLayout&      a_grids,
                  const EBISLayout&             a_ebisl,
                  const Box&                    a_domain,
                  const Real&                   a_dx,
                  bool                          a_doingRZ,
                  const IntVectSet&             a_invalidIVS)
{
  Box debugBoxCoar(f_debugIVCoar, f_debugIVCoar);
  Box debugBoxFine = refine(debugBoxCoar, 2);
  Real localSum= 0.0;
  Real cellVolXY = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      cellVolXY *= a_dx;
    }
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& box = a_grids.get(dit());
      IntVectSet ivsIter(box);
      ivsIter -= a_invalidIVS;
      const EBISBox ebisBox = a_ebisl[dit()];
      for (VoFIterator vofit(ivsIter, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
        {
          Real compVol = 1.0;
          Real kappa = ebisBox.volFrac(vofit());
          compVol = kappa*cellVolXY;

          Real dataVal = a_data[dit()](vofit(), 0);
          localSum += dataVal*compVol;
        }
    }

  Vector<Real> sumVec;
  int baseproc = 0;
  gather(sumVec, localSum, baseproc);
  Real massTot = 0;
  if (procID() == baseproc)
    {
      CH_assert(sumVec.size() == numProc());
      for (int ivec = 0; ivec < numProc(); ivec++)
        {
          massTot += sumVec[ivec];
        }
    }
  broadcast(massTot, baseproc);

  return massTot;
}
/***************/
Real sumMass(const Vector<LevelData<EBCellFAB> *>& a_data,
             const Vector<DisjointBoxLayout>&      a_grids,
             const Vector<EBISLayout>&             a_ebisl,
             const Box&                            a_coarsestDomain,
             const Vector<int>&                    a_refRat,
             const Real&                           a_dxCoar,
             bool                                  a_doingRZ)
{
  Real massTot = 0.0;
  Box curDomain = a_coarsestDomain;
  Real curDx = a_dxCoar;
  for (int ilev = 0; ilev < a_data.size(); ilev++)
    {
      f_ilev = ilev;
      IntVectSet invalidIVS;
      if ((ilev + 1) < (a_data.size()))
        {
          //medi grids are valid where there is no fine grid
          for (LayoutIterator lit=a_grids[ilev+1].layoutIterator();lit.ok();++lit)
            {
              Box cedBox = coarsen(a_grids[ilev+1].get(lit()), a_refRat[ilev]);
              invalidIVS |= cedBox;
            }
        }

      Real massLevel = sumMassLevel(*a_data[ilev], a_grids[ilev], a_ebisl[ilev],
                                    curDomain, curDx, a_doingRZ, invalidIVS);
      massTot += massLevel;

      curDomain.refine(a_refRat[ilev]);
      curDx /= a_refRat[ilev];
    }

  return massTot;
}

/***************/
int testConservation(const Vector<DisjointBoxLayout>& a_grids,
                     const Box&                       a_coarsestDomain,
                     const Vector<int>&               a_refRat,
                     const Real&                      a_dxCoar,
                     bool                             a_doingRZ)
{
  int eekflag = 0;
  int nlevels = a_grids.size();
  int nghost = 4;
  int nvar = 1;
  IntVect ivghost = IntVect::Zero;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  Vector<EBISLayout> ebisl(a_grids.size());
  Vector<LevelData<EBCellFAB>* > data(a_grids.size(), NULL);
  Box curDomain = a_coarsestDomain;
  for (int ilev = 0; ilev < nlevels ; ilev++)
    {
      ebisPtr->fillEBISLayout(ebisl[ilev], a_grids[ilev], curDomain, nghost);
      curDomain.refine(a_refRat[ilev]);

      EBCellFactory fact(ebisl[ilev]);
      data[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], nvar, ivghost, fact);
      EBLevelDataOps::setToZero(*data[ilev]);
    }

  //set data to something
  setDenseData(data, a_grids, ebisl, a_coarsestDomain,
               a_refRat, a_dxCoar, a_doingRZ);

  //find the total mass before over all
  //sum mass excluding data covered by finer levels
  Real sumBeforeAverage = sumMass(data, a_grids, ebisl, a_coarsestDomain,
                                  a_refRat, a_dxCoar, a_doingRZ);

  //average down data
  Interval interv(0,0);
  for (int ilev = nlevels-1; ilev > 0; ilev--)
    {
      Box coarDomain = a_coarsestDomain;
      for (int jlev = 0; jlev < ilev-1; jlev++)
        {
          coarDomain.refine(a_refRat[jlev]);
        }
      Real dxFine = a_dxCoar;
      for (int jlev = 0; jlev < ilev; jlev++)
        {
          dxFine /= a_refRat[jlev];
        }
      EBCoarseAverage averageOp(a_grids[ilev  ],
                                a_grids[ilev-1],
                                ebisl[ilev  ],
                                ebisl[ilev-1],
                                coarDomain,
                                a_refRat[ilev-1],
                                nvar,
                                ebisPtr);
      averageOp.average(*data[ilev-1],
                        *data[ilev  ],
                        interv);
    }

  //sum at level 0 (with nothing invalid) should now  be the same as sum over all levels
  //excluding stuff covered by coarser levels
  IntVectSet invalidIVS;
  f_ilev = 0;
  Real sumAfterAverage = sumMassLevel(*data[0], a_grids[0], ebisl[0],
                                      a_coarsestDomain, a_dxCoar, a_doingRZ, invalidIVS);

  if (Abs(sumBeforeAverage - sumAfterAverage) > 1.0e-6)
    eekflag = 42;

  for (int ilev = 0; ilev < nlevels ; ilev++)
    {
      delete data[ilev];
      data[ilev] = NULL;
    }
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
            const Vector<int>& a_refRat)
{
  Vector<Vector<Box> > newMeshes;
  ParmParse pp;
  int eekflag= 0;
  int topLevel = a_refRat.size() - 2;
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

      BRMeshRefine mesher(a_domainCoarsest, a_refRat,
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
          finerDomain.refine(a_refRat[ilev]);
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

  return eekflag;
}
/**********/
int makeGeometry(Box& a_domainFinest,
                 Vector<int>& a_refRat,
                 Real& a_dxCoar)
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

  int numLevels;
  pp.get("max_level", numLevels);
  numLevels++;
  pp.getarr("ref_ratio",a_refRat,0, numLevels);

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
  if (dx < 0)
    {
      pout() << " probhi lower than problo " << endl;
      return(-2);
    }
  a_dxCoar = dx;
  for (int ilev = 0; ilev < numLevels-1; ilev++)
    {
      a_dxCoar *= a_refRat[ilev];
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
int
main(int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int eekflag = 0;
  //begin forever present scoping trick
  {

    const char* in_file = "aveconserve.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    Box domainCoarsest, domainFinest;
    Vector<DisjointBoxLayout> grids;
    Vector<int> refRat;
    //make the gometry.  this makes the first Chombo_EBIS
    //and defines it using a geometryservice
    Real dxCoar;

    eekflag =  makeGeometry(domainFinest,  refRat, dxCoar);
    if (eekflag != 0)
      {
        pout() << "fullredisttest: problem in makegeom" << endl;
        return eekflag;
      }
    int numlevels = refRat.size();
    domainCoarsest = domainFinest;
    for (int ilev = 0; ilev < numlevels-1; ilev++)
      {
        domainCoarsest.coarsen(refRat[ilev]);
      }
    eekflag = makeLayouts(grids, domainCoarsest,refRat);
    if (eekflag != 0)
      {
        pout() << "fullredisttest: problem in makelayouts" << endl;
        return eekflag;
      }

    eekflag = testConservation(grids, domainCoarsest, refRat, dxCoar, false);
    if (eekflag != 0)
      {
        pout() << "aveConserveTest: problem in testConservation" << endl;
        return eekflag;
      }
    if (eekflag != 0) return eekflag;
  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  pout() << "average conservation test passed" << endl;
  return eekflag;
}
/***************/
