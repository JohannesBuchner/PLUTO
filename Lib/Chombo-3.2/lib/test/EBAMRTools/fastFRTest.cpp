#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves mon oct 29, 2001

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
#include "DebugDump.H"
#include "SlabService.H"
#include "EBDebugDump.H"
#include "EBCoarseAverage.H"
#include "EBFastFR.H"
#include "UsingNamespace.H"

/*****************/
/*****************/
int fluxRegTest(void);
/******************/
/******************/
int
makeLayouts(DisjointBoxLayout& a_dblCoar,
            DisjointBoxLayout& a_dblFine,
            const Box& a_domainCoar,
            const Box& a_domainFine);
/***************/
/***************/
int makeGeometry(Box& domain);
/***************/
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& nghost);

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int eekflag=0;
  {//begin forever present scoping trick

    const char* in_file = "fluxreg.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    eekflag = fluxRegTest();
    if (eekflag == 0)
      pout() << "ebfluxregister test passed." << endl;
    else
      pout() << "ebfluxregister test FAILED with error code "
           << eekflag << endl;

  }//end omnipresent scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}

/*****************/
/*****************/

int fluxRegTest()
{
  int nref = 2;
  int eekflag = 0;
  ParmParse pp;
  int iverbose;
  pp.get("verbose", iverbose);
  bool verbose = (iverbose==1);

  Box domainCoar, domainFine;
  eekflag = makeGeometry(domainFine);
  if (eekflag != 0)
    return eekflag;
  domainCoar = coarsen(domainFine, nref);

  DisjointBoxLayout dblFine,dblCoar;
  eekflag = makeLayouts(dblCoar, dblFine, domainCoar, domainFine);

  Interval interv(0,0);
  if (verbose)
    {
      pout()  << "dblCoarCell is" << endl;
      dumpDBL(&dblCoar);
      pout()    << "dblFineCell is" << endl;
      dumpDBL(&dblFine);
    }
  EBISLayout ebislFine, ebislCoar;
  int nghost = 3;
  makeEBISL(ebislFine, dblFine, domainFine, nghost);
  makeEBISL(ebislCoar, dblCoar, domainCoar, nghost);

  EBCellFactory factFine(ebislFine);
  EBCellFactory factCoar(ebislCoar);
  IntVect ivghost = IntVect::Unit;

  LevelData<EBCellFAB> fineData(dblFine, 1,ivghost, factFine);
  LevelData<EBCellFAB> coarData(dblCoar, 1,ivghost, factCoar);

  LevelData<EBCellFAB> extraDense(dblCoar, 1,ivghost, factCoar);

  //set data and flux registers to zero
  for (DataIterator coarIt = coarData.dataIterator();
      coarIt.ok(); ++coarIt)
      coarData[coarIt()].setVal(0.);
  for (DataIterator fineIt = fineData.dataIterator();
      fineIt.ok(); ++fineIt)
    fineData[fineIt()].setVal(0.);

  ProblemDomain probDomainFine(domainFine);
  ProblemDomain probDomainCoar(domainCoar);
  EBLevelGrid eblgFine(dblFine, ebislFine, probDomainFine);
  EBLevelGrid eblgCoar(dblCoar, ebislCoar, probDomainCoar);
  EBFastFR fluxReg(eblgFine, eblgCoar, nref, 1);

  fluxReg.setToZero();
  //loop through directions
  Real scale = 1.0;
  Real fluxVal = 4.77;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //increment and decrement
      //flux registers with equal size fluxes
      for (DataIterator coarIt = coarData.dataIterator();
          coarIt.ok(); ++coarIt)
        {
          const Box&  boxCoar = dblCoar.get(coarIt());
          const EBISBox& ebisBox = ebislCoar[coarIt()];
          IntVectSet ivsBC(boxCoar);
          EBFaceFAB edgeFluxReg(ebisBox, boxCoar, idir, 1);
          edgeFluxReg.setVal(fluxVal);
          for (SideIterator sit; sit.ok(); ++sit)
            fluxReg.incrementCoarseBoth(edgeFluxReg, scale,
                                        coarIt(), interv, idir, sit());
        }
      for (DataIterator fineIt = fineData.dataIterator();
          fineIt.ok(); ++fineIt)
        {
          const Box&  boxFine = dblFine.get(fineIt());
          const EBISBox& ebisBox = ebislFine[fineIt()];
          IntVectSet ivsBF(boxFine);
          EBFaceFAB edgeFluxReg(ebisBox, boxFine, idir, 1);
          edgeFluxReg.setVal(fluxVal);

          for (SideIterator sit; sit.ok(); ++sit)
            fluxReg.incrementFineBoth(edgeFluxReg, scale,
                                      fineIt(),  interv, idir, sit());
        }
    }

  //reflux what ought to be zero into zero and the result should be zero
  //except where the coarse-fine boundary gets crossed by the embedded
  //boundary.  That should get fixed by the extramass thing.
  fluxReg.reflux(coarData, interv, scale);

  for (DataIterator coarIt = coarData.dataIterator();
      coarIt.ok(); ++coarIt)
    extraDense[coarIt()].setVal(0.);

  // now add extra density to soltuion
  //in the end the solution should return to zero
  fluxReg.incrementDensityArray(extraDense, interv, scale);
  for (DataIterator coarIt = coarData.dataIterator();
      coarIt.ok(); ++coarIt)
    coarData[coarIt()] += extraDense[coarIt()];

  DataIterator datIt = coarData.dataIterator();
  for (datIt.reset(); datIt.ok(); ++datIt)
    {
      const EBCellFAB& data = coarData[datIt()];
      IntVectSet ivsBox(dblCoar.get(datIt()));
      const EBISBox& ebisBox = ebislCoar[datIt()];
      Real rmax = 0.;
      Real rmin = 0.;
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          rmax = Max(rmax, Abs(data(vof, 0)));
          rmin = Min(rmax, Abs(data(vof, 0)));

#ifdef CH_USE_FLOAT
          Real tolerance = 1.0e-6;
#else
          Real tolerance = 1.0e-10;
#endif
          if ((rmax > tolerance)||(rmin > tolerance))
            {
              pout() << "fastFRTest failed the test at  "
                     << " vof = " << vof << endl;
              pout() << "   rmax: "      << rmax      << ", or"
                     <<   " rmin: "      << rmin      << " >"
                     <<   " tolerance: " << tolerance << ")" << endl;
              eekflag = 42;
              return eekflag;
            }
        }
    }

  return eekflag;
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
makeLayouts(DisjointBoxLayout& a_dblCoar,
            DisjointBoxLayout& a_dblFine,
            const Box& a_domainCoar,
            const Box& a_domainFine)
{
  //set up mesh refine object
  int eekflag= 0;
  Vector<Box> vboxCoar;
  Vector<Box> vboxFine;
#if 1
  ParmParse pp;
  int maxsize;
  pp.get("maxboxsize",maxsize);
  int bufferSize = 1;
  int blockFactor = 2;
  Real fillrat = 0.75;
  Vector<int> refRat(2,2);
  BRMeshRefine meshRefObj(a_domainCoar, refRat, fillrat,
                          blockFactor, bufferSize, maxsize);

  Vector<Vector<Box> > oldMeshes(2);
  oldMeshes[0] = Vector<Box>(1, a_domainCoar);
  oldMeshes[1] = Vector<Box>(1, a_domainFine);

  //set up coarse tags
  int nc = a_domainCoar.size(0);
  int nmi = nc/2;//16
  int nqu = nc/4;//8
  int ntf = (nc*3)/4;  //24
  int nte = (nc*3)/8; //12
  int nfe = (nc*5)/8; //20
#if (CH_SPACEDIM ==2)
  Box boxf1(IntVect(0, nqu), IntVect(nmi-1,ntf-1));
  Box boxf2(IntVect(nmi,nte), IntVect(ntf-1,nfe-1));
  Box boxf3(IntVect(nqu,0  ), IntVect(nfe-1,nqu-1));
  Box boxf4(IntVect(nfe,nqu), IntVect(nc -1,nte-1));
  //  Box boxf5(IntVect(nqu,nqu), IntVect(ntf-1,ntf-1));
#else
  Box boxf1(IntVect(0, nqu,nqu), IntVect(nmi-1,ntf-1,ntf-1));
  Box boxf2(IntVect(nmi,nte,nte), IntVect(ntf-1,nfe-1,nfe-1));
  Box boxf3(IntVect(nqu,0,0  ), IntVect(nfe-1,nqu-1,nqu-1));
  Box boxf4(IntVect(nfe,nqu,nqu), IntVect(nc -1,nte-1,nte-1));
  //  Box boxf5(IntVect(nqu,nqu), IntVect(ntf-1,ntf-1))
#endif
  IntVectSet tags;
  tags |= boxf1;
  tags |= boxf2;
  tags |= boxf3;
  tags |= boxf4;

  //  tags |= boxf5;
  int baseLevel = 0;
  int topLevel = 0;
  Vector<Vector<Box> > newMeshes;
  meshRefObj.regrid(newMeshes, tags, baseLevel,
                    topLevel, oldMeshes);

  vboxCoar = newMeshes[0];
  vboxFine = newMeshes[1];
#else
  vboxCoar = Vector<Box>(1, a_domainCoar);
  Box shrunkFine = a_domainFine;
  int nc = Max(a_domainCoar.size(0)/4, 2);
  shrunkFine.grow(-nc);
  vboxFine = Vector<Box>(1, shrunkFine);
#endif

  Vector<int>  procAssignFine, procAssignCoar;
  eekflag = LoadBalance(procAssignFine,vboxFine);
  if (eekflag != 0) return eekflag;
  eekflag = LoadBalance(procAssignCoar,vboxCoar);
  if (eekflag != 0) return eekflag;

  a_dblFine.define(vboxFine, procAssignFine);
  a_dblCoar.define(vboxCoar, procAssignCoar);
  return eekflag;
}
/**********/
/**********/
int makeGeometry(Box& a_domain)
{
  int eekflag =  0;
  //parse input file
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Vector<int> n_cell(SpaceDim);
  pp.getarr("n_cell",n_cell,0,SpaceDim);
  Real probhi;
  pp.get("prob_hi", probhi);
  Real dx = probhi/n_cell[0];

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

  int iverbose;
  pp.get("verbose", iverbose);
  bool verbose = (iverbose==1);

  int whichgeom;
  pp.get("which_geom",whichgeom);
  if (whichgeom == 0)
    {
      //allregular
      if (verbose)
        pout() << "all regular geometry" << endl;
      AllRegularService regserv;
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, dx, regserv, 2048);
    }
  else if (whichgeom == 1)
    {
      if (verbose)
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
      ebisPtr->define(a_domain, origin, dx, workshop, 2048);
    }
  else if (whichgeom == 2)
    {
      if (verbose)
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
      ebisPtr->define(a_domain, origin, dx, slab,2048);
    }
  else
    {
      //bogus which_geom
      pout() << " bogus which_geom input = " << whichgeom;
      eekflag = 33;
    }
  return eekflag;
}
/*****************/
/*****************/
