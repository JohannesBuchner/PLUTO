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
#include "DebugDump.H"
#include "SlabService.H"
#include "EBDebugDump.H"
#include "EBCoarseAverage.H"
#include "EBLevelGrid.H"
#include "EBFastFR.H"

#include "UsingNamespace.H"

/*****************/
/*****************/
int fluxRegTest(void);
int nonEBFluxRegTest(void);
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

    eekflag = nonEBFluxRegTest();
    if (eekflag == 0)
      pout() << "noneb fluxregister test passed." << endl;
    else
      pout() << "noneb fluxregister test FAILED with error code "
           << eekflag << endl;

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

  Box domainCoar, domainFine;
  eekflag = makeGeometry(domainFine);
  if (eekflag != 0)
    return eekflag;
  domainCoar = coarsen(domainFine, nref);

  DisjointBoxLayout dblFine,dblCoar;
  eekflag = makeLayouts(dblCoar, dblFine, domainCoar, domainFine);

  Interval interv(0,0);

  EBISLayout ebislFine, ebislCoar;
  int nghost = 3;
  makeEBISL(ebislFine, dblFine, domainFine, nghost);
  makeEBISL(ebislCoar, dblCoar, domainCoar, nghost);

  EBCellFactory factFine(ebislFine);
  EBCellFactory factCoar(ebislCoar);
  IntVect ivghost = IntVect::Unit;

  LevelData<EBCellFAB> fineData(dblFine, 1,ivghost, factFine);
  LevelData<EBCellFAB> coarData(dblCoar, 1,ivghost, factCoar);

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
  //increment and decrement
  //flux registers with equal size fluxes
  for (DataIterator coarIt = coarData.dataIterator();
      coarIt.ok(); ++coarIt)
    {
      const Box&  boxCoar = dblCoar.get(coarIt());
      const EBISBox& ebisBox = ebislCoar[coarIt()];
      for (int idir = 0; idir < SpaceDim; idir++)
        {

          EBFaceFAB edgeFluxReg(ebisBox, boxCoar, idir, 1);
          edgeFluxReg.setVal(fluxVal);
          for (SideIterator sit; sit.ok(); ++sit)
            fluxReg.incrementCoarseBoth(edgeFluxReg, scale,
                                        coarIt(), interv, idir, sit());
        }
    }

  for (DataIterator fineIt = fineData.dataIterator();
      fineIt.ok(); ++fineIt)
    {
      const Box&  boxFine = dblFine.get(fineIt());
      const EBISBox& ebisBox = ebislFine[fineIt()];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB edgeFluxReg(ebisBox, boxFine, idir, 1);
          edgeFluxReg.setVal(fluxVal);
          for (SideIterator sit; sit.ok(); ++sit)
            fluxReg.incrementFineBoth(edgeFluxReg, scale,
                                      fineIt(),  interv, idir, sit());
        }
    }

  fluxReg.reflux(coarData, interv, scale);

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
          if ((rmax > 1.0e-10)||(rmin > 1.0e-10))
            {
              pout() << "EBFastFR failed the all Regular Conservation test at "
                     << " vof = " << vof <<endl;
              eekflag = 42;
              return eekflag;
            }
        }
    }

  return eekflag;
}

int nonEBFluxRegTest()
{
  int nref = 2;
  int eekflag = 0;
  ParmParse pp;

  Box domainCoar, domainFine;
  eekflag = makeGeometry(domainFine);
  if (eekflag != 0)
    return eekflag;
  domainCoar = coarsen(domainFine, nref);

  DisjointBoxLayout dblFine,dblCoar;
  eekflag = makeLayouts(dblCoar, dblFine, domainCoar, domainFine);

  Interval interv(0,0);
  IntVect ivghost = IntVect::Unit;

  LevelData<FArrayBox> fineData(dblFine, 1,ivghost);
  LevelData<FArrayBox> coarData(dblCoar, 1,ivghost);

  //set data and flux registers to zero
  for (DataIterator coarIt = coarData.dataIterator();
      coarIt.ok(); ++coarIt)
    coarData[coarIt()].setVal(0.);
  for (DataIterator fineIt = fineData.dataIterator();
      fineIt.ok(); ++fineIt)
    fineData[fineIt()].setVal(0.);

  ProblemDomain probDomainFine(domainFine);
  ProblemDomain probDomainCoar(domainCoar);
  LevelFluxRegister fluxReg(dblFine, dblCoar, domainFine, nref, 1);

  fluxReg.setToZero();
  //loop through directions
  Real scale = 1.0;
  Real fluxVal = 4.77;
  //increment and decrement
  //flux registers with equal size fluxes
  for (DataIterator coarIt = coarData.dataIterator();
      coarIt.ok(); ++coarIt)
    {
      const Box&  boxCoar = dblCoar.get(coarIt());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Box edgeBox = surroundingNodes(boxCoar, idir);
          FArrayBox edgeFluxReg(edgeBox, 1);
          edgeFluxReg.setVal(fluxVal);
          fluxReg.incrementCoarse(edgeFluxReg, scale,
                                  coarIt(), interv, interv, idir);
        }
    }

  for (DataIterator fineIt = fineData.dataIterator();
      fineIt.ok(); ++fineIt)
    {
      const Box&  boxFine = dblFine.get(fineIt());
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Box edgeBox = surroundingNodes(boxFine, idir);
          FArrayBox edgeFluxReg(edgeBox, 1);
          edgeFluxReg.setVal(fluxVal);
          for (SideIterator sit; sit.ok(); ++sit)
            {
              fluxReg.incrementFine(edgeFluxReg, scale,
                                    fineIt(),  interv,
                                    interv, idir, sit());

            }
        }
    }

  fluxReg.reflux(coarData, interv, interv, scale);

  DataIterator datIt = coarData.dataIterator();
  for (datIt.reset(); datIt.ok(); ++datIt)
    {
      const FArrayBox& data = coarData[datIt()];
      Real rmax = 0.;
      Real rmin = 0.;
      for (BoxIterator boxit(dblCoar.get(datIt())); boxit.ok(); ++boxit)
        {
          const IntVect& vof = boxit();
          rmax = Max(rmax, Abs(data(vof, 0)));
          rmin = Min(rmax, Abs(data(vof, 0)));
          if ((rmax > 1.0e-10)||(rmin > 1.0e-10))
            {
              pout() << "EBFastFR failed the all Regular Conservation test at "
                     << " vof = " << vof <<endl;
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

  //  pout() << "all regular geometry" << endl;
  AllRegularService regserv;
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain, origin, dx, regserv, 2048);

  return eekflag;
}
/*****************/
/*****************/
