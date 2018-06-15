#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

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
#include "EBCoarseAverage.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"

#include "UsingNamespace.H"

/******************/
/******************/
int
makeLayout(DisjointBoxLayout& dbl1,
           const Box& domain);
/***************/
/***************/
int makeGeometry(Box& domain,
                 Real& dx);
/***************/
/***************/
int makeEBISL(EBISLayout& a_ebisl,
              const DisjointBoxLayout& a_grids,
              const Box& a_domain,
              const int& nghost);

/***************/
/***************/
int getError(LevelData<EBCellFAB>& a_errorCoar,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const Real& dxFine,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar,
             const Real& dxCoar);
/***************/
/***************/
int compareError(const LevelData<EBCellFAB>& a_fineError,
                 const EBISLayout& a_ebislFine,
                 const DisjointBoxLayout& a_gridsFine,
                 const Box& a_domainFine,
                 const LevelData<EBCellFAB>& a_coarError,
                 const EBISLayout& a_ebislCoar,
                 const DisjointBoxLayout& a_gridsCoar,
                 const Box& a_domainCoar);

/***************/
/***************/
Real exactFunc(const IntVect& a_iv,
               const Real& a_dx);
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
    const char* in_file = "ramp.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    //make the gometry.  this makes the first Chombo_EBIS
    Box domainFine, domainMedi, domainCoar;
    Real dxFine, dxMedi, dxCoar;
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domainFine,  dxFine);
    CH_assert(eekflag == 0);

    domainMedi = coarsen(domainFine, 2);
    domainCoar = coarsen(domainMedi, 2);
    dxMedi = 2.0*dxFine;
    dxCoar = 2.0*dxMedi;
    //make grids
    DisjointBoxLayout gridsFine, gridsMedi, gridsCoar;
    eekflag = makeLayout(gridsFine, domainFine);
    CH_assert(eekflag == 0);
    coarsen(gridsMedi, gridsFine, 2);
    coarsen(gridsCoar, gridsMedi, 2);

    ///create ebislayout
    int nghost = 2;
    EBISLayout ebislFine, ebislMedi, ebislCoar;
    eekflag = makeEBISL(ebislFine, gridsFine, domainFine, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislMedi, gridsMedi, domainMedi, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, nghost);
    if (eekflag != 0) return eekflag;

    int nvar = 1;
    EBCellFactory ebcellfactMedi(ebislMedi);
    EBCellFactory ebcellfactCoar(ebislCoar);
    LevelData<EBCellFAB> errorMedi(gridsMedi, nvar,
                                   IntVect::Zero,ebcellfactMedi);
    LevelData<EBCellFAB> errorCoar(gridsCoar, nvar,
                                   IntVect::Zero,ebcellfactCoar);

    eekflag = getError(errorMedi,
                       ebislFine, gridsFine, domainFine, dxFine,
                       ebislMedi, gridsMedi, domainMedi, dxMedi);
    if (eekflag != 0) return eekflag;

    eekflag = getError(errorCoar,
                       ebislMedi, gridsMedi, domainMedi, dxMedi,
                       ebislCoar, gridsCoar, domainCoar, dxCoar);
    if (eekflag != 0) return eekflag;

    eekflag =  compareError(errorMedi, ebislMedi, gridsMedi, domainMedi,
                            errorCoar, ebislCoar, gridsCoar, domainCoar);
    if (eekflag != 0) return eekflag;

  }//end scoping trick
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}
/***************/
/***************/
int getError(LevelData<EBCellFAB>& a_errorCoar,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const Real& a_dxFine,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar,
             const Real& a_dxCoar)
{
  int eekflag = 0;
  int nvar = 1;
  int nref = 2;
  EBCellFactory ebcellfactFine(a_ebislFine);
  EBCellFactory ebcellfactCoar(a_ebislCoar);
  LevelData<EBCellFAB> phiFine(a_gridsFine, nvar,
                               IntVect::Zero,ebcellfactFine);
  LevelData<EBCellFAB> phiCoar(a_gridsCoar, nvar,
                               IntVect::Zero,ebcellfactCoar);

  //fill phi fine and phiCExact
  //put phiCexact into a_errorCoar
  for (DataIterator dit = a_gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_gridsFine.get(dit()));
      EBCellFAB& phiFineFAB = phiFine[dit()];
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), a_dxFine);
          phiFineFAB(vof, 0) = rightAns;
        }

      ivsBox = IntVectSet(a_gridsCoar.get(dit()));
      EBCellFAB& phiCoarFAB = a_errorCoar[dit()];
      for (VoFIterator vofit(ivsBox, a_ebislCoar[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), a_dxCoar);
          phiCoarFAB(vof, 0) = rightAns;
        }
    }

  //average phiFine onto phiC;
  EBCoarseAverage ebAveOp(a_gridsFine, a_gridsCoar,
                          a_ebislFine, a_ebislCoar,
                          a_domainCoar, nref, nvar, Chombo_EBIS::instance());

  Interval zeroiv(0,0);
  ebAveOp.average(phiCoar, phiFine, zeroiv);
  //error = phiC - phiCExact
  for (DataIterator dit = a_gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      EBCellFAB& errorFAB = a_errorCoar[dit()];
      EBCellFAB& phiCoarFAB = phiCoar[dit()];

      errorFAB -= phiCoarFAB;
    }

  return eekflag;
}
/***************/
/***************/
int compareError(const LevelData<EBCellFAB>& a_errorFine,
                 const EBISLayout& a_ebislFine,
                 const DisjointBoxLayout& a_gridsFine,
                 const Box& a_domainFine,
                 const LevelData<EBCellFAB>& a_errorCoar,
                 const EBISLayout& a_ebislCoar,
                 const DisjointBoxLayout& a_gridsCoar,
                 const Box& a_domainCoar)
{
  EBCellFactory factFine(a_ebislFine);
  EBCellFactory factCoar(a_ebislCoar);
  int eekflag = 0;
  for (int itype = 0; itype < 3; itype++)
    {
      EBNormType::NormMode normtype;
      if (itype == 0)
        {
          normtype = EBNormType::OverBoth;
          pout() << endl << " Using all uncovered cells." << endl  << endl;
        }
      else if (itype == 1)
        {
          normtype = EBNormType::OverOnlyRegular;
          pout() << endl << " Using only regular cells." << endl << endl;
        }
      else
        {
          normtype = EBNormType::OverOnlyIrregular;
          pout() << endl << " Using only irregular cells." << endl << endl;
        }
      int comp = 0;
      for (int inorm = 0; inorm < 3; inorm++)
        {

          if (inorm == 0)
            {
              pout() << endl << " Using max norm." << endl;
            }
          else
            {
              pout() << endl << " Using L-" << inorm << "norm." << endl;
            }
          Real coarnorm = EBArith::norm(a_errorCoar,
                                        a_gridsCoar, a_ebislCoar,
                                        comp, inorm, normtype);
          Real finenorm = EBArith::norm(a_errorFine,
                                        a_gridsFine, a_ebislFine,
                                        comp, inorm, normtype);
          pout() << "Coarse Error Norm = " << coarnorm << endl;
          pout() << "Fine   Error Norm = " << finenorm << endl;
          if (Abs(finenorm) > 1.0e-8)
            {
              Real order = log(coarnorm/finenorm)/log((Real)2.0);
              pout() << "Order of scheme = " << order << endl;
            }
        }
    }

  return eekflag;
}

/***************/
/***************/
Real exactFunc(const IntVect& a_iv, const Real& a_dx)
{
  Real retval;
  Real probHi;
  ParmParse pp;
  pp.get("prob_hi",probHi);
  RealVect xloc;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xloc[idir] = a_dx*(Real(a_iv[idir]) + 0.5);
    }
  Real x = xloc[0]/probHi;
  Real y = xloc[1]/probHi;
  retval = 1.0 + x + y + x*x + y*y;
  //  retval = xloc[0]*xloc[0];
  return retval;
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
int makeGeometry(Box& a_domain,
                 Real& a_dx)
{
  int eekflag =  0;
  //parse input file
  ParmParse pp;
  RealVect origin = RealVect::Zero;
  Vector<int> n_cell(SpaceDim);
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

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  a_dx = (prob_hi-prob_lo[0])/n_cell[0];
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
      ebisPtr->define(a_domain, origin, a_dx, regserv);
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
      vectDx *= a_dx;

      GeometryShop workshop(ramp,0,vectDx);
      //this generates the new EBIS
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, workshop);
    }
  else
    {
      //bogus which_geom
      pout() << " bogus which_geom input = " << whichgeom;
      eekflag = 33;
    }

  return eekflag;
}
/***************/
/***************/
