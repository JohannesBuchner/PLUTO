#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves mon aug 27, 2001

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
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "EBPWLFineInterp.H"
#include "DebugDump.H"
#include "SlabService.H"
#include "EBAMRIO.H"

#include "UsingNamespace.H"

/******************/
/******************/
int
makeLayoutFull(DisjointBoxLayout& dbl1,
               const Box& domain);
int
makeLayoutPart(DisjointBoxLayout& dbl1,
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
             const Real& dxCoar,
             const int& a_refRatio);
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
void outputError(const LevelData<EBCellFAB>&   a_errorFine,
                 const LevelData<EBCellFAB>&   a_errorCoar,
                 const DisjointBoxLayout&      a_gridsFine,
                 const DisjointBoxLayout&      a_gridsCoar,
                 const Box&                    a_domainFine,
                 const Box&                    a_domainCoar,
                 const string& a_fileFine,
                 const string& a_fileCoar);
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
  {//begin forever present scoping trick

    const char* in_file = "pwlinterp.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    //make the gometry.  this makes the first Chombo_EBIS
    Box domainFineFull, domainCoarFull;
    Box domainFinePart, domainCoarPart;
    Real dxFineFull, dxCoarFull;
    Real dxFinePart, dxCoarPart;
    //and defines it using a geometryservice
    eekflag =  makeGeometry(domainFinePart,  dxFinePart);
    CH_assert(eekflag == 0);
    int nlevels = 2;
    Vector<int> refRatio(nlevels);
    pp.getarr("ref_ratio", refRatio,0,nlevels);

    domainCoarPart = coarsen(domainFinePart, 2);
    domainFineFull = coarsen(domainFinePart, refRatio[0]);
    domainCoarFull = coarsen(domainFineFull, 2);

    dxCoarPart =         2.0*dxFinePart;
    dxFineFull = refRatio[0]*dxFinePart;
    dxCoarFull =         2.0*dxFineFull;

    //make grids
    DisjointBoxLayout gridsFinePart, gridsCoarPart;
    DisjointBoxLayout gridsFineFull, gridsCoarFull;

    eekflag = makeLayoutPart(gridsCoarPart, domainCoarFull);
    eekflag = makeLayoutFull(gridsCoarFull, domainCoarFull);
    CH_assert(eekflag == 0);
    refine(gridsFinePart, gridsCoarPart, 2);
    refine(gridsFineFull, gridsCoarFull, 2);

    ///create ebislayout
    int nghost = 0;
    EBISLayout ebislFinePart, ebislCoarPart;
    EBISLayout ebislFineFull, ebislCoarFull;

    eekflag = makeEBISL(ebislFinePart, gridsFinePart, domainFinePart, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislCoarPart, gridsCoarPart, domainCoarPart, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislFineFull, gridsFineFull, domainFineFull, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislCoarFull, gridsCoarFull, domainCoarFull, nghost);
    if (eekflag != 0) return eekflag;

    int nvar = 1;
    EBCellFactory ebcellfactFinePart(ebislFinePart);
    EBCellFactory ebcellfactCoarPart(ebislCoarPart);
    LevelData<EBCellFAB> errorFinePart(gridsFinePart, nvar,
                                   IntVect::Zero,ebcellfactFinePart);
    LevelData<EBCellFAB> errorCoarPart(gridsCoarPart, nvar,
                                   IntVect::Zero,ebcellfactCoarPart);

    eekflag = getError(errorFinePart,
                       ebislFinePart, gridsFinePart, domainFinePart, dxFinePart,
                       ebislFineFull, gridsFineFull, domainFineFull, dxFineFull,
                       refRatio[0]);
    if (eekflag != 0) return eekflag;

    eekflag = getError(errorCoarPart,
                       ebislCoarPart, gridsCoarPart, domainCoarPart, dxCoarPart,
                       ebislCoarFull, gridsCoarFull, domainCoarFull, dxCoarFull,
                       refRatio[0]);
    if (eekflag != 0) return eekflag;

    eekflag =  compareError(errorFinePart, ebislFinePart, gridsFinePart, domainFinePart,
                            errorCoarPart, ebislCoarPart, gridsCoarPart, domainCoarPart);
    if (eekflag != 0) return eekflag;


#if CH_SPACEDIM==2
    string fileFine("pltFineError.2d.hdf5");
    string fileCoar("pltCoarError.2d.hdf5");
#else
    string fileFine("pltFineError.3d.hdf5");
    string fileCoar("pltCoarError.3d.hdf5");
#endif

    outputError(errorFinePart,
                errorCoarPart,
                gridsFinePart,
                gridsCoarPart,
                domainFinePart,
                domainCoarPart,
                fileFine,
                fileCoar);

  }//end omnipresent scoping trick
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}
/***************/
/***************/
int getError(LevelData<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const Real& a_dxFine,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar,
             const Real& a_dxCoar,
             const int& a_refRatio)
{
  int eekflag = 0;
  int nvar = 1;

  EBCellFactory ebcellfactFine(a_ebislFine);
  EBCellFactory ebcellfactCoar(a_ebislCoar);
  LevelData<EBCellFAB> phiFine(a_gridsFine, nvar,
                               IntVect::Zero,ebcellfactFine);
  LevelData<EBCellFAB> phiCoar(a_gridsCoar, nvar,
                               IntVect::Zero,ebcellfactCoar);

  //fill phi fine and phiCExact
  //put phiFexact into a_errorFine
  for (DataIterator dit = a_gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_gridsFine.get(dit()));
      phiFine[dit()].setVal(0.0);
      EBCellFAB& phiFineFAB = a_errorFine[dit()];
      phiFineFAB.setCoveredCellVal(0.0,0);
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), a_dxFine);
          phiFineFAB(vof, 0) = rightAns;
        }
    }

  for (DataIterator dit = a_gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_gridsCoar.get(dit()));
      EBCellFAB& phiCoarFAB = phiCoar[dit()];
      phiCoarFAB.setCoveredCellVal(0.0,0);
      for (VoFIterator vofit(ivsBox, a_ebislCoar[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), a_dxCoar);
          phiCoarFAB(vof, 0) = rightAns;
        }
    }

  EBPWLFineInterp interpOp(a_gridsFine, a_gridsCoar,
                           a_ebislFine, a_ebislCoar,
                           a_domainCoar, a_refRatio, nvar);
  Interval zeroiv(0,0);
  interpOp.interpolate(phiFine, phiCoar, zeroiv);
  //error = phiC - phiCExact
  for (DataIterator dit = a_gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      EBCellFAB& errorFAB = a_errorFine[dit()];
      EBCellFAB& phiFineFAB = phiFine[dit()];

      errorFAB -= phiFineFAB;
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
              Real order = log(coarnorm/finenorm)/log(2.0);
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
  Real z = 0.0;
  if (SpaceDim==3)
    {
      z = xloc[2]/probHi;
    }
  //retval = 1.0 + x*x + y*y + z*z + x*x*x + y*y*y + z*z*z;
  retval = cos(20*x)*cos(20*y)*cos(20*z);
  //  if (x>0.74)
  //     {
  //       retval = 1;
  //     }
  //   else
  //     {
  //       retval = 0;
  //     }
  //retval = 1.0 + x*x + y*y;
  //  retval = 2.*xloc[0];
  //retval  = x + y;
  //retval  = x*x;
  //  retval  = x;
  //retval  = 1.0;
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
  const CH_XD::EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());
  ebisPtr->fillEBISLayout(a_ebisl, a_grids, a_domain, a_nghost);
  return 0;
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
makeLayoutPart(DisjointBoxLayout& a_dbl,
               const Box& a_domainCoarsest)
{
  //  Vector<int> vecRefRat(2, 2);
  ParmParse pp;
  int eekflag= 0;
  int blockFactor, bufferSize, maxSize;
  Real fillRat;
  pp.get("block_factor", blockFactor);
  pp.get("buffer_size", bufferSize);
  pp.get("maxboxsize", maxSize);
  pp.get("fill_ratio", fillRat);
  int nlevels =  2;
  Vector<int> vecRefRat(nlevels);
  pp.getarr("ref_ratio", vecRefRat,0,nlevels);
  BRMeshRefine mesher(a_domainCoarsest, vecRefRat,
                      fillRat, blockFactor, bufferSize, maxSize);

  int topLevel  =  0;
  int baseLevel =  0;
  //tags at base level
  IntVectSet tags;
  eekflag = makeTags(tags, a_domainCoarsest);
  if (eekflag < 0) return eekflag;

  Vector<Vector<Box> > oldMeshes(nlevels);
  oldMeshes[0] = Vector<Box>(1, a_domainCoarsest);
  Box finerDomain = a_domainCoarsest;
  for (int ilev = 1; ilev < nlevels; ilev++)
    {
      finerDomain.refine(vecRefRat[ilev]);
      oldMeshes[ilev] = Vector<Box>(1, finerDomain);
    }
  Vector<Vector<Box> > newMeshes;
  mesher.regrid(newMeshes, tags, baseLevel, topLevel, oldMeshes);

  Vector<int> procAssign;
  eekflag = LoadBalance(procAssign, newMeshes[1]);
  if (eekflag != 0) return eekflag;
  a_dbl.define(newMeshes[1], procAssign);
  int iverbose;
  pp.get("verbose", iverbose);
  if (iverbose == 1)
    {
      pout() << "the grids coarser domain= " << a_domainCoarsest
             << " = " << a_dbl << endl;
    }

  return eekflag;
}
/***************/
/***************/
int
makeLayoutFull(DisjointBoxLayout& a_dbl,
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

  int nlevels = 2;
  Vector<int> refRatio(nlevels);
  pp.getarr("ref_ratio", refRatio,0,nlevels);

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

  //now coarsen a_domain so it's the finest domain (level 1, fine hierarchy)
  a_domain.refine(refRatio[0]*2);

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  a_dx = (prob_hi-prob_lo[0])/n_cell[0];
  a_dx /= 2*refRatio[0];//so it's dx on the finest domain (level 1, fine hierarchy)

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
      CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
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
      CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_domain, origin, a_dx, workshop);
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
      CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      RealVect origin = RealVect::Zero;
      ebisPtr->define(a_domain, origin, a_dx, slab);
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
void outputError(const LevelData<EBCellFAB>&   a_errorFine,
                 const LevelData<EBCellFAB>&   a_errorCoar,
                 const DisjointBoxLayout&      a_gridsFine,
                 const DisjointBoxLayout&      a_gridsCoar,
                 const Box&                    a_domainFine,
                 const Box&                    a_domainCoar,
                 const string& a_fileFine,
                 const string& a_fileCoar)
{
#ifdef CH_USE_HDF5
  int nvar = 1;
  Vector<string> names(1, string("var0"));
  bool replaceCovered = true;
  Vector<Real> coveredValues(nvar, 0.0);
  //values that don't matter in output file
  Real dxFine = 1.0;
  Real dxCoar = 1.0;
  Real time   = 1.0;
  Real dt     = 1.0;

  ParmParse pp;

  writeEBHDF5(a_fileFine, a_gridsFine, a_errorFine, names,a_domainFine, dxFine, dt, time,replaceCovered, coveredValues);
  writeEBHDF5(a_fileCoar, a_gridsCoar, a_errorCoar, names,a_domainCoar, dxCoar, dt, time,replaceCovered, coveredValues);
#endif
}/***************/
/***************/
