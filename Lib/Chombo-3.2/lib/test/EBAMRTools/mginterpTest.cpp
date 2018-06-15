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
#include "SphereIF.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelDataOps.H"
#include "EBMGInterp.H"
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
  retval = cos(10*x)*cos(10*y)*cos(10*z);
  //retval = x + y;
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

int getError(LevelData<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const Real& a_dxFine,
             int a_doLinear)
{
  int eekflag = 0;
  int nvar = 1;
  IntVect giv = IntVect::Unit;
  DisjointBoxLayout gridsCoar;
  coarsen(gridsCoar, a_gridsFine, 2);
  Box domainCoar = coarsen(a_domainFine, 2);
  EBISLayout  ebislCoar;
  Real dxCoar = 2*a_dxFine;

  eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, 2);

  EBCellFactory ebcellfactFine(a_ebislFine);
  EBCellFactory ebcellfactCoar( ebislCoar);
  LevelData<EBCellFAB> phiFine(a_gridsFine, nvar,giv ,ebcellfactFine);
  LevelData<EBCellFAB> phiCoar(  gridsCoar, nvar,giv ,ebcellfactCoar);

  //set phi to zero because this is all incremental.
  EBLevelDataOps::setVal(phiFine, 0.0);

  //fill phi fine and phiCExact
  //put phiFexact into a_errorFine
  for (DataIterator dit = a_gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      IntVectSet ivsBox(a_gridsFine.get(dit()));
      phiFine[dit()].setVal(0.0);
      EBCellFAB& errorFAB = a_errorFine[dit()];
      errorFAB.setCoveredCellVal(0.0,0);
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), a_dxFine);
          errorFAB(vof, 0) = rightAns;
        }
    }

  for (DataIterator dit = gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      IntVectSet ivsBox(gridsCoar.get(dit()));
      EBCellFAB& phiCoarFAB = phiCoar[dit()];
      phiCoarFAB.setCoveredCellVal(0.0,0);
      for (VoFIterator vofit(ivsBox, ebislCoar[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), dxCoar);
          phiCoarFAB(vof, 0) = rightAns;
        }
    }

  int refRat = 2;
  EBMGInterp interpOp(a_gridsFine, gridsCoar,
                      a_ebislFine, ebislCoar,
                      domainCoar, refRat, nvar,
                      Chombo_EBIS::instance(),
                      giv,  true, true);
  Interval zeroiv(0,0);
  if (a_doLinear)
    interpOp.pwlInterp(phiFine, phiCoar, zeroiv);
  else
    interpOp.pwcInterp(phiFine, phiCoar, zeroiv);

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
  else if (whichgeom == 5)
    {
      pout() << "sphere geometry" << endl;
      vector<Real> sphere_center(SpaceDim);
      pp.getarr("sphere_center",sphere_center, 0, SpaceDim);
      RealVect sphereCenter;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          sphereCenter[idir] = sphere_center[idir];
        }
      Real sphereRadius;
      pp.get("sphere_radius", sphereRadius);
      bool fluidinside;
      pp.get("fluid_inside", fluidinside);
      SphereIF sphereIF(sphereRadius, sphereCenter, fluidinside);
      int verbosity = 0;
      GeometryShop workshop(sphereIF,verbosity, a_dx*RealVect::Unit);
      //this generates the new EBIS
      CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      RealVect origin = RealVect::Zero;
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
int
main(int argc, char** argv)

{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int eekflag = 0;
  {//begin forever present scoping trick

    const char* in_file = "mginterp.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    //make the gometry.  this makes the first Chombo_EBIS
    Box domainFine, domainCoar;
    Real dxFine, dxCoar;
    eekflag =  makeGeometry(domainFine, dxFine);
    //and defines it using a geometryservice
    CH_assert(eekflag == 0);
    dxCoar = 2*dxFine;
    domainCoar = coarsen(domainFine, 2);

    DisjointBoxLayout gridsFine, gridsCoar;
    eekflag = makeLayoutFull(gridsFine, domainFine);

    coarsen(gridsCoar, gridsFine, 2);

    ///create ebislayout
    int nghost = 2;
    EBISLayout ebislFine, ebislCoar;

    eekflag = makeEBISL(ebislFine, gridsFine, domainFine, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, nghost);
    if (eekflag != 0) return eekflag;

    int nvar = 1;
    EBCellFactory ebcellfactFine(ebislFine);
    EBCellFactory ebcellfactCoar(ebislCoar);
    LevelData<EBCellFAB> errorFine(gridsFine, nvar, IntVect::Zero,ebcellfactFine);
    LevelData<EBCellFAB> errorCoar(gridsCoar, nvar, IntVect::Zero,ebcellfactCoar);

    for (int doLinear = 0; doLinear < 2; doLinear++)
      {
        if (doLinear == 0)
          {
            pout() << "doing piecewise constant interpolation" << endl;
          }
        else
          {
            pout() << "doing piecewise linear interpolation" << endl;
          }
        eekflag = getError(errorFine, ebislFine, gridsFine, domainFine, dxFine, doLinear);
        if (eekflag != 0) return eekflag;

        eekflag = getError(errorCoar, ebislCoar, gridsCoar, domainCoar, dxCoar, doLinear);
        if (eekflag != 0) return eekflag;

        eekflag =  compareError(errorFine, ebislFine, gridsFine, domainFine,
                                errorCoar, ebislCoar, gridsCoar, domainCoar);
        if (eekflag != 0) return eekflag;


#if CH_SPACEDIM==2
        string fileFine("pltFineError.pwc.2d.hdf5");
        string fileCoar("pltCoarError.pwc.2d.hdf5");
        if (doLinear == 1)
          {
            fileFine = string("pltFineError.pwl.2d.hdf5");
            fileCoar = string("pltCoarError.pwl.2d.hdf5");
          }

#else
        string fileFine("pltFineError.3d.hdf5");
        string fileCoar("pltCoarError.3d.hdf5");
        if (doLinear == 1)
          {
            fileFine = string("pltFineError.pwl.3d.hdf5");
            fileCoar = string("pltCoarError.pwl.3d.hdf5");
          }
#endif

        outputError(errorFine,
                    errorCoar,
                    gridsFine,
                    gridsCoar,
                    domainFine,
                    domainCoar,
                    fileFine,
                    fileCoar);
      }
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
