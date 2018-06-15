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
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "SphereIF.H"
#include "MultiSphereIF.H"
#include "EBQuadCFInterp.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBAMRIO.H"
#include "EBLevelGrid.H"
#include "EBCoarsen.H"

#include "UsingNamespace.H"

/******************/
/******************/
int
makeLayout(DisjointBoxLayout& dbl1,
           const Box& domain);
/***************/
/***************/
int makeGeometry(Box& domain,
                 RealVect& dx);
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
             const RealVect& dxFine);
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
Real exactFunc(const IntVect&  a_iv,
               const RealVect& a_dx);

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

    const char* in_file = "halfquad.inputs";
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    //make the gometry.  this makes the first Chombo_EBIS
    Box domainFine, domainCoar;
    RealVect dxFine, dxCoar;

    //and defines it using a geometryservice
    eekflag =  makeGeometry(domainFine,  dxFine);
    CH_assert(eekflag == 0);

    domainCoar = coarsen(domainFine, 2);
    dxCoar = 2.0*dxFine;
    //make grids
    DisjointBoxLayout gridsFine, gridsCoar;
    eekflag = makeLayout(gridsFine, domainFine);
    CH_assert(eekflag == 0);
    coarsen(gridsCoar, gridsFine, 2);

    ///create ebislayout
    int nghost = 1;
    EBISLayout ebislFine, ebislCoar;
    eekflag = makeEBISL(ebislFine, gridsFine, domainFine, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, nghost);
    if (eekflag != 0) return eekflag;

    int nvar = 1;
    EBCellFactory ebcellfactFine(ebislFine);
    EBCellFactory ebcellfactCoar(ebislCoar);
    IntVect ivghost = nghost*IntVect::Unit;
    LevelData<EBCellFAB> errorFine(gridsFine, nvar,
                                   ivghost,ebcellfactFine);
    LevelData<EBCellFAB> errorCoar(gridsCoar, nvar,
                                   ivghost,ebcellfactCoar);

    eekflag = getError(errorCoar,ebislCoar, gridsCoar, domainCoar, dxCoar);
    if (eekflag != 0) return eekflag;

    eekflag = getError(errorFine,ebislFine, gridsFine, domainFine, dxFine);
    if (eekflag != 0) return eekflag;

    eekflag =  compareError(errorFine, ebislFine, gridsFine, domainFine,
                            errorCoar, ebislCoar, gridsCoar, domainCoar);
    if (eekflag != 0) return eekflag;

#if CH_SPACEDIM==2
    string fileFine("pltFineError.2d.hdf5");
    string fileCoar("pltCoarError.2d.hdf5");
#else
    string fileFine("pltFineError.3d.hdf5");
    string fileCoar("pltCoarError.3d.hdf5");
#endif

    outputError(errorFine,
                errorCoar,
                gridsFine,
                gridsCoar,
                domainFine,
                domainCoar,
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
int compareError(const LevelData<EBCellFAB>& a_errorFine,
                 const EBISLayout&           a_ebislFine,
                 const DisjointBoxLayout&    a_gridsFine,
                 const Box&                  a_domainFine,
                 const LevelData<EBCellFAB>& a_errorCoar,
                 const EBISLayout&           a_ebislCoar,
                 const DisjointBoxLayout&    a_gridsCoar,
                 const Box&                  a_domainCoar)
{
  //the norms here must go over ghost cells
  //too or the test is really silly
  BoxLayout blNormFine, blNormCoar;
  blNormFine.deepCopy(a_gridsFine);
  blNormCoar.deepCopy(a_gridsCoar);
  blNormFine.grow(1);
  blNormCoar.grow(1);
  blNormFine&= a_domainFine;
  blNormCoar&= a_domainCoar;
  blNormFine.close();
  blNormCoar.close();

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
                                        blNormCoar, a_ebislCoar,
                                        comp, inorm, normtype);
          Real finenorm = EBArith::norm(a_errorFine,
                                        blNormFine, a_ebislFine,
                                        comp, inorm, normtype);
          pout() << "Coarse Error Norm = " << coarnorm << endl;
          pout() << "Fine   Error Norm = " << finenorm << endl;
          if (Abs(finenorm) > 1.0e-12)
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
Real exactFunc(const IntVect& a_iv,
               const RealVect& a_dx)
{
  Real retval;
  Real probHi;
  ParmParse pp;
  pp.get("prob_hi",probHi);
  RealVect xloc,x,x2,x3;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xloc[idir] = a_dx[idir]*(Real(a_iv[idir]) + 0.5);
      x[idir] = xloc[idir]/probHi;
      x2[idir] = pow(x[idir],2);
      x3[idir] = pow(x[idir],3);
    }

  retval = 1.0;//constant termx
  //for (int idir = 0; idir < 1; idir++)
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      retval += x[idir];//linear terms
      retval += x2[idir];//quadratic terms
      retval += x3[idir];//cubic terms
      //      retval += 100.0*sin(x[idir]);//sin terms
      for (int jdir = 0; jdir < SpaceDim; jdir++)
        {
          if (idir!=jdir)
            {
              retval += x[ idir]*x[ jdir];//mixed linear    - linear    terms
              retval += x[ idir]*x2[jdir];//mixed linear    - quadratic terms
              retval += x2[idir]*x2[jdir];//mixed quadratic - quadratic terms
            }
        }
    }
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
makeLayout(DisjointBoxLayout& a_dbl,
           const Box& a_domainFine)
{
  //set up mesh refine object
  ParmParse pp;
  int eekflag= 0;
  int maxsize;
  pp.get("maxboxsize",maxsize);
  int bufferSize = 2;
  int blockFactor = 8;
  Real fillrat = 0.7;
  Box domainCoar = coarsen(a_domainFine, 2);
  Vector<int> refRat(2,2);
  BRMeshRefine meshRefObj(domainCoar, refRat, fillrat,
                          blockFactor, bufferSize, maxsize);

  Vector<Vector<Box> > oldMeshes(2);
  oldMeshes[0] = Vector<Box>(1,   domainCoar);
  oldMeshes[1] = Vector<Box>(1, a_domainFine);

  //set up coarse tags
  int nc = domainCoar.size(0);
  Box halfBox = domainCoar;
  halfBox.growHi(0, -nc/2);
  IntVectSet tags(halfBox);


  int baseLevel = 0;
  int topLevel = 0;
  Vector<Vector<Box> > newMeshes;
  meshRefObj.regrid(newMeshes, tags, baseLevel,
                    topLevel, oldMeshes);

  const Vector<Box>& vbox = newMeshes[1];
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vbox);
  if (eekflag != 0) return eekflag;
  a_dbl.define(vbox, procAssign);
  return eekflag;
}
/**********/
int makeGeometry(Box& a_domain,
                 RealVect& a_dx)
{
  int eekflag =  0;
  //parse input file
  ParmParse pp;
  int verbosity = 0;
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
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_dx[idir] = (prob_hi-prob_lo[idir])/n_cell[idir];
      origin[idir] = prob_lo[idir];
    }
  int whichgeom;
  pp.get("which_geom",whichgeom);
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  if (whichgeom == 0)
    {
      //allregular
      pout() << "all regular geometry" << endl;
      AllRegularService regserv;
      ebisPtr->define(a_domain, origin, a_dx[0], regserv);
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

      GeometryShop workshop(ramp,verbosity,a_dx);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, a_dx[0], workshop);
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

      bool     insideRegular = false;
      SphereIF implicit(sphereRadius,sphereCenter,insideRegular);

      GeometryShop workshop(implicit,verbosity,a_dx);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, a_dx[0], workshop);
    }
  else if (whichgeom == 6)
    {
      pout() << "multisphere geometry" << endl;
      int numSpheres;
      pp.get("num_spheres", numSpheres);
      Vector<Real>     radius(numSpheres);
      Vector<RealVect> center(numSpheres);
      for (int isphere = 0; isphere < numSpheres; isphere++)
        {
          char radiusString[80];
          char centerString[80];
          sprintf(radiusString, "sphere_radius_%d", isphere);
          sprintf(centerString, "sphere_center_%d", isphere);
          vector<Real> sphere_center(SpaceDim);
          Real sphereRadius;
          pp.get(radiusString, sphereRadius);
          pp.getarr(centerString,sphere_center, 0, SpaceDim);
          RealVect sphereCenter;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              sphereCenter[idir] = sphere_center[idir];
            }
          center[isphere] = sphereCenter;
          radius[isphere] = sphereRadius;
        }
      bool inside = false;
      MultiSphereIF spheres(radius, center, inside);
      GeometryShop workshop(spheres,verbosity,a_dx);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, a_dx[0], workshop);
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
int getError(LevelData<EBCellFAB>& a_errorFine,
             const EBISLayout& a_ebislFine,
             const DisjointBoxLayout& a_gridsFine,
             const Box& a_domainFine,
             const RealVect& a_dxFine)
{
  ParmParse pp;
  int eekflag = 0;
  int nvar = 1;
  int nghost = 1;
  int nref = 2;
  pp.get("ref_ratio",nref);
  int maxsize;
  pp.get("maxboxsize",maxsize);
  //generate coarse dx, domain
  Box domainCoar = coarsen(a_domainFine, nref);

  RealVect dxCoar = nref*a_dxFine;
  //make the coarse grids completely cover the domain
  Vector<Box> vbox;
  domainSplit(domainCoar, vbox,  maxsize);
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vbox);
  DisjointBoxLayout gridsCoar(vbox, procAssign);
  //make coarse ebisl
  EBISLayout ebislCoar;
  eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, nghost);
  if (eekflag != 0) return eekflag;
  // pout() << "getError dom coar = " << domainCoar << endl;;
  // pout() << "getError grids coar = " << gridsCoar << endl;;
  // pout() << "getError dom fine = " << a_domainFine << endl;;
  // pout() << "getError grids fine = " << a_gridsFine << endl;;

  //create data at both refinemenets
  EBCellFactory ebcellfactFine(a_ebislFine);
  EBCellFactory ebcellfactCoar(ebislCoar);
  LevelData<EBCellFAB> phiFine(a_gridsFine, nvar,
                               IntVect::Unit,ebcellfactFine);
  LevelData<EBCellFAB> phiCoar(gridsCoar, nvar,
                               IntVect::Unit,ebcellfactCoar);

  //fill phi fine and phiCExact
  //put phiFexact into a_errorFine and phiFine
  //this should make the error at the interior = 0
  for (DataIterator dit = a_gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      Box grownBox = grow(a_gridsFine.get(dit()), 1);
      grownBox &= a_domainFine;
      IntVectSet ivsBox(grownBox);
      phiFine[dit()].setVal(0.0);
      a_errorFine[dit()].setVal(0.0);
      EBCellFAB& phiFineFAB = phiFine[dit()];
      EBCellFAB& errFineFAB = a_errorFine[dit()];
      phiFineFAB.setCoveredCellVal(0.0,0);
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), a_dxFine);
          phiFineFAB(vof, 0) = rightAns;
          errFineFAB(vof, 0) = rightAns;
        }
    }

  for (DataIterator dit = gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      Box grownBox = grow(gridsCoar.get(dit()), 1);
      grownBox &= domainCoar;
      IntVectSet ivsBox(grownBox);
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


  LayoutData<IntVectSet> cfivs;
  EBArith::defineCFIVS(cfivs, a_gridsFine, a_domainFine);
  EBQuadCFInterp interpOp(a_gridsFine,  gridsCoar,
                          a_ebislFine,  ebislCoar,
                          domainCoar, nref, nvar,
                          cfivs, Chombo_EBIS::instance(), true);

  Interval interv(0,0);
  interpOp.interpolate(phiFine, phiCoar, interv);

  //error = phiF - phiFExact
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
  IntVect gv  = a_errorFine.ghostVect();

  ParmParse pp;

  writeEBHDF5(a_fileFine, a_gridsFine, a_errorFine, names,a_domainFine, dxFine, dt, time,replaceCovered, coveredValues,gv);
  writeEBHDF5(a_fileCoar, a_gridsCoar, a_errorCoar, names,a_domainCoar, dxCoar, dt, time,replaceCovered, coveredValues,gv);
#endif
}/***************/
/***************/
