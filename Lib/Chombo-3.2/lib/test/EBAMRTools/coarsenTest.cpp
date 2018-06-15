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
#include "EBCoarsen.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "PlaneIF.H"
#include "SphereIF.H"
#include "UnionIF.H"
#include "CH_Attach.H"
#include "ComplementIF.H"
#include "TransformIF.H"
#include "TiltedCylinderIF.H"
#include "IntersectionIF.H"
#include "EBAMRIO.H"
#include "DebugDump.H"
#include "EBLevelDataOps.H"

#include "UsingNamespace.H"

/******************/
/******************/
int
makeLayoutFine(DisjointBoxLayout& dbl1,
               const Box& domain);
/******************/
/******************/
int
makeLayoutCoar(DisjointBoxLayout& dbl1,
               const Box& domain);
/***************/
/***************/
int makeGeometry(const Box& domain,
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
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar,
             const RealVect& dxCoar,
             const int& a_nref);
/***************/
/***************/
int compareError(const LevelData<EBCellFAB>& a_fineError,
                 const EBISLayout& a_ebislFine,
                 const DisjointBoxLayout& a_gridsFine,
                 const Box& a_domainFine,
                 const LevelData<EBCellFAB>& a_coarError,
                 const EBISLayout& a_ebislCoar,
                 const DisjointBoxLayout& a_gridsCoar,
                 const Box& a_domainCoar,
                 const int& a_nref);
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
Real exactFunc(const IntVect&  a_iv,
               const RealVect& a_dx);
/***************/
/***************/
BaseIF* makeChamber(const Real& radius,
                    const Real& thick,
                    const Real& offset,
                    const Real& height);

BaseIF* makePlate(const Real& height,
                  const Real& thick,
                  const Real& radius,
                  const int&  doHoles,
                  const Real& holeRadius,
                  const Real& holeSpace);

BaseIF* makeVanes(const int&      num,
                  const Real&     thick,
                  const RealVect& normal,
                  const Real&     innerRadius,
                  const Real&     outerRadius,
                  const Real&     offset,
                  const Real&     height);

BaseIF* makeVane(const Real&     thick,
                 const RealVect& normal,
                 const Real&     innerRadius,
                 const Real&     outerRadius,
                 const Real&     offset,
                 const Real&     height,
                 const Real&     angle);

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
    const char* in_file = "coarsen.inputs";
    //    registerDebugger();
//     char* in_file = argv[1];
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    //define the geometry object.
    //need the new and delete because of
    //strong construction
    //make the gometry.  this makes the first Chombo_EBIS
    Box domainFine,  domainCoar;
    RealVect dxFine,  dxCoar;
    int refRatio;
    pp.get("ref_ratio", refRatio);
    //and defines it using a geometryservice

    //set things up
    {
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
      domainCoar.setSmall(lo);
      domainCoar.setBig(hi);
      domainFine  = refine(domainCoar , 2);
    }
    {//make the geometry even finer
      Box domainFinest = refine(domainFine , refRatio);
      eekflag =  makeGeometry(domainFinest,  dxFine);
      CH_assert(eekflag == 0);
    }
    dxCoar = dxFine*2.0;

    //make grids, starting with the finest and coarsen from there
    DisjointBoxLayout gridsFine,gridsCoar;
    eekflag = makeLayoutCoar(gridsFine, domainFine);//this covers the whole domain
    CH_assert(eekflag == 0);
    coarsen(gridsCoar, gridsFine, refRatio);

    ///create ebislayout
    int nghost = 4;
    EBISLayout ebislFine, ebislCoar;
    eekflag = makeEBISL(ebislFine, gridsFine, domainFine, nghost);
    if (eekflag != 0) return eekflag;
    eekflag = makeEBISL(ebislCoar, gridsCoar, domainCoar, nghost);
    if (eekflag != 0) return eekflag;

    int nvar = 1;
    EBCellFactory ebcellfactFine(ebislFine);
    EBCellFactory ebcellfactCoar(ebislCoar);
    LevelData<EBCellFAB> errorFine(gridsFine, nvar,
                                   4*IntVect::Unit,ebcellfactFine);
    LevelData<EBCellFAB> errorCoar(gridsCoar, nvar,
                                   4*IntVect::Unit,ebcellfactCoar);

    eekflag = getError(errorFine,
                       ebislFine, gridsFine, domainFine, dxFine,
                       refRatio);
    if (eekflag != 0) return eekflag;

    eekflag = getError(errorCoar,
                       ebislCoar, gridsCoar, domainCoar, dxCoar,
                       refRatio);
    if (eekflag != 0) return eekflag;

    eekflag =  compareError(errorFine, ebislFine, gridsFine, domainFine,
                            errorCoar, ebislCoar, gridsCoar, domainCoar,
                            refRatio);
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

  }//end scoping trick
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}
/***************/
/***************/
int getError(LevelData<EBCellFAB>& a_errorCoar,
             const EBISLayout& a_ebislCoar,
             const DisjointBoxLayout& a_gridsCoar,
             const Box& a_domainCoar,
             const RealVect& a_dxCoar,
             const int& a_nref)
{
  int eekflag = 0;
  //make fine stuff
  RealVect dxFine = a_dxCoar/a_nref;
  DisjointBoxLayout gridsFine;
  Box domainFine = refine(a_domainCoar,a_nref);
  eekflag = makeLayoutFine(gridsFine, domainFine);//this does not need to cover the whole domain
  CH_assert(eekflag == 0);
  EBISLayout ebislFine;
  int nghost = 4;
  eekflag = makeEBISL(ebislFine, gridsFine, domainFine, nghost);
  if (eekflag != 0) return eekflag;

  int nvar = 1;
  EBCellFactory ebcellfactFine(ebislFine);
  EBCellFactory ebcellfactCoar(a_ebislCoar);
  LevelData<EBCellFAB> phiFine(gridsFine, nvar,
                               4*IntVect::Unit,ebcellfactFine);
  LevelData<EBCellFAB> phiCoar(a_gridsCoar, nvar,
                               4*IntVect::Unit,ebcellfactCoar);

  EBLevelDataOps::setVal(phiFine,1.e20);//bogus value
  EBLevelDataOps::setVal(phiCoar,1.e20);//bogus value
  //fill phi fine and phiCExact
  //put phiCexact into a_errorCoar
  for (DataIterator dit = gridsFine.dataIterator();
      dit.ok(); ++dit)
    {
      IntVectSet ivsBox(gridsFine.get(dit()));
      EBCellFAB& phiFineFAB = phiFine[dit()];
      for (VoFIterator vofit(ivsBox, ebislFine[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), dxFine);
          phiFineFAB(vof, 0) = rightAns;
        }
    }
  phiFine.exchange();
  for (DataIterator dit = a_gridsCoar.dataIterator();
      dit.ok(); ++dit)
    {
      IntVectSet ivsBox = IntVectSet(a_gridsCoar.get(dit()));
      EBCellFAB& errCoarFAB = a_errorCoar[dit()];
      EBCellFAB& phiCoarFAB = phiCoar[dit()];
      for (VoFIterator vofit(ivsBox, a_ebislCoar[dit()].getEBGraph());
          vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real rightAns = exactFunc(vof.gridIndex(), a_dxCoar);
          errCoarFAB(vof, 0) = rightAns;
          phiCoarFAB(vof, 0) = rightAns;
        }
    }
  //coarsen phiFine onto phiC;
  EBLevelGrid eblgFine(gridsFine,ebislFine,domainFine);
  EBLevelGrid eblgCoar(a_gridsCoar,a_ebislCoar,a_domainCoar);
  EBCoarsen ebCoarsenOp(eblgFine, eblgCoar,
                        a_nref, nvar);

  Interval zeroiv(0,0);
  ebCoarsenOp.coarsenFine(phiCoar, phiFine, zeroiv);
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
                 const Box& a_domainCoar,
                 const int& a_nref)
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
          if (Abs(finenorm) > 1.0e-12)
            {
              Real order = log(coarnorm/finenorm)/log(Real(a_nref));
              pout() << "Order of scheme = " << order << endl;
            }
        }
    }

  return eekflag;
}

/***************/
/***************/
Real exactFunc(const IntVect& a_iv, const RealVect& a_dx)
{
  Real retval;
  Real probHi;
  ParmParse pp;
  pp.get("prob_hi",probHi);
  RealVect xloc;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xloc[idir] = a_dx[idir]*(Real(a_iv[idir]) + 0.5);
    }

  Real x = xloc[0]/probHi;
  Real y = xloc[1]/probHi;
  Real z = 0.0;
  if (SpaceDim==3)
    {
      z = xloc[2]/probHi;
    }
  //retval = 1.0 + x*x + y*y + z*z + x*x*x + y*y*y + z*z*z;
  retval = cos(2*x)*cos(2*y)*cos(2*z);
  //  if (x>0.74)
  //     {
  //       retval = 1;
  //     }
  //   else
  //     {
  //       retval = 0;
  //     }
  //  retval = x*x*x;
  //retval = x*x + y*y + z*z;
  //  retval = 2.*xloc[0];
  //retval  = x + y;
  //retval  = x*x;
  //retval  = y;
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
makeLayoutFine(DisjointBoxLayout& a_dblFine,
               const Box& a_domainFine)
{
  //set up mesh refine object
  ParmParse pp;
  int eekflag= 0;
  int maxsize;
  pp.get("maxboxsize",maxsize);
  int bufferSize = 1;
  int blockFactor = 2;
  Real fillrat = 1.0;
  Box domainCoar = coarsen(a_domainFine, 2);
  Vector<int> refRat(2,2);
  BRMeshRefine meshRefObj(domainCoar, refRat, fillrat,
                          blockFactor, bufferSize, maxsize);

  Vector<Vector<Box> > oldMeshes(2);
  oldMeshes[0] = Vector<Box>(1,   domainCoar);
  oldMeshes[1] = Vector<Box>(1, a_domainFine);

  //set up coarse tags
  int nc = domainCoar.size(0);
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
#else
  Box boxf1(IntVect(0, nqu,nqu), IntVect(nmi-1,ntf-1,ntf-1));
  Box boxf2(IntVect(nmi,nte,nte), IntVect(ntf-1,nfe-1,nfe-1));
  Box boxf3(IntVect(nqu,0,0  ), IntVect(nfe-1,nqu-1,nqu-1));
  Box boxf4(IntVect(nfe,nqu,nqu), IntVect(nc -1,nte-1,nte-1));
#endif
  IntVectSet tags;
  tags |= boxf1;
  tags |= boxf2;
  tags |= boxf3;
  tags |= boxf4;

  int baseLevel = 0;
  int topLevel = 0;
  Vector<Vector<Box> > newMeshes;
  meshRefObj.regrid(newMeshes, tags, baseLevel,
                    topLevel, oldMeshes);

  const Vector<Box>& vbox = newMeshes[1];
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vbox);
  if (eekflag != 0) return eekflag;
  a_dblFine.define(vbox, procAssign);
  return eekflag;
}
/***************/
/***************/
int
makeLayoutCoar(DisjointBoxLayout& a_dblCoar,
               const Box& a_domainCoar)
{
  ParmParse pp;
  int eekflag= 0;
  int ipieces;
  ipieces = Max(ipieces, 1);
  int maxsize;
  pp.get("maxboxsize",maxsize);
  Vector<Box> vbox(1, a_domainCoar);
  domainSplit(a_domainCoar, vbox,  maxsize);
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
  a_dblCoar.define(vbox, procAssign);
  return eekflag;
}
/**********/
int makeGeometry(const Box& a_finestDomain,
                 RealVect& a_finestDx)
{
  int eekflag =  0;
  //parse input file
  ParmParse pp;
  RealVect origin = RealVect::Zero;

  Vector<Real> prob_lo(SpaceDim, 1.0);
  Real prob_hi;
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_finestDx[idir] = (prob_hi-prob_lo[idir])/a_finestDomain.size(idir);
      origin[idir] = prob_lo[idir];
    }
  int verbosity = 0;
  int whichgeom;
  pp.get("which_geom",whichgeom);
  CH_XD::EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  if (whichgeom == 0)
    {
      //allregular
      pout() << "all regular geometry" << endl;
      AllRegularService regserv;
      ebisPtr->define(a_finestDomain, origin, a_finestDx[0], regserv);
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

      GeometryShop workshop(ramp,verbosity,a_finestDx);
      //this generates the new EBIS
      ebisPtr->define(a_finestDomain, origin, a_finestDx[0], workshop);
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
      GeometryShop workshop(implicit,verbosity,a_finestDx);
      //this generates the new EBIS
      ebisPtr->define(a_finestDomain, origin, a_finestDx[0], workshop);
    }
  else if (whichgeom == 18)
    {
      pout() << "Low swirl burner geometry" << endl;
      //        AttachDebugger();

      Box domain;

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          origin[idir] = prob_lo[idir];
        }

      Real outerRadius;
      pp.get("outer_radius",outerRadius);

      Real outerThick;
      pp.get("outer_thick",outerThick);

      Real outerHeight;
      pp.get("outer_height",outerHeight);

      Real outerOffset = ((prob_hi - prob_lo[0]) - outerHeight) / 2.0 + prob_lo[0];

      Real innerRadius;
      pp.get("inner_radius",innerRadius);

      Real innerThick;
      pp.get("inner_thick",innerThick);

      Real innerOffset = 0.0;
      innerOffset += outerOffset;

      Real innerHeight;
      pp.get("inner_height",innerHeight);

      Real plateHeight;
      pp.get("plate_height",plateHeight);
      plateHeight += outerOffset;

      Real plateThick;
      pp.get("plate_thick",plateThick);

      int doHoles;
      pp.get("do_holes",doHoles);

      Real holeRadius;
      pp.get("hole_radius",holeRadius);

      Real holeSpace;
      pp.get("hole_space",holeSpace);

      int vaneNum;
      pp.get("vane_num",vaneNum);

      Real vaneThick;
      pp.get("vane_thick",vaneThick);

      RealVect vaneNorm;

      Vector<Real> vectVaneNorm;
      pp.getarr("vane_norm",vectVaneNorm,0,SpaceDim);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          vaneNorm[idir] = vectVaneNorm[idir];
        }

      Real vaneOffset;
      pp.get("vane_offset",vaneOffset);

      Real vaneHeight = innerHeight - 2*vaneOffset;

      vaneOffset += outerOffset;

      // Make the outer chamber
      BaseIF* outerChamber = makeChamber(outerRadius,outerThick,
                                         outerOffset,outerHeight);

      // Make the inner chamber
      BaseIF* innerChamber = makeChamber(innerRadius,innerThick,
                                         innerOffset,innerHeight);

      // Make the inner plate with holes
      BaseIF* holyPlate = makePlate(plateHeight,plateThick,innerRadius,
                                    doHoles,holeRadius,holeSpace);

      // Make the vanes
      BaseIF* vanes = makeVanes(vaneNum,vaneThick,vaneNorm,innerRadius,outerRadius,
                                vaneOffset,vaneHeight);

      // Union all the pieces together
      Vector<BaseIF*> pieces;

      pieces.push_back(outerChamber);
      pieces.push_back(innerChamber);
      pieces.push_back(holyPlate);
      pieces.push_back(vanes);

      UnionIF swirl(pieces);
      ComplementIF outside(swirl,true);

      GeometryShop workshop(outside,verbosity,a_finestDx);

      // This generates the new EBIS
      EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
      ebisPtr->define(a_finestDomain, origin, a_finestDx[0], workshop);

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
}
/***************/
BaseIF* makeChamber(const Real& radius,
                    const Real& thick,
                    const Real& offset,
                    const Real& height)
{
  RealVect zero(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  bool inside = true;

  Vector<BaseIF*> pieces;

  // Create a chamber
  TiltedCylinderIF chamberOut(radius + thick/2.0,xAxis,zero, inside);
  TiltedCylinderIF chamberIn (radius - thick/2.0,xAxis,zero,!inside);

  IntersectionIF infiniteChamber(chamberIn,chamberOut);

  pieces.push_back(&infiniteChamber);

  RealVect normal1(D_DECL(1.0,0.0,0.0));
  RealVect point1(D_DECL(offset,0.0,0.0));
  PlaneIF plane1(normal1,point1,inside);

  pieces.push_back(&plane1);

  RealVect normal2(D_DECL(-1.0,0.0,0.0));
  RealVect point2(D_DECL(offset+height,0.0,0.0));
  PlaneIF plane2(normal2,point2,inside);

  pieces.push_back(&plane2);

  IntersectionIF* chamber = new IntersectionIF(pieces);

  return chamber;
}

BaseIF* makePlate(const Real& height,
                  const Real& thick,
                  const Real& radius,
                  const int&  doHoles,
                  const Real& holeRadius,
                  const Real& holeSpace)
{
  RealVect zero(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  bool inside = true;

  // Create the plate without holes
  Vector<BaseIF*> pieces;

  RealVect normal1(D_DECL(1.0,0.0,0.0));
  RealVect point1(D_DECL(height,0.0,0.0));
  PlaneIF plane1(normal1,point1,inside);

  pieces.push_back(&plane1);

  RealVect normal2(D_DECL(-1.0,0.0,0.0));
  RealVect point2(D_DECL(height+thick,0.0,0.0));
  PlaneIF plane2(normal2,point2,inside);

  pieces.push_back(&plane2);

  TiltedCylinderIF middle(radius,xAxis,zero,inside);

  pieces.push_back(&middle);

  IntersectionIF plate(pieces);

  // Make the drills
  Vector<BaseIF*> drillBits;

  // Compute how many drills are needed in each direciton - 2*num+1 -
  // conservatively
  int num = (int)((radius - holeRadius) / holeSpace + 1.0);

  if (doHoles != 0)
  {
    for (int i = -num; i <= num; i++)
    {
      for (int j = -num; j <= num; j++)
      {
        RealVect center(D_DECL(0.0,i*holeSpace,j*holeSpace));
        TiltedCylinderIF* drill = new TiltedCylinderIF(holeRadius,xAxis,center,inside);

        drillBits.push_back(drill);
      }
    }
  }

  UnionIF drills(drillBits);
  ComplementIF notDrills(drills,true);

  // Drill the plate
  IntersectionIF* holyPlate = new IntersectionIF(plate,notDrills);

  return holyPlate;
}

BaseIF* makeVanes(const int&      num,
                  const Real&     thick,
                  const RealVect& normal,
                  const Real&     innerRadius,
                  const Real&     outerRadius,
                  const Real&     offset,
                  const Real&     height)
{
  RealVect zeroVect(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));

  Vector<BaseIF*> eachVane;

  for (int i = 0; i < num; i++)
  {
    Real angle = i*2.*M_PI/num;

    BaseIF* oneVane = makeVane(thick,normal,innerRadius,outerRadius,offset,height,angle);

    eachVane.push_back(oneVane);
  }

  UnionIF* allVanes = new UnionIF(eachVane);

  return allVanes;
}

BaseIF* makeVane(const Real&     thick,
                 const RealVect& normal,
                 const Real&     innerRadius,
                 const Real&     outerRadius,
                 const Real&     offset,
                 const Real&     height,
                 const Real&     angle)
{
  RealVect zeroVect(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  bool inside = true;

  Vector<BaseIF*> vaneParts;

  Real sinTheta = sin(angle);
#if CH_SPACEDIM == 3
  Real cosTheta = cos(angle);

  // Each side of the vane (infinite)
  // rotate the normal around x-axis
  RealVect normal1(D_DECL(normal[0],cosTheta*normal[1]-sinTheta*normal[2],sinTheta*normal[1]+cosTheta*normal[2]));
  // rotate point on top of vane around x-axis
  RealVect point(D_DECL(offset+height/2.0,-thick/2.0,0.0));
  RealVect point1(D_DECL(point[0],cosTheta*point[1]-sinTheta*point[2],sinTheta*point[1]+cosTheta*point[2]));
  PlaneIF plane1(normal1,point1,inside);

  vaneParts.push_back(&plane1);

  RealVect normal2(-normal1);
  // rotate point on bottom (-point[2] of vane around x-axis
  RealVect point2(D_DECL(point[0],-cosTheta*point[1]-sinTheta*point[2],-sinTheta*point[1]+cosTheta*point[2]));
  PlaneIF plane2(normal2,point2,inside);

  vaneParts.push_back(&plane2);
#endif

  // Make sure we only get something to the right of the origin
  RealVect normal3(D_DECL(0.0,-sinTheta,cosTheta));
  RealVect point3(D_DECL(0.0,0.0,0.0));
  PlaneIF plane3(normal3,point3,inside);

  vaneParts.push_back(&plane3);

  // Cut off the top and bottom
  RealVect normal4(D_DECL(1.0,0.0,0.0));
  RealVect point4(D_DECL(offset,0.0,0.0));
  PlaneIF plane4(normal4,point4,inside);

  vaneParts.push_back(&plane4);

  RealVect normal5(D_DECL(-1.0,0.0,0.0));
  RealVect point5(D_DECL(offset+height,0.0,0.0));
  PlaneIF plane5(normal5,point5,inside);

  vaneParts.push_back(&plane5);

  // The outside of the inner cylinder
  TiltedCylinderIF inner(innerRadius,xAxis,zeroVect,!inside);

  vaneParts.push_back(&inner);

  // The inside of the outer cylinder
  TiltedCylinderIF outer(outerRadius,xAxis,zeroVect,inside);

  vaneParts.push_back(&outer);

  IntersectionIF* vane = new IntersectionIF(vaneParts);

  return vane;
}
/***************/
