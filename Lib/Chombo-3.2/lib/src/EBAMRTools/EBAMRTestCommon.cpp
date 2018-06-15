#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <list>
#include <set>
using std::cout;

#include "parstream.H"

#include "EBAMRTestCommon.H"
#include "EBDebugOut.H"
#include "NamespaceHeader.H"
/***************/
int 
EBAMRTestCommon::
compareError(const LevelData<EBCellFAB>& a_errorFine,
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
      for (int icomp = 0; icomp < a_errorCoar.nComp(); icomp++)
        {
          pout() << endl << " Component = " <<  icomp << endl;
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
                                            icomp, inorm, normtype);
              Real finenorm = EBArith::norm(a_errorFine,
                                            blNormFine, a_ebislFine,
                                            icomp, inorm, normtype);
              pout() << "Coarse Error Norm = " << coarnorm << endl;
              pout() << "Fine   Error Norm = " << finenorm << endl;
              if (Abs(finenorm) > 1.0e-12)
                {
                  Real order = log(coarnorm/finenorm)/log(2.0);
                  pout() << "Order of scheme = " << order << endl;
                }
            }
        }
    }

  return eekflag;
}
/***************/
int 
EBAMRTestCommon::
makeEBISL(EBISLayout& a_ebisl,
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
EBAMRTestCommon::
makeHalfLayouts(DisjointBoxLayout& a_gridsFine,
                EBISLayout&        a_ebislFine,
                DisjointBoxLayout& a_gridsCoar,
                EBISLayout&        a_ebislCoar,
                const Box& a_domainFine,
                const Box& a_domainCoar,
                int a_refrat,
                int a_halfDir
                )
{
  //set up mesh refine object
  ParmParse pp;
  int eekflag= 0;
  int maxsize;
  pp.get("maxboxsize",maxsize);
  int bufferSize = 2;
  int blockFactor = 8;
  Real fillrat = 0.7;

  Vector<int> refRat(2,2);
  BRMeshRefine meshRefObj(a_domainCoar, refRat, fillrat,
                          blockFactor, bufferSize, maxsize);

  Vector<Vector<Box> > oldMeshes(2);
  oldMeshes[0] = Vector<Box>(1, a_domainCoar);
  oldMeshes[1] = Vector<Box>(1, a_domainFine);

  //set up coarse tags
  int nc = a_domainCoar.size(a_halfDir);
  Box halfBox = a_domainCoar;
  halfBox.growHi(a_halfDir, -nc/2);
  IntVectSet tags(halfBox);


  int baseLevel = 0;
  int topLevel = 0;
  Vector<Vector<Box> > newMeshes;
  meshRefObj.regrid(newMeshes, tags, baseLevel,
                    topLevel, oldMeshes);

  const Vector<Box>& vboxFine = newMeshes[1];
  Vector<int>  procAssign;
  eekflag = LoadBalance(procAssign,vboxFine);
  if (eekflag != 0) return eekflag;
  a_gridsFine.define(vboxFine, procAssign);


  //make the coarse grids completely cover the domain
  Vector<Box> vboxCoar;
  domainSplit(a_domainCoar, vboxCoar,  maxsize);
  eekflag = LoadBalance(procAssign,vboxCoar);
  
  a_gridsCoar.define(vboxCoar, procAssign);

  ///create ebislayout
  int nghost = 4;
  eekflag = makeEBISL(a_ebislFine, a_gridsFine, a_domainFine, nghost);
  if (eekflag != 0) return eekflag;
  eekflag = makeEBISL(a_ebislCoar, a_gridsCoar, a_domainCoar, nghost);
  if (eekflag != 0) return eekflag;
  return eekflag;
}
/**********/
int 
EBAMRTestCommon::
makeGeometry(Box& a_domain,
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

  Real prob_hi;
  pp.get("prob_hi",prob_hi);
  a_dx = RealVect::Unit*(prob_hi/n_cell[0]);
  origin  = RealVect::Zero;

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
int 
EBAMRTestCommon::
checkForZero(const LevelData<EBCellFAB>& a_errorVelo,
             const DisjointBoxLayout&    a_gridsFine,
             const EBISLayout&           a_ebislFine,
             const Box&                  a_domainFine,
             string a_funcname)
{
  int eekflag = 0;
  Real eps = 1.0e-8;
#ifdef CH_USE_FLOAT
  eps = 1.0e-3;
#endif

  for (DataIterator dit = a_gridsFine.dataIterator(); dit.ok(); ++dit)
    {
      Box grownBox = a_gridsFine.get(dit());
      grownBox.grow(1);
      grownBox &= a_domainFine;
      IntVectSet ivsBox(grownBox);
      for (VoFIterator vofit(ivsBox, a_ebislFine[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          int ihere = 0;
          if (vof.gridIndex() == EBDebugPoint::s_ivd)
            {
              ihere = 1;
            }
          for (int ivar = 0; ivar < a_errorVelo.nComp(); ivar++)
            {
              Real errorIn = a_errorVelo[dit()](vof, ivar);
              if (Abs(errorIn) > eps)
                {
                  pout() << "check for zero error too big for test " << a_funcname << endl;
                  pout() << "ivar = " << ivar <<  endl;
                  pout() << "vof = " << vof.gridIndex() << " error = " << errorIn << endl;
                  return -1;
                }
            }
        }

    }
  return eekflag;
}
/***************/
void 
EBAMRTestCommon::
outputError(const LevelData<EBCellFAB>&   a_errorFine,
            const LevelData<EBCellFAB>&   a_errorCoar,
            const DisjointBoxLayout&      a_gridsFine,
            const DisjointBoxLayout&      a_gridsCoar,
            const Box&                    a_domainFine,
            const Box&                    a_domainCoar,
            const string& a_fileFine,
            const string& a_fileCoar)
{
#ifdef CH_USE_HDF5
  int nvar = a_errorFine.nComp();
  Vector<string> names(nvar, string("fake_name"));

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
}

#include "NamespaceFooter.H"
