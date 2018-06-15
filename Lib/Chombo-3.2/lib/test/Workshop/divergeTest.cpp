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
#include "VoFIterator.H"
#include "AllRegularService.H"
#include "GeometryShop.H"
#include "PlaneIF.H"
#include "SphereIF.H"
#include "RhodoneaIF.H"
#include "UnionIF.H"
#include "ComplementIF.H"
#include "TransformIF.H"
#include "TiltedCylinderIF.H"
#include "IntersectionIF.H"
#include "EBArith.H"
#include "PolyGeom.H"
#include "EBDebugDump.H"
#include "DebugDump.H"
#include "UsingNamespace.H"

typedef Real (*FluxFunc)(const RealVect& location);

Real constFunc(const RealVect& location)
{
  return 1.0;
}

Real linearFunc(const RealVect& location)
{
  Real retval = location[0];
//  for (int idir = 0; idir < SpaceDim; idir++)
//    retval += location[idir];

  return retval;
}




int makeLayout(DisjointBoxLayout& a_dbl,
               const Box& a_domain)
{
  ParmParse pp;
  int eekflag = 0;
  int ipieces;
  ipieces = Max(ipieces, 1);
  int maxsize;
  Vector<Box> vbox(1, a_domain);

  pp.get("maxboxsize",maxsize);

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

int makeGeometry(Box& a_domain)
{
  Real dx;
  RealVect origin;
  int eekflag =  0;
  //parse input file
  ParmParse pp;

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
  dx = (prob_hi-prob_lo[0])/n_cell[0];
  RealVect dxVect = dx*RealVect::Unit;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = prob_lo[idir];
    }
  int verbosity = 0;
  int whichgeom;
  pp.get("which_geom",whichgeom);
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  if (whichgeom == 0)
    {
      //allregular
      pout() << "all regular geometry" << endl;
      AllRegularService regserv;
      ebisPtr->define(a_domain, origin, dx, regserv);
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

      GeometryShop workshop(ramp,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
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
      GeometryShop workshop(implicit,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
    }
  else if (whichgeom == 13)
    {
      pout() << "rhodonea geometry" << endl;
      vector<Real> tmp(SpaceDim);
      pp.getarr("rhodonea_center", tmp, 0, SpaceDim);
      RealVect rhodoneaCenter;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          rhodoneaCenter[idir] = tmp[idir];
        }
      Real innerRadius;
      pp.get("inner_radius", innerRadius);

      Real outerRadius;
      pp.get("outer_radius", outerRadius);

      int frequency;
      pp.get("frequency", frequency);

      bool     insideRegular = false;

      RhodoneaIF implicit(innerRadius, outerRadius, frequency,
                          rhodoneaCenter, insideRegular);

      GeometryShop workshop(implicit,verbosity,dxVect);
      //this generates the new EBIS
      ebisPtr->define(a_domain, origin, dx, workshop);
    }
  else
    {
      //bogus which_geom
      pout() << " bogus which_geom input = " << whichgeom;
      eekflag = 33;
    }

  return eekflag;
}

Real
divergence(const FluxFunc  a_func,
           const EBISBox&  a_ebisBox,
           const VolIndex& a_vof,
           const Real&     a_dx)
{
  RealVect problo = RealVect::Zero;
  Real retval = 0;
  RealVect boundOff  = a_ebisBox.bndryCentroid(a_vof);
  boundOff *= a_dx;

  RealVect centLoc = EBArith::getVofLocation(a_vof, a_dx*RealVect::Unit, problo);
  centLoc += boundOff;

  Real bndryArea = a_ebisBox.bndryArea(a_vof);
  Real bndryFlux = a_func(centLoc);
  RealVect normal = a_ebisBox.normal(a_vof);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real bndryContrib = -bndryFlux*bndryArea*normal[idir];
      Real openContrib = 0;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          Real rsign = isign;
          Vector<FaceIndex> faces = a_ebisBox.getFaces(a_vof, idir, sit());
          for (int iface = 0; iface < faces.size(); ++iface)
            {
              RealVect faceLoc = EBArith::getFaceLocation(faces[iface], a_dx*RealVect::Unit, problo);
              RealVect faceOff = a_ebisBox.centroid(faces[iface]);
              faceOff[idir] = 0.0;
              faceOff *= a_dx;
              faceLoc += faceOff;

              Real areaFrac = a_ebisBox.areaFrac(faces[iface]);
              Real faceFlux = a_func(faceLoc);
              openContrib += rsign*faceFlux*areaFrac;
            }
        }
      retval += openContrib + bndryContrib;

    }
  retval /= a_dx;
  return retval;
}
/***************/
int checkLevel(const EBISLayout& a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box& a_domain,
               const Real& a_dx)
{
  int eekflag =  0;
#ifdef CH_USE_FLOAT
  Real tolerance = 1.0e-2;
#else
  Real tolerance = 1.0e-8;
#endif
  for (DataIterator dit  = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      IntVectSet ivs = a_ebisl[dit()].getIrregIVS(grid);
      for (VoFIterator vofit(ivs, a_ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          {
            Real constDivergence = divergence(constFunc, a_ebisl[dit()], vofit(), a_dx);
            if (Abs(constDivergence) > tolerance)
              {
                pout() << "failed divergence of constant test at vof " << vofit().gridIndex() << endl;
                return -1;
              }
          }

          {
            Real kappa = a_ebisl[dit()].volFrac(vofit());
            Real linearDivergence = divergence(linearFunc, a_ebisl[dit()], vofit(), a_dx);
            Real rightAns = kappa;
            if (Abs(linearDivergence - rightAns) > tolerance)
              {
                pout() << "failed divergence of linear function test at vof " << vofit().gridIndex() << endl;
                return -2;
              }
          }
        }
    }

  return eekflag;
}
int checkAllLevels(const Box& a_domain)
{
  int eekflag = 0;
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());
  Box levelDomain = a_domain;
  Real levelDx = 1.0/a_domain.size(0);
  for (int ilev = 0; ilev < ebisPtr->numLevels(); ilev++)
    {
      CH_assert(!levelDomain.isEmpty());
      Vector<Box> vbox(1,levelDomain);
      Vector<int> proc(1, 0);
      DisjointBoxLayout levelDBL(vbox, proc);
      EBISLayout levelEBISL;
      ebisPtr->fillEBISLayout(levelEBISL, levelDBL, levelDomain, 1);
      eekflag = checkLevel(   levelEBISL, levelDBL, levelDomain, levelDx);
      if (eekflag != 0)
        {
          pout() << "problem at domain = " << levelDomain << ", ilev = " << ilev << endl;
          return eekflag;
        }
      levelDomain.coarsen(2);
      levelDx *= 2;
    }
  return 0;
}
int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  // begin forever present scoping trick
  {
    const char* in_file;
    if (argc < 2 || argv[1][0] == '-')
      {
        in_file = "ramp.inputs";
      }
    else
      {
        in_file = argv[1];
      }
    //parse input file
    ParmParse pp(0,NULL,NULL,in_file);
    Box domain;
    int eekflag = 0;
    eekflag =  makeGeometry(domain);
    if (eekflag != 0)
      {
        pout() << "non zero eek detected = " << eekflag << endl;
        MayDay::Error("problem in makeGeometry");
      }

    eekflag = checkAllLevels(domain);
    if (eekflag != 0)
      {
        pout() << "checkAllLevels: eek = " << eekflag << endl;
        MayDay::Error("problem in checkEBISL");
      }
    pout() << "diverge test passed" << endl;
  } // end scoping trick

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}



