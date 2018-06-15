#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stdio.h>

#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DebugDump.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "PolyGeom.H"
#include "GeometryShop.H"

#include "EBFABView.H"

#include "PlaneIF.H"
#include "SphereIF.H"
#include "TiltedCylinderIF.H"
#include "IntersectionIF.H"
#include "UnionIF.H"
#include "ComplementIF.H"
#include "TransformIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

/***************/
// Define an EBIS.
/***************/
BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file = "cappedCylinder.inputs";

    if (argc > 1)
    {
      in_file = argv[1];
    }

    // Parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    Box domain;
    RealVect origin;
    Real dx;

    BaseIF* implicit;

    // Make geometry
    implicit = makeGeometry(domain,origin,dx);

    // Make grids
    DisjointBoxLayout grids;
    makeLayout(grids,domain);

    // Create ebislayout
    int nghost = 0;
    EBISLayout ebisl;
    makeEBISL(ebisl,grids,domain,nghost);

    // Make a LevelData
    int nComps = 1;

    IntVect ghost = IntVect::Unit;
    ghost *= nghost;

    RefCountedPtr<DataFactory<EBCellFAB> > rcpFactory(new EBCellFactory(ebisl));
    LevelData<EBCellFAB> level(grids,nComps,ghost,*rcpFactory);

    // Put some data in the data holders
    fillData(level,origin,dx,*implicit);

    // Done with this object
    delete implicit;

    // Write the data and the EB out
    const char* basename = "cappedCylinder";

    char name[1000];
    sprintf(name,"%s%dd.hdf5",basename,SpaceDim);

#ifdef CH_USE_HDF5
    writeEBLevelname(&level,name);
#endif
  } // End scoping trick

  // Clean up
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();

  CH_TIMER_REPORT();

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}

BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx)
{
  RealVect direction;
  RealVect center;
  Real     radius;
  Real     length;
  int      capping;
  bool     insideRegular;

  // parse input file
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
      pout() << "Bogus number of cells input = " << n_cell[ivec];
      exit(1);
    }

    hi[ivec] = n_cell[ivec] - 1;
  }

  a_domain.setSmall(lo);
  a_domain.setBig(hi);

  Vector<Real> prob_lo(SpaceDim,1.0);
  Real prob_hi;

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);

  a_dx = (prob_hi-prob_lo[0])/n_cell[0];

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_origin[idir] = prob_lo[idir];
  }

  // ParmParse doesn't get RealVects, so work-around with Vector<Real>
  Vector<Real> vectorDirection;
  pp.getarr("direction",vectorDirection,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    direction[idir] = vectorDirection[idir];
  }

  Vector<Real> vectorCenter;
  pp.getarr("center",vectorCenter,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    center[idir] = vectorCenter[idir];
  }

  pp.get("radius",radius);
  pp.get("length",length);

  // capping == 0 -> flat caps
  // capping != 0 -> spherical caps
  pp.get("capping",capping);

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  // Start with an infinite cylinder whose axis is the x-axis
  RealVect zero(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  bool inside = true;
  TiltedCylinderIF cylinder(radius,xAxis,zero,inside);

  // Create to planes to make it finite (via intersections)
  RealVect normal1(D_DECL(1.0,0.0,0.0));
  RealVect point1(D_DECL(-length/2.0,0.0,0.0));
  PlaneIF plane1(normal1,point1,inside);

  RealVect normal2(D_DECL(-1.0,0.0,0.0));
  RealVect point2(D_DECL(length/2.0,0.0,0.0));
  PlaneIF plane2(normal2,point2,inside);

  IntersectionIF interval(plane1,plane2);
  IntersectionIF finiteCylinder(cylinder,interval);

  TransformIF* transform;

  // Add spherical caps if needed
  if (capping == 0)
  {
    transform = new TransformIF(finiteCylinder);
  }
  else
  {
    SphereIF sphere1(radius,point1,inside);
    SphereIF sphere2(radius,point2,inside);

    UnionIF spheres(sphere1,sphere2);
    UnionIF cappedCylinder(finiteCylinder,spheres);

    transform = new TransformIF(cappedCylinder);
  }

  // Rotate the x-axis to the cylinder axis
  transform->rotate(xAxis,direction);

  // Translate to the center
  transform->translate(center);

  // Complement if necessary
  ComplementIF insideOut(*transform,!insideRegular);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(insideOut,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return insideOut.newImplicitFunction();
}
