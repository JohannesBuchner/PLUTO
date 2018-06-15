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
#include "TorusIF.H"
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
    const char* in_file = "borromeanRings.inputs";

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
    const char* basename = "borromeanRings";

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
  RealVect center;
  Real     endRadius;
  Real     linkRadius;
  Real     length;
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
  Vector<Real> vectorCenter;
  pp.getarr("center",vectorCenter,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    center[idir] = vectorCenter[idir];
  }

  pp.get("endRadius",endRadius);
  pp.get("linkRadius",linkRadius);
  pp.get("length",length);

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  // Define the three axes (in 3D)
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  RealVect yAxis(D_DECL(0.0,1.0,0.0));
  RealVect zAxis(D_DECL(0.0,0.0,1.0));

  // Define the origin
  RealVect zero(D_DECL(0.0,0.0,0.0));

  // Stay to the inside until the end
  bool inside = true;

  // Start with an infinite cylinder whose axis is the x-axis
  TiltedCylinderIF cylinder(linkRadius,xAxis,zero,inside);

  // Create to planes to make it finite (via intersections)
  RealVect normal1(D_DECL(1.0,0.0,0.0));
  RealVect point1(D_DECL(-length/2.0,0.0,0.0));
  PlaneIF plane1(normal1,point1,inside);

  RealVect normal2(D_DECL(-1.0,0.0,0.0));
  RealVect point2(D_DECL(length/2.0,0.0,0.0));
  PlaneIF plane2(normal2,point2,inside);

  // Cut the infinite cylinder into a finite cylinder
  IntersectionIF interval(plane1,plane2);
  IntersectionIF finiteCylinder(cylinder,interval);

  // Add spherical caps
  SphereIF sphere1(linkRadius,point1,inside);
  SphereIF sphere2(linkRadius,point2,inside);

  UnionIF twoSpheres1(sphere1,sphere2);
  UnionIF cappedCylinder(finiteCylinder,twoSpheres1);

  // Create two spherically capped cylinders and separate them in the y
  // direction
  TransformIF straight1(cappedCylinder);
  RealVect trans1(D_DECL(0.0, endRadius,0.0));
  straight1.translate(trans1);

  TransformIF straight2(cappedCylinder);
  RealVect trans2(D_DECL(0.0,-endRadius,0.0));
  straight2.translate(trans2);

  // Start with a torus for the end of the links
  TorusIF torus(endRadius,linkRadius,point1,inside);

  // Cut half of it away
  PlaneIF plane3(normal2,point1,inside);
  IntersectionIF halfTorus(torus,plane3);

  // Add spherically caps
  SphereIF sphere3(linkRadius,point1,inside);

  TransformIF transSphere1(sphere3);
  transSphere1.translate(trans1);

  TransformIF transSphere2(sphere3);
  transSphere2.translate(trans2);

  UnionIF twoSpheres2(transSphere1,transSphere2);
  UnionIF cappedTorus1(halfTorus,twoSpheres2);

  // Make a rotated copy
  TransformIF cappedTorus2(cappedTorus1);
  cappedTorus2.rotate(M_PI,zero,zAxis);

  // Put the link along the x axis together
  UnionIF straights(straight1,straight2);
  UnionIF cappedTori(cappedTorus1,cappedTorus2);
  UnionIF linkX(straights,cappedTori);

  // Create and orient the link along the y axis
  TransformIF linkY(linkX);
  linkY.rotate(M_PI/2.0,zero,zAxis);
  linkY.rotate(M_PI/2.0,zero,yAxis);

  // Create and orient the link along the y axis
  TransformIF linkZ(linkY);
  linkZ.rotate(M_PI/2.0,zero,xAxis);
  linkZ.rotate(M_PI/2.0,zero,zAxis);

  // Put everything together
  UnionIF twoRings(linkX,linkY);
  UnionIF borromeanRings(twoRings,linkZ);

  TransformIF finalCenter(borromeanRings);
  finalCenter.translate(center);

  // Complement if necessary
  ComplementIF insideOut(finalCenter,!insideRegular);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(insideOut,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return insideOut.newImplicitFunction();
}
