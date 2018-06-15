#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <cmath>

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

#include "SphereIF.H"
#include "TiltedCylinderIF.H"
#include "UnionIF.H"
#include "IntersectionIF.H"
#include "ComplementIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

/***************/
// Define an EBIS.
/***************/
BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx);

/***************/
// Get the distance from a point to a line
/***************/
Real getDistance(const RealVect& a_point,
                 const RealVect& a_lineDirection,
                 const RealVect& a_linePoint);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file = "channelWithSpheres.inputs";

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
    const char* basename = "channelWithSpheres";

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
  RealVect cylDirection;
  RealVect cylPoint;
  Real     cylRadius;

  int      sphNumber;
  Real     sphMinRadius;
  Real     sphMaxRadius;

  bool     insideRegular;

  // Initialize random number generator
  srand48(time(0));

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
  pp.getarr("cylDirection",vectorDirection,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    cylDirection[idir] = vectorDirection[idir];
  }

  Vector<Real> vectorPoint;
  pp.getarr("cylPoint",vectorPoint,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    cylPoint[idir] = vectorPoint[idir];
  }

  pp.get("cylRadius",cylRadius);

  pp.get("sphNumber",sphNumber);

  pp.get("sphMinRadius",sphMinRadius);
  pp.get("sphMaxRadius",sphMaxRadius);

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  // Stay inside until the end
  bool inside = true;

  // The vector of spheres
  Vector<BaseIF*> spheres;
  spheres.resize(sphNumber);

  // Place random spheres inside the cylinder
  int i = 0;

  while (i < sphNumber)
  {
    RealVect center;
    Real     radius;

    // Get a random center and radius for a sphere
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      center[idir] = (prob_hi - prob_lo[idir])*drand48() + prob_lo[idir];
    }

    radius = (sphMaxRadius - sphMinRadius)*drand48() + sphMinRadius;

    // Get the distance from the center of the sphere to the axis of the
    // cylinder
    Real dist = getDistance(center,cylDirection,cylPoint);

    // If the sphere is inside cylinder, use it
    if (dist + radius <= 0.95 * cylRadius)
    {
      spheres[i] = new SphereIF(radius,center,inside);
      i++;
    }
  }

  // Take the union of the insides of the spheres
  UnionIF insideSpheres(spheres);

  // Complement to get the outside of the spheres
  ComplementIF outsideSpheres(insideSpheres,true);

  // Define the cylinder
  TiltedCylinderIF cylinder(cylRadius,cylDirection,cylPoint,inside);

  // Intersect the inside of the cylinder with the outside of the spheres
  IntersectionIF implicit(cylinder,outsideSpheres);

  // Complement if necessary
  ComplementIF insideOut(implicit,!insideRegular);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(insideOut,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return insideOut.newImplicitFunction();
}

Real getDistance(const RealVect& a_point,
                 const RealVect& a_lineDirection,
                 const RealVect& a_linePoint)
{
  RealVect diff = a_linePoint;
  diff -= a_point;

  Real normDiff2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    normDiff2 += diff[idir] * diff[idir];
  }

  Real normDirection2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    normDirection2 += a_lineDirection[idir] * a_lineDirection[idir];
  }

  Real dotDiffDirection2 = 0.0;
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    dotDiffDirection2 += diff[idir] * a_lineDirection[idir];
  }
  dotDiffDirection2 = dotDiffDirection2 * dotDiffDirection2;

  return sqrt(normDiff2 - dotDiffDirection2 / normDirection2);
}
