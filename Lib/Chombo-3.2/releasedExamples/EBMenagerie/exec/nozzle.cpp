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

#include "Vector.H"
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
#include "ComplementIF.H"
#include "IntersectionIF.H"
#include "UnionIF.H"
#include "LatheIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

/***************/
// Define an EBIS.
/***************/
BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx);

UnionIF* makeCrossSection(const Vector<Vector<RealVect> >& a_polygons);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file;

    if (SpaceDim == 2)
    {
      in_file = "nozzle2d.inputs";
    }
    else if (SpaceDim == 3)
    {
      in_file = "nozzle3d.inputs";
    }
    else
    {
      MayDay::Error("SpaceDim must be 2 or 3");
    }

    if (argc > 1 && strcmp(argv[1],"nozzle.inputs") != 0)
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

    createEBDistributionFiles();

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
    const char* basename = "nozzle";

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
  bool insideRegular;

  BaseIF* retval;

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

  Vector<Real> prob_lo(SpaceDim,0.0);
  Vector<Real> prob_hi(SpaceDim,1.0);

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.getarr("prob_hi",prob_hi,0,SpaceDim);

  a_dx = (prob_hi[0]-prob_lo[0])/n_cell[0];

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_origin[idir] = prob_lo[idir];
  }

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  // Data for polygons making up nozzle
  Vector<Vector<RealVect> > polygons;

  // Nothing initially
  polygons.resize(0);

  // For building each polygon
  Vector<RealVect> polygon;
  RealVect point(RealVect::Zero);

  // This is how far the polygons extend - the physical domain should be small
  // than this (20 meters)
  Real far = 20000.0;

  // Add body (two pieces) which union to one

  // Piece 1 - add the points and then save the list as a polygon
  polygon.resize(0);

  point[0] =       far; point[1] =  5.64012000;
  polygon.push_back(point);

  point[0] = 1.2499200; point[1] =  5.64012000;
  polygon.push_back(point);

  point[0] = 1.2499200; point[1] = -0.95678300;
  polygon.push_back(point);

  point[0] =       far; point[1] = -0.95678300;
  polygon.push_back(point);

  polygons.push_back(polygon);

  // Piece 2 - add the points and then save the list as a polygon
  polygon.resize(0);

  point[0] =  1.2499200; point[1] =  5.64012000;
  polygon.push_back(point);

  point[0] =  1.1058800; point[1] =  5.64012000;
  polygon.push_back(point);

  point[0] =  0.3966470; point[1] =  1.62330000;
  polygon.push_back(point);

  point[0] =  0.3964730; point[1] =  0.42332500;
  polygon.push_back(point);

  point[0] =  0.7493570; point[1] =  0.00176610;
  polygon.push_back(point);

  point[0] =  1.2499200; point[1] = -0.00777410;
  polygon.push_back(point);

  polygons.push_back(polygon);

  // Add poppet - add the points and then save the list as a polygon
  polygon.resize(0);

  point[0] =  0.8528400; point[1] = -far;
  polygon.push_back(point);

  point[0] =  0.8528400; point[1] = -0.92157200;
  polygon.push_back(point);

  point[0] =  0.8521970; point[1] = -0.32181000;
  polygon.push_back(point);

  point[0] =  0.0000000; point[1] =  0.69393192;
  polygon.push_back(point);

  point[0] =  0.0000000; point[1] = -far;
  polygon.push_back(point);

  polygons.push_back(polygon);

  // Make the vector of (convex) polygons (vectors of points) into a union
  // of convex polygons, each made from the intersection of a set of half
  // planes/spaces - all represented by implicit functions.
  UnionIF* crossSection = makeCrossSection(polygons);

  if (SpaceDim == 2)
  {
    // In 2D use "as is"

    // Complement if necessary
    ComplementIF insideOut(*crossSection,!insideRegular);
    retval = insideOut.newImplicitFunction();
  }
  else
  {
    // In 3D rotate about the z-axis and complement if necessary
    LatheIF implicit(*crossSection,insideRegular);
    retval = implicit.newImplicitFunction();
  }

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(*retval,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return retval;
}

UnionIF* makeCrossSection(const Vector<Vector<RealVect> >& a_polygons)
{
  // The final result
  UnionIF* retval;

  // Get the number of polygons and make this inside of the domain
  // the inside of the polygons
  int numPolys = a_polygons.size();
  bool inside = true;

  // A list of all the polygons as implicit functions
  Vector<BaseIF*> polytopes;
  polytopes.resize(0);

  // Process each polygon
  for (int p = 0; p < numPolys; p++)
  {
    // All the half planes/spaces used to make a polygon
    Vector<BaseIF*> planes;
    planes.resize(0);

    // Get the current polygon (as a vector of points)
    const Vector<RealVect>& polygon = a_polygons[p];

    // Get the number of points in the polygon
    int numPts = polygon.size();

    // Process each pair of points
    for (int n = 0; n < numPts; n++)
    {
      // The normal and point is space used to specify each half plane/space
      RealVect normal(RealVect::Zero);
      RealVect point;

      // Set the normal remembering that the last point connects to the first
      // point.
      normal[0] = -(polygon[(n+1) % numPts][1] - polygon[n][1]);
      normal[1] =  (polygon[(n+1) % numPts][0] - polygon[n][0]);

      point = polygon[n];

      // Generate the appropriate half plane/space (as an implicit function)
      PlaneIF* plane;
      plane = new PlaneIF(normal,point,inside);

      // Save the result
      planes.push_back(plane);
    }

    // Intersect all the half planes/spaces to create an implicit function
    // that represents the polygon
    IntersectionIF* polygonIF = new IntersectionIF(planes);

    polytopes.push_back(polygonIF);
  }

  // Union all the polygon implicit functions to get the implicit function
  // returned
  retval = new UnionIF(polytopes);

  return retval;
}
