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
#include "SphereIF.H"
#include "PolynomialIF.H"
#include "ComplementIF.H"
#include "IntersectionIF.H"
#include "UnionIF.H"
#include "TransformIF.H"
#include "LatheIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

/***************/
// Define an EBIS.
/***************/
BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx);

BaseIF* makeCrossSection();

BaseIF* makeOuterChamber();

BaseIF* makeInlet();

BaseIF* makeRectangle(const Real & a_xSize,
                      const Real & a_ySize,
                      const Real & a_yOffset);

BaseIF* makeTrapezoid(const Real & a_xSizeBottom,
                      const Real & a_xSizeTop,
                      const Real & a_ySize,
                      const Real & a_yOffset);

BaseIF* makeConvexPolygon(const Vector<RealVect> & a_corners);

BaseIF* makeCurve(const Real & a_xSizeBottom,
                  const Real & a_xSizeTop,
                  const Real & a_ySize,
                  const Real & a_yOffset);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file = "turbine.inputs";

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
    const char* basename = "turbine";

    char name[1000];
    sprintf(name,"%s%dd.hdf5",basename,SpaceDim);

#ifdef CH_USE_HDF5
    writeEBLevelname(&level,name);
#endif
  } // End scoping trick

  // Clean up
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->clear();

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif

  return 0;
}

BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     Real&     a_dx)
{
  bool insideRegular = true;

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

  // Generate the cross-section of the turbine (2D) which can be rotated about
  // the y-axis to produce the 3D geometry
  BaseIF* crossSection = makeCrossSection();

  if (SpaceDim == 2)
  {
    // In 2D, use the cross-section
    retval = crossSection;
  }
  else
  if (SpaceDim == 3)
  {
    // In 3D, rotate the cross-section about the y-axis
    LatheIF rotated(*crossSection,insideRegular);
    TransformIF implicit(rotated);

    // Nominally, "LatheIF" makes the axis of rotation end up being the z-axis
    // but to make things consistent with the 2D case we rotate the z-axis
    // to the y-axis
    RealVect zAxis = BASISREALV(2);
    RealVect yAxis = BASISREALV(1);

    implicit.rotate(zAxis,yAxis);

    retval = implicit.newImplicitFunction();
  }
  else
  {
    MayDay::Error("Turbine geometry only defined in 2D and 3D");
  }

  // Generate the EBIndexSpace
  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(*retval,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return retval;
}

BaseIF* makeCrossSection()
{
  // Generate the geometry for the chamber and the inlet where the fluid is
  // inside the chamber and outside the inlet
  BaseIF* outer = makeOuterChamber();
  BaseIF* inlet = makeInlet();

  // Intersect the fluid inside the chamber with the fluid outside the inlet
  IntersectionIF* retval = new IntersectionIF(*outer,*inlet);

  return retval;
}

BaseIF* makeOuterChamber()
{
  // Below the cross-section (2D) is made by combining rectangles, trapezoids,
  // and a shape curved on the left and right.  When these are rotated about
  // the y-axis in 3D (see above), the rectangles become cylinders, the
  // trapezoids become truncated cones, and the curved shape becomes a curved
  // solid of rotation.

  // A small number used to extend the overall chamber slightly outside the
  // computational domain
  Real epsilon = 1.0e-10;

  // These parameters define the inlet chamber
  Real inletChamberDiameter1 =  50.1;
  Real inletChamberDiameter2 =  75.9;

  Real inletChamberOffset    =   0.0;
  Real inletChamberLength1   =  16.9;
  Real inletChamberLength2   = 104.0;
  Real inletChamberLength3   =  16.9;

  Real inletChamberLength    = inletChamberLength1 +
                               inletChamberLength2 +
                               inletChamberLength3;

  // The start of the inlet chamber is a rectangle in 2D
  BaseIF* inletChamber1 = makeRectangle(inletChamberDiameter1,
                                        inletChamberLength1 + epsilon,
                                        inletChamberOffset  - epsilon);

  // Then it expands making a trapezoid in 2D
  BaseIF* inletChamber2 = makeTrapezoid(inletChamberDiameter1,
                                        inletChamberDiameter2,
                                        inletChamberLength2,
                                        inletChamberOffset  +
                                        inletChamberLength1);

  // And finally continues as a rectangle in 2D
  BaseIF* inletChamber3 = makeRectangle(inletChamberDiameter2,
                                        inletChamberLength3,
                                        inletChamberOffset  +
                                        inletChamberLength1 +
                                        inletChamberLength2);

  // Union these together to make the inlet chamber
  UnionIF inletChamber12(*inletChamber1,*inletChamber2);
  UnionIF inletChamber  (inletChamber12,*inletChamber3);

  // These parameters define the main chamber
  Real mainChamberDiameter = 101.6;

  Real mainChamberOffset   = inletChamberOffset +
                             inletChamberLength;
  Real mainChamberLength   = 420.5;

  // The entire main chamber is a rectangle in 2D
  BaseIF* mainChamber = makeRectangle(mainChamberDiameter,
                                      mainChamberLength,
                                      mainChamberOffset);

  // These parameters define the exit chamber (see "makeCurve" for the
  // details)
  Real exitChamberDiameter1 = mainChamberDiameter;
  Real exitChamberDiameter2 = 50.1;

  Real exitChamberOffset = mainChamberOffset +
                           mainChamberLength;
  Real exitChamberLength = 153.6;

  // The exit chamber is curved on the left and right in 2D
  BaseIF* exitChamber = makeCurve(exitChamberDiameter1,
                                  exitChamberDiameter2,
                                  exitChamberLength + epsilon,
                                  exitChamberOffset);

  // Combine all the chambers together to form one chamber
  UnionIF startChambers(inletChamber,*mainChamber);
  UnionIF chamber(startChambers,*exitChamber);

  return chamber.newImplicitFunction();
}

BaseIF* makeInlet()
{
  // Below the cross-section (2D) is made by combining rectangles and
  // trapezoids.  When these are rotated about the y-axis in 3D (see above),
  // the rectangles become cylinders and the trapezoids become truncated
  // cones.

  // A small number used to extend the overall chamber slightly outside the
  // computational domain
  Real epsilon = 1.0e-10;

  // These parameters define the inlet
  Real inletDiameter1 =   8.4;
  Real inletDiameter2 =  27.6;
  Real inletDiameter3 =  63.6;

  Real inletOffset    =   0.0;
  Real inletLength1   =  67.8;
  Real inletLength2   =  18.6;
  Real inletLength3   =  48.8;
  Real inletLength4   =   4.7;

  // The start of the inlet is a rectangle in 2D
  BaseIF* inletSection1 = makeRectangle(inletDiameter1,
                                        inletLength1 + epsilon,
                                        inletOffset  - epsilon);

  // Which expands forming a trapezoid in 2D
  BaseIF* inletSection2 = makeTrapezoid(inletDiameter1,
                                        inletDiameter2,
                                        inletLength2,
                                        inletOffset  +
                                        inletLength1);

  // Then there is a continuation rectangle in 2D
  BaseIF* inletSection3 = makeRectangle(inletDiameter2,
                                        inletLength3,
                                        inletOffset  +
                                        inletLength1 +
                                        inletLength2);

  // And, finally, a larger rectangle a the top in 2D
  BaseIF* inletSection4 = makeRectangle(inletDiameter3,
                                        inletLength4,
                                        inletOffset  +
                                        inletLength1 +
                                        inletLength2 +
                                        inletLength3);

  // Put together all the pieces of the inlet
  UnionIF inletSection12(*inletSection1,*inletSection2);
  UnionIF inletSection34(*inletSection3,*inletSection4);
  UnionIF inlet(inletSection12,inletSection34);

  // Define the fluid to be outside the inlet
  return new ComplementIF(inlet,true);
}

BaseIF* makeRectangle(const Real & a_xSize,
                      const Real & a_ySize,
                      const Real & a_yOffset)
{
  // A rectangle is simply a trapezoid where the bottom and top are equal in
  // size
  return makeTrapezoid(a_xSize,
                       a_xSize,
                       a_ySize,
                       a_yOffset);
}

BaseIF* makeTrapezoid(const Real & a_xSizeBottom,
                      const Real & a_xSizeTop,
                      const Real & a_ySize,
                      const Real & a_yOffset)
{
  // Define the trapezoid via its corners as a polygon
  Vector<RealVect> corners;

  // Traverse the corners counter-clockwise to guarantee the fluid is inside
  // the polygon
  RealVect p1(D_DECL(-a_xSizeBottom/2.0,a_yOffset          ,0.0));
  RealVect p2(D_DECL( a_xSizeBottom/2.0,a_yOffset          ,0.0));
  RealVect p3(D_DECL( a_xSizeTop   /2.0,a_yOffset + a_ySize,0.0));
  RealVect p4(D_DECL(-a_xSizeTop   /2.0,a_yOffset + a_ySize,0.0));

  corners.push_back(p1);
  corners.push_back(p2);
  corners.push_back(p3);
  corners.push_back(p4);

  // Make the polygon
  return makeConvexPolygon(corners);
}

BaseIF* makeConvexPolygon(const Vector<RealVect> & a_corners)
{
  // All the half planes/spaces used to make a polygon
  Vector<BaseIF*> planes;
  planes.resize(0);

  // Get the number of corners in the polygon
  int numPts = a_corners.size();

  // This inside is in the direction of the normal
  bool inside = true;

  // Process each pair of points
  for (int n = 0; n < numPts; n++)
  {
    // The normal and point is space used to specify each half plane/space
    RealVect normal(RealVect::Zero);
    RealVect point;

    // Set the normal remembering that the last point connects to the first
    // point.
    normal[0] = -(a_corners[(n+1) % numPts][1] - a_corners[n][1]);
    normal[1] =  (a_corners[(n+1) % numPts][0] - a_corners[n][0]);

    point = a_corners[n];

    // Generate the appropriate half plane/space (as an implicit function)
    PlaneIF* plane;
    plane = new PlaneIF(normal,point,inside);

    // Save the result
    planes.push_back(plane);
  }

  // Intersect all the half planes/spaces to create an implicit function
  // that represents the polygon
  IntersectionIF* polygonIF = new IntersectionIF(planes);

  // Return the implicit function representing the polygon
  return polygonIF;
}

BaseIF* makeCurve(const Real & a_diameterBottom,
                  const Real & a_diameterTop,
                  const Real & a_length,
                  const Real & a_offset)
{
  // Define a curve to be used as the left and right side of the exit chamber.
  // I start with the 2D, cubic equation:  x = 2y^3 - 3y^2 + 1 which goes
  // through (1,0) and (0,1) and has an extremum at each of those points.  I
  // then scale and tranlate it to match the diameters and lengths given.

  // Define the initial equation
  Vector<PolyTerm> curveTerms;

  PolyTerm curTerm;

  curTerm.coef = 1.0;
  curTerm.powers = IntVect(D_DECL(1,0,0));
  curveTerms.push_back(curTerm);

  curTerm.coef = -2.0;
  curTerm.powers = IntVect(D_DECL(0,3,0));
  curveTerms.push_back(curTerm);

  curTerm.coef = 3.0;
  curTerm.powers = IntVect(D_DECL(0,2,0));
  curveTerms.push_back(curTerm);

  curTerm.coef = -1.0;
  curTerm.powers = IntVect(D_DECL(0,0,0));
  curveTerms.push_back(curTerm);

  // Make this into an implicit function
  PolynomialIF curvePolyUnit(curveTerms,true);

  // Scale and translate the function to be in the needed location
  TransformIF side1(curvePolyUnit);

  side1.scale(RealVect(D_DECL(a_diameterBottom/2.0 - a_diameterTop/2.0,
                              a_length,
                              1.0)));
  side1.translate(RealVect(D_DECL(a_diameterTop/2.0,
                                  a_offset,
                                  0.0)));

  // One side is just the other side reflected in x
  TransformIF side2(side1);
  side2.scale(RealVect(D_DECL(-1.0,1.0,1.0)));

  // Combine the the sides
  IntersectionIF sides(side1,side2);

  // Create the bottom and top
  RealVect normalBottom = BASISREALV(1);
  RealVect pointBottom(D_DECL(0.0,a_offset,0.0));

  PlaneIF bottom(normalBottom,pointBottom,true);

  RealVect normalTop = -BASISREALV(1);
  RealVect pointTop(D_DECL(0.0,a_offset + a_length,0.0));

  PlaneIF top(normalTop,pointTop,true);

  // Combine the sides, bottom, and top
  IntersectionIF bottomTop(bottom,top);
  IntersectionIF curve(sides,bottomTop);

  return curve.newImplicitFunction();
}
