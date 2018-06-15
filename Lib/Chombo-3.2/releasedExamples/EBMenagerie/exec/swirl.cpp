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
                 const Real&     height);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file = "swirl.inputs";

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
    const char* basename = "swirl";

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

  UnionIF insideSwirl(pieces);
  ComplementIF outsideSwirl(insideSwirl,true);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(outsideSwirl,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return outsideSwirl.newImplicitFunction();
}

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
  RealVect zero(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));

  BaseIF* oneVane = makeVane(thick,normal,innerRadius,outerRadius,offset,height);

  Vector<BaseIF*> eachVane;

  for (int i = 0; i < num; i++)
  {
    Real angle = i * 2*M_PI / num;

    TransformIF* currentVane = new TransformIF(*oneVane);
    currentVane->rotate(angle,zero,xAxis);

    eachVane.push_back(currentVane);
  }

  UnionIF* allVanes = new UnionIF(eachVane);

  return allVanes;
}

BaseIF* makeVane(const Real&     thick,
                 const RealVect& normal,
                 const Real&     innerRadius,
                 const Real&     outerRadius,
                 const Real&     offset,
                 const Real&     height)
{
  RealVect zero(D_DECL(0.0,0.0,0.0));
  RealVect xAxis(D_DECL(1.0,0.0,0.0));
  bool inside = true;

  Vector<BaseIF*> vaneParts;

  // Each side of the vane (infinite)
  RealVect normal1(normal);
  RealVect point1(D_DECL(offset+height/2.0,-thick/2.0,0.0));
  PlaneIF plane1(normal1,point1,inside);

  vaneParts.push_back(&plane1);

  RealVect normal2(-normal);
  RealVect point2(D_DECL(offset+height/2.0,thick/2.0,0.0));
  PlaneIF plane2(normal2,point2,inside);

  vaneParts.push_back(&plane2);

  // Make sure we only get something to the right of the origin
  RealVect normal3(D_DECL(0.0,0.0,1.0));
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
  TiltedCylinderIF inner(innerRadius,xAxis,zero,!inside);

  vaneParts.push_back(&inner);

  // The inside of the outer cylinder
  TiltedCylinderIF outer(outerRadius,xAxis,zero,inside);

  vaneParts.push_back(&outer);

  IntersectionIF* vane = new IntersectionIF(vaneParts);

  return vane;
}
