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

#include "SphereIF.H"
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

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file = "twoSpheres.inputs";

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
    const char* basename = "twoSpheres";

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
  RealVect center1;
  Real     radius1;

  RealVect center2;
  Real     radius2;

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
  pp.getarr("center1",vectorCenter,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    center1[idir] = vectorCenter[idir];
  }

  pp.get("radius1",radius1);

  pp.getarr("center2",vectorCenter,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    center2[idir] = vectorCenter[idir];
  }

  pp.get("radius2",radius2);

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  // Stay inside until the end
  bool inside = true;

  // Sphere 1 is the given "radius1" and "center1"
  SphereIF sphere1(radius1,center1,inside);

  // Sphere 2 is the given "radius2" and "center2"
  SphereIF sphere2(radius2,center2,inside);

  // Take the union of the spheres
  UnionIF implicit(sphere1,sphere2);

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
