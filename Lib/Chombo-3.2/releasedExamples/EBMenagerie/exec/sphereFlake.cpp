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

#include "BaseIF.H"
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

/***************/
// Define a sphereFlake
/***************/
void makeChildSpheres(Vector<BaseIF *>& a_spheres,
                      int&              a_iter,
                      int&              a_i,
                      int&              a_numIters,
                      Real&             a_parentRadius,
                      Real&             a_shrinkRadius,
                      RealVect&         a_parentCenter,
                      Real&             a_distance);

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
      in_file = "sphereFlake2d.inputs";
    }
    else
    if (SpaceDim == 3)
    {
      in_file = "sphereFlake3d.inputs";
    }
    else
    {
      MayDay::Error("SpaceDim must be 2 or 3");
    }

    if (argc > 1 && strcmp(argv[1],"sphereFlake.inputs") != 0)
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
    const char* basename = "sphereFlake";

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
  Real     radius;
  Real     distance;
  int      numIters;
  Real     shrinkRadius;

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

  pp.get("radius",radius);
  pp.get("distance",distance);
  pp.get("numIters",numIters);
  pp.get("shrinkRadius",shrinkRadius);

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  // Stay inside until the end
  bool inside = true;

  //figure out how many spheres we will have (= sum( (SpaceDim*2)^iter) )
  int sphNumber = 0;
  for (int iter = 0; iter <= numIters; iter++)
    {
      int power = 1;
      for (int jter = 0; jter < iter; jter++)
        {
          power *= SpaceDim*2;//couldn't find a pow() for int's...
        }
      sphNumber += power;
    }

  // The vector of spheres
  Vector<BaseIF*> spheres;
  spheres.resize(sphNumber);

  //make the parent sphere
  spheres[0] = new SphereIF(radius,center,inside);

  //now make the child spheres
  int iter = 0;
  int i = 1;//this is a sphere counter
  makeChildSpheres(spheres,
                   iter,
                   i,
                   numIters,
                   radius,
                   shrinkRadius,
                   center,
                   distance);

  // Take the union of the insides of the spheres
  UnionIF sphereFlake(spheres);

  // Complement if necessary
  ComplementIF insideOut(sphereFlake,!insideRegular);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(insideOut,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return insideOut.newImplicitFunction();
}

void makeChildSpheres(Vector<BaseIF *>& a_spheres,
                      int&              a_iter,
                      int&              a_i,
                      int&              a_numIters,
                      Real&             a_parentRadius,
                      Real&             a_shrinkRadius,
                      RealVect&         a_parentCenter,
                      Real&             a_distance)
{
  Real childRadius = a_parentRadius*a_shrinkRadius;
  bool inside = true;
  if (a_iter < a_numIters)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              int lohi = sign(sit());

              RealVect center = a_parentCenter;
              //center[idir] += lohi*(a_parentRadius*(1.0 + a_shrinkRadius));
              center[idir] += lohi*(a_parentRadius + childRadius*a_distance);

              // add a child sphere to the list
              // child Sphere is the given "childRadius" and "centerChild"
              a_spheres[a_i] = new SphereIF(childRadius,center,inside);
              a_i += 1;

              //recursive call for children
              {
                a_iter += 1;
                makeChildSpheres(a_spheres,
                                 a_iter,
                                 a_i,
                                 a_numIters,
                                 childRadius,
                                 a_shrinkRadius,
                                 center,
                                 a_distance);
                a_iter -= 1;
              }
            }
        }
    }
}
