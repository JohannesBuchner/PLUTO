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

#include "GaussianIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

/***************/
// Define an EBIS.
/***************/
BaseIF* makeGeometry(Box&      a_domain,
                     RealVect& a_origin,
                     RealVect& a_dx);

int main (int argc, char** argv)
{
#ifdef CH_MPI
  MPI_Init(&argc,&argv);
#endif

  // Begin forever present scoping trick
  {
    const char* in_file = "gaussian.inputs";

    if (argc > 1)
    {
      in_file = argv[1];
    }

    // Parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    Box domain;
    RealVect origin;
    RealVect dx;

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
    const char* basename = "gaussian";

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
                     RealVect& a_dx)
{
  RealVect point;
  RealVect normal;
  bool     convex;

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
  Vector<Real> prob_hi(SpaceDim,1.0);

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.getarr("prob_hi",prob_hi,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_dx[idir] = (prob_hi[idir]-prob_lo[idir])/n_cell[idir];
      a_origin[idir] = prob_lo[idir];
    }

  // ParmParse doesn't get RealVects, so work-around with Vector<Real>
  Vector<Real> vectorCenter;
  Vector<Real> vectorSig2;
  pp.getarr("origin",vectorCenter,0,SpaceDim);
  pp.getarr("sig2"  ,vectorSig2  ,0,SpaceDim);
  RealVect origin;
  RealVect sig2;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      origin[idir] = vectorCenter[idir];
      sig2[  idir] = vectorSig2[idir];
    }

  Real height;
  pp.get("height",height);

  int upDir;
  pp.get("upDir",upDir);

  pp.get("convex",convex);

  GaussianIF implicit(origin,height,sig2,upDir,convex);

  GeometryShop workshop(implicit,0,a_dx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx[0],workshop);

  return implicit.newImplicitFunction();
}
