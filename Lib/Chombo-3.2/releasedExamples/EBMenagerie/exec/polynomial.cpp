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

#include "PolynomialIF.H"

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
    const char* in_file = "polynomial.inputs";

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
    const char* basename = "polynomial";

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

  // Parm Parse doesn't get bools, so work-around with int
  int intInsideRegular;
  pp.get("insideRegular",intInsideRegular);

  if (intInsideRegular != 0) insideRegular = true;
  if (intInsideRegular == 0) insideRegular = false;

  // Read the maximum total power for the polynomial
  int totalMaxPower;
  pp.get("power",totalMaxPower);

  // Generate all the powers of x, y, and z (3D) where the power of x varies
  // most rapidly, followed by y and z (3D) and the sum of the powers is
  // always less than totalMaxPower.
  Vector<IntVect> powers;
  int total = 0;

  // State of the iteration and the maximum power for each dimension
  IntVect state = IntVect::Zero;
  IntVect maxPower;

  // The final dimension can make it to the maximum power
  maxPower[SpaceDim-1] = totalMaxPower;

  // Loop over all possibilities
  while (state[SpaceDim-1] <= maxPower[SpaceDim-1])
  {
    // Determine the maximum power each dimension can currently have
    int curMaxPower = totalMaxPower;
    for (int i = SpaceDim-1; i >= 0; i--)
    {
      maxPower[i] = curMaxPower;
      curMaxPower -= state[i];
    }

    // Make room of the next set of powers
    total++;
    powers.resize(total);

    // The current state gives the powers in each dimension
    powers[total-1] = state;

    // Increment the x power - fastiest varying
    state[0]++;

    // If greater that the maximum it can be...
    if (state[0] > maxPower[0])
    {
      // For all dimensions...
      for (int i = 0; i < SpaceDim-1; i++)
      {
        // Zero the current power and increment the next dimension's power
        state[i] = 0;
        state[i+1]++;

        // If the next dimension's power is okay go on
        if (state[i+1] <= maxPower[i+1])
        {
          break;
        }
      }
    }
  }

  // Allocate space for the polynomial coefficients and read them in
  Vector<Real> coefs;
  coefs.resize(total);

  pp.getarr("coefs",coefs,0,total);

  Vector<PolyTerm> polynomial;
  polynomial.resize(total);

  // Set up the polynomial
  for (int i = 0; i < total; i++)
  {
    polynomial[i].coef   = coefs[i];
    polynomial[i].powers = powers[i];
  }

  PolynomialIF implicit(polynomial,insideRegular);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_dx;

  GeometryShop workshop(implicit,0,vectDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  ebisPtr->define(a_domain,a_origin,a_dx,workshop);

  return implicit.newImplicitFunction();
}
