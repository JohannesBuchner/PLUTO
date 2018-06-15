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

#include "DEMIF.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

/***************/
// Define an EBIS.
/***************/
BaseIF* makeGeometry(Vector<Box>&      a_domain,
                     RealVect&         a_origin,
                     Vector<RealVect>& a_dx,
                     Vector<int>&      a_refRatio);

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
        in_file = "digitalElevationModel2d.inputs";
      }
    else if (SpaceDim == 3)
      {
        in_file = "digitalElevationModel3d.inputs";
      }
    else
      {
        MayDay::Error("SpaceDim must be 2 or 3");
      }

    if (argc > 1 && strcmp(argv[1],"digitalElevationModel.inputs") != 0)
    {
      in_file = argv[1];
    }

    // Parse input file
    ParmParse pp(0,NULL,NULL,in_file);

    Vector<Box> domain;
    RealVect origin;
    Vector<RealVect> dx;
    Vector<int> refRatio;

    BaseIF* implicit;

    // Make geometry
    implicit = makeGeometry(domain,origin,dx,refRatio);

    const int numLevels = domain.size();

    // Make grids
    Vector<DisjointBoxLayout> grids;
    makeLayouts(grids,refRatio,dx,origin,domain);

    // Create ebislayout
    int nghost = 1;
    Vector<EBISLayout> ebisl;
    makeEBISL(ebisl,grids,domain,nghost);

    // Make a LevelData
    int nComps = 1;

    IntVect ghost = IntVect::Unit;
    ghost *= nghost;

    Vector<LevelData<EBCellFAB>* > data(numLevels);
    for (int ilev = 0; ilev < numLevels; ilev++)
      {
        EBCellFactory ebcellfact(ebisl[ilev]);
        data[ilev] = new LevelData<EBCellFAB>(grids[ilev],
                                              nComps,
                                              ghost,
                                              ebcellfact);
      }

    for (int ilev = 0;ilev<numLevels;ilev++)
      {
        // Put some data in the data holders
        fillData(*data[ilev],origin,dx[ilev],*implicit,ghost);
      }

    // Done with this object
    delete implicit;

    // Write the data and the EB out
    const char* basename = "digitalElevationModel";

    char name[1000];
    sprintf(name,"%s%dd.hdf5",basename,SpaceDim);
#ifdef CH_USE_HDF5
    Vector<string> names(nComps);

    for (int i = 0; i < nComps; i++)
      {
        names[i] = "Level-Set";
      }
    CH_assert(nComps==1);

    const bool replaceCovered = false;
    Vector<Real> coveredValues(nComps,0.0);

    const Real dt = 1.e99;
    const Real time = 1.e99;
    writeEBHDF5(name,
                grids,
                data,
                names,
                domain[0],
                dx[0][0],
                dt,
                time,
                refRatio,
                numLevels,
                replaceCovered,
                coveredValues,
                ghost);
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

BaseIF* makeGeometry(Vector<Box>&      a_domain,
                     RealVect&         a_origin,
                     Vector<RealVect>& a_dx,
                     Vector<int>&      a_refRatio)
{
  pout() << "Begin making geometry " << endl;

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

  int numLevels;
  pp.get("numLevels",numLevels);
  a_domain.resize(numLevels);
  a_dx.resize(numLevels);
  a_refRatio.resize(numLevels);
  pp.getarr("ref_ratio",a_refRatio,0,numLevels);

  a_domain[0].setSmall(lo);
  a_domain[0].setBig(hi);

  Vector<Real> prob_lo(SpaceDim,1.0);
  Vector<Real> prob_hi(SpaceDim,1.0);

  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.getarr("prob_hi",prob_hi,0,SpaceDim);

  IntVect ncell;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_dx[0][idir] = (prob_hi[idir]-prob_lo[idir])/n_cell[idir];
      a_origin[idir] = prob_lo[idir];
      ncell[idir] = n_cell[idir];
    }
  pout() << "coarsest dx = " << a_dx[0][0] << endl;
  pout() << "coarsest dy = " << a_dx[0][1] << endl;
#if CH_SPACEDIM ==3
  pout() << "coarsest dz = " << a_dx[0][2] << endl;
#endif

  pout() << "dxAspect = " << a_dx[0][0]/a_dx[0][SpaceDim-1] << endl;
  pout() << "dyAspect = " << a_dx[0][1]/a_dx[0][SpaceDim-1] << endl;
#if CH_SPACEDIM ==3
  pout() << "dzAspect = " << a_dx[0][2]/a_dx[0][SpaceDim-1] << endl;
#endif

  //refine down to finer levels
  for (int ilev = 1;ilev<numLevels;ilev++)
    {
      a_domain[ilev] = a_domain[ilev-1];
      a_domain[ilev].refine(a_refRatio[ilev]);

      a_dx[ilev] = a_dx[ilev-1];
      a_dx[ilev] /= Real(a_refRatio[ilev]);
    }
  RealVect finestDx     = a_dx[numLevels-1];
  Box      finestDomain = a_domain[numLevels-1];

  //interpType = 1 -> bilinear interpolation
  //interpType = 2 -> biquadratic interpolation
  //interpType = 3 -> bicubic interpolation type1 (this performs cubic interpolation in 1D, matching all value)
  //interpType = 4 -> bicubic interpolation type2 (this performs cubic interpolation in 1D, matching values and derivatives)
  int interpType;
  pp.get("interpType",interpType);

  //bottomBuffer is space added below the bathymetry,
  //  (the distance from the deepest spot to the domain box)
  Real bottomBuffer;
  pp.get("bottomBuffer",bottomBuffer);

  //verticalScale is used for testing anisotropic vs isotropic geometry, if verticalScale=1.0 then this is the true geometry;
  //   if verticalScale is 100 then all elevations are multiplied by 100...
  Real verticalScale;
  pp.get("verticalScale",verticalScale);

  //highGround is the elevation given for nodata points with all land neighbors
  //  (useful for higher order interpolation)
  Real highGround;
  pp.get("highGround",highGround);
  Real truncElev = 1.e99;

  //get the file name for the Digital Elevation Model
  std::string demFile;
  pp.get("DEM_file",demFile);

  DEMIF implicit(ncell,
                 interpType,
                 finestDx,
                 demFile,
                 bottomBuffer,
                 truncElev,
                 highGround,
                 verticalScale);

  int verbosity = 0;
  GeometryShop workshop(implicit,
                        verbosity,
                        finestDx);

  // This generates the new EBIS
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();

  int loadEBISFile;
  pp.get("loadEBISFile",loadEBISFile);
  if (loadEBISFile==1)
    {
      //get ebindexspace
      pout() << " recreating  geometry from file: dem.ebis " << endl;
      //define ebis from file input
#ifdef CH_USE_HDF5
      HDF5Handle handleIn("dem.ebis", HDF5Handle::OPEN_RDONLY);
      ebisPtr->define(handleIn);
      handleIn.close();
#endif
    }
  else
    {
      ebisPtr->define(finestDomain,a_origin,finestDx[0],workshop);
    }

  //write EBIndexSpace to a file
  {
    //get ebindexspace
    pout() << " writing  geometry to file: dem.ebis" << endl;
#ifdef CH_USE_HDF5
    HDF5Handle handle("dem.ebis", HDF5Handle::CREATE);
    ebisPtr->write(handle);
    handle.close();
#endif
  }

  return implicit.newImplicitFunction();
}
