#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Purpose:
//  Test the infrastucture for periodic boundary conditions
//
// Usage:
//  <program-name> [-q|-v] ...
//
//  where:
//    -q means run quietly (only pass/fail messages printed)
//    -v means run verbosely (all messages printed)
//    -writePlots means write hdf5 plotfile outputs
//    -noPlots means no hdf5 plotfile output (default)
//
//    ... all non-option arguments are ignored (Chombo convention)
//
//  Default is `-v'
//
//  Unknown options are treated as errors.
//

// Include files:
#include <cstdio>

#include "parstream.H"
#include "LoadBalance.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H"
#include "DebugDump.H"

#include "FineInterp.H"
#include "QuadCFInterp.H"
#include "ExtrapFillPatch.H"
#include "AMRIO.H"
#include "UsingNamespace.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

void initData(LevelData<FArrayBox>& a_data, const Real a_dx);

/// Global variables for handling output

static const char *pgmname = "testPeriodic";
static const char *indent = "   ";
//static const char *indent2 = "      " ;

static bool verbose;
static bool writePlots = false;

/// Code:

int
main(int argc ,char *argv[] )
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  int status = 0;
  // Make this test pass automatically if DIM > 3, to avoid FAIL message.
  if (SpaceDim <= 3)
  {
    parseTestOptions(argc, argv);

    // establish periodic domain -- first multiply periodic in all directions
    int baseDomainSize = 8;
    int numGhost = 2;
    IntVect ghostVect(numGhost*IntVect::Unit);

    Box baseDomBox(IntVect::Zero, (baseDomainSize-1)*IntVect::Unit);
    ProblemDomain baseDomain(baseDomBox);

    // set periodic in all directions
    for (int dir=0; dir<SpaceDim; dir++)
    {
      baseDomain.setPeriodic(dir, true);
    }

    {
      // QuadCFInterp test
      int nRef = 2;
      ProblemDomain fineDomain(baseDomain);
      fineDomain.refine(nRef);

      const Box domainBox = baseDomain.domainBox();
      Vector<Box> crseBoxes(1, domainBox);
      Vector<int> crseProcAssign(1,0);
      DisjointBoxLayout crseGrids(crseBoxes, crseProcAssign, baseDomain);

      // note the lack of ghost cells on the base level
      LevelData<FArrayBox> crseData(crseGrids, 1, IntVect::Zero);

      Real dxCrse = 1.0/(baseDomainSize);
      Real dxFine = dxCrse/nRef;
      initData(crseData, dxCrse);

#if (CH_SPACEDIM == 1)
      Vector<Box> fineBoxes(2);
      Vector<int> fineProcAssign(2);

      fineBoxes[0] = Box(IntVect::Zero, 5*IntVect::Unit);
      fineBoxes[1] = Box(IntVect(D_DECL6(6,10,10,6,10,0)), IntVect(D_DECL6(11,15,15,11,15,15)));
#else
      // try multiple boxes for fine level
      Vector<Box> fineBoxes(4);
      Vector<int> fineProcAssign(4);

      fineBoxes[0] = Box(IntVect::Zero, 9*IntVect::Unit);
      fineBoxes[1] = Box(IntVect(D_DECL6(6,10,10,6,10,0)), IntVect(D_DECL6(11,15,15,11,15,15)));
      // this box should fail the disjointness test
      //fineBoxes[1] = Box(IntVect(D_DECL(7,10,10)), 17*IntVect::Unit);
      fineBoxes[2] = Box(IntVect(D_DECL6(12,2,0,12,2,0)), IntVect(D_DECL6(15,7,7,15,7,7)));
      fineBoxes[3] = Box(IntVect(D_DECL6(12,12,0,12,12,0)), 15*IntVect::Unit);
#endif

      int loadbalancestatus = LoadBalance(fineProcAssign, fineBoxes);
      CH_assert (loadbalancestatus == 0);

      DisjointBoxLayout fineGrids;
      fineGrids.define(fineBoxes, fineProcAssign, fineDomain);
      LevelData<FArrayBox> fineData(fineGrids, 1, ghostVect);

      // set these to a bogus value to start with
      DataIterator fineDit = fineData.dataIterator();
      for (fineDit.begin(); fineDit.ok(); ++fineDit)
        {
          fineData[fineDit()].setVal(1.0e9);
        }

      FineInterp interpolator(fineGrids, 1, nRef, fineDomain);

      // turn on fine-data averaging even though it's unnecessary
      // since this is a test (so we'd like to know if it's broken)
      bool averageFineData = true;
      interpolator.interpToFine(fineData, crseData, averageFineData);

      // now fill in ghost cells
      QuadCFInterp filler(fineGrids, &crseGrids, dxFine, nRef, 1,
              fineDomain);

      filler.coarseFineInterp(fineData, crseData);
    }

    {
      // fineInterp near boundary test
      int numLev = 3;
      Vector<int> nRefVect(numLev);
      nRefVect[0] = 4;
      nRefVect[1] = 2;
      nRefVect[2] = -1;
      Vector<Real> dxVect(numLev);
      Vector<ProblemDomain> levelDomains(numLev);
      levelDomains[0] = Box(IntVect::Zero, 15*IntVect::Unit);
      dxVect[0] = 1.0/(levelDomains[0].domainBox().size(0));
      for (int lev=1; lev<numLev; lev++)
        {
          levelDomains[lev] = levelDomains[lev-1];
          levelDomains[lev].refine(nRefVect[lev-1]);
          dxVect[lev] = dxVect[lev-1]/nRefVect[lev-1];
        }

      Vector<DisjointBoxLayout> levelGrids(numLev);
      Vector<LevelData<FArrayBox>* > levelData(numLev, NULL);

      const Box domainBox = levelDomains[0].domainBox();
      Vector<Box> crseBoxes(1, domainBox);
      Vector<int> crseProcAssign(1,0);
      DisjointBoxLayout level0Grids(crseBoxes, crseProcAssign,
                                    levelDomains[0]);
      levelGrids[0] = level0Grids;

      // note the lack of ghost cells on the base level
      levelData[0] = new LevelData<FArrayBox>(level0Grids,
                                              1, IntVect::Zero);

      initData(*levelData[0], dxVect[0]);

      Vector<Box> level1Boxes(1, levelDomains[1].domainBox());
      int loadbalancestatus = LoadBalance(crseProcAssign, level1Boxes);
      CH_assert (loadbalancestatus == 0);

      IntVect newGhostVect(4*IntVect::Unit);
      DisjointBoxLayout fineGrids;
      fineGrids.define(level1Boxes, crseProcAssign, levelDomains[1]);
      levelGrids[1] = fineGrids;
      levelData[1] = new LevelData<FArrayBox>(fineGrids, 1,
                                              newGhostVect);

      Vector<Box> level2Boxes(1);
      //level2Boxes[0] = Box(20*IntVect::Unit, 119*IntVect::Unit);
      level2Boxes[0] = Box(100*IntVect::Unit, 127*IntVect::Unit);
      DisjointBoxLayout finestGrids(level2Boxes, crseProcAssign,
                                    levelDomains[2]);
      levelGrids[2] = finestGrids;
      levelData[2] = new LevelData<FArrayBox>(finestGrids, 1,
                                              newGhostVect);

      for (int lev=1; lev<numLev; lev++)
        {
          LevelData<FArrayBox>& fineData = *levelData[lev];
          LevelData<FArrayBox>& crseData = *levelData[lev-1];

          // set these to a bogus value to start with
          DataIterator fineDit = fineData.dataIterator();
          for (fineDit.begin(); fineDit.ok(); ++fineDit)
            {
              fineData[fineDit()].setVal(1.0e9);
            }

          FineInterp interpolator(levelGrids[lev], 1, nRefVect[lev-1],
                                  levelDomains[lev]);

          // turn on fine-data averaging even though it's unnecessary
          // since this is a test (so we'd like to know if it's broken)
          bool averageFineData = true;
          interpolator.interpToFine(fineData, crseData, averageFineData);

          pout() << "grids at level = " << lev << endl;
          dumpDBL(&levelGrids[lev]);
          pout() << "grids at level = " << lev-1 << endl;
          dumpDBL(&levelGrids[lev-1]);

          // now fill in ghost cells
          QuadCFInterp filler(levelGrids[lev], &levelGrids[lev-1],
                              dxVect[lev],
                              nRefVect[lev-1], 1,
                              levelDomains[lev]);

          filler.coarseFineInterp(fineData, crseData);
        }

#ifdef CH_USE_HDF5
      // dump data
      if (writePlots)
        {
          string fname = "interp.hdf5";
          Vector<string> varNames(1,"phi");
          Real dt = 0.0;
          Real time= 0.0;
          WriteAMRHierarchyHDF5(fname, levelGrids,
                                levelData, varNames,
                                levelDomains[0].domainBox(),
                                dxVect[0], dt, time,
                                nRefVect, numLev);
        }
#endif

      // clean up memory
      for (int lev=0; lev<levelData.size(); lev++)
        {
          if (levelData[lev] != NULL)
            {
              delete levelData[lev];
              levelData[lev] = NULL;
            }
        }
    }

    // now test partially periodic case (also coarsen back to original domain
    // no such thing as "partially periodic" for 1D; for 1D, test moving
    // fine-grid boxes from low-end of domain to high-end of domain

    {

      baseDomain.setPeriodic(0,false);
      // try similar tests to the previous case
      int nRef = 2;
      ProblemDomain fineDomain(baseDomain);
      fineDomain.refine(nRef);

      const Box domainBox = baseDomain.domainBox();
      Vector<Box> crseBoxes(1, domainBox);
      Vector<int> crseProcAssign(1,0);
      DisjointBoxLayout crseGrids(crseBoxes, crseProcAssign, baseDomain);

      // note the lack of ghost cells on the base level
      LevelData<FArrayBox> crseData(crseGrids, 1, IntVect::Zero);

      Real dxCrse = 1.0/(baseDomainSize);
      Real dxFine = dxCrse/nRef;
      initData(crseData, dxCrse);

#if (CH_SPACEDIM == 1)
      Vector<Box> fineBoxes(2);
      Vector<int> fineProcAssign(2);

      fineBoxes[0] = Box(12*IntVect::Unit, 15*IntVect::Unit);
      fineBoxes[1] = Box(IntVect(D_DECL6(6,10,10,6,10,10)), IntVect(D_DECL6(11,15,15,11,15,15)));
#else
      // try multiple boxes for fine level
      // note grids are different in this case to allow for all periodic checking
      Vector<Box> fineBoxes(3);
      Vector<int> fineProcAssign(3);

      fineBoxes[0] = Box(IntVect::Zero, 9*IntVect::Unit);
      fineBoxes[1] = Box(IntVect(D_DECL6(6,10,10,6,10,10)), IntVect(D_DECL6(11,15,15,11,15,15)));
      // this box should fail the disjointness test
      //fineBoxes[1] = Box(IntVect(D_DECL(7,10,10)), 17*IntVect::Unit);
      fineBoxes[2] = Box(IntVect(D_DECL6(12,2,0,12,2,0)), IntVect(D_DECL6(15,7,7,15,7,7)));
      //fineBoxes[3] = Box(IntVect(D_DECL(12,12,0)), 15*IntVect::Unit);
#endif

      int loadbalancestatus = LoadBalance(fineProcAssign, fineBoxes);
      CH_assert (loadbalancestatus == 0);

      DisjointBoxLayout fineGrids;
      fineGrids.define(fineBoxes, fineProcAssign, fineDomain);
      LevelData<FArrayBox> fineData(fineGrids, 1, ghostVect);

      // set these to a bogus value to start with
      DataIterator fineDit = fineData.dataIterator();
      for (fineDit.begin(); fineDit.ok(); ++fineDit)
        {
          fineData[fineDit()].setVal(1.0e9);
        }

      FineInterp interpolator(fineGrids, 1, nRef, fineDomain);

      // turn on fine-data averaging even though it's unnecessary
      // since this is a test (so we'd like to know if it's broken)
      bool averageFineData = true;
      interpolator.interpToFine(fineData, crseData, averageFineData);

      // now fill in first layer of ghost cells
      QuadCFInterp filler(fineGrids, &crseGrids, dxFine, nRef, 1,
                          fineDomain);

      filler.coarseFineInterp(fineData, crseData);

      // now fill in remaining ghost cell with extrapolation
      Interval extrapInterval(1,1);
      ExtrapFillPatch extrapFiller(fineGrids, fineDomain,
                                   extrapInterval);

      for (int dir=0; dir<SpaceDim; dir++)
        {
          extrapFiller.fillExtrap(fineData,dir,0,1);
        }

      fineData.exchange(fineData.interval());

    }

  } // if (SpaceDim <= 3)
  pout() << indent << pgmname << ": "
         << ( (status == 0) ? "passed all tests" : "failed at least one test,")
         << endl;

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status ;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if ( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-writePlots" ,3 ) == 0 )
            {
              writePlots = true;
              // argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-noPlots" ,3 ) == 0 )
            {
              writePlots = false;
              // argv[i] = "" ;
            }
          else
            {
              break ;
            }
        }
    }
  return ;
}

void
initData(LevelData<FArrayBox>& a_data, const Real a_dx)
{

  DataIterator dit = a_data.dataIterator();
  const DisjointBoxLayout& interiorBoxes = a_data.getBoxes();

  for (dit.begin(); dit.ok(); ++dit)
    {
      // first set to a bogus value which will persist in ghost cells
      // after initialization
      a_data[dit()].setVal(1.0e9);

      // this will be slow, but who cares?
      FArrayBox& localData = a_data[dit()];
      BoxIterator boxIt(interiorBoxes[dit()]);
      Real localVal;
      for (boxIt.begin(); boxIt.ok(); ++boxIt)
        {
          const IntVect& loc = boxIt();
          // linear profile
          //localVal = a_dx*(D_TERM(loc[0]+0.5,
          //                      + loc[1]+0.5,
          //                      + loc[2]+0.5));

          // quadratic profile
          localVal = a_dx*a_dx*(D_TERM6((loc[0]+0.5)*(loc[0]+0.5),
                                        + (loc[1]+0.5)*(loc[1]+0.5),
                                        + (loc[2]+0.5)*(loc[2]+0.5),
                                        + (loc[3]+0.5)*(loc[3]+0.5),
                                        + (loc[4]+0.5)*(loc[4]+0.5),
                                        + (loc[5]+0.5)*(loc[5]+0.5)));
          localData(loc) = localVal;
        }

    }
}
