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
//    -write means write out the boxes; otherwise will read in
//           boxes from a file and compare with boxes generated
//           in tests.
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
#include "FluxBox.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H"
#include "DebugDump.H"

#include "LevelFluxRegister.H"

#ifdef CH_MPI
#include "mpi.h"
#endif

#include "UsingNamespace.H"

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

void initData(LevelData<FArrayBox>& a_data, const Real a_dx);

void initEdgeData(LevelData<FluxBox>& a_data, const Real a_dx);

/// Global variables for handling output

static const char *pgmname = "testPeriodic";
static const char *indent  = "   "    ;
// static const char *indent2 = "      " ;

static bool verbose = true ;

/// Code:

int
main(int argc ,char *argv[] )
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions(argc, argv);

  // establish periodic domain -- first multiply periodic in all directions
  int status = 0;
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
    // doubly periodic flux register test
    int nRef = 2;
    ProblemDomain fineDomain(baseDomain);
    fineDomain.refine(nRef);

    const Box domainBox = baseDomain.domainBox();
    //Vector<Box> crseBoxes(1, domainBox);
    int numBoxes = 4;
    if (SpaceDim == 1)
      {
        numBoxes = 2;
      }
    Vector<Box> crseBoxes(numBoxes);
    crseBoxes[0] = Box(IntVect::Zero, 4*IntVect::Unit);
    crseBoxes[1] = Box(5*IntVect::Unit, 7*IntVect::Unit);
    if (SpaceDim > 1)
      {
        crseBoxes[2] = Box(IntVect(D_DECL6(0,5,0,0,0,0)),
                           IntVect(D_DECL6(4,7,7,7,7,7)));
        crseBoxes[3] = Box(IntVect(D_DECL6(5,0,0,0,0,0)),
                           IntVect(D_DECL6(7,4,7,7,7,7)));
      }


    Vector<int> crseProcAssign(crseBoxes.size(),0);
    DisjointBoxLayout crseGrids(crseBoxes, crseProcAssign, baseDomain);

    // note the lack of ghost cells on the base level
    LevelData<FArrayBox> crseData(crseGrids, 1, IntVect::Zero);
    LevelData<FluxBox> crseFlux(crseGrids, 1, IntVect::Zero);

    Real dxCrse = 1.0/(baseDomainSize);
    Real dxFine = dxCrse/nRef;

    // try multiple boxes for fine level
    int numFineBoxes = 4;
    if (SpaceDim == 1)
      {
        numFineBoxes = 2;
      }
    Vector<Box> fineBoxes(numFineBoxes);
    Vector<int> fineProcAssign(numFineBoxes);

    fineBoxes[0] = Box(IntVect::Zero, 9*IntVect::Unit);
    fineBoxes[1] = Box(IntVect(D_DECL6(12,12,0,0,0,0)),
                       15*IntVect::Unit);
    if (SpaceDim > 1)
      {
        fineBoxes[2] = Box(IntVect(D_DECL6(6,10,10,10,10,10)),
                           IntVect(D_DECL6(11,15,15,15,15,15)));
        // this box should fail the disjointness test
        //fineBoxes[1] = Box(IntVect(D_DECL(7,10,10)), 17*IntVect::Unit);
        fineBoxes[3] = Box(IntVect(D_DECL6(12,2,0,0,0,0)),
                           IntVect(D_DECL6(15,7,7,7,7,7)));
      }

    int loadbalancestatus = LoadBalance(fineProcAssign, fineBoxes);
    CH_assert (loadbalancestatus == 0);

    DisjointBoxLayout fineGrids;
    fineGrids.define(fineBoxes, fineProcAssign, fineDomain);
    LevelData<FluxBox> fineFlux(fineGrids, 1, IntVect::Zero);

    // set coarse fluxes to 1, fine fluxes to 2 to start
    // coarse-grid cell-centered data is set to 0 to see effect of reflux
    // this tests reflux locations...
    DataIterator fineDit = fineFlux.dataIterator();
    for (fineDit.begin(); fineDit.ok(); ++fineDit)
      {
        fineFlux[fineDit()].setVal(2.0);
      }

    DataIterator crseDit = crseData.dataIterator();
    for (crseDit.begin(); crseDit.ok(); ++crseDit)
      {
        crseFlux[crseDit()].setVal(1.0);
        crseData[crseDit()].setVal(0.0);
      }

    LevelFluxRegister testFR(fineGrids, crseGrids,
                             fineDomain, nRef, 1);

    testFR.setToZero();

    // increment coarse fluxes
    Real crseScale = 1.0;
    Interval comps = crseData.interval();
    for (crseDit.begin(); crseDit.ok(); ++crseDit)
      {
        FluxBox& thisCrseFlux = crseFlux[crseDit()];
        for (int dir=0; dir<SpaceDim; ++dir)
          {
            FArrayBox& thisCrseFluxDir = thisCrseFlux[dir];
            testFR.incrementCoarse(thisCrseFluxDir,
                                   crseScale,
                                   crseDit(),
                                   comps, comps, dir);
          }
      }

    // increment fine fluxes
    Real fineScale = 1.0;
    for (fineDit.begin(); fineDit.ok(); ++fineDit)
      {
        FluxBox& thisFineFlux = fineFlux[fineDit()];
        for (int dir=0; dir<SpaceDim; ++dir)
          {
            FArrayBox& thisFineFluxDir = thisFineFlux[dir];

            testFR.incrementFine(thisFineFluxDir,
                                 fineScale, fineDit(),
                                 comps, comps, dir,
                                 Side::Lo);

            testFR.incrementFine(thisFineFluxDir,
                                 fineScale, fineDit(),
                                 comps, comps, dir,
                                 Side::Hi);
          }
      }

    // now apply reflux correction
    Real refluxScale = 1.0;
    testFR.reflux(crseData, refluxScale);

    // now test accuracy (of locations) -- redo test with linear
    // profiles.  all flux contributions should cancel out, leaving
    // no refluxing contributions.

    initEdgeData(fineFlux, dxFine);
    initEdgeData(crseFlux, dxCrse);

    for (crseDit.begin(); crseDit.ok(); ++crseDit)
      {
        crseData[crseDit()].setVal(0.0);
      }

    testFR.setToZero();

    // increment coarse fluxes
        // increment coarse fluxes
    crseScale = 1.0;
    for (crseDit.begin(); crseDit.ok(); ++crseDit)
      {
        FluxBox& thisCrseFlux = crseFlux[crseDit()];
        for (int dir=0; dir<SpaceDim; ++dir)
          {
            FArrayBox& thisCrseFluxDir = thisCrseFlux[dir];
            testFR.incrementCoarse(thisCrseFluxDir,
                                   crseScale,
                                   crseDit(),
                                   comps, comps, dir);
          }
      }

    // increment fine fluxes
    fineScale = 1.0;
    for (fineDit.begin(); fineDit.ok(); ++fineDit)
      {
        FluxBox& thisFineFlux = fineFlux[fineDit()];
        for (int dir=0; dir<SpaceDim; ++dir)
          {
            FArrayBox& thisFineFluxDir = thisFineFlux[dir];

            testFR.incrementFine(thisFineFluxDir,
                                 fineScale, fineDit(),
                                 comps, comps, dir,
                                 Side::Lo);

            testFR.incrementFine(thisFineFluxDir,
                                 fineScale, fineDit(),
                                 comps, comps, dir,
                                 Side::Hi);
          }
      }

    // now apply reflux correction
    refluxScale = 1.0;
    testFR.reflux(crseData, refluxScale);
  }

  {
    // now test partially periodic case (also coarsen back to original domain
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
    LevelData<FluxBox> crseFlux(crseGrids, 1, IntVect::Zero);

    // Real dxCrse = 1.0/(baseDomainSize);
    // Real dxFine = dxCrse/nRef;

    // try multiple boxes for fine level
    // note grids are different in this case to allow for all periodic checking
    int numFineBoxes = 3;
    if (SpaceDim == 1)
      {
        numFineBoxes = 2;
      }
    Vector<Box> fineBoxes(numFineBoxes);
    Vector<int> fineProcAssign(numFineBoxes);

    fineBoxes[0] = Box(IntVect::Zero, 9*IntVect::Unit);
    fineBoxes[1] = Box(IntVect(D_DECL6(12,2,0,0,0,0)),
                       IntVect(D_DECL6(15,7,7,7,7,7)));
    if (SpaceDim > 1)
      {
        fineBoxes[2] = Box(IntVect(D_DECL6(6,10,10,10,10,10)),
                           IntVect(D_DECL6(11,15,15,15,15,15)));
      }
    // this box should fail the disjointness test
    //fineBoxes[1] = Box(IntVect(D_DECL(7,10,10)), 17*IntVect::Unit);
    //fineBoxes[3] = Box(IntVect(D_DECL(12,12,0)), 15*IntVect::Unit);

    int loadbalancestatus = LoadBalance(fineProcAssign, fineBoxes);
    CH_assert (loadbalancestatus == 0);

    DisjointBoxLayout fineGrids;
    fineGrids.define(fineBoxes, fineProcAssign, fineDomain);
    LevelData<FluxBox> fineFlux(fineGrids, 1, IntVect::Zero);

    // set coarse fluxes to 1, fine fluxes to 2 to start
    // coarse-grid cell-centered data is set to 0 to see effect of reflux
    DataIterator fineDit = fineFlux.dataIterator();
    for (fineDit.begin(); fineDit.ok(); ++fineDit)
      {
        fineFlux[fineDit()].setVal(2.0);
      }

    DataIterator crseDit = crseData.dataIterator();
    for (crseDit.begin(); crseDit.ok(); ++crseDit)
      {
        crseFlux[crseDit()].setVal(1.0);
        crseData[crseDit()].setVal(0.0);
      }

    LevelFluxRegister testFR(fineGrids, crseGrids,
                             fineDomain, nRef, 1);

    testFR.setToZero();

    // increment coarse fluxes
    Real crseScale = 1.0;
    Interval comps = crseData.interval();
    for (crseDit.begin(); crseDit.ok(); ++crseDit)
      {
        FluxBox& thisCrseFlux = crseFlux[crseDit()];
        for (int dir=0; dir<SpaceDim; ++dir)
          {
            FArrayBox& thisCrseFluxDir = thisCrseFlux[dir];
            testFR.incrementCoarse(thisCrseFluxDir,
                                   crseScale,
                                   crseDit(),
                                   comps, comps, dir);
          }
      }

    // increment fine fluxes
    Real fineScale = 1.0;
    for (fineDit.begin(); fineDit.ok(); ++fineDit)
      {
        FluxBox& thisFineFlux = fineFlux[fineDit()];
        for (int dir=0; dir<SpaceDim; ++dir)
          {
            FArrayBox& thisFineFluxDir = thisFineFlux[dir];

            testFR.incrementFine(thisFineFluxDir,
                                 fineScale, fineDit(),
                                 comps, comps, dir,
                                 Side::Lo);

            testFR.incrementFine(thisFineFluxDir,
                                 fineScale, fineDit(),
                                 comps, comps, dir,
                                 Side::Hi);
          }
      }

    // now apply reflux correction
    Real refluxScale = 1.0;
    testFR.reflux(crseData, refluxScale);
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  pout() << indent << pgmname << ": "
     << ( (status == 0) ? "passed all tests" : "failed at least one test,")
     << " status = " << status
     << endl;

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
                                        + (loc[5]+0.5)*(loc[5]+0.5) ) );
          localData(loc) = localVal;
        }

    }
}

void
initEdgeData(LevelData<FluxBox>& a_data, const Real a_dx)
{
  DataIterator dit = a_data.dataIterator();
  const DisjointBoxLayout& interiorBoxes = a_data.getBoxes();

  for (dit.begin(); dit.ok(); ++dit)
    {
      // first set to a bogus value which will persist in ghost cells
      // after initialization
      a_data[dit()].setVal(1.0e9);

      for (int dir=0; dir<SpaceDim; dir++)
    {
      // this will be slow, but who cares?
      FArrayBox& localData = a_data[dit()][dir];
      Box edgeBox = interiorBoxes[dit()];
      edgeBox.surroundingNodes(dir);
      BoxIterator boxIt(edgeBox);
      Real localVal;
      for (boxIt.begin(); boxIt.ok(); ++boxIt)
        {
          const IntVect& loc = boxIt();
          Real realLoc[SpaceDim];
          D_TERM6(realLoc[0] = a_dx*(loc[0]+0.5);,
                  realLoc[1] = a_dx*(loc[1]+0.5);,
                  realLoc[2] = a_dx*(loc[2]+0.5);,
                  realLoc[3] = a_dx*(loc[3]+0.5);,
                  realLoc[4] = a_dx*(loc[4]+0.5);,
                  realLoc[5] = a_dx*(loc[5]+0.5); )

          // now correct for edge centering
          realLoc[dir] = realLoc[dir] -0.5*a_dx;

          // linear profile
          localVal = D_TERM6(realLoc[0],
                             +realLoc[1],
                             +realLoc[2],
                             +realLoc[3],
                             +realLoc[4],
                             +realLoc[5] );

          // quadratic profile
          //localVal = a_dx*a_dx*(D_TERM((loc[0]+0.5)*(loc[0]+0.5),
          //                   + (loc[1]+0.5)*(loc[1]+0.5),
          //                   + (loc[2]+0.5)*(loc[2]+0.5)));
          localData(loc) = localVal;
        }

    }
    }
}
