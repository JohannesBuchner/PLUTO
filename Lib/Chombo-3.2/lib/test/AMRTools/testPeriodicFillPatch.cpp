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
#include "LevelData.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H"
#include "AMRIO.H"
#include "CONSTANTS.H"

#include "FABView.H"

#include "FineInterp.H"
#include "PiecewiseLinearFillPatch.H"

#ifdef CH_MPI
#include "mpi.h"
#endif
#include "UsingNamespace.H"

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

void initData(LevelData<FArrayBox>& a_data, const Real a_dx);

Real getSolutionVal(const RealVect& loc)
{
  // need a function which is linear, but also periodic. This
  // works, as long as we don't have a box which spans the centerline
  Real retVal = 0;
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (loc[dir] <0.5)
        {
          retVal += (dir+1)*loc[dir];
        }
      else
        {
          retVal += (dir+1)*(loc[dir]-1.0);
        }
    }
  return retVal;
}

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

  int status = 0;
  // Make this test pass automatically if DIM > 3, to avoid FAIL message.
  if (SpaceDim <= 3)
  {
    // establish periodic domain -- first multiply periodic in all directions
    int baseDomainSize = 16;
    int numGhost = 2;
    IntVect ghostVect(numGhost*IntVect::Unit);

    Box baseDomBox(IntVect::Zero, (baseDomainSize-1)*IntVect::Unit);
    ProblemDomain baseDomain(baseDomBox);

    // set periodic in all directions
    for (int dir=0; dir<SpaceDim; dir++)
    {
      baseDomain.setPeriodic(dir, true);
    }

#if 0
    // set up physical BC class for partially periodic tests
    DomainGhostBC physBC;
    for (int dir=0; dir<SpaceDim; dir++)
    {
      SideIterator sit;
      for (sit.begin(); sit.ok(); ++sit)
        {
          // it should be easy to see if Neumann bcs are doing the right thing
          NeumannBC thisBC(dir,sit());
          physBC.setBoxGhostBC(thisBC);
        }
    }
#endif

    {
      // FineInterp test
      int nRef = 2;
      int fineDomainSize = nRef*baseDomainSize;
      ProblemDomain fineDomain(baseDomain);
      fineDomain.refine(nRef);

      const Box domainBox = baseDomain.domainBox();
      Vector<Box> crseBoxes(1, domainBox);
      Vector<int> crseProcAssign(1,0);
      DisjointBoxLayout crseGrids(crseBoxes, crseProcAssign, baseDomain);

      // note the lack of ghost cells on the base level
      LevelData<FArrayBox> crseData(crseGrids, 1, IntVect::Zero);

      Real dx = 1.0/(baseDomainSize);
      initData(crseData, dx);

      // try multiple boxes for fine level
      int numBoxes = 3;
      Vector<Box> fineBoxes(numBoxes);
      Vector<int> fineProcAssign(numBoxes);

      for (int n=0; n<numBoxes; n++)
        {
          int halfDomain = fineDomainSize/2;
          int quarterDomain = fineDomainSize/4;
          if (n == 0)
            {
              fineBoxes[n] = Box(IntVect::Zero, (quarterDomain-1)*IntVect::Unit);
            }
          else if (n == 1)
            {
              fineBoxes[n] = Box((3*quarterDomain)*IntVect::Unit,
                                 (fineDomainSize-1)*IntVect::Unit);
            }
          else if (n == 2)
            {
              // this box should fail the disjointness test
              //fineBoxes[1] = Box(IntVect(D_DECL(7,10,10)), 17*IntVect::Unit);
              fineBoxes[2] = Box((halfDomain+2*nRef)*IntVect::Unit,
                                 (3*quarterDomain -1)*IntVect::Unit);
            }
        } // end loop over boxes for setup.

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

      interpolator.interpToFine(fineData, crseData);

      // now fill in ghost cells
      PiecewiseLinearFillPatch filler(fineGrids, crseGrids, 1, baseDomain,
                                      nRef, numGhost);


      filler.fillInterp(fineData, crseData, crseData, 1.0, 0, 0, 1);

      // finally, do an exchange
      fineData.exchange();

      // now test to see if we got the right answer:
      Real tolerance = 1.0e-10;
      Real fineDx = dx/nRef;

      DataIterator dit = fineGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& thisData = fineData[dit];
          // note that this will loop over ghost cells as well
          BoxIterator bit(thisData.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc(iv);
              loc += half*RealVect::Unit;
              loc *= fineDx;
              Real exactVal = getSolutionVal(loc);
              Real error = exactVal - thisData(iv,0);
              if (abs(error) > tolerance)
                {
                  // check to see if this is an interior or ghost cell
                  if (fineGrids[dit].contains(iv))
                    {
                      // interior cell
                      if (verbose)
                        {
                          pout() << " failed doubly-periodic fineInterp test at"
                                 << iv << ", exactVal = " << exactVal
                                 << ", error = " << error <<endl;
                        }
                      status += 1;
                    }
                  else
                    {
                      // ghost cell
                      if (verbose)
                        {
                          pout() << " failed doubly-periodic PWL test at"
                                 << iv << ", exactVal = " << exactVal
                                 << ", error = " << error <<endl;
                        }
                      status += 1000;
                    }
                } // end if we fail error tolerance test
            } // end loop over cells in box
        } // end loop over grids

    }

#if 0
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

      Real dx = 1.0/(baseDomainSize);
      initData(crseData, dx);

      // try multiple boxes for fine level
      Vector<Box> fineBoxes(3);
      Vector<int> fineProcAssign(3);

      fineBoxes[0] = Box(IntVect::Zero, 9*IntVect::Unit);
      fineBoxes[1] = Box(IntVect(D_DECL(6,10,10)), 15*IntVect::Unit);
      // this box should fail the disjointness test
      //fineBoxes[1] = Box(IntVect(D_DECL(7,10,10)), 17*IntVect::Unit);
      fineBoxes[2] = Box(IntVect(D_DECL(12,0,0)), IntVect(D_DECL(15,7,7)));

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

      interpolator.interpToFine(fineData, crseData);

      // now fill in ghost cells
      PiecewiseLinearFillPatch filler(fineGrids, crseGrids, 1, baseDomain,
                                      nRef, numGhost);


      filler.fillInterp(fineData, crseData, crseData, 1.0, 0, 0, 1);

      fineData.exchange(fineData.interval());

    }
#endif
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
          const IntVect& iv = boxIt();
          RealVect loc(iv);
          loc += half*RealVect::Unit;
          loc *= a_dx;
          localVal = getSolutionVal(loc);
          localData(iv) = localVal;
        }

    }
}
