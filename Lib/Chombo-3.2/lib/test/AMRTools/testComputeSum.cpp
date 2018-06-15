#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <vector>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "computeSum.H"
#include "FABView.H"
#include "DebugDump.H"

#include "UsingNamespace.H"


/// Prototypes:
int
testComputeSum();

void
parseTestOptions(int argc ,char* argv[]) ;

void
initData(LevelData<FArrayBox>& data,
         Real dx);

/// Global variables for handling output:
static const char* pgmname = "testComputeSum" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;
Real tolerance = 1.0e-10;

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testComputeSum() ;
  pout() << indent << pgmname ;
  if ( status == 0 )
    {
      pout() << " passed." << endl;
    }
  else
    {
      pout() << " failed with result code " << status << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status ;
}

int
testComputeSum()
{
  int status = 0;


  int numCells = 32;
  Real domainSize = 2.0;
  IntVect loVect = IntVect::Zero;
  IntVect hiVect = (numCells-1)*IntVect::Unit;
  Box domainBox(loVect, hiVect);
  ProblemDomain baseDomain(domainBox);

  int maxBoxSize = numCells/2;

  Vector<Box> vectBoxes;
  domainSplit(baseDomain, vectBoxes, maxBoxSize, 1);
  Vector<int> procAssign(vectBoxes.size(), 0);
  LoadBalance(procAssign, vectBoxes);
  DisjointBoxLayout level0Grids(vectBoxes, procAssign, baseDomain);
  Real dx0 = domainSize/numCells;


  // first test -- single level
  {
    int numLevels = 1;
    Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
    IntVect ghostVect = IntVect::Unit;
    phi[0] = new LevelData<FArrayBox>(level0Grids, 1, ghostVect);

    initData(*phi[0], dx0);

    Vector<int> nRef(numLevels, 4);
    Real sum = computeSum(phi, nRef, dx0);

    Real exactVal = D_TERM6(domainSize, *domainSize, *domainSize,
                            *domainSize, *domainSize, *domainSize);

    if (abs(sum - exactVal) > tolerance)
      {
        if (verbose)
          {
            pout() << "Single-level computeSum: expected" << exactVal
                   << ", computed sum = " << sum << endl;
          }
        status += 1;
      }

    // clean up storage
    delete phi[0];
    phi[0] = NULL;

  }



  // second test -- single level with maxLevel > 0
  {
    int numLevels = 3;
    Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
    IntVect ghostVect = IntVect::Unit;
    // levels 1 and 2 are undefined
    phi[0] = new LevelData<FArrayBox>(level0Grids, 1, ghostVect);

    initData(*phi[0], dx0);

    Vector<int> nRef(numLevels, 4);
    Real sum = computeSum(phi, nRef, dx0);

    Real exactVal = D_TERM6(domainSize, *domainSize, *domainSize,
                            *domainSize, *domainSize, *domainSize);

    if (abs(sum - exactVal) > tolerance)
      {
        if (verbose)
          {
            pout() << "single level, maxLevel > 0 computeSum: expected"
                   << exactVal
                   << ", computed sum = " << sum << endl;
          }
        status += 10;
      }

    // clean up storage
    delete phi[0];
    phi[0] = NULL;

  }


  // third test -- multiple levels
  {
    int numLevels = 3;
    Vector<int> nRef(numLevels, 4);

    Vector<LevelData<FArrayBox>* > phi(numLevels, NULL);
    IntVect ghostVect = IntVect::Unit;
    // level 0
    phi[0] = new LevelData<FArrayBox>(level0Grids, 1, ghostVect);

    initData(*phi[0], dx0);

    // level 1:
    {
      Vector<Box> level1Boxes;

      // do this in a bit of a silly way
      int quarterDomain = numCells / 4;
      int threeQuarterDomain = 3*numCells/4;
      IntVect loVect(quarterDomain);
      IntVect hiVect(IntVect::Unit*threeQuarterDomain);
      hiVect -= IntVect::Unit;
      Box spanBox(loVect, hiVect);
      int maxBox1 = quarterDomain;

      ProblemDomain level1Domain(baseDomain);
      level1Domain.refine(nRef[0]);

      ProblemDomain spanDomain(spanBox);
      domainSplit(spanDomain, level1Boxes, maxBox1, 1);
      Vector<int> level1ProcAssign(level1Boxes.size(), 0);
      LoadBalance(level1ProcAssign, level1Boxes);
      DisjointBoxLayout level1Grids(level1Boxes, level1ProcAssign,
                                    level1Domain);

      phi[1] = new LevelData<FArrayBox>(level1Grids, 1, IntVect::Unit);

      Real dx1 = dx0/nRef[0];

      initData(*phi[1], dx1);

    }

    Real sum = computeSum(phi, nRef, dx0);

    Real exactVal = D_TERM6(domainSize, *domainSize, *domainSize,
                            *domainSize, *domainSize, *domainSize);

    if (abs(sum - exactVal) > tolerance)
      {
        if (verbose)
          {
            pout() << "single level, maxLevel > 0 computeSum: expected"
                   << exactVal
                   << ", computed sum = " << sum << endl;
          }
        status += 10;
      }

    // clean up storage
    delete phi[0];
    phi[0] = NULL;

    delete phi[1];
    phi[1] = NULL;

  }

  return status;
}

void
initData(LevelData<FArrayBox>& a_data,
         Real a_dx)
{
  // for now, set phi to 1 everywhere
  Real phiVal = 1.0;

  DataIterator dit = a_data.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_data[dit].setVal(phiVal);
    }
}

///
// Parse the standard test options (-p -q) out of the command line.
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
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // argv[i] = "" ;
      }
    }
  }
  return ;
}
