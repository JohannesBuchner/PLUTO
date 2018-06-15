#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include <cstring>

#include "REAL.H"
#include "Vector.H"
#include "CH_Timer.H"
#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "MayDay.H"
#include "ProblemDomain.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "LevelData.H"

#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:
void
parseTestOptions( int argc ,char* argv[] );

int testExchange(void);

/// Global variables for handling output:
static const char *pgmname = "interiorExchangeTest";
static const char *indent2 = "      ";
static bool verbose = true;

/// Code:
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions(argc,argv);

  if ( verbose ) pout() << indent2 << "Beginning " << pgmname << " ..." << endl;

  ///
  // Run the tests
  ///
  int ret = testExchange();
  if (ret == 0)
    {
      pout() << "exchange test passed" << endl;
    }
  else
    {
      pout() << "exchange test failed with code" << ret << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return ret;
}


int testExchange(void)
{

  int domsize = 32;
  int maxbox = 16;
  Box bigBox = Box(IntVect::Zero, (domsize-1)*IntVect::Unit);
  ProblemDomain baseLevelDomain(bigBox);
  ProblemDomain dom2;
  dom2.define(baseLevelDomain);
  Vector<Box> baseLevelBoxes;
  domainSplit(baseLevelDomain, baseLevelBoxes, maxbox, maxbox);

  const int numLevels = 1;
  Vector<int> ranks;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      LoadBalance(ranks, baseLevelBoxes);
    }

  DisjointBoxLayout grids(baseLevelBoxes, ranks);
  int nghost = 4;
  LevelData< BaseFab<int> > data(grids, SpaceDim, nghost*IntVect::Unit);

  //set the data on the boxes according to which quadrant it is in
  int midpt = domsize/2;
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      for (BoxIterator bit(grids[dit()]); bit.ok(); ++bit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              int  val = 0;
              if (bit()[idir] < midpt)
                {
                  val = -1;
                }
              else
                {
                  val = 1;
                }
              data[dit()](bit(), idir) = val;
            }
        }
    }
  //exchange the data to fill the ghost cells
  data.exchange();
  //check the answer
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      Box grownBox = data[dit()].box();
      grownBox &= baseLevelDomain;
      for (BoxIterator bit(grownBox); bit.ok(); ++bit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              int  val = 0;
              if (bit()[idir] < midpt)
                {
                  val = -1;
                }
              else
                {
                  val = 1;
                }
              if (data[dit()](bit(), idir) != val)
                {
                  pout() << "value at " << bit() << "looks wrong" << endl;
                  return -42;
                }
            }
        }
    }
  return 0;
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
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
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
