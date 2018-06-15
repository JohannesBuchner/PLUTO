#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves, Monday, April 23, 2007

#include <cmath>
#include "Misc.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "TimedDataIterator.H"
#include "UsingNamespace.H"

/// Global variables for test output:

static const char *pgmname = "dblTimeTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

/**
   dblTimeTest returns:
    0: all tests passed
 */
extern int
dblTest(void);

/// Code:
int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  bool passed = true;
  int icode = 0;

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  icode = dblTest();
  if (icode != 0)
    {
      pout() << indent << pgmname
           << " failed with error code " << icode
           << endl;
      passed = false;
    }

  if (passed)
    {
      pout() << indent << pgmname
           << " passed all tests"
           << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}

int
dblTest()
{
  DisjointBoxLayout loLayoutFine;
#if (CH_SPACEDIM == 1)
  // no longer a bogus spacedim
  Box b1(IntVect(7), IntVect(7));
  Box b2(IntVect(3), IntVect(3));
  Box b3(IntVect(9), IntVect(9));
  Box b4(IntVect(-1), IntVect(-1));
#elif (CH_SPACEDIM == 2)
  Box b1(IntVect(7,6), IntVect(7,9));
  Box b2(IntVect(3,0), IntVect(3,3));
  Box b3(IntVect(9,4), IntVect(9,5));
  Box b4(IntVect(-1,4), IntVect(-1,12));
#elif (CH_SPACEDIM == 3)
  Box b1(IntVect(7,6,0), IntVect(7,9,0));
  Box b2(IntVect(3,0,0), IntVect(3,3,0));
  Box b3(IntVect(9,4,0), IntVect(9,5,0));
  Box b4(IntVect(-1,4,0), IntVect(-1,12,0));
#elif (CH_SPACEDIM == 4)
  Box b1(IntVect(7,6,0,7), IntVect(7,9,0,7));
  Box b2(IntVect(3,0,0,3), IntVect(3,3,0,3));
  Box b3(IntVect(9,4,0,9), IntVect(9,5,0,9));
  Box b4(IntVect(-1,4,0,-1), IntVect(-1,12,0,-1));
#elif (CH_SPACEDIM == 5)
  Box b1(IntVect(7,6,0,7,6), IntVect(7,9,0,7,9));
  Box b2(IntVect(3,0,0,3,0), IntVect(3,3,0,3,3));
  Box b3(IntVect(9,4,0,9,4), IntVect(9,5,0,9,5));
  Box b4(IntVect(-1,4,0,-1,4), IntVect(-1,12,0,-1,12));
#elif (CH_SPACEDIM == 6)
  Box b1(IntVect(7,6,0,7,6,0), IntVect(7,9,0,7,9,0));
  Box b2(IntVect(3,0,0,3,0,0), IntVect(3,3,0,3,3,0));
  Box b3(IntVect(9,4,0,9,4,0), IntVect(9,5,0,9,5,0));
  Box b4(IntVect(-1,4,0,-1,4,0), IntVect(-1,12,0,-1,12,0));
#else
  bogus spacedim;
#endif
  loLayoutFine.addBox(b1,0);
  loLayoutFine.addBox(b2,0);
  loLayoutFine.addBox(b3,0);
  loLayoutFine.addBox(b4,0);
  loLayoutFine.close();
  int nsleep = 3;
  int totsleepProc =0;
  LevelData<FArrayBox> loFabFine(loLayoutFine,1);
  TimedDataIterator ditFine =loFabFine.timedDataIterator();
  ditFine.enableTime();
  ditFine.clearTime();
  for (ditFine.reset(); ditFine.ok(); ++ditFine)
    {
      sleep(nsleep);
      totsleepProc += nsleep;
    }
  ditFine.disableTime();
  Vector<unsigned long long> allTime = ditFine.getTime();
  Real totTimeMeasured = 0;
  if (verbose)
    {
      for (int ivec = 0; ivec < allTime.size(); ivec++)
        {
          pout() << "allTime[" << ivec << "]=" << allTime[ivec] << endl;
          totTimeMeasured += allTime[ivec];
        }
    }
  Real eps = 0.01;
  if (((totTimeMeasured - totsleepProc)/totsleepProc) > eps)
    {
      return -1;
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
