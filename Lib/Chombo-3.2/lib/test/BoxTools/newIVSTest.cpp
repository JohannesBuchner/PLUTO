#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include "IntVect.H"
#include "Box.H"
#include "IntVectSet.H"
#include "DebugOut.H"
#include "UsingNamespace.H"

/// Global variables for test output:

static const char *pgmname = "newIVSTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

/**
   testIVS returns:
    0: all tests passed
 */
int testIVS(void);

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  bool passed = true;
  int icode = 0;

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  icode = testIVS();
  if (icode != 0)
    {
      pout() << indent << pgmname
           << ": failed with error code " << icode
           << endl;
      passed = false;
    }

  if (passed)
    pout() << indent << pgmname
           << ": passed all tests"
           << endl;
#ifdef CH_MPI
  MPI_Finalize();
#endif

  return icode;
}

int testIVS()
{
  int eekflag = 0;
#if (CH_SPACEDIM == 1)
  IntVect lo(0);
  IntVect hi(1);
  IntVect low(0);
  IntVect mid(2);
  IntVect high(4);
#elif CH_SPACEDIM==2
  IntVect lo(0, 0);
  IntVect hi(0, 1);
  IntVect low(0, 0);
  IntVect mid(0, 2);
  IntVect high(0, 4);
#elif CH_SPACEDIM==3
  IntVect lo(0, 0, 0);
  IntVect hi(0, 1, 0);
  IntVect low(0, 0, 0);
  IntVect mid(0, 2, 0);
  IntVect high(0, 4, 0);
#elif CH_SPACEDIM==4
  IntVect lo(0, 0, 0, 0);
  IntVect hi(0, 1, 0, 1);
  IntVect low(0, 0, 0, 0);
  IntVect mid(0, 2, 0, 2);
  IntVect high(0, 4, 0, 4);
#elif CH_SPACEDIM==5
  IntVect lo(0, 0, 0, 0, 0);
  IntVect hi(0, 1, 0, 1, 0);
  IntVect low(0, 0, 0, 0, 0);
  IntVect mid(0, 2, 0, 2, 0);
  IntVect high(0, 4, 0, 4, 0);
#elif CH_SPACEDIM==6
  IntVect lo(0, 0, 0, 0, 0, 0);
  IntVect hi(0, 1, 0, 1, 0, 1);
  IntVect low(0, 0, 0, 0, 0, 0);
  IntVect mid(0, 2, 0, 2, 0, 2);
  IntVect high(0, 4, 0, 4, 0, 4);
#else
  bogus_ch_spacedim();
#endif

  Box box(lo,hi);
  IntVectSet ivs0(box);
  ivs0.convert();
  for (IVSIterator ivsit(ivs0); ivsit.ok(); ++ivsit)
    {
      if (!box.contains(ivsit()))
        {
          eekflag = 4;
          return eekflag;
        }
    }
  IntVectSet ivs1(Box(low, high));
  ivs1 -= Box(mid, high);

  IntVectSet ivs2(ivs1);
  for (IVSIterator ivsit(ivs2); ivsit.ok(); ++ivsit)
    {
      if (!ivs1.contains(ivsit()))
        {
          eekflag = 2;
          return eekflag;
        }
    }

  if (verbose)
    {
      pout() << "testIVS" << endl;
      pout() << "ivs1:" << endl;
      dumpIVS(&ivs1);
      pout() << "ivs2:" << endl;
      dumpIVS(&ivs2);
    }
  IntVectSet ivs3 = ivs1 | ivs2;
  for (IVSIterator ivsit(ivs2); ivsit.ok(); ++ivsit)
    {
      if (!ivs1.contains(ivsit()))
        {
          eekflag = 6;
          return eekflag;
        }
    }

  if (verbose)
    {
      pout() << "ivs1 union with ivs2" << endl;
      dumpIVS(&ivs3);
    }

  for (IVSIterator ivsit(ivs3); ivsit.ok(); ++ivsit)
    {
      if (!ivs1.contains(ivsit()))
        {
          eekflag = 1;
          return eekflag;
        }
    }

  return eekflag;
}
/*****************/
///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
/*****************/
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
