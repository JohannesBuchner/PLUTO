#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Tues, Oct 5 1999

#include <cmath>
#include "Misc.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "UsingNamespace.H"

/// Global variables for test output:

static const char *pgmname = "copyTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

/**
   copyTest returns:
    0: all tests passed
 */
extern int
copyTest(void);

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

  icode = copyTest();
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
copyTest()
{
  IntVect ivlo = IntVect::Zero;
  IntVect ivhi = 7*IntVect::Unit;
  Box box(ivlo, ivhi);
  int ncomp = 5;
  FArrayBox fab(box, ncomp);
  Vector<Real> exactVal(ncomp);
  for (int i = 0; i < ncomp; i++)
    {
      exactVal[i] = Real(i) + 1.2345;
      fab.setVal(exactVal[i],i);
    }
  int dstvar = 3;
  int srcvar = 2;

  //now copy from one comp to another
  Interval srcI(srcvar, srcvar);
  Interval dstI(dstvar, dstvar);
  fab.copy( fab.box(), dstI, fab.box(), fab, srcI );

  int icode = 0;
  Real tol = 1.0e-10;
  for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      Real fabVal = fab(bit(), dstvar);
      if (Abs(fabVal-exactVal[srcvar]) > tol)
        {
          icode = 1;
          return(icode);
        }
    }


  return(icode);
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
