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

static const char *pgmname = "ivsfabTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

/// Prototypes:

void parseTestOptions( int argc ,char* argv[] ) ;
int  sum(const IntVect & a_iv)
{
  return D_TERM6(     a_iv[0],
                 +  2*a_iv[1],
                 +  4*a_iv[2],
                 +  8*a_iv[3],
                 + 16*a_iv[4],
                 + 32*a_iv[5]);
}

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
  int errors = 0 ;
  Box box(IntVect::Zero, IntVect::Unit);
  IntVectSet  ivs(box);
  int ncomp = 3;
  IVSFAB<Real> fab(ivs, ncomp);

  fab.setVerbose(verbose);

  // test that the IVS is the same
  IntVectSet ivstest(fab.getIVS()) ;
  if ( ! ( ivs == ivstest ) )
  {
    if ( verbose ) pout() << "error: " << pgmname << ": failed IntVectSet consistency test" << endl;
    ++errors ;
  }
  // test that numIvs() returns the expected result
  if ( fab.numIvs() != box.volume() )
  {
    if ( verbose ) pout() << "error: " << pgmname << ": failed numIvs() test" << endl;
    ++errors ;
  }
  // test that setVal changes all the elements
  fab.setVal(7.5);
  for ( BoxIterator bit(box) ; bit.ok() ; ++bit )
  {
    for (int icomp = 0; icomp < ncomp; icomp++)
    {
      if ( fab(bit(),icomp) != 7.5 )
      {
        if ( verbose )
          pout() << "error: " << pgmname << ": value at (" << bit() << "," << icomp
               << ") is not what setVal() set" << endl;
        ++errors;
      }
    }
  }

  //[NOTE: no way to test the output of this, so it effectively only checks for abort.]
  if ( verbose ) dumpIVSFAB(&fab);

  // visual check that indexing+iterator doesnt fail
  Real val = 0.0;
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
  {
    for (int icomp = 0; icomp < ncomp; icomp++)
    {
          fab(ivsit(), icomp) = val;
          val += 1.0;
    }
  }
  //[NOTE: see above]
  if ( verbose ) dumpIVSFAB(&fab);

  // test that indexing produces the expected results
  // by putting values based on the IntVect into the fab
  // and taking them back out
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
  {
    for (int icomp = 0; icomp < ncomp; icomp++)
    {
      fab(ivsit(), icomp) = 10*sum(ivsit()) + icomp;
    }
  }
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
  {
    for (int icomp = 0; icomp < ncomp; icomp++)
    {
      val = 10*sum(ivsit()) + icomp;
      if ( fab(ivsit(), icomp) != val )
      {
        if ( verbose )
          pout() << "error: " << pgmname << ": value at (" << ivsit() << "," << icomp
               << ") is " << fab(ivsit(), icomp) << ") but should be " << val << endl ;
        ++errors;
      }
    }
  }

  return errors;
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
