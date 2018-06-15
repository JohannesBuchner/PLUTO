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
using std::endl ;
#include "parstream.H"
#include "RealVect.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testRealVect();

/// Global variables for handling output:
static const char *pgmname = "testRealVect" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testRealVect() ;

  if ( ret == 0 )
    {
      pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

#ifdef CH_USE_FLOAT
const Real eps = 0.75e-6 ;
#else
const Real eps = 1e-14 ;
#endif

int
testRealVect()
{
  int errors = 0 ;

  { // operator RealVect / Real
    bool err = false ;
    RealVect testrv( D_DECL6( 10.0 ,20.0 ,30.0, 40.0, 50.0, 60.0 ) );
    RealVect testdiv = testrv / 5.0 ;
    for (int dir=0; dir<SpaceDim; dir++)
      {
        Real relerror = Abs ( 1.0 - (testdiv[dir] / (2.0 * (Real)(dir+1) ) ) );
        err |= ( relerror > eps ) ;
      }

    if ( err )
      {
        pout() << indent << pgmname
             << ": RealVect failed the operator/(RV,R) test" << endl;
        pout() << indent2
             << "info: expected (" D_TERM6( << 2.0 ,<< 4.0 ,<< 6.0,
                                            << 8.0 ,<< 10.0 ,<< 12.0 )
             << "), got " << testdiv << endl;
        ++errors ;
      }
    else
      {
        if ( verbose )
          {
            pout() << indent << pgmname
                 << ": RealVect passed the operator/(RV,R) test" << endl;
          }
      }
  }

  { // operator Real / RealVect
    bool err = false ;
    RealVect testrv( D_DECL6( 10.0 ,20.0 ,30.0, 40.0, 50.0, 60.0 ) );
    RealVect testdiv = 5.0 / testrv ;
    for (int dir=0; dir<SpaceDim; dir++)
      {
        Real relerror = Abs ( 1.0 - ( testdiv[dir] * 2.0 * (Real)( dir+1 ) ) );
        err |= ( relerror > eps ) ;
      }
    if ( err )
      {
        pout() << indent << pgmname
             << ": RealVect failed the operator/(R,RV) test" << endl;
        pout() << indent2
             << "info: expected (" D_TERM6( << 0.5 ,<< 0.25 ,<< 1.0/6.0,
                                            << 0.125 ,<< 0.1 ,<< 1.0/12.0 )
             << "), got " << testdiv << endl;
        ++errors ;
      }
    else
      {
        if ( verbose )
        {
          pout() << indent << pgmname
               << ": RealVect passed the operator/(R,RV) test" << endl;
        }
      }
  }

  { // operator RealVect / RealVect
    bool err = false ;
    RealVect testrv1( D_DECL6( 10.0, 20.0, 30.0, 40.0, 50.0, 60.0 ) );
    RealVect testrv2( D_DECL6(  5.0 , 5.0 , 5.0, 5.0, 5.0, 5.0 ) );
    const RealVect testdiv = testrv1 / testrv2 ;
    for (int dir=0; dir<SpaceDim; dir++)
      {
        Real relerror = Abs (1.0 - (testdiv[dir] / (2.0 * (Real)(dir+1) ) ) );
        err |= ( relerror > eps ) ;
      }
    if ( err )
      {
        pout() << indent << pgmname
             << ": RealVect failed the operator/(RV,RV) test" << endl;
        pout() << indent2
             << "info: expected (" D_TERM6( << 2.0 ,<< 4.0 ,<< 6.0,
                                            << 8.0 ,<< 10.0 ,<< 12.0 )
             << "), got " << testdiv << endl;
        ++errors ;
      }
    else
      {
        if ( verbose )
        {
          pout() << indent << pgmname
               << ": RealVect passed the operator/(RV,RV) test" << endl;
        }
      }
  }

  return errors ;
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
