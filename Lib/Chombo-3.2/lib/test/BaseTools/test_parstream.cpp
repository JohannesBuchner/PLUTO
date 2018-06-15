#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/// Test code for Chombo/BoxTools/parstream.*

// There are 3 functions to test:
//    pout(), setPoutBaseName() and poutFileName().
//
// It's difficult to test pout() in serial since it should always write to
// std::cout.  Too bad a program can't redirect its own output.  In
// parallel, pout() writes to a file so it can be tested.
//
// The other 2 functions relate to parallel filenames and have should
// do nothing significant in serial, but the insignificant behavior
// can be tested.
//
// Summary of testing:
//  pout() tests:
//   1) pout() stream is in good state
//   2) output to pout() doesn't abort
//   3) pout() stream is in good state after setPoutBaseName()
//   4) output to pout() after setPoutBaseName() doesn't abort
//   5) pout() stream is in good state after poutFileName()
//   6) output to pout() after poutFileName() doesn't abort
// [serial-only tests]
//   7) poutFileName() returns string "cout"
//   8) pout() is std::cout
//   9) pout() is std::cout after setPoutBaseName("foo")
//  10) pout() is std::cout after setPoutBaseName("cout")
//  11) pout() is std::cout after poutFileName()
// [parallel-only tests]
//   8) poutFileName() returns default of "pout.<procnum>" on
//       proc <procnum>
//   9) output to pout() writes to "pout.<procnum>" file on proc <procnum>
//  10) setPoutBaseName("pout") doesn't abort
//  11) output to pout() after setPoutBaseName("pout") doesn't abort
//      and appends to pout.<procnum> on proc <procnum>
//  12) pout() is in good state after output
//  13) poutFileName() after setPoutBaseName("pout") returns
//      "pout.<procnum>" on proc <procnum>
//  14) output to pout() after poutFileName() doesn't abort and
//      appends to pout.<procnum> on proc <procnum>
//  15) setPoutBaseName("test") doesn't abort
//  16) output to pout() after setPoutBaseName("test") writes to
//      test.<procnum> on proc <procnum> and doesn't change
//      pout.<procnum>
//  17) poutFileName() after setPoutBaseName("test") returns
//      test.<procnum> on proc <procnum>
//  18) setPoutBaseName("test2") doesn't abort
//  19) poutFileName() returns "test2.<procnum> on proc <procnum>
//  20) output to pout() doesn't abort and leaves stream in good state
//      and writes to "test2.<procnum>" and doesn't write to "test.<procnum>"
//      or to "pout.<procnum>
//
// Usage:
//   The main program handles several standard command-line options.
//     -d       runs the program with debug output enabled (sets 'debug' var)
//     -h       print the help message and exit ('usage' global var)
//     -q       runs the program quietly (sets 'quiet'  var)
//     -v       runs the program verbosely (sets 'verbose' global var)
//     -E <eps> set threshold for flt.pt. tests ('epsilon' global var)
//
// History:
//  13Mar06 <dbs> initial design and coding based on test_harness.skel
//  <dbs> == David Serafini, LBL/ANAG <dbs@hpcrdm.lbl.gov>
//
//////////////////////////////////////////////////////////////////

// header file for this test
#include "test_parstream.H"

// header file for the class to test
#include "parstream.H"
#include <fstream>
#include "UsingBaseNamespace.H"

using std::ifstream;

// This generates the main program to test the class.
MAKE_TEST( pout )

////////////////////////////////////////////////////////////////

// some global variables for convenience

#ifdef CH_MPI
static std::string DefaultBaseName         = "pout" ;
static std::string ExpectedDefaultFileName = "pout.0" ;
#else
static std::string ExpectedDefaultFileName = "cout" ;
#endif
static std::string TestBaseName1 = "test1" ;
static std::string TestBaseName2 = "test2" ;
static std::string TestLine1 = "***_this_is_test_line_1_***" ;
static std::string TestLine2 = "***_this_is_test_line_2_***" ;

// This function does the actual testing.
//[NOTE: MPI_Initialize() has already been called so we can't test the
//       behaviour of the pout functions in the before-MPI-initialized case.]

void Test::Testpout( int argc, char *argv[] )
{
  { // some tests for the accessor and manipulator functions before
    // pout() is invoked.
    Test( "accessors and manipulators" );
    std::string default_filename = poutFileName() ;
    if ( procID() == 0 )
      Fail( default_filename != ExpectedDefaultFileName ,false );
    setPoutBaseName( TestBaseName1 );
    std::string filename_1 = poutFileName() ;
#ifdef CH_MPI
    // other procs besides 0 and 1 have different filenames, but don't bother checking them
    if ( procID() <= 1 )
    {
      std::string expected_filename = TestBaseName1 ;
      expected_filename += ( procID() == 0 ) ? ".0" : ".1" ;
      Fail( filename_1 != expected_filename ,true );
      if ( verbose && filename_1 != expected_filename )
      {
        std::cout << std::endl << "fail[p" << procID() << "] getting pout filename: expected ["
                  << expected_filename << "], got [" << filename_1 << "]" << std::endl;
      }
    }
    else
    {
      // only need to know if proc==0,1 fails
      Fail( false , true );
    }
#else
    // in serial, the filename never changes
    Fail( filename_1 != ExpectedDefaultFileName ,true );
    if ( verbose && filename_1 != ExpectedDefaultFileName )
    {
      std::cout << std::endl << "fail: getting pout filename: expected ["
                << ExpectedDefaultFileName << "], got [" << filename_1 << "]" << std::endl;
    }
#endif
  }

  { // test that the stream is ok
#ifndef CH_MPI
    Test( "serial pout() is std::cout" );
    Fail( pout() != std::cout );
#endif
    Test( "that pout() is initialized in good state" );
    Fail( !pout() ,false );
    Fail( ! pout().good() ,false );
    Fail( pout().fail() ,false );
    Fail( pout().bad() ,true );
  }

  { // test that stream works
    Test( "that pout() works" );
    pout() << TestLine1 << std::endl ;
    // same tests as above -- not proof that the write worked, but better than nothing
    Fail( !pout() ,false );
    Fail( ! pout().good() ,false );
    Fail( pout().fail() ,false );
    Fail( pout().bad() ,false );
    // extra test
    Fail( pout().eof() ,true );
  }
#ifdef CH_MPI
  { // extra tests for parallel: check the file contents
    std::string filename1 = poutFileName() ;
    pout().flush() ;
    ifstream pin( filename1.c_str() );
    Test( "that the pout file has the right name and contents" );
    Fail( !pin ,false );
    std::string line1 ; pin >> line1 ;
    Fail( line1 != TestLine1 ,false );
    if ( verbose && line1 != TestLine1 )
    {
      std::cout << std::endl << "fail[p" << procID() << "] reading pout file: expected [" << TestLine1
                << "], got [" << line1 << "]" << std::endl ;
    }
    // that should be the only data in the file
    std::string line2; pin >> line2 ;
    Fail( ! pin.eof() ,true );
    if ( verbose && ! pin.eof() )
    {
      std::cout << std::endl << "fail[p" << procID() << "] expected eof() but didn't get it.  Next line is ["
                << line2 << "]" << std::endl ;
    }
  }
  { // test that setting a new name writes to a new file and doesn't change the old file
    Test( "changing the pout base name" );
    setPoutBaseName( TestBaseName2 );
    std::string filename2 = poutFileName() ;
    if ( debug )
    {
      std::cout << std::endl << "dbg: new filename is [" << filename2 << "] ... " ;
    }
    int compare = strncmp( filename2.c_str() ,TestBaseName2.c_str() ,TestBaseName2.size() ) ;
    if ( debug )
    {
      std::cout << std::endl << "dbg: strcmp(filename,basename) returned " << compare << " ... " ;
    }
    Fail( compare != 0 ,false ) ;
    pout() << TestLine1 << std::endl ;
    pout() << TestLine2 << std::endl ;
    pout().flush() ;
    // switch filenames again to close the file
    setPoutBaseName( DefaultBaseName );
    pout() << DefaultBaseName << std::endl ;
    // read the previous output file
    ifstream pin( filename2.c_str() );
    Fail( !pin ,false );
    if ( verbose && !pin )
    {
      std::cout << "fail[p" << procID() << "] opening file [" << filename2 << "] for reading ... " ;
    }
    std::string line1,line2 ;
    pin >> line1 ;
    pin >> line2 ;
    if ( debug )
    {
      std::cout << std::endl << "dbg: line1,2 are [" << line1 << "],[" << line2 << "] ... " ;
    }
    Fail( line1 != TestLine1 ,false ) ;
    Fail( line2 != TestLine2 ,true ) ;
  }
#endif
}

// This function prints the "usage" message for the test program.
void Test::Usage( )
{
  std::cout << "usage: " << Name << " [-dhqv]" << std::endl ;
}
