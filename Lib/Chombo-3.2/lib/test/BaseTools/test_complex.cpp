#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/// This tests how std::complex is implemented to make it easier to
/// figure out how to write code for Chombo's Complex::re() and ::im()
/// functions.  Since these need to be as fast as possible, they are
/// compiler dependent and can be whatever ugly hack gets the job done.
/// D.B.Serafini   Nov2004

// Usage:
//  <program-name> [-hqv] ...
//
//  where:
//    -h means print the help message and exit
//    -q means run quietly (only pass/fail messages printed)
//    -v means run verbosely (all messages printed)
//    ... all non-option arguments are ignored (Chombo convention)
//
//  Default is `-v'
//
//  Unknown options are treated as errors.

#include <iostream>
using std::endl;
#include "parstream.H"
#include "CH_Complex.H"
#include <cstring> //for strncmp()

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "MayDay.H"
#include "UsingBaseNamespace.H"

/// Global variables for handling output

static const char *pgmname = "testComplex" ;
static const char *indent = "   " ,*indent2 = "      " ;
static const char *usage = "[-hqv]" ;
static bool verbose = true ;

/// functions
extern "C" void check_complex_datarep( double * x );
void parseTestOptions( int argc ,char* argv[] );
void print_usage()
{
  pout() << "usage: " << pgmname << " " << usage << endl;
  MayDay::Abort();
}

////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int status = 0 ;
  parseTestOptions( argc ,argv ) ;
  if ( verbose ) pout() << indent << "Beginning " << pgmname << " ..." << endl ;

  std::complex<double> c(123.4567,765.4321);
  double r1,i1,r2,i2 ;
  r1 = c.real() ; i1 = c.imag() ;
  if ( verbose ) pout() << indent2 << "before: c = " << c << endl ;

  // see if the data is stored in the expected way
  check_complex_datarep( (double*)&c );

  size_t csize = sizeof(c) ,dsize=2*sizeof(r1) ;
  if ( csize != dsize )
  {
    pout() << "error: " << pgmname
         << ": sizeof std::complex [" << csize << "] is not twice the sizeof double ["
         << dsize << "]" << endl;
    ++status;
  }

  if ( verbose ) pout() << indent2 << "after: c = " << c << endl ;
  r2 = c.real() ; i2 = c.imag() ;
  if ( r1 != i2 || i1 != r2 )
  {
    pout() << "error: " << pgmname
         << ": the real,imag parts of std::complex do not have the usual format" << endl;
    ++status ;
  }

  Complex nullConstruct;
  Complex compConstruct = Complex(1.234, 5.678);

  pout() << indent << pgmname << ": "
       << ( (status == 0) ? "passed all tests," : "failed at least one test,")
       << endl;

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status ;
}

////////////////////////////////////////////////////////////////

// see if the complex is stored in memory in two consecutive
// words, real part first
void check_complex_datarep( double * x )
{
  double tmp = x[0] ;
  x[0] = x[1] ;
  x[1] = tmp ;
}

///
// Parse the standard test options (-v -q) out of the command line.
// The non-standard option -w means write an output file instead of reading input.
// Stop parsing when a non-option argument is found.
///
void parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
  {
    if ( argv[i][0] == '-' )
    {
      //if it is an option
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
      else if ( strncmp( argv[i] ,"-h" ,3 ) == 0 )
      {
        print_usage();
      }
      else
      {
        break ;
      }
    }
  }
  return ;
}
