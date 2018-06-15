#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// testParmParse.cpp

// tests ParmParse.  using the input file inputs.testParmParse you
// should get the output in 2d:
//
//   testing query functions
//   foo1 = 1
//   foo2 = 2.5
//   pp.countname(foo2_name_stds) = 2
//   foo3 = 3
//   foo4stdv[0] = 4
//   foo4stdv[1] = 5
//   foo5_stds = 'hello'
//   foo6stdvs[0] = 'man'
//   foo6stdvs[1] = 'plan'
//   foo7stdv[0] = 7
//   foo7stdv[1] = 8
//
//   testing get functions
//   foo1 = 1
//   foo2 = 2.5
//   pp.countname(foo2_name_stds) = 2
//   foo3 = 3
//   foo4stdv[0] = 4
//   foo4stdv[1] = 5
//   foo5_stds = 'hello'
//   foo6stdvs[0] = 'man'
//   foo6stdvs[1] = 'plan'
//   foo7stdv[0] = 7
//   foo7stdv[1] = 8
//
// and in 3d you should get:
//
//   testing query functions
//   foo1 = 1
//   foo2 = 2.5
//   pp.countname(foo2_name_stds) = 2
//   foo3 = 3
//   foo4stdv[0] = 4
//   foo4stdv[1] = 5
//   foo4stdv[2] = 6
//   foo5_stds = 'hello'
//   foo6stdvs[0] = 'man'
//   foo6stdvs[1] = 'plan'
//   foo6stdvs[2] = 'canal'
//   foo7stdv[0] = 7
//   foo7stdv[1] = 8
//   foo7stdv[2] = 9
//
//   testing get functions
//   foo1 = 1
//   foo2 = 2.5
//   pp.countname(foo2_name_stds) = 2
//   foo3 = 3
//   foo4stdv[0] = 4
//   foo4stdv[1] = 5
//   foo4stdv[2] = 6
//   foo5_stds = 'hello'
//   foo6stdvs[0] = 'man'
//   foo6stdvs[1] = 'plan'
//   foo6stdvs[2] = 'canal'
//   foo7stdv[0] = 7
//   foo7stdv[1] = 8
//   foo7stdv[2] = 9
//

#include <iostream>
#include <iomanip>
#include <string>

//#include "SPACE.H"
#include "REAL.H"
#include "parstream.H"
#include "ParmParse.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingBaseNamespace.H"

/// Prototypes:

using std::endl;
using std::boolalpha;

// since this is now in BaseTools, need to define our own SpaceDim

int SpaceDim = CH_SPACEDIM;

void
parseTestOptions(int argc ,char* argv[]) ;

int
testParmParse(int argc, char* argv[]);

void
testQuery();

void
testGet();

void
testPrefix();

/// Global variables for handling output:

static const char *pgmname = "testParmParse" ;
static const char *indent = "   " ,*indent2 = "      " ;
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

  int ret = testParmParse(argc, argv);
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

char dummy[1];

int
testParmParse(int argc, char* argv[])
{
  dummy[0] = 0;

  // use the first argument that isn't an option as the input filename
  // and take it out of the argument list
  const char* in_file = NULL;
  for ( int i=1 ; i<argc ; ++i )
    {
      if ( strlen(argv[i]) > 0 && argv[i][0] != '-' )
      {
        in_file = argv[i] ;
        argv[i] = dummy;
        break ;
      }
    }

  if ( in_file == NULL )
    {
      in_file = "inputs.testParmParse";
    }

  ParmParse pp(argc-1, argv+1, NULL, in_file);

  if (verbose)
  {
    pout() << indent2 << pgmname
         << ": testing query functions" << endl;
    pout() << indent2 << "ParmParse dumpTable() before testQuery():" << endl ;
    pp.dumpTable( pout() );
  }

  testQuery();

  if (verbose)
  {
    pout() << indent2 << pgmname
         << ": testing get functions" << endl;
    pout() << indent2 << "ParmParse dumpTable() before testGet():" << endl ;
    pp.dumpTable( pout() );
  }

  testGet();


  testPrefix();

  // if there were any errors, the program would have aborted
  pout() << indent << pgmname
       << ": ParmParse tests passed." << endl ;

  return 0 ;
}

void
testQuery()
{
    ParmParse pp;

    int foo1 = 0;
    pp.query("foo1", foo1);
    if (verbose) pout() << indent2 << pgmname << ": query int foo1 = " << foo1 << endl;

    Real foo2 = 0;
    pp.query("foo2", foo2);
    if (verbose) pout() << indent2 << pgmname << ": query Real foo2 = " << foo2 << endl;

    std::string foo2_name_stds("foo2");
    int num_foo2_stds = pp.countname(foo2_name_stds);
    if (verbose) pout() << indent2 << pgmname << ": std::string pp.countname(foo2) = " << num_foo2_stds << endl;

    Real foo3 = 0;
    pp.query("foo3", foo3);
    if (verbose) pout() << indent2 << pgmname << ": query Real foo3 = " << foo3 << endl;

    std::vector<Real> foo4stdv(SpaceDim);
    pp.queryarr("foo4", foo4stdv, 0, SpaceDim);
    if (verbose)
      {
        pout() << indent2 << pgmname << ": queryarr std::vector<Real> foo4 = " << endl ;
        for (int d = 0; d < SpaceDim; ++d)
          {
            pout() << indent2 << "[" << d << "] = " << foo4stdv[d] << endl;
          }
      }

    std::string foo5_stds;
    pp.query("foo5", foo5_stds);
    if (verbose) pout() << indent2 << pgmname << ": query std::string foo5 = '" << foo5_stds << "'" << endl;

    std::vector<std::string> foo6stdvs;
    pp.queryarr("foo6", foo6stdvs, 0, SpaceDim);
    if (verbose)
      {
        pout() << indent2 << pgmname << ": queryarr std::vector<std::string> foo6 =" << endl ;
        for (int d = 0; d < SpaceDim; ++d)
          {
            pout() << indent2 << "[" << d << "] = '" << foo6stdvs[d] << "'" << endl;
          }
      }

    std::vector<int> foo7stdv(SpaceDim);
    pp.queryarr("foo7", foo7stdv, 0, SpaceDim);
    if (verbose)
      {
        pout() << indent2 << pgmname << ": queryarr std::vector<int> foo7 =" << endl ;
        for (int d = 0; d < SpaceDim; ++d)
          {
            pout() << indent2 << "[" << d << "] = " << foo7stdv[d] << endl;
          }
      }

    pout() << boolalpha ; // output symbolic representation of boolean (true,false) instead of (1,0)
    bool foobool1=false,foobool2=false,foobool3=true,foobool4=true,foobool5=false,foobool6=true;
    pp.query("foobool1",foobool1);
    pp.query("foobool2",foobool2);
    pp.query("foobool3",foobool3);
    pp.query("foobool4",foobool4);
    pp.query("foobool5",foobool5);
    pp.query("foobool6",foobool6);
    if (verbose)
      {
        pout() << indent2 << pgmname << ": query foobool1: input=true,  expected result true,  result = " << foobool1 << endl ;
        pout() << indent2 << pgmname << ": query foobool2: input=1,     expected result true,  result = " << foobool2 << endl ;
        pout() << indent2 << pgmname << ": query foobool3: input=false, expected result false, result = " << foobool3 << endl ;
        pout() << indent2 << pgmname << ": query foobool4: input=0,     expected result false, result = " << foobool4 << endl ;
        pout() << indent2 << pgmname << ": query foobool5: input=TRUE,  expected result true,  result = " << foobool5 << endl ;
        pout() << indent2 << pgmname << ": query foobool6: input=FALSE, expected result false, result = " << foobool6 << endl ;
      }
    if ( !( foobool1 && foobool2 && ! foobool3 && ! foobool4 && foobool5 && ! foobool6 ) )
    {
      MayDay::Error( "bool query test failed" ,1 );
    }
}

void
testGet()
{
    ParmParse pp;

    int foo1 = 0;
    pp.get("foo1", foo1);
    if (verbose) pout() << indent2 << pgmname << ": get int foo1 = " << foo1 << endl;

    Real foo2 = 0;
    pp.get("foo2", foo2);
    if (verbose) pout() << indent2 << pgmname << ": get Real foo2 = " << foo2 << endl;

    std::string foo2_name_stds("foo2");
    int num_foo2_stds = pp.countname(foo2_name_stds);
    if (verbose) pout() << indent2 << pgmname << ": std::string pp.countname(foo2) = " << num_foo2_stds << endl;

    Real foo3 = 0;
    pp.get("foo3", foo3);
    if (verbose) pout() << indent2 << pgmname << ": get Real foo3 = " << foo3 << endl;

    std::vector<Real> foo4stdv(SpaceDim);
    pp.getarr("foo4", foo4stdv, 0, SpaceDim);
    if (verbose)
      {
        pout() << indent2 << pgmname << ": queryarr std::vector<Real> foo4 = " << endl ;
        for (int d = 0; d < SpaceDim; ++d)
          {
            pout() << indent2 << "[" << d << "] = " << foo4stdv[d] << endl;
          }
      }

    std::string foo5_stds;
    pp.get("foo5", foo5_stds);
    if (verbose) pout() << indent2 << pgmname << ": get std::string foo5 = '" << foo5_stds << "'" << endl;

    std::vector<std::string> foo6stdvs;
    pp.getarr("foo6", foo6stdvs, 0, SpaceDim);
    if (verbose)
      {
        pout() << indent2 << pgmname << ": getarr std::vector<std::string> foo6 =" << endl ;
        for (int d = 0; d < SpaceDim; ++d)
          {
            pout() << indent2 << "[" << d << "] = '" << foo6stdvs[d] << "'" << endl;
          }
      }

    std::vector<int> foo7stdv(SpaceDim);
    pp.getarr("foo7", foo7stdv, 0, SpaceDim);
    if (verbose)
      {
        pout() << indent2 << pgmname << ": getarr std::vector<int> foo7 =" << endl ;
        for (int d = 0; d < SpaceDim; ++d)
          {
            pout() << indent2 << "[" << d << "] = " << foo7stdv[d] << endl;
          }
      }
}


void
testPrefix()
{
    ParmParse ppsub("sub");

    const char* subPrefix = ppsub.prefix();
    std::string prefixString(subPrefix);

    Real subOne = 0;
    ppsub.query("sub1", subOne);


    Real subFooOne = 0;
    ppsub.query("foo1", subFooOne);

    if (verbose)
      {
        pout() << indent2 << pgmname << ": prefix() test:  expected result 'sub',  result = " << subPrefix << endl ;
        pout() << indent2 << pgmname << ": query sub.sub1:  expected result 5.0,  result = " << subOne << endl ;
        pout() << indent2 << pgmname << ": query sub.foo1:  expected result 100.0,  result = " << subFooOne << endl ;
      }
    if ( !( (prefixString == "sub") && (subOne == 5.0) && (subFooOne == 100.0) ) )
    {
      MayDay::Error( "Prefix test failed" ,1 );
    }
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
