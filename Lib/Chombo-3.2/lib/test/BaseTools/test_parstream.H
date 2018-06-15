#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TEST_PARSTREAM_H_
#define _TEST_PARSTREAM_H_

// -*- Mode: C++; Modified: "Tue 14 Mar 2006 16:01:09 by dbs"; -*-

///  Test harness for unit and regression tests for a class.

/**
 * To use this harness, copy the skeleton file "test_harness.skel"
 * into a C++ source file and follow the instructions therein.
 *
 * The class to be tested should have the declaration:\\
 *
 *    friend void Test::TestCLASSNAME(int argc, char * argv[]);\\
 *
 * where CLASSNAME is the name of the class to test.
 *
 * The test code file should also define the function:
 *
 *    void Test::usage();\\
 *
 * to print the standard usage message.  The harness handles calling this
 * function when the "-h" command line option is specified, and terminates
 * the program with exit code = 101.  The skeleton file provides an example
 * implementation of this function.
 *
 * The implementation of Test::TestCLASSNAME should use code like: \\
 *
 *      Test("string to describe test");\\
 *      ... code to do the test ...\\
 *      Fail(boolean_value_true_if_test_failed);\\
 *
 * for each test.  Alternatively, the Fail() function can take an optional
 * 2nd argument, a bool, which allows multiple conditions to be checked for
 * a single Test(), producing an appropriate message,  and is used as follows:
 *
 *      Test("string to describe test");\\
 *      ... code to do the test ...\\
 *      Fail(boolean_value_true_if_test_failed, false); ... \\
 *      Fail(another_boolean_value_true_if_test_failed, false); ...\\
 *      Fail(last_boolean_value_true_if_test_failed, true);\\
 *
 * The harness contains the code for the main program.  This program will
 * parse the standard command line options, run Test::TestCLASSNAME(), count
 * the number of successful and failed tests, report results at the end and
 * exit with the exit code set to the percentage of failures as an int
 * (i.e. 0--100).
 *
 * The Test::TestCLASSNAME() function may use the following global variables:\\
 *   bool debug    //if true, print debug output        (default:false)\\
 *   bool quiet    //if true, don't print normal output (default:false)\\
 *   bool verbose  //if true, print extra output        (default:false)\\
 *   const char *Name  //base name of the executable program from argv[0]
 *   Real epsilon  //threshold for testing floating pt. values (def:1e-6)\\
 *
 * and functions:\\
 *   Real reldiff(a,b) //returns |(a-b)/max(b)| unless a or b =0, IWC returns 1.
 *   Real absdiff(a,b) //returns |a-b|
 *
 * The variables are set by the command-line options.  See "Usage" below.
 *
 * The Test::TestCLASSNAME function may reparse the command line options
 * if it uses other options in addition to those supported by the harness.
 *
 * Usage:
 *   The main program handles several standard command-line options.\\
 *     -d       runs the program with debug output enabled (sets 'debug' var)\\
 *     -h       print the help message and exit (Test::Usage() function)\\
 *     -q       runs the program quietly (sets 'quiet' global var)\\
 *     -v       runs the program verbosely (sets 'verbose' global var)\\
 *     -E <eps> set threshold for flt.pt. tests ('epsilon' global var)\\
 */

#ifdef CH_MPI
#define MAKE_USAGE static void usage_function(void) \
  {                                                 \
  Test::Usage();                                    \
  MPI_Finalize();                                   \
  exit(101);                                        \
  }
#else
#define MAKE_USAGE static void usage_function(void) \
  {                                                 \
  Test::Usage();                                    \
  exit(101);                                        \
  }
#endif

// This is the macro that the test code file should use.

#define MAKE_TEST(CLASSNAME) \
namespace Test                                         \
  {                                                    \
    void Test ## CLASSNAME(int argc, char *argv[]);  \
    void Usage(void);                                \
  }                                                    \
  static void test_function(int argc, char *argv[])  \
  {                                                    \
    cout << " Testing class " << # CLASSNAME << endl; \
    Test::Test ## CLASSNAME(argc, argv);            \
  }                                                    \
MAKE_USAGE

#include <iostream>
#include <unistd.h>     //for getopt()
#include <stdlib.h>     //for exit(), atof()
#include <string.h>     //for strrchr()
#include "REAL.H"       //for Real
#include "SPMD.H"       //for parallel stuff
//#include "DebugOut.H"   //for view*() debugger routines
#include "UsingBaseNamespace.H"

using std::cout;
using std::cerr;
using std::endl;

/// macro for the test harness
#define Test(s)                                                       \
  ++tcount;                                                           \
  if (!quiet)                                                         \
  {                                                                   \
    cout << "  testing " << s << " ... ";                             \
  }                                                                   \
  if (tcount1>0)                                                    \
  {                                                                   \
    cout << endl << "WARNING: previous Test() not finished." << endl; \
    tcount1=0;                                                        \
  }

/// prototypes for the test harness functions
void AbortAll(bool b);            //if b, abort program
void Fail(bool b);               //if b, test failed
void Fail(bool b, bool last);     //multiple conditions; if any b, test failed
static void test_function(int argc, char* argv[]); //wrapper for test-specific code
static void usage_function(void);                  //<ditto>
static const char * ch_basename(const char * const argv0);

Real reldiff(Real a, Real b)
{
  if (a == b)
  {
    return (Real)0;
  }

  if (a == (Real)0 || b == (Real)0)
  {
    return (Real)1;
  }

  if (a > b)
  {
    return (a-b)/a;
  }

  return (b-a)/b;
}

Real absdiff(Real a, Real b)
{
  return a>b ? a-b : b-a;
}
////////////////////////////////////////////////////////////////

// global vars for test harness

static int tcount=0,         //number of tests
           tpassed=0,        //number of tests passed
           tfailed=0,        //number of tests failed
           tcount1=0;        //number of tests in current set

Real epsilon=0.000001;    //error threshold for flt.pt. tests

bool quiet=false;   //if true, don't print pass/fail msgs (-q option)
bool debug=false;   //if true, print debug output (-d)
bool verbose=false; //if true, print extra output (-v)
const char *Name;   //name of this program as called (without path)


////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  opterr = 0; //disable stderr msgs from getopt()
#ifdef CH_MPI
  MPI_Init(&argc ,&argv);
  {
#endif

    // set the program name variable
    Name = ch_basename(argv[0]);

    // Parse and remove standard command line options [dhqvE],
    // leave the rest in argv and adjust argc accordingly
    {
#define ARG_IS(str) strncmp(argv[ic], str ,strlen(str)+1) == 0
#define DBG                                                           \
      if (debug)                                                      \
      {                                                               \
        cout << "debug: main: option parser found canonical option [" \
             << argv[ic] << "]" <<endl;                               \
      }

#define DBGA                                                          \
      if (debug)                                                      \
      {                                                               \
        cout << "debug: main: option parser found canonical option [" \
             << argv[ic] << "]=|" << argv[ic+1] << "|" <<endl;        \
      }

      int iic = 1;
      for (int ic = 1; ic < argc; ++ic)
      {
        if (ARG_IS("-d"))
        {
          DBG;
          debug = true;
          quiet = false;
        }
        else if (ARG_IS("-h"))
        {
          DBG;
          usage_function();
        }
        else if (ARG_IS("-q"))
        {
          DBG;
          quiet = true;
          debug = false;
          verbose = false;
        }
        else if (ARG_IS("-v"))
        {
          DBG;
          verbose = true;
          quiet = false;
        }
        else if (ARG_IS("-E"))
        {
          DBGA;
          epsilon = atof(argv[++ic]);
        }
        else
        {
          // not a canonical option
          if (debug || verbose)
            cout << "debug: main: option parser is skipping command-line option ["
                 << argv[ic] << "]" << endl;
          // move the option down, if necessary,
          // and move the target index up one
          if (iic < ic) argv[iic] = argv[ic];
          ++iic;
        }
      }
      // if any args were removed, adjust arg count
      if (iic < argc) argc = iic;

#undef ARG_IS
#undef DBG
#undef DBGA
    }

    // this runs the tests
    test_function(argc ,argv);

    // final summary message
    //[NOTE: this prints even with -q; make it go away by redirecting stdout to /dev/null.]
    cout << " Completed " << tcount << " tests, ";
    if (tfailed > 0)
    {
      cout << tfailed << " failed." << endl;
    }
    else
    {
      cout << "all passed." << endl;
    }

#ifdef CH_MPI
  }
  MPI_Finalize();
#endif

  if (tcount > 0) return 100*tfailed / tcount;
  else             return 0;
}


////////////////////////////////////////////////////////////////

// test harness functions

void AbortAll(bool b)
{
  Fail(b);
  if (b)
  {
    cout << " Aborting tests.  Failed " << tfailed << " of " << tcount << " tests." << endl;
    exit(100*tfailed / tcount);
  }
}

void Fail(bool b)
{
  if (b)
  {
    ++tfailed;
    if (!quiet) cout << "failed." << endl;
  }
  else
  {
    ++tpassed;
    if (!quiet) cout << "passed." << endl;
  }
  return;
}

void Fail(bool b, bool last)
{
  static int tfailed1=0;
  ++tcount1;
  if (b)
  {
    ++tfailed;
    ++tfailed1;
  }
  else
  {
    ++tpassed;
  }
  if (last)
  {
    if (tfailed1 > 0)
    {
      if (!quiet) cout << "failed " << tfailed1 << " of " << tcount1 << "." << endl;
    }
    else
    {
      if (!quiet) cout << tcount1 << " passed." << endl;
    }
    tcount1 = 0;
    tfailed1 = 0;
  }
  else
  {
    // increment total test count only if this isn't a singleton test
    ++tcount;
  }
  return;
}


static const char * ch_basename(const char * const argv0)
{
#ifdef _WIN32
  char * base = argv0;
  char * c;
  for (c = argv0; *c != '\0'; ++c)
  {
    if (*c == '/') base = c + 1;
    if (*c == '\\') base = c + 1;
  }
  return base;
#else
  // Get the part of filename after the directory path.
  // If no '/', then return the whole string, else return
  // the part after the slash.
  const char * name = strrchr(argv0, '/');
  return (name == NULL) ? argv0 : ++name;
}
#endif

#endif
