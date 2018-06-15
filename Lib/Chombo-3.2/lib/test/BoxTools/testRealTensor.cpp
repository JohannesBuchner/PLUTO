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
#include <iostream>
using std::endl;

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "parstream.H"
#include "RealTensor.H"

#include "UsingNamespace.H"

void parseTestOptions( int argc ,char* argv[]);

int testRealTensor();

/// Global variables for handling output:
static const char *pgmname = "testRealTensor";
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

/// Code:

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv );

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl;

  ///
  // Run the tests
  ///
  int ret = testRealTensor();

  if ( ret == 0 )
    {
      pout() << indent << pgmname << " passed all tests" << endl;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

int testRealTensor()
{
  int errors = 0;

  // Check the determinant of the identity matrix
  {
    bool err = false;

#if CH_SPACEDIM == 1
    RealTensor ident(1.0);
#elif CH_SPACEDIM == 2
    RealTensor ident(1.0, 0.0,
                     0.0, 1.0);
#elif CH_SPACEDIM == 3
    RealTensor ident(1.0, 0.0, 0.0,
                     0.0, 1.0, 0.0,
                     0.0, 0.0, 1.0);
#elif CH_SPACEDIM == 4
    RealTensor ident(1.0, 0.0, 0.0, 0.0,
                     0.0, 1.0, 0.0, 0.0,
                     0.0, 0.0, 1.0, 0.0,
                     0.0, 0.0, 0.0, 1.0);
#elif CH_SPACEDIM == 5
    RealTensor ident(1.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 1.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 1.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 1.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 1.0);
#elif CH_SPACEDIM == 6
    RealTensor ident(1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
#endif

    Real det = ident.det();

    if (det != 1.0)
    {
      err = true;
    }

    if (err)
      {
        pout() << indent << pgmname
             << ": RealTensor failed the determinant of the identity test" << endl;
        pout() << indent2
             << "info: expected 1, "
             << "got " << det << endl;
        ++errors;
      }
    else
      {
        if (verbose)
          {
            pout() << indent << pgmname
               << ": RealTensor passed the determinant of the identity test" << endl;
          }
      }
  }

  return errors;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void parseTestOptions( int argc ,char* argv[])
{
  for (int i = 1; i < argc; ++i)
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i] ,"-v" ,3 ) == 0)
            {
              verbose = true;
            }
          else
          if (strncmp( argv[i] ,"-q" ,3 ) == 0)
            {
              verbose = false;
            }
          else
            {
              break;
            }
        }
    }
}
