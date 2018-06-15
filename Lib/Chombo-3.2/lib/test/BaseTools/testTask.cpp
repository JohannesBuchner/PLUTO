#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <cstring>

#include "parstream.H"
#include "SPMD.H"

#include "UsingBaseNamespace.H"

/*
 *  Test program for doing task-based parallelism with Chombo inside a global MPI communicator
 *  that is split.
 */


/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;


/// Global variables for handling output:
static const char *pgmname = "testTask" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;
static bool throwException = false;

/// Code:

int
main(int argc, char* argv[])
{
  int rank = 0;
  int size = 1;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &size);
  MPI_Comm localComm;
  MPI_Comm_split(Chombo_MPI::comm, rank, 0, &localComm);
  // split will be suprisingly slow
  Chombo_MPI::comm = localComm;
  char baseN[30];
  sprintf(baseN, "pout.%d", rank);
  setPoutBaseName(std::string(baseN));
#endif
  parseTestOptions( argc ,argv ) ;
  int err = 0;
  pout() << "my size is " << size << std::endl;
  pout() << "my original rank is " << rank << " my new rank is " << procID() << std::endl;

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl ;


  if ( verbose )
    {
      if (err == 0)
        {
          pout() << indent << "testTask test passed." << std::endl;
        }
      else
        {
          pout() << indent << "testTask FAILED!!! "<<err<<std::endl;
        }
    }

#ifdef CH_MPI
  MPI_Comm_free(&localComm);
  MPI_Finalize();
#endif

  return err;
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
      else if ( strncmp( argv[i],"-e", 3) == 0 )
        {
          throwException = true;
        }

      else
      {
        break ;
      }
    }
  }
  return ;
}
