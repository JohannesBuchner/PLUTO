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

#include "CH_Attach.H"
#include "MayDay.H"
#include "parstream.H"

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testRCP();

/// Global variables for handling output:
static const char *pgmname = "testCH_Attach" ;
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
    pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl ;

  int err = registerDebugger();

#ifdef CH_MPI
//   MPI_Comm otherComm;
//   MPI_Comm_rank(otherComm, &err); //using unitialized communicator
#endif

  if ( verbose )
    {
      if (err == 0)
        {
          pout() << indent << "testCH_Attach test passed." << std::endl;
        }
      else
        {
          pout() << indent << "testchCH_Attach FAILED!!! "<<err<<std::endl;
        }
    }

#ifdef CH_MPI
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
      else
      {
        break ;
      }
    }
  }
  return ;
}
