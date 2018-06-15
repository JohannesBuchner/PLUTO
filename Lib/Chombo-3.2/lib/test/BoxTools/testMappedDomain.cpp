#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/// Global variables for test output:
static const char *pgmname = "testMappedDomain" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

#if 0 // No MappedDomain.H anywhere in sight!

#include "MappedDomain.H"
#include "SPMD.H"
#include "parstream.H"
#include "UsingNamespace.H"


const int destProc = 0;

int testMappedDomain();

int main(int argc, char* argv[])
{
  using std::ostream;
  bool passed = true;
  int icode = 0;

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif


  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << std::endl ;

  icode = testMappedDomain();

  if (procID() == destProc)
  {
      if (icode != 0)
      {
        pout() << indent << pgmname
               << ": failed with error code " << icode
               << std::endl;
          passed = false;
      }
      if (passed)
      {
        pout() << indent << pgmname
               << ": passed all tests"
               << std::endl;
      }
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return icode;
}


int testMappedDomain()
{
  int boxsize = 8;
  int boxcount   = 6;
  int blocksize= boxsize*boxcount;

  IntVect low = 32*IntVect::Unit;
  IntVect hi  = (32+(blocksize)-1)*IntVect::Unit;

  Box block1(low, hi);
  Box block2=block1;
  block2.shift(0, blocksize);
  Box block3=block1;
  block3.shift(0,-blocksize);

  return 0;
}

#else  // 0

#include <iostream>
int main()
{
    pout() << "testMappedDomain dummy implementation\n";
}
#endif
