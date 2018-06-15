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
#include <cmath>
#include "Box.H"
#include "SPMD.H"
#include "REAL.H"
#include "MayDay.H"
#include "IntVectSet.H"
#include "BRMeshRefine.H"
#include "parstream.H"
#include <iostream>
#include <fstream>
#include <string>
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

using std::endl;

/// Global variables for test output:
static const char *pgmname = "domainSplitTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

const int destProc = 0;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;
/**
   domainSplitTest returns:
   0: all tests passed
   1: boxes coming out of domain split are not in
   the original domain.
   2: boxes coming out of domain split have a side
   that is too long
   3: domainsplit returned an error code
   4: boxes do not cover  domain
  */
int domainSplitTest(void);

int main(int argc, char* argv[])
{
  using std::ostream;
  bool passed = true;
  int icode = 0;

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  icode = domainSplitTest();

  if (procID() == destProc)
  {
      if (icode != 0)
      {
          pout() << indent << pgmname
               << ": failed with error code " << icode
               << endl;
          passed = false;
      }
      if (passed)
      {
          pout() << indent << pgmname
               << ": passed all tests"
               << endl;
      }
  }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return icode;
}
int domainSplitTest()
{
  int retflag = 0;
  int domlen = 100;
  int maxlen = 25;
  if (SpaceDim > 3)
    {
      domlen = 50;
      maxlen = 13;
    }

  int blockf = 2;

  Box domain(IntVect::Zero, (domlen-1)*IntVect::Unit);

  Vector<Box> vbox;

  domainSplit(domain, vbox, maxlen, blockf);

  // useful to switch between dense and tree IVS's...
  //IntVectSet::setMaxDense(1000000000);
  IntVectSet ivsDom(domain);
  for (int ibox = 0; ibox < vbox.size(); ibox++)
    {
      const Box& boxloc = vbox[ibox];
      if (!domain.contains(boxloc))
        {
          if (verbose)
            pout() << "error: " << boxloc << "not contained in " << domain  << endl;
          retflag += 1;
        }
      if ( verbose )
        pout() << "subtracting box " << boxloc << " from ivs" << endl;
      ivsDom -= boxloc;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (boxloc.size(idir) > maxlen)
            {
              if (verbose)
                pout() << "error: " << boxloc << "is too big in " << idir
                     << " direction (maxlen == " << maxlen << ")"  << endl;
              retflag += 2;
            }
        }
    }
  if (!ivsDom.isEmpty())
    {
      if (verbose)
        pout() << "error: domain not covered " << endl;
      retflag += 4;
    }
  if (verbose)
    {
      pout() << "domain = " << domain << endl;
      pout() << "max box size = " << maxlen << endl;
      pout() << "boxes = " << endl;
      for (int ibox = 0; ibox < vbox.size(); ibox++)
        {
          const Box& boxloc = vbox[ibox];
          pout() << boxloc << endl;
        }
      pout() << "box sizes = " << endl;
      for (int ibox = 0; ibox < vbox.size(); ibox++)
        {
          const Box& boxloc = vbox[ibox];
          for (int idir = 0; idir < SpaceDim; idir ++)
            pout() << boxloc.size(idir)  << "   ";
          pout() << endl;
        }
    }

  // test to see if granularity is screwing things up
  {
    Box testDomainBox(IntVect::Zero, IntVect(D_DECL6(1, 63, 15,
                                                     15,15,15)));
    int ncellmax = 16;
    int block_factor = 2;
    Vector<Box> boxes;

    domainSplit(testDomainBox, boxes, ncellmax, block_factor);

    // check results
    bool passedDomainSplitTest = true;
    if (boxes.size() == 0 )
      {
        // test fails
        passedDomainSplitTest = false;
      }

    // check to be sure that boxes cover domain and no box is
    // greater than ncellmax
    IntVectSet testIVS(testDomainBox);
    for (int n=0; n<boxes.size(); n++)
      {
        testIVS -= boxes[n];
        if (!testDomainBox.contains(boxes[n]))
          passedDomainSplitTest = false;

        if (boxes[n].longside() > ncellmax)
          passedDomainSplitTest =false;
      }

    if (verbose)
      pout() << indent << pgmname << " domainSplit thin domain test: "
             << (( passedDomainSplitTest ) ? "passed" : "failed")
             << endl ;
    if (!passedDomainSplitTest) retflag += 8;
  }

  return retflag;
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
