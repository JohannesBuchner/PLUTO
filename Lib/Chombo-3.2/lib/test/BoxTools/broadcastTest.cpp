#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include "parstream.H"
#include "Misc.H"
#include "SPMD.H"
#include "REAL.H"
#include "Box.H"
#include "MayDay.H"
#include <iostream>
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Global variables for test output:
static const char *pgmname = "broadcastTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

using std::endl;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;
/**
  broadcastTest returns:
  0: all tests passed
  1: integer value wrong
  2: real value wrong
  3: box value wrong
  4: vector<int> size wrong
  5: vector<int> value wrong
  6: vector<real> size wrong
  7: vector<real> value wrong
  8: vector<box> size wrong
  9: vector<box> value wrong
  */
int broadcastTest(void);

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
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  icode = broadcastTest();
  if (icode != 0)
    {
      pout() << indent << pgmname
           << ": failed with error code " << icode
           << endl;
      passed = false;
    }

  if (passed)
    pout() << indent << pgmname
         << ": passed all tests"
         << endl;

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return icode;
}
int intVal(int input)
{
  int ifactor = 42;
  return(ifactor*(input+1));
}
Real realVal(Real input)
{
  Real rfactor = 42.0;
  return(rfactor*(input+1.0));
}
Box boxVal(int input)
{
  unsigned int icareful = abs(input);
  IntVect hi = icareful*IntVect::Unit;
  return(Box(IntVect::Zero, hi));
}
/**
  broadcastTest returns:
  0: all tests passed
  1: integer value wrong
  2: real value wrong
  3: box value wrong
  4: vector<int> size wrong
  5: vector<int> value wrong
  6: vector<real> size wrong
  7: vector<real> value wrong
  8: vector<box> size wrong
  9: vector<box> value wrong
  */
int broadcastTest()
{
  int retflag = 0;
  int srcProc = 0;

  if (retflag == 0)
    {
      if (verbose)
        pout() << "testing integer broadcast  " << endl;
      //make a local int on ever processor
      // with a different number
      int ilocl = intVal(procID());
      int icorr = intVal(srcProc);

      broadcast(ilocl, srcProc);

      if (ilocl != icorr)
        {
          if (verbose)
            pout() << "broadcastTest: int  val = " << ilocl
                 << " is wrong; right answer = " << icorr << endl;
          retflag = 1;
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout() << "testing real broadcast  " << endl;
      Real rlocl = realVal(procID());
      Real rcorr = realVal(srcProc);

      broadcast(rlocl, srcProc);

      if (rlocl != rcorr)
        {
          if (verbose)
            pout() << "broadcastTest: real val = " << rlocl
                 << " is wrong; right answer = " << rcorr << endl;
          retflag = 2;
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout() << "testing box broadcast  " << endl;
      Box rlocl = boxVal(procID());
      Box rcorr = boxVal(srcProc);

      broadcast(rlocl, srcProc);

      if (rlocl != rcorr)
        {
          if (verbose)
            pout() << "broadcastTest: box val = " << rlocl
                 << " is wrong- right answer =" << rcorr << endl;
          retflag = 3;
        }
    }
  int loclVsize = 10;
  int corrVsize = 7;
  if (retflag == 0)
    {
      if (verbose)
        pout() << "testing vector<int> broadcast  " << endl;

      Vector<int> ilocl;
      if (procID() != srcProc)
        ilocl.resize(loclVsize);
      else
        ilocl.resize(corrVsize);

      for (int ivec = 0; ivec < ilocl.size(); ivec++)
        ilocl[ivec] = intVal(ivec*procID());

      Vector<int> icorr(corrVsize);
      for (int ivec = 0; ivec < corrVsize; ivec++)
        icorr[ivec] = intVal(ivec*srcProc);

      broadcast(ilocl, srcProc);

      if (ilocl.size() != corrVsize)
        {
          if (verbose)
            pout() << "broadcasttest: vector<int> size = " << ilocl.size() <<
              "  is wrong; correct size = " << corrVsize << endl;
          retflag = 4;
        }
      if (retflag == 0)
        {
          for (int ivec = 0; ivec < corrVsize; ivec++)
            {
              if (ilocl[ivec] != icorr[ivec])
                {
                  if (verbose)
                    pout() << "broadcastTest: int  val at"
                         << ivec << " = "
                         << ilocl[ivec]
                         << " is wrong; right answer = "
                         << icorr[ivec] << endl;
                  retflag = 5;
                }
            }
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout() << "testing vector<Real> broadcast  " << endl;

      Vector<Real> rlocl;
      if (procID() != srcProc)
        rlocl.resize(loclVsize);
      else
        rlocl.resize(corrVsize);

      for (int ivec = 0; ivec < rlocl.size(); ivec++)
        rlocl[ivec] = realVal(ivec*procID());

      Vector<Real> rcorr(corrVsize);
      for (int ivec = 0; ivec < corrVsize; ivec++)
        rcorr[ivec] = realVal(ivec*srcProc);

      broadcast(rlocl, srcProc);

      if (rlocl.size() != corrVsize)
        {
          if (verbose)
            pout() << "broadcasttest: vector<Real> size = "
                 << rlocl.size() <<
              "  is wrong; correct size = " << corrVsize << endl;
          retflag = 6;
        }
      if (retflag == 0)
        {
          for (int ivec = 0; ivec < corrVsize; ivec++)
            {
              if (rlocl[ivec] != rcorr[ivec])
                {
                  if (verbose)
                    pout() << "broadcastTest: Real  val at"
                         << ivec << " = "
                         << rlocl[ivec]
                         << " is wrong; right answer = "
                         << rcorr[ivec] << endl;
                  retflag = 7;
                }
            }
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout() << "testing vector<Box> broadcast  " << endl;

      Vector<Box> blocl(loclVsize);
      if (procID() != srcProc)
        blocl.resize(loclVsize);
      else
        blocl.resize(corrVsize);

      for (int ivec = 0; ivec < blocl.size(); ivec++)
        blocl[ivec] = boxVal(ivec*procID());

      Vector<Box> bcorr(corrVsize);
      for (int ivec = 0; ivec < corrVsize; ivec++)
        bcorr[ivec] = boxVal(ivec*srcProc);

      broadcast(blocl, srcProc);

      if (blocl.size() != corrVsize)
        {
          if (verbose)
            pout() << "broadcasttest: vector<Box> size = " << blocl.size() <<
              "  is wrong; correct size = " << corrVsize << endl;
          retflag = 8;
        }
      if (retflag == 0)
        {
          for (int ivec = 0; ivec < corrVsize; ivec++)
            {
              if (blocl[ivec] != bcorr[ivec])
                {
                  if (verbose)
                    pout() << "broadcastTest: Box  val at"
                         << ivec << " = "
                         << blocl[ivec]
                         << " is wrong; right answer = "
                         << bcorr[ivec] << endl;
                  retflag = 9;
                }
            }
        }
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
