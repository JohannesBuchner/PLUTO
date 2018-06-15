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
static const char *pgmname = "gatherTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

const int destProc = 0;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;
/**
  gatherTest returns:
  0: all tests passed
  1: size of vector<int> wrong
  2: values in vector<int> wrong
  3: size of vector<real> wrong
  4: values in vector<real> wrong
  5: size of vector<vector<int>> wrong
  6: values in vector<vector<int>> wrong
  7: size of vector<vector<real>> wrong
  8: values in vector<vector<real>> wrong
  9: size of vector<box> wrong
  10: values in vector<box> wrong
  11: size of vector<vector<box>> wrong
  12: values in vector<vector<box>> wrong
  */
int gatherTest(void);

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

  icode = gatherTest();

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
int gatherTest()
{
  int retflag = 0;
  int nProcess = numProc();
  if (retflag == 0)
    {
      if (verbose)
        pout ()  << indent2 << "testing integer gather " << endl;
      //make a local int on ever processor
      // with a different number
      int loclInt = intVal(procID());
      if (verbose)
        pout ()  << indent2 << "loclInt = " << loclInt << endl;
      Vector<int> allInts;

      gather(allInts, loclInt, destProc);

      if (procID() == destProc)
        {
          if (allInts.size() != nProcess)
            {
              if (verbose)
                pout ()  << indent2 << "gatherTest: integer vector size wrong" << endl;
              retflag = 1;
            }
          for (int iproc = 0; iproc < nProcess; iproc++)
            {
              if (retflag != 0) break;
              int correctInt = intVal(iproc);
              if (verbose)
                {
                  pout ()  << indent2 << "correctInt = " << correctInt << endl;
                  pout ()  << indent2 << "allInts[" << iproc << "] = "
                       <<  allInts[iproc] << endl;
                }
              if (allInts[iproc] != correctInt)
                {
                  if (verbose)
                    {
                      pout ()  << indent2 << "gatherTest: integer vector wrong at "
                           << iproc << " = " << allInts[iproc] << endl;
                      pout ()  << indent2 << "correct value is" << correctInt << endl;
                    }
                  retflag = 2;
                }
            }
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout ()  << indent2 << "testing real gather " << endl;
      //make a local Real on ever processor
      // with a different number
      Real loclReal = realVal(procID());
      Vector<Real> allReals;

      gather(allReals, loclReal, destProc);

      if (procID() == destProc)
        {
          if (allReals.size() != nProcess)
            {
              if (verbose)
                pout ()  << indent2 << "gatherTest: Real vector size wrong" << endl;
              retflag = 3;
            }
          for (int iproc = 0; iproc < nProcess; iproc++)
            {
              if (retflag != 0) break;
              Real correctReal = realVal(iproc);
              if (allReals[iproc] != correctReal)
                {
                  if (verbose)
                    {
                      pout ()  << indent2 << "gatherTest: Real vector wrong at "
                           << iproc << " = " << allReals[iproc] << endl;
                      pout ()  << indent2 << "correct value is" << correctReal << endl;
                    }
                  retflag = 4;
                }
            }
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout ()  << indent2 << "testing vector<int> gather " << endl;

      int nelem = 10;
      Vector<int> loclVInt(nelem);
      for (int ivec = 0; ivec < nelem; ivec++)
        loclVInt[ivec] = intVal(ivec*procID());
      Vector< Vector<int> > allVInts;

      gather(allVInts, loclVInt, destProc);

      if (procID() == destProc)
        {
          if (allVInts.size() != nProcess)
            {
              if (verbose)
                pout ()  << indent2 << "gatherTest: vector< vector<int> > size wrong"
                     << endl;
              retflag = 5;
            }
          for (int iproc = 0; iproc < nProcess; iproc++)
            {
              if (retflag != 0) break;

              for (int ivec = 0; ivec < nelem; ivec++)
                {
                  int correctInt = intVal(ivec*iproc);
                  if (allVInts[iproc][ivec] != correctInt)
                    {
                      if (verbose)
                        {
                          pout ()  << indent2 << "gatherTest: integer vector wrong at [][]"
                               << iproc << ivec <<
                            " = " << allVInts[iproc][ivec] << endl;
                          pout ()  << indent2 << "correct value is" << correctInt << endl;
                        }
                      retflag = 6;
                    }
                }
            }
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout ()  << indent2 << "testing vector<Real> gather " << endl;
      int nelem = 10;
      Vector<Real> loclVReal(nelem);
      for (int ivec = 0; ivec < nelem; ivec++)
        loclVReal[ivec] = realVal(ivec*procID());
      Vector< Vector<Real> > allVReals;

      gather(allVReals, loclVReal, destProc);

      if (procID() == destProc)
        {
          if (allVReals.size() != nProcess)
            {
              if (verbose)
                pout ()  << indent2 << "gatherTest: vector< vector<Real> > size wrong"
                     << endl;
              retflag = 7;
            }
          for (int iproc = 0; iproc < nProcess; iproc++)
            {
              if (retflag != 0) break;

              for (int ivec = 0; ivec < nelem; ivec++)
                {
                  Real correctReal = realVal(ivec*iproc);
                  if (allVReals[iproc][ivec] != correctReal)
                    {
                      if (verbose)
                        {
                          pout ()  << indent2 << "gatherTest: Real vector wrong at [][]"
                               << iproc << ivec <<
                            " = " << allVReals[iproc][ivec] << endl;
                          pout ()  << indent2 << "correct value is" << correctReal << endl;
                        }
                      retflag = 8;
                    }
                }
            }
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout ()  << indent2 << "testing Box gather " << endl;
      //make a local int on ever processor
      // with a different number
      Box loclBox = boxVal(procID());
      if (verbose)
        pout ()  << indent2 << "loclBox = " << loclBox << endl;
      Vector<Box> allBox;

      gather(allBox, loclBox, destProc);

      if (procID() == destProc)
        {
          if (allBox.size() != nProcess)
            {
              if (verbose)
                pout ()  << indent2 << "gatherTest: box vector size wrong" << endl;
              retflag = 9;
            }
          for (int iproc = 0; iproc < nProcess; iproc++)
            {
              if (retflag != 0) break;
              Box correctBox = boxVal(iproc);
              if (verbose)
                {
                  pout ()  << indent2 << "correctBox = " << correctBox << endl;
                  pout ()  << indent2 << "allBox[" << iproc << "] = "
                       <<  allBox[iproc] << endl;
                }
              if (allBox[iproc] != correctBox)
                {
                  if (verbose)
                    {
                      pout ()  << indent2 << "gatherTest: Box vector wrong at "
                           << iproc << " = " << allBox[iproc] << endl;
                      pout ()  << indent2 << "correct value is" << correctBox << endl;
                    }
                  retflag = 10;
                }
            }
        }
    }
  if (retflag == 0)
    {
      if (verbose)
        pout ()  << indent2 << "testing vector<Box> gather " << endl;
      int nelem = 10;
      Vector<Box> loclVBox(nelem);
      for (int ivec = 0; ivec < nelem; ivec++)
        loclVBox[ivec] = boxVal(ivec*procID());
      Vector< Vector<Box> > allVBoxs;

      gather(allVBoxs, loclVBox, destProc);

      if (procID() == destProc)
        {
          if (allVBoxs.size() != nProcess)
            {
              if (verbose)
                pout ()  << indent2 << "gatherTest: vector< vector<Box> > size wrong"
                     << endl;
              retflag = 11;
            }
          for (int iproc = 0; iproc < nProcess; iproc++)
            {
              if (retflag != 0) break;

              for (int ivec = 0; ivec < nelem; ivec++)
                {
                  Box correctBox = boxVal(ivec*iproc);
                  if (allVBoxs[iproc][ivec] != correctBox)
                    {
                      if (verbose)
                        {
                          pout ()  << indent2 << "gatherTest: Box vector wrong at [][]"
                               << iproc << ivec <<
                            " = " << allVBoxs[iproc][ivec] << endl;
                          pout ()  << indent2 << "correct value is" << correctBox << endl;
                        }
                      retflag = 12;
                    }
                }
            }
        }
    }
  if (retflag == 0)
    {
      if (verbose)
          pout () << indent2 << "testing IntVectSet gather " << endl;
      int nelem = 10;
      Vector<Box> loclVBox(nelem);
      for (int ivec = 0; ivec < nelem; ivec++)
      {
        const IntVect iv = (ivec+nelem*procID()) * IntVect::Unit;
        loclVBox[ivec] = Box(iv,iv);
      }
      IntVectSet ivs;
      for (int i = 0; i < nelem; ++i)
      {
        ivs |= loclVBox[i];
      }
      Vector<IntVectSet> all_ivses;

      gather(all_ivses, ivs, destProc);

      if (procID() == destProc)
        {
          bool pass = true;
          if (all_ivses.size() != nProcess)
            {
                pass = false;
              if (verbose)
                pout ()  << indent2 << "gatherTest: vector< IntVectSet > size wrong"
                     << endl;
              retflag = 13;
            }
          for (int iproc = 0; iproc < nProcess; iproc++)
            {
              if (retflag != 0) break;

              for (int ivec = 0; ivec < nelem; ivec++)
                {

                    Vector<Box> correct_vbox(nelem);
                    for (int ivec = 0; ivec < nelem; ivec++)
                    {
                        const IntVect iv = (ivec+nelem*iproc) * IntVect::Unit;
                        correct_vbox[ivec] = Box(iv,iv);
                    }
                    IntVectSet correct_ivs;
                    for (int i = 0; i < nelem; ++i)
                    {
                        correct_ivs |= correct_vbox[i];
                    }
                    IntVectSet diff1(all_ivses[iproc]);
                    diff1 -= correct_ivs;

                    IntVectSet diff2(correct_ivs);
                    diff2 -= all_ivses[iproc];

                    if ( !diff1.isEmpty() || !diff2.isEmpty() )
                    {
                        pass = false;
                      if (verbose)
                        {
                          pout ()  << indent2 << "gatherTest: IntVectSet wrong at iproc = "
                                << iproc << endl;
                          pout () << indent2 << "received IntVectSet: " << endl;
                          all_ivses[iproc].printBoxes(pout ());
                          pout () << endl;
                          pout ()  << indent2 << "correct value is: " << endl;
                          correct_ivs.printBoxes(pout ());
                          pout () << endl;
                        }
                      retflag = 14;
                    }
                }
            }
            if (verbose)
              {
                  if (pass)
                  {
                      pout () << indent2 << "IntVectSet test passed." << endl;
                  }
                  else
                  {
                      pout () << indent2 << "IntVectSet test FAILED!!!" << endl;
                  }
              }
        }
    }
#ifdef CH_MPI
  //make sure all processors have the same return flag
  int rootProc = destProc;
  MPI_Bcast(&retflag, 1, MPI_INT, rootProc, Chombo_MPI::comm);
#endif
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
