#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// This code reads a plotfile and writes another.
// The output file does not have ghost cells

#include <iostream>
#include <strstream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "LayoutIterator.H"

// One more function for MPI
void dumpmemoryatexit();

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  // setChomboMPIErrorHandler();
  MPI_Barrier(Chombo_MPI::comm);  // Barrier #1
#endif

  // ------------------------------------------------
  // parse command line and input file
  // ------------------------------------------------
  // infile must be first
  if (argc < 2)
    {
      cerr << "  need inputs file" << endl;
      abort();
    }

  char* in_file = argv[1];

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  ostrstream outFile;

  outFile << "ng." << in_file << ends;

  // declare memory
  Vector<LevelData<FArrayBox>* > inData;
  Vector<string> inVars; // inData variable names
  Vector<DisjointBoxLayout> inGrids;
  Box inDomain;
  Real inDx, inDt, inTime;
  Vector<int> inRefRatio;
  int inNumLevels;
  string inFileName(in_file);

  ReadAMRHierarchyHDF5(inFileName,
                       inGrids,
                       inData,
                       inVars,
                       inDomain,
                       inDx,
                       inDt,
                       inTime,
                       inRefRatio,
                       inNumLevels);

  Vector<LevelData<FArrayBox>* > outData(inNumLevels);

  // allocate output data -- same domain as input
  for (int level = 0; level < inNumLevels; level++)
    {
      outData[level] = new LevelData<FArrayBox>(inGrids[level],
                                                inVars.size(), IntVect::Zero);
      //copy data to non-ghosted data-structure
      Interval inInterval(0,inVars.size());
      inData[level]->copyTo(*outData[level]);
    }

  WriteAMRHierarchyHDF5(outFile.str(),
                        inGrids,
                        outData,
                        inVars,
                        inDomain,
                        inDx,
                        inDt,
                        inTime,
                        inRefRatio,
                        inNumLevels);

  // clean up memory
  for (int level = 0; level < inNumLevels; level++)
    {
      delete inData[level];
      delete outData[level];
    }

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main
