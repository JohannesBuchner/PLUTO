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
#include "BoxLayout.H"
#include "LayoutIterator.H"
#include "PiecewiseLinearFillPatch.H"

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
  if ((argc < 3) || (argc > 4) )
    {
      cerr << "  usage: addGhosts<*>.ex in_file out_file <nGhost>" << endl;
      abort();
    }

  const char* in_file = argv[1];
  const char* out_file = argv[2];

  // default ghost cells is 1
  int nGhost = 1;
  if (argc == 4)
    {
      nGhost = atoi(argv[3]);
    }

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

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

  IntVect ghostVect(nGhost*IntVect::Unit);

  // use ProblemDomain even though we don't know periodicity info
  ProblemDomain levelDomain(inDomain);

  // allocate output data -- same domain as input
  for (int level = 0; level < inNumLevels; level++)
    {
      outData[level] = new LevelData<FArrayBox>(inGrids[level],
                                                inVars.size(), ghostVect);
      //copy data2 to ghosted data-structure
      Interval inInterval(0,inVars.size());
      inData[level]->copyTo(*outData[level]);

      // now fill data from coarser levels if appropriate
      if (level > 0)
        {

          PiecewiseLinearFillPatch patcher(inGrids[level],
                                           inGrids[level-1],
                                           inVars.size(),
                                           levelDomain,
                                           inRefRatio[level-1],
                                           nGhost);
          Real timeInterpCoeff = 0.0;
          patcher.fillInterp(*outData[level],
                             *outData[level-1],
                             *outData[level-1],
                             timeInterpCoeff,
                             0, 0, outData[level]->nComp());

        }
      // finally, refine problem domain
      if (level < inNumLevels-1)
        {
          levelDomain.refine(inRefRatio[level]);
        }


    }

  string outFileName(out_file);
  WriteAMRHierarchyHDF5(outFileName,
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
