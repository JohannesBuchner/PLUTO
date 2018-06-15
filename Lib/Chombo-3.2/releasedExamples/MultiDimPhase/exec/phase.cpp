#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//===========================================================================
// driver.cpp
//
//===========================================================================
#include <iostream>

// 2d includes
#define CH_SPACEDIM 2
#include "ParmParse.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "amrPhase.H"
#include "CH_Attach.H"

#include "UsingNamespace.H"

//using namespace Chombo::D2;


const unsigned int nComps = 2;

//===========================================================================
// simple multi-dimensional phase-space integrator
//
//===========================================================================
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  //registerDebugger();

  { // Begin nested scope

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank=0;
    number_procs=1;
#endif

    if (argc < 2)
      { std::cerr << " usage: " << argv[0] << " <input_file>\n"; exit(0); }
    char* in_file = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,in_file);

    amrPHASE amrObject;

    // set up initial grids, initialize data, etc.
    amrObject.initialize();

    ParmParse pp2("main");
    int maxStep;
    Real maxTime;
    //Real startTime;
    pp2.get("maxTime", maxTime);
    pp2.get("maxStep", maxStep);

    amrObject.run(maxTime, maxStep);

  }  // end nested scope


#ifdef CH_MPI
  MPI_Finalize();
#endif

  return 0;
}
