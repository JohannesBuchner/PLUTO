#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include  <iostream>
#include "LevelData.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Vector.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "computeSum.H"
#include "ProblemDomain.H"
#include "BCFunc.H"
#include "AMRNodeOp.H"
#include "AMRMultiGrid.H"
#include "BoxIterator.H"
#include "BiCGStabSolver.H"
#include "NodePoissonUtilities.H"
#include "NodeAMRIO.H"
#include <fstream>


/*****/
void runSolver()
{
  PoissonParameters params;
  getPoissonParameters(params);

  ParmParse pp2;

  Vector<DisjointBoxLayout> vectGrids;
  setGrids(vectGrids,  params);

  int numlevels = params.numLevels;
  for (int ilev = 0; ilev < vectGrids.size(); ilev++)
    {
      pout() << "grid at level " << ilev << " = " << vectGrids[ilev] << endl;
    }

  Vector<LevelData<NodeFArrayBox>* > phi(numlevels, NULL);
  Vector<LevelData<NodeFArrayBox>* > rhs(numlevels, NULL);

  int ncomp  = 1;
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      //everybody gets ghost cells
      phi[ilev] = new LevelData<NodeFArrayBox>(vectGrids[ilev], ncomp,  2*IntVect::Unit);
      rhs[ilev] = new LevelData<NodeFArrayBox>(vectGrids[ilev], ncomp,  2*IntVect::Unit);
      LevelData<NodeFArrayBox>& phifabs= *phi[ilev];
      LevelData<NodeFArrayBox>& rhsfabs= *rhs[ilev];
      for (DataIterator dit = phifabs.dataIterator(); dit.ok(); ++dit)
        {
          phifabs[dit()].setVal(0.);
          rhsfabs[dit()].setVal(1.);
        }
    }


  AMRMultiGrid<LevelData<NodeFArrayBox> > solver;
  BiCGStabSolver<LevelData<NodeFArrayBox> >   bottomSolver;
  bottomSolver.m_verbosity = 0;
  nodeDefineSolver(solver, vectGrids, bottomSolver, params);
  int lbase;
  pp2.get("lbase", lbase);

  solver.solve(phi, rhs, numlevels-1, lbase);
  //  solver.relaxOnly(phi, rhs, numlevels-1, lbase);

  nodeOutputData(phi, vectGrids, params.coarsestDomain, params.refRatio,
                 params.coarsestDx, numlevels,  string("phi.hdf5"), string("phi"));

  nodeOutputData(rhs, vectGrids, params.coarsestDomain, params.refRatio,
                 params.coarsestDx, numlevels,  string("rhs.hdf5"), string("rhs"));

  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      delete phi[ilev];
      delete rhs[ilev];
    }
}
/*****/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //scoping trick
  {
    if (argc < 2)
      {
        cerr<< " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }
    char* in_file = argv[1];
    ParmParse  pp(argc-2,argv+2,NULL,in_file);

    int iterations;
    pp.get("iterations", iterations);
    for (int iiter = 0; iiter < iterations; iiter++)
      {
        runSolver();
      }

  }//end scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return(0);
}

