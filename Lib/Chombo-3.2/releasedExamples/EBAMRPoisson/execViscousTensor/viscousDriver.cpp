#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "PoissonUtilities.H"

#include "EBFABView.H"
#include "EBDebugDump.H"

#include "EBViscousTensorOp.H"
#include "EBViscousTensorOpFactory.H"
#include "EBLevelDataOps.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLevelDataOps.H"
#include "EBSimpleSolver.H"
#include "BiCGStabSolver.H"
#include "EBEllipticLoadBalance.H"
#include "EBLoadBalance.H"
#include "LoadBalance.H"
#include "EBLevelGrid.H"
#include "CH_Timer.H"
#include "memusage.H"
#include "memtrack.H"

/********/
void solve(const PoissonParameters&  a_params)
{
  int nvar = SpaceDim;
  Vector<DisjointBoxLayout> grids;
  Vector<EBISLayout>        ebisl;
  getAllIrregRefinedLayouts(grids, ebisl, a_params);


  //define  data
  Vector<LevelData<EBCellFAB>* > phi(a_params.numLevels);
  Vector<LevelData<EBCellFAB>* > rhs(a_params.numLevels);
  ParmParse pp;
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      EBCellFactory factory(ebisl[ilev]);
      phi[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, a_params.ghostPhi, factory);
      rhs[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, a_params.ghostRHS, factory);

      //for now just set phi to zero and rhs to -1.
      EBLevelDataOps::setVal(*phi[ilev], 0.0);
      EBLevelDataOps::setVal(*rhs[ilev], 1.0);
    }

  //create the solver
  AMRMultiGrid<LevelData<EBCellFAB> > solver;
  pout() << "defining  solver" << endl;

  BiCGStabSolver<LevelData<EBCellFAB> > newBottomSolver;
  newBottomSolver.m_verbosity = 0;
  //EBSimpleSolver newBottomSolver;
  //  newBottomSolver.setNumSmooths(32);
  defineViscousTensorSolver(solver, grids, ebisl, newBottomSolver, a_params);


  pout() << "solving " << endl;
  //solve the equation
  solver.init(phi, rhs, a_params.maxLevel, 0);

  solver.solveNoInit(phi, rhs, a_params.maxLevel, 0);


#ifdef CH_USE_HDF5
  bool fileOut;
  pp.get("do_file_output", fileOut);
  if (fileOut)
    {
      pout() << "outputting the answer to file" << endl;
//       UnfreedMemory();
      //output the answer
      char charstr[100];
      sprintf(charstr, "vel.%dd.hdf5", SpaceDim);
      string filename(charstr);
      Vector<string> names(SpaceDim, string("vel"));
      bool replaceCovered;
      Vector<Real> coveredValues(SpaceDim, -1.0);
      Real dumReal =  1.0;
      writeEBHDF5(filename, grids, phi, names,
                  a_params.coarsestDomain.domainBox(), dumReal, dumReal, dumReal,
                  a_params.refRatio, a_params.numLevels,
                  replaceCovered, coveredValues);
      char charstrrhs[100];
      sprintf(charstrrhs, "rhs.%dd.hdf5", SpaceDim);
      string filenamerhs(charstrrhs);
      Vector<string> namesrhs(SpaceDim, string("rhs"));
      writeEBHDF5(filenamerhs, grids, rhs, namesrhs,
                  a_params.coarsestDomain.domainBox(), dumReal, dumReal, dumReal,
                  a_params.refRatio, a_params.numLevels,
                  replaceCovered, coveredValues);
    }
#endif
  //clean up memory
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      delete phi[ilev] ;
      delete rhs[ilev] ;
    }
}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters params;

    //read params from file
    getPoissonParameters(params);

    //define geometry from given params
    definePoissonGeometry(params);



    //solve the stinking problem and output everything
    solve(params);

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  }
  // End scoping trick

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
