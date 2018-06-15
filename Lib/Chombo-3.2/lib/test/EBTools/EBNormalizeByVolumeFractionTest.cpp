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
#include "DebugDump.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBNormalizeByVolumeFraction.H"
#include "GeometryShop.H"
#include "SphereIF.H"
#include "EBAMRIO.H"
#include "EBFABView.H"

#include "UsingNamespace.H"

//----------------------------------------------------------------------------
int
main(int argc, char* argv[])
{
   if (SpaceDim == 1)
   {
      pout() << "Not testing in 1D." << endl;
      exit(0);
   }
#ifdef CH_MPI
   MPI_Init(&argc, &argv);
#endif
   {
      // Set up a domain with a spherical (circular in 2D) cavity in the center.
      // The domain has sides of length 2, and the cavity is located at (1,1,1)
      // with a radius of 0.5.
      Real L = 2.0, R = 0.5;
      RealVect center = IntVect::Unit, origin = IntVect::Zero;
      int ebMaxSize = 1024, ebMaxCoarsen = -1;
      int N = 32; // Number of cells on a side.
      RealVect dx = (L / N) * RealVect::Unit;
      int verbosity = 0;
      IntVect lo = IntVect::Zero, hi = (N - 1) * IntVect::Unit;
      Box box(lo, hi);
      SphereIF sphere(R, center, false);
      GeometryShop workshop(sphere, verbosity, dx);
      EBIndexSpace& indexSpace = *Chombo_EBIS::instance();
      indexSpace.define(box, origin, dx[0], workshop, ebMaxSize, ebMaxCoarsen);
      const ProblemDomain& domain = indexSpace.getBox(0);
      DisjointBoxLayout grids = indexSpace.levelGrids(0);
      EBISLayout indexSpaceLayout;
      indexSpace.fillEBISLayout(indexSpaceLayout, grids, domain, 4); // 4 ghost cells in a stencil.

      // Set up a field that is 1 on the regular cells and kappa on the
      // irregular cells.
      // NOTE: The following line gets an "invalid operator" error.
      // LevelData<EBCellFAB> phi(grids, 1);
      EBCellFactory fact(indexSpaceLayout);
      LevelData<EBCellFAB> phi(grids, 1, 4*IntVect::Unit, fact);
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
      {
         const EBISBox& box = indexSpaceLayout[dit()];
         EBCellFAB& phiFAB = phi[dit()];
         phiFAB.setVal(1.0); // Set all cell values to 1.

         // Now go over the irregular cells and set them to the
         // volume fraction.
         const IntVectSet& irregCells = box.getIrregIVS(grids[dit()]);
         for (VoFIterator vit(irregCells, box.getEBGraph()); vit.ok(); ++vit)
         {
            VolIndex vi = vit();
            phiFAB(vi, 0) = box.volFrac(vi);
         }
      }

      // Now normalize phi by the volume fractions.
      EBLevelGrid levelGrid(grids, indexSpaceLayout, domain);
      EBNormalizeByVolumeFraction normalize(levelGrid);
      normalize(phi);

      // Now verify that phi == 1 on the irregular cells.
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
      {
         const EBISBox& box = indexSpaceLayout[dit()];
         EBCellFAB& phiFAB = phi[dit()];
         const IntVectSet& irregCells = box.getIrregIVS(grids[dit()]);
         for (VoFIterator vit(irregCells, box.getEBGraph()); vit.ok(); ++vit)
         {
            VolIndex vi = vit();
            Real val = phiFAB(vi,0);
            Real tol= 1.0e-15;
#ifdef CH_USE_FLOAT
            tol = 1.0e-6;
#endif
            if (Abs(val - 1.0) > tol)
            {
               pout() << "FAIL: phi != 1 on irregular cell!" << endl;
               exit(-1);
            }
         }
      }
   }
   pout() << "EBNormalize test passed" << endl;
#ifdef CH_MPI
   MPI_Finalize();
#endif
}
//----------------------------------------------------------------------------
