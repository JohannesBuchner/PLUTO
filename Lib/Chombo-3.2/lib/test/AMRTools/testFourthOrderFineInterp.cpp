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
#include <iomanip>
#include <sstream>
#include <string>
using std::string;
#include  <iostream>
#include "parstream.H"
#include "CONSTANTS.H"
#include "FABView.H"
#include "DebugDump.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "ParmParse.H"
#include "LoadBalance.H"
#include "FourthOrderFineInterp.H"
// for parallel debugging
#include "CH_Attach.H"
using std::cerr;
using std::endl;
#include "CH_Timer.H"

#include "UsingNamespace.H"

/**
 * Test FourthOrderFineInterp.
 * petermc, 4 Oct 2012
 */

// ---------------------------------------------------------
#ifdef CH_MPI
void reduceReal(Real&           a_val,
                const MPI_Op&   a_mpiOp)
{
  Real recv;
  int resultMPI = MPI_Allreduce(&a_val, &recv,
                                1, MPI_CH_REAL,
                                a_mpiOp, Chombo_MPI::comm);

  if (resultMPI != MPI_SUCCESS)
    { //bark!!!
      MayDay::Error("sorry, but I had a communication error on reduceReal");
    }
  a_val = recv;
}
#endif

// ---------------------------------------------------------
void
setHierarchy(const IntVect& a_lenCoarse,
             const IntVect& a_refineVect,
             Vector<Box>& a_coarboxes,
             Vector<Box>& a_fineboxes,
             Box& a_domc,
             Box& a_domf)
{
  a_domc = Box(IntVect::Zero, a_lenCoarse - IntVect::Unit);

  a_coarboxes.resize(0);

  // If CH_SPACEDIM >= 2: coarboxes[0:2]

  // 1 0
  // 1 2

  // If CH_SPACEDIM == 1: coarBoxes[0:1]

  // 0 1

  int ichop = a_domc.size(0)/2;
  Box bc1 = a_domc;
  Box bc2 = bc1.chop(0, ichop);
#if (CH_SPACEDIM > 1)
  Box bc3 = bc2.chop(1, ichop);
  a_coarboxes.push_back(bc3);
#endif
  a_coarboxes.push_back(bc1);
  a_coarboxes.push_back(bc2);

  IntVect lenFine = a_refineVect * a_lenCoarse;
  a_domf = Box(IntVect::Zero, lenFine - IntVect::Unit);

  IntVect u = lenFine / 8;

#if (CH_SPACEDIM ==1)

  // a_fineboxes[0:2]

  // 0 0 * 1 1 * 2 2

  a_fineboxes.resize(3);
  a_fineboxes[0].define(IntVect(0),
                        IntVect(2*u[0]) - IntVect::Unit);
  a_fineboxes[1].define(IntVect(3*u[0]),
                        IntVect(5*u[0]) - IntVect::Unit);
  a_fineboxes[2].define(IntVect(6*u[0]),
                        IntVect(8*u[0]) - IntVect::Unit);

#elif (CH_SPACEDIM ==2)

  // a_fineboxes[0:3]

  // * * * * * * * *
  // * * * * * * * *
  // 0 0 0 0 * * * *
  // 0 0 0 0 1 1 * *
  // 0 0 0 0 1 1 * *
  // 0 0 0 0 * 3 3 3
  // * * 2 2 2 * * *
  // * * 2 2 2 * * *

  a_fineboxes.resize(4);
  a_fineboxes[0].define(IntVect(0, 2*u[1]),
                        IntVect(4*u[0], 6*u[1]) - IntVect::Unit);
  a_fineboxes[1].define(IntVect(4*u[0], 3*u[1]),
                        IntVect(6*u[0], 5*u[1]) - IntVect::Unit);
  a_fineboxes[2].define(IntVect(2*u[0], 0),
                        IntVect(5*u[0], 2*u[1]) - IntVect::Unit);
  a_fineboxes[3].define(IntVect(5*u[0], 2*u[1]),
                        IntVect(8*u[0], 3*u[1]) - IntVect::Unit);

#elif (CH_SPACEDIM==3)

  // a_fineboxes[0:3] extended to cover everything in the third dimension

  // * * * * * * * *
  // * * * * * * * *
  // 0 0 0 0 * * * *
  // 0 0 0 0 1 1 * *
  // 0 0 0 0 1 1 * *
  // 0 0 0 0 * 3 3 3
  // * * 2 2 2 * * *
  // * * 2 2 2 * * *

  a_fineboxes.resize(4);
  a_fineboxes[0].define(IntVect(0, 2*u[1], 0),
                        IntVect(4*u[0], 6*u[1], 8*u[2]) - IntVect::Unit);
  a_fineboxes[1].define(IntVect(4*u[0], 3*u[1], 0),
                        IntVect(6*u[0], 5*u[1], 8*u[2]) - IntVect::Unit);
  a_fineboxes[2].define(IntVect(2*u[0], 0, 0),
                        IntVect(5*u[0], 2*u[1], 8*u[2]) - IntVect::Unit);
  a_fineboxes[3].define(IntVect(5*u[0], 2*u[1], 0),
                        IntVect(8*u[0], 3*u[1], 8*u[2]) - IntVect::Unit);

#else
  MayDay::Error("testFourthOrderFineInterp does not work with DIM > 3");
  // error bogus spacedim;
#endif

  return ;
}


// ---------------------------------------------------------
Real avgFunVal(const RealVect& a_cellLo,
               const RealVect& a_cellHi,
               const RealVect& a_XYZ)
{
  // Return average of function
  // f(x, y, z) = sin(2*pi*x/X) * cos(4*pi*y/Y) * sin(6*pi*z/Z)
  // over the rectangle a_cellLo : a_cellHi.

  // Average of sin(2*pi*x/X) over [a, b] is
  // (-X/(2*pi)) / (b-a) * (cos(2*pi*b/X) - cos(2*pi*a/X)).

  // Average of cos(4*pi*y/Y) over [a, b] is
  // (X/(4*pi)) / (b-a) * (sin(4*pi*b/Y) - sin(4*pi*a/Y)).

  // Average of sin(6*pi*z/Z) over [a, b] is
  // (-Z/(6*pi)) / (b-a) * (cos(6*pi*b/Z) - cos(6*pi*a/Z)).

  // Real Pi = 4.0 * atan(1.0);
  Real fac, a, b;
  Real retval = 1.;

  fac = 2. * Pi / a_XYZ[0];
  a = a_cellLo[0];
  b = a_cellHi[0];
  retval *= ( (-1./fac) / (b-a) ) * (cos(fac*b) - cos(fac*a));
    
#if CH_SPACEDIM >= 2
  fac = 4. * Pi / a_XYZ[1];
  a = a_cellLo[1];
  b = a_cellHi[1];
  retval *= ( (1./fac) / (b-a) ) * (sin(fac*b) - sin(fac*a));
#endif

#if CH_SPACEDIM >= 3
  fac = 6. * Pi / a_XYZ[2];
  a = a_cellLo[2];
  b = a_cellHi[2];
  retval *= ( (-1./fac) / (b-a) ) * (cos(fac*b) - cos(fac*a));
#endif

  return retval;
}


// ---------------------------------------------------------
int main(int argc, char* argv[])
{
  int status = 0; // number of errors detected.
  // Do nothing if DIM > 3.
#if (CH_SPACEDIM <= 3)
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  //scoping trick
  {
    // test parameters
    const int testOrder = 4; // expected order of convergence

    // A test is considered as a failure
    //  if its convergence rate is smaller than below.
    Real targetConvergeRate = testOrder*0.9;

    // number of ghost cells being interpolated
    int numGhost = 0;
    // This does NOT fill in fine ghost cells; if you want it to do so,
    // you need to define an IntVectSet as is done in FourthOrderFillPatch.

    // refinement ratio between coarse and fine levels
    int refRatio = 4;

    Vector<Interval> fixedDimsAll;
    // Always include empty interval, meaning NO fixed dimensions.
    fixedDimsAll.push_back(Interval());
    if (SpaceDim == 3)
      {
        fixedDimsAll.push_back(Interval(SpaceDim-1, SpaceDim-1));
      }

    // real domain has length 1. in every dimension
    RealVect physLength = RealVect::Unit;

    for (int ifixed = 0; ifixed < fixedDimsAll.size(); ifixed++)
      {
        const Interval& fixedDims = fixedDimsAll[ifixed];
        int nfixed = fixedDims.size();

        IntVect interpUnit = IntVect::Unit;
        IntVect refineVect = refRatio * IntVect::Unit;
        for (int dirf = fixedDims.begin(); dirf <= fixedDims.end(); dirf++)
          {
            interpUnit[dirf] = 0;
            refineVect[dirf] = 1;
          }

        // No ghost cells in fixed dimensions.
        IntVect ghostVect = numGhost * interpUnit;

#ifdef CH_USE_FLOAT
        // Single precision doesn't get us very far.
#if (CH_SPACEDIM == 1)
        int domainLengthMin = 16;
        int domainLengthMax = 32;
#else
        int domainLengthMin = 32;
        int domainLengthMax = 64;
#endif
#endif

#ifdef CH_USE_DOUBLE
        int domainLengthMin = 32;
#if (CH_SPACEDIM >= 3)
        int domainLengthMax = 64;
        if (nfixed > 0) domainLengthMax = 128;
#else
        int domainLengthMax = 512;
#endif
#endif

        Vector<int> domainLengths;
        int len = domainLengthMin;
        while (len <= domainLengthMax)
          {
            domainLengths.push_back(len);
            len *= 2;
          }
        int nGrids = domainLengths.size();

        if (nGrids > 0)
          {
            pout() << endl
                   << "Testing FourthOrderFineInterp DIM=" << SpaceDim;
            if (nfixed > 0)
              {
                pout() << " fixing " << fixedDims.begin()
                       << ":" << fixedDims.end();
              }
            pout() << endl;
            pout() << "size  "
                   << "max diff   rate";
            pout() << endl;
          }

        Real diffMaxCoarser;
        for (int iGrid = 0; iGrid < nGrids; iGrid++)
          {
            int domainLength = domainLengths[iGrid];
            IntVect lenCoarse = domainLength * IntVect::Unit;

            RealVect dxCoarseVect = physLength / RealVect(lenCoarse);
            RealVect dxFineVect = dxCoarseVect / RealVect(refineVect);

            Box coarDomain, fineDomain;
            Vector<Box> coarBoxes, fineBoxes;
            setHierarchy(lenCoarse,
                         refineVect,
                         coarBoxes,
                         fineBoxes,
                         coarDomain,
                         fineDomain);

            mortonOrdering(coarBoxes);
            mortonOrdering(fineBoxes);

            Vector<int> procCoar(coarBoxes.size());
            Vector<int> procFine(fineBoxes.size());
            LoadBalance(procCoar, coarBoxes);
            LoadBalance(procFine, fineBoxes);

            DisjointBoxLayout dblCoar(coarBoxes, procCoar);
            DisjointBoxLayout dblFine(fineBoxes, procFine);

            ProblemDomain probCoar(coarDomain);
            ProblemDomain probFine(fineDomain);

            FourthOrderFineInterp interpolator;
            int ncomp = 1;
            // Set coarseGhostsFill = ceil(numGhosts / refRatio).
            int coarseGhostsFill = numGhost / refRatio;
            if (coarseGhostsFill * refRatio < numGhost) coarseGhostsFill++;

            interpolator.define(dblFine, ncomp, refRatio,
                                probFine, coarseGhostsFill, fixedDims);

            LevelData<FArrayBox> exactCoarse(dblCoar, ncomp);

            // Set exact value of function on coarse level.
            for (DataIterator dit = exactCoarse.dataIterator(); dit.ok(); ++dit)
              {
                CH_TIME("computing exactCoarse on FAB");
                FArrayBox& exactCoarseFab = exactCoarse[dit];
                const Box& bx = exactCoarseFab.box();
                // Get average of function on each cell.
                for (BoxIterator bit(bx); bit.ok(); ++bit)
                  {
                    const IntVect& iv = bit();
                    RealVect cellLo = RealVect(iv) * dxCoarseVect;
                    RealVect cellHi = cellLo + dxCoarseVect;
                    exactCoarseFab(iv, 0) = avgFunVal(cellLo, cellHi, physLength);
                  }
              }
            
            /*
              Define data holders on fine level.
            */
            Interval intvlExact(0, ncomp-1);
            Interval intvlCalc(ncomp, 2*ncomp-1);
            Interval intvlDiff(2*ncomp, 3*ncomp-1);
            LevelData<FArrayBox> allFine;
            LevelData<FArrayBox> exactFine, calcFine, diffFine;
            { CH_TIME("allocate allFine");
              allFine.define(dblFine, 3*ncomp, ghostVect);
              
              aliasLevelData(exactFine, &allFine, intvlExact);
              aliasLevelData(calcFine, &allFine, intvlCalc);
              aliasLevelData(diffFine, &allFine, intvlDiff);
            }
            
            /*
              Fill in data holders on fine level.
            */

            // Set exact value of function on fine level.
            Real exactMax = 0.;
            for (DataIterator dit = exactFine.dataIterator(); dit.ok(); ++dit)
              {
                CH_TIME("computing exactFine on FAB");
                FArrayBox& exactFineFab = exactFine[dit];
                const Box& bx = exactFineFab.box();
                // Get average of function on each cell.
                for (BoxIterator bit(bx); bit.ok(); ++bit)
                  {
                    const IntVect& iv = bit();
                    RealVect cellLo = RealVect(iv) * dxFineVect;
                    RealVect cellHi = cellLo + dxFineVect;
                    exactFineFab(iv, 0) = avgFunVal(cellLo, cellHi, physLength);
                  }
                Real exactFabMax = exactFineFab.norm(0);
                if (exactFabMax > exactMax) exactMax = exactFabMax;
              }
#ifdef CH_MPI
            reduceReal(exactMax, MPI_MAX);
#endif

            // Fill in calcFine with function values at valid cells.
            interpolator.interpToFine(calcFine, exactCoarse);

            /*
              Got results.  Now take differences,
              diffFine = calcFine - exactFine.
            */
            Real diffMax = 0.;
            for (DataIterator dit = diffFine.dataIterator(); dit.ok(); ++dit)
              { CH_TIME("calculating diffFine");
                FArrayBox& diffFineFab = diffFine[dit];
                const FArrayBox& exactFineFab = exactFine[dit];
                const FArrayBox& calcFineFab = calcFine[dit];
                
                diffFineFab.setVal(0.);
                diffFineFab.copy(exactFineFab);
                diffFineFab.minus(calcFineFab);

                Real diffFabMax = diffFineFab.norm(0);
                if (diffFabMax > diffMax) diffMax = diffFabMax;
                // dummy statement in order to get around gdb bug
                int dummy_unused = 0; dummy_unused = 0;
              }
#ifdef CH_MPI
            reduceReal(diffMax, MPI_MAX);
#endif

            pout() << setw(4) << domainLength;
            pout() << "  " << scientific << setprecision(4)
                   << diffMax;
            if (iGrid > 0)
              {
                Real ratio = diffMaxCoarser / diffMax;
                Real rate = log(ratio) / log(2.0);
                pout () << " " << fixed << setprecision(2) << rate;
                if (rate < targetConvergeRate)
                  {
                    status += 1;
                  }
              }
            pout() << endl;
            
            diffMaxCoarser = diffMax;
          } // end loop over grid sizes
      } // end loop over which dimensions are fixed
    if (status==0)
      {
        pout() <<  "All tests passed!\n";
      }
    else
      {
        pout() <<  status << " tests failed!\n";
      }
  } // end scoping trick
#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
#endif

  return status;
}
