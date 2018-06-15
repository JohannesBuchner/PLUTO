#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// this code is the chombo version of the compare utility
// which takes two plotfiles (one fine "exact" solution, and
// one coarser "computed" solution) and computes L1, L2, and
// Max "errors".  Norms are computed only on valid regions
// of each AMR grid.  assumes that exact solution is a single
// fine grid.

#include <iostream>
#include <strstream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "LayoutIterator.H"
#include "LoHiSide.H"
#include "GhostBC.H"
#include "FineInterp.H"

#include "computeSum.H"
#include "computeNorm.H"

#include "HOExtrapBC.H"
#include "AverageHO_F.H"

Vector<string> noAverageVars;
Vector<string> sumVars;

void init(string&         a_exactRoot,
          string&         a_computedRoot,
          string&         a_errorRoot,
          int&            a_intFieldSize,
          int&            a_numCrseStart,
          Vector<string>& a_errorVars,
          int&            a_numCrseFinish,
          int&            a_crseStep,
          int&            a_crseMult,
          bool&           a_doPlots,
          bool&           a_isTimeDep,
          bool&           a_isSameSize,
          Real&           a_bogusValue,
          bool&           a_computeRelativeError,
          bool&           a_removeMean,
          bool&           a_doGhostCells,
          bool&           a_useUnitDomain,
          bool&           a_HOaverage,
          bool&           a_verbose);

void computeAMRError(Vector<LevelData<FArrayBox>* >&       a_error,
                     const Vector<string>&                 a_errorVars,
                     const Vector<LevelData<FArrayBox>* >& a_computedSoln,
                     const Vector<string>&                 a_computedVars,
                     const Vector<DisjointBoxLayout>&      a_computedGrids,
                     const Real                            a_computedDx,
                     const Vector<int>&                    a_computedRefRatio,
                     const Vector<LevelData<FArrayBox>* >& a_exactSoln,
                     const Vector<string>&                 a_exactVars,
                     const Real                            a_exactDx,
                     Real                                  a_bogus_value,
                     bool                                  a_HOaverage,
                     bool                                  a_computeRelativeError);

void computeSameSizeError(Vector<LevelData<FArrayBox>* >&       a_error,
                          const Vector<string>&                 a_errorVars,
                          const Vector<LevelData<FArrayBox>* >& a_computedSoln,
                          const Vector<string>&                 a_computedVars,
                          const Vector<DisjointBoxLayout>&      a_computedGrids,
                          const Real                            a_dx,
                          const Vector<int>&                    a_refRatio,
                          const Vector<LevelData<FArrayBox>* >& a_exactSoln,
                          const Vector<string>&                 a_exactVars,
                          Real                                  a_bogus_value,
                          bool                                  a_computeRelativeError,
                          bool                                  a_doGhostCells=false);

void constructErrorNames(Vector<string>&       a_errorNames,
                         const Vector<string>& a_errorVars);

void constructPlotFileName(ostrstream&   a_fileName,
                           const string& a_fileRoot,
                           const int     a_intFieldWidth,
                           const int     a_step);

bool nonAverageVar(const string& a_errorName);

bool sumVar(const string& a_errorName);

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
  ParmParse pp(argc-2, argv+2, NULL, in_file);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  pout().setf(ios::scientific);
  pout().precision(4);

  // set defaults
  bool isSameSize = false;
  bool doGhostCells = false;
  bool doPlots = true;
  bool isTimeDep = false;
  bool verbose = false;
  bool computeRelativeError = false;
  bool removeMean = false;
  bool useUnitDomain = false;
  bool HOaverage = false;

  // this is the initial value for the error
  Real bogusValue = 0.0;

  string exactRoot, computedRoot, errorRoot;

  Vector<string> errorVars;

  int numCrseFinish, numCrseStart, crseStep, crseMult, intFieldSize;

  init(exactRoot, computedRoot, errorRoot, intFieldSize,
       numCrseStart, errorVars, numCrseFinish, crseStep,
       crseMult, doPlots, isTimeDep, isSameSize, bogusValue,
       computeRelativeError, removeMean, doGhostCells, useUnitDomain,
       HOaverage,verbose);


  int nStep, exactStep;

  if (!isTimeDep)
  {
    crseStep = 1;
    numCrseFinish = numCrseStart;
  }

  for (nStep = numCrseStart; nStep <= numCrseFinish; nStep += crseStep)
  {
    if (verbose)
    {
      pout() << "starting step " << nStep << endl;
    }

    exactStep = nStep*crseMult;

    ostrstream exactFile;
    ostrstream computedFile;
    ostrstream errorFile;

    exactFile.fill('0');
    computedFile.fill('0');
    errorFile.fill('0');

    if (isTimeDep)
    {
      constructPlotFileName(exactFile, exactRoot, intFieldSize, exactStep);
      constructPlotFileName(computedFile, computedRoot, intFieldSize, nStep);
      constructPlotFileName(errorFile, errorRoot, intFieldSize, nStep);

    }
    else
    {
      // if not time dependent, file roots are really filenames
      exactFile << exactRoot << ends;
      computedFile << computedRoot << ends;
      errorFile << errorRoot << ends;
    }

    pout() << "exact Filename = " << exactFile.str() << endl;
    pout() << "computed Filename = " << computedFile.str() << endl;
    if (doPlots)
    {
      pout() << "error Filename = " << errorFile.str() << endl;
    }

    // declare memory
    Vector<LevelData<FArrayBox>* > exactSoln;
    Vector<string> exactVars; // exact solution variable names
    Vector<DisjointBoxLayout> exactGrids;
    Box exactDomain;
    Real exactDx, exactDt, exactTime;
    Vector<int> exactRefRatio;
    int exactNumLevels;
    IntVect ghostVect = IntVect::Unit;
    string exactFileName(exactFile.str());

    // get exact solution
    if (verbose)
    {
      pout() << "read exact solution..." << endl;
    }

    ReadAMRHierarchyHDF5(exactFileName,
                         exactGrids,
                         exactSoln,
                         exactVars,
                         exactDomain,
                         exactDx,
                         exactDt,
                         exactTime,
                         exactRefRatio,
                         exactNumLevels);

    if (verbose)
    {
      pout () << "done reading exact soln" << endl;
    }

    // we assume that exact soln is single-grid if we're doing averaging
    if (!isSameSize)
    {
      CH_assert(exactNumLevels == 1);
    }

    Vector<LevelData<FArrayBox>* > computedSoln;
    Vector<string> computedVars; // computed soln variable names
    Vector<DisjointBoxLayout> computedGrids;
    Box computedDomain;
    Real computedDx, computedDt, computedTime;
    Vector<int> computedRefRatio;
    int computedNumLevels;
    string computedFileName(computedFile.str());

    //ghostVect = IntVect::Zero;
    // now read in computed solution
    if (verbose)
    {
      pout() << "read computed solution..." << endl;
    }

    ReadAMRHierarchyHDF5(computedFileName,
                         computedGrids,
                         computedSoln,
                         computedVars,
                         computedDomain,
                         computedDx,
                         computedDt,
                         computedTime,
                         computedRefRatio,
                         computedNumLevels);

    if (verbose)
    {
      pout() << "done reading computed solution" << endl;
    }

    // reality check
    if ((computedDomain != exactDomain) && isSameSize)
    {
      MayDay::Error("Incompatible exact and computed domains for sameSize comparison");
    }

    int numExact    = exactVars.size();
    int numComputed = computedVars.size();

    int numError = errorVars.size();
    // If no errorVars were specified
    if (numError == 0)
    {
      // Set errorVars to the intersection of exactVars and computedVars
      // This numVars^2 method should be changed to something more efficient
      for (int iExact = 0; iExact < numExact; iExact++)
      {
        for (int iComp = 0; iComp < numComputed; iComp++)
        {
          if (exactVars[iExact] == computedVars[iComp])
          {
            errorVars.push_back(exactVars[iExact]);
            break;
          }
        }
      }

      numError = errorVars.size();
    }
    else
    {
      // if errorVars were specified, then do a quick check that
      // they're present in both exactVars and computedVars
      for (int errVarNo = 0; errVarNo<errorVars.size(); ++errVarNo)
        {
          bool foundComputed = false;
          for (int i=0; i<numComputed; i++)
            {
              if (errorVars[errVarNo] == computedVars[i])
                {
                  foundComputed = true;
                }
            } // end loop over exact variables
          if (!foundComputed)
            {
              pout() << "errorVar " << errorVars[errVarNo]
                     << " not found in computed solution!"
                     << endl;
              MayDay::Error();
            }

          bool foundExact = false;
          for (int i=0; i<numExact; i++)
            {
              if (errorVars[errVarNo] == exactVars[i])
                {
                  foundExact = true;
                }
            } // end loop over exact variables
          if (!foundExact)
            {
              pout() << "errorVar " << errorVars[errVarNo]
                     << " not found in exact solution!"
                     << endl;
              MayDay::Error();
            }

        } // end loop over errorVars
    } // end if errorVars was specified in the inputs file


    Vector<string> errorNames;
    errorNames.resize(numError);

    constructErrorNames(errorNames, errorVars);

    Vector<LevelData<FArrayBox>* > error(computedNumLevels);

    // allocate error -- same domain as computed solution
    for (int level = 0; level < computedNumLevels; level++)
    {
      error[level] = new LevelData<FArrayBox>(computedGrids[level],
                                              numError, ghostVect);
    }

    if (!isSameSize)
    {
      if (verbose)
      {
        pout () << "compute AMR error..." << endl;
      }

      computeAMRError(error,
                      errorVars,
                      computedSoln,
                      computedVars,
                      computedGrids,
                      computedDx,
                      computedRefRatio,
                      exactSoln,
                      exactVars,
                      exactDx,
                      bogusValue,
                      HOaverage,
                      computeRelativeError);

      if (verbose)
      {
        pout() << "done computing AMR error" << endl;
      }
    }
    else
    {
      // first make sure refRatios are the same
      for (int lev = 0; lev < computedRefRatio.size() - 1; lev++)
      {
        CH_assert(computedRefRatio[lev] == exactRefRatio[lev]);
      }

      CH_assert(exactDx == computedDx);

      if (verbose)
      {
        pout () << "compute sameSize error..." << endl;
      }

      computeSameSizeError(error,
                           errorVars,
                           computedSoln,
                           computedVars,
                           computedGrids,
                           computedDx,
                           computedRefRatio,
                           exactSoln,
                           exactVars,
                           bogusValue,
                           computeRelativeError,
                           doGhostCells);

      if (verbose)
      {
        pout() << "done computing sameSize error" << endl;
      }
    }

    Vector<Real> mean(numError);

    // remove mean
    if (removeMean)
    {
      int lBase = 0;
      int numLevels = error.size();

      for (int err = 0; err < numError; err++)
      {
        Interval errComps(err, err);

        if (verbose)
        {
          pout() << "compute mean for component " << err << endl;
        }

        Real volume;
        mean[err] = computeSum(volume,
                               error,
                               computedRefRatio,
                               computedDx,
                               errComps,
                               lBase);

        mean[err] /= volume;

        if (verbose)
        {
          pout() << "comp = "   << err
                 << ", mean = " << mean
                 << ", vol = "  << volume
                 << endl;
        }

        for (int level = 0; level < numLevels; level++)
        {
          LevelData<FArrayBox>& thisLevelError = *error[level];

          const DisjointBoxLayout levelGrids = error[level]->getBoxes();

          DataIterator dit = levelGrids.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
          {
            const Box thisBox = levelGrids[dit()];
            FArrayBox& thisFabError = thisLevelError[dit()];

            thisFabError.plus(-mean[err],err);
          } // end loop over errors

          thisLevelError.exchange();
        } // end loop over levels
      } // end loop over errors
    }

    // now compute norms

    pout() << "error, step " << nStep << ": L1, L2, Max, sum" << endl;
    for (int err = 0; err < numError; err++)
    {
      Interval errComps(err, err);
      Real L0, L1, L2, sum;

      if (verbose)
      {
        pout() << "compute error for component " << err << endl;
      }

      int lBase = 0;
      int normType = 0;
      Real normDx = computedDx;
      if (useUnitDomain)
        {
          normDx = 1.0/computedDomain.size(0);
        }
      L0 = computeNorm(error,
                       computedRefRatio,
                       normDx,
                       errComps,
                       normType,
                       lBase);

      normType = 1;
      L1 = computeNorm(error,
                       computedRefRatio,
                       normDx,
                       errComps,
                       normType,
                       lBase);

      normType = 2;
      L2 = computeNorm(error,
                       computedRefRatio,
                       normDx,
                       errComps,
                       normType,
                       lBase);
      
      sum = computeSum(error,
                       computedRefRatio,
                       normDx,
                       errComps,
                       lBase);


      pout() << errorNames[err] << ": "
             << L1 << ", "
             << L2 << ", "
             << L0 << ", "
             << sum;

      if (removeMean)
      {
        pout() << " (" << mean[err] << ")";
      }

      pout() << endl;
    } // end loop over errors

    if (doPlots)
    {
      if (verbose)
      {
        pout() << "begin writing hdf5 file..." << endl;
      }

      WriteAMRHierarchyHDF5(errorFile.str(),
                            computedGrids,
                            error,
                            errorNames,
                            computedDomain,
                            computedDx,
                            computedDt,
                            computedTime,
                            computedRefRatio,
                            computedNumLevels);

      if (verbose)
      {
        pout() << "done writing hdf5 file" << endl;
      }
    }

    // clean up memory
    for (int level = 0; level < exactNumLevels; level++)
    {
      if (exactSoln[level] != NULL)
      {
        delete exactSoln[level];
        exactSoln[level] = NULL;
      }
    }

    for (int level = 0; level < computedNumLevels; level++)
    {
      if (computedSoln[level] != NULL)
      {
        delete computedSoln[level];
        computedSoln[level] = NULL;
      }

      if (error[level] != NULL)
      {
        delete error[level];
        error[level] = NULL;
      }
    }
  } // end loop over timesteps

#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main

// this function averages down the fine solution to the valid
// regions of the computed solution, then subtracts ir from
// the computed solution.  (error is exact-computed)
void computeAMRError(Vector<LevelData<FArrayBox>* >&       a_error,
                     const Vector<string>&                 a_errorVars,
                     const Vector<LevelData<FArrayBox>* >& a_computedSoln,
                     const Vector<string>&                 a_computedVars,
                     const Vector<DisjointBoxLayout>&      a_computedGrids,
                     const Real                            a_computedDx,
                     const Vector<int>&                    a_computedRefRatio,
                     const Vector<LevelData<FArrayBox>* >& a_exactSoln,
                     const Vector<string>&                 a_exactVars,
                     const Real                            a_exactDx,
                     Real                                  a_bogus_value,
                     bool                                  a_HOaverage,
                     bool                                  a_computeRelativeError)
{
  int numLevels = a_computedSoln.size();

  CH_assert(a_exactSoln.size() == 1);
  CH_assert(a_error.size() == numLevels);
  CH_assert(a_exactDx <= a_computedDx);
  CH_assert(a_computedRefRatio.size() >= numLevels - 1);

  if (a_exactDx == a_computedDx)
  {
    cerr << "Exact dx and computed dx are equal." << endl;
  }

  // check whether input file selects "sum all variables"
  bool sumAll = false;
  ParmParse pp;
  pp.query("sumAll",sumAll);
  
  // const DisjointBoxLayout& exactGrids = a_exactSoln[0]->getBoxes();

  Real dxLevel = a_computedDx;

  // do a bit of sleight-of-hand in the case where there are no
  // ghost cells in the exact solution -- allocate a temporary which
  // _has_ ghost cells, and do a copyTo
  LevelData<FArrayBox>* exactSolnPtr = NULL;
  bool allocatedMemory = false;
  if (a_exactSoln[0]->ghostVect() == IntVect::Zero)
    {
      exactSolnPtr = new LevelData<FArrayBox>(a_exactSoln[0]->getBoxes(),
                                              a_exactSoln[0]->nComp(),
                                              IntVect::Unit);
      a_exactSoln[0]->copyTo(*exactSolnPtr);

      allocatedMemory = true;
     }
   else
     {
       // if there are ghost cells, we can use the exactSoln as-is
       exactSolnPtr = a_exactSoln[0];
     }
   LevelData<FArrayBox>& exactSolnRef = *exactSolnPtr;

   // first need to set boundary conditions on exactsoln
   // this is for the Laplacian which is needed in AverageHO
   DataIterator ditFine = exactSolnRef.dataIterator();
   DomainGhostBC exactBC;
   Interval exactComps(0, a_exactVars.size() - 1);
   for (int dir = 0; dir < SpaceDim; dir++)
   {
     SideIterator sit;
     for (sit.reset(); sit.ok(); ++sit)
     {
       // use HO extrapolation at physical boundaries
       HOExtrapBC thisBC(dir, sit(), exactComps);
       exactBC.setBoxGhostBC(thisBC);
     }
   }

   for (ditFine.begin(); ditFine.ok(); ++ditFine)
   {
     FArrayBox& thisFineSoln = exactSolnRef[ditFine()];
     const Box& fineBox = exactSolnRef.getBoxes()[ditFine()];
     exactBC.applyInhomogeneousBCs(thisFineSoln, fineBox, a_exactDx);
   }
   exactSolnRef.exchange(exactComps);

   // outer loop is over levels
   for (int level = 0; level < numLevels; level++)
   {
     LevelData<FArrayBox>& thisLevelError = *a_error[level];
     LevelData<FArrayBox>& thisLevelComputed = *a_computedSoln[level];

     // compute refinement ratio between solution at this level
     // and exact solution
     Real nRefTemp = (dxLevel / a_exactDx);
     int nRefExact = (int) nRefTemp;

     // this is to do rounding properly if necessary
     if (nRefTemp - nRefExact > 0.5) nRefExact += 1;

     // make sure it's not zero
     if (nRefExact == 0) nRefExact =1;

     const DisjointBoxLayout levelGrids = a_error[level]->getBoxes();
     const DisjointBoxLayout fineGrids = a_exactSoln[0]->getBoxes();
     DisjointBoxLayout coarsenedFineGrids;

     // petermc, 14 Jan 2014: Replace this because fineGrids might
     // not be coarsenable by nRefExact.
     // coarsen(coarsenedFineGrids, fineGrids, nRefExact);
     int nCoarsenExact = nRefExact;
     while ( !fineGrids.coarsenable(nCoarsenExact) && (nCoarsenExact > 0) )
       {
         // Divide nCoarsenExact by 2 until fineGrids is coarsenable by it.
         nCoarsenExact /= 2;
       }
     if (nCoarsenExact == 0)
       {
         nCoarsenExact = 1;
       }
     coarsen(coarsenedFineGrids, fineGrids, nCoarsenExact);

     int numExact = a_exactVars.size();
     LevelData<FArrayBox> averagedExact(coarsenedFineGrids, numExact);

     Box fineRefBox(IntVect::Zero, (nCoarsenExact-1)*IntVect::Unit);

     // average fine solution down to coarsened-fine level
     // loop over grids and do HO averaging down
     //DataIterator crseExactDit = coarsenedFineGrids.dataIterator();
     for (ditFine.reset(); ditFine.ok(); ++ditFine)
       {
         const Box fineBox = exactSolnRef.getBoxes()[ditFine()];
         FArrayBox fineTemp(fineBox, 1);
         Box coarsenedFineBox(fineBox);
         coarsenedFineBox.coarsen(nCoarsenExact);
         if (a_exactDx < a_computedDx)
           {
             // loop over components
             for (int comp = 0; comp < numExact; comp++)
               {
                 Box coarseBox(coarsenedFineGrids.get(ditFine()));
                 coarseBox &= coarsenedFineBox;
                 
                 if (!coarseBox.isEmpty())
                   {
                     // for now, this is a quick and dirty way to avoid
                     // stepping out of bounds if there are no ghost cells.
                     // LapBox will be the box over which we are
                     // able to compute the Laplacian.
                     Box LapBox = exactSolnRef[ditFine()].box();
                     LapBox.grow(-1);
                     LapBox &= fineBox;
                     fineTemp.setVal(0.0);
                     int doHO = 0;
                     if (a_HOaverage)
                       { 
                         doHO = 1;
                       }
                     
                     // average by default
                     int doAverage = 1;
                     
                     if (sumAll || sumVar(a_exactVars[comp]))
                       {
                         doAverage = 0;
                       }
                     
                     // average or sum, based on booleans
                     FORT_AVERAGEHO(CHF_FRA1(averagedExact[ditFine], comp),
                                    CHF_CONST_FRA1(exactSolnRef[ditFine], comp),
                                    CHF_FRA1(fineTemp, 0),
                                    CHF_BOX(coarseBox),
                                    CHF_BOX(LapBox),
                                    CHF_CONST_INT(nCoarsenExact),
                                    CHF_BOX(fineRefBox),
                                    CHF_INT(doHO),
                                    CHF_INT(doAverage));
                   } // end if crseBox not empty
               } // end loop over comps
           }
         else
           {
             // if cell sizes are the same, then copy
             averagedExact[ditFine].copy(exactSolnRef[ditFine]);
           }
         
       } // end loop over exact solution boxes
     
     int nRefineComputed = nRefExact / nCoarsenExact;
     LevelData<FArrayBox>* thisLevelComputedRefinedPtr = &thisLevelComputed;
     LevelData<FArrayBox>* thisLevelErrorRefinedPtr = &thisLevelError;
     if (nRefineComputed > 1)
       {
         // Do piecewise constant interpolation (replication) by nRefineComputed
         // on thisLevelComputed.
         DisjointBoxLayout levelRefinedGrids;
         refine(levelRefinedGrids, levelGrids, nRefineComputed);
         int nCompComputed = thisLevelComputed.nComp();
         IntVect ghostVectComputed = nRefineComputed * thisLevelComputed.ghostVect();
         thisLevelComputedRefinedPtr =
           new LevelData<FArrayBox>(levelRefinedGrids, nCompComputed, ghostVectComputed);
         ProblemDomain levelDomain = levelRefinedGrids.physDomain();
         FineInterp interpolator(levelRefinedGrids, nCompComputed, nRefineComputed,
                                 levelDomain);
         interpolator.pwcinterpToFine(*thisLevelComputedRefinedPtr, thisLevelComputed);

         int nCompErr = thisLevelError.nComp();
         IntVect ghostVectErr = nRefineComputed * thisLevelError.ghostVect();
         thisLevelErrorRefinedPtr =
           new LevelData<FArrayBox>(levelRefinedGrids, nCompErr, ghostVectErr);
       }

     // initialize error to 0
     // also initialize error to a bogus value
     DataIterator levelDit = thisLevelError.dataIterator();
     for (levelDit.begin(); levelDit.ok(); ++levelDit)
       {
         (*thisLevelErrorRefinedPtr)[levelDit].setVal(a_bogus_value);
       }

     Box refComputedBox(IntVect::Zero, (nRefineComputed-1)*IntVect::Unit);
     // loop over variables
     for (int nErr = 0; nErr < a_errorVars.size(); nErr++)
       {
         string thisErrVar = a_errorVars[nErr];
         bool done = false;

         // first loop over exact variables
         for (int exactComp = 0; exactComp < a_exactVars.size(); exactComp++)
           {
             string thisExactVar = a_exactVars[exactComp];
             // check if this exact variable is "the one"
             if ((thisExactVar == thisErrVar) || nonAverageVar(thisErrVar))
               {
                 int computedComp = 0;
                 // now loop over computed variables
                 while (!done && (computedComp < a_computedVars.size()))
                   {
                     if (a_computedVars[computedComp] == thisErrVar)
                       {
                         if (!nonAverageVar(thisErrVar))
                           {
                             // copy averaged exact solution -> error
                             // and then subtract computed solution
                             Interval exactInterval(exactComp, exactComp);
                             Interval errorInterval(nErr, nErr);
                             averagedExact.copyTo(exactInterval, *thisLevelErrorRefinedPtr,
                                                  errorInterval);
                           }
                         
                         DataIterator levelDit = thisLevelError.dataIterator();
                         for (levelDit.reset(); levelDit.ok(); ++levelDit)
                           {
                             FArrayBox& thisComputedRefined = (*thisLevelComputedRefinedPtr)[levelDit];
                             FArrayBox& thisErrorRefined = (*thisLevelErrorRefinedPtr)[levelDit];
                             if (a_computeRelativeError)
                               {
                                 // do this a little strangely -- relative
                                 // error is one - computed/exact.
                                 thisErrorRefined.divide(thisComputedRefined, computedComp, nErr, 1);
                                 thisErrorRefined.invert(-1.0, nErr, 1);
                                 thisErrorRefined.plus(1.0, nErr, 1);
                               }
                             else
                               {
                                 thisErrorRefined.minus(thisComputedRefined, computedComp, nErr, 1);
                               }
                             if (nRefineComputed > 1)
                               {
                                 FArrayBox& thisError = thisLevelError[levelDit];
                                 Box coarseBox = thisError.box();
                                 // Average thisErrorRefined to thisError.
                                 int doHO = 0;
                                 if (a_HOaverage)
                                   { 
                                     doHO = 1;
                                   }
                                 CH_assert(doHO == 0);

                                 // for now, this is a quick and dirty way to avoid
                                 // stepping out of bounds if there are no ghost cells.
                                 // LapBox will be the box over which we are
                                 // able to compute the Laplacian.
                                 Box LapBox = thisErrorRefined.box();
                                 LapBox.grow(-1);
                                 // LapBox &= fineBox;
                                 FArrayBox fineTemp(thisErrorRefined.box(), 1);
                                 fineTemp.setVal(0.0);
                                 
                                 // average by default
                                 int doAverage = 1;
                                 // average or sum, based on booleans
                                 FORT_AVERAGEHO(CHF_FRA1(thisError, nErr),
                                                CHF_CONST_FRA1(thisErrorRefined, nErr),
                                                CHF_FRA1(fineTemp, 0),
                                                CHF_BOX(coarseBox),
                                                CHF_BOX(LapBox),
                                                CHF_CONST_INT(nRefineComputed),
                                                CHF_BOX(refComputedBox),
                                                CHF_INT(doHO),
                                                CHF_INT(doAverage));
                                 
                               }
                           } // end loop over coarse grids
                         
                         done = true;
                       } // if computedVar is a_errorVar
                     
                     computedComp += 1;
                   } // end loop over a_computedVars
                 
                 if (!done)
                   {
                     pout() << "Variable " << thisErrVar  << " not found!!!" << endl;
                     MayDay::Error();
                   }
               } // end if this exactVar is correct
           } // end loop over exact variables
       } // end loop over errors
     
     if (nRefineComputed > 1)
       {
         delete thisLevelComputedRefinedPtr;
         delete thisLevelErrorRefinedPtr;
       }
     
     // now need to set covered regions to 0
     if (level < numLevels - 1)
       {
         // will need to loop over all boxes in finer level, not just
         // those on this processor...
         const BoxLayout& finerGrids = a_computedSoln[level + 1]->boxLayout();
         LayoutIterator fineLit = finerGrids.layoutIterator();
         
         // outer loop over this level's grids, since there are fewer of them
         DataIterator levelDit = thisLevelError.dataIterator();
         for (levelDit.reset(); levelDit.ok(); ++levelDit)
           {
             const Box& coarseBox = levelGrids[levelDit()];
             FArrayBox& thisError = thisLevelError[levelDit()];
             int numError = thisError.nComp();
             
             for (fineLit.reset(); fineLit.ok(); ++fineLit)
               {
                 Box fineBox(finerGrids[fineLit()]);
                 // now coarsen box down to this level
                 fineBox.coarsen(a_computedRefRatio[level]);
                 // if coarsened fine box intersects error's box, set
                 // overlap to 0
                 fineBox &= coarseBox;
                 if (!fineBox.isEmpty())
                   {
                     thisError.setVal(0.0, fineBox, 0, numError);
                   }
               } // end loop over finer-level grids
           } // end loop over this-level grids
         
         // this is a good place to update dx as well
         dxLevel = dxLevel / a_computedRefRatio[level];
       } // end if there is a finer level
     
     thisLevelError.exchange();
   } // end loop over levels
   
   // clean up if we need to
   if (allocatedMemory)
     {
       delete exactSolnPtr;
       exactSolnPtr = NULL;
     }
}


// this function works on two solutions on equivalent grids.
// It subtracts the computed solution from the exact solution
// if a_doGhostCells == true, then does box-by-box comparison,
// including ghost cells (boxes must be the same for each).
// Otherwise, only does this for valid cells, but boxes don't
// need to be the same.
void computeSameSizeError(Vector<LevelData<FArrayBox>* >&       a_error,
                          const Vector<string>&                 a_errorVars,
                          const Vector<LevelData<FArrayBox>* >& a_computedSoln,
                          const Vector<string>&                 a_computedVars,
                          const Vector<DisjointBoxLayout>&      a_computedGrids,
                          const Real                            a_dx,
                          const Vector<int>&                    a_refRatio,
                          const Vector<LevelData<FArrayBox>* >& a_exactSoln,
                          const Vector<string>&                 a_exactVars,
                          Real                                  a_bogus_value,
                          bool                                  a_computeRelativeError,
                          bool                                  a_doGhostCells)

{
  int numLevels = a_computedSoln.size();

  CH_assert(a_exactSoln.size() == numLevels);
  CH_assert(a_error.size() == numLevels);
  CH_assert(a_refRatio.size() >= numLevels - 1);

  Real dxLevel = a_dx;

  // outer loop is over levels
  for (int level = 0; level < numLevels; level++)
  {
    LevelData<FArrayBox>& thisLevelError = *a_error[level];
    LevelData<FArrayBox>& thisLevelComputed = *a_computedSoln[level];
    LevelData<FArrayBox>& thisLevelExact = *a_exactSoln[level];

    const DisjointBoxLayout levelGrids = thisLevelComputed.getBoxes();
    const DisjointBoxLayout exactGrids = thisLevelExact.getBoxes();

    DataIterator levelDit = levelGrids.dataIterator();
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      // initialize error to a bogus value
      thisLevelError[levelDit()].setVal(a_bogus_value);
    }

    // loop over variables
    for (int nErr = 0; nErr < a_errorVars.size(); nErr++)
    {
      string thisErrVar = a_errorVars[nErr];
      bool done = false;

      // this is where things differ between the ghost-cell
      // and non-ghost-cell approach.
      if (a_doGhostCells)
      {
        // this is the older approach to things --
        // do everything grid-by-grid

        // first loop over exact variables
        for (int exactComp = 0; exactComp < a_exactVars.size(); exactComp++)
        {
          string thisExactVar = a_exactVars[exactComp];
          // check if this exact variable is "the one"
          if (thisExactVar == thisErrVar)
          {
            int computedComp = 0;

            // now loop over computed variables
            while (!done && (computedComp < a_computedVars.size()))
            {
              if (a_computedVars[computedComp] == thisErrVar)
              {
                // copy exact solution -> error
                // and then subtract computed solution
                DataIterator exactDit = thisLevelExact.dataIterator();
                for (levelDit.reset(); levelDit.ok(); ++levelDit)
                {
                  FArrayBox& thisComputed = thisLevelComputed[levelDit()];
                  FArrayBox& thisError = thisLevelError[levelDit()];
                  const Box& thisBox = levelGrids[levelDit()];

                  for (exactDit.begin(); exactDit.ok(); ++exactDit)
                  {
                    if (thisBox.contains(exactGrids[exactDit()]))
                    {
                      thisError.copy(thisLevelExact[exactDit()],
                                     exactComp, nErr, 1);
                    } // end if exact and computed boxes match
                  } // end loop over exact grids

                  if (a_computeRelativeError)
                  {
                    // do this a little strangely -- relative
                    // error is one - computed/exact.
                    thisError.divide(thisComputed, computedComp, nErr, 1);
                    thisError.invert(-1.0, nErr, 1);
                    thisError.plus(1.0, nErr, 1);
                  }
                  else
                  {
                    thisError.minus(thisComputed, computedComp, nErr, 1);
                  }
                } // end loop over grids

                done = true;
              } // if a_computedVar is a_errorVar

              computedComp += 1;
            } // end loop over a_computedVars

            if (!done)
            {
              pout() << "Variable " << thisErrVar  << " not found!!!" << endl;
              MayDay::Error();
            }
          } // end if this exactVar is correct
        }  // end loop over exact variables
      }
      else
        // non-ghost cell case; this is simpler:
      {
        // first loop over exact variables and copy into error
        for (int exactComp=0; exactComp<a_exactVars.size(); ++exactComp)
        {
          string thisExactVar = a_exactVars[exactComp];

          // check if this exact variable is "the one"
          if (thisExactVar == thisErrVar)
          {
            // copy exact solution -> error
            Interval exactInterval(exactComp, exactComp);
            Interval errInterval(nErr, nErr);
            thisLevelExact.copyTo(exactInterval,
                                  thisLevelError,
                                  errInterval);
            done = true;
          } // end if this exact var is the error var
        } // end loop over exact comps

        if (!done)
        {
          pout() << "Variable " << thisErrVar
                 << " not found in exact solution!!!" << endl;
          MayDay::Error();
        }

        done = false;
        int computedComp = 0;
        // now loop over computed variables and subtract computed solution
        while (!done && (computedComp < a_computedVars.size()))
        {
          if (a_computedVars[computedComp] == thisErrVar)
          {
            for (levelDit.reset(); levelDit.ok(); ++levelDit)
            {
              FArrayBox& thisComputed = thisLevelComputed[levelDit()];
              FArrayBox& thisError = thisLevelError[levelDit()];

              thisError.minus(thisComputed, computedComp, nErr, 1);
            } // end loop over computed/error grids

            done = true;
          } // if a_computedVar is a_errorVar

          computedComp += 1;
        } // end loop over a_computedVars

        if (!done)
        {
          pout() << "Variable " << thisErrVar  << " not found!!!" << endl;
          MayDay::Error();
        }
      } // end non-ghost-cell case
    } // end loop over errors

    // now need to set covered regions to 0
    if (level < numLevels - 1)
    {
      // will need to loop over all boxes in finer level, not just
      // those on this processor...
      const BoxLayout& finerGrids = a_computedSoln[level + 1]->boxLayout();
      LayoutIterator fineLit = finerGrids.layoutIterator();

      // outer loop over this level's grids, since there are fewer of them
      DataIterator levelDit = thisLevelError.dataIterator();
      for (levelDit.reset(); levelDit.ok(); ++levelDit)
      {
        const Box& coarseBox = levelGrids[levelDit()];
        FArrayBox& thisError = thisLevelError[levelDit()];
        int numError = thisError.nComp();

        for (fineLit.reset(); fineLit.ok(); ++fineLit)
        {
          Box fineBox(finerGrids[fineLit()]);
          // now coarsen box down to this level
          fineBox.coarsen(a_refRatio[level]);
          // if coarsened fine box intersects error's box, set
          // overlap to 0
          fineBox &= coarseBox;
          if (!fineBox.isEmpty())
          {
            thisError.setVal(0.0, fineBox, 0, numError);
          }
        } // end loop over finer-level grids
      } // end loop over this-level grids

      // this is a good place to update dx as well
      dxLevel = dxLevel / a_refRatio[level];
    } // end if there is a finer level

    // finally, if we're not doing ghost cells, do an exchange just
    // to "prettify" the output
    if (!a_doGhostCells)
    {
      thisLevelError.exchange(thisLevelError.interval());
    }
  } // end loop over levels
}

void constructErrorNames(Vector<string>&       a_errorNames,
                         const Vector<string>& a_errorVars)
{
  CH_assert(a_errorNames.size() == a_errorVars.size());

  // for now, don't do anything fancy -- just copy
  for (int i = 0; i < a_errorVars.size(); i++)
  {
    a_errorNames[i] = a_errorVars[i];
  }
}

void init(string&         a_exactRoot,
          string&         a_computedRoot,
          string&         a_errorRoot,
          int&            a_intFieldSize,
          int&            a_numCrseStart,
          Vector<string>& a_errorVars,
          int&            a_numCrseFinish,
          int&            a_crseStep,
          int&            a_crseMult,
          bool&           a_doPlots,
          bool&           a_isTimeDep,
          bool&           a_isSameSize,
          Real&           a_bogusValue,
          bool&           a_computeRelativeError,
          bool&           a_removeMean,
          bool&           a_doGhostCells,
          bool&           a_useUnitDomain,
          bool&           a_HOaverage,
          bool&           a_verbose)
{
  ParmParse ppCompare("compare");

  int temp = a_doPlots;
  ppCompare.query("doPlots", temp);
  a_doPlots = (temp == 1);

  temp = a_verbose;
  ppCompare.query("verbose", temp);
  a_verbose = (temp == 1);

  ppCompare.query("bogus_value", a_bogusValue);

  temp = a_HOaverage;
  ppCompare.query("HOaverage", temp);
  a_HOaverage = (temp == 1);

  temp = a_computeRelativeError;
  ppCompare.query("computeRelativeError", temp);
  a_computeRelativeError = (temp == 1);

  temp = a_removeMean;
  ppCompare.query("removeMean", temp);
  a_removeMean = (temp == 1);

  temp = a_isSameSize;
  ppCompare.query("sameSize", temp);
  a_isSameSize = (temp == 1);

  // ghost cells thing only really matters for sameSize
  if (a_isSameSize)
  {
    temp = a_doGhostCells;
    ppCompare.query("doGhostCells", temp);
    a_doGhostCells = (temp == 1);
  }

  temp = a_useUnitDomain;
  ppCompare.query("useUnitDomain", temp);
  a_useUnitDomain = (temp == 1);

  temp = a_isTimeDep;
  ppCompare.query("isTimeDep", temp);
  a_isTimeDep = (temp == 1);

  ppCompare.get("exactRoot", a_exactRoot);
  ppCompare.get("computedRoot", a_computedRoot);
  if (a_doPlots)
  {
    ppCompare.get("errorRoot", a_errorRoot);
  }

  a_numCrseStart = 0;
  ppCompare.query("numCrseStart", a_numCrseStart);
  pout() << "numCrseStart " << a_numCrseStart << endl;

  if (a_isTimeDep)
  {
    ppCompare.get("numCrseFinish", a_numCrseFinish);
    a_crseStep = 1;
    ppCompare.query("crseStep", a_crseStep);
    pout() << "numCrseFinish " << a_numCrseFinish <<
      " crseStep " << a_crseStep << endl;
    ppCompare.get("mult", a_crseMult);
    pout() << "mult = " << a_crseMult << endl;
    a_intFieldSize = 4;
    ppCompare.query("intFieldWidth", a_intFieldSize);
  }

  int nErr = ppCompare.countval("error_var");
  if (nErr > 0)
  {
    a_errorVars.resize(nErr);
    ppCompare.getarr("error_var", a_errorVars, 0, nErr);
  }

  nErr = ppCompare.countval("no_average_var");
  noAverageVars.resize(nErr);
  if (nErr > 0)
  {
    ppCompare.getarr("no_average_var", noAverageVars, 0, nErr);
  }

  nErr = ppCompare.countval("sum_var");
  sumVars.resize(nErr);
  if (nErr > 0)
  {
    ppCompare.getarr("sum_var", sumVars, 0, nErr);
  }

}

void constructPlotFileName(ostrstream&   a_fileName,
                           const string& a_fileRoot,
                           const int     a_intFieldWidth,
                           const int     a_step)
{
  // this is kinda klugy, but what are ya gonna do?
  a_fileName << a_fileRoot
             << setw(a_intFieldWidth) << a_step
             << "."
             << setw(1) << CH_SPACEDIM
             << "d.hdf5" << ends;
}


// this function returns true if a_errorName is on the list of
// error variables which are not compared against an averaged fine
// solution
bool nonAverageVar(const string& a_errorName)
{
  bool match = false;
  for (int n = 0; n < noAverageVars.size(); n++)
  {
    if (a_errorName == noAverageVars[n])
    {
      match = true;
    }
  }

  return match;
}


// this function returns true if a_errorName is on the list of
// error variables which are compared against sums of the fine-level
// solution rather than averages of the fine solution
bool sumVar(const string& a_errorName)
{
  bool match = false;
  for (int n = 0; n < sumVars.size(); n++)
  {
    if (a_errorName == sumVars[n])
    {
      match = true;
    }
  }

  return match;
}
