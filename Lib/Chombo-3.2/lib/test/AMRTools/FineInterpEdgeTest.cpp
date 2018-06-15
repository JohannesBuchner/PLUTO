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
#include <fstream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::string;
using std::ifstream;
using std::ios;

#include "FineInterpFace.H"
#include "PiecewiseLinearFillPatchFace.H"
#include "ProblemDomain.H"
#include "parstream.H"
#include "BRMeshRefine.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
    //#define ABRFILES
#ifdef ABRFILES
//#include "AMRBoxesAndRanksIO.H"
#endif

#include "FABView.H"

#include "UsingNamespace.H"

bool use_ABRFile = false;

enum probtypes
{
  linear,
  sinusoidal,
  NUM_TYPES
};

bool dump_plots = true;

int probtype = linear;
//int probtype = sinusoidal;

Real errorEps = 1.0e-10;

void
initializeInterpTest(FArrayBox& edgePhi,
                     Real dx,
                     int dir)
{
  Real Pi = 4.0*atan(1.0);

  RealVect offset = 0.5*RealVect::Unit;
  offset[dir] = 0.0;
  RealVect shift = RealVect::Zero;
  //shift[0] = 0.5;

  for (int var = 0; var < edgePhi.nComp(); var++)
    {
      BoxIterator bit(edgePhi.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc(iv);
          loc += offset;
          loc *= dx;
          loc += shift;

          Real thisPhi;
          if (probtype == linear)
            {
              thisPhi = D_TERM6(loc[0], +loc[1], +loc[2],
                                +loc[3],+loc[4],+loc[5]);
            }

          else if (probtype == sinusoidal)
            {
              loc *= 2*Pi;
              thisPhi = D_TERM6(sin(loc[0]), +sin(loc[1]), +sin(loc[2]),
                                +sin(loc[3]), +sin(loc[4]), +sin(loc[5]));
            }
          edgePhi(iv, var) = thisPhi;

        }
    }

}

int
main(int argc, char* argv[])
{

#ifdef CH_MPI
  // Start MPI
  MPI_Init(&argc,&argv);
#ifdef CH_AIX
  //  H5dont_atexit();
#endif
  // setChomboMPIErrorHandler();
#endif

    // just to make things simple
    int verbosity = 4;

    int status = 0;

    if (verbosity > 1)
      {
        pout() << "beginning FineInterpEdgeTest" << endl;
      }

  // put everything into a local scope, currently
  // this only works for 2D and 3D. (will expand eventually)
  if ((SpaceDim ==2) || (SpaceDim == 3))
  {
    int rank, number_procs;

#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank = 0;
    number_procs = 1;
#endif

    if ((rank == 0) && (verbosity > 2))
      {
        pout() << " number_procs = " << number_procs << endl;
      }


    string grids_file;
    if (SpaceDim == 2)
      {
        grids_file = "interpEdgeTestGrids.2d.dat";
      }
    else if (SpaceDim == 3)
      {
        grids_file = "interpEdgeTestGrids.3d.dat";
      }

    Vector<Vector<Box> > amrBoxes;
    Vector<DisjointBoxLayout> amrGrids;
    Vector<Vector<int> > procIDs;
    Vector<LevelData<FluxBox>* > amrData;
    Vector<LevelData<FluxBox>* > amrError;
    Vector<ProblemDomain> amrDomains;
    Vector<int> ref_ratios;

    bool is_periodic[SpaceDim];
    for (int dir=0; dir<SpaceDim; dir++)
      {
        //is_periodic[dir] = false;
        is_periodic[dir] = true;
      }

#ifdef ABRFILES
    if (use_ABRFile)
      {
        string fname;
        if (SpaceDim == 3)
          {
            fname = "interpTest.3d.abr";
          }
        else if (SpaceDim == 2)
          {
            fname = "interpTest.2d.abr";
          }
        else
          {
            MayDay::Error("undefined abrfile name");
          }

        Box baseDomainBox;
        readABRfile(amrBoxes,
                    procIDs,
                    ref_ratios,
                    baseDomainBox,
                    fname,
                    false);

        // now set things up...
        int numLevels = amrBoxes.size();
        amrGrids.resize(numLevels);
        amrData.resize(numLevels, NULL);
        amrError.resize(numLevels, NULL);
        amrDomains.resize(numLevels);

        amrDomains[0].define(baseDomainBox,
                             is_periodic);

        for (int lev=0; lev<numLevels; lev++)
          {
            if (lev > 0)
              {
                amrDomains[lev] = amrDomains[lev-1];
                amrDomains[lev].refine(ref_ratios[lev-1]);
              }

            bool doLoadBalance = true;
            if (doLoadBalance)
              {
                LoadBalance(procIDs[lev], amrBoxes[lev]);
              }
            amrGrids[lev].define(amrBoxes[lev], procIDs[lev],
                                 amrDomains[lev]);

          }
      }
    else
#endif // end if using ABR files
      {
        int max_level = 5;
        int max_grid_size = 128;
        //int num_cells = 64;
        int num_cells = 32;

        ProblemDomain prob_domain(IntVect::Zero,
                                  IntVect((num_cells-1)*IntVect::Unit),
                                  is_periodic);

#ifdef MPI
        MPI_Barrier(Chombo_MPI::comm);
        if (procID() == 0)
          {
#endif

            ifstream is(grids_file.c_str(), ios::in);
            if (is.fail())
              {
                pout() << "cannot open grids file " << grids_file << endl;
                MayDay::Error();
              }
#ifdef MPI
            if (verbosity >= 3)
              {
                pout() << "procID: " << procID()
                       << "  opening gridfile \n" << endl;
              }
#endif

            // format of file -- number of levels, then for
            // each level starting with level 1, number of
            // grids on level, list of boxes
            int in_numLevels;
            is >> in_numLevels;
            CH_assert (in_numLevels <= max_level+1);
            while (is.get() != '\n');
            amrBoxes.resize(in_numLevels);
            procIDs.resize(in_numLevels);

            amrGrids.resize(in_numLevels);
            amrData.resize(in_numLevels, NULL);
            amrError.resize(in_numLevels, NULL);

            amrDomains.resize(in_numLevels);
            // for the moment, nref = 2 everywhere
            //ref_ratios.resize(in_numLevels, 2);
            ref_ratios.resize(in_numLevels, 4);


            // check to see if coarsest level needs to be broken up
            domainSplit(prob_domain,amrBoxes[0], max_grid_size,4);
            // for now, just put all boxes on proc 0
            procIDs[0].resize(amrBoxes[0].size());
            LoadBalance(procIDs[0], amrBoxes[0]);

            amrDomains[0] = prob_domain;
            amrGrids[0].define(amrBoxes[0], procIDs[0], amrDomains[0]);

            if (verbosity >= 3)
              {
                pout() << "level 0: ";
                for (int n=0; n<amrBoxes[0].size(); n++)
                  {
                    pout () << amrBoxes[0][0] << endl;
                  }
              }
            // now loop over levels, starting with level 1
            int ngrid;
            for (int lev=1; lev<in_numLevels; lev++)
              {
                amrDomains[lev] = amrDomains[lev-1];
                amrDomains[lev].refine(ref_ratios[lev-1]);

                is >> ngrid;
                if (verbosity >= 3)
                  {
                    pout() << "level " << lev << " numGrids = "
                           << ngrid << endl;
                    pout() << "Grids: ";
                  }
                while (is.get() != '\n');
                amrBoxes[lev].resize(ngrid);
                procIDs[lev].resize(ngrid);

                for (int i=0; i<ngrid; i++)
                  {
                    Box bx;
                    is >> bx;
                    while (is.get() != '\n');
                    // quick check on box size
                    Box bxRef(bx);
                    // not really sure why i was doing this (holdover from
                    // legacy code)
                    if (bxRef.longside() > max_grid_size)
                      {
                        pout() << "Grid " << bx << " too large" << endl;
                        MayDay::Error();
                      }
                    if (verbosity >= 3)
                      {
                        pout() << bx << endl;
                      }
                    amrBoxes[lev][i] = bx;
                  } // end loop over boxes on this level


                // distribute boxes and define disjointBoxLayout
                LoadBalance(procIDs[lev], amrBoxes[lev]);
                amrGrids[lev].define(amrBoxes[lev], procIDs[lev], amrDomains[lev]);
              } // end loop over levels

          } // end if old-style grid inputs file

#ifdef MPI
      } // end if procID = 0
    //broadcast (amrBoxes, uniqueProc(SerialTask::compute));
    MPI_Barrier(Chombo_MPI::comm);
    broadcast (amrGrids, 0);
    MPI_Barrier(Chombo_MPI::comm);
#endif
    //int nghost = 4;
    // mirror Francesco by using 8 ghost cells
    int nghost = 8;
    IntVect ghostVect(nghost*IntVect::Unit);

    // now intialize data
    for (int lev=0; lev<amrGrids.size(); lev++)
      {
        amrData[lev] = new LevelData<FluxBox>(amrGrids[lev], 1,ghostVect);
        amrError[lev] = new LevelData<FluxBox>(amrGrids[lev],1,ghostVect);
      }

    // now initialize data on coarsest level
    DataIterator dit0 = amrData[0]->dataIterator();
    Real dx0 = 1.0/amrDomains[0].domainBox().size(0);
    LevelData<FluxBox>& level0Data = *amrData[0];

    for (dit0.reset(); dit0.ok(); ++dit0)
      {
        FluxBox& thisData = level0Data[dit0()];
        for (int dir=0; dir<SpaceDim; dir++)
          {

            FArrayBox& thisFab = thisData[dir];
            initializeInterpTest(thisFab,
                                 dx0,
                                 dir);
          }
      }

    // now loop over levels and interp
    for (int lev=1; lev<amrGrids.size(); lev++)
      {
        DataIterator levelDit = amrData[lev]->dataIterator();
        for (levelDit.begin(); levelDit.ok(); ++levelDit)
          {
            FluxBox& thisFlux = (*amrData[lev])[levelDit];
            for (int dir=0; dir<SpaceDim; dir++)
              {
                Real testVal = 1000+dir;
                thisFlux[dir].setVal(testVal);
              }
          }
      
        FineInterpFace levelInterp(amrGrids[lev],
                                   amrData[lev]->nComp(),
                                   ref_ratios[lev-1],
                                   amrDomains[lev]);


        // turn on fine-data averaging even though it's unnecessary
        // since this is a test (so we'd like to know if it's broken)
        bool averageFineData = true;
        levelInterp.interpToFine(*amrData[lev], *amrData[lev-1], 
                                 averageFineData);

        PiecewiseLinearFillPatchFace filpatcher(amrGrids[lev],
                                                amrGrids[lev-1],
                                                amrData[lev]->nComp(),
                                                amrDomains[lev-1],
                                                ref_ratios[lev-1],
                                                nghost);

        Real timeInterpCoeff = 1;
        filpatcher.fillInterp(*amrData[lev],
                              *amrData[lev-1],
                              *amrData[lev-1],
                              timeInterpCoeff,
                              0,0,amrData[lev]->nComp());

        Interval dataComps = amrData[lev]->interval();
        amrData[lev]->exchange(dataComps);

      }

    // now compute error
    // (no need to compute error on level 0)
    for (int lev=1; lev<amrData.size(); lev++)
      {
        Real dx = 1.0/amrDomains[lev].domainBox().size(0);
        DataIterator levelDit = amrError[lev]->dataIterator();
        LevelData<FluxBox>& levelError = *amrError[lev];
        if (verbosity > 2)
          {
            pout() << "Level = " << lev << endl;
          }
        for (int faceDir=0; faceDir<SpaceDim; faceDir++)
          {
            int levelPatchNo = 0;
            if (verbosity > 2)
              {
                pout () << "    faceDir == " << faceDir << endl;
              }
            for (levelDit.reset(); levelDit.ok(); ++levelDit)
              {
                initializeInterpTest(levelError[levelDit()][faceDir],
                                     dx,
                                     faceDir);

                FArrayBox& thisError = levelError[levelDit()][faceDir];
                FArrayBox& thisData = (*amrData[lev])[levelDit()][faceDir];

                thisError -= thisData;

                if (verbosity > 3)
                  {
                    IntVect testPoint(D_DECL6(128,44,0,0,0,0));
                    if (faceDir == 2 && thisData.box().contains(testPoint))
                      {
                        pout () << "patchNo: " << levelPatchNo
                                << ", testFace value = "
                                << thisData(testPoint,0) << endl;
                      }
                  }

                //ArrayView(&thisError);
                Box domainValidBox = thisError.box();
                domainValidBox &= amrDomains[lev];

                // if linear + periodic doesn't work at edges, so
                // screen out cells which are corrupted
                if (probtype == linear)
                  {
                    if (amrDomains[lev].isPeriodic() )
                      {
                        Box domainBox = amrDomains[lev].domainBox();
                        for (int dir=0; dir<SpaceDim; dir++)
                          {
                            if (amrDomains[lev].isPeriodic(dir))
                              {
                                domainBox.grow(dir,-(ref_ratios[lev-1]));
                              }
                          }
                        domainBox.surroundingNodes(faceDir);
                        domainValidBox &= domainBox;
                      }
                  }


                Real L1Err = thisError.norm(domainValidBox,1,0,1);
                Real MaxErr = thisError.norm(domainValidBox,0,0,1);

                if (abs(L1Err) > errorEps)
                  {
                    if (verbosity > 2)
                      {
                        pout() << "Failed L1(error) test" << endl;
                      }
                    ++status;
                  }

                if (abs(MaxErr) > errorEps)
                  {
                    if (verbosity > 2)
                      {
                        pout() << "Failed Max(error) test" << endl;
                      }
                    ++status;
                  }

                if (verbosity > 2)
                  {
                    pout() << "     patchNo = " << levelPatchNo
                           << " L1 error: " << L1Err
                           << " Max error: " << MaxErr << endl;
                  }
                levelPatchNo++;
              }
          }
      }

    pout() << "exiting" << endl;

    if (verbosity > 1)
      {
        if (status == 0)
          {
            pout() << "FineInterpEdgeTest PASSED all tests" << endl;
          }
        else
          {
            pout() << "FineInterpEdgeTest FAILED at least one test" << endl;            }
      }

    // clean up memory
    for (int lev=0; lev<amrGrids.size(); lev++)
      {
        delete(amrData[lev]);
        delete(amrError[lev]);
      }

  } //end local scope
  else
  {
    MayDay::Warning("FineInterpEdgeTest not implemented for this Dimensionality");
  }


#ifdef CH_MPI
  // End MPI
  MPI_Finalize();
#endif

  return status;
}







