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
using std::ifstream;
using std::ios;

#include "LevelFluxRegisterEdge.H"
#include "parstream.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "CoarseAverageFace.H"
#include "EdgeDataBox.H"
#include "AMRIO.H"
#include "BoxIterator.H"
#include "computeNorm.H"
#include "CONSTANTS.H"

#include "UsingNamespace.H"

enum rhsTypes
{
  quadratic = 0,
  sinusoidal,
  constantX,
  NUM_RHS_TYPES
};

int phiType = sinusoidal;
//int phiType = constantX;

#if defined(CH_USE_DOUBLE)
Real errorEps = 1.0e-10;
#elif defined(CH_USE_FLOAT)
// need much-relaxed tolerances for single-precision
Real errorEps = 5.0e-3;
#endif

bool writePlots = true;

void
initializeRefluxTest(FArrayBox& a_edgeData,
                     const Real a_dx,
                     const int a_dir)
{
  // this replaces the fortran subroutine in refluxEdgeTest.ChF

  RealVect offSet(RealVect::Zero);
  offSet[a_dir] = 0.5;

  // shift sinusoid so that it's not constant/zero at periodic boundaries
  RealVect xShift(RealVect::Unit);
  xShift *= 0.10;

  BoxIterator bit(a_edgeData.box());
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      for (int var=0; var<a_edgeData.nComp(); var++)
        {
          RealVect x(iv);
          x += offSet;
          x *= a_dx;

          Real thisPhi;
          if (phiType == quadratic)
            {
              thisPhi = x[1]*(x[1]-1.0)*x[0]*x[0];
            }
          else if (phiType == sinusoidal)
            {
              x += xShift;
              thisPhi = 1.0+D_TERM6(sin(2*Pi*x[0]*2.0),
                                    *cos(2*Pi*x[1]*2.0),
                                    *sin(2*Pi*x[2]*2.0),
                                    *sin(2*Pi*x[3]*2.0),
                                    *cos(2*Pi*x[4]*2.0),
                                    *sin(2*Pi*x[5]*2.0));

            }
          else if (phiType == constantX)
            {
              Real constVal = 1.0;
              if (a_dir == 0)
                {
                  thisPhi = constVal;
                }
              else
                {
                  thisPhi = 0.0;
                }
            }
          else
            {
              MayDay::Error("bad phiType");
            }

          /*
               thisPhi = thisPhi + y**2;
c               thisPhi = thisPhi + z**2]
c  this puts a flux in the x-direction
c               thisPhi = (1-CHF_ID(dir,1))*(1-CHF_ID(dir,2))*thisPhi
c   this puts a flux in the y-direction
               thisPhi = (1-CHF_ID(dir,0))*(1-CHF_ID(dir,2))*thisPhi
c   this puts a flux in the z-direction
c               thisPhi = (1-CHF_ID(dir,0))*(1-CHF_ID(dir,1))*thisPhi

          */

          a_edgeData(iv, var) = thisPhi;
        } // end loop over variables
    } // end loop over box

}

//  computes curl of face-centered data
void computeCurlEdge(FArrayBox& a_curl,
                     const FArrayBox& a_data,
                     const Box& a_gridBox,
                     const Real a_dx,
                     const int a_curlComp,
                     Real a_edgeComp)
{

  Real oneOnDx = 1.0/a_dx;
  int dir1;
  int dir2;

  if ((SpaceDim == 2) || (a_curlComp == 2))
    {
      dir1 = 0;
      dir2 = 1;
    }
  else if (a_curlComp == 0)
    {
      dir1= 1;
      dir2= 2;
    }
  else if (a_curlComp == 1)
    {
      dir1=2;
      dir2=0;
    }
  else
    {
      MayDay::Error("ComputeCurl -- bad direction");
    }


  // ::: since this is edge-centered, difference is between
  // ::: i+1 and i
  IntVect ii1(IntVect::Zero);
  ii1[dir1] = 1;

  IntVect ii2(IntVect::Zero);
  ii2[dir2] = 1;

  //Real oneOnDx = 1.0/a_dx;

  BoxIterator bit(a_gridBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      Real thisCurl;

      if (a_edgeComp == dir1)
        {
          IntVect iv2 = iv + ii2;
          thisCurl =  oneOnDx*(a_data(iv2, 0) - a_data(iv,0));

          a_curl(iv) -= thisCurl;
        }
      else if (a_edgeComp == dir2)
        {
          IntVect iv1 = iv + ii1;
          thisCurl = oneOnDx*(a_data(iv1, 0) - a_data(iv,0));
          a_curl(iv,0) += thisCurl;
        }
    } // end loop over cells
}




// ---------------------------------------------------------
// helper functions (lifted from Divergence.cpp in order to make this test
// self-sufficient)

void
singleBoxDivergence(const FArrayBox& uEdgeDir,
                    FArrayBox& a_div,
                    const Box& cellBox,
                    Real a_dx,
                    int a_dir)
{

  IntVect offset = IntVect::Zero;
  offset[a_dir] = 1;

  Real oneOnDx = 1.0/a_dx;

  for (int n=0; n<a_div.nComp(); n++)
    {
      // this is going to be really slow, but it's a test, so
      // it shouldn't be that bad...
      BoxIterator cellBit(cellBox);
      for (cellBit.begin(); cellBit.ok(); ++cellBit)
        {
          IntVect cellIV = cellBit();
          IntVect loFace = cellIV;
          IntVect hiFace = loFace + offset;

          a_div(cellIV, n) += oneOnDx*(uEdgeDir(hiFace,n) - uEdgeDir(loFace,n));

        } // end loop over cells
    } // end loop over components


}


void levelDivergenceMAC(LevelData<FArrayBox>& a_div,
                        const LevelData<FluxBox>& a_uEdge,
                        const Real a_dx)
{

  // silly way to do this until i figure out a better
  // way to make this dimensionally-independent
  CH_assert (a_uEdge.nComp() >= a_div.nComp());


  DataIterator dit = a_div.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      a_div[dit()].setVal(0.0);

      const FluxBox& thisFluxBox = a_uEdge[dit()];
      Box cellBox(thisFluxBox.box());
      // just to be sure we don't accidentally trash memory...
      cellBox &= a_div[dit()].box();

      // now loop over coordinate directions and add to divergence
      for (int dir=0; dir<SpaceDim; dir++)
        {
          const FArrayBox& uEdgeDir = thisFluxBox[dir];

          // this would normally be the Divergence fortran

          singleBoxDivergence(uEdgeDir,
                              a_div[dit()],
                              cellBox,
                              a_dx,
                              dir);
        }

    }

}



int
main(int argc, char* argv[])
{

  int status = 0;

#ifdef CH_MPI
  // Start MPI
  MPI_Init(&argc,&argv);
#ifdef CH_AIX
  //  H5dont_atexit();
#endif
  // setChomboMPIErrorHandler();
#endif

  int rank, number_procs;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
  rank = 0;
  number_procs = 1;
#endif

  if (rank == 0)
  {
    pout() << " number_procs = " << number_procs << endl;
  }

  // for now, this test only works in 3D
  if (SpaceDim == 3)
    {
      // just to make things simple
      int verbosity = 3;

      bool read_grids_from_file = false;

      string grids_file = "refluxTestGrids.dat";
      Vector<Vector<Box> > amrBoxes;
      Vector<DisjointBoxLayout> amrGrids;
      Vector<Vector<int> > procIDs;
      // data is edge-centered (as opposed to face centered)
      // this is node-centered in 2d
      Vector<LevelData<EdgeDataBox>* > amrData;
      Vector<LevelData<FArrayBox>* > amrDivergence;
      Vector<LevelData<FluxBox>* > amrCurl;
      Vector<ProblemDomain> amrDomains;
      Vector<int> ref_ratios;
      Vector<Real> amrDx;
      Vector<LevelFluxRegisterEdge*> amrFR;

      int max_level = 1;
      //int max_grid_size = 64;
      int max_grid_size = 16;
      //  int num_cells = 64;
      int num_cells = 32;
      //int num_cells = 16;

      bool is_periodic[SpaceDim];
      for (int dim=0; dim<SpaceDim; dim++)
        is_periodic[dim] = true;

      ProblemDomain prob_domain(IntVect::Zero,
                                IntVect((num_cells-1)*IntVect::Unit),
                                is_periodic);

#ifdef MPI
      MPI_Barrier(Chombo_MPI::comm);
      if (procI() == 0)
        {
#endif
          int numLevels = max_level+1;
          if (read_grids_from_file)
            {
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

              numLevels = in_numLevels;
              amrBoxes.resize(numLevels);
              procIDs.resize(numLevels);
              amrDomains.resize(numLevels);
              amrDomains[0] = prob_domain;
              // for the moment, nref = 2 everywhere
              int nRef = 4;
              ref_ratios.resize(numLevels, nRef);

              // now loop over levels, starting with level 1 and read grids
              // from file
              int ngrid;
              for (int lev=1; lev<numLevels; lev++)
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
                }
            } // end if reading grids from file
          else
            {
              // define grids here
              amrBoxes.resize(numLevels);
              procIDs.resize(numLevels);
              amrDomains.resize(numLevels);
              amrDomains[0] = prob_domain;

              int nRef = 4;
              ref_ratios.resize(numLevels, nRef);

              // for now, a single box offset from the center
              // assume a single refined level
              CH_assert(numLevels == 2);
              for (int lev=1; lev<numLevels; lev++)
                {
                  amrDomains[lev] = amrDomains[lev-1];
                  amrDomains[lev].refine(ref_ratios[lev-1]);

                  //int ngrid = 1;
                  //int ngrid = 2;
                  //int ngrid = 3;
                  //int ngrid = 4;
                  int ngrid = 5;
                  if (verbosity >= 3)
                    {
                      pout() << "level " << lev << " numGrids = " << ngrid << endl;
                      pout() << "Grids: ";
                    }
                  amrBoxes[lev].resize(ngrid);
                  procIDs[lev].resize(ngrid);

                  const Box& domainBox = amrDomains[lev].domainBox();
                  const IntVect& domLo = domainBox.smallEnd();
                  const IntVect& domHi = domainBox.bigEnd();
                  const IntVect& domSize = domainBox.size();

                  IntVect halfDomain = domSize;
                  halfDomain /= 2;
                  IntVect quarterDomain = domSize;
                  quarterDomain /= 4;
                  IntVect eighthDomain = domSize;
                  eighthDomain /= 8;


                  for (int i=0; i<ngrid; i++)
                    {
                      Box bx;
                      //#if 0
                      if (i==2)
                        {
                          IntVect lo = domLo;
                          lo += halfDomain;
                          IntVect hi = lo;
                          hi += quarterDomain;
                          hi -= IntVect::Unit;
                          bx.define(lo,hi);
                        }
                      else if (i == 3)
                        {
                          IntVect lo = domLo;
                          lo += halfDomain;
                          lo[0] -= eighthDomain[0];
                          IntVect hi = lo;
                          hi += eighthDomain;
                          hi -= IntVect::Unit;
                          bx.define(lo,hi);
                        }
                      else if (i == 0)
                        {
                          IntVect lo = domLo;
                          IntVect hi = lo;
                          hi += quarterDomain;
                          hi -= IntVect::Unit;
                          bx.define(lo,hi);

                          // put "do-nothing" shifts in explicitly to quiet
                          // compiler warnings about empty macro arguments
                          D_TERM6(bx.shift(0,0);/* do nothing in the x-dir */,
                                  bx.shift(1,quarterDomain[1]);,
                                  bx.shift(2,quarterDomain[2]);,
                                  bx.shift(3,0);,/* nothing in extra dirs */
                                  bx.shift(4,0);,/* nothing in extra dirs */
                                  bx.shift(5,0);) /* nothing in extra dirs */

                          D_TERM6(bx.shift(0,0);/* do nothing in the x-dir */,
                                  bx.shift(1,eighthDomain[1]);,
                                  bx.shift(2,eighthDomain[2]);,
                                  bx.shift(3,0);,/* nothing in extra dirs */
                                  bx.shift(4,0);,/* nothing in extra dirs */
                                  bx.shift(5,0);)/* nothing in extra dirs */
                         }
                      //else if (i == 3)
                      else if (i == 1)
                        {
                          IntVect lo = domHi;
                          IntVect hi = domHi;
                          lo -= eighthDomain;
                          //lo -= quarterDomain;
                          lo += IntVect::Unit;
                          bx.define(lo,hi);
#if 0
                          D_TERM(bx.shift(0,0);/* do nothing in the x-direction */ ,
                                 bx.shift(1,-3*quarterDomain[1]/2);,
                                 bx.shift(2,-3*quarterDomain[2]/2);)
#endif
                          // no-op shifts are to avoid empty macro
                          //argument compiler warnings
                          D_TERM6(bx.shift(0,0);/* do nothing in the x-dir */ ,
                                  bx.shift(1,-quarterDomain[1]);,
                                  bx.shift(2,-quarterDomain[2]);,
                                  bx.shift(3,0);,/* do nothing in extra dirs */
                                  bx.shift(4,0);,/* do nothing in extra dirs */
                                  bx.shift(5,0);)/* do nothing in extra dirs */

                          D_TERM6(bx.shift(0,0);/* do nothing in the x-dir */ ,
                                  bx.shift(1,-eighthDomain[1]);,
                                  bx.shift(2,-eighthDomain[2]);,
                                  bx.shift(3,0);,/* do nothing in extra dirs */
                                  bx.shift(4,0);,/* do nothing in extra dirs */
                                  bx.shift(5,0);)/* do nothing in extra dirs */

                          //bx.shift(0,-7*eighthDomain[0]);
                          //bx.shift(0,-3*quarterDomain[0]);


                  }
                else if (i == 4)
                  {
                    IntVect lo = domLo;
                    IntVect hi = eighthDomain;
                    hi -= IntVect::Unit;

                    bx.define(lo,hi);
                    D_TERM6(bx.shift(0,0);/* do nothing in the x-dir */,
                            bx.shift(1,7*eighthDomain[1]);,
                            bx.shift(2,7*eighthDomain[2]);,
                            bx.shift(3,0);,/* do nothing in extra dirs */
                            bx.shift(4,0);,/* do nothing in extra dirs */
                            bx.shift(5,0);)/* do nothing in extra dirs */

                  }
                if (verbosity >= 3)
                  {
                    pout() << bx << endl;
                  }
                amrBoxes[lev][i] = bx;
              } // end loop over boxes on this level


          } //end loop over levels
      } // end if we're not reading grids from file

      amrData.resize(numLevels, NULL);
      amrCurl.resize(numLevels,NULL);
      amrDivergence.resize(numLevels, NULL);
      amrGrids.resize(numLevels);

      amrDx.resize(numLevels,-1.0);
      amrFR.resize(numLevels,NULL);

      // check to see if coarsest level needs to be broken up
      domainSplit(prob_domain,amrBoxes[0], max_grid_size,4);

      procIDs[0].resize(amrBoxes[0].size());
      LoadBalance(procIDs[0], amrBoxes[0]);

      amrGrids[0].define(amrBoxes[0], procIDs[0], amrDomains[0]);

      if (verbosity >= 3)
      {
        pout() << "level 0: ";
        for (int n=0; n<amrBoxes[0].size(); n++)
        {
          pout () << amrBoxes[0][n] << endl;
        }
      }

      for (int lev=1; lev<numLevels; lev++)
        {
          // distribute boxes and define disjointBoxLayout
          LoadBalance(procIDs[lev], amrBoxes[lev]);
          amrGrids[lev].define(amrBoxes[lev], procIDs[lev], amrDomains[lev]);
        } // end loop over levels

#ifdef MPI
  } // end if procID = 0
    //broadcast (amrBoxes, uniqueProc(SerialTask::compute));
    MPI_Barrier(Chombo_MPI::comm);
    broadcast (amrGrids, 0);
    MPI_Barrier(Chombo_MPI::comm);
#endif

    // now initialize data
    for (int lev=0; lev<amrGrids.size(); lev++)
    {
      amrData[lev] = new LevelData<EdgeDataBox>(amrGrids[lev], 1);
      amrCurl[lev] = new LevelData<FluxBox>(amrGrids[lev],1, IntVect::Unit);
      amrDivergence[lev] = new LevelData<FArrayBox>(amrGrids[lev],1);
      if (lev>0)
      {
        // initialize flux register
        amrFR[lev] = new LevelFluxRegisterEdge(amrGrids[lev],
                                               amrGrids[lev-1],
                                               amrDomains[lev],
                                               ref_ratios[lev-1],1);
      }


      // loop over grids on this level and initialize edge-centered data
      DataIterator ditLev = amrData[lev]->dataIterator();
      amrDx[lev] = 1.0/amrDomains[lev].domainBox().size(0);
      for (ditLev.reset(); ditLev.ok(); ++ditLev)
      {
        EdgeDataBox& thisData = (*amrData[lev])[ditLev()];
        for (int dir=0; dir<SpaceDim; ++dir)
        {
          initializeRefluxTest(thisData[dir],
                               amrDx[lev],
                               dir);
#if 0
          FArrayBox & thisFab = thisData[dir];
          if (lev == 0)
          {
            thisFab.setVal(0.0);
          }
          else if (dir==0)
          {
            thisFab.setVal(1.0);
          }
          else
          {
            thisFab.setVal(0.0);
          }
#endif
        }
      }

    }


    // compute level-operator curls and their divergences
    for (int lev=0; lev<amrData.size(); lev++)
    {
      const DisjointBoxLayout& levelGrids = amrGrids[lev];
      DataIterator ditLev = amrData[lev]->dataIterator();
      for (ditLev.reset(); ditLev.ok(); ++ditLev)
      {
        FluxBox& thisCurl = (*amrCurl[lev])[ditLev()];
        thisCurl.setVal(0.0);
        const Box& thisGridBox = levelGrids[ditLev()];
        EdgeDataBox& thisData = (*amrData[lev])[ditLev()];
        for (int curlComp=0; curlComp<SpaceDim; curlComp++)
        {
          for (int edgeComp=0; edgeComp<SpaceDim; edgeComp++)
          {
            if (edgeComp != curlComp)
            {
              FArrayBox& thisCurlFab = thisCurl[curlComp];
              FArrayBox& thisDataFab = thisData[edgeComp];
              Box faceBox = thisGridBox;
              faceBox.surroundingNodes(curlComp);
              computeCurlEdge(thisCurlFab,
                              thisDataFab,
                              faceBox,
                              amrDx[lev],
                              curlComp,
                              edgeComp);
            }
          }
        }
      } // end loop over grids for curl computation

      // compute level divergence of curl
      levelDivergenceMAC((*amrDivergence[lev]), (*amrCurl[lev]),
                         amrDx[lev]);
    } // end loop over levels

#ifdef CH_USE_HDF5
    Real bogusValue = 1.0e8;
#endif
    Vector<string> varNames(1,"levelDiv");

    if (writePlots)
      {
#ifdef CH_USE_HDF5
        // dump out level divergence
        string fname1 = "levelDiv.3d.hdf5";
        WriteAMRHierarchyHDF5(fname1, amrGrids, amrDivergence,
                              varNames, amrDomains[0].domainBox(), amrDx[0],
                              bogusValue, bogusValue, ref_ratios,
                              amrDivergence.size());
#endif
      }

    // compute and report norms of divergence
    if (verbosity >= 2)
      {
        Real L1Div, L2Div, maxDiv;
        Interval normInterval(0,0);
        int normType = 1;
        L1Div = computeNorm(amrDivergence,
                            ref_ratios,
                            amrDx[0],
                            normInterval,
                            normType);

        normType = 2;
        L2Div = computeNorm(amrDivergence,
                            ref_ratios,
                            amrDx[0],
                            normInterval,
                            normType);

        normType = 0;
        maxDiv = computeNorm(amrDivergence,
                             ref_ratios,
                             amrDx[0],
                             normInterval,
                             normType);

        pout() << "level-only divergence: L1(div) = " << L1Div
               << ", L2(div) = " << L2Div
               << ", max(div) = " << maxDiv << endl;

      }




    // now load flux registers and average down face-centered values
    for (int fineLev=amrData.size()-1; fineLev>0; fineLev--)
    {
      // first increment with coarse data
      LevelData<EdgeDataBox>& coarseData = *amrData[fineLev-1];
      LevelData<EdgeDataBox>& fineData = *amrData[fineLev];
      Interval dataComps(0,0);

      amrFR[fineLev]->setToZero();

      DataIterator crseDit = coarseData.dataIterator();
      for (crseDit.reset(); crseDit.ok(); ++crseDit)
      {

        //for (int faceDir=0; faceDir<SpaceDim; faceDir++) {
        Real scale = 1.0;
        for (int edgeDir=0; edgeDir<SpaceDim; edgeDir++)
        {
          // if edgeDir and faceDir are the same, there
          // should be appropriate faces
          //if (edgeDir != faceDir) {
          //#if 0
          amrFR[fineLev]->incrementCoarse(coarseData[crseDit()][edgeDir],
                                          scale, crseDit(), dataComps,
                                          dataComps);
          //#endif
          //}
          //}
        }
      }

      //#if 0
      // now increment with fine data
      DataIterator fineDit = fineData.dataIterator();
      for (fineDit.reset(); fineDit.ok(); ++fineDit)
      {
        Real scale = 1.0;

        for (int edgeDir=0; edgeDir<SpaceDim; edgeDir++)
          {

            //#if 0
            amrFR[fineLev]->incrementFine(fineData[fineDit()][edgeDir],
                                          scale, fineDit(),
                                          dataComps, dataComps);
            //#endif
          }
      }



      //#endif
      // now do averaging down
      //CoarseAverageEdge avgDown(amrGrids[fineLev],1, ref_ratios[fineLev-1]);
      CoarseAverageFace avgDown(amrGrids[fineLev],1, ref_ratios[fineLev-1]);

      avgDown.averageToCoarse((*amrCurl[fineLev-1]), *amrCurl[fineLev]);

      // compute divergence before refluxing
      for (int lev=0; lev<amrData.size(); lev++)
      {
        LevelData<FluxBox>* fineDataPtr = NULL;
        Real fineDx = 0;
        if ((lev+1) < amrData.size())
        {
          fineDataPtr = amrCurl[lev+1];
          fineDx = amrDx[lev+1];
        }

        levelDivergenceMAC((*amrDivergence[lev]), (*amrCurl[lev]),
                           amrDx[lev]);
      }

      if (writePlots)
        {
#ifdef CH_USE_HDF5
          string fname2 = "prefluxDiv.3d.hdf5";
          WriteAMRHierarchyHDF5(fname2, amrGrids, amrDivergence,
                                varNames, amrDomains[0].domainBox(), amrDx[0],
                                bogusValue, bogusValue, ref_ratios,
                                amrDivergence.size());
#endif
        }

    // compute and report norms of divergence
    if (verbosity >= 2)
      {
        Real L1Div, L2Div, maxDiv;
        Interval normInterval(0,0);
        int normType = 1;
        L1Div = computeNorm(amrDivergence,
                            ref_ratios,
                            amrDx[0],
                            normInterval,
                            normType);

        normType = 2;
        L2Div = computeNorm(amrDivergence,
                            ref_ratios,
                            amrDx[0],
                            normInterval,
                            normType);

        normType = 0;
        maxDiv = computeNorm(amrDivergence,
                             ref_ratios,
                             amrDx[0],
                             normInterval,
                             normType);

        pout() << "pre-reflux divergence: L1(div) = " << L1Div
               << ", L2(div) = " << L2Div
               << ", max(div) = " << maxDiv << endl;

      }


      // finally, do refluxing
      Real scale = 1.0/amrDx[fineLev-1];

#if 0
      // debugging tool to isolate refluxing effect
      LevelData<FluxBox>& crseCurl = *amrCurl[fineLev-1];
      for (crseDit.reset(); crseDit.ok(); ++crseDit)
      {
        crseCurl[crseDit()].setVal(0.0);
      }
#endif

      bool doReflux = true;
      if (doReflux)
      {
        amrFR[fineLev]->refluxCurl(*amrCurl[fineLev-1], scale);
      }

      // do this again, just to be sure
      avgDown.averageToCoarse((*amrCurl[fineLev-1]), *amrCurl[fineLev]);
    }


    // now that we've done that, compute composite divergence
    // of new composite curl

    for (int lev=0; lev<amrData.size(); lev++)
    {
      LevelData<FluxBox>* fineDataPtr = NULL;
      Real fineDx = 0;
      if ((lev+1) < amrData.size())
      {
        fineDataPtr = amrCurl[lev+1];
        fineDx = amrDx[lev+1];
      }

      /*
      compDivergenceMAC(*amrDivergence[lev], *amrCurl[lev], fineDataPtr,
                        amrDx[lev], &fineDx, ref_ratios[lev],
                        amrDomains[lev]);
      */

      levelDivergenceMAC((*amrDivergence[lev]), (*amrCurl[lev]),
                         amrDx[lev]);
    }

    if (writePlots)
      {
#ifdef CH_USE_HDF5
        string fname2 = "compDiv.3d.hdf5";
        WriteAMRHierarchyHDF5(fname2, amrGrids, amrDivergence,
                              varNames, amrDomains[0].domainBox(), amrDx[0],
                              bogusValue, bogusValue, ref_ratios,
                              amrDivergence.size());
#endif
      }

    // compute and report norms of divergence
    if (verbosity >= 2)
      {
        Real L1Div, L2Div, maxDiv;
        Interval normInterval(0,0);
        int normType = 1;
        L1Div = computeNorm(amrDivergence,
                            ref_ratios,
                            amrDx[0],
                            normInterval,
                            normType);

        normType = 2;
        L2Div = computeNorm(amrDivergence,
                            ref_ratios,
                            amrDx[0],
                            normInterval,
                            normType);

        normType = 0;
        maxDiv = computeNorm(amrDivergence,
                             ref_ratios,
                             amrDx[0],
                             normInterval,
                             normType);

        if (verbosity > 2)
          {
            pout() << "post-reflux divergence: L1(div) = " << L1Div
                   << ", L2(div) = " << L2Div
                   << ", max(div) = " << maxDiv << endl;
          }

        if (L1Div > errorEps)
          {
            status += 1;
          }

        if (L2Div > errorEps)
          {
            status += 1;
          }

        if (maxDiv > errorEps)
          {
            status += 1;
          }
      }


    // finally, just dump out effect of refluxing
    for (int fineLev=amrData.size()-1; fineLev>0; fineLev--)
    {
      // do refluxing
      Real scale = 1.0/amrDx[fineLev-1];


      // debugging tool to isolate refluxing effect
      LevelData<FluxBox>& crseCurl = *amrCurl[fineLev-1];
      DataIterator crseDit= crseCurl.dataIterator();
      for (crseDit.reset(); crseDit.ok(); ++crseDit)
      {
        crseCurl[crseDit()].setVal(0.0);
      }

      amrFR[fineLev]->refluxCurl(*amrCurl[fineLev-1], scale);

    }

    for (int lev=0; lev<amrData.size(); lev++)
    {
      LevelData<FluxBox>* fineDataPtr = NULL;
      Real fineDx = 0;
      if ((lev+1) < amrData.size())
      {
        fineDataPtr = amrCurl[lev+1];
        fineDx = amrDx[lev+1];
      }

      /*      compDivergenceMAC(*amrDivergence[lev], *amrCurl[lev], fineDataPtr,
                        amrDx[lev], &fineDx, ref_ratios[lev],
                        amrDomains[lev]);
      */
      levelDivergenceMAC(*amrDivergence[lev], *amrCurl[lev], amrDx[lev]);
    }

    if (writePlots)
      {
#ifdef CH_USE_HDF5
        string fname3 = "divReflux.3d.hdf5";
        WriteAMRHierarchyHDF5(fname3, amrGrids, amrDivergence,
                              varNames, amrDomains[0].domainBox(), amrDx[0],
                              bogusValue, bogusValue, ref_ratios,
                              amrDivergence.size());
#endif
      }


    // now clear data
    for (int lev=0; lev<amrGrids.size(); lev++)
      {
        delete amrData[lev];
        delete amrCurl[lev];
        delete amrDivergence[lev];
        if (amrFR[lev] != NULL)
          {
            delete amrFR[lev];
          }
      }

    if (verbosity > 0)
      {
        if (status == 0)
          {
            pout() << "refluxEdgeTest -- PASSED all tests!" << endl;
          }
        else
          {
            pout() << "refluxEdgeTest -- FAILED at least one test!" << endl;
          }
      }

    } // end if SpaceDim == 3
  else
    {
      MayDay::Warning("refluxEdgeTest not defined for this dimensionality");
    }
#ifdef CH_MPI
  // End MPI
    MPI_Finalize();
#endif

    return status;
}






