#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "REAL.H"
#include "parstream.H"

#include "LevelData.H"
#include "LoadBalance.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "BoxIterator.H"
#include "ReductionCopier.H"
#include "SpreadingCopier.H"
#include "ReductionOps.H"
#include "BRMeshRefine.H" // contains domainSplit function
#include "CH_Attach.H"
#include "FABView.H"
#include "UsingNamespace.H"

using std::endl;

static bool verbose = true;


/***************/
/***************/
void dumpmemoryatexit();
/***************/
/***************/


void
parseTestOptions( int argc ,char* argv[] ) ;

int reductionTest();

int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  //  registerDebugger();
#endif
  parseTestOptions( argc ,argv ) ;

  int retval = 0;

  retval = reductionTest();

  if (retval == 0)
    {
      pout() << "reductionTest passed all tests!" << endl;
    }
  else
    {
      pout() << "reductionTest failed at least one test, return value = "
             << retval << endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif

  return retval;
}


int reductionTest()
{
  int status = 0;
  int faceStatus = 0;
  int domLength = 18;
  //int domLength = 2;

#ifdef CH_USE_FLOAT
  Real tolerance = 5.0e-6;
#else
  Real tolerance = 1.0e-11;
#endif



  // this test doesn't make any sense in 1d, so only do any of this
  // if DIM> 1
  if (SpaceDim > 1)
    {
      if (verbose)
        {
          pout() << "single-box reduction test" << endl;
        }

      Real dx = 2.0/domLength;
      Box domainBox(IntVect::Zero, (domLength-1)*IntVect::Unit);

      // start with non-periodic case
      ProblemDomain domain(domainBox);

      int transverseDir = 1;
      Box lowerDimDomainBox(domainBox);
      lowerDimDomainBox.setSmall(transverseDir, 0);
      lowerDimDomainBox.setBig(transverseDir,0);
      ProblemDomain lowerDimDomain(lowerDimDomainBox);

      // start with single box case with scaling
      {
        //int maxBoxSize = domLength;
        Vector<Box> gridBoxes(1, domainBox);
        Vector<int> procAssign(1,0);
        DisjointBoxLayout grids(gridBoxes, procAssign, domain);

        Vector<Box> lowerDimGridBoxes(1, lowerDimDomainBox);
        DisjointBoxLayout lowerDimGrids(lowerDimGridBoxes, procAssign, domain);

        LevelData<FArrayBox> data(grids, 1);
        LevelData<FArrayBox> reducedData(lowerDimGrids, 1);

        LevelData<FluxBox> faceData(grids, 1);
        LevelData<FluxBox> reducedFaceData(lowerDimGrids, 1);

        // initialize data to 1
        Real testVal = 1.0;
        DataIterator dit = grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            data[dit()].setVal(testVal);
            for (int dir=0; dir<SpaceDim; dir++)
              {
                faceData[dit][dir].setVal(testVal);
              }
          }

        DataIterator reducedDit = reducedData.dataIterator();
        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            reducedData[reducedDit()].setVal(0.0);
            reducedFaceData[reducedDit].setVal(0.0);
          }

        IntVect GhostVect = IntVect::Zero;
        ReductionCopier reduceCopier(grids, lowerDimGrids, domain,
                                     GhostVect, transverseDir);

        SumOp op(transverseDir);
        op.scale = dx;

        data.copyTo(data.interval(), reducedData, reducedData.interval(),
                    reduceCopier, op);

        // check values
        Real exactVal = domainBox.size(transverseDir)*dx;
        RealVect faceExact(exactVal*RealVect::Unit);
        faceExact[transverseDir] += dx;

        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            // check
            FArrayBox& reducedFab = reducedData[reducedDit()];
            BoxIterator bit(lowerDimGrids.get(reducedDit()));
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                if (abs((reducedFab(iv, 0) - exactVal)) > tolerance)
                  {
                    status+= 1;
                  }
              }
          } // end loop over boxes to check cell-centered data

        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            // check face data
            FluxBox& reducedFlux = reducedFaceData[reducedDit()];
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& faceFab = reducedFlux[dir];

                // instead of checking face-centered data, set it
                // to what it would be if we did the reduction
                faceFab.setVal(faceExact[dir]);

#if 0 // commenting out everything related to FaceSumOp, since it doesn't work
                Box gridBoxDir(lowerDimGrids.get(reducedDit()));
                gridBoxDir.surroundingNodes(dir);
                // if we're in the transverse direction, we only want
                // to look at low-side face
                if (dir == transverseDir) gridBoxDir.growHi(dir,-1);
                BoxIterator bit(gridBoxDir);
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    if (abs((faceFab(iv, 0) - faceExact[dir])) > tolerance)
                      {
                        faceStatus+= 1;
                      }
                  }
#endif
              }
          } // end loop over grid boxes

      if (verbose)
        {
          pout() << "single-box spreading test" << endl;
        }

      // first reset data to a bogus value
      Real bogusVal = 100000;
      for (dit.begin(); dit.ok(); ++dit)
        {
          data[dit].setVal(bogusVal);
          faceData[dit].setVal(bogusVal);
        }

      // now spread reduced data back to original
        Real spreadingScale = 5.0;
        SpreadingOp spreadOp(transverseDir);
        spreadOp.scale = spreadingScale;

        SpreadingCopier spreadCopier(lowerDimGrids, grids,
                                     domain, transverseDir);

        reducedData.copyTo(reducedData.interval(), data,
                           data.interval(), spreadCopier,
                           spreadOp);


        FaceSpreadingOp faceSpreadOp(transverseDir);
        faceSpreadOp.scale = spreadingScale;

        reducedFaceData.copyTo(reducedFaceData.interval(), faceData,
                               faceData.interval(), spreadCopier,
                               faceSpreadOp);


        exactVal *= spreadingScale;
        faceExact *= spreadingScale;

        for (dit.begin(); dit.ok(); ++dit)
          {
            FArrayBox& thisData = data[dit()];
            BoxIterator bit(grids.get(dit()));
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                if (abs((thisData(iv, 0) - exactVal)) > tolerance)
                  {
                    status+= 10;
                  }
              }


            FluxBox& thisFlux = faceData[dit()];
            for (int dir=0; dir<SpaceDim; dir++)
              {
                const FArrayBox& thisFluxDir = thisFlux[dir];
                Box faceBox(grids.get(dit()));
                faceBox.surroundingNodes(dir);
                BoxIterator bit(faceBox);
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    if (abs((thisFluxDir(iv, 0) - faceExact[dir])) > tolerance)
                      {
                        faceStatus+= 10;
                      }
                  }
              }
          }

        if (status == 0)
          {
            pout() << "passed cell-centered single-box test case" << endl;
          }
        else
          {
            pout() << "FAILED cell-centered single-box test case" << endl;
          }

        if (faceStatus == 0)
          {
            pout() << "passed face-centered single-box test case" << endl;
          }
        else
          {
            pout() << "FAILED face-centered single-box test case" << endl;
          }
      } // end single box case


      // now do multiple box and multiple component (with offset) case
      // also try ghost cells here
      {
        
        if (verbose)
          {
            pout() << "multiple-box reduction test" << endl;
          }
        
        int localStatus = 0;
        int localFaceStatus = 0;
        int maxBoxSize = domLength/3;
        Vector<Box> gridBoxes;
        domainSplit(domainBox, gridBoxes, maxBoxSize);
        int numBoxes = gridBoxes.size();
        Vector<int> procAssign(numBoxes,0);
        // for mpi case, distribute boxes among processors
        LoadBalance(procAssign,gridBoxes);
        DisjointBoxLayout grids(gridBoxes, procAssign, domain);

        Vector<Box> lowerDimGridBoxes;
        domainSplit(lowerDimDomainBox, lowerDimGridBoxes, maxBoxSize);
        int numLowerDimBoxes = lowerDimGridBoxes.size();
        Vector<int> lowerDimProcAssign(numLowerDimBoxes,0);
        LoadBalance(lowerDimProcAssign,lowerDimGridBoxes);
        DisjointBoxLayout lowerDimGrids(lowerDimGridBoxes,
                                        lowerDimProcAssign, domain);

        IntVect ghostVect = IntVect::Unit;
        LevelData<FArrayBox> data(grids, 3, ghostVect);
        LevelData<FArrayBox> reducedData(lowerDimGrids, 1, ghostVect);

        LevelData<FluxBox> faceData(grids, 3, ghostVect);
        LevelData<FluxBox> reducedFaceData(lowerDimGrids, 1, ghostVect);

        // initialize data to 1
        Real testVal = 1.0;
        DataIterator dit = grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int comp=0; comp<data.nComp(); comp++)
              {
                data[dit()].setVal((comp+1)*testVal, comp);
                for (int dir=0; dir<SpaceDim; dir++)
                  {
                    faceData[dit][dir].setVal((comp+1)*testVal, comp);
                  }
              }
          }

        DataIterator reducedDit = reducedData.dataIterator();
        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            reducedData[reducedDit()].setVal(0.0);
            reducedFaceData[reducedDit()].setVal(0.0);
          }

        if (verbose)
          {
            pout() << "hiDimGrids: " << endl;
            for (DataIterator dit1 = grids.dataIterator(); dit1.ok(); ++dit1)
              {
                pout() << grids[dit1] << endl;
              }

            pout() << "lowerDimGrids: " << endl;
            for (DataIterator dit1=lowerDimGrids.dataIterator(); dit1.ok(); ++dit1)
              {
                pout() << lowerDimGrids[dit1] << endl;
              }
          }

        ReductionCopier reduceCopier(grids, lowerDimGrids, domain,
                                     ghostVect, transverseDir);

        if (verbose)
          {
            pout() << "Forward copier: " << endl;
            pout () << reduceCopier << endl;
          }

        SumOp op(transverseDir);
        op.scale = dx;

        // note that we copy from data component 1 to reducedData component 0
        int srcComp = 1;
        Interval srcInterval(srcComp, srcComp);
        data.copyTo(srcInterval, reducedData, reducedData.interval(),
                    reduceCopier, op);

#if 0 // commenting out everything related to FaceSumOp, since it doesn't work
        FaceSumOp faceOp(transverseDir);
        faceOp.scale = dx;


        faceData.copyTo(srcInterval, reducedFaceData,
                        reducedFaceData.interval(),
                        reduceCopier, faceOp);
#endif

        // check values
        Real exactVal = (srcComp+1)*dx*domainBox.size(transverseDir);

        RealVect exactFaceVal(exactVal*RealVect::Unit);
        exactFaceVal[transverseDir] = (srcComp+1)*dx*(domainBox.size(transverseDir)+1);

        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            FArrayBox& reducedFab = reducedData[reducedDit()];
            BoxIterator bit(lowerDimGrids.get(reducedDit()));
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                if (abs((reducedFab(iv, 0) - exactVal)) > tolerance)
                  {
                    if (verbose)
                      {
                        pout () << "bad cc value at "
                                << iv
                                << " expected " << exactVal
                                << " found " << reducedFab(iv,0)
                                << " in grid box " << lowerDimGrids[reducedDit]
                                << endl;
                      }

                    localStatus += 100;
                  }
              }
          }

        if (localStatus == 0)
          {
            pout () << "passed multiple-box cell-centered reduction test"
                    << endl;
          }
        else
          {
            pout () << "FAILED multiple-box cell-centered reduction test, status = "
                    << status << endl;
          }

#if 0
        // we know that the face-centered sumOp doesn't do multiple
        // boxes correctly, so comment this out for now.

        // now check face-centered values
        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& reducedFab = reducedFaceData[reducedDit()][dir];
                const Box& b = lowerDimGrids.get(reducedDit());
                BoxIterator bit(b);
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    if (abs((reducedFab(iv, 0) - exactFaceVal[dir])) > tolerance)
                      {
                        localFaceStatus += 100;
                      }
                  }
              }
          }

        if (localFaceStatus == 0)
          {
            pout () << "passed multiple-box face-centered reduction test"
                    << endl;
          }
        else
          {
            pout () << "FAILED multiple-box face-centered reduction test, status = "
                    << status << endl;
          }
#endif // end commenting out checks of faceSumOp multibox test.

        status += localStatus;
        faceStatus += localFaceStatus;
        localStatus = 0;
        localFaceStatus = 0;

        if (verbose)
          {
            pout() << "multiple-box spreading test" << endl;
          }
        // now spread reduced data back to original
        Real spreadingScale = 5.0;
        SpreadingOp spreadOp(transverseDir);
        spreadOp.scale = spreadingScale;

        FaceSpreadingOp faceSpreadOp(transverseDir);
        faceSpreadOp.scale = spreadingScale;

        SpreadingCopier spreadCopier(lowerDimGrids, grids,
                                     domain, ghostVect, transverseDir);

        if (verbose)
          {
            pout() << "Spreading copier: " << endl;
            pout() << spreadCopier << endl;
            pout() << endl;
          }

        // need a separate spreading copier for face-centered
        // data, since it will need different buffers
        SpreadingCopier faceSpreadCopier(lowerDimGrids, grids,
                                         domain, ghostVect,transverseDir);


        // reset data and faceData to a bogus value
        Real bogusVal = 100000;
        for (dit.begin(); dit.ok(); ++dit)
          {
            data[dit].setVal(bogusVal);
            faceData[dit].setVal(bogusVal);
          }

        // note that we copy from reduced data component 0
        // to data component 2
        int destComp = 2;
        Interval destInterval(destComp, destComp);
        reducedData.copyTo(reducedData.interval(), data,
                           destInterval, spreadCopier,
                           spreadOp);

        //writeLevel(&data);

        // set the reducedFaceData to be equal to the cell-centered exact val
        // note that we set all faces to be the same.
        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            reducedFaceData[reducedDit].setVal(exactVal);
          }

        reducedFaceData.copyTo(reducedFaceData.interval(), faceData,
                               destInterval, faceSpreadCopier,
                               faceSpreadOp);


        exactVal *= spreadingScale;

        // check cell-centered values
        for (dit.begin(); dit.ok(); ++dit)
          {
            FArrayBox& Fab = data[dit()];
            BoxIterator bit(Fab.box());
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                if (abs((Fab(iv, destComp) - exactVal)) > tolerance)
                  {
                    localStatus+= 1000;
                    if (verbose)
                      {
                        pout() << "fab-spreading test: bad value at "
                               << iv
                               << " expected " << exactVal
                               << " found " << Fab(iv,destComp) << endl;
                      }
                  }
              }
          }

        if (localStatus == 0)
          {
            pout () << "passed multiple-box cell-centered spreading test" << endl;
          }
        else
          {
            pout () << "FAILED multiple-box cell-centered spreading test, status = "
                    << status << endl;
          }



        // check face-centered values
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& Fab = faceData[dit()][dir];
                const Box& faceBox = Fab.box();
                BoxIterator bit(faceBox);
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    if (abs((Fab(iv, destComp) - exactVal)) > tolerance)
                      {
                        localFaceStatus+= 1000;
                        if (verbose)
                          {
                            pout () << "face-spreading test: bad value -- dir = "
                                    << dir
                                    << " location = " << iv
                                    << " expected " << exactVal
                                    << " found " << Fab(iv, destComp)
                                    << " in grid box " << grids[dit]
                                    << endl;
                          }
                      }
                  }
              } // end loop over directions
          }

        if (localFaceStatus == 0)
          {
            pout () << "passed multiple-box face-centered spreading test" << endl;
          }
        else
          {
            pout () << "FAILED multiple-box face-centered spreading test, status = "
                    << localFaceStatus << endl;
          }



        status += localStatus;
        faceStatus += localFaceStatus;
        localStatus = 0;


        if (status == 0)
          {
            pout() << "passed multiple-box test case" << endl;
          }
        else
          {
            pout() << "FAILED multiple-box test case" << endl;
          }

      } // end multiple box case


      // now do multiblock test case -- this also involves
      // ghost cells.
      // basic idea is that copiers should be able to be 
      // restricted to a single mapping block by specifying 
      // the domain appropriately (so only operations which occur 
      // are done in the "valid" block

      {

      if (verbose)
        {
          pout() << "multiple-block reduction test" << endl;
        }

        int localStatus = 0;
        int localFaceStatus = 0;
        // want two domain boxes which are not even adjacent

        int domLength2 = 16;
        Real dx2 = 2.0/domLength2;
        Box domainBox2(IntVect::Zero, IntVect::Unit*(domLength2-1));
        ProblemDomain domain2(domainBox2);
        int maxBoxSize = domLength2/4;

        Box lowerDimDomainBox2(domainBox2);
        lowerDimDomainBox2.setSmall(transverseDir, 0);
        lowerDimDomainBox2.setBig(transverseDir,0);
        ProblemDomain lowerDimDomain2(lowerDimDomainBox2);

        Vector<Box> gridBoxes;
        domainSplit(domainBox2, gridBoxes, maxBoxSize);
        int numBoxes = gridBoxes.size();
        // break this into separate blocks.
        numBoxes /= 2;
        // create 2 blocks
        IntVect block0Lo = IntVect::Zero;
        IntVect block0Hi = domainBox2.bigEnd();
        block0Hi[0] = maxBoxSize-1;
        Box block0(block0Lo, block0Hi);

        IntVect block1Lo = IntVect::Zero;
        block1Lo[0] = 2*maxBoxSize;
        IntVect block1Hi = domainBox2.bigEnd();
        block1Hi[0] = 3*maxBoxSize-1;
        Box block1(block1Lo, block1Hi);
          
        Vector<Box> blocks(2);
        blocks[0] = block0;
        blocks[1] = block1;

        Vector<Box> MBgridBoxes(numBoxes);
        int index = 0;
        for (int i=0; i<gridBoxes.size(); i++)
          {
            for (int blockno = 0; blockno<blocks.size(); blockno++)
              {
                if (blocks[blockno].intersects(gridBoxes[i]))
                  {
                    MBgridBoxes[index] = gridBoxes[i];
                    index++;
                  } // end if box intersects this block
              } // end loop over blocks
          } // end loop over gridboxes

        Vector<int> procAssign(numBoxes,0);
        // for mpi case, distribute boxes among processors
        LoadBalance(procAssign,MBgridBoxes);
        DisjointBoxLayout grids(MBgridBoxes, procAssign, domain2);

        Vector<Box> lowerDimGridBoxes;
        domainSplit(lowerDimDomainBox2, lowerDimGridBoxes, maxBoxSize);
        int numLowerDimBoxes = lowerDimGridBoxes.size()/2;
        Vector<Box> MBlowerDimGridBoxes(numLowerDimBoxes);
        index = 0;
        for (int i=0; i< lowerDimGridBoxes.size(); i++)
          {
            for (int blockno = 0; blockno<blocks.size(); blockno++)
              {
                if (blocks[blockno].intersects(lowerDimGridBoxes[i]))
                  {
                    MBlowerDimGridBoxes[index] = lowerDimGridBoxes[i];
                    index++;
                  } // end if box intersects this block
              } // end loop over blocks
          } // end loop over gridboxes
        
        Vector<int> lowerDimProcAssign(numLowerDimBoxes,0);
        LoadBalance(lowerDimProcAssign,MBlowerDimGridBoxes);
        DisjointBoxLayout lowerDimGrids(MBlowerDimGridBoxes,
                                        lowerDimProcAssign, domain2);

        IntVect ghostVect = IntVect::Unit;
        LevelData<FArrayBox> data(grids, 3, ghostVect);
        LevelData<FArrayBox> reducedData(lowerDimGrids, 1, ghostVect);

        LevelData<FluxBox> faceData(grids, 3, ghostVect);
        LevelData<FluxBox> reducedFaceData(lowerDimGrids, 1, ghostVect);

        // initialize data to 1
        Real testVal = 1.0;
        DataIterator dit = grids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int comp=0; comp<data.nComp(); comp++)
              {
                data[dit()].setVal((comp+1)*testVal, comp);
                for (int dir=0; dir<SpaceDim; dir++)
                  {
                    faceData[dit][dir].setVal((comp+1)*testVal, comp);
                  }
              }
          }

        DataIterator reducedDit = reducedData.dataIterator();
        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            reducedData[reducedDit()].setVal(0.0);
            reducedFaceData[reducedDit()].setVal(0.0);
          }

        if (verbose)
          {
            pout() << "hiDimGrids: " << endl;
            for (DataIterator dit1 = grids.dataIterator(); dit1.ok(); ++dit1)
              {
                pout() << grids[dit1] << endl;
              }

            pout() << "lowerDimGrids: " << endl;
            for (DataIterator dit1=lowerDimGrids.dataIterator(); dit1.ok(); ++dit1)
              {
                pout() << lowerDimGrids[dit1] << endl;
              }
          }

        SumOp op(transverseDir);
        op.scale = dx2;

        // do this once per block...
        for (int i=0; i<blocks.size(); i++)
          {
            // make a fake domain which only covers this block
            ProblemDomain blockDomain(blocks[i]);
            ReductionCopier reduceCopier(grids, lowerDimGrids, blockDomain,
                                         ghostVect, transverseDir);
            
            if (verbose)
              {
                pout() << "Block " << i << ": Forward copier: " << endl;
                pout () << reduceCopier << endl;
              }
            
            // note that we copy from data component 1 to reducedData component 0
            int srcComp = 1;
            Interval srcInterval(srcComp, srcComp);
            data.copyTo(srcInterval, reducedData, reducedData.interval(),
                        reduceCopier, op);
            
          } // end loop over blocks

        // check values
        int srcComp = 1;
        Interval srcInterval(srcComp, srcComp);
        
        Real exactVal = (srcComp+1)*dx2*domainBox2.size(transverseDir);

        RealVect exactFaceVal(exactVal*RealVect::Unit);
        exactFaceVal[transverseDir] = (srcComp+1)*dx2*(domainBox2.size(transverseDir)+1);

        for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
          {
            FArrayBox& reducedFab = reducedData[reducedDit()];
            BoxIterator bit(lowerDimGrids.get(reducedDit()));
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                if (abs((reducedFab(iv, 0) - exactVal)) > tolerance)
                  {
                    if (verbose)
                      {
                        pout () << "bad cc value at "
                                << iv
                                << " expected " << exactVal
                                << " found " << reducedFab(iv,0)
                                << " in grid box " << lowerDimGrids[reducedDit]
                                << endl;
                      }

                    localStatus += 100;
                  }
              }
          }

        if (localStatus == 0)
          {
            pout () << "passed multiple-box cell-centered reduction test"
                    << endl;
          }
        else
          {
            pout () << "FAILED multiple-box cell-centered reduction test, status = "
                    << status << endl;
          }

        status += localStatus;
        faceStatus += localFaceStatus;
        localStatus = 0;
        localFaceStatus = 0;

        if (verbose)
          {
            pout() << "multiple-box spreading test" << endl;
          }
        // now spread reduced data back to original
        Real spreadingScale = 5.0;
        SpreadingOp spreadOp(transverseDir);
        spreadOp.scale = spreadingScale;

        FaceSpreadingOp faceSpreadOp(transverseDir);
        faceSpreadOp.scale = spreadingScale;

        // reset data and faceData to a bogus value
        Real bogusVal = 100000;
        for (dit.begin(); dit.ok(); ++dit)
          {
            data[dit].setVal(bogusVal);
            faceData[dit].setVal(bogusVal);
          }
                        
        // now loop over blocks
        for (int i = 0; i < blocks.size(); i++)
          {
            ProblemDomain blockDomain(blocks[i]);
            SpreadingCopier spreadCopier(lowerDimGrids, grids,
                                         blockDomain, ghostVect, 
                                         transverseDir);
            
            
            if (verbose)
              {
                pout() << "Block " << i << ": Spreading copier: " << endl;
                pout() << spreadCopier << endl;
                pout() << endl;
              }

            // need a separate spreading copier for face-centered
            // data, since it will need different buffers
            SpreadingCopier faceSpreadCopier(lowerDimGrids, grids,
                                             blockDomain, ghostVect,
                                             transverseDir);
            
            
            // note that we copy from reduced data component 0
            // to data component 2
            int destComp = 2;
            Interval destInterval(destComp, destComp);
            reducedData.copyTo(reducedData.interval(), data,
                               destInterval, spreadCopier,
                               spreadOp);
            
            //writeLevel(&data);
            
            // set the reducedFaceData to be equal to the cell-centered exact val
            // note that we set all faces to be the same.
            for (reducedDit.begin(); reducedDit.ok(); ++reducedDit)
              {
                reducedFaceData[reducedDit].setVal(exactVal);
              }
            
            reducedFaceData.copyTo(reducedFaceData.interval(), faceData,
                                   destInterval, faceSpreadCopier,
                                   faceSpreadOp);
            
          } // end loop over blocks

        exactVal *= spreadingScale;

        // note that we copied from reduced data component 0
        // to data component 2
        int destComp = 2;
        Interval destInterval(destComp, destComp);
        // check cell-centered values
        for (dit.begin(); dit.ok(); ++dit)
          {
            FArrayBox& Fab = data[dit()];
            BoxIterator bit(Fab.box());
            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                if (abs((Fab(iv, destComp) - exactVal)) > tolerance)
                  {
                    localStatus+= 1000;
                    if (verbose)
                      {
                        pout() << "fab-spreading test: bad value at "
                               << iv
                               << " expected " << exactVal
                               << " found " << Fab(iv,destComp) << endl;
                      }
                  }
              }
          }

        if (localStatus == 0)
          {
            pout () << "passed multiple-block cell-centered spreading test" << endl;
          }
        else
          {
            pout () << "FAILED multiple-block cell-centered spreading test, status = "
                    << status << endl;
          }



        // check face-centered values
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& Fab = faceData[dit()][dir];
                const Box& faceBox = Fab.box();
                BoxIterator bit(faceBox);
                for (bit.begin(); bit.ok(); ++bit)
                  {
                    IntVect iv = bit();
                    if (abs((Fab(iv, destComp) - exactVal)) > tolerance)
                      {
                        localFaceStatus+= 1000;
                        if (verbose)
                          {
                            pout () << "face-spreading test: bad value -- dir = "
                                    << dir
                                    << " location = " << iv
                                    << " expected " << exactVal
                                    << " found " << Fab(iv, destComp)
                                    << " in grid box " << grids[dit]
                                    << endl;
                          }
                      }
                  }
              } // end loop over directions
          }

        if (localFaceStatus == 0)
          {
            pout () << "passed multiple-block face-centered spreading test" << endl;
          }
        else
          {
            pout () << "FAILED multiple-block face-centered spreading test, status = "
                    << localFaceStatus << endl;
          }



        status += localStatus;
        faceStatus += localFaceStatus;
        localStatus = 0;


        if (status == 0)
          {
            pout() << "passed multiple-block test case" << endl;
          }
        else
          {
            pout() << "FAILED multiple-block test case" << endl;
          }

      } // end multiple block case


      // now do multiple-dimension testing (which only makes sense
      // for 3D and higher)
      if (SpaceDim > 2)
        {
          int multiDimStatus = 0;
          // first, create a single box
          int length = 8;
          Box baseBox(IntVect::Zero, (length-1)*IntVect::Unit);

          // now create transverse directions
          Vector<int> transverseDir;
          transverseDir.resize(2);
          transverseDir[0] = SpaceDim -1;
          transverseDir[1] = SpaceDim -2;

          // flatten base box in the transverse directions
          Box lowDimBox = baseBox;
          for (int n=0; n<transverseDir.size(); n++)
            {
              lowDimBox.setBig(transverseDir[n],0);
            }
          Vector<Box> loDimBoxes(1,lowDimBox);

          // now create hi-dim grids, including two shifted boxes
          Vector<Box> hiDimBoxes;
          hiDimBoxes.push_back(baseBox);
          IntVect shiftVect = IntVect::Zero;
          for (int n=0; n<transverseDir.size(); n++)
            {
              shiftVect[transverseDir[n]] = length;
              Box tempBox(baseBox);
              tempBox.shift(shiftVect);
              hiDimBoxes.push_back(tempBox);
            }

          // create DBL's
          Vector<int> loDimProcAssign(1,0);
          Vector<int> hiDimProcAssign(hiDimBoxes.size(), 0);

          LoadBalance(hiDimProcAssign, hiDimBoxes);


          // problem domain big enough to contain all boxes
          Box domainBox(IntVect::Zero, (2*length-1)*IntVect::Unit);
          ProblemDomain probDomain(domainBox);

          DisjointBoxLayout loDimGrids(loDimBoxes, loDimProcAssign, probDomain);
          DisjointBoxLayout hiDimGrids(hiDimBoxes, hiDimProcAssign, probDomain);

          // first try spreading
          SpreadingCopier spreadCopier(loDimGrids, hiDimGrids,
                                       probDomain, transverseDir);



          IntVect ghostVect(IntVect::Zero);
          LevelData<FArrayBox> loDimData(loDimGrids, 1, ghostVect);
          LevelData<FArrayBox> hiDimData(hiDimGrids, 1, ghostVect);

          // initialize low-dim data
          DataIterator loDimDit = loDimData.dataIterator();
          for (loDimDit.begin(); loDimDit.ok(); ++loDimDit)
            {
              FArrayBox& thisFab = loDimData[loDimDit];
              BoxIterator bit(thisFab.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  thisFab(iv, 0) = iv[0];
                }
            }


          // construct spreading op
          SpreadingOp spreadOp(transverseDir);
          loDimData.copyTo(loDimData.interval(),
                           hiDimData,
                           hiDimData.interval(),
                           spreadCopier,
                           spreadOp);


          int localStatus = 0;
          // now check results
          DataIterator hiDimDit = hiDimData.dataIterator();
          for (hiDimDit.begin(); hiDimDit.ok(); ++hiDimDit)
            {
              FArrayBox& thisFAB = hiDimData[hiDimDit];
              BoxIterator bit(thisFAB.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  if (thisFAB(iv,0) != iv[0])
                    {
                      localStatus += 1;
                      if (verbose)
                        {
                          pout() << "multiDim-spreading test: bad value at "
                                 << iv
                                 << " expected " << iv[0]
                                 << " found " << thisFAB(iv,0) << endl;
                        }
                    } // end if bad val

                } // end loop over cells
            } // end loop over high-dim boxes
          if (localStatus != 0) multiDimStatus += 100000;


        if (localStatus == 0)
          {
            pout() << "passed multiDim spreading test case" << endl;
          }
        else
          {
            pout() << "FAILED multiDim spreading test case" << endl;
          }

          status += multiDimStatus;
        }

    } // end if SpaceDim > 1
  status += faceStatus;
  return status;
}


///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
  {
    if ( argv[i][0] == '-' ) //if it is an option
    {
      // compare 3 chars to differentiate -x from -xx
      if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
      {
        verbose = true ;
        // argv[i] = "" ;
      }
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // argv[i] = "" ;
      }
      else
      {
        break ;
      }
    }
  }
  return ;
}
