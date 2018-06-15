#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Misc.H"   //Abs()
#include "LevelFluxRegister.H"
#include "LoadBalance.H"
#include "DebugDump.H"
#include "LayoutIterator.H"
#include "AMRIO.H"
#include "UsingNamespace.H"

/*****************/
/*****************/
void
setDefaults( int& nref,
             Vector<Box>& coarboxes,
             Vector<Box>& fineboxes,
             Box& domc, Box& domf);

int fluxRegTest();
void parseTestOptions(int argc ,char* argv[]) ;
/*****************/
/*****************/

/// Global variables for handling output:
static const char *pgmname = "newFluxRegisterTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = false;

/*****************/
/*****************/

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int eekflag = fluxRegTest();
  if (eekflag == 0)
    {
      pout() << indent << pgmname
           << ": test passed." << endl;
    }
  else
    {
      pout() << indent << pgmname
           << ": test FAILED with error code "
           << eekflag << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return eekflag;
}

/*****************/
/*****************/

int fluxRegTest()
{
#ifdef CH_USE_HDF5
  writeLevel(NULL);
#endif
  int retflag = 0;

  int nref;
  Vector<Box> fineboxes;
  Vector<Box> coarboxes;
  Box domf, domc;

  //set coarse and fine grid boxes
  setDefaults(nref, coarboxes,  fineboxes, domc, domf);

  {
    // first do nonperiodic test

    Interval interv(0,0);

    //set up coarse and fine grids
    DisjointBoxLayout dblFineCell,dblCoarCell;
    Vector<int> procAssignCoar(coarboxes.size(), 0);
    Vector<int> procAssignFine(fineboxes.size(), 0);
    LoadBalance(procAssignCoar, coarboxes);
    LoadBalance(procAssignFine, fineboxes);

    dblCoarCell.define(coarboxes, procAssignCoar);
    dblFineCell.define(fineboxes, procAssignFine);

    dblCoarCell.close();
    dblFineCell.close();

    LevelData<FArrayBox> coarData(dblCoarCell, 1);
    LevelData<FArrayBox> fineData(dblFineCell, 1);

    DataIterator coarIt = coarData.dataIterator();
    DataIterator fineIt = fineData.dataIterator();

    LevelFluxRegister fluxReg(dblFineCell,
                              dblCoarCell,
                              domf, nref,
                              1);

    //set data and flux registers to zero
    for (coarIt.reset(); coarIt.ok(); ++coarIt)
      coarData[coarIt()].setVal(0.);

    fluxReg.setToZero();

    //increment and decrement
    //flux registers with equal size fluxes
    Real scale = 1.0;
    Real fluxVal = 4.77;
    for (coarIt.reset(); coarIt.ok(); ++coarIt)
      {
        const Box&  cellBoxCoar = dblCoarCell.get(coarIt());
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            Box edgeBoxCoar = surroundingNodes(cellBoxCoar, idir);
            FArrayBox edgeFlux(edgeBoxCoar,1);
            edgeFlux.setVal(fluxVal);
            DataIndex dataIndGlo = coarIt();
            fluxReg.incrementCoarse(edgeFlux, scale,
                                    dataIndGlo, interv, interv, idir);
          }
      }

    for (fineIt.reset(); fineIt.ok(); ++fineIt)
      {
        const Box&  cellBoxFine = dblFineCell.get(fineIt());
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            Box edgeBoxFine = surroundingNodes(cellBoxFine, idir);
            FArrayBox edgeFlux(edgeBoxFine,1);
            edgeFlux.setVal(fluxVal);
            SideIterator sit;
            DataIndex dataIndGlo = fineIt();
            for (sit.reset(); sit.ok(); ++sit)
              {
                fluxReg.incrementFine(edgeFlux, scale,
                                      dataIndGlo,  interv, interv, idir, sit());
              }
          }
      }

    //reflux what ought to be zero into zero and the result should be zero
    fluxReg.reflux(coarData, scale);

    DataIterator datIt = coarData.dataIterator();
    for (datIt.reset(); datIt.ok(); ++datIt)
      {
        const FArrayBox& data = coarData[datIt()];
        Real rmax = Abs(data.max());
        Real rmin = Abs(data.min());
        if ((rmax > 1.0e-10)||(rmin > 1.0e-10))
          {
            pout() << indent << pgmname
                 << ": fluxRegister failed the nonperiodic conservation test = " << endl;
            retflag = 1;
          }
      }
  }
  // end non-periodic test

  // now do the same thing all over again, this time with a periodic domain
  {
    ProblemDomain coarseDomain(domc);
    ProblemDomain fineDomain(domf);
    for (int dir=0; dir<SpaceDim; dir++)
      {
        coarseDomain.setPeriodic(dir, true);
        fineDomain.setPeriodic(dir, true);
      }

    Interval interv(0,0);

    //set up coarse and fine grids
    DisjointBoxLayout dblFineCell,dblCoarCell;
    Vector<int> procAssignCoar(coarboxes.size(), 0);
    Vector<int> procAssignFine(fineboxes.size(), 0);
    LoadBalance(procAssignCoar, coarboxes);
    LoadBalance(procAssignFine, fineboxes);

    dblCoarCell.define(coarboxes, procAssignCoar, coarseDomain);
    dblFineCell.define(fineboxes, procAssignFine, fineDomain);

    dblCoarCell.close();
    dblFineCell.close();

    LevelData<FArrayBox> coarData(dblCoarCell, 1);
    LevelData<FArrayBox> fineData(dblFineCell, 1);

    DataIterator coarIt = coarData.dataIterator();
    DataIterator fineIt = fineData.dataIterator();

    LevelFluxRegister fluxReg(dblFineCell,
                              dblCoarCell,
                              fineDomain, nref,
                              1);

    //set data and flux registers to zero
    for (coarIt.reset(); coarIt.ok(); ++coarIt)
      coarData[coarIt()].setVal(0.);

    fluxReg.setToZero();

    //increment and decrement
    //flux registers with equal size fluxes
    Real scale = 1.0;
    Real fluxVal = 4.77;
    for (coarIt.reset(); coarIt.ok(); ++coarIt)
      {
        const Box&  cellBoxCoar = dblCoarCell.get(coarIt());
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            Box edgeBoxCoar = surroundingNodes(cellBoxCoar, idir);
            FArrayBox edgeFlux(edgeBoxCoar,1);
            edgeFlux.setVal(fluxVal);
            DataIndex dataIndGlo = coarIt();
            fluxReg.incrementCoarse(edgeFlux, scale,
                                    dataIndGlo, interv, interv, idir);
          }
      }

    for (fineIt.reset(); fineIt.ok(); ++fineIt)
      {
        const Box&  cellBoxFine = dblFineCell.get(fineIt());
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            Box edgeBoxFine = surroundingNodes(cellBoxFine, idir);
            FArrayBox edgeFlux(edgeBoxFine,1);
            edgeFlux.setVal(fluxVal);
            SideIterator sit;
            DataIndex dataIndGlo = fineIt();
            for (sit.reset(); sit.ok(); ++sit)
              {
                fluxReg.incrementFine(edgeFlux, scale,
                                      dataIndGlo,  interv, interv, idir, sit());
              }
          }
      }

    //reflux what ought to be zero into zero and the result should be zero
    fluxReg.reflux(coarData, scale);

    DataIterator datIt = coarData.dataIterator();
    for (datIt.reset(); datIt.ok(); ++datIt)
      {
        const FArrayBox& data = coarData[datIt()];
        Real rmax = Abs(data.max());
        Real rmin = Abs(data.min());
        if ((rmax > 1.0e-10)||(rmin > 1.0e-10))
          {
            pout() << indent << pgmname
                 << ": fluxRegister failed the periodic conservation test " << endl;
            retflag += 2;
          }
      }
  } // end periodic test

  return retflag;
}

/*****************/
void
setDefaults( int& nref,
             Vector<Box>& coarboxes,
             Vector<Box>& fineboxes,
             Box& domc, Box& domf)
{
  nref = 2;
  int nc = 16;
  int nc16 = nc/2;
  int nc8 = nc/4;
  int nc24 = 3*nc/4;
  int nc12 = 3*nc/8;
  int nc20 = 5*nc/8;

  IntVect ivclo = IntVect::Zero;

  fineboxes.resize(0);

  IntVect ivchi = (nc-1) * IntVect::Unit;

#if (CH_SPACEDIM ==1)
  fineboxes.resize(3);
  fineboxes[0].define(IntVect(0), IntVect(nc8-1));
  fineboxes[1].define(IntVect(nc12), IntVect(nc20-1));
  fineboxes[2].define(IntVect(nc24), ivchi);

  // this is to prevent unused variable warning in 1d
  int temp =nc16;
  nc16 = temp;

#else
  Box boxf1(IntVect(D_DECL6(0 ,nc8,0, 0,0,0)),
            IntVect(D_DECL6(nc16-1,nc24,nc-1, nc-1,nc-1,nc-1)));
  Box boxf2(IntVect(D_DECL6(nc16,nc12,0, 0,0,0)),
            IntVect(D_DECL6(nc24-1,nc20-1,nc-1, nc-1,nc-1,nc-1)));
  Box boxf3(IntVect(D_DECL6(nc8 ,0 ,0, 0,0,0)),
            IntVect(D_DECL6(nc20-1,nc8-1,nc-1, nc-1,nc-1,nc-1)));
  Box boxf4(IntVect(D_DECL6(nc20,nc8,0, 0,0,0)),
            IntVect(D_DECL6(nc-1,nc12-1,nc-1, nc-1,nc-1,nc-1)));

  fineboxes.resize(4);
  fineboxes[0] = boxf1;
  fineboxes[1] = boxf2;
  fineboxes[2] = boxf3;
  fineboxes[3] = boxf4;
#endif

  for (int i=0; i<fineboxes.size(); ++i)
    {
      fineboxes[i].refine(nref);
    }

  domc = Box(ivclo, ivchi);
  int ichop = domc.size(0)/2;
  Box bc1 = domc;
  Box bc2 = bc1.chop(0, ichop);
#if (CH_SPACEDIM > 1)
  Box bc3 = bc2.chop(1, ichop);
  coarboxes.push_back(bc3);
#endif
  coarboxes.push_back(bc1);
  coarboxes.push_back(bc2);


  domf = refine(domc,nref);
}

///
// Parse the standard test options (-v -q) out of the command line.
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
        }
    }
  return ;
}
