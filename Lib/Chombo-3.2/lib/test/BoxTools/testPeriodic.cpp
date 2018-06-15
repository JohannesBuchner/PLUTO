#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Purpose:
//  Test the infrastucture for periodic boundary conditions
//
// Usage:
//  <program-name> [-q|-v] ...
//
//  where:
//    -q means run quietly (only pass/fail messages printed)
//    -v means run verbosely (all messages printed)
//    -write means write out the boxes; otherwise will read in
//           boxes from a file and compare with boxes generated
//           in tests.
//    ... all non-option arguments are ignored (Chombo convention)
//
//  Default is `-v'
//
//  Unknown options are treated as errors.
//

// Include files:
#include <cstdio>

#include "parstream.H"
#include "LoadBalance.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H"


#ifdef CH_MPI
#include <mpi.h>
#endif

#include "UsingNamespace.H"

//////////////////////////////////////////////////////////////
using std::endl;

void parseTestOptions(int argc, char* argv[]);

void initData(LevelData<FArrayBox>& a_data, const Real a_dx);

Real dataVal(const IntVect& a_loc, int a_comp);

//
// For testing LevelData::apply()...a function, and a functor:
//
void leveldataApplyFunc(const Box& box, int comp, FArrayBox& fab)
{
  fab.setVal( 12.34, box, 0 );
}
class LevelDataApplyFunctor : public LevelData<FArrayBox>::ApplyFunctor
{
public:
  LevelDataApplyFunctor( Real x ) : m_x(x)
  {
  }

  virtual void operator()( const Box& box, int comp, FArrayBox& fab ) const
  {
    fab.setVal( m_x, box, comp );
  }

private:
  Real m_x;
};

/// Global variables for handling output

static const char *pgmname = "testPeriodic";
static const char *indent = "   ";
//static const char *indent2 = "      " ;

static bool verbose = false ;

/// Code:

int
main(int argc ,char *argv[] )
{

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions(argc, argv);

  // establish periodic domain -- first multiply periodic in all directions
  int status = 0;
  int baseDomainSize = 8;
  int numGhost = 2;
  IntVect ghostVect(numGhost*IntVect::Unit);

  Box baseDomBox(IntVect::Zero, (baseDomainSize-1)*IntVect::Unit);
  ProblemDomain baseDomain(baseDomBox);

  // set periodic in all directions
  for (int dir=0; dir<SpaceDim; dir++)
    {
      baseDomain.setPeriodic(dir, true);
    }

  // quick check
  //  should be (3^SpaceDim) -1 vectors
  int numVect =0;
  if (verbose) pout() << "ShiftIterator shift vectors:" << endl;
  ShiftIterator sit1 = baseDomain.shiftIterator();
  for (sit1.begin(); sit1.ok(); ++sit1)
    {
      numVect++;
      if (verbose) pout() << "   " << sit1() << endl;
    }
  Real exact =  pow(3.0,SpaceDim) -1;
  Real eps = 0.0001;
  if ((exact - numVect) > eps)
    {
      //fail
      if (verbose)
        {
          pout() << "failed ShiftIterator count test " << endl;
        }
      status += 1;
    }

  ShiftIterator sit2 = baseDomain.shiftIterator();

  {
    // first try a single box on this level
    if (verbose) pout() <<  "Single box test -- " << endl;

    Vector<Box> levelBoxes(1, baseDomBox);
    Vector<int> procAssign(1, 0);
    DisjointBoxLayout levelGrids(levelBoxes,procAssign, baseDomain);

    LevelData<FArrayBox> tester(levelGrids, 1, ghostVect);

    Real dx = 1.0/baseDomainSize;

    initData(tester, dx);

    // test exchange in this case -- should fill ghost cells with
    // periodically copied data
    if (verbose) pout() << "    begin exchange..." << endl;

    Interval comps = tester.interval();
    tester.exchange(comps);
    if (verbose) pout() << "    done exchange..." << endl;
    DataIterator dit = tester.dataIterator();

    for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& testFab = tester[dit()];
        // now check to be sure that this passed
        for (int dir=0; dir<SpaceDim; dir++)
          {
            if (verbose) pout() << "check lo cells in dir " << dir << endl;
            Box loGhostBox = adjCellLo(baseDomBox, dir, numGhost);
            BoxIterator loBit(loGhostBox);
            for (loBit.begin(); loBit.ok(); ++loBit)
              {
                const IntVect iv = loBit();
                IntVect validLoc = iv;
                validLoc.shift(dir, baseDomBox.size(dir));
                for (int comp=0; comp<tester.nComp(); comp++)
                  {
                    Real validVal = dataVal(validLoc, comp);
                    if (testFab(iv, comp) != validVal)
                      {
                        if (verbose)
                          {
                            pout() << "failed single-box exchange test at "
                                   << iv << endl;
                          }
                        status += 10;
                      }
                  }
              } // end loop over lo ghost cells

            if (verbose) pout() << "check hi cells in dir " << dir << endl;
            Box hiGhostBox = adjCellHi(baseDomBox, dir, numGhost);
            BoxIterator hiBit(hiGhostBox);
            for (hiBit.begin(); hiBit.ok(); ++hiBit)
              {
                const IntVect iv = hiBit();
                IntVect validLoc = iv;
                validLoc.shift(dir, -baseDomBox.size(dir));
                for (int comp=0; comp<tester.nComp(); comp++)
                  {
                    Real validVal = dataVal(validLoc, comp);
                    if (testFab(iv, comp) != validVal)
                      {
                        if (verbose)
                          {
                            pout() << "failed single-box exchange test at "
                                   << iv << endl;
                          }
                        status += 100;
                      }
                  }
              } // end loop over hi ghost cells
          } // end loop over directions
      } // end loop over boxes
    if (verbose) pout() << "done single-box test" << endl;
  } // end single-box test

  // now try refining things, this time look at 2 boxes
  baseDomain.refine(2);
  {
    if (verbose) pout() << "2-box test..." << endl;

    const Box domainBox = baseDomain.domainBox();
    Vector<Box> levelBoxes1(1, domainBox);
    Vector<int> procAssign1(1,0);
    DisjointBoxLayout levelGrids1(levelBoxes1, procAssign1, baseDomain);

    LevelData<FArrayBox> tester0(levelGrids1, 1, IntVect::Zero);
    LevelData<FArrayBox> tester1(levelGrids1, 1, ghostVect);

    Real dx = 1.0/(2.0*baseDomainSize);
    initData(tester1, dx);
    initData(tester0, dx);

    // first test exchange for single box
    Interval comps = tester1.interval();
    Copier ex;
    ex.exchangeDefine(levelGrids1, ghostVect);
    pout() << "copier = " << ex << endl;
    tester1.exchange(comps, ex);

    // now try multiple boxes
    Vector<Box> levelBoxes2(3);
    Vector<int> procAssign2(3);

    if (SpaceDim == 1)
      {
        // do this a little differently in 1d
        levelBoxes2[0] = Box(IntVect::Zero, 9*IntVect::Unit);
        levelBoxes2[1] = Box(10*IntVect::Unit, 11*IntVect::Unit);
        levelBoxes2[2] = Box(14*IntVect::Unit, 15*IntVect::Unit);
      }
    else
      {
        levelBoxes2[0] = Box(IntVect::Zero, 9*IntVect::Unit);
        levelBoxes2[1] = Box(IntVect(D_DECL6(7,10,10,10,10,10)),
                             15*IntVect::Unit);
        // this box should fail the disjointness test
        //levelBoxes2[1] = Box(IntVect(D_DECL6(7,10,10,10,10,10)),17*IntVect::Unit);
        levelBoxes2[2] = Box(IntVect(D_DECL6(11,0,0,0,0,0)),
                             IntVect(D_DECL6(15,7,7,7,7,7)));
      }

    int loadbalancestatus = LoadBalance(procAssign2, levelBoxes2);
    CH_assert (loadbalancestatus == 0);

    DisjointBoxLayout levelGrids2;
    levelGrids2.define(levelBoxes2, procAssign2, baseDomain);
    LevelData<FArrayBox> tester2(levelGrids2, 1, ghostVect);

    if (verbose) pout() << "    begin copyTo test" << endl;
    tester0.copyTo(comps, tester2, comps);
    if (verbose) pout() << "    done copyTo, now test values" << endl;
    // check that periodic copies worked OK:

    DataIterator dit2 = tester2.dataIterator();
    for (dit2.begin(); dit2.ok(); ++dit2)
      {
        const Box& thisBox = tester2[dit2()].box();
        if (!domainBox.contains(thisBox))
          {
            FArrayBox& destFab = tester2[dit2];
            for (int dir=0; dir<SpaceDim; dir++)
              {
                if (verbose)
                  {
                    pout() << "      check lo cells in dir "
                           << dir << endl;
                  }

                Box loGhost = adjCellLo(domainBox, dir, numGhost);
                loGhost &= thisBox;
                if (!loGhost.isEmpty())
                {
                  BoxIterator loBit(loGhost);
                  for (loBit.begin(); loBit.ok(); ++loBit)
                    {
                      IntVect iv = loBit();
                      IntVect srcIV = iv;
                      srcIV.shift(dir, domainBox.size(dir));
                      for (int comp=0; comp<destFab.nComp(); comp++)
                        {
                          Real validVal = dataVal(srcIV, comp);
                          if (validVal != destFab(iv, comp))
                            {
                              if (verbose)
                                {
                                  pout() << "failed multi-box copy test at "
                                         << iv << endl;
                                }
                              status += 1000;
                            } // end if we fail test
                        } // end loop over components
                    } // end loop over ghost cells outside domain
                } // end if there are lo-end ghost cells

                if (verbose)
                  {
                    pout() << "      check hi cells in dir "
                           << dir << endl;
                  }

                Box hiGhost = adjCellHi(domainBox, dir, numGhost);
                hiGhost &= thisBox;
                if (!hiGhost.isEmpty())
                  {
                    BoxIterator hiBit(hiGhost);
                    for (hiBit.begin(); hiBit.ok(); ++hiBit)
                      {
                        IntVect iv = hiBit();
                        IntVect srcIV = iv;
                        srcIV.shift(dir, -domainBox.size(dir));
                        for (int comp=0; comp<destFab.nComp(); comp++)
                          {
                            Real validVal = dataVal(srcIV, comp);
                            if (validVal != destFab(iv, comp))
                              {
                                if (verbose)
                                  {
                                    pout() << "failed multi-box copy test at "
                                           << iv << endl;
                                  }
                                status += 10000;
                              } // end if we fail test
                          } // end loop over components
                      } // end loop over high-end ghost cells
                  } // end if there are hi-end ghost cells

              } // end loop over directions
          } // end if there are ghost cells
      } // end loop over grids in destination
    if (verbose) pout() << "done 2-box test" << endl;
  }

  // this last test makes no sense in 1D
  if (SpaceDim > 1)
    {
      if (verbose) pout() << "partially-periodic case" << endl;
      // now test partially periodic case (also coarsen back to original domain)
      baseDomain.coarsen(2);
      baseDomain.setPeriodic(0,false);

      // try similar tests to the second case, only coarsened
      const Box domainBox = baseDomain.domainBox();
      Vector<Box> levelBoxes1(1, domainBox);
      Vector<int> procAssign1(1,0);
      DisjointBoxLayout levelGrids1(levelBoxes1, procAssign1, baseDomain);

      LevelData<FArrayBox> tester0(levelGrids1, 1, IntVect::Zero);
      LevelData<FArrayBox> tester1(levelGrids1, 1, ghostVect);

      Real dx = 1.0/(baseDomainSize);
      initData(tester1, dx);
      initData(tester0, dx);

      if (verbose) pout() << "    partially-periodic exchange..." << endl;
      // first test exchange for single box
      Interval comps = tester1.interval();
      tester1.exchange(comps);

      if (verbose) pout() << "    ....  done" << endl;

      if (verbose) pout() << "   multi-box partially periodic case" << endl;
      // now try multiple boxes
      Vector<Box> levelBoxes2(3);
      Vector<int> procAssign2(3);

      levelBoxes2[0] = Box(IntVect::Zero, 4*IntVect::Unit);
      levelBoxes2[1] = Box(IntVect(D_DECL6(3,5,5,5,5,5)), 7*IntVect::Unit);
      // this box should fail the disjointness test
      //levelBoxes2[1] = Box(IntVect(D_DECL6(3,5,5)), 10*IntVect::Unit);
      levelBoxes2[2] = Box(IntVect(D_DECL6(6,0,0,0,0,0)),
                           IntVect(D_DECL6(7,1,1,1,1,1)));

      int loadbalancestatus = LoadBalance(procAssign2, levelBoxes2);
      CH_assert (loadbalancestatus == 0);

      DisjointBoxLayout levelGrids2(levelBoxes2, procAssign2, baseDomain);
      LevelData<FArrayBox> tester2(levelGrids2, 1, ghostVect);

      tester0.copyTo(comps, tester2, comps);

      if (verbose) pout() << "done partially-periodic case" << endl;

      // We don't care about the results of these -- just want to see if they compile:
      if ( argc == 12345 )
      {
        tester0.apply( leveldataApplyFunc );
        tester0.apply( LevelDataApplyFunctor(3.14159) );
      }
    }

  // I'm thinking all pout()'s should be before MPI_Finalize... (ndk)
  pout() << indent << pgmname << ": "
         << ( (status == 0) ? "passed all tests" : "failed at least one test,")
         << endl;


#ifdef CH_MPI
  MPI_Finalize();
#endif

  return status ;
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
          else
            {
              break ;
            }
        }
    }
  return ;
}

void
initData(LevelData<FArrayBox>& a_data, const Real a_dx)
{

  DataIterator dit = a_data.dataIterator();
  const DisjointBoxLayout& interiorBoxes = a_data.getBoxes();

  for (dit.begin(); dit.ok(); ++dit)
    {
      // first set to a bogus value which will persist in ghost cells
      // after initialization
      a_data[dit()].setVal(1.0e9);

      // this will be slow, but who cares?
      FArrayBox& localData = a_data[dit()];
      BoxIterator boxIt(interiorBoxes[dit()]);
      Real localVal;
      for (boxIt.begin(); boxIt.ok(); ++boxIt)
        {
          const IntVect& loc = boxIt();
          for (int comp=0; comp<localData.nComp(); comp++)
            {
              localVal = dataVal(loc, comp);
              localData(loc, comp) = localVal;
            }
        }

    }
}

Real dataVal(const IntVect& a_loc, int a_comp)
{
  Real localVal = a_comp*(D_TERM6(0.1*a_loc[0]+0.5,
                                  +a_loc[1]+0.5,
                                  +100*a_loc[2]+0.5,
                                  +10*a_loc[3]+0.5,
                                  +1000*a_loc[4]+0.5,
                                  +10000*a_loc[5]+0.5));

  return localVal;
}
