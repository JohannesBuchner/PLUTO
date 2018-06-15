#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//#if defined(NDEBUG)
//int  main()
//{
//  return -1;
//}
//#else

#include <cmath>
#include <algorithm>

#include "CH_HDF5.H"
#include "BRMeshRefine.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LoadBalance.H"
#include "DebugDump.H"
#include "UGIO.H"
#include "TestCommon.H"
#include "UsingNamespace.H"


void
parseTestOptions( int argc ,char* argv[] ) ;

int
test();

/// Global variables for handling output:
static const char *pgmname = "HDF5data" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  //  InitWriteFAB();
  parseTestOptions( argc ,argv ) ;

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int icode = test();
  if (icode != 0)
    {
      pout() << indent << pgmname <<" failed"<<endl;
    }
  else
    {
      pout() << indent << pgmname <<" passed"<<endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}

struct values
{
  static void setVal1(const Box& box, int comps, BaseFab<int>& fab)
    {
      int center = (box.smallEnd()[0] + box.bigEnd()[0])/2;
      fab.setVal(center, 0);
      fab.setVal(2*center, 1);
      fab.setVal(3*center, 2);

    }
  static void setVal2(const Box& box, int comps, BaseFab<int>& fab)
    {
      if ( SpaceDim == 1 )
        {
          fab.setVal(173);
        }
      else
        {
          int center = (box.smallEnd()[1] + box.bigEnd()[1])/2;
          fab.setVal(center);
        }
    }

  static void setVal3(const Box& box, int comps, FArrayBox& fab)
    {
      if ( SpaceDim < 3 )
        {
          fab.setVal(311);
        }
      else
        {
          int center = (box.smallEnd()[2] + box.bigEnd()[2])/2;
          for (int c=0; c<comps; ++c)
                    {
                      for (BoxIterator bit(box); bit.ok(); ++bit)
                            {
                              fab(bit(), c) = c + center + bit()[0];
                            }
                    }
        }
    }
};

struct
{
  IntVect offset;
  Real    value;
} typedef Moment;


// returns 0 on all tests passed.

int test()
{
#ifdef CH_USE_HDF5

  int error;
  HDF5Handle testFile;

  CH_assert(!testFile.isOpen());

  error = testFile.open("data.h5", HDF5Handle::CREATE);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "File creation failed "<<error<<endl;
      return error;
    }

  CH_assert(testFile.isOpen());

  Box domain(IntVect::Zero, 20*IntVect::Unit);

  DisjointBoxLayout plan1, plan2;
  {
    IntVectSet tags;

    IntVect center = 10*IntVect::Unit;

    setCircleTags(tags, 6, 1, center); //circle/sphere

    buildDisjointBoxLayout(plan1, tags, domain);

    tags.makeEmpty();

    setCircleTags(tags, 5, 2, center);

    buildDisjointBoxLayout(plan2, tags, domain);
  }
  if ( verbose )
  {
    pout() << "plan1: " << procID() << "...." << plan1 << endl;
    pout() << "plan2: " << procID() << "...." << plan2 << endl;
  }

  //test LayoutData<Real> specialization
  LayoutData<Real>  specialReal(plan1);


  LayoutData<Moment>   vlPlan(plan1);


  LevelData<BaseFab<int> > level1(plan1, 3, IntVect::Unit);

  LevelData<BaseFab<int> > level2;

  level2.define(level1);

  level2.define(plan2, 1);

  for (DataIterator i(level2.dataIterator()); i.ok(); ++i)
  {
    level2[i()].setVal(2);
  }

  level1.apply(values::setVal1);

  level2.apply(values::setVal2);

  HDF5HeaderData set1;
  Real dx=0.004;
  Box b1(IntVect(D_DECL6(1,2,1,1,2,1)), IntVect(D_DECL6(4,4,4,4,4,4)));
  Box b2(IntVect(D_DECL6(5,2,1,5,2,1)), IntVect(D_DECL6(12,4,4,12,4,4)));
  int currentStep = 2332;

  set1.m_string["name"] = "set1";
  set1.m_real["dx"] = dx;
  set1.m_int["currentStep"] = currentStep;
  set1.m_intvect["some intvect or other"] = b1.smallEnd();
  set1.m_box["b1"] = b1;
  set1.m_box["b2"] = b2;

  testFile.setGroupToLevel(1);

  error = write(testFile, plan1);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box write failed "<<error<<endl;
      return error;
    }

  error = write(testFile, level1, "level1 state vector");
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "BoxLayoutData 1 write failed "<<error<<endl;
      return error;
    }

  testFile.setGroupToLevel(2);

  error = write(testFile, plan2);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box2 write failed "<<error<<endl;
      return error;
    }

  error = write(testFile, level2, "level2 state vector");
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "BoxLayoutData 2 write failed "<<error<<endl;
      return error;
    }

  LevelData<FArrayBox> state(plan2, 3);

  state.apply(values::setVal3);
  testFile.setGroupToLevel(0);
  set1.writeToFile(testFile);

  error = write(testFile, plan2);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box2 write failed "<<error<<endl;
      return error;
    }

  testFile.setGroup("/");

  error = writeLevel(testFile, 0, state, 2, 1, 0.001, b2, 2);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "BoxLayoutData 2 write failed "<<error<<endl;
      return error;
    }

  set1.writeToFile(testFile);
  set1.writeToFile(testFile);

  testFile.close();

 CH_assert(!testFile.isOpen());

  // test the utility functions ReadUGHDF5 and WriteUGHDF5

  WriteUGHDF5("UGIO.hdf5", plan2, state, domain);

  ReadUGHDF5("UGIO.hdf5", plan2, state, domain);
  //========================================================================
  //
  //  now, read this data back in
  //
  //========================================================================

  BoxLayoutData<BaseFab<int> > readlevel1, readlevel2;

  error = testFile.open("data.h5", HDF5Handle::OPEN_RDONLY);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "File open failed "<<error<<endl;
      return error;
    }

  testFile.setGroupToLevel(2);
  Vector<Box> boxes;
  error = read(testFile, boxes);

  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box read failed "<<error<<endl;
      return error;
    }
  boxes.sort();
  Vector<int> assign;
  error = LoadBalance(assign, boxes);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "BoxLayout LoadBalance failed "<<error<<endl;
      return error;
    }
  BoxLayout readplan2(boxes, assign);
  readplan2.close();
  error = read(testFile, readlevel2, "level2 state vector", readplan2);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "BoxLayoutData<BaseFab<int>> read failed "<<error<<endl;
      return error;
    }

  testFile.setGroupToLevel(1);
  error = read(testFile, boxes);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box read failed "<<error<<endl;
      return error;
    }

  error = LoadBalance(assign, boxes);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "BoxLayout LoadBalance failed "<<error<<endl;
      return error;
    }
  BoxLayout readplan1(boxes, assign);
  readplan1.close();
  if ( verbose )
  {
    pout() << "readplan1: " << procID() << "...." << readplan1 << endl;
    pout() << "readplan2: " << procID() << "...." << readplan2 << endl;
  }
  error = read(testFile, readlevel1, "level1 state vector", readplan1);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "BoxLayoutData<BaseFab<int>> read failed "<<error<<endl;
      return error;
    }

  if ( verbose )
    pout() << plan1<<readplan1<<endl;

  // real test of IO, make sure the data is the same coming and going

  DataIterator l1  = level1.dataIterator();
  DataIterator rl1 = readlevel1.dataIterator();
  DataIterator l2  = level2.dataIterator();
  DataIterator rl2 = readlevel2.dataIterator();

  if (level1.boxLayout().size() != readlevel1.boxLayout().size())
    {
      if ( verbose )
        pout() << indent2 << "level1.size() != readl1.size() read failed "<<error<<endl;
      return 1;
    }
  if (level2.boxLayout().size() != readlevel2.boxLayout().size())
    {
      if ( verbose )
        pout() << indent2 << "level2.size() != readl2.size() read failed "<<error<<endl;
      return 1;
    }

  // we can assume that BoxLayout IO is tested in HDF5boxIO
  BaseFab<int>* before, *after;
  for (; l1.ok(); ++l1, ++rl1)
    {

      before = &(level1[l1()]); after = &(readlevel1[rl1()]);
      for (int c=0; c<before->nComp(); ++c)
        {
          for (BoxIterator it(level1.box(l1())); it.ok(); ++it)
            {
              if ((*before)(it(), c) != (*after)(it(), c))
                {
                  if ( verbose )
                    pout() << indent2 << "l1 != readl1 read failed "<<error<<endl;
                  return 2;
                }
            }
        }
    }

  for (; l2.ok(); ++l2, ++rl2)
    {

      before = &(level2[l2()]); after = &(readlevel2[rl2()]);
      for (int c=0; c<before->nComp(); ++c)
        {
          for (BoxIterator it(level2.box(l2())); it.ok(); ++it)
            {
              if ((*before)(it(), c) != (*after)(it(), c))
                {
                  if ( verbose )
                    pout() << indent2 << "level2 != readlevel2 read failed "<<error<<endl;
                  return 3;
                }
            }
        }
    }

  LevelData<FArrayBox> readState;
  Real dt, time;
  int refRatio;

  testFile.setGroup("/");

  error = readLevel(testFile, 0, readState, dx, dt, time, b2, refRatio);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "readLevel failed "<<error<<endl;
      return error;
    }

#ifndef CH_MPI
  // OK, now try to read one FArrayBox at a time
  // problem with DataIterator and running the out-of-core in parallel, so
  // have to think about that for now  BVS.
  FArrayBox readFAB;
  Interval  interval(1,2);
  testFile.setGroup("/");
  int index=0;
  for (DataIterator dit(state.dataIterator()); dit.ok(); ++index,++dit)
        {
          FArrayBox& fab = state[dit()];
          readFArrayBox(testFile, readFAB, 0, index, interval);
          for (BoxIterator it(state.box(dit())) ; it.ok() ; ++it)
                {
                  if (readFAB(it(), 0) != fab(it(), 1))
                        {
                          if ( verbose )
                                pout() << indent2 << "state != after for out-of-core "<<error<<endl;
                          return 3;
                        }
                }
        }

#endif

  testFile.close();

 CH_assert(!testFile.isOpen());

#endif // CH_USE_HDF5

  return 0;
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


//#endif
