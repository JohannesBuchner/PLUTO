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
using std::sqrt;
#include <iostream>
using std::endl;

#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "BRMeshRefine.H"
#include "IntVectSet.H"
#include "LoadBalance.H"
#include "TestCommon.H"
#include "UsingNamespace.H"


void
parseTestOptions( int argc ,char* argv[] ) ;

int
test();

/// Global variables for handling output:
static const char *pgmname = "HDF5boxIO" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
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

// returns 0 on all tests passed.

int test()
{
#ifdef CH_USE_HDF5
  int error;
  HDF5Handle testFile;

  CH_assert(!testFile.isOpen());

  error = testFile.open("boxIO.h5", HDF5Handle::CREATE);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "File creation failed "<<error<<endl;
      return error;
    }

 CH_assert(testFile.isOpen());

  DisjointBoxLayout plan1, plan2;
  {
    Box domain(IntVect::Zero, 20*IntVect::Unit);

    IntVectSet tags;

    IntVect center = 10*IntVect::Unit;

    setCircleTags(tags, 6, 3, center); //circle/sphere

    buildDisjointBoxLayout(plan1, tags, domain);

    tags.makeEmpty();

    setCircleTags(tags, 5, 2, center);

    buildDisjointBoxLayout(plan2, tags, domain);
  }

  testFile.setGroupToLevel(0);
  error = write(testFile, plan1);
   if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box write failed "<<error<<endl;
      return error;
    }

  testFile.setGroupToLevel(1);

  error = write(testFile, plan2);
    if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box write failed "<<error<<endl;
      return error;
    }

  testFile.close();

 CH_assert(!testFile.isOpen());

  //====================================================================

  DisjointBoxLayout readplan1, readplan2;

  testFile.open("boxIO.h5", HDF5Handle::OPEN_RDONLY);

  testFile.setGroupToLevel(1);
  Vector<Box> boxes;
  Vector<int> assignments;
  error = read(testFile, boxes);
  if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box read failed "<<error<<endl;
      return error;
    }
  LoadBalance(assignments, boxes);
  readplan2.define(boxes, assignments);

  testFile.setGroupToLevel(0);

  error = read(testFile, boxes);
   if (error != 0)
    {
      if ( verbose )
        pout() << indent2 << "box read failed "<<error<<endl;
      return error;
    }
   LoadBalance(assignments, boxes);
   readplan1.define(boxes, assignments);

   // now the test:
   LayoutIterator p1 = plan1.layoutIterator();
   LayoutIterator rp1 = readplan1.layoutIterator();
   LayoutIterator p2 = plan2.layoutIterator();
   LayoutIterator rp2 = readplan2.layoutIterator();

   if (plan1.size() != readplan1.size())
     {
       if ( verbose )
         pout() << indent2 << "plan1 size different on read "<<endl;
       return 1;
     }

   for (;p1.ok(); ++p1, ++rp1)
     {
       if (plan1.get(p1()) != readplan1.get(rp1()))
         {
           if ( verbose )
             {
               pout() << indent2 << "plan1 != readplan1  on read \n";
               pout()<<plan1<<"\n\n"<<readplan1<<endl;
             }
           return 2;
         }

     }

   if (plan2.size() != readplan2.size())
     {
       if ( verbose )
         pout() << indent2 << "plan2 size different on read "<<endl;
       return 1;
     }

   for (;p2.ok(); ++p2, ++rp2)
     {
       if (plan2.get(p2()) != readplan2.get(rp2()))
         {
           if ( verbose )
             {
               pout() << indent2 << "plan2 != readplan2  on read \n";
               pout()<<plan2<<"\n\n"<<readplan2<<endl;
             }
           return 2;

         }
     }

   testFile.close();
#endif
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

