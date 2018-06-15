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
#include "Box.H"
#include "Vector.H"
#include "parstream.H"
#include "TreeIntVectSet.H"
#include "BoxIterator.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"
using std::endl;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testTreeIntVectSet();

/// Global variables for handling output:
static const char *pgmname = "testTreeIntVectSet" ;
static const char *indent2 = "      " ;
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
  int ret = testTreeIntVectSet() ;
  if (ret == 0)
    pout() << indent2 << pgmname << " passed." << endl ;
  else
    pout() << indent2 << pgmname << " failed." << endl ;
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

int
testTreeIntVectSet()
{
  {
    Box b1(IntVect(D_DECL6(13,5,-1, -1, -1, -1)),
           31 * IntVect::Unit);
    Box b2(32 * IntVect::Unit,
           53 * IntVect::Unit);
    Box b3(14 * IntVect::Unit,
           53 * IntVect::Unit);

    /*
    TreeIntVectSet ivs(b);

    ivs -= b2;
    ivs -= b1;

    for (TreeIntVectSetIterator it(ivs); it.ok(); ++it)
      {
        if (b2.contains(it()) || b1.contains(it()))
          {
            if ( verbose )
              {
                pout() << indent2 << "subtraction left in bad data: ";
                pout() << pgmname << endl ;
              }
            return 1;
          }

        if (!b3.contains(it()))
          {
            if ( verbose )
              {
                pout() << indent2 << "IntVect beyond domain: ";
                pout() << pgmname << endl ;
              }
            return 1;
          }
      }
    */

    int numpts = b2.numPts()+b3.numPts() + b1.numPts() -(b2&b3).numPts() -
      (b2&b1).numPts() - (b1&b3).numPts();

    TreeIntVectSet ivs2(b3);

    ivs2 |= b2;

    ivs2 |= b1;

    pout() <<endl;
    int count=0;

    TreeIntVectSet otherSet(ivs2);

    TreeIntVectSetIterator it(otherSet), it2(ivs2);

    for (it.begin(); it.ok(); ++it)
      {
        if (!(b2.contains(it()) || b1.contains(it()) || b3.contains(it())))
          {
            if ( verbose )
              {
                pout() << indent2 << "union point falls outside region: ";
                pout() << pgmname << endl ;
              }
            return 2;
          }
        count++;
      }
    if (count != numpts)
      {
        pout() << indent2 << "count != numpts from union ";
        pout() << pgmname << endl ;
        return 2;
      }

    IntVect shifter(D_DECL6(1,-2,3,1,-2,3));
    ivs2.shift(shifter);
    b1.shift(shifter);
    b2.shift(shifter);
    b3.shift(shifter);
    for (it2.begin(); it2.ok(); ++it2)
      {
        if (!(b2.contains(it2()) || b1.contains(it2()) || b3.contains(it2())))
          {
            if ( verbose )
              {
                pout() << indent2 << "union point falls outside region: shift";
                pout() << pgmname << endl ;
              }
            return 2;
          }
      }

    b1.refine(4);
    b2.refine(4);
    b3.refine(4);
    ivs2.refine(4);
    numpts = b2.numPts()+b3.numPts() + b1.numPts() -(b2&b3).numPts() -
      (b2&b1).numPts() - (b1&b3).numPts();
    count = 0;
    for (it2.begin(); it2.ok(); ++it2)
      {
        if (!(b2.contains(it2()) || b1.contains(it2()) || b3.contains(it2())))
          {
            if ( verbose )
              {
                pout() << indent2 << "union point falls outside region:refine ";
                pout() << pgmname << endl ;
              }
            return 2;
          }
        count++;
      }
    if (count != numpts)
      {
        pout() << indent2 << "count != numpts from union ";
        pout() << pgmname << endl ;
        return 2;
      }

  }

  {
    IntVect a(D_DECL6(-2,-1,0,-2,-1,0)), b(D_DECL6(0,1,1,0,1,1)),
      c(D_DECL6(14,3,4,14,3,4)), d(D_DECL6(16,6,5,16,6,5));
    Box b1(a,c);
    Box b2(b,d);
    Box b3(b,c);
    TreeIntVectSet ivs;
    ivs|= b1;
    ivs&= b2;

    long count = 0;
    for (TreeIntVectSetIterator it(ivs); it.ok(); ++it)
      {
        if (!b3.contains(it()))
          {
            if ( verbose )
              {
                pout() << indent2 << "intersection point falls outside region:";
                pout() << pgmname << endl ;
              }
            return 2;
          }
        count++;
      }
    if (count != b3.numPts())
      {
        pout() << indent2 << "count != numpts from intersection ";
        pout() << pgmname << endl ;
        return 2;
      }

    BoxIterator bit(b1);
    for (bit.begin(); bit.ok(); ++bit)
      {
        if (b2.contains(bit()))
          {
            if (!ivs.contains(bit()))
              {
                if ( verbose )
                  {
                    pout() << indent2 << "intersection point missing";
                    pout() << pgmname << endl ;
                  }
              }
          }
        else
          {
            if (ivs.contains(bit()))
              {
                if ( verbose )
                  {
                    pout() << indent2 << "extraneous intersection point";
                    pout() << pgmname << endl ;
                  }
              }
          }
      }
  }
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
