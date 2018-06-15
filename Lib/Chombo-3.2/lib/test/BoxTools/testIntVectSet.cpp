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

using std::endl;

#include "REAL.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "BoxIterator.H"
#include "TreeIntVectSet.H"
#include "ProblemDomain.H"
#include "MayDay.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testIntVectSet();

/// Global variables for handling output:
static const char *pgmname = "testIntVectSet" ;
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
  int ret = testIntVectSet() ;

  if ( ret == 0 )
    {
      pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed one or more tests" << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

int
testIntVectSet()
{
  bool eekflag = false;

  {
    // dan's silly test
    Box testBox(IntVect::Zero, 3*IntVect::Unit);
    ProblemDomain testDom(testBox);
    for (int dir=0; dir<SpaceDim; dir++)
      {
        testDom.setPeriodic(dir,true);
      }

    int nestingRadius = 1;
    IntVectSet testIVS;
    testIVS |= testBox;
    testIVS.nestingRegion(nestingRadius, testDom);

    TreeIntVectSet tivs;
    Box set(3*IntVect::Unit, 4*IntVect::Unit);
    tivs |= set;

    for (BoxIterator bit(set); bit.ok(); ++bit)
      {
        IntVect image = bit();
        if (!testDom.image(image))
        {
          MayDay::Error(" A point managed to fall outside fully periodic");
        }
        if (!testBox.contains(image))
        {
          MayDay::Error("image point falls outside domain box");
        }
      }
    if (tivs.numPts() != D_TERM6(2,*2,*2,*2,*2,*2))
    {
      MayDay::Error("points went missing on ProblemDomain intersection");
    }

    if ( verbose )
          pout() << indent << pgmname
               << ": IntVectSet passed periodic ProblemDomain test"
               << endl;

  }
  {
    IntVectSet ivs;
    IntVect a(-1*IntVect::Unit), b(0*IntVect::Unit);
    ivs |= a;
    if (!ivs.contains(a))
      {
        pout() << indent << pgmname
             << ": failed -1 -1 -1 add test"
             << endl;
        eekflag = true;
      }
    if (ivs.numPts() != 1)
      {
        pout() << indent << pgmname
             << ": failed single IntVect count test"
             << endl;
        eekflag = true;
      }
    ivs-=a;
    if (ivs.contains(a))
      {
        pout() << indent << pgmname
             << ": failed [-1 -1 -1] remove test"
             << endl;
        eekflag = true;
      }
    if (ivs.numPts() != 0 || !ivs.isEmpty())
      {
        pout() << indent << pgmname
             << ": failed single IntVect removal count test"
             << endl;
        eekflag = true;
      }
    ivs|=a;
    ivs|=b;
    if (!ivs.contains(a) || !ivs.contains(b))
      {
        pout() << indent << pgmname
             << ": failed [-1 -1 -1] [0,0,0] add test"
             << endl;
        eekflag = true;
      }
    if (ivs.numPts() != 2)
      {
        pout() << indent << pgmname
             << ": failed two IntVect count test"
             << endl;
        eekflag = true;
      }
  }

#ifndef __IBMCPP__  // not working on AIX with xlC compiler
   // this does nothing, but makes sure dumpTree is avaiable for debugging work
  dumpTree(NULL);
#endif

  Box b1(12 * IntVect::Unit,
         31 * IntVect::Unit);
  Box b2(32 * IntVect::Unit,
         53 * IntVect::Unit);

  {
    // add, remove box.  success = empty IntVectSet.
    IntVectSet ivs;
    ivs |= b1;
    ivs |= b2;

    ivs -= b1;
    ivs -= b2;

    if (ivs.isEmpty())
      {
        if ( verbose )
          pout() << indent << pgmname
               << ": IntVectSet passed the add/remove Box test"
               << endl;
      }
    else
      {
        pout() << indent << pgmname
             << ": IntVectSet FAILED the add/remove Box test"
             << endl;
        eekflag = true;
      }
  }

  // intersect with box, remove box.  success = empty IntVectSet.
  {
    IntVectSet ivs;
    ivs |= b1;
    ivs |= b2;

    ivs &= b1;
    ivs -= b1;

    if (ivs.isEmpty())
      {
        if ( verbose )
          pout() << indent << pgmname
               << ": IntVectSet passed the intersect/remove Box test"
               << endl;
      }
    else
      {
        pout() << indent << pgmname
             << ": IntVectSet FAILED the intersect/remove Box test"
             << endl;
        eekflag = true;
      }
  }

  // add box.  success = every IntVect is in the IVS.
  {
    IntVectSet ivs;
    ivs |= b1;
    BoxIterator bit(b1);
    Vector<IntVect> iv_fail;
    for (bit.begin(); bit.ok(); bit.next())
      {
        IntVect iv = bit();
        if (!ivs.contains(iv))
          {
            iv_fail.push_back(iv);
          }
      }
    if (iv_fail.size() == 0)
      {
        if ( verbose )
          pout() << indent << pgmname
               << ": IntVectSet passed the every IntVect in the Box test"
               << endl;
      }
    else
      {
        pout() << indent << pgmname
             << ": IntVectSet FAILED the every IntVect in the Box test"
             << endl;
        eekflag = true;
      }
  }

  // intersect disjoint IntVectSets.  success = empty IntVectSet.
  {
    IntVectSet ivs1(b1);
    IntVectSet ivs2(b2);
    IntVectSet ivs_intersect = ivs1 & ivs2;

    if (ivs_intersect.isEmpty())
      {
        if ( verbose )
          pout() << indent << pgmname
               << ": IntVectSet passed the disjoint intersection test"
               << endl;
      }
    else
      {
        pout() << indent << pgmname
             << ": IntVectSet FAILED the disjoint intersection test"
             << endl;
        eekflag = true;
      }
  }
  // remove  one intvect from box.
  {
    IntVectSet ivs(b1);
    IntVect iv = (b1.smallEnd() + b1.bigEnd() ) / 2;
    ivs -= iv;
    BoxIterator bit(b1);
    Vector<IntVect> iv_fail;
    for (bit.begin(); bit.ok(); bit.next())
      {
        if (iv == bit())
        {
        }
        else if (!ivs.contains(bit()))
          {
            iv_fail.push_back(bit());
          }
      }
    if (iv_fail.size() == 0 && !ivs.contains(iv))
      {
        if ( verbose )
          pout() << indent << pgmname
               << ": IntVectSet passed the remove IntVect test"
               << endl;
      }
    else
      {
        pout() << indent << pgmname
             << ": IntVectSet FAILED the remove IntVect test"
             << endl;
        eekflag = true;
      }
  }

  // remove and add one intvect from box.  success = entire box is in IntVectSet.
  {
    IntVectSet ivs(b1);
    IntVect iv = (b1.smallEnd() + b1.bigEnd() ) / 2;
    ivs -= iv;
    ivs |= iv;
    BoxIterator bit(b1);
    Vector<IntVect> iv_fail;
    for (bit.begin(); bit.ok(); bit.next())
      {
        IntVect iv = bit();
        if (!ivs.contains(iv))
          {
            iv_fail.push_back(iv);
          }
      }
    if (iv_fail.size() == 0)
      {
        if ( verbose )
          pout() << indent << pgmname
               << ": IntVectSet passed the remove/add IntVect test"
               << endl;
      }
    else
      {
        pout() << indent << pgmname
             << ": IntVectSet FAILED the remove/add IntVect test"
             << endl;
        eekflag = true;
      }
  }

  {
    Box b1(IntVect(D_DECL6(13,5,-1, 13, 5, -1)),
           31 * IntVect::Unit);
    Box b2(32 * IntVect::Unit,
           53 * IntVect::Unit);
    Box b3(14 * IntVect::Unit,
           53 * IntVect::Unit);

    int numpts = b2.numPts()+b3.numPts() + b1.numPts() -(b2&b3).numPts() -
      (b2&b1).numPts() - (b1&b3).numPts();

    IntVectSet ivs2(b3);

    ivs2 |= b2;

    ivs2 |= b1;

    int count=0;

    IntVectSet otherSet(ivs2);

    IVSIterator it(otherSet), it2(ivs2);

    for (it.begin(); it.ok(); ++it)
      {
        if (!(b2.contains(it()) || b1.contains(it()) || b3.contains(it())))
          {
            if ( verbose )
              {
                pout() << indent2 << "union point falls outside region: ";
                pout() << pgmname << endl ;
              }
            eekflag = true;
          }
        count++;
      }
    if (count != numpts)
      {
        pout() << indent2 << "count != numpts from union ";
        pout() << pgmname << endl ;
        eekflag = true;
      }
    else
      {
        if ( verbose )
          pout() << indent << pgmname
               << ": IntVectSet passed the union domain test"
               << endl;
      }

    IntVect shifter(D_DECL6(1,-2,3,1,-2,3));
    ivs2.shift(shifter);
    b1.shift(shifter);
    b2.shift(shifter);
    b3.shift(shifter);
    bool fail = false ;
    for (it2.begin(); it2.ok(); ++it2)
      {
        if (!(b2.contains(it2()) || b1.contains(it2()) || b3.contains(it2())))
          {
            if ( verbose )
              {
                pout() << indent2 << "union point falls outside region: shift";
                pout() << pgmname << endl ;
              }
            eekflag = true;
            fail = true ;
          }
      }
    if ( !fail && verbose )
      {
        pout() << indent << pgmname
             << ": IntVectSet passed the shift domain test"
             << endl;
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
            eekflag = true;
          }
        count++;
      }
    if (count != numpts)
      {
        pout() << indent2 << "count != numpts from union ";
        pout() << pgmname << endl ;
        eekflag = true;
      }
    else
      {
        if ( verbose )
          {
            pout() << indent << pgmname
                 << ": IntVectSet passed the refine domain test"
                 << endl;
          }
      }
  }

  {
    IntVect a(D_DECL6(-2,-1,0,-2,-1,0)),
      b(D_DECL6(0,1,1,0,1,1)),
      c(D_DECL6(14,3,4,14,3,4)),
      d(D_DECL6(16,6,5,16,6,5));
    Box b1(a,c);
    Box b2(b,d);
    Box b3(b,c);
    IntVectSet ivs;
    ivs|= b1;
    ivs&= b2;

    long count = 0;
    for (IVSIterator it(ivs); it.ok(); ++it)
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
        eekflag = true;
      }
    else
      {
        if ( verbose )
          {
            pout() << indent << pgmname
                 << ": IntVectSet passed the intersection domain test"
                 << endl;
          }
      }

    bool fail = false ;
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
                eekflag = true;
                fail = true ;

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
                eekflag = true;
                fail = true ;

              }
          }
      }
    if ( !fail && verbose )
      {
        pout() << indent << pgmname
             << ": IntVectSet passed the intersection tests"
             << endl;
      }
  }

  {
    IntVect a(D_DECL6(-2,-1,0,-2,-1,0)),
      b(D_DECL6(0,1,1,0,1,1)),
      c(D_DECL6(14,3,4,14,3,4)),
      d(D_DECL6(16,6,5,16,6,5));
    Box b1(a,c);
    Box b2(b,d);
    Box b3(b,c);
    TreeIntVectSet ivs, ivs2;
    ivs2  |= b1;
    ivs   |= b2;

    ivs&=ivs2;

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
        eekflag = true;
      }
    else if (  verbose )
      {
        pout() << indent << pgmname
             << ": IntVectSet passed the TreeIntVectSet intersection domain test"
             << endl;
      }

    BoxIterator bit(b1);
    bool fail = false ;
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
                eekflag = true;
                fail = true ;

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
                eekflag = true;
                fail = true ;

              }
          }
      }
    if ( !fail && verbose )
      {
        pout() << indent << pgmname
             << ": IntVectSet passed the TreeIntVectSet intersection tests"
             << endl;
      }

    {
      int s = ivs.linearSize();
      char* buffer = (char*)malloc(s);
      ivs.linearOut(buffer);
      TreeIntVectSet ivsIn;
      ivsIn.linearIn(buffer);

      if (!(ivs == ivsIn))
        {
          if ( verbose )
            {
              pout() << indent2 << "TreeIntVectSet linearization failed";
              pout() << pgmname << endl ;
            }
          eekflag = true;

        }
      else if ( verbose )
        {
          pout() << indent << pgmname
               << ": IntVectSet passed the TreeIntVectSet linearization test"
               << endl;
        }

      free(buffer);
    }
    {
      b3.shift(0,4);
      ivs |= b3;
      int s = ivs.linearSize();
      char* buffer = (char*)malloc(s);
      ivs.linearOut(buffer);
      TreeIntVectSet ivsIn;
      ivsIn.linearIn(buffer);

      if (!(ivs == ivsIn))
        {
          if ( verbose )
            {
              pout() << indent2 << "TreeIntVectSet linearization failed";
              pout() << pgmname << endl ;
            }
          eekflag = true;

        }
      else if ( verbose )
        {
          pout() << indent << pgmname
               << ": IntVectSet passed the TreeIntVectSet union linearization test"
               << endl;
        }

      free(buffer);
    }

    {
      TreeIntVectSet emptyIVS;
      ivs = emptyIVS;

      int s = ivs.linearSize();
      char* buffer = (char*)malloc(s);
      ivs.linearOut(buffer);
      TreeIntVectSet ivsIn;
      ivsIn.linearIn(buffer);

      if (!(ivs == ivsIn))
        {
          if ( verbose )
            {
              pout() << indent2 << "TreeIntVectSet linearization failed";
              pout() << pgmname << endl ;
            }
          eekflag = true;

        }
      else if ( verbose )
        {
          pout() << indent << pgmname
               << ": IntVectSet passed the empty TreeIntVectSet linearization test"
               << endl;
        }

      free(buffer);
    }
  }

  // test IntVectSet::grow
  {
    bool fail = false;
    Box smallBox(IntVect::Zero, 2*IntVect::Unit);
    int numGrow = 2;
    Box grownBox(smallBox);
    grownBox.grow(numGrow);

    {
      // this should be a dense IVS
      IntVectSet ivs(smallBox);
      ivs.grow(numGrow);

      // success -- ivs == grownBox
      if (ivs.minBox() != grownBox)
        {
          if ( verbose )
            {
              pout() << indent2 << "DenseIntVectSet grow test failed minBox test"
                   << pgmname << endl ;
            }
          eekflag = true;
          fail = true;
        }

      BoxIterator bit(grownBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          if (!ivs.contains(bit()))
            {
              if ( verbose )
                {
                  pout() << indent2 << "fail DenseIntVectSet grow test "
                       << pgmname << endl ;
                }
              eekflag = true;
              fail = true;
            }
        }

      if (verbose && !fail)
        {
          pout() << indent << pgmname
               << " passed DenseIntVectSet::grow test" << endl;
        }
    }


    {
      // this should be a tree IVS
      IntVectSet ivs;
      ivs |=  smallBox;
      ivs.grow(numGrow);

      // success -- ivs == grownBox
      if (ivs.minBox() != grownBox)
        {
          if ( verbose )
            {
              pout() << indent2 << "TreeIntVectSet grow test failed minBox test"
                   << pgmname << endl ;
            }
          eekflag = true;
          fail = true;
        }

      BoxIterator bit(grownBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          if (!ivs.contains(bit()))
            {
              if ( verbose )
                {
                  pout() << indent2 << "fail TreeIntVectSet grow test "
                       << pgmname << endl ;
                }
              eekflag = true;
              fail = true;
            }
        }

      if (verbose && !fail)
        {
          pout() << indent << pgmname
               << " passed TreeIntVectSet::grow test" << endl;
        }
    }


  }


  // make sure no strange destructor behaviour floats around....
  {
    Box b1(12 * IntVect::Unit,
           31 * IntVect::Unit);
    Box b2(32 * IntVect::Unit,
           53 * IntVect::Unit);
    Box b3(14 * IntVect::Unit,
           51 * IntVect::Unit);

    IntVectSet mask1(b3);

    mask1 -= b2;
    mask1 -= b1;

    IntVectSet mask3(mask1);
    IntVectSet mask4 = mask3; //ensure copy and assignment of dense works

    for (IVSIterator it(mask3), it4(mask4); it.ok(); ++it, ++it4)
      {
        if (b2.contains(it()) || b1.contains(it()))
          {
            if ( verbose )
              {
                pout() << indent2 << "subtraction left in bad data: ";
                pout() << pgmname << endl ;
              }
            eekflag = true;
          }

        if (!b3.contains(it()))
          {
            if ( verbose )
              {
                pout() << indent2 << "IntVect beyond domain: ";
                pout() << pgmname << endl ;
              }
            eekflag = true;
          }
        if (it() != it4())
          {
            if ( verbose )
              {
                pout() << indent2 << "IntVectSet copy/assign problem: ";
                pout() << pgmname << endl ;
              }
            eekflag = true;
          }
      }

    IntVectSet mask2(b3);
    mask2 &= b1;

    for (IVSIterator it(mask2); it.ok(); ++it)
      {
        if (!(b3.contains(it()) && b1.contains(it())))
          {
            if ( verbose )
              {
                pout() << indent2 << "intersection contains extra data: ";
                pout() << pgmname << endl ;
              }
            eekflag = true;
          }
      }
    for (BoxIterator bit(b1&b3); bit.ok(); ++bit)
      {
        if (!mask2.contains(bit()))
          {
            if ( verbose )
              {
                pout() << indent2 << "intersection missed data: ";
                pout() << pgmname << endl ;
              }
            eekflag = true;
          }
      }
    {
      int s = mask3.linearSize();
      char* buffer = (char*)malloc(s);
      mask3.linearOut(buffer);
      IntVectSet ivsIn;
      ivsIn.linearIn(buffer);

      if (!(mask3 == ivsIn))
        {
          if ( verbose )
            {
              pout() << indent2 << "DenseIntVectSet linearization failed";
              pout() << pgmname << endl ;
            }
          eekflag = true;

        }

      free(buffer);
    }
    {
      mask3.makeEmpty();
      int s = mask3.linearSize();
      char* buffer = (char*)malloc(s);
      mask3.linearOut(buffer);
      IntVectSet ivsIn;
      ivsIn.linearIn(buffer);

      if (!(mask3 == ivsIn))
        {
          if ( verbose )
            {
              pout() << indent2 << "DenseIntVectSet linearization failed";
              pout() << pgmname << endl ;
            }
          eekflag = true;

        }

      free(buffer);
    }

  }

  {
    Box b1(IntVect(D_DECL6(13,5,-1,13,5,-1)),
           31 * IntVect::Unit);
    Box b2(32 * IntVect::Unit,
           53 * IntVect::Unit);
    Box b3(14 * IntVect::Unit,
           51 * IntVect::Unit);

    IntVectSet ivs2;

    ivs2 |= b3;

    ivs2 |= b2;

    ivs2 |= b1;

    b1.coarsen(2); b2.coarsen(2); b3.coarsen(2);
    ivs2.coarsen(2);

    int count = 0;
    int numpts = b2.numPts()+b3.numPts() + b1.numPts() -(b2&b3).numPts() -
      (b2&b1).numPts() - (b1&b3).numPts();

    IVSIterator it(ivs2);
    for (it.begin(); it.ok(); ++it)
      {
        if (!(b2.contains(it()) || b1.contains(it()) || b3.contains(it())))
          {
            if ( verbose )
              {
                pout() << indent2 << "box.coarsen != ivs.coarsen: point falls outside region: ";
                pout() << pgmname << endl ;
              }
            eekflag = true;
          }
        count++;
      }
    if (count != numpts)
      {
        pout() << indent2 << "box.coarsen != ivs.coarsen: count != numpts from union ";
        pout() << pgmname << endl ;
        eekflag = true;
      }

  }
  {
    Box b1(IntVect(D_DECL6(13,5,-1,13,5,-1)),
           31 * IntVect::Unit);
    Box b2(32 * IntVect::Unit,
           53 * IntVect::Unit);
    Box b3(14 * IntVect::Unit,
           51 * IntVect::Unit);

    TreeIntVectSet ivs2, ivs3, ivs4;

    ivs4.clear();
    ivs4.clear();

    ivs2 |= b3;

    ivs2 |= b2;

    ivs3 |= b1;

        ivs2 |= ivs3;

    b1.coarsen(2); b2.coarsen(2); b3.coarsen(2);

    ivs2.coarsen(2);

    int count = 0;
    int numpts = b2.numPts()+b3.numPts() + b1.numPts() -(b2&b3).numPts() -
      (b2&b1).numPts() - (b1&b3).numPts();

    TreeIntVectSetIterator it(ivs2);
    for (it.begin(); it.ok(); ++it)
      {
        if (!(b2.contains(it()) || b1.contains(it()) || b3.contains(it())))
          {
            if ( verbose )
              {
                pout() << indent2 << "box.coarsen != tivs.coarsen: point falls outside region: ";
                pout() << pgmname << endl ;
              }
            eekflag = true;
          }
        count++;
      }
    if (count != numpts)
      {
        pout() << indent2 << "box.coarsen != tivs.coarsen: count != numpts from union ";
        pout() << pgmname << endl ;
        eekflag = true;
      }
    ivs2.clear();
    ivs2.clear();
  }
  {
    Box b1(IntVect(D_DECL6(5,5,5,5,5,5)),
           10 * IntVect::Unit);
    Box b2(6 * IntVect::Unit,
           7 * IntVect::Unit);
    Box b3(IntVect(D_DECL6(6,8,8,6,8,8)),
           IntVect(D_DECL6(7,9,8,7,9,8)));
    Box b4(IntVect(D_DECL6(8,6,5,8,6,5)),
           IntVect(D_DECL6(9,7,5,9,7,5)));
    Box b5(IntVect(D_DECL6(8,8,8,8,8,8)),
           IntVect(D_DECL6(9,9,9,9,9,9)));

    TreeIntVectSet ivs1, ivs2;

    ivs1 |= b1;

    ivs2 |= b2;
    ivs2 |= b3;
    ivs2 |= b4;
    ivs2 |= b5;

    TreeIntVectSet ivs3(ivs2);

    ivs2 &= ivs1;
    TreeIntVectSetIterator it(ivs3);
    for (it.begin(); it.ok(); ++it)
      {
        if (!ivs2.contains(it()))
        {
          MayDay::Error ("TreeIntVectSet::&= operator broken");
        }
      }
  }
  {
    /*   TreeIntVectSet ivs1;
    ivs1 |= IntVect::Zero;
    ivs1.coarsen(4);
    // pout() << "ivs1:\n" << ivs1 << std::endl;
    TreeIntVectSet ivs2;
    ivs2 |= 32*IntVect::Unit;
    ivs2.coarsen(4);
    //pout() << "ivs2:\n" << ivs2 << std::endl;

    ivs1 |= ivs2; // should contain just (0,0) and (8,8) afterwards
    if (ivs1.numPts() != 2)
      {
        pout() <<"failed simple IntVectSet coarsen and union\n";
        eekflag = true;
      }
    if (!ivs1.contains(IntVect::Zero) || !ivs1.contains(8*IntVect::Unit))
      {
        pout() <<"failed simple IntVectSet coarsen and union\n";
        eekflag = true;
      }
    */
  }
    {
     // this test appears to be breaking for Donna
     IntVectSet testIVS;
     Box testIVSBox(IntVect::Zero, 31*IntVect::Unit);
     // force a tree IVS
     testIVS |= testIVSBox;

     // create periodic domain box (shouldn't be an issue, but
     // am duplicating the test)
     Box domBox(IntVect::Zero, 63*IntVect::Unit);

     ProblemDomain thisDom(domBox);
     for (int dir=0; dir<SpaceDim; dir++)
       {
         thisDom.setPeriodic(dir,true);
       }

     int bufferSize = 1;
     testIVS.nestingRegion(bufferSize, thisDom);

     // check results
     Box nestedBox(IntVect::Unit, 30*IntVect::Unit);

     if (testIVS.numPts() != nestedBox.numPts())
       {
         eekflag = true;
       }
     testIVS -= nestedBox;
     if (!testIVS.isEmpty())
       {
         eekflag = true;
       }
   }

  if (eekflag) return 1;
  else        return 0;
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
