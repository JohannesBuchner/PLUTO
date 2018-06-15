#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>

#include "REAL.H"
#include "Box.H"
#include "parstream.H"
#include "BoxIterator.H"
#include "RefCountedPtr.H"
#include "Vector.H"
#include "CH_Timer.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testBox();

/// Global variables for handling output:
static const char *pgmname = "testBox" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
  CH_TIMERS("main()");
  CH_TIMER("wholething", t1);
  CH_START(t1);

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  Vector<int> vi(10);
  RefCountedPtr< Vector<int> > rcp( new Vector<int>( vi ) );

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testBox() ;
#ifdef CH_MPI
  MPI_Finalize();
#endif

  CH_STOP(t1);
  CH_TIMER_REPORT();

  return ret;
}

void
breakpointHook()
{
}

int
testBox()
{
  {
    IntVect lo(IntVect::Unit), hi(IntVect::Unit*3);
    hi[0]-=1;
    Box b(lo, hi);
    b.refine(2);
    b.coarsen(2);
    b.surroundingNodes();
    b.enclosedCells();

    Box b2(b);
    b2.shift(IntVect::Unit);
    b &= b2;
  }
  int return_code = 0;

#if (CH_SPACEDIM == 1)
  const Box b1(IntVect(  0 ), IntVect(  7 ));
  const long numPtsB1 = 8;
  const Box b2(IntVect( 12 ), IntVect( 15 ));
  const long numPtsB2 = 4;

#if !CH_USE_64
  const long n = 1 << 15;
#else
  const long n = 1 << 30;
#endif
  const Box bBig(IntVect( 0 ), IntVect( n-1 ));
  const long numPtsBBig = n;
#elif (CH_SPACEDIM == 2)
  const Box b1(IntVect(  0,  0 ), IntVect(  7,  7 ));
  const long numPtsB1 = 8*8;
  const Box b2(IntVect( 12,  4 ), IntVect( 15, 15 ));
  const long numPtsB2 = 4*12;

#if !CH_USE_64
  const long n = 1 << 15;
#else
  const long n = 1 << 30;
#endif
  const Box bBig(IntVect( 0, 0 ), IntVect( n-1, n-1 ));
  const long numPtsBBig = n*n;
#elif (CH_SPACEDIM == 3)
  const Box b1(IntVect(  0,  0,  0 ), IntVect(  7,  7,  7 ));
  const long numPtsB1 = 8*8*8;
  const Box b2(IntVect( 12, 12,  4 ), IntVect( 15, 15, 15 ));
  const long numPtsB2 = 4*4*12;

#if !CH_USE_64
  const long n = 1 << 10;
#else
  const long n = 1 << 20;
#endif
  const Box bBig(IntVect( 0, 0, 0 ), IntVect( n-1, n-1, n-1 ));
  const long numPtsBBig = n*n*n;
#elif (CH_SPACEDIM == 4)
  const Box b1(IntVect(  0,  0, 0,  0 ), IntVect(  7,  7,  7, 7 ));
  const long numPtsB1 = 8*8*8*8;
  const Box b2(IntVect( 12, 12, 12, 4 ), IntVect( 15, 15, 15, 15 ));
  const long numPtsB2 = 4*4*4*12;

#if !CH_USE_64
  const long n = 1 <<  7;
#else
  const long n = 1 << 15;
#endif
  const Box bBig(IntVect( 0, 0, 0, 0 ), IntVect( n-1, n-1, n-1, n-1 ));
  const long numPtsBBig = n*n*n*n;
#elif (CH_SPACEDIM == 5)
  const Box b1(IntVect(  0,  0,  0,  0, 0 ), IntVect(  7,  7,  7,  7,  7));
  const long numPtsB1 = 8*8*8*8*8;
  const Box b2(IntVect( 12, 12, 12, 12, 4 ), IntVect( 15, 15, 15, 15, 15));
  const long numPtsB2 = 4*4*4*4*12;

#if !CH_USE_64
  const long n = 1 <<  6;
#else
  const long n = 1 << 12;
#endif
  const Box bBig(IntVect( 0, 0, 0, 0, 0 ), IntVect( n-1, n-1, n-1, n-1, n-1 ));
  const long numPtsBBig = n*n*n*n*n;
#elif (CH_SPACEDIM == 6)
  const Box b1(IntVect( 0,  0,  0,  0,  0, 0), IntVect( 7,  7,  7,  7,  7, 7));
  const long numPtsB1 = 8*8*8*8*8*8;
  const Box b2(IntVect(12, 12, 12, 12, 12, 4), IntVect(15, 15, 15, 15, 15,15));
  const long numPtsB2 = 4*4*4*4*4*12;

#if !CH_USE_64
  const long n = 1 <<  5;
#else
  const long n = 1 << 10;
#endif
  const Box bBig(IntVect( 0, 0, 0, 0, 0, 0 ), IntVect( n-1, n-1, n-1, n-1, n-1, n-1 ));
  const long numPtsBBig = n*n*n*n*n*n;
#else
#error invalid CH_SPACEDIM.  Must be 1, 2, 3, 4, 5, or 6.
#endif

  const long numPtsB1Computed = b1.numPts();
  const long numPtsB2Computed = b2.numPts();
  const long numPtsBBigComputed = bBig.numPts();

  if (verbose)
  {
    pout() << indent2 << "num points in box 1: " << numPtsB1Computed
                      << " (" << numPtsB1 << ")" << endl;
    pout() << indent2 << "num points in box 2: " << numPtsB2Computed
                      << " (" << numPtsB2 << ")" << endl;
    pout() << indent2 << "num points in box big: " << numPtsBBigComputed
                      << " (" << numPtsBBig << ")" << endl;
  }

  if (numPtsB1 != numPtsB1Computed)
  {
    if (verbose)
    {
      pout() << indent2 << "num points in box 1 FAILED!!!" << endl;
      pout() << indent2 << "num points in box 1: " << numPtsB1 << " (actual)"
                        << " != " << numPtsB1Computed << " (computed)" << endl;
    }
    return_code = -1;
  }

  if (numPtsB2 != numPtsB2Computed)
  {
    if (verbose)
    {
      pout() << indent2 << "num points in box 2 FAILED!!!" << endl;
      pout() << indent2 << "num points in box 2: " << numPtsB2 << " (actual)"
                        << " != " << numPtsB2Computed << " (computed)" << endl;
    }
    return_code = -1;
  }

  if (numPtsBBig != numPtsBBigComputed)
  {
    if (verbose)
    {
      pout() << indent2 << "num points in box big FAILED!!!" << endl;
      pout() << indent2 << "num points in box big: " << numPtsBBig << " (actual)"
                        << " != " << numPtsBBigComputed << " (computed)" << endl;
    }
    return_code = -1;
  }

//
// null construction tests
  const Box b_empty;
  if (verbose)
  {
    pout() << indent2 << "test box 1: " << b1 << endl;
    pout() << indent2 << "test box 2: " << b2 << endl;
    pout() << indent2 << "empty null box: " << b_empty << endl;
  }

  const bool b_empty_isempty = b_empty.isEmpty();
  if (!b_empty_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "null constructor isEmpty test FAILED!!!" << endl;
      pout() << indent2 << "b_empty_isempty = " << b_empty_isempty << endl;
    }
    return_code = -1;
  }
  const long b_empty_numpts = b_empty.numPts();
  if (b_empty_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "null constructor numPts test FAILED!!!" << endl;
      pout() << indent2 << "b_empty_numpts = " << b_empty_numpts << endl;
    }
    return_code = -1;
  }
//
// contains tests
  const bool b1_contains_empty = b1.contains(b_empty);
  if (!b1_contains_empty)
  {
    if (verbose)
    {
      pout() << indent2 << "contains empty test 1 FAILED!!!" << endl;
      return_code = -1;
    }
  }
//
  const bool b2_contains_empty = b2.contains(b_empty);
  if (!b2_contains_empty)
  {
    if (verbose)
    {
      pout() << indent2 << "contains empty test 2 FAILED!!!" << endl;
      return_code = -1;
    }
  }
//
  const bool empty_contains_iv0 = b_empty.contains(IntVect::Zero);
  if (empty_contains_iv0)
  {
    if (verbose)
    {
      pout() << indent2 << "empty contains iv test 1 FAILED!!!" << endl;
    }
    return_code = -1;
  }
//
  const bool empty_contains_iv1 = b_empty.contains(IntVect::Unit);
  if (empty_contains_iv1)
  {
    if (verbose)
    {
      pout() << indent2 << "empty contains iv test 2 FAILED!!!" << endl;
    }
    return_code = -1;
  }
//
  const bool empty_contains_empty = b_empty.contains(b_empty);
  if (!empty_contains_empty)
  {
    if (verbose)
    {
      pout() << indent2 << "empty contains empty test 1 FAILED!!!" << endl;
    }
    return_code = -1;
  }
//
// conversion to node centered tests
  Box beconv1 = b_empty;
  beconv1.convert(IndexType::TheNodeType());

  const bool beconv1_isempty = beconv1.isEmpty();
  if (!beconv1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "conversion to node isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "beconv1_isempty = " << beconv1_isempty << endl;
    }
    return_code = -1;
  }
  const long beconv1_numpts = beconv1.numPts();
  if (beconv1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "conversion to node numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "beconv1_numpts = " << beconv1_numpts << endl;
    }
    return_code = -1;
  }
//
  Box beconv2 = b_empty;
  beconv2.convert(IntVect::Unit);

  const bool beconv2_isempty = beconv2.isEmpty();
  if (!beconv2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "conversion to node isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "beconv2_isempty = " << beconv2_isempty << endl;
    }
    return_code = -1;
  }
  const long beconv2_numpts = beconv2.numPts();
  if (beconv2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "conversion to node numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "beconv2_numpts = " << beconv2_numpts << endl;
    }
    return_code = -1;
  }
//
  Box beconv3 = b_empty;
  beconv3.convert(0,IndexType::NODE);

  const bool beconv3_isempty = beconv3.isEmpty();
  if (!beconv3_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "conversion to node isEmpty test 3 FAILED!!!" << endl;
      pout() << indent2 << "beconv3_isempty = " << beconv3_isempty << endl;
    }
    return_code = -1;
  }
  const long beconv3_numpts = beconv3.numPts();
  if (beconv3_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "conversion to node numPts test 3 FAILED!!!" << endl;
      pout() << indent2 << "beconv3_numpts = " << beconv3_numpts << endl;
    }
    return_code = -1;
  }
//
// surrounding nodes tests
  const Box besn1 = surroundingNodes(b_empty);

  const bool besn1_isempty = besn1.isEmpty();
  if (!besn1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "surrounding nodes isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "besn1_isempty = " << besn1_isempty << endl;
    }
    return_code = -1;
  }
  const long besn1_numpts = besn1.numPts();
  if (besn1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "surrounding nodes numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "besn1_numpts = " << besn1_numpts << endl;
    }
    return_code = -1;
  }
//
  Box besn2(b_empty);
  besn2.surroundingNodes();

  const bool besn2_isempty = besn2.isEmpty();
  if (!besn2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "surrounding nodes isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "besn2_isempty = " << besn2_isempty << endl;
    }
    return_code = -1;
  }
  const long besn2_numpts = besn2.numPts();
  if (besn2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "surrounding nodes numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "besn2_numpts = " << besn2_numpts << endl;
    }
    return_code = -1;
  }
//
// enclosed cell tests
  const Box beec1 = enclosedCells(besn1);

  const bool beec1_isempty = beec1.isEmpty();
  if (!beec1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "enclosed cells isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "beec1_isempty = " << beec1_isempty << endl;
    }
    return_code = -1;
  }
  const long beec1_numpts = beec1.numPts();
  if (beec1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "enclosed cells numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "beec1_numpts = " << beec1_numpts << endl;
    }
    return_code = -1;
  }
//
  Box beec2(besn2);
  beec2.enclosedCells();

  const bool beec2_isempty = beec2.isEmpty();
  if (!beec2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "enclosed cells isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "beec2_isempty = " << beec2_isempty << endl;
    }
    return_code = -1;
  }
  const long beec2_numpts = beec2.numPts();
  if (beec2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "enclosed cells numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "beec2_numpts = " << beec2_numpts << endl;
    }
    return_code = -1;
  }
//
// intersection tests
  const Box bi1 = b1 & b2;

  const bool bi1_isempty = bi1.isEmpty();
  if (!bi1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "intersection isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "bi1_isempty = " << bi1_isempty << endl;
    }
    return_code = -1;
  }
  const long bi1_numpts = bi1.numPts();
  if (bi1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "intersection numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "bi1_numpts = " << bi1_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bi2(b1);
  bi2 &= b2;

  const bool bi2_isempty = bi2.isEmpty();
  if (!bi2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "intersection isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "bi2_isempty = " << bi2_isempty << endl;
    }
    return_code = -1;
  }
  const long bi2_numpts = bi2.numPts();
  if (bi2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "intersection numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "bi2_numpts = " << bi2_numpts << endl;
    }
    return_code = -1;
  }
//
// grow tests
  const Box bg1 = grow(b1,-4);

  const bool bg1_isempty = bg1.isEmpty();
  if (!bg1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "bg1_isempty = " << bg1_isempty << endl;
    }
    return_code = -1;
  }
  const long bg1_numpts = bg1.numPts();
  if (bg1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "bg1_numpts = " << bg1_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bg2(b1);
  bg2.grow(-4);

  const bool bg2_isempty = bg2.isEmpty();
  if (!bg2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "bg2_isempty = " << bg2_isempty << endl;
    }
    return_code = -1;
  }
  const long bg2_numpts = bg2.numPts();
  if (bg2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "bg2_numpts = " << bg2_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box bg3 = grow(b1, -4*BASISV(0));

  const bool bg3_isempty = bg3.isEmpty();
  if (!bg3_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 3 FAILED!!!" << endl;
      pout() << indent2 << "bg3_isempty = " << bg3_isempty << endl;
    }
    return_code = -1;
  }
  const long bg3_numpts = bg3.numPts();
  if (bg3_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 3 FAILED!!!" << endl;
      pout() << indent2 << "bg3_numpts = " << bg3_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bg4(b1);
  bg4.grow(-4*BASISV(0));

  const bool bg4_isempty = bg4.isEmpty();
  if (!bg4_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 4 FAILED!!!" << endl;
      pout() << indent2 << "bg4_isempty = " << bg4_isempty << endl;
    }
    return_code = -1;
  }
  const long bg4_numpts = bg4.numPts();
  if (bg4_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 4 FAILED!!!" << endl;
      pout() << indent2 << "bg4_numpts = " << bg4_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bg5(b1);
  bg5.growLo(0,-8);

  const bool bg5_isempty = bg5.isEmpty();
  if (!bg5_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 5 FAILED!!!" << endl;
      pout() << indent2 << "bg5_isempty = " << bg5_isempty << endl;
    }
    return_code = -1;
  }
  const long bg5_numpts = bg5.numPts();
  if (bg5_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 5 FAILED!!!" << endl;
      pout() << indent2 << "bg5_numpts = " << bg5_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bg6(b1);
  bg6.growHi(0,-8);

  const bool bg6_isempty = bg6.isEmpty();
  if (!bg6_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 6 FAILED!!!" << endl;
      pout() << indent2 << "bg6_isempty = " << bg6_isempty << endl;
    }
    return_code = -1;
  }
  const long bg6_numpts = bg6.numPts();
  if (bg6_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 6 FAILED!!!" << endl;
      pout() << indent2 << "bg6_numpts = " << bg6_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box bg7 = grow(b_empty,1);

  const bool bg7_isempty = bg7.isEmpty();
  if (!bg7_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 7 FAILED!!!" << endl;
      pout() << indent2 << "bg7_isempty = " << bg7_isempty << endl;
    }
    return_code = -1;
  }
  const long bg7_numpts = bg7.numPts();
  if (bg7_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 7 FAILED!!!" << endl;
      pout() << indent2 << "bg7_numpts = " << bg7_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bg8(b_empty);
  bg8.grow(1);

  const bool bg8_isempty = bg8.isEmpty();
  if (!bg8_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 8 FAILED!!!" << endl;
      pout() << indent2 << "bg8_isempty = " << bg8_isempty << endl;
    }
    return_code = -1;
  }
  const long bg8_numpts = bg8.numPts();
  if (bg8_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 8 FAILED!!!" << endl;
      pout() << indent2 << "bg8_numpts = " << bg8_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box bg9 = grow(b_empty, BASISV(0));

  const bool bg9_isempty = bg9.isEmpty();
  if (!bg9_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 9 FAILED!!!" << endl;
      pout() << indent2 << "bg9_isempty = " << bg9_isempty << endl;
    }
    return_code = -1;
  }
  const long bg9_numpts = bg9.numPts();
  if (bg9_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 9 FAILED!!!" << endl;
      pout() << indent2 << "bg9_numpts = " << bg9_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bg10(b_empty);
  bg10.grow(BASISV(0));

  const bool bg10_isempty = bg10.isEmpty();
  if (!bg10_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 10 FAILED!!!" << endl;
      pout() << indent2 << "bg10_isempty = " << bg10_isempty << endl;
    }
    return_code = -1;
  }
  const long bg10_numpts = bg10.numPts();
  if (bg10_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 10 FAILED!!!" << endl;
      pout() << indent2 << "bg10_numpts = " << bg10_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bg11(b_empty);
  bg11.growLo(0,1);

  const bool bg11_isempty = bg11.isEmpty();
  if (!bg11_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 11 FAILED!!!" << endl;
      pout() << indent2 << "bg11_isempty = " << bg11_isempty << endl;
    }
    return_code = -1;
  }
  const long bg11_numpts = bg11.numPts();
  if (bg11_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 11 FAILED!!!" << endl;
      pout() << indent2 << "bg11_numpts = " << bg11_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bg12(b_empty);
  bg12.growHi(0,1);

  const bool bg12_isempty = bg12.isEmpty();
  if (!bg12_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "grow isEmpty test 12 FAILED!!!" << endl;
      pout() << indent2 << "bg12_isempty = " << bg12_isempty << endl;
    }
    return_code = -1;
  }
  const long bg12_numpts = bg12.numPts();
  if (bg12_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "grow numPts test 12 FAILED!!!" << endl;
      pout() << indent2 << "bg12_numpts = " << bg12_numpts << endl;
    }
    return_code = -1;
  }
//
// coarsen tests
  const Box bc1 = coarsen(b_empty,2);

  const bool bc1_isempty = bc1.isEmpty();
  if (!bc1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "coarsen isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "bc1_isempty = " << bc1_isempty << endl;
    }
    return_code = -1;
  }
  const long bc1_numpts = bc1.numPts();
  if (bc1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "coarsen numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "bc1_numpts = " << bc1_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box bc2 = coarsen(b_empty,2*BASISV(0));

  const bool bc2_isempty = bc2.isEmpty();
  if (!bc2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "coarsen isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "bc2_isempty = " << bc2_isempty << endl;
    }
    return_code = -1;
  }
  const long bc2_numpts = bc2.numPts();
  if (bc2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "coarsen numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "bc2_numpts = " << bc2_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bc3(b_empty);
  bc3.coarsen(2);

  const bool bc3_isempty = bc3.isEmpty();
  if (!bc3_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "coarsen isEmpty test 3 FAILED!!!" << endl;
      pout() << indent2 << "bc3_isempty = " << bc3_isempty << endl;
    }
    return_code = -1;
  }
  const long bc3_numpts = bc3.numPts();
  if (bc3_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "coarsen numPts test 3 FAILED!!!" << endl;
      pout() << indent2 << "bc3_numpts = " << bc3_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bc4(b_empty);
  bc4.coarsen(2*BASISV(0));

  const bool bc4_isempty = bc4.isEmpty();
  if (!bc4_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "coarsen isEmpty test 4 FAILED!!!" << endl;
      pout() << indent2 << "bc4_isempty = " << bc4_isempty << endl;
    }
    return_code = -1;
  }
  const long bc4_numpts = bc4.numPts();
  if (bc4_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "coarsen numPts test 4 FAILED!!!" << endl;
      pout() << indent2 << "bc4_numpts = " << bc4_numpts << endl;
    }
    return_code = -1;
  }
//
// refine tests
  const Box br1 = refine(b_empty, 2);

  const bool br1_isempty = br1.isEmpty();
  if (!br1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "refine isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "br1_isempty = " << br1_isempty << endl;
    }
    return_code = -1;
  }
  const long br1_numpts = br1.numPts();
  if (br1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "refine numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "br1_numpts = " << br1_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box br2 = refine(b_empty, 2*BASISV(0));

  const bool br2_isempty = br2.isEmpty();
  if (!br2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "refine isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "br2_isempty = " << br2_isempty << endl;
    }
    return_code = -1;
  }
  const long br2_numpts = br2.numPts();
  if (br2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "refine numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "br2_numpts = " << br2_numpts << endl;
    }
    return_code = -1;
  }
//
  Box br3(b_empty);
  br3.refine(2);

  const bool br3_isempty = br3.isEmpty();
  if (!br3_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "refine isEmpty test 3 FAILED!!!" << endl;
      pout() << indent2 << "br3_isempty = " << br3_isempty << endl;
    }
    return_code = -1;
  }
  const long br3_numpts = br3.numPts();
  if (br3_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "refine numPts test 3 FAILED!!!" << endl;
      pout() << indent2 << "br3_numpts = " << br3_numpts << endl;
    }
    return_code = -1;
  }
//
  Box br4(b_empty);
  br4.refine(2*BASISV(0));

  const bool br4_isempty = br4.isEmpty();
  if (!br4_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "refine isEmpty test 4 FAILED!!!" << endl;
      pout() << indent2 << "br4_isempty = " << br4_isempty << endl;
    }
    return_code = -1;
  }
  const long br4_numpts = br4.numPts();
  if (br4_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "refine numPts test 4 FAILED!!!" << endl;
      pout() << indent2 << "br4_numpts = " << br4_numpts << endl;
    }
    return_code = -1;
  }
//
// shift tests
  Box bsh1(b_empty);
  bsh1.shift(0,1);
  const bool bsh1_isempty = bsh1.isEmpty();
  if (!bsh1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "shift isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "bsh1_isempty = " << bsh1_isempty << endl;
    }
    return_code = -1;
  }
  const long bsh1_numpts = bsh1.numPts();
  if (bsh1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "shift numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "bsh1_numpts = " << bsh1_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bsh2(b_empty);
  bsh2.shift(BASISV(0));
  const bool bsh2_isempty = bsh2.isEmpty();
  if (!bsh2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "shift isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "bsh2_isempty = " << bsh2_isempty << endl;
    }
    return_code = -1;
  }
  const long bsh2_numpts = bsh2.numPts();
  if (bsh2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "shift numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "bsh2_numpts = " << bsh2_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bsh3(b_empty);
  bsh3.shiftHalf(0,1);
  const bool bsh3_isempty = bsh3.isEmpty();
  if (!bsh3_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "shift isEmpty test 3 FAILED!!!" << endl;
      pout() << indent2 << "bsh3_isempty = " << bsh3_isempty << endl;
    }
    return_code = -1;
  }
  const long bsh3_numpts = bsh3.numPts();
  if (bsh3_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "shift numPts test 3 FAILED!!!" << endl;
      pout() << indent2 << "bsh3_numpts = " << bsh3_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bsh4(b_empty);
  bsh4.shiftHalf(BASISV(0));
  const bool bsh4_isempty = bsh4.isEmpty();
  if (!bsh4_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "shift isEmpty test 4 FAILED!!!" << endl;
      pout() << indent2 << "bsh4_isempty = " << bsh4_isempty << endl;
    }
    return_code = -1;
  }
  const long bsh4_numpts = bsh4.numPts();
  if (bsh4_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "shift numPts test 4 FAILED!!!" << endl;
      pout() << indent2 << "bsh4_numpts = " << bsh4_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box bsh5 = b_empty + BASISV(0);
  const bool bsh5_isempty = bsh5.isEmpty();
  if (!bsh5_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "shift isEmpty test 5 FAILED!!!" << endl;
      pout() << indent2 << "bsh5_isempty = " << bsh5_isempty << endl;
    }
    return_code = -1;
  }
  const long bsh5_numpts = bsh5.numPts();
  if (bsh5_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "shift numPts test 5 FAILED!!!" << endl;
      pout() << indent2 << "bsh5_numpts = " << bsh5_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box bsh6 = b_empty - BASISV(0);
  const bool bsh6_isempty = bsh6.isEmpty();
  if (!bsh6_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "shift isEmpty test 6 FAILED!!!" << endl;
      pout() << indent2 << "bsh6_isempty = " << bsh6_isempty << endl;
    }
    return_code = -1;
  }
  const long bsh6_numpts = bsh6.numPts();
  if (bsh6_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "shift numPts test 6 FAILED!!!" << endl;
      pout() << indent2 << "bsh6_numpts = " << bsh6_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bsh7(b_empty);
  bsh7 += BASISV(0);
  const bool bsh7_isempty = bsh7.isEmpty();
  if (!bsh7_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "shift isEmpty test 7 FAILED!!!" << endl;
      pout() << indent2 << "bsh7_isempty = " << bsh7_isempty << endl;
    }
    return_code = -1;
  }
  const long bsh7_numpts = bsh7.numPts();
  if (bsh7_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "shift numPts test 7 FAILED!!!" << endl;
      pout() << indent2 << "bsh7_numpts = " << bsh7_numpts << endl;
    }
    return_code = -1;
  }
//
  Box bsh8(b_empty);
  bsh8 -= BASISV(0);
  const bool bsh8_isempty = bsh8.isEmpty();
  if (!bsh8_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "shift isEmpty test 8 FAILED!!!" << endl;
      pout() << indent2 << "bsh8_isempty = " << bsh8_isempty << endl;
    }
    return_code = -1;
  }
  const long bsh8_numpts = bsh8.numPts();
  if (bsh8_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "shift numPts test 8 FAILED!!!" << endl;
      pout() << indent2 << "bsh8_numpts = " << bsh8_numpts << endl;
    }
    return_code = -1;
  }
//
// bdry tests
  const Box bb1 = bdryLo(b_empty,0,1);
  const bool bb1_isempty = bb1.isEmpty();
  if (!bb1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "bdry isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "bb1_isempty = " << bb1_isempty << endl;
    }
    return_code = -1;
  }
  const long bb1_numpts = bb1.numPts();
  if (bb1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "bdry numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "bb1_numpts = " << bb1_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box bb2 = bdryHi(b_empty,0,1);
  const bool bb2_isempty = bb2.isEmpty();
  if (!bb2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "bdry isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "bb2_isempty = " << bb2_isempty << endl;
    }
    return_code = -1;
  }
  const long bb2_numpts = bb2.numPts();
  if (bb2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "bdry numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "bb2_numpts = " << bb2_numpts << endl;
    }
    return_code = -1;
  }
//
// adjCell tests
  const Box bac1 = adjCellLo(b_empty,0,1);
  const bool bac1_isempty = bac1.isEmpty();
  if (!bac1_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "adjCell isEmpty test 1 FAILED!!!" << endl;
      pout() << indent2 << "bac1_isempty = " << bac1_isempty << endl;
    }
    return_code = -1;
  }
  const long bac1_numpts = bac1.numPts();
  if (bac1_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "adjCell numPts test 1 FAILED!!!" << endl;
      pout() << indent2 << "bac1_numpts = " << bac1_numpts << endl;
    }
    return_code = -1;
  }
//
  const Box bac2 = adjCellHi(b_empty,0,1);
  const bool bac2_isempty = bac2.isEmpty();
  if (!bac2_isempty)
  {
    if (verbose)
    {
      pout() << indent2 << "adjCell isEmpty test 2 FAILED!!!" << endl;
      pout() << indent2 << "bac2_isempty = " << bac2_isempty << endl;
    }
    return_code = -1;
  }
  const long bac2_numpts = bac2.numPts();
  if (bac2_numpts != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "adjCell numPts test 2 FAILED!!!" << endl;
      pout() << indent2 << "bac2_numpts = " << bac2_numpts << endl;
    }
    return_code = -1;
  }
//
// minbox tests
  const Box bm1 = minBox(b_empty, b1);
  const bool bm1_eq_b1 = bm1 == b1;
  if (!bm1_eq_b1)
  {
    if (verbose)
    {
      pout()  << indent2 << "minbox test 1 FAILED!!!" << endl;
    }
    return_code = -1;
  }
//
  const Box bm2 = minBox(b1,b_empty);
  const bool bm2_eq_b1 = bm2 == b1;
  if (!bm2_eq_b1)
  {
    if (verbose)
    {
      pout()  << indent2 << "minbox test 2 FAILED!!!" << endl;
    }
    return_code = -1;
  }
//
  Box bm3(b1);
  bm3.minBox(b_empty);
  const bool bm3_eq_b1 = bm3 == b1;
  if (!bm3_eq_b1)
  {
    if (verbose)
    {
      pout() << indent2 << "minbox test 3 FAILED!!!" << endl;
    }
    return_code = -1;
  }
//
  Box bm4(b_empty);
  bm4.minBox(b1);
  const bool bm4_eq_b1 = bm4 == b1;
  if (!bm4_eq_b1)
  {
    if (verbose)
    {
      pout() << indent2 << "minbox test 4 FAILED!!!" << endl;
    }
    return_code = -1;
  }
//
// lessthan tests
  const bool lt1 = b_empty < b1;
  if (lt1)
  {
    if (verbose)
    {
      pout() << indent2 << "lessthan test 1 FAILED!!!" << endl;
    }
    return_code = -1;
  }
  const bool lt2 = b1 < b_empty;
  if (!lt2)
  {
    if (verbose)
    {
      pout() << indent2 << "lessthan test 2 FAILED!!!" << endl;
    }
    return_code = -1;
  }
//
// size tests
  const IntVect size_empty = b_empty.size();
  if (size_empty != IntVect::Zero)
  {
    if (verbose)
    {
      pout() << indent2 << "empty size test 1 FAILED!!!" << endl;
      pout() << indent2 << "size_empty = " << size_empty << endl;
    }
    return_code = -1;
  }
//
  const int size_empty_0 = b_empty.size(0);
  if (size_empty_0 != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "empty size test 2 FAILED!!!" << endl;
      pout() << indent2 << "size_empty_0 = " << size_empty_0 << endl;
    }
    return_code = -1;
  }
//
  const IntVect size_i1 = bi1.size();
  if (size_i1 != IntVect::Zero)
  {
    if (verbose)
    {
      pout() << indent2 << "empty size test 3 FAILED!!!" << endl;
      pout() << indent2 << "size_i1 = " << size_i1 << endl;
    }
    return_code = -1;
  }
//
  const int size_i1_0 = bi1.size(0);
  if (size_i1_0 != 0)
  {
    if (verbose)
    {
      pout() << indent2 << "empty size test 4 FAILED!!!" << endl;
      pout() << indent2 << "size_i1_0 = " << size_i1_0 << endl;
    }
    return_code = -1;
  }
//
  if (return_code == 0)
  {
    pout() << indent << "box tests passed." << endl;
  }
  else
  {
    pout() << indent << "box tests FAILED!!!" << endl;
  }
  breakpointHook();

  return return_code;
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
