#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Tues, Oct 5 1999

#include <cmath>

#include "Misc.H"
#include "FArrayBox.H"
#include "BoxIterator.H"
#include "DebugOut.H"

#ifdef PGF90
// uncomment this line for pgf90...
//#include "fabTestF_F.H"
#endif

#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:
void
parseTestOptions( int argc ,char* argv[] ) ;

/**
   fabTest returns:
     0: all tests passed
    -1: setVal seems to be broken
    -2: operator*= seems to be broken
    -3: operator+= seems to be broken
    -4: operator-= seems to be broken
    -5: operator/= seems to be broken
    -6: L_\infty norm is broken
    -7: L_1 norm is broken
    -8: L_2 norm is broken
    -9: dot product is broken
   -10: linearIn, linearOut, size() is broken
   -11: axby seems to be broken (added by JNJ, 1/21/2010)
 */
extern int
fabTest(void);

/// Global variables for handling output:
static const char *pgmname = "fabTest" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

#ifdef PGF90
// define main loop here.
static int zz = 0;
extern "C" void pghpf_init(int *);
int __argc_save;
char **__argv_save;
#endif

/// Code:
int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

#ifdef PGF90
  __argc_save = argc;
  __argv_save = argv;

  pghpf_init(&zz);
#endif

  bool passed = true;
  int icode = 0;

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

#ifdef PGF90
  //  FORT_F90FABTEST(CHF_INT(icode));
#endif

  icode = fabTest();
  if (icode != 0)
    {
      pout() << indent << pgmname
           << ": FArrayBox test failed with error code " << icode
           << endl;
      passed = false;
    }

  if (passed)
      pout() << indent << pgmname
           << ": FArrayBox passed all tests"
           << endl;

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}

int fabTest()
{
  const int nComps = 3;

  //this is the error code that will get returned if a test fails
  int icode = 0;

  int ihigh = 8;
  Real oldInit;
  Real newInit = 9.87654321e+10;
#ifdef CH_USE_FLOAT
  Real eps = 1.0e-6;
#else
  Real eps = 1.0e-12;
#endif
  //test 1.  make a box of a certain size.  set it to one.  see if
  //add, subtract, multiply divide
  IntVect hiVect = (ihigh-1)*IntVect::Unit;
  IntVect loVect(D_DECL6(2,3,4,2,3,4));
  Box B(loVect, hiVect);
  FArrayBox ourfabs[5];
  ourfabs[0].define(B,nComps);

  ///initialization test
#ifdef CH_USE_SETVAL
  if (((ourfabs[0].max() > BaseFabRealSetVal+eps)||
      (ourfabs[0].max() < BaseFabRealSetVal-eps)) ||
     ((ourfabs[0].min() > BaseFabRealSetVal+eps)||
      (ourfabs[0].min() < BaseFabRealSetVal-eps)))
    {
      pout() << indent << pgmname
           << ": fab::initial setval 1 might have failed"
           << endl ;
      icode = -1;
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::initial setval 1 passed"
           << endl ;
    }
#else
  if (((ourfabs[0].max() > 0.0+eps)||
      (ourfabs[0].max() < 0.0-eps)) ||
     ((ourfabs[0].min() > 0.0+eps)||
      (ourfabs[0].min() < 0.0-eps)))
    {
      pout() << indent << pgmname
           << ": fab::initial calloc 1 might have failed"
           << endl ;
      icode = -2;
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::initial calloc 1 passed"
           << endl ;
    }
#endif

  oldInit = BaseFabRealSetVal;
  BaseFabRealSetVal = newInit;
  ourfabs[1].define(B,nComps);

#ifdef CH_USE_SETVAL
  if (((ourfabs[1].max() > newInit+eps)||
      (ourfabs[1].max() < newInit-eps)) ||
     ((ourfabs[1].min() > newInit+eps)||
      (ourfabs[1].min() < newInit-eps)))
    {
      pout() << indent << pgmname
           << ": fab::initial setval 2 might have failed"
           << endl ;
      icode = -3;
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::initial setval 2 passed"
           << endl ;
    }
#else
  if (((ourfabs[1].max() > 0.0+eps)||
      (ourfabs[1].max() < 0.0-eps)) ||
     ((ourfabs[1].min() > 0.0+eps)||
      (ourfabs[1].min() < 0.0-eps)))
    {
      pout() << indent << pgmname
           << ": fab::initial calloc 2 might have failed"
           << endl ;
      icode = -4;
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::initial calloc 2 passed"
           << endl ;
    }
#endif

  BaseFabRealSetVal = oldInit;

  ///setVal test
  ourfabs[0].setVal(1.);

  if (((ourfabs[0].max() > 1.0+eps)||
      (ourfabs[0].max() < 1.0-eps)) ||
     ((ourfabs[0].min() > 1.0+eps)||
      (ourfabs[0].min() < 1.0-eps)))
    {
      pout() << indent << pgmname
           << ": fab::setval might have failed"
           << endl ;
      icode = -5;
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::setval passed"
           << endl ;
    }


  //
  // Test minIndex() and maxIndex()
  //
  IntVect theindex = B.smallEnd()+2;
  Real saveme = ourfabs[0](theindex, 0 );
  ourfabs[0].setVal( -100.0, Box(theindex, theindex), 0 );
  IntVect answer = ourfabs[0].minIndex(ourfabs[0].box(), 0);
  if ( answer != theindex )
  {
    pout() << "minIndex() might have failed (returned " << answer
         << ", should be " << theindex << ")" << endl;
    icode = -17;
  }

  ourfabs[0].setVal( +100.0, Box(theindex, theindex), 0 );
  answer = ourfabs[0].maxIndex(0);
  if ( answer != theindex )
  {
    pout() << "maxIndex() might have failed (returned " << answer
         << ", should be " << theindex << ")" << endl;
    icode = -18;
  }

  ourfabs[0].setVal( saveme, Box(theindex, theindex), 0 );


  ///arithmetic operator tests
  ourfabs[1].copy(ourfabs[0]);
  ourfabs[1] *= -2.0;
  if (((ourfabs[1].max() > -2.0+eps)||
      (ourfabs[1].max() < -2.0-eps)) ||
     ((ourfabs[1].min() > -2.0+eps)||
      (ourfabs[1].min() < -2.0-eps)))
    {
      pout() << indent << pgmname
           << ": fab::operator*= might have failed"
           << endl ;
      icode = -6;
      return(icode);
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::operator*= passed"
           << endl ;
    }

  ourfabs[2].define(B,nComps);
  ourfabs[2].copy(ourfabs[0]);
  ourfabs[2] -= ourfabs[1];

  if (((ourfabs[2].max() > 3.0+eps)||
      (ourfabs[2].max() < 3.0-eps)) ||
     ((ourfabs[2].min() > 3.0+eps)||
      (ourfabs[2].min() < 3.0-eps)))
    {
      pout() << indent << pgmname
           << ": fab::operator+= might have failed"
           << endl ;
      icode = -7;
      return(icode);
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::operator+= passed"
           << endl ;
    }

  if (((ourfabs[1].max() > -2.0+eps)||
      (ourfabs[1].max() < -2.0-eps)) ||
     ((ourfabs[1].min() > -2.0+eps)||
      (ourfabs[1].min() < -2.0-eps)))
    {
      pout() << indent << pgmname
           << ": fab::operator-= might have failed"
           << endl ;
      icode = -8;
      return(icode);
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::operator-= passed"
           << endl ;
    }

  ourfabs[3].define(B,nComps);
  ourfabs[3].setVal(8.0);
  ourfabs[3] /= ourfabs[1];
  ourfabs[3] *= -1.0;
  if (((ourfabs[3].max() > 4.0+eps)||
      (ourfabs[3].max() < 4.0-eps)) ||
     ((ourfabs[3].min() > 4.0+eps)||
      (ourfabs[3].min() < 4.0-eps)))
    {
      pout() << indent << pgmname
           << ": fab::operator/= might have failed"
           << endl ;
      icode = -9;
      return(icode);
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::operator/= passed"
           << endl ;
    }

  ///norm tests
  int normStatus = 0;
  if (Abs(ourfabs[0].norm(0,0,1) - 1.0) > eps) normStatus += 1;

  if (Abs(ourfabs[1].norm(0,0,1) - 2.0) > eps) normStatus += 10;
  if (Abs(ourfabs[2].norm(0,0,1) - 3.0) > eps) normStatus += 100;
  if (Abs(ourfabs[3].norm(0,0,1) - 4.0) > eps) normStatus += 1000;
  if (normStatus != 0)
    {
      pout() << indent << pgmname
           << ": fab::norm(0) might have failed"
           << endl ;
      icode = -10;
      return(icode);
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::norm(0) passed"
           << endl ;
    }

  Real rnorm = B.numPts();

  normStatus = 0;
  if (Abs(ourfabs[0].norm(1,0,3)/rnorm - 3*1.0) > eps) normStatus += 1;
  if (Abs(ourfabs[1].norm(1,0,3)/rnorm - 3*2.0) > eps) normStatus += 10;
  if (Abs(ourfabs[2].norm(1,1,2)/rnorm - 2*3.0) > eps) normStatus += 100;
  if (Abs(ourfabs[3].norm(1,2,1)/rnorm - 1*4.0) > eps) normStatus += 1000;
  if (normStatus != 0)
    {
      pout() << indent << pgmname
           << ": fab::norm(1) might have failed"
           << endl ;
      icode = -11;
      return(icode);
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::norm(1) passed"
           << endl ;
    }

  // norm on restricted box.
  Box restrictedBox( B );
  restrictedBox.grow( -IntVect::Unit );
  if (Abs(ourfabs[0].norm(restrictedBox,1,1,2) - 2*restrictedBox.numPts()) > eps)
  {
    pout() << indent << pgmname
         << ": fab::norm(1) on restricted box might have failed"
         << endl ;
    icode = -11;
    return(icode);
  }

  rnorm = sqrt(rnorm);
  normStatus = 0;
  if (Abs(ourfabs[0].norm(2,0,1)/rnorm - 1.0) > eps) normStatus += 1;
  if (Abs(ourfabs[1].norm(2,0,1)/rnorm - 2.0) > eps) normStatus += 10;
  if (Abs(ourfabs[2].norm(2,0,1)/rnorm - 3.0) > eps) normStatus += 100;
  if (Abs(ourfabs[3].norm(2,0,1)/rnorm - 4.0) > eps) normStatus += 1000;
  if (normStatus != 0)
    {
      pout() << indent << pgmname
           << ": fab::norm(2) might have failed"
           << endl ;
      icode = -12;
      return(icode);
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::norm(2) passed"
           << endl ;
    }


  // dot product tests
  Real numPts = B.numPts();
  if ((Abs(ourfabs[0].dotProduct(ourfabs[0])/numPts -  1.0*nComps) > eps) ||
     (Abs(ourfabs[0].dotProduct(ourfabs[1])/numPts +  2.0*nComps) > eps) ||
     (Abs(ourfabs[0].dotProduct(ourfabs[2])/numPts -  3.0*nComps) > eps) ||
     (Abs(ourfabs[0].dotProduct(ourfabs[3])/numPts -  4.0*nComps) > eps) ||
     (Abs(ourfabs[1].dotProduct(ourfabs[0])/numPts +  2.0*nComps) > eps) ||
     (Abs(ourfabs[1].dotProduct(ourfabs[1])/numPts -  4.0*nComps) > eps) ||
     (Abs(ourfabs[1].dotProduct(ourfabs[2])/numPts +  6.0*nComps) > eps) ||
     (Abs(ourfabs[1].dotProduct(ourfabs[3])/numPts +  8.0*nComps) > eps) ||
     (Abs(ourfabs[2].dotProduct(ourfabs[0])/numPts -  3.0*nComps) > eps) ||
     (Abs(ourfabs[2].dotProduct(ourfabs[1])/numPts +  6.0*nComps) > eps) ||
     (Abs(ourfabs[2].dotProduct(ourfabs[2])/numPts -  9.0*nComps) > eps) ||
     (Abs(ourfabs[2].dotProduct(ourfabs[3])/numPts - 12.0*nComps) > eps) ||
     (Abs(ourfabs[3].dotProduct(ourfabs[0])/numPts -  4.0*nComps) > eps) ||
     (Abs(ourfabs[3].dotProduct(ourfabs[1])/numPts +  8.0*nComps) > eps) ||
     (Abs(ourfabs[3].dotProduct(ourfabs[2])/numPts - 12.0*nComps) > eps) ||
     (Abs(ourfabs[3].dotProduct(ourfabs[3])/numPts - 16.0*nComps) > eps))
  {
    pout() << "ourfabs[0].dotProduct(ourfabs[0])/numPts = (should=3.0) "
         << ourfabs[0].dotProduct(ourfabs[0])/numPts << '\n';
    pout() << "ourfabs[1].dotProduct(ourfabs[1])/numPts = (should=-6.0) "
         << ourfabs[1].dotProduct(ourfabs[0])/numPts << '\n';

    pout() << indent << pgmname
         << ": fab::dotProduct might have failed"
         << endl ;
    icode = -13;
    return(icode);
  }
  else if (verbose)
  {
    pout() << indent2 << pgmname
         << ": fab::dotProduct 1 passed"
         << endl ;
  }

  IntVect hiVect2 = (ihigh-2)*IntVect::Unit;
  IntVect loVect2(D_DECL6(0,1,2,0,1,2));
  Box B2(loVect2, hiVect2);

  FArrayBox fabFive (B2,nComps);

  ourfabs[4].define(B,nComps);
  ourfabs[4].setVal((Real)  1.0,0);
  ourfabs[4].setVal((Real) 10.0,1);
  ourfabs[4].setVal((Real)100.0,2);

  fabFive .setVal((Real)  2.0,0);
  fabFive .setVal((Real) 20.0,1);
  fabFive .setVal((Real)200.0,2);

  Box B3 = B & B2;
  numPts = B3.numPts();

  if ( (Abs(ourfabs[4].dotProduct(fabFive )/numPts - 20202.0) > eps) ||
      (Abs(fabFive.dotProduct(ourfabs[4])/numPts - 20202.0) > eps))
  {
    pout() << indent << pgmname
         << ": fab::dotProduct might have failed"
         << endl ;
    icode = -14;
    return(icode);
  }
  else if (verbose)
  {
    pout() << indent2 << pgmname
         << ": fab::dotProduct 2 passed"
         << endl ;
  }

  FArrayBox iotaFab( Box(IntVect::Zero, IntVect(D_DECL6(1,2,3,1,2,3))), 2 );
  //FArrayBox iotaFab( Box(IntVect::Zero, IntVect(D_DECL6(1,0,1,1,1,1))), 2 );
  for ( int comp=0;comp<2;++comp )
  {
      for ( int i=0;i<iotaFab.box().numPts();++i )
      {
        iotaFab.dataPtr()[i + comp*iotaFab.box().numPts()] = i + comp*10;
      }
  }
  //dumpFAB( &iotaFab );


  //
  // MaskLT, maskLE, etc.
  //
#if CH_SPACEDIM > 1
  // Cop-out for until I code in the correct results for 1D

  BaseFab<int> mask;

  iotaFab.maskLT( mask, 4, 0);
  int maskCount = std::count( mask.dataPtr(), mask.dataPtr()+mask.box().numPts(), 1 );
  if ( maskCount != 4 )
  {
    if (verbose)
    {
      pout() << "maskLT failed: returned " << maskCount << ", should be 4\n";
    }
    icode = -20;
  }
  iotaFab.maskLT( mask, 14, 1);
  maskCount = std::count( mask.dataPtr(), mask.dataPtr()+mask.box().numPts(), 1 );
  if ( maskCount != 4 )
  {
    if (verbose)
    {
      pout() << "maskLT failed: returned " << maskCount << ", should be 4\n";
    }
    icode = -20;
  }

  iotaFab.maskLE( mask, 4, 0);
  maskCount = std::count( mask.dataPtr(), mask.dataPtr()+mask.box().numPts(), 1 );
  if ( maskCount != 5 )
  {
    if (verbose)
    {
      pout() << "maskLE failed: returned " << maskCount << ", should be 5\n";
    }
    icode = -20;
  }
  iotaFab.maskLE( mask, 14, 1);
  maskCount = std::count( mask.dataPtr(), mask.dataPtr()+mask.box().numPts(), 1 );
  if ( maskCount != 5 )
  {
    if (verbose)
    {
      pout() << "maskLE failed: returned " << maskCount << ", should be 5\n";
    }
    icode = -20;
  }

  iotaFab.maskEQ( mask, 4, 0);
  maskCount = std::count( mask.dataPtr(), mask.dataPtr()+mask.box().numPts(), 1 );
  if ( maskCount != 1 )
  {
    if (verbose)
    {
      pout() << "maskEQ failed: returned " << maskCount << ", should be 1\n";
    }
    icode = -20;
  }
  iotaFab.maskEQ( mask, 14, 1);
  maskCount = std::count( mask.dataPtr(), mask.dataPtr()+mask.box().numPts(), 1 );
  if ( maskCount != 1 )
  {
    if (verbose)
    {
      pout() << "maskEQ failed: returned " << maskCount << ", should be 1\n";
    }
    icode = -20;
  }
#endif // CH_SPACEDIM > 1

  //
  // Test sum()
  //
  {
    Real sumAnswer = iotaFab.sum( 1,1 );
    int numPts = iotaFab.box().numPts();
    Real correctAnswer = 10*numPts + numPts*(numPts-1)/2;
    if ( sumAnswer != correctAnswer )
    {
      if (verbose)
      {
        pout() << "sum() failed: returned " << sumAnswer << ", should be "
             << correctAnswer << '\n';
      }
      icode = -21;
    }
  }


  //
  // Test invert()
  //
  {
    FArrayBox inverseFab( iotaFab.box(), iotaFab.nComp() );
    inverseFab.copy( iotaFab );
    inverseFab += 1.0;
    inverseFab.invert( 1.0 );

    Real sumAnswer = inverseFab.sum( 1,1 );
    Real correctAnswer;
    if     ( CH_SPACEDIM == 1 ) correctAnswer = 0.174242;
    else if ( CH_SPACEDIM == 2 ) correctAnswer = 0.451761;
    else if ( CH_SPACEDIM == 3 ) correctAnswer = 1.18924;
    else if ( CH_SPACEDIM == 4 ) correctAnswer = 1.717286;
    else if ( CH_SPACEDIM == 5 ) correctAnswer = 2.68844325;
    else if ( CH_SPACEDIM == 6 ) correctAnswer = 4.0224202;
    else MayDay::Error("undefined dimension for inverse test");
    if ( Abs( sumAnswer - correctAnswer ) > 0.00001 )
    {
      if (verbose)
      {
        pout() << "inverse() failed: returned " << sumAnswer << ", correct="
             << correctAnswer << '\n';
      }
      icode = -22;
    }
  }

  //===alias testing==============================================
  {
    Interval comps(1,2);
    FArrayBox aliasFive(comps, ourfabs[4]);
    FArrayBox aliasSix;
    aliasSix.define(comps, fabFive);

    if ((aliasFive.dotProduct(aliasSix)/numPts - 20200.0 > eps) ||
        (aliasSix .dotProduct(aliasFive)/numPts - 20200.0 > eps))
      {
        pout() << indent << pgmname
             << ": fab aliasing might have failed"
             << endl ;
        icode = -15;
        return(icode);
      }
    else if (verbose)
      {
        pout() << indent2 << pgmname
             << ": fab aliasing passed"
             << endl ;
      }
  } // aliased FArrayBox's drop out of scope here. Test that original data unharmed.

  if (ourfabs[4].dotProduct(fabFive )/numPts - 20202.0 > eps)
    {
      pout() << indent << pgmname
           << ": fab aliasing might have had side-effect on original"
           << endl ;
      icode = -16;
      return(icode);
    }

  // ======linearization tests=====================================
  Box B4(loVect+2*IntVect::Unit, hiVect);
  CH_assert(!B4.isEmpty());

  BaseFab<Real> intfab1(B, 2), intfab2(B, 3), intfab3(B4, 1);

  intfab1.setVal(-11);
  intfab2.setVal(2);
  intfab3.setVal(3);

  Interval interval1(0, 1);
  Interval interval2(0, 0);
  Interval interval3(2, 2);
  int sz = intfab1.size(B4, interval1);
  void* buf = malloc(sz);
  CH_assert(buf != NULL);
  intfab1.linearOut(buf, B4, interval1);
  intfab2.linearIn(buf, B4, interval1);

  for (BoxIterator bit(B); bit.ok(); ++bit)
    { // check that non-modified data actually is non-modifed
      if (!(intfab1(bit(), 0) == -11 &&
           intfab1(bit(), 1) == -11 &&
           intfab2(bit(), 2) == 2))
        {
          pout() << indent << pgmname
               << ": BaseFab linearization test #1 might have failed."
               << endl ;
          return -10;
        }
    }
  for (BoxIterator bit(B4); bit.ok(); ++bit)
    { //check that modified data is properly modified.
      if (!(intfab2(bit(), 0) == -11 &&
           intfab2(bit(), 1) == -11))
        {
          pout() << indent << pgmname
               << ": BaseFab linearization test #2 might have failed"
               << endl ;
          return -10;
        }
    }

  int sz2 = intfab3.size(B4, interval2);
  //reuse buffer here if possible
  if (sz2 > sz)
    {
      free(buf);
      buf = malloc(sz2);
    }
  CH_assert(buf != NULL);
  intfab3.linearOut(buf, B4, interval2);
  intfab2.linearIn(buf, B4, interval3);

  for (BoxIterator bit(B4); bit.ok(); ++bit)
    {
      if (!(intfab3(bit(), 0) == 3 &&
           intfab2(bit(), 2) == 3))
        {
          pout() << indent << pgmname
               << ": BaseFab linearization test #3 might have failed"
               << endl ;
          return -10;
        }
    }

  free(buf);

  if ( verbose ) pout() << indent2 << pgmname
                     << ": fab::linearIn/Out passed"
                     << endl ;

  /// axby method tests
  ourfabs[2].axby(ourfabs[0], ourfabs[1], 1.0, 2.0);
  if (((ourfabs[2].max() > (ourfabs[0].max() + 2.0*ourfabs[1].max() + eps))||
      (ourfabs[2].max() < (ourfabs[0].max() + 2.0*ourfabs[1].max() - eps))) ||
     ((ourfabs[2].min() > (ourfabs[0].min() + 2.0*ourfabs[1].min() + eps)) ||
      (ourfabs[2].min() < (ourfabs[0].min() + 2.0*ourfabs[1].min() - eps))))
    {
      pout() << indent << pgmname
           << ": fab::axby might have failed"
           << endl ;
      icode = -11;
      return(icode);
    }
  else if (verbose)
    {
      pout() << indent2 << pgmname
           << ": fab::axby passed"
           << endl ;
    }

  //=======end linearization tests =========================================
  return(icode);
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
