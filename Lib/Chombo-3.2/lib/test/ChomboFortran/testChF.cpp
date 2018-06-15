#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// This program tests ChomboFortran.

// Note: the standard_io tests are not done by default because
//       they can't automatically determine if the result is
//       correct and they might fail by core-dumping.
// Note: the MayDay tests are not done by default because they
//       are _intended_ to core dump if they succeed.

// Usage: a.out [-mqsv] [-E num] [args...]*
//
// Where:
//   -m selects the mayday tests
//   -q selects `quiet' execution: only pass/fail messages are printed
//   -s selects the standard_io tests
//   -v selects `verbose' execution: all messages are printed
//   -E specifies the threshold value (aka epsilon) to use for
//      comparing floating point numbers (should be <1)
//   args... are ignored.
//
// Defaults:  -q
//
// Returns:
//  The exit code of the program is the result code from the first
//  test routine that fails.

#include <cmath> //for pow()
#include <cstring>
#include <iostream>
using std::endl;
#include <memory>
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "parstream.H"
#include "SPMD.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "BoxIterator.H"
#include "CH_Complex.H"
#include "FORT_PROTO.H"
//[NOTE: this is a hand-generated header for a .F file, not a .ChF. <dbs>]
#include "test_fm2.H"

#include "UsingNamespace.H"

#include "ChFSubs_F.H"
#include "ChFIOF_F.H"
#include "test_fm_F.H"

/// Prototypes:
void
parseTestOptions(bool *option_s ,bool *option_m ,int argc ,char* argv[]) ;

/// Global variables for handling output:
static const char* pgmname = "testChF" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;
#ifdef CH_USE_DOUBLE
static Real Epsilon = 1.0e-15 ;
#else
static Real Epsilon = 1.0e-6 ;
#endif

int
main(int argc ,char* argv[])
{
  int status = 0 ;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  {//scoping trick to make sure everything is destructed before MPI_FINALIZE is called

    bool test_stdio = false ; //dont do this test by default
    bool test_mayday = false ; //dont do this by default either

    int proc = procID() ;

    parseTestOptions( &test_stdio ,&test_mayday ,argc ,argv ) ;

    if ( verbose )
      pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

    { // Test File input/output in Chombo Fortran
      int i1 = 23 ;
      Real r1 = 23.45 ;
      int status2;

      FORT_CHF_FIO( CHF_CONST_INT( i1 ) ,
                    CHF_CONST_REAL( r1 ),
                    CHF_CONST_INT( proc ),
                    CHF_INT(status2) ) ;

      if ( status2 == 0 )
      {
        if ( verbose )
          pout() << indent << pgmname << " passed the file I/O tests." << endl ;
      }
      else
      {
        pout() << indent << pgmname << " failed the file I/O tests with return code " << status2 << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }

    { // Test standard input/output in ChF
      if ( test_stdio )
      {
        int i1 = 23 ;
        Real r1 = 23.45 ;
        pout() << indent << pgmname << ": running std I/O tests.  You should see some output" << endl ;
        int status2;
        FORT_CHF_STDIO( CHF_CONST_INT( i1 ) ,CHF_CONST_REAL( r1 ), CHF_INT(status2) ) ;
        if ( status2 == 0 )
        {
          if ( verbose )
            pout() << indent << pgmname << " passed the std I/O tests." << endl ;
        }
        else
        {
          pout() << indent << pgmname << " failed the std I/O tests with return code " << status << endl ;
          if ( status == 0 ) status = status2 ;
        }
      }
    }

    { // Test ChomboFortran argument list macros for Vector<>
      // Test conditions: (from ChFSubs.ChF)
      //  vectorI length 1, value 12345678, modified to 87654321
      //  vectorCI length 2, values 23456789,98765432, not modified
      //  vectorR length 3, values 1.1,22.22,333.333, modified to 1/<value>
      //  vectorCR length 4, values 1+1/2^20,1/<value[0]>,1+1/2^40,1/<value[2]>
      //  NOTE: vectorCR[2:3] will be 1 when PRECISION=FLOAT

      int status2 = 0 ;

      const int I1 = 12345678, I1r = 87654321, I20 = 23456789, I21 = 98765432 ;
      const Real R10 = 1.1, R11 = 22.22, R12 = 333.333, R10r = 1.0/R10, R11r = 1.0/R11, R12r = 1.0/R12 ;
      const Real TwoTo20 = 1<<20 ;
      const Real R20 = 1.0+(1.0/TwoTo20), R21 = 1.0/R20 ;
      const Real R22 = 1.0+(1.0/(TwoTo20*TwoTo20)), R23 = 1.0/R22 ;
      Complex C1(R20,R23);
      Vector<int>  vectorI(1)  ; vectorI[0]  = I1  ;
      Vector<int>  vectorCI(2) ; vectorCI[0] = I20 ; vectorCI[1] = I21 ;
      Vector<Real> vectorR(3)  ; vectorR[0]  = R10 ; vectorR[1]  = R11 ; vectorR[2] = R12 ;
      Vector<Real> vectorCR(4) ; vectorCR[0] = R20 ; vectorCR[1] = R21 ; vectorCR[2] = R22 ; vectorCR[3] = R23 ;

      Vector<Complex> vectorCplx(2, C1);

      // test chfptr class
      chfptr p( vectorI.size()-1 );
      if ( p.i != 0 )
      {
       if ( verbose )
         pout() << indent2 << "error: chfptr.i returned " << p.i << " instead of 0" << endl ;
        ++status2 ;
      }
      if ( (const int *)p != &p.i )
      {
        if ( verbose )
          pout() << indent2 << "error: (const int * const)chfptr returned " << (const int *)p
               << "instead of " << &p.i << endl ;
        ++status2 ;
      }

      FORT_TESTVIRC_1( CHF_VI( vectorI ), CHF_CONST_VI( vectorCI ),
                       CHF_VR( vectorR ), CHF_CONST_VR( vectorCR ), CHF_VC( vectorCplx ),
                      CHF_CONST_REAL( Epsilon ),CHF_INT( status2 ) );
      // check return values
      if ( vectorI[0] != I1r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned int arg 1 = " << vectorI[0]
               << " instead of " << I1r << endl;
        ++status2 ;
      }
      if ( vectorCI[0] != I20 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned const int arg 2 [0] = " << vectorCI[0]
               << " instead of " << I20 << endl ;
        ++status2 ;
      }
      if ( vectorCI[1] != I21 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned const int arg 2 [1] = " << vectorCI[1]
               << " instead of " << I21 << endl ;
        ++status2 ;
      }
      if ( vectorR[0] != R10r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned Real arg 3 [0] = " << vectorR[0]
               << " instead of " << R10r << endl;
        ++status2 ;
      }
      if ( vectorR[1] != R11r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned Real arg 3 [1] = " << vectorR[1]
               << " instead of " << R11r << endl;
        ++status2 ;
      }
      if ( vectorR[2] != R12r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned Real arg 3 [2] = " << vectorR[2]
               << " instead of " << R12r << endl;
        ++status2 ;
      }
      if ( vectorCR[0] != R20 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned const Real arg 4 [0] = " << vectorCR[0]
               << " instead of " << R20 << endl;
        ++status2 ;
      }
      if ( vectorCR[1] != R21 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned const Real arg 4 [1] = " << vectorCR[1]
               << " instead of " << R21 << endl;
        ++status2 ;
      }
      if ( vectorCR[2] != R22 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned const Real arg 4 [2] = " << vectorCR[2]
               << " instead of " << R22 << endl;
        ++status2 ;
      }
      if ( vectorCR[3] != R23 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTVIR_1 returned const Real arg 4 [3] = " << vectorCR[3]
               << " instead of " << R23 << endl;
        ++status2 ;
      }
      // done with tests
      if ( status2 == 0 )
      {
        if ( verbose )
          pout() << indent << pgmname << " passed the ChF Vector<> tests." << endl ;
      }
      else
      {
        if ( verbose )
          pout() << indent << pgmname << " failed " << status2 << " of the ChF Vector<> tests." << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }

    { // Test ChomboFortran 1D array macros
      // Test conditions: (from ChFSubs.ChF)
      //  arrayI length 1, value 12345678, modified to 87654321
      //  arrayCI length 2, values 23456789,98765432, not modified
      //  arrayR length 3, values 1.1,22.22,333.333, modified to 1/<value>
      //  arrayCR length 4, values 1+1/2^20,1/<value[0]>,1+1/2^40,1/<value[2]>
      //  NOTE: arrayCR[2:3] will be 1 when PRECISION=FLOAT

      int status2 = 0 ;

      const int I1 = 12345678, I1r = 87654321, I20 = 23456789, I21 = 98765432 ;
      const Real R10 = 1.1, R11 = 22.22, R12 = 333.333, R10r = 1.0/R10, R11r = 1.0/R11, R12r = 1.0/R12 ;
      const Real TwoTo20 = 1<<20 ;
      const Real R20 = 1.0+(1.0/TwoTo20), R21 = 1.0/R20 ;
      const Real R22 = 1.0+(1.0/(TwoTo20*TwoTo20)), R23 = 1.0/R22 ;

      int arrayI[1]  ; arrayI[0]  = I1  ;
      const int arrayCI[2] =
      {
        I20,I21
      };
      Real arrayR[3]  ; arrayR[0]  = R10 ; arrayR[1]  = R11 ; arrayR[2]  = R12 ;
      const Real arrayCR[4] =
      {
        R20 ,R21 ,R22 ,R23
      };

      FORT_TESTIR1D_1( CHF_I1D( arrayI,1 ), CHF_CONST_I1D( arrayCI,2 ),
                       CHF_R1D( arrayR,3 ), CHF_CONST_R1D( arrayCR,4 ),
                       CHF_CONST_REAL( Epsilon ),CHF_INT( status2 ) );
      // check return values
      if ( arrayI[0] != I1r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned int arg 1 = " << arrayI[0]
               << " instead of " << I1r << endl;
        ++status2 ;
      }
      if ( arrayCI[0] != I20 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned const int arg 2 [0] = " << arrayCI[0]
               << " instead of " << I20 << endl ;
        ++status2 ;
      }
      if ( arrayCI[1] != I21 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned const int arg 2 [1] = " << arrayCI[1]
               << " instead of " << I21 << endl ;
        ++status2 ;
      }
      if ( arrayR[0] != R10r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned Real arg 3 [0] = " << arrayR[0]
               << " instead of " << R10r << endl;
        ++status2 ;
      }
      if ( arrayR[1] != R11r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned Real arg 3 [1] = " << arrayR[1]
               << " instead of " << R11r << endl;
        ++status2 ;
      }
      if ( arrayR[2] != R12r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned Real arg 3 [2] = " << arrayR[2]
               << " instead of " << R12r << endl;
        ++status2 ;
      }
      if ( arrayCR[0] != R20 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned const Real arg 4 [0] = " << arrayCR[0]
               << " instead of " << R20 << endl;
        ++status2 ;
      }
      if ( arrayCR[1] != R21 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned const Real arg 4 [1] = " << arrayCR[1]
               << " instead of " << R21 << endl;
        ++status2 ;
      }
      if ( arrayCR[2] != R22 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned const Real arg 4 [2] = " << arrayCR[2]
               << " instead of " << R22 << endl;
        ++status2 ;
      }
      if ( arrayCR[3] != R23 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTIR1D_1 returned const Real arg 4 [3] = " << arrayCR[3]
               << " instead of " << R23 << endl;
        ++status2 ;
      }
      // done with tests
      if ( status2 == 0 )
      {
        if ( verbose )
          pout() << indent << pgmname << " passed the ChF *1D tests." << endl ;
      }
      else
      {
        if ( verbose )
          pout() << indent << pgmname << " failed " << status2 << " of the ChF *1D tests." << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }

    { // Test calling ChomboFortran from ChomboFortran, using the 1D array macros
      //
      //[NOTE: all the code here is the same as above for TESTIR1D except the
      //       subr name and the error msgs. <dbs>]
      //
      // Test conditions: (from ChFSubs.ChF)
      //  arrayI length 1, value 12345678, modified to 87654321
      //  arrayCI length 2, values 23456789,98765432, not modified
      //  arrayR length 3, values 1.1,22.22,333.333, modified to 1/<value>
      //  arrayCR length 4, values 1+1/2^20,1/<value[0]>,1+1/2^40,1/<value[2]>
      //  NOTE: arrayCR[2:3] will be 1 when PRECISION=FLOAT

//      const int IVerbose = verbose ? 1 : 0 ;
      int status2 = 0 ;

      const int I1 = 12345678, I1r = 87654321, I20 = 23456789, I21 = 98765432 ;
      const Real R10 = 1.1, R11 = 22.22, R12 = 333.333, R10r = 1.0/R10, R11r = 1.0/R11, R12r = 1.0/R12 ;
      const Real TwoTo20 = 1<<20 ;
      const Real R20 = 1.0+(1.0/TwoTo20), R21 = 1.0/R20 ;
      const Real R22 = 1.0+(1.0/(TwoTo20*TwoTo20)), R23 = 1.0/R22 ;

      int arrayI[2]  ; arrayI[0]  = I1  ;
      const int arrayCI[2] =
      {
        I20,I21
      };
      Real arrayR[3]  ; arrayR[0]  = R10 ; arrayR[1]  = R11 ; arrayR[2]  = R12 ;
      const Real arrayCR[4] =
      {
        R20 ,R21 ,R22 ,R23
      };

      FORT_TESTCALL( CHF_I1D( arrayI,2 ), CHF_CONST_I1D( arrayCI,2 ),
                     CHF_R1D( arrayR,3 ), CHF_CONST_R1D( arrayCR,4 ),
                     CHF_CONST_REAL( Epsilon ),CHF_INT( status2 ) );
      // check return values
      if ( arrayI[0] != I1r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned int arg 1 = " << arrayI[0]
               << " instead of " << I1r << endl;
        ++status2 ;
      }
      if ( arrayCI[0] != I20 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned const int arg 2 [0] = " << arrayCI[0]
               << " instead of " << I20 << endl ;
        ++status2 ;
      }
      if ( arrayCI[1] != I21 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned const int arg 2 [1] = " << arrayCI[1]
               << " instead of " << I21 << endl ;
        ++status2 ;
      }
      if ( arrayR[0] != R10r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned Real arg 3 [0] = " << arrayR[0]
               << " instead of " << R10r << endl;
        ++status2 ;
      }
      if ( arrayR[1] != R11r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned Real arg 3 [1] = " << arrayR[1]
               << " instead of " << R11r << endl;
        ++status2 ;
      }
      if ( arrayR[2] != R12r )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned Real arg 3 [2] = " << arrayR[2]
               << " instead of " << R12r << endl;
        ++status2 ;
      }
      if ( arrayCR[0] != R20 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned const Real arg 4 [0] = " << arrayCR[0]
               << " instead of " << R20 << endl;
        ++status2 ;
      }
      if ( arrayCR[1] != R21 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned const Real arg 4 [1] = " << arrayCR[1]
               << " instead of " << R21 << endl;
        ++status2 ;
      }
      if ( arrayCR[2] != R22 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned const Real arg 4 [2] = " << arrayCR[2]
               << " instead of " << R22 << endl;
        ++status2 ;
      }
      if ( arrayCR[3] != R23 )
      {
        if ( verbose )
          pout() << indent2 << "error: TESTCALL returned const Real arg 4 [3] = " << arrayCR[3]
               << " instead of " << R23 << endl;
        ++status2 ;
      }
      // done with tests
      if ( status2 == 0 )
      {
        if ( verbose )
          pout() << indent << pgmname << " passed the ChF calls ChF tests." << endl ;
      }
      else
      {
        if ( verbose )
          pout() << indent << pgmname << " failed " << status2 << " of the ChF calls ChF tests." << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }

    { // Test ChomboFortran CHF_DSELECT macro
      const int IVerbose = verbose ? 1 : 0 ;
      int status2 = 0 ;
      FORT_TESTDSEL_1( CHF_CONST_INT( SpaceDim ),
                       CHF_CONST_INT( IVerbose ), CHF_INT( status2 ) );
      // done with tests
      if ( status2 == 0 )
      {
        if ( verbose )
          pout() << indent << pgmname << " passed the CHF_DSELECT tests." << endl ;
      }
      else
      {
        if ( verbose )
          pout() << indent << pgmname << " failed " << status2 << " of the CHF_DSELECT tests." << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }

    { // Test ChomboFortran CHF_MULTIDO macro
      Box          ibox(IntVect::Zero,IntVect::Unit*3) ;
      BaseFab<int> ifab(ibox,1) ;
      const int ivol = (int)ibox.volume() ;
      const int IVerbose = verbose ? 1 : 0 ;
      int status2 = 0 ;
      ifab.setVal(1) ;
      FORT_TESTMDO_1( CHF_FIA1( ifab,0 ) ,CHF_BOX( ibox ), CHF_CONST_INT(ivol),
                      CHF_CONST_INT( IVerbose ), CHF_INT( status2 ) );
      if ( status2 == 0 )
      {
        if ( verbose )
          pout() << indent << pgmname << " passed the CHF_MULTIDO tests." << endl ;
      }
      else
      {
        if ( verbose )
          pout() << indent << pgmname << " failed " << status2 << " of the CHF_MULTIDO tests." << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }



    { // Test ChomboFortran CHF_AUTOMULTIDO macro
      Box          ibox(IntVect::Zero,IntVect::Unit*3) ;
      // need offset box to test offset indexing
      Box offsetBox = ibox;
      offsetBox.shift(IntVect::Unit);
      BaseFab<int> ifab(ibox,1) ;
      const int ivol = (int)ibox.volume() ;
      const int IVerbose = verbose ? 1 : 0 ;
      int status2 = 0 ;
      ifab.setVal(1) ;
      FORT_TESTAMDO_1(CHF_FIA1( ifab,0 ) ,CHF_BOX( ibox ),
                      CHF_CONST_INT(ivol), CHF_BOX( offsetBox),
                      CHF_CONST_INT( IVerbose ), CHF_INT( status2 ) );
      if ( status2 == 0 )
      {
        if ( verbose )
          pout() << indent << pgmname << " passed the CHF_AUTOMULTIDO tests." << endl ;
      }
      else
      {
        if ( verbose )
          pout() << indent << pgmname << " failed " << status2 << " of the CHF_AUTOMULTIDO tests." << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }


    { // Test CHF_ID array
      int status2 = 0 ;
      FORT_TEST_CHFID(CHF_INT(status2));

      if ( status2 == 0 )
      {
        if ( verbose )
          pout() << indent << pgmname << " passed the CHF_ID tests." << endl ;
      }
      else
      {
        if ( verbose )
          pout() << indent << pgmname << " failed " << status2 << " of the CHF_ID tests." << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }

    { // Test calling MayDay routines (testing only one is enough).
      //[NOTE: this causes the program to fail, so it's not done by default (use -m).]
      if ( test_mayday )
      {
        // this shouldn't return, but if it does, something went wrong
        pout() << indent << pgmname << ": running MayDay tests. An Error or Abort should happen." << endl ;
        int status2;
        FORT_CHF_MAYDAY(CHF_INT(status2) ) ;
        pout() << indent << pgmname << " failed the MayDay tests with return code " << status << endl ;
        if ( status == 0 ) status = status2 ;
      }
    }

    { // Regression test for parsing end-of-line comments.
      // If this compiles ok, then probably it works.
      int status2 ;
      //[NOTE: this aborts on MacOS with gcc v4.0.  Probably a compiler bug. <dbs>]
      FORT_TEST_PARSE1( CHF_INT( status2 ) ) ;
      if ( verbose )
      {
        pout() << indent << pgmname ;
        if ( status2 == 0 ) pout() << " passed " ;
        else               pout() << " failed " << status2 << " of " ;
        pout() << "the parse tests." << endl ;
      }
      if ( status == 0 ) status = status2 ;
    }

    { // Regression test for statement label bug found by Francesco Minati in gcc v3.x
      // If this compiles ok, then it works.
      int i ,status2 = 0 ;
      i = 1 ; FORT_TEST_FM1(CHF_INT(i)) ;
      if ( i < 100 ) ++status2 ;
      i = 2 ; FORT_TEST_FM1(CHF_INT(i)) ;
      if ( i < 100 ) ++status2 ;
      i = 3 ; FORT_TEST_FM1(CHF_INT(i)) ;
      if ( i < 100 ) ++status2 ;
      i = 4 ; FORT_TEST_FM1(CHF_INT(i)) ;
      if ( i < 100 ) ++status2 ;
      i = 5 ; FORT_TEST_FM1(CHF_INT(i)) ;
      if ( i < 100 ) ++status2 ;
      if ( verbose )
      {
        pout() << indent << pgmname ;
        if ( status2 == 0 ) pout() << " passed " ;
        else               pout() << " failed " << status2 << " of " ;
        pout() << "the FM1 tests." << endl ;
      }
      if ( status == 0 ) status = status2 ;
    }

    { // Regression test for leading [Tab] bug found by Francesco.
      int status2;
      FORT_TEST_FM2( CHF_INT(status2) ) ;
      if ( verbose )
      {
        pout() << indent << pgmname ;
        if ( status2 == 0 ) pout() << " passed " ;
        else               pout() << " failed " << status2 << " of " ;
        pout() << "the FM2 tests." << endl ;
      }
      if ( status == 0 ) status = status2 ;
    }

    { // Test shifting macros.  Define a box centered at (0,0) and loop over
      // uniform refined data, adding it to underlying coarse data at
      // calculated coarse indices.  If it works correctly, the coarse data
      // should also be uniform (try it without shifting to see what would
      // otherise happend).

      int status2 = 0;
      int status3 = 0;
      const int refRatio = 2;
      Box crBox (-IntVect::Unit, IntVect::Unit);
      Box fnBox = refine(crBox, refRatio);
      const int crCellSum = D_TERM6( refRatio,  // Expected test result
                                    *refRatio,
                                    *refRatio,
                                    *refRatio,
                                    *refRatio,
                                    *refRatio);
      const IntVect crShiftToZero = crBox.smallEnd();
      const IntVect fnShiftToZero = scale(crShiftToZero, refRatio);

      // Integer
      {
        BaseFab<int> crFIA(crBox, 1);
        BaseFab<int> fnFIA(fnBox, 1);
        fnFIA.setVal(1);

        // FIA
        crFIA.setVal(0);
        FORT_TEST_SHIFT_FIA(CHF_FIA_SHIFT(crFIA, crShiftToZero),
                            CHF_BOX_SHIFT(fnBox, fnShiftToZero),
                            CHF_CONST_FIA_SHIFT(fnFIA, fnShiftToZero),
                            CHF_CONST_INT(refRatio));
        status3 = 0;
        for (BoxIterator bit(crBox); bit.ok(); ++bit)
          {
            status3 += (crFIA(bit(), 0) != crCellSum);
          }
        if (status3)
          {
            ++status2;
            if (verbose)
              {
                pout() << indent << pgmname
                       << " failed the CHF_FIA_SHIFT test." << endl;
              }
          }

        // FIA1
        crFIA.setVal(0);
        FORT_TEST_SHIFT_FIA1(CHF_FIA1_SHIFT(crFIA, 0, crShiftToZero),
                             CHF_BOX_SHIFT(fnBox, fnShiftToZero),
                             CHF_CONST_FIA1_SHIFT(fnFIA, 0, fnShiftToZero),
                             CHF_CONST_INT(refRatio));
        status3 = 0;
        for (BoxIterator bit(crBox); bit.ok(); ++bit)
          {
            status3 += (crFIA(bit(), 0) != crCellSum);
          }
        if (status3)
          {
            ++status2;
            if (verbose)
              {
                pout() << indent << pgmname
                       << " failed the CHF_FIA1_SHIFT test." << endl;
              }
          }
      }

      // Real
      {
        FArrayBox crFRA(crBox, 1);
        FArrayBox fnFRA(fnBox, 1);
        fnFRA.setVal(1.);

        // FRA
        crFRA.setVal(0.);
        FORT_TEST_SHIFT_FRA(CHF_FRA_SHIFT(crFRA, crShiftToZero),
                            CHF_BOX_SHIFT(fnBox, fnShiftToZero),
                            CHF_CONST_FRA_SHIFT(fnFRA, fnShiftToZero),
                            CHF_CONST_INT(refRatio));
        status3 = 0;
        for (BoxIterator bit(crBox); bit.ok(); ++bit)
          {
            status3 += (crFRA(bit(), 0) != (Real)crCellSum);
          }
        if (status3)
          {
            ++status2;
            if (verbose)
              {
                pout() << indent << pgmname
                       << " failed the CHF_FRA_SHIFT test." << endl;
              }
          }

        // FRA1
        crFRA.setVal(0);
        FORT_TEST_SHIFT_FRA1(CHF_FRA1_SHIFT(crFRA, 0, crShiftToZero),
                             CHF_BOX_SHIFT(fnBox, fnShiftToZero),
                             CHF_CONST_FRA1_SHIFT(fnFRA, 0, fnShiftToZero),
                             CHF_CONST_INT(refRatio));
        status3 = 0;
        for (BoxIterator bit(crBox); bit.ok(); ++bit)
          {
            status3 += (crFRA(bit(), 0) != (Real)crCellSum);
          }
        if (status3)
          {
            ++status2;
            if (verbose)
              {
                pout() << indent << pgmname
                       << " failed the CHF_FRA1_SHIFT test." << endl;
              }
          }
      }

      if (status2 == 0)
        {
          if ( verbose )
            pout() << indent << pgmname << " passed the CHF_SHIFT tests."
                   << endl;
        }
      else
        {
          if (verbose)
            {
              pout() << indent << pgmname << " failed " << status2 << " of the "
                "CHF_SHIFT tests." << endl;
            }
          if (status == 0) status = status2;
        }
    }
  }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status ;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Also parse the `-s' option to set the first argument
///
void
parseTestOptions(bool *option_s, bool *option_m,int argc ,char* argv[] )
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
          else if ( strncmp( argv[i] ,"-s" ,3 ) == 0 )
            {
              *option_s = true ;
            }
          else if ( strncmp( argv[i] ,"-m" ,3 ) == 0 )
            {
              *option_m = true ;
            }
          else if ( strncmp( argv[i] ,"-E" ,3 ) == 0 )
            {
              Epsilon = atof( argv[++i] );
            }
        }
    }
  return ;
}
