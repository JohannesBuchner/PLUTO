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
#include "CH_Timer.H"
#include "BaseFabMacros.H"
#include "FArrayBox.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testBaseFabMacros();

// this is a FArrayBox-derived class so that we can test the "This"
// macros
class MacroTesterFab : public FArrayBox
{
public:
  MacroTesterFab(const Box& a_box, int a_nComp);

  virtual ~MacroTesterFab()
  {
  }

  // returns status after tests; 0 means everything passed
  int testMacros(bool& a_verbose);

};


/// Global variables for handling output:
static const char *pgmname = "testBaseFabMacros" ;
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

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int status = testBaseFabMacros() ;

  // I'm thinking all pout()'s should be before MPI_Finalize... (ndk)
  pout() << indent << pgmname << ": "
         << ( (status == 0) ? "passed all tests" : "failed at least one test,")
         << endl;


#ifdef CH_MPI
  MPI_Finalize();
#endif

  CH_STOP(t1);
  CH_TIMER_REPORT();


  return status;
}

int
testBaseFabMacros()
{
  int return_code = 0;

  int boxHi = 7;
  IntVect loVect = IntVect::Zero;
  IntVect hiVect = boxHi*IntVect::Unit;
  Box box1(loVect,hiVect);

  FArrayBox fab1(box1, 1);
  fab1.setVal(2.0);

  IntVect offset = 3*IntVect::Unit;
  loVect += offset;
  hiVect += offset;
  Box box2(loVect, hiVect);
  FArrayBox fab2(box2,1);
  fab2.setVal(1.0);

  // this code implements all of the macros in
  // BaseFabMacros.H which don't require a "this" (in other words,
  // all of the ones which can be called outside of a BaseFab-derived
  // object). I've left the "this" macros commented out here just to
  // document the ones we're _not_ testing.
  // We may want to come back at some point and expand coverage to the
  // "this" macros.
  //
  // At the moment, all we're really
  // testing for here is that they compile and run to completion
  // we should come back at some point and add substantive testing
  // as the need (or perception of need) arises.

  ForAllXBPencil(Real, fab1, box1, 0, 1)
    {
      Real* pencil = xR;
      for (int i=0; i<thisLen; i++)
        {
          pencil[i] += 1.0;
        }
    }   EndForPencil

  int testVal = 3.0;
  BoxIterator bit(box1);
  for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (fab1(iv, 0) != testVal)
        {
          if (verbose)
            {
              pout() << "  ForAllXBPencil test FAILED at " << iv << endl;
            }
          return_code += 1;
        }
    }


  ForAllXBNNnoindx(Real, fab1, box1, 0, 1)
    {
      fab1R += 1;
    } EndFor

  ++testVal;
  for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (fab1(iv, 0) != testVal)
        {
          if (verbose)
            {
              pout() << "  ForAllXBNNnoindx test FAILED at " << iv << endl;
            }
          return_code += 1;
        }
    }


  ForAllXBNN(Real, fab1, box1, 0, 1)
    {
      fab1R += 1;
    } EndFor

  ++testVal;
  for (bit.reset(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (fab1(iv, 0) != testVal)
        {
          if (verbose)
            {
              pout() << "  ForAllXBNN test FAILED at " << iv << endl;
            }
          return_code += 1;
        }
    }



  int sum = 0;
  int numPts = box1.numPts();
  // const version needs to be done a bit differently --
  // sum over values in fab
  ForAllXCBNN(Real, fab1, box1, 0, 1)
    {
      sum += fab1R;
    } EndFor;

  if (sum != numPts*testVal)
    {
      if (verbose)
        {
          pout() << "  ForAllXCBNN test FAILED!" << endl;
        }
      return_code += 1;
    }

  // now test "This" macros
  MacroTesterFab testerFab(box1, 1);

  int fabStatus = testerFab.testMacros(verbose);
  if (fabStatus != 0)
    {
      if (verbose)
        {
          pout() << "  Fab-internal macro test FAILED!" << endl;
        }
      return_code += 100;
    }

#if 0
   // (DFM 1/25/10) -- I'm not sure what this macro is doing, and we don't
   // appear to use it anywhere (it's not defined in 4D, for example),
   // so I'm deferring implementing a test for it.
   // We may want to consider removing it from BaseFabMacros.H
   ForAllRevXBNYCBNNN(Real, fab1, box1, 0, fab2, box2, 0, 1, 0)
 {

 } EndFor;
#endif

  return return_code;
}

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

MacroTesterFab::MacroTesterFab(const Box& a_box, int a_nComp)
{
  FArrayBox::define(a_box, a_nComp);
}

int
MacroTesterFab::testMacros(bool& a_verbose)
{
  int status = 0;

  const Box& bx = box();
  int numPts = bx.numPts();

  FArrayBox fab1(bx,1);
  fab1.setVal(1.0);


  FArrayBox fab2(bx,1);
  fab2.setVal(2.0);

  // do a setval for the first test...
  ForAllThisBNN(Real,bx, 0, 1)
    {
      thisR = 1.0;
    } EndFor;

  BoxIterator bit(bx);
  Real testVal = 1.0;
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (this->operator()(iv,0) != testVal)
        {
          if (a_verbose)
            {
              pout() << "  ForAllThisBNN test FAILED at " << iv << "!" << endl;
            }
          status += 1;
        }
    }


  Real thisSum = 0.0;
  ForAllThisCBNN(Real,bx, 0, 1)
    {
      thisSum += thisR;
    } EndFor;

  Real tolerance = 1.0e-8;
  thisSum -= numPts;
  if (thisSum < 0) thisSum *= -1.0;
  if ( thisSum > tolerance)
    {
      if (a_verbose)
        {
          pout() << "  ForAllThisCBNN test FAILED! " << endl;
        }
      status += 1;
    }

  ForAllThisBNNXC(Real, bx, 0, 1, fab1, 0)
    {
      thisR += fab1R;
    } EndForTX;

  testVal += 1.0;
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (this->operator()(iv,0) != testVal)
        {
          if (a_verbose)
            {
              pout() << "  ForAllThisBNNXC test FAILED at "
                     << iv << "!" << endl;
            }
          status += 1;
        }
    }

  thisSum = 0.0;
  ForAllThisCBNNXC(Real, bx, 0, 1, fab1, 0)
    {
      thisSum += (thisR + fab1R);
    } EndForTX;

  testVal += 1.0;
  thisSum -= testVal*numPts;
  if (thisSum < 0) thisSum *= -1.0;
  if ( thisSum > tolerance)
    {
      if (a_verbose)
        {
          pout() << "  ForAllThisCBNNXC test FAILED! " << endl;
        }
      status += 1;
    }



  ForAllThisBNNXCBN(Real, bx, 0, 1, fab1, bx,0)
    {
      thisR += fab1R;
    } EndForTX;

  // already incremented in the last test
  //testVal += 1.0;
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (this->operator()(iv,0) != testVal)
        {
          if (a_verbose)
            {
              pout() << "  ForAllThisBNNXCBN test FAILED at "
                     << iv << "!" << endl;
            }
          status += 1;
        }
    }


  ForAllThisBNNXCBNYCBN(Real,bx,0,1,fab1,bx,0,fab2,bx,0)
    {
      thisR += (fab1R + fab2R);
    } EndForTX;

  testVal += 3.0;
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (this->operator()(iv,0) != testVal)
        {
          if (a_verbose)
            {
              pout() << "  ForAllThisBNNXCBNYCBN test FAILED at "
                     << iv << "!" << endl;
            }
          status += 1;
        }
    }

  thisSum = 0.0;
  ForAllThisCPencil(Real, bx, 0, 1)
    {
      const Real* pencil = &thisR;
      for (int i=0; i<thisLen; i++)
        {
          thisSum += pencil[i];
        }
    }   EndForPencil

  thisSum -= testVal*numPts;
  if (thisSum < 0) thisSum *= -1.0;
  if ( thisSum > tolerance)
    {
      if (a_verbose)
        {
          pout() << "  ForAllThisCPencil test FAILED! " << endl;
        }
      status += 1;
    }

  ForAllThisPencil(Real, bx, 0, 1)
    {
      Real* pencil = &thisR;
      for (int i=0; i<thisLen; i++)
        {
          pencil[i] += 1;
        }
    }   EndForPencil


  testVal += 1.0;
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (this->operator()(iv,0) != testVal)
        {
          if (a_verbose)
            {
              pout() << "  ForAllThisPencil test FAILED at "
                     << iv << "!" << endl;
            }
          status += 1;
        }
    }


  return status;
}
