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
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "CH_Attach.H"
#include "UsingNamespace.H"

/// Global variables for test output:

static const char *pgmname = "copyTest" ;
static const char *indent = "   " ,*indent2 = "      " ;
static bool verbose = true ;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

/**
   copyTest returns:
    0: all tests passed
 */
extern int
copyTest(void);

/// Code:
int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  bool passed = true;
  int icode = 0;
  //registerDebugger();
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  icode = copyTest();
  if (icode != 0)
    {
      pout() << indent << pgmname
           << " failed with error code " << icode
           << endl;
      passed = false;
    }

  if (passed)
    {
      pout() << indent << pgmname
           << " passed all tests"
           << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}

int
copyTest()
{
  DisjointBoxLayout loLayoutFine, loLayoutCoar;
#if (CH_SPACEDIM == 1)
  // no longer a bogus spacedim
  Box b1(IntVect(7), IntVect(7));
  Box b2(IntVect(3), IntVect(3));
  Box b3(IntVect(9), IntVect(9));
  Box b4(IntVect(-1), IntVect(-1));
#elif (CH_SPACEDIM == 2)
  Box b1(IntVect(7,6), IntVect(7,9));
  Box b2(IntVect(3,0), IntVect(3,3));
  Box b3(IntVect(9,4), IntVect(9,5));
  Box b4(IntVect(-1,4), IntVect(-1,12));
#elif (CH_SPACEDIM == 3)
  Box b1(IntVect(7,6,0), IntVect(7,9,0));
  Box b2(IntVect(3,0,0), IntVect(3,3,0));
  Box b3(IntVect(9,4,0), IntVect(9,5,0));
  Box b4(IntVect(-1,4,0), IntVect(-1,12,0));
#elif (CH_SPACEDIM == 4)
  Box b1(IntVect(7,6,0,0), IntVect(7,9,0,0));
  Box b2(IntVect(3,0,0,0), IntVect(3,3,0,0));
  Box b3(IntVect(9,4,0,0), IntVect(9,5,0,0));
  Box b4(IntVect(-1,4,0,0), IntVect(-1,12,0,0));
#elif (CH_SPACEDIM == 5)
  Box b1(IntVect(7,6,0,0,0), IntVect(7,9,0,0,0));
  Box b2(IntVect(3,0,0,0,0), IntVect(3,3,0,0,0));
  Box b3(IntVect(9,4,0,0,0), IntVect(9,5,0,0,0));
  Box b4(IntVect(-1,4,0,0,0), IntVect(-1,12,0,0,0));
#elif (CH_SPACEDIM == 6)
  Box b1(IntVect(7,6,0,0,0,0), IntVect(7,9,0,0,0,0));
  Box b2(IntVect(3,0,0,0,0,0), IntVect(3,3,0,0,0,0));
  Box b3(IntVect(9,4,0,0,0,0), IntVect(9,5,0,0,0,0));
  Box b4(IntVect(-1,4,0,0,0,0), IntVect(-1,12,0,0,0,0));
#else
  bogus spacedim;
#endif
  Vector<Box> boxesCoar(3);
  Vector<int> procsCoar(3, 0);
  Vector<Box> boxesFine(4);
  Vector<int> procsFine(4, 0);
  boxesCoar[0] = b1;
  boxesCoar[1] = b2;
  boxesCoar[2] = b3;

  boxesFine[0] = b1;
  boxesFine[1] = b2;
  boxesFine[2] = b3;
  boxesFine[3] = b4;

  loLayoutCoar.define(boxesCoar, procsCoar);
  loLayoutFine.define(boxesFine, procsFine);

  LevelData<FArrayBox> loFabCoar(loLayoutCoar,1);
  LevelData<FArrayBox> loFabFine(loLayoutFine,1);

  DataIterator ditCoar =loFabCoar.dataIterator();

  Real fluxval = 4.77;
  int icode = 0;
  for (ditCoar.reset(); ditCoar.ok(); ++ditCoar)
    {
      FArrayBox& fab = loFabCoar[ditCoar];
      Real val;
      const Box& b = loLayoutCoar.get(ditCoar);
      if (b == b2)
        val = -fluxval;
      else
        val = 0.0;

      fab.setVal(val);
      val = 1.0;
    }

  DataIterator ditFine =loFabFine.dataIterator();
  for (ditFine.reset(); ditFine.ok(); ++ditFine)
    {
      FArrayBox& fab = loFabFine[ditFine()];
      Real val = -fluxval/2.0;
      fab.setVal(val);
    }
  Copier copierCoarToFine(loLayoutCoar, loLayoutFine, loFabFine.ghostVect());
  pout() << "layout = " << loLayoutCoar << endl;
  pout() << "copier = " << copierCoarToFine << endl;
  loFabCoar.copyTo( loFabCoar.interval(),
                    loFabFine,
                    loFabFine.interval(), copierCoarToFine);

#ifdef CH_USE_FLOAT
  Real eps = 1.0e-5;
#else
  Real eps = 1.0e-10;
#endif

  double normAnswer;
  if     ( CH_SPACEDIM == 1 ) normAnswer = 4.77;
  else if ( CH_SPACEDIM == 2 ) normAnswer = 9.54;
  else if ( CH_SPACEDIM == 3 ) normAnswer = 9.54;
  else if ( CH_SPACEDIM == 4 ) normAnswer = 9.54;
  else if ( CH_SPACEDIM == 5 ) normAnswer = 9.54;
  else if ( CH_SPACEDIM == 6 ) normAnswer = 9.54;
  else CH_assert( (CH_SPACEDIM>0) && (CH_SPACEDIM<7) );
  if ( Abs( norm( loFabCoar, Interval(0,0), 2 ) - normAnswer ) > eps )
  {
    if (verbose)
      {
        pout() << "norm(BoxLayoutData) is wrong, norm error ("
             << Abs( norm( loFabCoar, Interval(0,0), 2 ) - normAnswer )
             << ") is more than tolerance ("
             << eps
             << ")"
             << endl;
      }
    icode = -3;
  }

  for (ditFine.reset(); ditFine.ok(); ++ditFine)
    {
      FArrayBox& fab = loFabFine[ditFine()];
      Real rmax = fab.max(0);
      Real rmin = fab.min(0);
      Real valcomp;
      if (loLayoutFine.get(ditFine()) == b2)
        valcomp = -fluxval;
      else if (loLayoutFine.get(ditFine()) == b4)
        valcomp = -fluxval/2.0;
      else
        valcomp = 0.0;

      if (Abs(rmax - valcomp) > eps)
        {
          if (verbose)
            pout() << indent2 << "max is wrong" << endl;
          icode =-1;
        }
      if (Abs(rmin - valcomp) > eps)
        {
          if (verbose)
            pout() << indent2 << "min is wrong" << endl;
          icode =-2;
        }
    }
  if (verbose)
    {
      pout() << indent2 << "LoFabCoar:" << endl;
      for (ditCoar.reset(); ditCoar.ok(); ++ditCoar)
        {
          FArrayBox& fab = loFabCoar[ditCoar()];
          BoxIterator bit(fab.box());
          for (bit.reset(); bit.ok(); ++bit)
            {
              pout() << indent2 << bit() << indent2 << fab(bit(),0) << endl;
            }
        }
      pout() << indent2 << "LoFabFine:" << endl;
      for (ditFine.reset(); ditFine.ok(); ++ditFine)
        {
          FArrayBox& fab = loFabFine[ditFine()];
          BoxIterator bit(fab.box());
          for (bit.reset(); bit.ok(); ++bit)
            {
              pout() << indent2 << bit() << indent2 << fab(bit(),0) << endl;
            }
        }
    }
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
