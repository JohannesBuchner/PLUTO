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
#include "Vector.H"
#include "CH_Timer.H"
#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "MayDay.H"
#include "ProblemDomain.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "LevelData.H"

#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:
void
parseTestOptions( int argc ,char* argv[] );

int testMU(void);

/// Global variables for handling output:
static const char *pgmname = "testMemoryUsage";
static const char *indent2 = "      ";
static bool verbose = true;

/// Code:
int main(int argc, char* argv[])
{
  CH_TIMERS("mymain");
  CH_TIMER("wholething", t1);
  CH_START(t1);

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions(argc,argv);

  if ( verbose ) pout() << indent2 << "Beginning " << pgmname << " ..." << endl;

  ///
  // Run the tests
  ///
  int ret = testMU();
#ifdef CH_MPI
  MPI_Finalize();
#endif

  double timeForCode=0;
  CH_STOPV(t1, timeForCode);

  return ret;
}

void printMemoryInfo(char *s)
{
#ifdef CH_USE_MEMORY_TRACKING
  Real residentSetSize=0.0;
  Real size=0.0;
  getMemoryUsageFromOS(residentSetSize, size);
  char temp[120];
  sprintf(temp, "%30s res=%10.4f MB  s=%10.4f MB\n", s, residentSetSize, size);

  string seperator = "############################################################################################################";
  pout() << seperator << endl;
  pout() << temp;
  ReportAllocatedMemory(pout());
#endif
}

int testMU(void)
{
  CH_TIME("testMU");
  int return_code = 0;

#ifdef CH_USE_MEMORY_TRACKING
  const int N3=1000;
  const int N5=100000;
  const int N6=1000000;

  printMemoryInfo("begin");

  Vector<Real> a(N6);
  for (int i=0; i<N6; i++) a[i] = i;

  printMemoryInfo("after Vector<Real> a(N6)");

  Vector<Real> b(N6);
  for (int i=0; i<N6; i++) b[i] = i;

  printMemoryInfo("after Vector<Real> b(N6)");

  Vector<IntVect> vectIV(N6);

  printMemoryInfo("after Vector<IntVect> vectIV(N6)");

  Vector<Box> vectBox(N6);
  for (int i=0; i<N6; i++) vectBox[i] = Box(IntVect::Zero, 20*IntVect::Unit);

  printMemoryInfo("after Vector<Box> vectBox(N6)");


  const int bf=8;
  const int mega=10*8-1;
  pout() << " mega=" << mega << " (mega+1)*(mega+1)*(mega+1)=" << (mega+1)*(mega+1)*(mega+1) << endl;
  Box bigBox = Box(IntVect::Zero, mega*IntVect::Unit);
  ProblemDomain baseLevelDomain(bigBox);
  Vector<Box> baseLevelBoxes;
  domainSplit(baseLevelDomain, baseLevelBoxes, bf, bf);
  pout() << "Num baseLevelBoxes=" << baseLevelBoxes.size() << endl;

  printMemoryInfo("after domainSplit");

  const int numLevels = 1;
  Vector<int> ranks;
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    LoadBalance(ranks, baseLevelBoxes);
  }

  printMemoryInfo("after LoadBalance");

  ProblemDomain pd(bigBox);
  Vector<DisjointBoxLayout> vectGrids(numLevels);
  for (int ilev = 0; ilev < numLevels; ilev++)
  {
    vectGrids[ilev] = DisjointBoxLayout(baseLevelBoxes, ranks, pd);
  }

  printMemoryInfo("after DisjointBoxLayouts");

  Vector<LevelData<FArrayBox>* > rhs(numLevels, NULL);

  const int ncomp = 1;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      rhs[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], ncomp, IntVect::Zero);
      LevelData<FArrayBox>& rhsfabs= *rhs[ilev];
      for (DataIterator dit = rhs[ilev]->dataIterator(); dit.ok(); ++dit)
        {
          rhsfabs[dit()].setVal(0.0);
        }
    }

  printMemoryInfo("after rhs allocation");

  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      delete rhs[ilev];
    }

  printMemoryInfo("after delete rhs");

#endif
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
