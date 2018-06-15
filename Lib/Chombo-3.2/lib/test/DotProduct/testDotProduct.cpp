#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FArrayBox.H"
#include "DotProdF_F.H"
#include "SPMD.H"
#include "parstream.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

using std::endl;

/// Prototypes:
int
testDotProduct();

void
parseTestOptions(int argc ,char* argv[]);

/// Global variables for handling output:
static const char* pgmname = "testDotProduct";
static const char* indent = "   ";
static const char* indent2 = "      ";
static bool verbose = true;

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions(argc ,argv);
  if (verbose)
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl;

  int status = testDotProduct();

  if (status == 0)
    pout() << indent << pgmname << " passed all tests." << endl;
  else
    pout() << indent << pgmname << " failed with return code " << status << endl;
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}

int
testDotProduct()
{
  const int boxSize = 3;
  int status = 0;

  // test non-empty box
  {

    Box box(IntVect::Zero ,boxSize * IntVect::Unit);

    // test single component, zero-valued fabs
    FArrayBox fab1(box,1) ,fab2(box,1);

    fab1.setVal((Real)0.0); fab2.setVal((Real)0.0);

    Real dot_prod = (Real)-9988.7766;
    int component = 0;

    FORT_DOTPRODUCT(CHF_REAL(dot_prod),
                     CHF_CONST_FRA(fab1) ,CHF_CONST_FRA(fab2),
                     CHF_BOX(box),
                     CHF_CONST_INT(component),
                     CHF_CONST_INT(component)
);

    if (dot_prod == (Real)0.0)
      {
        if (verbose)
        {
          pout() << indent << pgmname << " passed zero-valued FAB test." << endl;
        }
      }
    else
      {
        pout() << indent << pgmname << " failed zero-valued FAB test." << endl;
        if (verbose)
        {
          pout() << indent2 << "dot_product result was " << dot_prod
                 << " but should have been 0.0" << endl;
        }
        status = -1;

      }

    // test single component, nonzero-valued fabs
    fab1.setVal((Real)1.0); fab2.setVal((Real)2.0);

    dot_prod = (Real)-9988.7766;

    component = 0;
    FORT_DOTPRODUCT(CHF_REAL(dot_prod),
                     CHF_CONST_FRA(fab1) ,CHF_CONST_FRA(fab2),
                     CHF_BOX(box),
                     CHF_CONST_INT(component),
                     CHF_CONST_INT(component)
);

    Real exact_dp = box.numPts() * (Real)2.0;

    if (dot_prod == exact_dp)
      {
        if (verbose)
        {
          pout() << indent << pgmname << " passed nonzero-valued FAB test." << endl;
        }
      }
    else
      {
        pout() << indent << pgmname << " failed nonzero-valued FAB test." << endl;
        if (verbose)
        {
          pout() << indent2 << "dot_product result was " << dot_prod
                 << " but should have been " << exact_dp << endl;
        }
        status = -1;
      }

    // test multicomponent, nonzero-valued fabs
    fab1.define(box,3); fab2.define(box,3);
    fab1.setVal((Real)1.0   ,0); fab2.setVal((Real)2.0   ,0);
    fab1.setVal((Real)10.0  ,1); fab2.setVal((Real)20.0  ,1);
    fab1.setVal((Real)100.0 ,2); fab2.setVal((Real)200.0 ,2);

    dot_prod = (Real)0.0;

    int startcomp = 0;
    int endcomp = fab1.nComp()-1;

    FORT_DOTPRODUCT(CHF_REAL(dot_prod),
                     CHF_CONST_FRA(fab1),
                     CHF_CONST_FRA(fab2),
                     CHF_BOX(box),
                     CHF_CONST_INT(startcomp),
                     CHF_CONST_INT(endcomp)
);

    exact_dp = (box.numPts() * ((Real)2.0 + (Real)200.0 + (Real)20000.0));

    if (dot_prod == exact_dp)
      {
        if (verbose)
        {
          pout() << indent << pgmname << " passed nonzero-valued, multicomponent FAB test." << endl;
        }
      }
    else
      {
        pout() << indent << pgmname << " failed nonzero-valued, multicomponent FAB test." << endl;
        if (verbose)
        {
          pout() << indent2 << "dot_product result was " << dot_prod
                 << " but should have been " << exact_dp << endl;
        }
        status = -1;
      }
  }

  return status;
}

///
// Parse the standard test options (-v -q) out of the command line.
///
void
parseTestOptions(int argc ,char* argv[])
{
  for (int i = 1; i < argc; ++i)
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp(argv[i] ,"-v" ,3) == 0)
            {
              verbose = true;
              // argv[i] = "";
            }
          else if (strncmp(argv[i] ,"-q" ,3) == 0)
            {
              verbose = false;
              // argv[i] = "";
            }
        }
    }
  return;
}
