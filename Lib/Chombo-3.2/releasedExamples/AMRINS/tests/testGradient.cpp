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
#include "LevelData.H"
#include "BoxIterator.H"
#include "BRMeshRefine.H"
#include "Gradient.H"
#include "SPMD.H"
#include "parstream.H"
#include "FABView.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

using std::cout;
using std::endl;

enum problemType
  {
    cubic = 0,
    sinusoidal,
    num_prob_types
  };

// for now, set this as a global variable
int probType = sinusoidal;

/// Prototypes:
int
testGradient();

void
parseTestOptions(int argc ,char* argv[]) ;

/// Global variables for handling output:
static const char* pgmname = "testDotProduct" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

void
initializePhi(LevelData<FArrayBox>& a_phi,
              const Real a_dx)
{
  RealVect offset(0.5*RealVect::Unit);

  Real Pi = 4.0*atan(1.0);
  DisjointBoxLayout grids = a_phi.getBoxes();
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisPhi = a_phi[dit];
      BoxIterator bit(grids[dit]);
      //BoxIterator bit(thisPhi.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect x(iv);
          x += offset;
          x *= a_dx;

          if (probType = cubic)
            {
              thisPhi(iv,0) = D_TERM(x[0]*x[0]*x[0],
                                     +x[1]*x[1]*x[1],
                                     +x[2]*x[2]*x[2]);
            }
          else if (probType = sinusoidal)
            {
              thisPhi(iv,0) = D_TERM(sin(2.0*Pi*x[0]),
                                     +sin(2.0*Pi*x[1]),
                                     +sin(2.0*Pi*x[2]));
            }
        }
    }
}


void
exactGradPhi(LevelData<FluxBox>& a_gradPhi,
             const Real a_dx)
{
  RealVect CCoffset(0.5*RealVect::Unit);

  Real Pi = 4.0*atan(1.0);
  DisjointBoxLayout grids = a_gradPhi.getBoxes();
  DataIterator dit = a_gradPhi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisGradPhi = a_gradPhi[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisGradPhiDir = thisGradPhi[dir];
          BoxIterator bit(thisGradPhiDir.box());

          RealVect faceOffset(CCoffset);
          faceOffset[dir] = 0.0;

          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect x(iv);
              x += faceOffset;
              x *= a_dx;

              for (int comp=0; comp<thisGradPhiDir.nComp(); comp++)
                {
                  if (probType == cubic)
                    {
                      thisGradPhiDir(iv,comp) = 3.0*x[comp]*x[comp];
                    }
                  else if (probType == sinusoidal)
                    {
                      thisGradPhiDir(iv,comp) = 2.0*Pi*cos(2.0*Pi*x[comp]);
                    }
                } // end loop over components
            } // end loop over faces
        } // end loop over face directions
    } // end loop over boxes
}

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testGradient();

  if ( status == 0 )
    pout() << indent << pgmname << " passed all tests." << endl ;
  else
    pout() << indent << pgmname << " failed with return code " << status << endl ;
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status ;
}

int
testGradient()
{
  int boxSize = 16 ;
  int domainSize = 32;

  int status = 0 ;

  // test internal boundary fix
  {
    Box domainBox(IntVect::Zero,
                  (domainSize-1)*IntVect::Unit);
    ProblemDomain probDomain(domainBox);
    // set to periodic domain
    for (int dir=0;dir<SpaceDim; dir++)
      {
        probDomain.setPeriodic(dir, true);
      }
    Vector<Box> splitBoxes;

    domainSplit(probDomain, splitBoxes, boxSize);

    Vector<int> procAssign(splitBoxes.size(), 0);

    DisjointBoxLayout splitGrids(splitBoxes, procAssign, probDomain);

    Vector<Box> unsplitBoxes(1, domainBox);
    Vector<int> unsplitProcAssign(1,0);

    DisjointBoxLayout unsplitGrids(unsplitBoxes, unsplitProcAssign,
                                   probDomain);

    Real dx = 1.0/domainSize;

    LevelData<FArrayBox> unsplitPhi(unsplitGrids, 1, IntVect::Unit);
    LevelData<FArrayBox> splitPhi(splitGrids, 1, IntVect::Unit);

    initializePhi(unsplitPhi, dx);
    initializePhi(splitPhi, dx);
    unsplitPhi.exchange();
    splitPhi.exchange();

    LevelData<FluxBox> unsplitGrad(unsplitGrids, SpaceDim, IntVect::Zero);
    LevelData<FluxBox> splitGrad(splitGrids, SpaceDim, IntVect::Zero);

    Gradient::levelGradientMAC(unsplitGrad, unsplitPhi, dx);
    Gradient::levelGradientMAC(splitGrad, splitPhi, dx);

    // now subtract one from the other; should be equal everywhere
    DataIterator singleDit = unsplitGrad.dataIterator();
    singleDit.begin();
    FluxBox& thisUnsplitGrad = unsplitGrad[singleDit];

    // now loop over split grids and subtract unsplit gradient
    Real maxDiff = 0.0;
    DataIterator splitDit = splitGrad.dataIterator();
    for (splitDit.begin(); splitDit.ok(); ++splitDit)
      {
        FluxBox& thisSplitGrad = splitGrad[splitDit];

        for (int dir=0; dir<SpaceDim; dir++)
          {
            FArrayBox& thisSplitGradDir = thisSplitGrad[dir];
            FArrayBox& thisUnsplitGradDir = thisUnsplitGrad[dir];

            thisSplitGradDir.minus(thisUnsplitGradDir,
                                   thisSplitGradDir.box(),
                                   0,0,SpaceDim);
            Real thisMaxDiff = thisSplitGradDir.norm(0,0,SpaceDim);
            if (abs(thisMaxDiff) > maxDiff) maxDiff = abs(thisMaxDiff);
          } // end loop over directions

      } // end loop over split grids


    Real errorTol = 1.0e-15;
    if (maxDiff > errorTol)
      {
        if (verbose)
          {
            pout() << "testGradient FAILED grid split independence test"
                   << " -- max diff = " << maxDiff
                   << endl;
          }
        status += 10;
      }
    else if (verbose)
      {
        pout() << "testGradient passed grid split independence test"
               << " -- Max(diff) = " << maxDiff
               << endl;
      }

  }


  { // test single-level convergence
    int numRefine = 3;

    Real oldMaxError[SpaceDim][SpaceDim];

    int thisDomainSize = boxSize;
    for (int n=0; n<numRefine; n++)
      {
        Box domainBox(IntVect::Zero,
                      (thisDomainSize-1)*IntVect::Unit);
        ProblemDomain probDomain(domainBox);
        // periodic domain for now
        for (int dir=0;dir<SpaceDim; dir++)
          {
            probDomain.setPeriodic(dir, true);
          }

        Vector<Box> splitBoxes;

        domainSplit(probDomain, splitBoxes, boxSize);

        Vector<int> procAssign(splitBoxes.size(), 0);

        DisjointBoxLayout grids(splitBoxes, procAssign,
                                probDomain);

        Real dx = 1.0/domainSize;

        LevelData<FArrayBox> phi(grids, 1, IntVect::Unit);

        initializePhi(phi, dx);
        phi.exchange();

        LevelData<FluxBox> grad(grids, SpaceDim, IntVect::Zero);
        LevelData<FluxBox> error(grids, SpaceDim, IntVect::Zero);

        Gradient::levelGradientMAC(grad, phi, dx);

        exactGradPhi(error, dx);

        // now subtract off computed gradPhi from exact to get error
        // also compute max(error) while we're at it.

      } // end loop over refinements
  } // end single-level convergence scope

  return status;
}

///
// Parse the standard test options (-v -q) out of the command line.
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
              argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              argv[i] = "" ;
            }
        }
    }
  return ;
}
