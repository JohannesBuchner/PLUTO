#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "1DPotentialSolve.H"

#define SAVE_SPACEDIM CH_SPACEDIM
#include "LevelData.H.multidim"
#include "FArrayBox.H.multidim"
#include "BCFunc.H.multidim"
#include "AMRPoissonOp.H.multidim"
#include "AMRMultiGrid.H.multidim"
#include "BiCGStabSolver.H.multidim"
#include "CoarseAverage.H.multidim"
#include "computeSum.H.multidim"

// At this point, CH_SPACEDIM=0; the .H.multidim files do that.

#include "Slicing.H.transdim"
#include "Injection.H.transdim"

// Declare the dimensionally-correct FORT_GRADCC:
//#define SAVE_SPACEDIM CH_SPACEDIM
#undef CH_SPACEDIM
#define CH_SPACEDIM 1
#include "FORT_PROTO.H"
#include "GradientF.1D_F.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM SAVE_SPACEDIM
#undef SAVE_SPACEDIM



using namespace Chombo;
using namespace CH_MultiDim;
using std::endl;

int doOneDimensionalPotentialSolve(Vector<D2::LevelData<D2::FArrayBox>* >& a_rhs,
                                   Vector<D2::LevelData<D2::FArrayBox>* >& a_gradient,
                                   Vector<D2::ProblemDomain>& a_amrDomains,
                                   Vector<Real>& a_amrDx,
                                   Vector<int>& a_amrRefRatios,
                                   int a_finest_level,
                                   int a_lBase,
                                   int a_numComps,
                                   int a_verbosity)
{

  if (a_verbosity >= 3)
    {
      pout() << "In 1dPotentialSolve" << endl;
    }
  int status =0;

  // first, take rhs to 1D
  Vector<D1::LevelData<D1::FArrayBox>* > oneDrhs(a_rhs.size(), NULL);
  Vector<D1::LevelData<D1::FArrayBox>* > oneDPotential(a_rhs.size(), NULL);
  Vector<D1::DisjointBoxLayout > oneDGrids(a_amrDomains.size());
  Vector<D1::ProblemDomain> oneDDomains(a_amrDomains.size());


  // this is a little silly, but necessary
  Vector<int> oneDRefRatios(a_amrRefRatios.size());
  for (int lev=0; lev<oneDRefRatios.size(); lev++)
    {
      oneDRefRatios[lev] = a_amrRefRatios[lev];
    }

  if (a_verbosity >= 3)
    {
      pout() << "entering AMR memory allocation..." << endl;
    }

  Vector<Real> oneDamrDx(a_amrDx.size());
  for (int lev=0; lev<= a_finest_level; lev++)
    {
      oneDamrDx[lev] = a_amrDx[lev];
      // slice in the y-direction
      D2::SliceSpec slice(1,0);

      //const D2::DisjointBoxLayout& LevelGrids2d = a_rhs[lev]->getBoxes();
      //sliceDisjointBoxLayout(oneDGrids[lev], LevelGrids2d,
      //                       slice);


      D1::Box oneDDomainBox;
      bool testBool;
      D2::ProblemDomain twodDomain= a_amrDomains[lev];
      D1::ProblemDomain onedDomain;
      onedDomain = sliceDomain(twodDomain,
                                slice, &testBool);
      oneDDomains[lev] = onedDomain;

      //      oneDDomains[lev].define(oneDDomainBox);
      //oneDDomains[lev].setPeriodic(0,a_amrDomains[lev].isPeriodic(0));

      oneDrhs[lev] = new D1::LevelData<D1::FArrayBox>;
      D1::LevelData<D1::FArrayBox>& oneDRhsLevel = *oneDrhs[lev];
      D2::LevelData<D2::FArrayBox>& twoDRhsLevel = *(a_rhs[lev]);

      // define oneDRhs here
      //int nComp = twoDRhsLevel.nComp();
      D1::IntVect oneDGhost;
      sliceIntVect(oneDGhost, twoDRhsLevel.ghostVect(), slice);
      //oneDRhsLevel.define(oneDGrids[lev], nComp, oneDGhost);

      sliceLevelData(oneDRhsLevel, twoDRhsLevel, slice);

      oneDGrids[lev] = oneDRhsLevel.getBoxes();
      // debug check
      if (a_verbosity > 4)
        {
          D1::pout () << oneDGrids[lev];
        }

      oneDPotential[lev] = new D1::LevelData<D1::FArrayBox>(oneDGrids[lev],
                                                            oneDRhsLevel.nComp(),
                                                            D1::IntVect::Unit);
    }

  if (a_verbosity >= 3)
    {
      pout() << "Memory allocated..." << endl;
    }

  // set up boundary conditions for potential solve; for now, hardwire
  // to neumann.
  D1::IntVect NeumBCtype = D1::IntVect::Zero;
  D1::RealVect bcVal = D1::RealVect::Zero;

  RefCountedPtr<D1::BCFunction> potentialBC = D1::ConstDiriNeumBC(NeumBCtype,
                                                                  bcVal,
                                                                  NeumBCtype,
                                                                  bcVal);

  if (a_verbosity >= 3)
    {
      pout() << "before sumRHS" << endl;
    }

  // since we're either doing Neumann or periodic BC's, rescale rhs to
  // ensure solvability
  Real sumRHS = D1::computeSum(oneDrhs,
                               oneDRefRatios,
                               oneDamrDx[0]);

  if (a_verbosity >= 3)
    {
      pout() << "after sumRHS" << endl;
    }

  // multiply by domain size to get integral of RHS
  sumRHS /= (oneDDomains[0].domainBox().size(0)*oneDamrDx[0]);

  for (int lev=0; lev<oneDrhs.size(); lev++)
    {
      D1::LevelData<D1::FArrayBox>& levelRHS = *oneDrhs[lev];
      D1::DataIterator levelDit = levelRHS.dataIterator();
      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          Real scale = sumRHS;
          levelRHS[levelDit] -= scale;
        }
    }

  if (a_verbosity >= 3)
    {
      pout() << "begin second sumRHS" << endl;
    }

  Real newSumRHS = D1::computeSum(oneDrhs, oneDRefRatios,
                                  oneDamrDx[0]);

  pout() << "old sum(rhs) = " << sumRHS
         << ",  new sum(rhs) after rescaling = " << newSumRHS
         << endl;




  Real alpha = 0.0;
  Real beta = 1.0;
  D1::AMRPoissonOpFactory oneDfactory;
  oneDfactory.define(oneDDomains[0],
                     oneDGrids,
                     oneDRefRatios,
                     oneDamrDx[0],
                     potentialBC,
                     alpha, beta);

  D1::BiCGStabSolver<D1::LevelData<D1::FArrayBox> > oneDbottomSolver;
  oneDbottomSolver.m_verbosity = a_verbosity - 2;


  // need to allocate AMRSolver for this
  // use poissop because it already has the correct boundary conditions
  D1::AMRMultiGrid<D1::LevelData<D1::FArrayBox> > potentialSolver;
  potentialSolver.define(oneDDomains[0],
                         oneDfactory,
                         &oneDbottomSolver,
                         a_finest_level+1);

  potentialSolver.m_verbosity = a_verbosity - 1;

  if (a_verbosity >= 4)
    {
      pout() << "starting AMR solve..." << endl;
    }

  // solve for potential
  bool zeroPhi = true;
  potentialSolver.solve(oneDPotential, oneDrhs,
                        a_finest_level, a_lBase, zeroPhi);

  if (a_verbosity >= 4)
    {
      pout() << "Out of AMR solve" << endl;
    }


  // average potential from finer-coarser, if necessary
  for (int lev=a_finest_level; lev>0; lev--)
    {
      D1::CoarseAverage averager(oneDGrids[lev],
                                 oneDGrids[lev-1],
                                 1,
                                 oneDRefRatios[lev-1]);

      averager.averageToCoarse(*oneDPotential[lev-1],
                               *oneDPotential[lev]);
    }

  // re-apply physical boundary conditions to potential
  // and then compute 1D gradient
  for (int lev=0; lev<=a_finest_level; lev++)
    {
      const D1::ProblemDomain& levelDomain = oneDDomains[lev];
      D1::LevelData<D1::FArrayBox>& levelPotential = *oneDPotential[lev];
      // store gradient in RHS space
      D1::LevelData<D1::FArrayBox>& levelGradient = *oneDrhs[lev];
      const D1::DisjointBoxLayout& levelGrids = levelGradient.getBoxes();
      bool isHomogeneous = false;

      D1::DataIterator dit = levelPotential.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          (*potentialBC)(levelPotential[dit()],
                         levelGrids[dit],
                         levelDomain,
                         a_amrDx[lev],
                         isHomogeneous);
        }

      levelPotential.exchange();

      // compute 1D gradients in x-direction
      // initially, just do CC gradient
      for (dit.begin(); dit.ok(); ++dit)
        {
          D1::Box gridBox = levelGrids[dit()];
          int gradDir = 0;
#define SAVE_SPACEDIM CH_SPACEDIM
#undef CH_SPACEDIM
#define CH_SPACEDIM 1

          FORT_GRADCC(CHF_FRA1(levelGradient[dit()],0),
                      CHF_CONST_FRA1(levelPotential[dit()],0),
                      CHF_BOX(gridBox),
                      CHF_CONST_REAL(a_amrDx[lev]),
                      CHF_INT(gradDir));

#undef CH_SPACEDIM
#define CH_SPACEDIM SAVE_SPACEDIM
#undef SAVE_SPACEDIM

          // really want the negative of the gradient
          levelGradient[dit] *= -1;
        }


      // this is where we can bring things back from 1D
      D2::SliceSpec slice(1,0);
      injectLevelData(*(a_gradient[lev]), levelGradient, slice);

    } // end loop over levels

  if (a_verbosity >= 3)
    {
      pout() << "Leaving 1DPotentialSolve" << endl;
    }

  return status;
}


