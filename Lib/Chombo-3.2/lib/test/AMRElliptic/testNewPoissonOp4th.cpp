#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <iomanip>
#include <sstream>

#include "BCFunc.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "FABView.H"
#include "MultiGrid.H"
#include "NewPoissonOp4.H"
#include "NewPoissonOp.H"
#include "RelaxSolver.H"

#include "UsingNamespace.H"


/**
 * Test fourth-order Multigrid Poisson solver:
 *    Laplace(u) = rhs
 *  with periodic boundary condition and
 *  exact u := \sum_(j=0)^(n_k) \sum_(d=0)^(D-1) sin(k_j x_d)
 *
 *  Qinghai Zhang MAY/20/2009
 */

const Real k[] =
{
  2*M_PI, 4*M_PI, 8*M_PI
};
//const Real k[] = { 2*M_PI};
const int nk = sizeof(k)/sizeof(k[0]);

void complexSin(Real* x, int* dir, Side::LoHiSide* side, Real* a_values)
{
  Real v = 0;
  for (int i=0; i<nk; i++)
    {
      v += D_TERM(sin(k[i]*x[0]), +sin(k[i]*x[1]), +sin(k[i]*x[2]));
    }
  a_values[0] = v;
}

// if forRHS=false, return \int u at x
// if forRHS=true, return \int Lap(u) at x
Real integralSin(const RealVect& x, bool forRHS)
{
  Real v = 0;
  for (int i=0; i<nk; i++)
    {
      v += (D_TERM(cos(k[i]*x[0]), +cos(k[i]*x[1]), +cos(k[i]*x[2])))
        * (forRHS ? k[i] : -1/k[i]);
    }
  return v;
}

// x is the cell center, dx is the cell width.
Real cellAvgSin(const RealVect& x, Real dx)
{
  return (integralSin(dx/2+x,false)-integralSin(-dx/2+x,false))/dx;
}

// cell avergage of RHS
Real cellAvgRHS(const RealVect& x, Real dx)
{
  return (integralSin(dx/2+x,true)-integralSin(-dx/2+x, true))/dx;
}

// periodic boundary condition for FArrayBox
void periodicSinBC(FArrayBox& a_state, const Box& a_valid,
                   const ProblemDomain& a_domain, Real a_dx, bool a_homo)
{
  const int nGhost = a_valid.smallEnd(0) - a_state.box().smallEnd(0);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Box toRegion = adjCellBox(a_valid, idir, sit(), nGhost);
          toRegion &= a_state.box();

          IntVect shift = IntVect::Zero;
          shift[idir] -= sign(sit())*a_valid.size(idir);
          for (BoxIterator bit(toRegion); bit.ok(); ++bit)
            {
              a_state(bit()) = a_state(bit()+shift);
            }
        } // end loop over sides
    } // end loop over directions
}

template <class T>
const std::string paddedStr(T r, std::streamsize n)
{
  std::ostringstream os;
  os << std::setw(n) << r;
  return os.str();
}

int main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif

  // test parameters
  const int testOrder = 4;  // order-of-convergence of the solver
#ifdef CH_USE_DOUBLE
  const int nGrids = 3;
#endif
  // single-precision numbers have a short mantissa
  //   and a convergence test on 3 grids is probably too much
#ifdef CH_USE_FLOAT
  const int nGrids = 2;
#endif
  const int nCells0 = 32;
  const RealVect xLo = RealVect::Zero;
  const Real xHi = 1.0;
  const int nGhosts = testOrder/2;
  const int resNT = 0;      // norm Type
  const int errNT = 0;
  const int maxCycles = 50;
  // A test is considered as a failure
  //  if its convergence rate is smaller than below.
  const Real targetConvergeRate = testOrder*0.9;
  // A test is considered as converged
  //  if the inital residual is reduced by the following factor.
#ifdef CH_USE_DOUBLE
  const Real residualConvergeRatio = 1.e-11;
#endif
  // single-precision numbers have a short mantissa,
  //  thus the ratio by which the residual can be reduced should be bigger.
#ifdef CH_USE_FLOAT
  const Real residualConvergeRatio = 2.e-5;
#endif

  // solver parameters
  const int nRelax = (SpaceDim>1) ? 2 : 3;          // m_pre=m_post
  // cycle types
  const std::string cycleStr[] =
  {
    "  V" , "FMG"
  };
  const int cycleType[] =
  {
    1, -1
  };
  const int nCycleType = 2;       // set to 1 for V-Cycle only.

  // test results holders
  int nCycles[nGrids];            // number of cycles to convergence.
  for (int i=0; i<nGrids; i++)    // initialize nCycles to maxCycles.
    nCycles[i] = maxCycles;
  Real resNorm[nGrids][maxCycles+1], errNorm[nGrids][maxCycles+1];
  Real convergeRate[nGrids-1][2];
  const Real log2r = 1.0/log(2.0);

  // The big nested loop for different cases.
  int status = 0;                      // number of errors detected.
  for (int j=0; j<nCycleType; j++)
  {
    pout() << "\n********************************************************\n"
           << "\nTesting MultiGrid::oneCycle(correction, residual)"
           << "; DIM=" << SpaceDim << std::endl
           << " cycle type = " << cycleStr[j]
           << "; m_pre = m_post = " << nRelax << std::endl;

    for (int iGrid=0; iGrid<nGrids; iGrid++)
    {
      int ref = 1;
      for (int i=0; i<iGrid; i++) ref*=2;
      const Real dx = xHi/nCells0/ref;
      const Box domain = refine(Box(IntVect::Zero, (nCells0-1)*IntVect::Unit),
                                ref);
      const Box ghostBox = grow(domain,nGhosts);

      pout() << "\n----------------------------------------------------\n";
      pout() << "nCells = " << nCells0*ref << " ; dx = " << dx << " \n";

      FArrayBox phi(ghostBox, 1);
      FArrayBox correction(ghostBox, 1);
      FArrayBox rhs(domain, 1);
      FArrayBox error(domain, 1);
      FArrayBox phiExact(domain, 1);
      FArrayBox residual(domain, 1);

      // set initial guess
      phi.setVal(0.0);
      // set RHS and the exact solution
      for (BoxIterator bit(domain); bit.ok(); ++bit)
        {
          const RealVect offset = bit()-domain.smallEnd();
          const RealVect x = xLo + dx*(0.5+offset);
          rhs(bit()) = cellAvgRHS( x, dx );
          phiExact(bit()) = cellAvgSin( x, dx );
        }

      // Initialize big objects
      MultiGrid<FArrayBox> solver;
      solver.m_numMG = 1;
      solver.m_bottom = nRelax;
      solver.m_pre = nRelax;
      solver.m_post = nRelax;
      solver.m_cycle = cycleType[j];
      // BiCGStab coupled with period BC will cause problems
      //  close to truncation error.
      //      BiCGStabSolver<FArrayBox> bottomSolver;
      RelaxSolver<FArrayBox> bottomSolver;
      bottomSolver.m_verbosity = 0;
      MGLevelOp<FArrayBox>* op = 0;
      if (testOrder==4)
      {
        NewPoissonOp4Factory opFactory;
        opFactory.define(dx*RealVect(IntVect::Unit), periodicSinBC);
        solver.define(opFactory, &bottomSolver, domain);
        op = opFactory.MGnewOp(domain,0);
      }
      else if (testOrder==2)
      {
        NewPoissonOpFactory opFactory;
        opFactory.define(dx*RealVect(IntVect::Unit), periodicSinBC);
        solver.define(opFactory, &bottomSolver, domain);
        op = opFactory.MGnewOp(domain,0);
      }

      // put the data into residual-correction form
      op->residual(residual, phi, rhs);
      resNorm[iGrid][0] = residual.norm(resNT);
      op->axby(error, phi, phiExact, 1, -1);
      errNorm[iGrid][0] = error.norm(errNT);
      solver.init(correction, residual);

      const std::streamsize columnSZ = 16;
      const std::string tab(" |");
      // The table header of results
      std::ostringstream oss[4];
      oss[0] << cycleStr[j] << "-Cycle N.O.";
      oss[1] << "Residual " << resNT << "-norm";
      oss[2] << "Error " << errNT << "-norm";
      oss[3] << "Residual-ratio";
      for (int i=0; i<4; i++)
        {
          pout() << paddedStr(oss[i].str(), columnSZ) << tab;
        }
      // The initial residual and error norms
      pout() << std::endl << paddedStr(0, columnSZ) << tab
             << paddedStr(resNorm[iGrid][0], columnSZ) << tab
             << paddedStr(errNorm[iGrid][0], columnSZ) << tab << std::endl;

      // Solve the problem using MultiGrid::oneCycle
      const Real tiny = resNorm[iGrid][0]*residualConvergeRatio;
      for (int i=1; i<=maxCycles; i++)
        {
          correction.setVal(0.0);
          solver.oneCycle(correction, residual);
          op->incr(phi, correction, 1);
          op->residual(residual, phi, rhs);
          resNorm[iGrid][i] = residual.norm(resNT);
          op->axby(error, phi, phiExact, 1, -1);
          // substract the averaged error for periodic BC.----
          error -= error.get(domain.smallEnd(),0);
          // -------------------------------------------------
          errNorm[iGrid][i] = error.norm(errNT);
          const Real ratio = resNorm[iGrid][i-1]/resNorm[iGrid][i];
          pout() << paddedStr(i, columnSZ) << tab
                 << paddedStr(resNorm[iGrid][i], columnSZ) << tab
                 << paddedStr(errNorm[iGrid][i], columnSZ) << tab
                 << paddedStr(ratio, columnSZ) << tab << std::endl;
          // stop the iteration if the residual is small enough.
          if (resNorm[iGrid][i]<tiny)
          {
            nCycles[iGrid] = i;
            break;
          }
        }
      delete op;
    } // end grid loop

    pout() << "\nConvergence Rates :\n";
    for (int i=0; i<nGrids-1; i++)
    {
      const Real ratio = errNorm[i][nCycles[i]]/errNorm[i+1][nCycles[i+1]];
      convergeRate[i][j] = log(ratio)*log2r;
      if (convergeRate[i][j] < targetConvergeRate)
      {
        status += 1;
      }
      pout() << "    " << convergeRate[i][j] << std::endl;
    }
  }// end cycle type

  if (status==0)
  {
    pout() <<  "All tests passed!\n";
  }
  else
  {
    pout() <<  status << " tests failed!\n";
  }

#ifdef CH_MPI
  MPI_Finalize ();
#endif

  return status;
}
