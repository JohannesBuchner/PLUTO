#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BCFunc.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "FABView.H"
#include "MultiGrid.H"
#include "NewPoissonOp.H"
#include "UsingNamespace.H"

/**
 *  Multigrid Poisson solver:
 *    Laplace(u) = rhs
 *
 *  A numerical test on page 64 of
 *   `A MULTIGRID TUTORIAL by Briggs, Henson & McCormick'
 *   is duplicated here.
 *  For 3D
 *  exact u := (x^2-x^4)*(y^2-y^4)*(z^2-z^4)
 *
 *  Qinghai Zhang Apr/22/2009
 */

// RHS
Real rhsFunc(const RealVect& a_x)
{
  Real r = 0.0;
  const RealVect x2 = a_x*a_x;
  const RealVect a = x2*(1-x2);
  const RealVect b = 2-12*x2;
  // cross products of coordinates prevents the use of D_TERM
  if (SpaceDim==1)
    r = b[0];
  else if (SpaceDim==2)
    r = b[0]*a[1] + a[0]*b[1];
  else if (SpaceDim==3)
    r = b[0]*a[1]*a[2] + a[0]*b[1]*a[2] + a[0]*a[1]*b[2];
  else
    MayDay::Error("Invalid Dimension!");
  return r;
}

// exact u
Real exactSolution(const RealVect& a_x)
{
  RealVect e = a_x*a_x;
  e *= 1. - e;
  return e.product();
}

// Constant Dirichlet Boundary condition
void constDiri(Real* pos, int* dir, Side::LoHiSide* side, Real* a_values)
{
  a_values[0] = exactSolution(RealVect(D_DECL(pos[0], pos[1], pos[2])));
}

void constDiriBC(FArrayBox& a_state, const Box& valid,
                 const ProblemDomain& a_domain,
                 Real a_dx, bool a_homogeneous)
{
  // changing the extrapolation order from 1 to 2
  //  would hurt the convergence in 1D.
  for (int i=0; i<SpaceDim; ++i)
    {
      DiriBC(a_state,valid,a_dx,a_homogeneous, constDiri, i, Side::Lo, 1);
      DiriBC(a_state,valid,a_dx,a_homogeneous, constDiri, i, Side::Hi, 1);
    }
}


int main(int argc, char* argv[])
{

#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif

  // test parameters
  const int nGrids = 3;
  const int nCells0 = 32;
  // xLo has to be zero in order for DiriBc to work.
  const RealVect xLo = RealVect::Zero;
  const Real xHi = 1.0;
  const Box box0(IntVect::Zero, (nCells0-1)*IntVect::Unit);
  const int nGhosts = 1;
  const int resNT = 2; // norm Type
  const int errNT = 0;
  // A test is considered as a failure
  // if its convergence rate is smaller than below.
  const Real targetConvergeRate = 1.75;

  // solver parameters
  // To converge within 10 V-Cycles in 1D,
  //  nRelax=3 is the minimum number of relaxations.
  const int nRelax = 3; // m_pre=m_post
  // cycle Type, 1 : V-Cycle; -1 : FMG-Cycle
  const int cycleType[2] =
  {
    1, -1
  };
  const std::string cycleStr[2] =
  {
    "   V" , " FMG"
  };

  // test results holder
  const int nCycles[2] =
  {
    9, 5
  };
  const int maxCycles = 10; // > max(nCycles)
  //  Real resNorm[nGrids][nCycles+1], errNorm[nGrids][nCycles+1];
  Real resNorm[nGrids][maxCycles], errNorm[nGrids][maxCycles];
  Real convergeRate[nGrids-1][2];
  const Real log2r = 1.0/log(2.0);

  // status records the number of errors detected.
  int status = 0;
  for (int j=0; j<2; j++)
  {
    pout() << "\n**************************************************\n"
           << "\nTesting MultiGrid::oneCycle(correction, residual)\n"
           << " cycle type = " << cycleStr[j]
           << "; m_pre = m_post = " << nRelax << "\n";

    for (int iGrid=0; iGrid<nGrids; iGrid++)
    {
      int ref = 1;
      for (int i=0; i<iGrid; i++)
        ref*=2;
      const Real dx = xHi/nCells0/ref;
      const Box domain = refine(box0,ref);
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
          rhs(bit()) = rhsFunc( x );
          phiExact(bit()) = exactSolution( x );
        }

      // Initialize big objects
      NewPoissonOpFactory opFactory;
      opFactory.define(dx*RealVect(IntVect::Unit), constDiriBC);
      MultiGrid<FArrayBox> solver;
      BiCGStabSolver<FArrayBox> bottomSolver;
      bottomSolver.m_verbosity = 0;
      MGLevelOp<FArrayBox>* op = opFactory.MGnewOp(domain,0);
      solver.m_numMG = 1;
      solver.m_bottom = 1;
      solver.m_pre = nRelax;
      solver.m_post = nRelax;
      solver.m_cycle = cycleType[j];
      solver.define(opFactory, &bottomSolver, domain);

      // put the data into residual-correction form
      op->residual(residual, phi, rhs);
      resNorm[iGrid][0] = residual.norm(resNT);
      op->axby(error, phi, phiExact, 1, -1);
      errNorm[iGrid][0] = error.norm(errNT);
      solver.init(correction, residual);

      // Solve the problem using MultiGrid::oneCycle
      for (int i=0; i<nCycles[j]; i++)
        {
          correction.setVal(0.0);
          solver.oneCycle(correction, residual);
          op->incr(phi, correction, 1);
          op->residual(residual, phi, rhs);
          resNorm[iGrid][i+1] = residual.norm(resNT);
          op->axby(error, phi, phiExact, 1, -1);
          errNorm[iGrid][i+1] = error.norm(errNT);
        }
      delete op;

      // output a table of results
      pout()<< cycleStr[j] << "-Cycle N.O. |  residual " << resNT
            << "-norm  |  Error " << errNT << "-norm  \n";
      for (int i=0; i<nCycles[j]+1; i++)
        {
          pout() << "         " << i << "      |    " << resNorm[iGrid][i]
                 << "        |    " << errNorm[iGrid][i]  << "\n";
        }
    } // end grid loop

    pout() << "\nConvergence Rate based on the error in the last cycle:\n";
    for (int i=0; i<nGrids-1; i++)
    {
      Real ratio = errNorm[i][nCycles[j]]/errNorm[i+1][nCycles[j]];
      convergeRate[i][j] = log(ratio)*log2r;
      if (convergeRate[i][j] < targetConvergeRate)
      {
        status += 1;
      }
      pout() << "    " << convergeRate[i][j] << "\n";
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
