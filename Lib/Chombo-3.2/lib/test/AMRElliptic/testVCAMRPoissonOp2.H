#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TESTVCAMRPOISSONOP2_H_
#define _TESTVCAMRPOISSONOP2_H_

// testVCAMRPoissonOp2 -- A test suite for the Revenge of the Variable-Coefficient
// AMR Poisson Operator (VCAMRPoissonOp2) of the form
// L(phi) = alpha A phi - beta div (B grad phi).

#include <iostream>
using std::endl;

#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "parstream.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "VCAMRPoissonOp2.H"
#include "AnalyticForms.H"
#include "BCFunc.H"
#include "DebugOut.H"
#include "AMRIO.H"
#include <cxxtest/TestSuite.h>
#include "UsingNamespace.H"

/// Application-specific global variables:
//#define NUM_RESOLUTIONS 1
#define NUM_RESOLUTIONS 3
#if CH_SPACEDIM == 2
//static int resolutions[NUM_RESOLUTIONS] = {64};
static int resolutions[NUM_RESOLUTIONS] =
{
  64, 128, 256
};
#else
//static int resolutions[NUM_RESOLUTIONS] = {32};
static int resolutions[NUM_RESOLUTIONS] =
{
  32, 64, 128
};
#endif
static int blockingFactor = 8;

// This test suite class is used to generate a C++ test harness with CxxTest.
class TestVCAMRPoissonOp2: public CxxTest::TestSuite
{
  public:

  // This sets up the grids used by the various tests below.
  void setUp()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      // Create boxes and problem domains on which these tests run.
      int N = resolutions[ires];
      m_dxs[ires] = 1.0 / (2*N);
      m_domains[ires] = Box(IntVect(D_DECL(-N,-N,-N)), IntVect(D_DECL(N-1, N-1, N-1)));
      m_problemDomains[ires] = ProblemDomain(m_domains[ires]);

      // Now initialize the disjoint box layouts.
      BRMeshRefine br;
      Box domain = m_domains[ires];
      domain.coarsen(blockingFactor);
      domain.refine(blockingFactor);
      domain.coarsen(blockingFactor);

      // Carve a funky shape in the domain.
      ProblemDomain junk(domain);
      IntVectSet pnd(domain);
      IntVectSet tags;
      for (BoxIterator bit(domain); bit.ok(); ++bit)
      {
        const IntVect& iv = bit();
        if (D_TERM(true, && iv[1]< 2*iv[0] && iv[1]>iv[0]/2, && iv[2] < domain.bigEnd(2)/2))
        {
          tags|=iv;
        }
      }
      Vector<Box> boxes;
      br.makeBoxes(boxes, tags, pnd, junk, 32/blockingFactor, 1);
      Vector<int> procs;
      LoadBalance(procs, boxes);
      for (int i=0; i< boxes.size(); ++i)
        boxes[i].refine(blockingFactor);
      m_dbls[ires] = new DisjointBoxLayout(boxes, procs);
    }
  }

  // This breaks everything down.
  void tearDown()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
      delete m_dbls[ires];
  }

  // Test the truncation error for the laplacian operator applied to a constant field
  // using Dirichlet boundary conditions.
  // Here, L(phi) = div grad(1) = 0.
  void testTruncErrorForLaplacianOfConstantWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> zero(new ConstantFunction(0.0, dx)),
                                  one(new ConstantFunction(1.0, dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(one));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *one, -1.0, *F0, BC,
                                  *one, *zero, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the laplacian operator applied to a constant field
  // using Neumann boundary conditions.
  // Here, L(phi) = div grad(1) = 0.
  void testTruncErrorForLaplacianOfConstantWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> zero(new ConstantFunction(0.0, dx)),
                                  one(new ConstantFunction(1.0, dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(one));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *one, -1.0, *F0, BC,
                                  *one, *zero, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the laplacian operator applied to a linear
  // field using Dirichlet boundary conditions.
  // Here, L(phi) = div grad(x + y + z) = 0.
  void testTruncErrorForLaplacianOfLinearFieldWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> zero(new ConstantFunction(0.0, dx)),
                                  linearField(new LinearFunction(RealVect::Unit, dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(linearField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *zero, -1.0, *F0, BC,
                                  *linearField, *zero, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the laplacian operator applied to a linear field
  // using Neumann boundary conditions.
  // Here, L(phi) = div grad(x + y + z) = 0.
  void testTruncErrorForLaplacianOfLinearFieldWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> zero(new ConstantFunction(0.0, dx)),
                                  linearField(new LinearFunction(RealVect::Unit, dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(linearField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *zero, -1.0, *F0, BC,
                                  *linearField, *zero, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the laplacian operator applied to a quadratic
  // field using Dirichlet boundary conditions.
  // Here, L(phi) = div grad(x^2 + y^2 + z^2) = 2*SpaceDim.
  void testTruncErrorForLaplacianOfQuadraticFieldWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> zero(new ConstantFunction(0.0, dx)),
                                  twoD(new ConstantFunction(2.0*SpaceDim, dx)),
                                  quadraticField(new SquareDistanceFunction(dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(quadraticField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *zero, -1.0, *F0, BC,
                                  *quadraticField, *twoD, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the laplacian operator applied to a quadratic field
  // using Neumann boundary conditions.
  // Here, L(phi) = div grad(x^2 + y^2 + z^2) = 2*SpaceDim.
  void testTruncErrorForLaplacianOfQuadraticFieldWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> zero(new ConstantFunction(0.0, dx)),
                                  twoD(new ConstantFunction(2.0*SpaceDim, dx)),
                                  quadraticField(new SquareDistanceFunction(dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(quadraticField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *zero, -1.0, *F0, BC,
                                  *quadraticField, *twoD, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the Helmholtz operator (I + div grad) applied
  // to a constant field using Dirichlet boundary conditions.
  // Here, L(phi) = (I + div grad)1 = 1.
  void testTruncErrorForHelmholtzOfConstantWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(one));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(1.0, *one, -1.0, *F0, BC,
                                  *one, *one, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the Helmholtz operator (I + div grad) applied
  // to a constant field using Neumann boundary conditions.
  // Here, L(phi) = (I + div grad)1 = 1.
  void testTruncErrorForHelmholtzOfConstantWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(one));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(1.0, *one, -1.0, *F0, BC,
                                  *one, *one, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the Helmholtz operator (I + div grad) applied
  // to a linear field using Dirichlet boundary conditions.
  // Here, L(phi) = (I + div grad)(x + y + z) = phi.
  void testTruncErrorForHelmholtzOfLinearFieldWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  linearField(new LinearFunction(RealVect::Unit, dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(linearField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(1.0, *one, -1.0, *F0, BC,
                                  *linearField, *linearField, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the Helmholtz operator (I + div grad) applied
  // to a constant field using Neumann boundary conditions.
  // Here, L(phi) = (I + div grad)(x + y + z) = phi.
  void testTruncErrorForHelmholtzOfLinearFieldWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  linearField(new LinearFunction(RealVect::Unit, dx));
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(linearField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(1.0, *one, -1.0, *F0, BC,
                                  *linearField, *linearField, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the Helmholtz operator (I + div grad) applied
  // to a quadratic field using Dirichlet boundary conditions.
  // Here, L(phi) = (I + div grad)(x^2 + y^2 + z^2) = phi + 2D.
  void testTruncErrorForHelmholtzOfQuadraticFieldWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  quadraticField(new SquareDistanceFunction(dx));

      // Compute the actual value of L(phi).
      LevelData<FArrayBox> answer(*m_dbls[ires], 1);
      answer.apply(*quadraticField);
      for (DataIterator dit = m_dbls[ires]->dataIterator(); dit.ok(); ++dit)
        answer[dit()] += 2.0 * SpaceDim;

      // Go to.
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(quadraticField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(1.0, *one, -1.0, *F0, BC,
                                  *quadraticField, answer, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the Helmholtz operator (I + div grad) applied
  // to a quadratic field using Dirichlet boundary conditions.
  // Here, L(phi) = (I + div grad)(x^2 + y^2 + z^2) = phi + 2D.
  void testTruncErrorForHelmholtzOfQuadraticFieldWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  quadraticField(new SquareDistanceFunction(dx));

      // Compute the actual value of L(phi).
      LevelData<FArrayBox> answer(*m_dbls[ires], 1);
      answer.apply(*quadraticField);
      for (DataIterator dit = m_dbls[ires]->dataIterator(); dit.ok(); ++dit)
        answer[dit()] += 2.0 * SpaceDim;

      // Go to.
      RefCountedPtr<FluxCoefficient> F0(new ConstantFluxCoefficient(1.0, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(quadraticField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(1.0, *one, -1.0, *F0, BC,
                                  *quadraticField, answer, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the "linear conductivity" operator div((x + y + z) grad)
  // applied to a constant field using Dirichlet boundary conditions.
  // Here, L(phi) = div((x + y + z)1 = 0.
  void testTruncErrorForLinearConductivityForConstantWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  zero(new ConstantFunction(0.0, dx));
      RefCountedPtr<FluxCoefficient> F(new LinearFluxCoefficient(RealVect::Unit, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(one));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *one, -1.0, *F, BC,
                                  *one, *zero, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the "linear conductivity" operator div((x + y + z) grad)
  // applied to a constant field using Neumann boundary conditions.
  // Here, L(phi) = div((x + y + z)1 = 0.
  void testTruncErrorForLinearConductivityForConstantWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  zero(new ConstantFunction(0.0, dx));
      RefCountedPtr<FluxCoefficient> F(new LinearFluxCoefficient(RealVect::Unit, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(one));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *one, -1.0, *F, BC,
                                  *one, *zero, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the "linear conductivity" operator div((x + y + z) grad)
  // applied to a linear field using Neumann boundary conditions.
  // Here, L(phi) = div[(x + y + z) grad(x + y + z)] = div([x, y, z]) = D.
  void testTruncErrorForLinearConductivityForLinearFieldWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  linearField(new LinearFunction(RealVect::Unit, dx)),
                                  D(new ConstantFunction(SpaceDim, dx));
      RefCountedPtr<FluxCoefficient> F(new LinearFluxCoefficient(RealVect::Unit, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(linearField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *one, -1.0, *F, BC,
                                  *linearField, *D, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the "linear conductivity" operator div((x + y + z) grad)
  // applied to a linear field using Neumann boundary conditions.
  // Here, L(phi) = div[(x + y + z) grad(x + y + z)] = div(x, y, z) = D.
  void testTruncErrorForLinearConductivityForLinearFieldWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  linearField(new LinearFunction(RealVect::Unit, dx)),
                                  D(new ConstantFunction(SpaceDim, dx));
      RefCountedPtr<FluxCoefficient> F(new LinearFluxCoefficient(RealVect::Unit, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(linearField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *one, -1.0, *F, BC,
                                  *linearField, *D, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the "linear conductivity" operator div((x + y + z) grad)
  // applied to a quadratic field using Dirichlet boundary conditions.
  // Here, L(phi) = div[(x + y + z) grad(x^2 + y^2 + z^2)] =
  //   div([2x(x + y + z), 2y(x + y + z), 2z(x + y + z)]) =
  //   4x + 2y + 2z + 2x + 4y + 2z + 2x + 2y + 4z =
  //   2 (1 + D) (x + y + z).
  void testTruncErrorForLinearConductivityForQuadraticFieldWithDirichletBCs()
  {
    Real two1PlusD = 2.0 * (1 + SpaceDim);
    RealVect answerCoeffs(D_DECL(two1PlusD, two1PlusD, two1PlusD));
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  quadraticField(new SquareDistanceFunction(dx)),
                                  answer(new LinearFunction(answerCoeffs, dx));
      RefCountedPtr<FluxCoefficient> F(new LinearFluxCoefficient(RealVect::Unit, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticDirichletBC(quadraticField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *one, -1.0, *F, BC,
                                  *quadraticField, *answer, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the "linear conductivity" operator div((x + y + z) grad)
  // applied to a quadratic field using Neumann boundary conditions.
  // Here, L(phi) = div[(x + y + z) grad(x^2 + y^2 + z^2)] =
  //   div([2x(x + y + z), 2y(x + y + z), 2z(x + y + z)]) =
  //   4x + 2y + 2z + 2x + 4y + 2z + 2x + 2y + 4z =
  //   2 (1 + D) (x + y + z).
  void testTruncErrorForLinearConductivityForQuadraticFieldWithNeumannBCs()
  {
    Real two1PlusD = 2.0 * (1 + SpaceDim);
    RealVect answerCoeffs(D_DECL(two1PlusD, two1PlusD, two1PlusD));
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  quadraticField(new SquareDistanceFunction(dx)),
                                  answer(new LinearFunction(answerCoeffs, dx));
      RefCountedPtr<FluxCoefficient> F(new LinearFluxCoefficient(RealVect::Unit, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(quadraticField));
      Real L1Norm, L2Norm, LInfNorm;
      this->doTruncationErrorTest(0.0, *one, -1.0, *F, BC,
                                  *quadraticField, *answer, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norm, L2Norm, LInfNorm);
      // This should be represented to machine precision.
      TS_ASSERT(L1Norm <= 1e-15);
      TS_ASSERT(L2Norm <= 1e-15);
      TS_ASSERT(LInfNorm <= 1e-15);
    }
  }

  // Test the truncation error for the "linear conductivity" operator div((x + y + z) grad)
  // applied to a quadratic field using Neumann boundary conditions.
  // Here, L(phi) = div[(x + y + z) grad(x^2 + y^2 + z^2)] =
  //   div([2x(x + y + z), 2y(x + y + z), 2z(x + y + z)]) =
  //   4x + 2y + 2z + 2x + 4y + 2z + 2x + 2y + 4z =
  //   2 (1 + D) (x + y + z).
  // This result can't be exactly reproduced by a second-order-accurate
  // operator, so we have to study the convergence of the truncation error.
  void ftestTruncErrorForLinearConductivityForQuadraticFieldWithNeumannBCs()
  {
    Real two1PlusD = 2.0 * (1 + SpaceDim);
    RealVect answerCoeffs(D_DECL(two1PlusD, two1PlusD, two1PlusD));
    Real L1Norms[NUM_RESOLUTIONS], L2Norms[NUM_RESOLUTIONS],
         LInfNorms[NUM_RESOLUTIONS];
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real dx = m_dxs[ires];
      RefCountedPtr<AnalyticForm> one(new ConstantFunction(1.0, dx)),
                                  quadraticField(new SquareDistanceFunction(dx)),
                                  answer(new LinearFunction(answerCoeffs, dx));
      RefCountedPtr<FluxCoefficient> F(new LinearFluxCoefficient(RealVect::Unit, dx));
      RefCountedPtr<BCFunction> BC(new AnalyticNeumannBC(quadraticField));
      this->doTruncationErrorTest(0.0, *one, -1.0, *F, BC,
                                  *quadraticField, *answer, *m_dbls[ires], m_problemDomains[ires], dx,
                                  L1Norms[ires], L2Norms[ires], LInfNorms[ires]);
    }

    // Obtain the rates of convergence from least-squares fits.
    Real L1Rate, L1ErrorBar;
    this->computeConvergenceRate(L1Norms, L1Rate, L1ErrorBar);
    TS_ASSERT(L1Rate >= 2.0);

    Real L2Rate, L2ErrorBar;
    this->computeConvergenceRate(L2Norms, L2Rate, L2ErrorBar);
    TS_ASSERT(L2Rate >= 1.9);

    Real LInfRate, LInfErrorBar;
    this->computeConvergenceRate(L2Norms, LInfRate, LInfErrorBar);
    TS_ASSERT(LInfRate >= 1.9);
  }

  private:

  // ------------------
  //  Helper functions
  // ------------------

  // This performs a generic truncation error test for a variable-coefficient
  // "Poisson" operator L(phi) = alpha A phi + beta div (B grad phi) with
  // the given analytic forms for A, B, and phi, and the expected value of
  // L(phi). Tolerances for the L1, L2, and Linf norms can also be specified.
  // You can also set a_plot to true if you want plots of phi and L(phi).
  void doTruncationErrorTest(Real a_alpha,
                             const LevelData<FArrayBox>::ApplyFunctor& a_Afunc,
                             Real a_beta,
                             const LevelData<FluxBox>::ApplyFunctor& a_Bfunc,
                             BCHolder a_BC,
                             const LevelData<FArrayBox>::ApplyFunctor& a_phi,
                             const LevelData<FArrayBox>& a_Lphi,
                             const DisjointBoxLayout& a_dbl,
                             const ProblemDomain& a_problemDomain,
                             Real a_dx,
                             Real& a_L1Norm,
                             Real& a_L2Norm,
                             Real& a_LInfNorm,
                             bool a_plot = false)
  {
    DataIterator dit(a_dbl);
    LevelData<FArrayBox> phi(a_dbl, 1, IntVect::Unit);
    LevelData<FArrayBox> lofphi(a_dbl, 1), error(a_dbl, 1);

    // Set up phi and the analytic form for L(phi).
    phi.apply(a_phi);

    // Set up A and B coefficients.
    RefCountedPtr<LevelData<FArrayBox> > A(new LevelData<FArrayBox>(a_dbl, 1));
    A->apply(a_Afunc);
    RefCountedPtr<LevelData<FluxBox> > B(new LevelData<FluxBox>(a_dbl, 1));
    B->apply(a_Bfunc);

    // Get an appropriate operator and apply it to phi.
    VCAMRPoissonOp2 amrop;
    amrop.define(a_dbl, a_dx, a_problemDomain, a_BC);
    amrop.setCoefs(A, B, a_alpha, a_beta);
    amrop.applyOp(lofphi, phi);

#ifdef CH_USE_HDF5
    if (a_plot)
    {
      // Plot stuff if we've been asked nicely.
      writeLevelname(&phi, "phi.hdf5");
      writeLevelname(&lofphi, "Lphi.hdf5");
      writeLevelname(&a_Lphi, "answer.hdf5");
    }
#endif

    // Compare the answers and compute the error norms.
    Real volume = D_TERM(a_dx, * a_dx, * a_dx);
    a_LInfNorm = 0.0;
    a_L1Norm = 0.0;
    a_L2Norm = 0.0;
    for (dit.begin(); dit.ok(); ++dit)
    {
      error[dit].axby(a_Lphi[dit()], lofphi[dit()], 1.0, -1.0);
      a_LInfNorm = Max(a_LInfNorm, error[dit].norm(0));
      a_L1Norm += error[dit].norm(1);
      Real L2 = error[dit].norm(2);
      a_L2Norm += L2*L2;
    }

    // Post-process the norms.
    a_L1Norm *= volume;
    a_L2Norm = volume * sqrt(a_L2Norm);

#ifdef CH_USE_HDF5
    // Dump the error if we're plotting things.
    if (a_plot)
      writeLevelname(&error, "error.hdf5");
#endif
  }

  // This performs a generic truncation error test for a variable-coefficient
  // "Poisson" operator L(phi) = alpha A phi + beta div (B grad phi) with
  // the given analytic forms for A, B, and phi, and the expected value of
  // L(phi). Tolerances for the L1, L2, and Linf norms can also be specified.
  // You can also set a_plot to true if you want plots of phi and L(phi).
  void doTruncationErrorTest(Real a_alpha,
                             const LevelData<FArrayBox>::ApplyFunctor& a_Afunc,
                             Real a_beta,
                             const LevelData<FluxBox>::ApplyFunctor& a_Bfunc,
                             BCHolder a_BC,
                             const LevelData<FArrayBox>::ApplyFunctor& a_phi,
                             const LevelData<FArrayBox>::ApplyFunctor& a_Lphi,
                             const DisjointBoxLayout& a_dbl,
                             const ProblemDomain& a_problemDomain,
                             Real a_dx,
                             Real& a_L1Norm,
                             Real& a_L2Norm,
                             Real& a_LInfNorm,
                             bool a_plot = false)
  {
    // Set up the analytic form for L(phi).
    LevelData<FArrayBox> ans(a_dbl, 1);
    ans.apply(a_Lphi);

    // Call the other version.
    this->doTruncationErrorTest(a_alpha, a_Afunc, a_beta, a_Bfunc, a_BC,
                                a_phi, ans, a_dbl, a_problemDomain, a_dx,
                                a_L1Norm, a_L2Norm, a_LInfNorm, a_plot);
  }

  // This method computes the convergence rate of a norm given the resolutions
  // at which its values are computed.
  void computeConvergenceRate(Real a_norm[],
                              Real& a_rate,
                              Real& a_errorBar)
  {
    Real Xs[NUM_RESOLUTIONS], Ys[NUM_RESOLUTIONS];
    Real sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
    Real epsilon = 1e-15;
    Real h0 = 1.0 / (2 * resolutions[0]);
    for (size_t ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Real h = 1.0 / (2 * resolutions[ires]);
      Xs[ires] = log(Abs(h/h0));
      Ys[ires] = log(Abs((a_norm[ires] + epsilon)/(a_norm[0] + epsilon)));
      sumX += Xs[ires];
      sumY += Ys[ires];
      sumXY += Xs[ires]*Ys[ires];
      sumXX += Xs[ires]*Xs[ires];
    }

    // Compute the coefficients for Y = A + B*X.
    Real n = 1.0 * NUM_RESOLUTIONS;
    Real Delta = n*sumXX - sumX*sumX;
    Real A = (sumXX*sumY - sumX*sumXY) / Delta;
    Real B = (n*sumXY - sumX*sumY) / Delta;

    // Compute the square error in the fit.
    Real err2 = 0.0;
    for (size_t ires = 0; ires < NUM_RESOLUTIONS; ++ires)
      err2 += pow(Ys[ires] - A - B*Xs[ires], 2.0);
    err2 /= (n - 2.0);

    // Compute the error in the B coefficient.
    Real errB = sqrt(n * err2 / Delta);

    // The convergence rate is B, and its uncertainty is the error.
    a_rate = B;
    a_errorBar = errB;
  }

  // ------
  //  Data
  // ------

  Box m_domains[NUM_RESOLUTIONS];
  ProblemDomain m_problemDomains[NUM_RESOLUTIONS];
  DisjointBoxLayout* m_dbls[NUM_RESOLUTIONS];
  Real m_dxs[NUM_RESOLUTIONS];
};

#endif
