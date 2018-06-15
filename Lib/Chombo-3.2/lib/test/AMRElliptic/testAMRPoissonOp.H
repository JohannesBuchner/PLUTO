#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TESTAMRPOISSONOP_H_
#define _TESTAMRPOISSONOP_H_

#include <iostream>
using std::endl;

#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "BoxIterator.H"
#include "FABView.H"

#include "NewPoissonOp.H"
#include "AMRPoissonOp.H"
#include "BCFunc.H"
#include <cxxtest/TestSuite.h>
#include "UsingNamespace.H"

/// Application-specific global variables:
#define NUM_RESOLUTIONS 3
#if CH_SPACEDIM == 2
static int resolutions[NUM_RESOLUTIONS] =
{
  64, 128, 256
};
#else
static int resolutions[NUM_RESOLUTIONS] =
{
  32, 64, 128
};
#endif
static int blockingFactor = 8;
static Real s_currentDx = 1.0;

//   u = x*x + y*y + z*z
//   du/dx = 2*x
//   du/dy = 2*y
//   du/dz = 2*z
//   Laplace(u) = 2*CH_SPACEDIM
extern "C"
{

  void Parabola_neum(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    switch (*dir)
      {
      case 0:
        a_values[0]=2*pos[0];
        return;
      case 1:
        a_values[0]=2*pos[1];
        return;
      case 2:
        a_values[0]=2*pos[2];
        return;
      default:
        MayDay::Error("no such dimension");
      };
  }
  void Parabola_diri(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    a_values[0] = D_TERM(pos[0]*pos[0],+pos[1]*pos[1],+pos[2]*pos[2]);
  }

  //this sets ghost cells in domain or out.
  //this makes zero sense outside of test environment
  void DirParabolaBC(FArrayBox& a_state,
                     const Box& valid,
                     const ProblemDomain& a_domain,
                     Real a_dx,
                     bool a_homogeneous)
  {
    for (int i=0; i<CH_SPACEDIM; ++i)
      {
        for (SideIterator sit; sit.ok(); ++sit)
          {
            DiriBC(a_state,
                   valid,
                   a_dx,
                   a_homogeneous,
                   Parabola_diri,
                   i,
                   sit(),2);
          }
      }
  }

  //this sets ghost cells in domain or out.
  //this makes zero sense outside of test environment
  void NeumParabolaBC(FArrayBox& a_state,
                      const Box& valid,
                      const ProblemDomain& a_domain,
                      Real a_dx,
                      bool a_homogeneous)
  {
    for (int i=0; i<CH_SPACEDIM; ++i)
      {
        for (SideIterator sit; sit.ok(); ++sit)
          {

            NeumBC(a_state,
                   valid,
                   a_dx,
                   a_homogeneous,
                   Parabola_neum,
                   i,
                   sit());
          }
      }
  }
}

static BCValueFunc pointFunc = Parabola_diri;

void parabola(const Box& box, int comps, FArrayBox& t)
{
  RealVect pos;
  int dir;
  Side::LoHiSide side;
  int num = 1;
  Real dx = s_currentDx;
  ForAllXBNN(Real,t, box, 0, comps)
  {
      num=nR;
      D_TERM(pos[0]=dx*(iR+0.5);, pos[1]=dx*(jR+0.5);, pos[2]=dx*(kR+0.5));
      pointFunc(&(pos[0]), &dir, &side, &tR);
  }EndFor;
}

// Test suite class.
class TestAMRPoissonOp: public CxxTest::TestSuite
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

  // Tests the poisson operator applied to a parabola
  // using Dirichlet boundary conditions.
  void testTruncationErrorForParabolaWithDirichletBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Box phiBox = m_domains[ires];
      ProblemDomain regularDomain = m_problemDomains[ires];
      Real dx = m_dxs[ires];
      phiBox.grow(1);

      FArrayBox phi(phiBox, 1);
      FArrayBox rhs(m_domains[ires], 1);
      FArrayBox lofphi(m_domains[ires], 1);

      rhs.setVal(2*CH_SPACEDIM);
      s_currentDx = dx;
      parabola(phiBox, 1, phi);

      RealVect pos(IntVect::Unit);
      pos *= dx;

      NewPoissonOp op;
      op.define(pos, regularDomain, DirParabolaBC);
      op.applyOp(lofphi, phi);

      lofphi -= 2*CH_SPACEDIM;
      Real norm0 = lofphi.norm(0);
      Real norm1 = lofphi.norm(1);
      Real norm2 = lofphi.norm(2);

      TS_ASSERT(norm0 <= 1e-6);
#if 0
    if (norm0 > 1E-6)
    {
      pout() << "applyOp for parabolic dirichlet function seems busted, max(L(phi))=" << norm0
        <<std::endl;
      pout() << "applyOp: 1-norm, 2-norm = " << norm1 << " , " << norm2 << std::endl ;
      return 1;
    }
    else if (verbose)
    {
      pout() << indent << "applyOp for parabolic dirichlet: norms(0,1,2) = "
        << norm0 << " , " << norm1 << " , " << norm2 << std::endl ;
    }
#endif
    }
  }

  // Tests the poisson operator applied to a parabola
  // using Neumann boundary conditions.
  void testTruncationErrorForParabolaWithNeumannBCs()
  {
    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      Box phiBox = m_domains[ires];
      ProblemDomain regularDomain = m_problemDomains[ires];
      Real dx = m_dxs[ires];
      phiBox.grow(1);

      FArrayBox phi(phiBox, 1);
      FArrayBox rhs(m_domains[ires], 1);
      FArrayBox lofphi(m_domains[ires], 1);

      rhs.setVal(2*CH_SPACEDIM);
      s_currentDx = dx;
      parabola(phiBox, 1, phi);

      RealVect pos(IntVect::Unit);
      pos *= dx;

      NewPoissonOp op;
      op.define(pos, regularDomain, NeumParabolaBC);
      op.applyOp(lofphi, phi);

      lofphi -= 2*CH_SPACEDIM;
      Real norm0 = lofphi.norm(0);
      Real norm1 = lofphi.norm(1);
      Real norm2 = lofphi.norm(2);
      TS_ASSERT(norm0 <= 1E-6);

#if 0
    if (norm0 > 1E-6)
    {
      pout()<< "applyOp for parabolic neumann function seems busted, max(L(phi))=" << norm0
        <<std::endl;
      pout() << "applyOp: 1-norm, 2-norm = " << norm1 << " , " << norm2 << std::endl ;
      return 1;
    }
    else if (verbose)
    {
      pout() << indent << "applyOp for parabolic neumann: norms(0,1,2) = "
        << norm0 << " , " << norm1 << " , " << norm2 << std::endl ;
    }
#endif
    }
  }

  void testAMRTruncationErrorForParabolaWithDirichletBCs()
  {
    struct setvalue
    {
      static void setFunc(const Box& box,
                          int comps, FArrayBox& t)
      {
        t.setVal(4.0);
      }
    };

    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      const DisjointBoxLayout& dbl = *m_dbls[ires];
      Real dx = m_dxs[ires];
      ProblemDomain regularDomain = m_problemDomains[ires];

      DataIterator dit(dbl);
      LevelData<FArrayBox> phi(dbl, 1, IntVect::Unit);
      LevelData<FArrayBox> rhs(dbl, 1);
      LevelData<FArrayBox> lofphi(dbl, 1);

      rhs.apply(setvalue::setFunc);
      s_currentDx = dx;
      phi.apply(parabola);

      RealVect pos(IntVect::Unit);
      pos *= dx;

      AMRPoissonOp amrop;
      amrop.define(dbl, pos[0], regularDomain, DirParabolaBC);
      amrop.applyOp(lofphi, phi);

      for (dit.begin(); dit.ok(); ++dit)
      {
        lofphi[dit] -= 2*CH_SPACEDIM;
        Real norm0 = lofphi[dit].norm(0);
        Real norm1 = lofphi[dit].norm(1);
        Real norm2 = lofphi[dit].norm(2);
        TS_ASSERT(norm0 <= 1e-6);
#if 0
      if (norm0 > 1E-6)
        {
          pout()<< "applyOp for AMR parabolic dirichlet function seems busted on Box " << lofphi[dit].box()
                << ", max(L(phi))=" << norm0 <<std::endl;
          pout() << "applyOp: 1-norm, 2-norm = " << norm1 << " , " << norm2 << std::endl ;
          return 1;
        }
      else if (verbose)
      {
        pout() << indent << "applyOp for AMR parabolic dirichlet on " << lofphi[dit].box()
               << ": norms(0,1,2) = " << norm0 << " , " << norm1 << " , " << norm2 << std::endl ;
      }
#endif
      }
    }
  }

  void testAMRTruncationErrorForParabolaWithNeumannBCs()
  {
    struct setvalue
    {
      static void setFunc(const Box& box,
                          int comps, FArrayBox& t)
      {
        t.setVal(4.0);
      }
    };

    for (int ires = 0; ires < NUM_RESOLUTIONS; ++ires)
    {
      const DisjointBoxLayout& dbl = *m_dbls[ires];
      Real dx = m_dxs[ires];
      ProblemDomain regularDomain = m_problemDomains[ires];

      DataIterator dit(dbl);
      LevelData<FArrayBox> phi(dbl, 1, IntVect::Unit);
      LevelData<FArrayBox> rhs(dbl, 1);
      LevelData<FArrayBox> lofphi(dbl, 1);

      rhs.apply(setvalue::setFunc);
      s_currentDx = dx;
      phi.apply(parabola);

      RealVect pos(IntVect::Unit);
      pos *= dx;

      AMRPoissonOp amrop;
      amrop.define(dbl, pos[0], regularDomain, NeumParabolaBC);
      amrop.applyOp(lofphi, phi);

      for (dit.begin(); dit.ok(); ++dit)
      {
        lofphi[dit] -= 2*CH_SPACEDIM;
        Real norm0 = lofphi[dit].norm(0);
        Real norm1 = lofphi[dit].norm(1);
        Real norm2 = lofphi[dit].norm(2);
        TS_ASSERT(norm0 <= 1E-6);
#if 0
      if (norm0 > 1E-6)
        {
          pout()<< "applyOp for AMR parabolic neumann function seems busted on Box " << lofphi[dit].box()
                << ", max(L(phi))=" << norm0 <<std::endl;
          pout() << "applyOp: 1-norm, 2-norm = " << norm1 << " , " << norm2 << std::endl ;
          return 1;
        }
      else if (verbose)
      {
        pout() << indent << "applyOp for AMR parabolic neumann on " << lofphi[dit].box()
               << ": norms(0,1,2) = " << norm0 << " , " << norm1 << " , " << norm2 << std::endl ;
      }
#endif
      }
    }
  }

  // Tests truncation error on the operator applied to a constant.
  void testConstantTruncationError()
  {
  }

  // Tests truncation error on the operator applied to a linearly-varying field.
  void testLinearTruncationError()
  {
  }

  // Tests truncation error on the operator applied to a quadratic field.
  void testQuadraticTruncationError()
  {
  }

  // Tests truncation error on the operator applied to a sinusoid.
  void testSinusoidTruncationError()
  {
  }

  private:

  Box m_domains[NUM_RESOLUTIONS];
  ProblemDomain m_problemDomains[NUM_RESOLUTIONS];
  DisjointBoxLayout* m_dbls[NUM_RESOLUTIONS];
  Real m_dxs[NUM_RESOLUTIONS];
};

#endif
