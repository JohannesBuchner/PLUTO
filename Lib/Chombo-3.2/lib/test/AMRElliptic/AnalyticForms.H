#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ANALYTICFORMS_H_
#define _ANALYTICFORMS_H_

#include "LevelData.H"
#include "BCFunc.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "BoxIterator.H"
#include "SPACE.H"
#include "UsingNamespace.H"
#include <iostream>
using namespace std;

// Here's a representation of an analytic form, f(x), that can be used to
// set up solutions and to evaluate Dirichlet and Neumann boundary conditions.
struct AnalyticForm: public LevelData<FArrayBox>::ApplyFunctor
{
  explicit AnalyticForm(Real a_dx)
    :
    m_dx(a_dx)
  {
  }

  // Gradient function: to be overridden.
  virtual void grad(const Box& a_box,
                    int a_comps,
                    FluxBox& a_FB) const = 0;

  Real m_dx;

  private:

  AnalyticForm();
  AnalyticForm(const AnalyticForm&);
  AnalyticForm& operator=(const AnalyticForm&);
};

// Here's a representation of a flux coefficient F(x).
struct FluxCoefficient: public LevelData<FluxBox>::ApplyFunctor
{
  explicit FluxCoefficient(Real a_dx)
    :
    m_dx(a_dx)
  {
  }

  Real m_dx;

  private:

  FluxCoefficient();
  FluxCoefficient(const FluxCoefficient&);
  FluxCoefficient& operator=(const FluxCoefficient&);
};

// Here's a representation of a Dirichlet boundary condition that assumes the
// given form on the boundary.
class AnalyticDirichletBC: public BCFunction
{
  public:

  AnalyticDirichletBC(RefCountedPtr<AnalyticForm>& a_form)
    :
    m_form(a_form)
  {
  }

  void operator()(FArrayBox& a_state,
                  const Box& a_valid,
                  const ProblemDomain& a_domain,
                  Real a_dx,
                  bool a_homogeneous)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_domain.isPeriodic(idir))
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();

          Box toRegion = adjCellBox(a_valid, idir, side, 1);
          toRegion &= a_state.box();

          // Evaluate the function on this box.
          FArrayBox f(toRegion, a_state.nComp());
          (*m_form)(toRegion, f.nComp(), f);

          for (BoxIterator bit(toRegion); bit.ok(); ++bit)
          {
            const IntVect& ivTo = bit();
            for (int icomp = 0; icomp < a_state.nComp(); icomp++)
            {
              Real ghostVal = (a_homogeneous) ? -a_state(ivTo, icomp) : f(ivTo, icomp);
              a_state(ivTo, icomp) = ghostVal;
            }
          }
        }
      }
    }
  }

  RefCountedPtr<AnalyticForm> m_form;
};

// Here's a representation of a Neumann boundary condition that assumes the
// given form on the boundary.
class AnalyticNeumannBC: public BCFunction
{
  public:

  AnalyticNeumannBC(RefCountedPtr<AnalyticForm>& a_form)
    :
    m_form(a_form)
  {
  }

  void operator()(FArrayBox& a_state,
                  const Box& a_valid,
                  const ProblemDomain& a_domain,
                  Real a_dx,
                  bool a_homogeneous)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_domain.isPeriodic(idir))
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();

          int isign = sign(side);

          Box toRegion = adjCellBox(a_valid, idir, side, 1);
          toRegion &= a_state.box();

          // Evaluate the function on this box.
          FluxBox gradF(toRegion, a_state.nComp());
          gradF.setVal(0.0);
          if (!a_homogeneous)
            m_form->grad(toRegion, gradF.nComp(), gradF);

          Box fromRegion = toRegion;
          fromRegion.shift(idir, -isign);
          a_state.copy(a_state, fromRegion, 0, toRegion, 0, a_state.nComp());

          for (BoxIterator bit(toRegion); bit.ok(); ++bit)
          {
            const IntVect& ivTo = bit();
            IntVect ivF = ivTo;
            if (isign < 0) ivF[idir] += 1; // Handle cell/face indexing gymnastics.

            for (int icomp = 0; icomp < a_state.nComp(); icomp++)
              a_state(ivTo, icomp) += Real(isign)*a_dx*gradF[idir](ivF, icomp);
            RealVect xF(D_DECL(a_dx * ivF[0], a_dx * ivF[1], a_dx * ivF[2]));
            RealVect xV(D_DECL(a_dx * (ivTo[0] + 0.5), a_dx * (ivTo[1] + 0.5), a_dx * (ivTo[2] + 0.5)));
          }
        }
      }
    }
  }

  RefCountedPtr<AnalyticForm> m_form;
};

// This is a representation of a constant function with a given value F0.
struct ConstantFunction: public AnalyticForm
{
  explicit ConstantFunction(Real a_F0, Real a_dx):
    AnalyticForm(a_dx),
    m_F0(a_F0)
  {
  }

  void operator()(const Box& a_box,
                  int a_comps,
                  FArrayBox& a_FAB) const
  {
    for (int icomp = 0; icomp < a_comps; ++icomp)
      a_FAB.setVal(m_F0, a_box, icomp);
  }

  void grad(const Box& a_box,
            int a_comps,
            FluxBox& a_FB) const
  {
    for (int idir = 0; idir < SpaceDim; ++idir)
      a_FB.setVal(0.0, a_box, idir, 0, a_comps);
  }

  Real m_F0;
};

// This is a representation of a constant function with a given value F0.
struct LinearFunction: public AnalyticForm
{
  explicit LinearFunction(RealVect a_coeffs, Real a_dx):
    AnalyticForm(a_dx),
    m_coeffs(a_coeffs)
  {
  }

  void operator()(const Box& a_box,
                  int a_comps,
                  FArrayBox& a_FAB) const
  {
    RealVect pos;
    ForAllXBNN(Real, a_FAB, a_box, 0, a_comps)
    {
      D_TERM(pos[0]=m_dx*(iR+0.5);, pos[1]=m_dx*(jR+0.5);, pos[2]=m_dx*(kR+0.5));
      D_TERM(a_FABR = m_coeffs[0]*pos[0], + m_coeffs[1]*pos[1], + m_coeffs[2]*pos[2]);
    }EndFor;
  }

  void grad(const Box& a_box,
            int a_comps,
            FluxBox& a_FB) const
  {
    for (int idir = 0; idir < SpaceDim; ++idir)
      a_FB.setVal(m_coeffs[idir], a_box, idir, 0, a_comps);
  }

  RealVect m_coeffs;
};

// This is a representation of the (quadratic) function that returns the square of
// the distance of the given point from the origin.
struct SquareDistanceFunction: public AnalyticForm
{
  explicit SquareDistanceFunction(Real a_dx):
    AnalyticForm(a_dx)
  {
  }

  void operator()(const Box& a_box,
                  int a_comps,
                  FArrayBox& a_FAB) const
  {
    RealVect pos;
    ForAllXBNN(Real, a_FAB, a_box, 0, a_comps)
    {
      D_TERM(pos[0]=m_dx*(iR+0.5);, pos[1]=m_dx*(jR+0.5);, pos[2]=m_dx*(kR+0.5));
      D_TERM(a_FABR = pos[0]*pos[0], + pos[1]*pos[1], + pos[2]*pos[2]);
    }EndFor;
  }

  void grad(const Box& a_box,
            int a_comps,
            FluxBox& a_FB) const
  {
    RealVect pos;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        for (BoxIterator bit(a_box); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          if (sit() == Side::Hi) iv[idir] += 1; // Handle cell/face indexing gymnastics.
          Real xi = m_dx*iv[idir];
          for (int icomp = 0; icomp < a_FB.nComp(); icomp++)
            a_FB[idir](iv, icomp) = 2.0 * xi;
        }
      }
    }
  }
};

// This is a representation of a constant flux function with a given value F0.
struct ConstantFluxCoefficient: public FluxCoefficient
{
  explicit ConstantFluxCoefficient(Real a_F0, Real a_dx):
    FluxCoefficient(a_dx),
    m_F0(a_F0)
  {
  }

  void operator()(const Box& a_box,
                  int a_comps,
                  FluxBox& a_FB) const
  {
    for (int idir = 0; idir < SpaceDim; ++idir)
      a_FB.setVal(m_F0, a_box, idir, 0, a_comps);
  }

  Real m_F0;
};

// This is a representation of a constant function with a given value F0.
struct LinearFluxCoefficient: public FluxCoefficient
{
  explicit LinearFluxCoefficient(RealVect a_coeffs, Real a_dx):
    FluxCoefficient(a_dx),
    m_coeffs(a_coeffs)
  {
  }

  void operator()(const Box& a_box,
                  int a_comps,
                  FluxBox& a_FB) const
  {
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        for (BoxIterator bit(a_box); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          if (sit() == Side::Hi) iv[idir] += 1; // Handle cell/face indexing gymnastics.
          // Remember that the center position of the face is aligned with the cell in the transverse coordinates!
//          RealVect I(D_DECL(iv[0], iv[1], iv[2]));
          RealVect I(D_DECL((idir == 0) ? iv[0] : (iv[0] + 0.5), (idir == 1) ? iv[1] : (iv[1] + 0.5), (idir == 2) ? iv[2] : (iv[2] + 0.5)));
          Real term = D_TERM(m_coeffs[0]*m_dx*I[0], + m_coeffs[1]*m_dx*I[1], + m_coeffs[2]*m_dx*I[2]);
          for (int icomp = 0; icomp < a_FB.nComp(); icomp++)
            a_FB[idir](iv, icomp) = term;
        }
      }
    }
  }

  RealVect m_coeffs;
};

#endif
