#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ROOTSOLVER_H_
#define _ROOTSOLVER_H_

/******************************************************************************/
/**
 * \file
 *
 * \brief Root solvers
 *
 *//*+*************************************************************************/

#include <cmath>
#include <algorithm>

#include "CH_assert.H"

#include "BaseNamespaceHeader.H"

namespace RootSolver
{

template <typename T> struct RootTr
{
};

// Default epsilon and tolerance for floats
template <> struct RootTr<float>
{
  enum
  {
    maxIter = 100
  };

  static float eps()
  {
    return 3.0E-7;
  }

  static float tolerance()
  {
    return 1.0e-6;
  }
};

// Default epsilon and tolerance for doubles
template <> struct RootTr<double>
{
  enum
  {
    maxIter = 100
  };

  static double eps()
  {
    return 3.0E-15;
  }

  static double tolerance()
  {
    return 1.0e-12;
  }
};

/*******************************************************************************
 */
///  Brent's root solver
/**
 *   \tparam T          Type for x and f(x) - must be floating point
 *   \tparam Func       Function object describing function to solve where
 *                      T f(x) = Func::operator()(const T& x)
 *
 *   \param[out] numIter
 *                      Number of iterations required for convergence to
 *                      specified tolerance.  If equal to MAXITER, the solution
 *                      is not within specified tolerance.
 *   \param[in]  f      Instance of function to solve
 *   \param[in]  aPt    Lower bound
 *   \param[in]  bPt    Upper bound
 *   \param[in]  tol    Tolerance for solve - essentially the spread of the
 *                      bracket.  This can be specified in absolute terms, or,
 *                      if given by an integer cast to T, the number of
 *                      significant digits to solve for in x.  The default
 *                      is given by the RootTr class.  Note that epsilon is
 *                      also considered for specifying the spread of the
 *                      bracket.
 *   \param[in]  MAXITER
 *                      Maximum iterations.  Default (100).  If reached,
 *                      a message is written to cerr but the program otherwise
 *                      completes
 *   \return            x where f(x) = 0
 *
 *   Example                                                           \verbatim
 *     #include <functional>
 *     #include "RootSolver.H"
 *     // Func is not allowed to be local until the C++0x standard
 *     struct Func : public std::unary_function<Real, Real>
 *     {
 *       Real operator()(const Real& a_x) const
 *         {
 *           return 5*std::pow(a_x, 5) - 3*std::pow(a_x, 3) + a_x;
 *         }
 *     };
 *     void foo()
 *     {
 *       int numIter;
 *       const Real xmin = -1.;
 *       const Real xmax =  1.;
 *       const Real x0 = RootSolver::Brent(numIter, Func(), xmin, xmax);
 *       if (numIter == RootTr<Real>::maxIter)
 *         {
 *           std::cout << "Uh oh\n";
 *         }
 *     }
 *                                                                  \endverbatim
 *//*+*************************************************************************/

template <typename T, typename Func>
T Brent(int&           numIter,
        const Func&    f,
        T              aPt,
        T              bPt,
        T              tol = RootTr<T>::tolerance(),
        const unsigned MAXITER = RootTr<T>::maxIter)
{
  CH_assert(tol >= 0.);
  CH_assert(MAXITER > 0);
  const T eps = RootTr<T>::eps();
  const int prec = (int)tol;
  if (((T)prec) == tol)
    {
      tol = std::pow(10., -std::abs(prec));
    }
  numIter = -1;

  //  Max allowed iterations and floating point precision
  unsigned i;
  T c, d, e;
  T p, q, r, s;

  T fa = f(aPt);
  T fb = f(bPt);

  //  Init these to be safe
  c = d = e = 0.0;

  if (fb*fa > 0)
    {
      MayDay::Abort("RootSolver::Brent: Root must be bracketed");
    }

  T fc = fb;

  for (i = 0; i < MAXITER; i++)
    {
      if (fb*fc > 0)
        {
          //  Rename a, b, c and adjust bounding interval d
          c = aPt;
          fc  = fa;
          d = bPt - aPt;
          e = d;
        }

      if (Abs(fc) < Abs(fb))
        {
          aPt = bPt;
          bPt = c;
          c   = aPt;
          fa  = fb;
          fb  = fc;
          fc  = fa;
        }

      //  Convergence check
      const T tol1  = 2.0 * eps * Abs(bPt) + 0.5 * tol;
      const T xm    = 0.5 * (c - bPt);

      if (Abs(xm) <= tol1 || fb == 0.0)
        {
          break;
        }

      if (Abs(e) >= tol1 && Abs(fa) > Abs(fb))
        {
          //  Attempt inverse quadratic interpolation
          s = fb / fa;
          if (aPt == c)
            {
              p = 2.0 * xm * s;
              q = 1.0 - s;
            }
          else
            {
              q = fa / fc;
              r = fb / fc;
              p = s * (2.0 * xm * q * (q-r) - (bPt-aPt) * (r-1.0));
              q = (q-1.0) * (r-1.0) * (s-1.0);
            }

          //  Check whether in bounds
          if (p > 0) q = -q;

          p = Abs(p);

          if (2.0 * p < std::min(((float)3.0)*xm*q-Abs(tol1*q),
                                 Abs(e*q)))
            {
              //  Accept interpolation
              e = d;
              d = p / q;
            }
          else
            {
              //  Interpolation failed, use bisection
              d = xm;
              e = d;
            }
        }
      else
        {
          //  Bounds decreasing too slowly, use bisection
          d = xm;
          e = d;
        }

      //  Move last best guess to a
      aPt = bPt;
      fa  = fb;

      //  Evaluate new trial root
      if (Abs(d) > tol1)
        {
          bPt = bPt + d;
        }
      else
        {
          if (xm < 0) bPt = bPt - tol1;
          else        bPt = bPt + tol1;
        }

      fb = f(bPt);
    }

  if (i >= MAXITER)
    {
      std::cerr  << "RootSolver::Brent: exceeded maximum iterations: "
                 << MAXITER << std::endl;
    }

  numIter = i;
  return bPt;
}

}  // End of namespace RootSolver

#include "BaseNamespaceFooter.H"

#endif
