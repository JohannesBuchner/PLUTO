#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// BVS, June 18, 2003

#ifndef _NOOPSOLVER_H_
#define _NOOPSOLVER_H_

#include "LinearSolver.H"
#include "NamespaceHeader.H"

template <class T>
class NoOpSolver : public LinearSolver<T>
{
public:
  NoOpSolver()
  {
  }

  virtual ~NoOpSolver()
  {
  }

  virtual void setHomogeneous(bool a_homogeneous)
  {
  }

  virtual void define(LinearOp<T>* a_factory, bool a_homogeneous)
  {
  }

  virtual void solve(T& a_phi, const T& a_rhs);

  virtual void setConvergenceMetrics(Real a_metric,
                                     Real a_tolerance)
  {
  }
};

template <class T>
void NoOpSolver<T>::solve(T& a_phi, const T& a_rhs)
{
  //as the name implies, this function does nothing.
}

#include "NamespaceFooter.H"
#endif /*_NOOPSOLVER_H_*/
