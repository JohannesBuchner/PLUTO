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

#ifndef _MERGESOLVER_H_
#define _MERGESOLVER_H_

#include "LinearSolver.H"
#include "parstream.H"
#include "BiCGStabSolver.H"
#include "LevelData.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

///
/**
   Elliptic solver using the BiCGStab algorithm that turns the bottom solve into a single grid and solves on one processor if it can.
 */

template <class T>
class MergeSolver : public BiCGStabSolver<LevelData<T> >
{
public:

  MergeSolver();

  virtual ~MergeSolver()
  {
  }

  ///
  /**
     reset whether the solver is homogeneous.
   */
  virtual void setHomogeneous(bool a_homogeneous)
  {
    BiCGStabSolver<LevelData<T> >::setHomogeneous(a_homogeneous);
  }

  virtual void define(LinearOp<LevelData<T> >* a_operator, bool a_homogeneous = false)
  {
    BiCGStabSolver<LevelData<T> >::define(a_operator, a_homogeneous);
  }

  virtual void setConvergenceMetrics(Real a_metric, Real a_tolerance)
  {
    BiCGStabSolver<LevelData<T> >::setConvergenceMetrics(a_metric, a_tolerance);
  }

  /// only overridden virtual function
  virtual void solve(LevelData<T>& a_phi, const LevelData<T>& a_rhs);

protected:

#ifdef CH_MPI
  MPI_Comm m_comm;
#endif
};

template <class T>
MergeSolver<T>::MergeSolver()
{
#ifdef CH_MPI
  MPI_Group proc0, cworld;
  int members[1] =
  {
    0
  };
  int err;
  err = MPI_Comm_group(Chombo_MPI::comm, &cworld);
  err = MPI_Group_incl(cworld , 1, members, &proc0);
  err = MPI_Comm_create(Chombo_MPI::comm, proc0, &m_comm);
#endif
}

template <class T>
void MergeSolver<T>::solve(LevelData<T>& a_phi, const LevelData<T>& a_rhs)
{
  CH_TIME("MergeSolver::solve");
  DisjointBoxLayout dbl = a_phi.disjointBoxLayout();
  Box region;
  int numPts=0;
  for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit)
    {
      region.minBox(dbl[lit]);
      numPts+=dbl[lit].numPts();
    }
  if (region.numPts() == numPts)
    {
      // bottomsolve can be done in just one box without any exchanges.
      Vector<Box> newBoxes(1, region);
      Vector<int> newProc(1,0);
      DisjointBoxLayout newDBL(newBoxes, newProc, dbl.physDomain());
      LevelData<T> phi(newDBL, a_phi.nComp(), a_phi.ghostVect());
      LevelData<T> rhs(newDBL, a_rhs.nComp(), a_rhs.ghostVect());
      a_phi.copyTo(phi);
      a_rhs.copyTo(rhs);
#ifdef CH_MPI
      MPI_Comm save=Chombo_MPI::comm;
      Chombo_MPI::comm = m_comm;
#endif
      if (procID() == 0)
        BiCGStabSolver<LevelData<T> >::solve(phi, rhs);
#ifdef CH_MPI
      Chombo_MPI::comm = save;
#endif
      phi.copyTo(a_phi);
    }
  else
    {
    BiCGStabSolver<LevelData<T> >::solve(a_phi, a_rhs);
  }
}

#include "NamespaceFooter.H"

#endif /*_MERGESOLVER_H_*/
