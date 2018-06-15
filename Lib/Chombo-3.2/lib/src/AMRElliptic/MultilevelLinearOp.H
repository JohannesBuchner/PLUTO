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

// We can assume that template class T has null construction.

#ifndef _MULTILEVELLINEAROP_H_
#define _MULTILEVELLINEAROP_H_

#include <cmath>

#include "RefCountedPtr.H"
#include "Vector.H"
#include "RealVect.H"
#include "LevelData.H"
#include "LinearSolver.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"

#include "NamespaceHeader.H"

///
/**
   Operator class for Linear solver for solving L(phi) = rhs for the case
   where <T> is a multilevel (AMR) data structure; as such, T is actually
   a Vector<T>.  The MultilevelLinearOp is also
   defined with a Vector<AMRLevelOp> which applies the operator in an AMR
   way.  For now, this is hardwired for T = LevelData<FArrayBox>, but we can
   hopefully get back redo this once I get the kinks worked out...
 */
template <class T>
class MultilevelLinearOp : public LinearOp< Vector<LevelData<T>* > >
{

public:
  /// default constructor (sets defaults for preconditioner)
  MultilevelLinearOp();

  /// define function -- opFactory should be able to define all required levels
  /**
  */
  virtual void define(const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_refRatios,
                      const Vector<ProblemDomain>& a_domains,
                      const Vector<RealVect>& a_vectDx,
                      RefCountedPtr<AMRLevelOpFactory<LevelData<T> > >& a_opFactory,
                      int a_lBase);

  ///
  virtual ~MultilevelLinearOp();

  /// compute residual using AMRLevelOps AMRResidual functionality
  /**
     a_lhs = L(a_phi) - a_rhs.
     If a_homogeneous is true,  evaluate the operator using
     homogeneous boundary conditions.
   */
  virtual void residual(  Vector<LevelData<T>* >& a_lhs,
                          const Vector<LevelData<T>* >& a_phi,
                          const Vector<LevelData<T>* >& a_rhs,
                          bool a_homogeneous = false);

  /// Apply preconditioner
  /**
     Given the current state of the residual the correction,
     apply your preconditioner to a_cor.
   */
  virtual void preCond(Vector<LevelData<T>* >& a_cor,
                       const Vector<LevelData<T>* >& a_residual);

  ///
  /**
     In the context of solving L(phi) = rhs, set a_lhs = L(a_phi).
     If a_homogeneous is true,
     evaluate the operator using homogeneous boundary conditions.
   */
  virtual void applyOp(Vector<LevelData<T>* >& a_lhs,
                       const Vector<LevelData<T>* >& a_phi,
                       bool a_homogeneous = false);

  ///
  /**
     Create data holder a_lhs that mirrors a_rhs.   You do not need
     to copy the data of a_rhs, just  make a holder the same size.
  */
  virtual void create(Vector<LevelData<T>* >& a_lhs,
                      const Vector<LevelData<T>* >& a_rhs);

  ///
  /**
     Clean up data holder before it goes out of scope.  This is necessary
     because create calls new.
  */
  virtual void clear(Vector<LevelData<T>* >& a_lhs);

  ///
  /**
     Set a_lhs  equal to a_rhs.
   */
  virtual void assign(Vector<LevelData<T>* >& a_lhs,
                      const Vector<LevelData<T>* >& a_rhs);

  ///
  /**
     Compute and return the dot product of a_1 and a_2.   In most
     contexts, this means return the sum over all data points of a_1*a_2.
   */
  virtual Real dotProduct(const Vector<LevelData<T>* >& a_1,
                          const Vector<LevelData<T>* >& a_2);

  ///
  /**
     Increment by scaled amount (a_lhs += a_scale*a_x).
   */
  virtual void incr  (Vector<LevelData<T>* >& a_lhs,
                      const Vector<LevelData<T>* >& a_x,
                      Real a_scale);

  ///
  /**
     Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
   */
  virtual void axby(Vector<LevelData<T>* >& a_lhs,
                    const Vector<LevelData<T>* >& a_x,
                    const Vector<LevelData<T>* >& a_y,
                    Real a_a,
                    Real a_b);

  ///
  /**
     Multiply the input by a given scale (a_lhs *= a_scale).
   */
  virtual void scale(Vector<LevelData<T>* >& a_lhs,
                     const Real& a_scale);

  ///
  /**
     Return the norm of  a_rhs.
     a_ord == 0  max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
   */
  virtual Real norm(const Vector<LevelData<T>* >& a_rhs,
                    int a_ord);

  ///
  /**
     Set a_lhs to zero.
   */
  virtual void setToZero(Vector<LevelData<T>* >& a_lhs);

  virtual void write(const Vector<LevelData<T>* >* a_data,
                     const char*                   a_filename);

  /// if true, preCond(...) function uses multigrid v-cycle(s)
  bool m_use_multigrid_preconditioner;

  /// number of mg v-cycles for each preCond call (if using mg preconditioner)
  int m_num_mg_iterations;

  /// number of smoothings in the multigrid preconditioner
  int m_num_mg_smooth;

  /// max depth for multigrid preconditioner (if using mg preconditioner)
  int m_preCondSolverDepth;

  ///
  Vector<DisjointBoxLayout> m_vectGrids;

  ///
  Vector<int> m_refRatios;

  ///
  Vector<ProblemDomain> m_domains;

  // cell spacing
  Vector<RealVect> m_vectDx;

  /// vector of level operators
  Vector<RefCountedPtr<AMRLevelOp<LevelData<T> > > > m_vectOperators;

  ///
  int m_lBase;

  /// Solver which is used by the preCond function
  /**
     if m_use_multigrid_preconditioner is true, this solver is
     used to apply a AMR multigrid v-cycle as the preconditioner
  */
  AMRMultiGrid<LevelData<T> > m_preCondSolver;

  /// bottom solver for m_preCondSolver
  LinearSolver<LevelData<T> >* m_precondBottomSolverPtr;

  /// has the preconditioner solver been initialized?
  bool m_isPrecondSolverInitialized;
};

#include "NamespaceFooter.H"

#include "MultilevelLinearOpI.H"
#endif /*_MULTILEVELLINEAROP_H_*/
