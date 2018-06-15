#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBSIMPLESOLVER_H_
#define _EBSIMPLESOLVER_H_

#include "LevelData.H"
#include "LinearSolver.H"
#include "MultiGrid.H"

#include "EBCellFAB.H"
#include "NamespaceHeader.H"

class EBSimpleSolver : public LinearSolver<LevelData<EBCellFAB> >
{
public:
  EBSimpleSolver();

  virtual ~EBSimpleSolver();

  virtual void setHomogeneous(bool a_homogeneous);

  virtual void define(LinearOp<LevelData<EBCellFAB> >* a_operator,
                      bool                             a_homogeneous);

  virtual void setNumSmooths(const int& a_numSmooths);

  virtual void solve(LevelData<EBCellFAB>&       a_phi,
                     const LevelData<EBCellFAB>& a_rhs);

private:
  bool m_isDefined;

  MGLevelOp<LevelData<EBCellFAB> >* m_operator;

  bool m_homogeneous;
  int  m_numSmooths;
};

#include "NamespaceFooter.H"
#endif
