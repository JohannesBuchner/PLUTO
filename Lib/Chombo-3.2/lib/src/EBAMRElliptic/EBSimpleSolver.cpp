#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBSimpleSolver.H"
#include "NamespaceHeader.H"

EBSimpleSolver::EBSimpleSolver()
{
  m_isDefined = false;

  m_operator = NULL;
  m_homogeneous = true;

  m_numSmooths = 4;
}

EBSimpleSolver::~EBSimpleSolver()
{
}

void EBSimpleSolver::setHomogeneous(bool a_homogeneous)
{
  m_homogeneous = a_homogeneous;
}

void EBSimpleSolver::define(LinearOp<LevelData<EBCellFAB> >* a_operator,
                            bool                             a_homogeneous)
{
  m_isDefined = true;

  m_operator = dynamic_cast <MGLevelOp<LevelData<EBCellFAB> >* > (a_operator);
  if (m_operator == NULL)
  {
    MayDay::Error("EBSimpleSolver::define - operator not an MGLevelOp");
  }

  m_homogeneous = a_homogeneous;
}

void EBSimpleSolver::setNumSmooths(const int& a_numSmooths)
{
  m_numSmooths = a_numSmooths;
}

void EBSimpleSolver::solve(LevelData<EBCellFAB>&       a_phi,
                           const LevelData<EBCellFAB>& a_rhs)
{
  CH_assert(m_isDefined);

  m_operator->relax(a_phi,a_rhs,m_numSmooths);
}
#include "NamespaceFooter.H"
