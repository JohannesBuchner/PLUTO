#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _OUTPUTAMRMULTIGRIDINSPECTOR_H_
#define _OUTPUTAMRMULTIGRIDINSPECTOR_H_

#include "AMRMultiGrid.H"
#include "NamespaceHeader.H"

//! \class OutputAMRMultiGridInspector
//! This subclass of AMRMultiGridInspector simply writes output files
//! containing residuals and corrections at each multigrid iteration.
template <class T>
class OutputAMRMultiGridInspector: public AMRMultiGridInspector<T>
{
  public:

  //! Create an AMRMultiGridInspector that inspects the given AMRMultiGrid solver and
  //! generates output files containing intermediate data.
  //! a_name The name of the inspector. This gets prepended to the output files.
  //! a_solver The AMRMultiGrid solver with which this inspector is associated.
  OutputAMRMultiGridInspector(const std::string& a_name,
                              AMRMultiGrid<T>& a_solver)
    :AMRMultiGridInspector<T>(),
     m_name(a_name),
     m_solver(&a_solver)
  {
  }

  //! Destructor.
  ~OutputAMRMultiGridInspector()
  {
  }

  // Overridden methods.
  void recordResiduals(const Vector<T*>& a_residuals,
                       int a_minLevel,
                       int a_maxLevel,
                       int a_iter)
  {
    char name[1024];
    snprintf(name, 1024, "%s.residual.iter.%d", m_name.c_str(), a_iter);
    std::string nameStr(name);
    m_solver->outputAMR(const_cast<Vector<T*>&>(a_residuals), nameStr,
                        a_maxLevel, a_minLevel);
  }

  void recordCorrections(const Vector<T*>& a_corrections,
                         int a_minLevel,
                         int a_maxLevel,
                         int a_iter)
  {
    char name[1024];
    snprintf(name, 1024, "%s.correction.iter.%d", m_name.c_str(), a_iter);
    std::string nameStr(name);
    m_solver->outputAMR(const_cast<Vector<T*>&>(a_corrections), nameStr,
                        a_maxLevel, a_minLevel);
  }

  private:

  std::string m_name;
  AMRMultiGrid<T>* m_solver;
};

#include "NamespaceFooter.H"
#endif
