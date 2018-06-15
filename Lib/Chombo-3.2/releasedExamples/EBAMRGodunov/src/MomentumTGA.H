#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MOMENTUMTGA_H_
#define _MOMENTUMTGA_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "EBFluxFAB.H"
#include "EBFluxRegister.H"
#include "ProblemDomain.H"
#include "BaseLevelTGA.H"
#include "EBViscousTensorOp.H"

#include "NamespaceHeader.H"

//! Momentum TGA adds another function to EBLevelTGA
//! that does not make sense in a broader context
class MomentumTGA: public EBLevelTGA
{

public:

  //! Create a new TGA level integrator.
  //! \param a_grids The disjoint box layout on which the level integrator is defined.
  //! \param a_refRat The refinement ratios for the boxes.
  //! \param a_level0Domain The coarsest grid level defining the problem domain.
  //! \param a_opFact A factory that generates level diffusion operators.
  //! \param a_solver A multigrid solver.
  MomentumTGA(const Vector<DisjointBoxLayout>&                               a_grids,
              const Vector<int>&                                             a_refRat,
              const ProblemDomain&                                           a_level0Domain,
              RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&      a_opFact,
              const RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >&     a_solver)
    :EBLevelTGA(a_grids, a_refRat, a_level0Domain, a_opFact,         a_solver)
  {
  }

  //! Destructor
  virtual ~MomentumTGA()
  {
  }

  /// compute volfrac(sigma dot grad U) from inputs
  void getShearStressDotGradU(LevelData<EBCellFAB>&       a_sigmaGradU,
                             const LevelData<EBCellFAB>& a_gradU,
                             int a_level)
  {
    EBViscousTensorOp* ebvto = dynamic_cast<EBViscousTensorOp*>(m_ops[a_level]);
    if (ebvto == NULL)
      {
        MayDay::Error("MomentumTGA::unrecognized operator (i need EBVTO)");
      }

    //flux gets multiplied by beta deep in the bowels here.
    resetSolverAlphaAndBeta(1.0, 1.0);
    ebvto->getShearStressDotGradU(a_sigmaGradU, a_gradU, a_level);
  }

  /// compute volfrac( div sigma u) from inputs
  void getKappaDivSigmaU(LevelData<EBCellFAB>&       a_divSigmaU,
                         const LevelData<EBCellFAB>& a_velocity,
                         const LevelData<EBCellFAB>* a_veloCoar,
                         int a_level)
  {
    EBViscousTensorOp* ebvto = dynamic_cast<EBViscousTensorOp*>(m_ops[a_level]);
    if (ebvto == NULL)
      {
        MayDay::Error("MomentumTGA::unrecognized operator (i need EBVTO)");
      }

    //flux gets multiplied by beta deep in the bowels here.
    resetSolverAlphaAndBeta(1.0, 1.0);
    ebvto->getKappaDivSigmaU(a_divSigmaU, a_velocity, a_veloCoar, a_level);
  }

};

//! MomentumBackwardEuler adds another function to EBLevelBackwardEuler
//! that does not make sense in a broader context
class MomentumBackwardEuler: public EBLevelBackwardEuler
{

public:

  //! Create a new TGA level integrator.
  //! \param a_grids The disjoint box layout on which the level integrator is defined.
  //! \param a_refRat The refinement ratios for the boxes.
  //! \param a_level0Domain The coarsest grid level defining the problem domain.
  //! \param a_opFact A factory that generates level diffusion operators.
  //! \param a_solver A multigrid solver.
  MomentumBackwardEuler(const Vector<DisjointBoxLayout>&                               a_grids,
                        const Vector<int>&                                             a_refRat,
                        const ProblemDomain&                                           a_level0Domain,
                        RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&      a_opFact,
                        const RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >&     a_solver)
    :EBLevelBackwardEuler(a_grids, a_refRat, a_level0Domain, a_opFact,         a_solver)
  {
  }

  //! Destructor
  virtual ~MomentumBackwardEuler()
  {
  }

  /// compute volfrac(sigma dot grad U) from inputs
  void getShearStressDotGradU(LevelData<EBCellFAB>&       a_sigmaGradU,
                             const LevelData<EBCellFAB>& a_gradU,
                             int a_level)
  {
    EBViscousTensorOp* ebvto = dynamic_cast<EBViscousTensorOp*>(m_ops[a_level]);
    if (ebvto == NULL)
      {
        MayDay::Error("MomentumBackwardEuler::unrecognized operator (i need EBVTO)");
      }

    resetSolverAlphaAndBeta(1.0, 1.0);
    ebvto->getShearStressDotGradU(a_sigmaGradU, a_gradU, a_level);
  }

  /// compute volfrac( div sigma u) from inputs
  void getKappaDivSigmaU(LevelData<EBCellFAB>&       a_divSigmaU,
                         const LevelData<EBCellFAB>& a_velocity,
                         const LevelData<EBCellFAB>* a_veloCoar,
                         int a_level)
  {
    EBViscousTensorOp* ebvto = dynamic_cast<EBViscousTensorOp*>(m_ops[a_level]);
    if (ebvto == NULL)
      {
        MayDay::Error("MomentumBackwardEuler::unrecognized operator (i need EBVTO)");
      }

    resetSolverAlphaAndBeta(1.0, 1.0);
    ebvto->getKappaDivSigmaU(a_divSigmaU, a_velocity, a_veloCoar, a_level);
  }

};

#include "NamespaceFooter.H"
#endif
