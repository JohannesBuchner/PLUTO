#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBCONSBACKWARDEULERINTEGRATOR_H_
#define _EBCONSBACKWARDEULERINTEGRATOR_H_

#include "AMR.H"
#include "AMRLevel.H"
#include "AMRMultiGrid.H"
#include "EBCellFAB.H"
#include "EBConductivityOp.H"
#include "RefCountedPtr.H"
#include "LevelData.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"

#include "NamespaceHeader.H"

// This is a one-off, single-level, conservative backward Euler time
// integrator for the heat equation with time-dependent mass density.
// NO AMR ALLOWED!  This class is not long for this world!
class EBConsBackwardEulerIntegrator
{
  public:

  // Construct this thing for the given operator.
  explicit EBConsBackwardEulerIntegrator(RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >& a_solver);

  // Destructor.
  virtual ~EBConsBackwardEulerIntegrator();

  // Integrate the heat equation to find a new temperature.
  void updateSolution(LevelData<EBCellFAB>& a_diffusionTerm,
                      const LevelData<EBCellFAB>& a_oldTemperature,
                      Real a_specificHeat,
                      const LevelData<EBCellFAB>& a_oldDensity,
                      const LevelData<EBCellFAB>& a_newDensity,
                      const LevelData<EBCellFAB>& a_source,
                      Real a_time,
                      Real a_timeStep,
                      bool a_zeroTemp);

  // Compute the diffusion term by integrating the temperature and then
  // backing out the time derivative of the energy.
  void computeDiffusion(LevelData<EBCellFAB>& a_diffusionTerm,
                        const LevelData<EBCellFAB>& a_oldTemperature,
                        Real a_specificHeat,
                        const LevelData<EBCellFAB>& a_oldDensity,
                        const LevelData<EBCellFAB>& a_newDensity,
                        const LevelData<EBCellFAB>& a_source,
                        Real a_time,
                        Real a_timeStep,
                        bool a_zeroTemp);

  private:

  // The multigrid solver.
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > m_solver;

  // The conductivity operator.
  EBConductivityOp* m_op;

  // The level grid on which this integrator is defined.
  EBLevelGrid m_eblg;

  // No-nos.
  EBConsBackwardEulerIntegrator();
  EBConsBackwardEulerIntegrator(const EBConsBackwardEulerIntegrator&);
  EBConsBackwardEulerIntegrator& operator=(const EBConsBackwardEulerIntegrator&);

};

#include "NamespaceFooter.H"
#endif
