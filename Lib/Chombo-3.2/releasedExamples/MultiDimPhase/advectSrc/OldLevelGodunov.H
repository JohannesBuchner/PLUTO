#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _OLDLEVELGODUNOV_H_
#define _OLDLEVELGODUNOV_H_

// Original design:                     Phil Colella
// Original implementation:             Brian Van Straalen
// Implementation date:                 April 12, 1999
// Rearchitected for unsplit algorithm: Phil Colella/Terry Ligocki
//                                      March 2002

#include "FArrayBox.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "ExtrapFillPatch.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "OldPatchGodunov.H"
#include "UsingNamespace.H"

/// Level Godunov - really a hyperbolic level integrator
/**
 */
class OldLevelGodunov
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  OldLevelGodunov();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~OldLevelGodunov();

  /// Actual constructor.
  /**
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const
     (coding standards).  a_refine is the refinement ratio between
     this level and the next coarser level.  For the coarsest level, an
     empty DisjointBoxLayout is passed in for coarserDisjointBoxLayout.
   */
  void define(const DisjointBoxLayout&     a_thisDisjointBoxLayout,
              const DisjointBoxLayout&     a_coarserDisjointBoxLayout,
              const ProblemDomain&         a_domain,
              const int&                   a_refineCoarse,
              const Real&                  a_dx,
              const OldPatchGodunov* const a_patchGodunovFactory,
              const bool&                  a_hasCoarser,
              const bool&                  a_hasFiner);

  /// Take one timestep for this grid level.
  /**
     For the coarsest level, UcoarseOld and UCoarseNew are empty
     LevelData<FArrayBox> objects.  Also, either UCoarseOld or
     UCoarseNew might be empty to indicate that t(nf) = t(nc) the
     one grid is at the current time and interpolation is not required
     for boundary condition generation.  U must be defined on the same
     DisjointBoxLayouts as were used in define().  a_flux is
     is used more generally to pass
     face-centered data in and out of the function. Coarse flux register
     is flux register with the next coarser level.  Fine flux register is
     the flux register with the next finer level.  To the finer level FR,
     this level is the coarse level.  To the coarser level FR, this level
     is the fine level.  If source terms do not exist, "a_S" should be null
     constructed and not defined (i.e. its define() should not be called).
   */
  Real step(LevelData<FArrayBox>&       a_U,
            LevelData<FArrayBox>        a_flux[CH_SPACEDIM],
            LevelFluxRegister&          a_finerFluxRegister,
            LevelFluxRegister&          a_coarserFluxRegister,
            const LevelData<FArrayBox>& a_S,
            const LevelData<FArrayBox>& a_UCoarseOld,
            const Real&                 a_TCoarseOld,
            const LevelData<FArrayBox>& a_UCoarseNew,
            const Real&                 a_TCoarseNew,
            const Real&                 a_time,
            const Real&                 a_dt);

  /// Get maximum wave speed
  /**
    */
  Real getMaxWaveSpeed(const LevelData<FArrayBox>& a_U);


protected:
  // define() has been called
  bool m_defined;

  // Box layout for this level
  DisjointBoxLayout m_grids;

  // Patch integrator
  OldPatchGodunov* m_patchGodunov;

  // Number of ghost cells need locally for this level
  int m_numGhost;

  // Interpolator for filling in ghost cells from the next coarser level
  PiecewiseLinearFillPatch m_patcher;

  // Grid spacing
  Real m_dx;

  // Problem domain - index space for this level
  ProblemDomain m_domain;

  // Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  // Number of conserved variables
  int m_numCons;

  // Number of conservative fluxes and other face centered quantities
  int m_numFluxes;

  // Flags for the existence of coarser and finer levels, respectively
  bool m_hasCoarser;
  bool m_hasFiner;


private:
  // Disallowed for all the usual reasons
  void operator=(const OldLevelGodunov& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  OldLevelGodunov(const OldLevelGodunov& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
