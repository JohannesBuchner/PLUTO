#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELADVECT_H_
#define _LEVELADVECT_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "Copier.H"
#include "PatchGodunov.H"
#include "AdvectPhysics.H"

#include "NamespaceHeader.H"

/// Advection integrator on a level
/**
 */
class LevelAdvect
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
  */
  LevelAdvect()
  {
    m_dx           = 0.0;
    m_refineCoarse = 0;
    m_isDefined    = false;
  }

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
  */
  ~LevelAdvect()
  {
  }

  /// Actual constructor.
  /**
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const
     (coding standards).
   */
  void define(/// advection physics class
              const AdvectPhysics&        a_gphys,
              /// box layout at this level
              const DisjointBoxLayout&    a_thisDisjointBoxLayout,
              /// box layout at next coarser level (or empty if this is coarsest level)
              const DisjointBoxLayout&    a_coarserDisjointBoxLayout,
              /// problem domain at this level
              const ProblemDomain&        a_domain,
              /// refinement ratio between this level and next coarser level
              const int&                  a_refineCoarse,
              /// whether to use limiting
              const bool&                 a_useLimiting,
              /// grid spacing at this level
              const Real&                 a_dx,
              /// whether there is a coarser level
              const bool&                 a_hasCoarser,
              /// whether there is a finer level
              const bool&                 a_hasFiner);

  /// Advance the solution by one timestep on this grid level.
  /**
   */
  Real step(/// conserved variables at this level, defined on a_thisDisjointBoxLayout in define(); gets updated in this routine
            LevelData<FArrayBox>&       a_U,
            /// flux register with next finer level
            LevelFluxRegister&          a_finerFluxRegister,
            /// flux register with next coarser level
            LevelFluxRegister&          a_coarserFluxRegister,
            /// advection velocity on faces
            LevelData<FluxBox>&         a_advectionVelocity,
            /// source term, or if no source term then this is null constructed and not defined
            const LevelData<FArrayBox>& a_S,
            /// conserved variables at coarser level at time of last coarser-level update, or empty if no coarser level; may be empty if interpolation not required
            const LevelData<FArrayBox>& a_UCoarseOld,
            /// time of last update at coarser level
            const Real&                 a_TCoarseOld,
            /// conserved variables at coarser level at time of next coarser-level update, or empty if no coarser level; may be empty if interpolation not required
            const LevelData<FArrayBox>& a_UCoarseNew,
            /// time of next update at coarser level
            const Real&                 a_TCoarseNew,
            /// current time
            const Real&                 a_time,
            /// time step
            const Real&                 a_dt);

  /// Convert velocity from face-centered to cell-centered
  /**
     In each direction, take average of normal component of velocity
     on the neighboring faces in that direction.
   */
  void averageVelToCC(/// cell-centered velocity
                      FArrayBox&                        a_normalVel,
                      /// face-centered velocity
                      const FluxBox  &                  a_advectionVel,
                      /// Box on which to return a_normalVel
                      const Box&                        a_box) const;

  /// Fill in ghost cells by exchange at this level and then by interpolation from coarser level (if any).
  void fillGhost(/// conserved variables at this level, with ghosts cells to be filled in
                 LevelData<FArrayBox>&       a_U,
                 /// conserved variables at coarser level at time of last coarser-level update, or empty if no coarser level; may be empty if interpolation not required
                 const LevelData<FArrayBox>& a_UCoarseOld,
                 /// time of last update at coarser level
                 const Real&                 a_TCoarseOld,
                 /// conserved variables at coarser level at time of next coarser-level update, or empty if no coarser level; may be empty if interpolation not required
                 const LevelData<FArrayBox>& a_UCoarseNew,
                 /// time of next update at coarser level
                 const Real&                 a_TCoarseNew,
                 /// time step
                 const Real&                 a_dt,
                 /// current time
                 const Real&                 a_time);

  /// Get maximum wave speed
  /**
    */
  Real getMaxWaveSpeed(/// conserved variables at this level
                       const LevelData<FArrayBox>& a_U,
                       /// face-centered velocities
                       LevelData<FluxBox>&         a_advectionVelocity);

protected:

  /// layout for this level
  DisjointBoxLayout m_grids;

  /// patch integrator
  PatchGodunov                     m_patchGodunov;

  /// physics class
  AdvectPhysics*                   m_advectPhysics;

  /// number of ghost cells need locally for this level
  int m_numGhost;

  /// exchange copier
  Copier m_exchangeCopier;

  /// interpolator for filling in ghost cells from the next coarser level
  PiecewiseLinearFillPatch m_patcher;

  /// grid spacing
  Real m_dx;

  /// problem domain - index space for this level
  ProblemDomain m_domain;

  /// refinement ratio between this level and the next coarser
  int m_refineCoarse;

  /// whether a coarser level exists
  bool m_hasCoarser;

  /// whether a finer level exists
  bool m_hasFiner;

  /// number of conserved variables (= 1)
  int m_numCons;//=1

  /// order of normal predictor
  int m_normalPredOrder;

  /// whether to use 4th-order slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  /// whether to do slope limiting in the primitive variables
  bool m_usePrimLimiting;

  /// whether to do slope limiting in the characteristic variables
  bool m_useCharLimiting;

  /// whether to do slope flattening - MUST BE USING 4th-order slopes
  bool m_useFlattening;

  /// whether to apply artificial viscosity of a set value
  bool m_useArtificialViscosity;

  /// artificial viscosity coefficient
  Real m_artificialViscosity;

  /// whether this object has been defined
  bool m_isDefined;

private:

  // Disallowed for all the usual reasons
  void operator=(const LevelAdvect&);
  LevelAdvect(const LevelAdvect&);
};

#include "NamespaceFooter.H"

#endif
