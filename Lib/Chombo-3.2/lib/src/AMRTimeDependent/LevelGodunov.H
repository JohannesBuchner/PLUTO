#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELGODUNOV_H_
#define _LEVELGODUNOV_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "PiecewiseLinearFillPatch.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "Copier.H"

#include "PatchGodunov.H"
#include "NamespaceHeader.H"

/// Level Godunov - really a hyperbolic level integrator
/**
 */
class LevelGodunov
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelGodunov();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelGodunov();

  /// Actual constructor.
  /**
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const
     (coding standards).
   */
  void define(/// box layout at this level
              const DisjointBoxLayout&    a_thisDisjointBoxLayout,
              /// box layout at next coarser level (or empty if this is coarsest level)
              const DisjointBoxLayout&    a_coarserDisjointBoxLayout,
              /// problem domain at this level
              const ProblemDomain&        a_domain,
              /// refinement ratio between this level and next coarser level
              const int&                  a_refineCoarse,
              /// grid spacing at this level
              const Real&                 a_dx,
              /// pointer to GodunovPhysics class
              const GodunovPhysics* const a_godunovFactory,
              /// order of the normal predictor: 1 for PLM, 2 for PPM
              const int&                  a_normalPredOrder,
              /// whether to use 4th-order slopes
              const bool&                 a_useFourthOrderSlopes,
              /// whether to apply slope limiting to primitive variables
              const bool&                 a_usePrimLimiting,
              /// whether to apply slope limiting to characteristic variables
              const bool&                 a_useCharLimiting,
              /// whether to apply slope flattening
              const bool&                 a_useFlattening,
              /// whether to apply artificial viscosity
              const bool&                 a_useArtificialViscosity,
              /// artificial viscosity coefficient
              const Real&                 a_artificialViscosity,
              /// whether there is a coarser level
              const bool&                 a_hasCoarser,
              /// whether there is a finer level
              const bool&                 a_hasFiner);

  /// Take one timestep for this grid level.
  Real step(/// conserved variables at this level, defined on a_thisDisjointBoxLayout in define()
            LevelData<FArrayBox>&       a_U,
            /// flux, for passing face-centered data in and out of the function
            LevelData<FArrayBox>        a_flux[CH_SPACEDIM],
            /// flux register with next finer level
            LevelFluxRegister&          a_finerFluxRegister,
            /// flux register with next coarser level
            LevelFluxRegister&          a_coarserFluxRegister,
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

  /// Compute the time-centered values of the primitive variables on cell faces.
  /**
      This API is used in cases where some operation over the whole
      level must be performed on the face-centered variables prior to
      the final difference update. Examples include incompressible
      flow and MHD, in which it is necessary to compute the projection
      of a face-centered vector field on its divergence-free part. To
      complete the differencing, it is necessary to call the member
      function computeUpdate().
  */
  void computeWHalf(/// primitive variables extrapolated to cell faces and a half time step
                    LayoutData<FluxBox>&        a_WHalf,
                    /// conserved variables on current grids
                    LevelData<FArrayBox>&       a_U,
                    /// source term on current grid
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

  /// Compute the increment in the conserved variables from face variables.
  /**
     Compute dU = dt*dU/dt, the change in the conserved variables over
     the time step.  Function returns the maximum stable time step.
   */
  Real computeUpdate(/// change in conserved variables over the time step
                     LevelData<FArrayBox>&       a_dU,
                     /// flux register with next finer level
                     LevelFluxRegister&          a_finerFluxRegister,
                     /// flux register with next coarser level
                     LevelFluxRegister&          a_coarserFluxRegister,
                     /// conserved variables at beginning of time step
                     const LevelData<FArrayBox>& a_U,
                     /// primitive variables extrapolated to cell faces and a half time step
                     const LayoutData<FluxBox>&  a_WHalf,
                     /// current time
                     const Real&                 a_time,
                     /// time step
                     const Real&                 a_dt);

  /// Return maximum wave speed
  /**
    */
  Real getMaxWaveSpeed(const LevelData<FArrayBox>& a_U);

  /// Get the GodunovPhysics pointer from the PatchGodunov member
  /**
      Returns the pointer to the GodunovPhysics object that is
      used in the PatchGodunov class.
  */
  GodunovPhysics* getGodunovPhysicsPtr();

  const GodunovPhysics* getGodunovPhysicsPtrConst() const ;

  /// Set whether to use high-order limiter.
  void highOrderLimiter(bool a_highOrderLimiter);

protected:
  // Box layout for this level
  DisjointBoxLayout m_grids;

  // Patch integrator
  LayoutData<PatchGodunov>    m_patchGodunov;

  // Number of ghost cells need locally for this level
  int m_numGhost;

  // Exchange copier
  Copier m_exchangeCopier;

  // Temporary storage space for conserved variables
  LevelData<FArrayBox> m_U;

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

  // Order of the normal predictor (1 -> PLM, 2-> PPM)
  int m_normalPredOrder;

  // Use 4th order slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  // Do slope limiting in the primitive or characteristic variables, respect.
  bool m_usePrimLimiting;
  bool m_useCharLimiting;

  // Do slope flattening - MUST BE USING 4th order slopes
  bool m_useFlattening;

  // Apply artificial viscosity of a set value
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;

  // Has this object been defined
  bool m_isDefined;


private:

  // Disallowed for all the usual reasons
  void operator=(const LevelGodunov&);
  LevelGodunov(const LevelGodunov&);
};

#include "NamespaceFooter.H"
#endif
