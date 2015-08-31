#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef LEVEL_GODUNOV_H
#define LEVEL_GODUNOV_H

#include "FluxBox.H"
#include "DisjointBoxLayout.H"
#include "ExtrapFillPatch.H"
#include "LevelFluxRegister.H"

#include "PatchPluto.H"
#include "NamespaceHeader.H"

/// Level Pluto - really a hyperbolic level integrator
/**
 */
class LevelPluto
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelPluto();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelPluto();

  /// Actual constructor.
  /**
     Inside the routine, we cast away const-ness on the data members
     for the assignment.  The arguments passed in are maintained const
     (coding standards).  a_refine is the refinement ratio between
     this level and the next coarser level.  For the coarsest level, an
     empty DisjointBoxLayout is passed in for coarserDisjointBoxLayout.
   */
  void define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
              const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
              const ProblemDomain&      a_domain,
              const int&                a_refineCoarse,
              const int&                a_level,
              const Real&               a_dx,
              const PatchPluto*  const  a_patchPlutoFactory,
              const bool&               a_hasCoarser,
              const bool&               a_hasFiner);

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
            LevelData<FArrayBox>&       a_split_tags,
            const LevelData<FArrayBox>& a_UCoarseOld,
            const Real&                 a_TCoarseOld,
            const LevelData<FArrayBox>& a_UCoarseNew,
            const Real&                 a_TCoarseNew,
            const Real&                 a_time,
            const Real&                 a_dt,
            const Real&                 a_cfl);

  void setGridLevel();

  /* Used by AMRLevelPluto to get the m_dV LevelData<FArrayBox> */
  #if GEOMETRY != CARTESIAN
   const LevelData<FArrayBox>& getdV() const;
  #endif

  Real getDlMin();
 
protected:
  // Box layout for this level
  DisjointBoxLayout m_grids;

  // Union of struct GRID for this level
  LayoutData<PatchGrid> m_structs_grid;

  // Storage for informations about cell volumes (1 = dV/m_dx^3, 2 = cylindrical radius)
  #if GEOMETRY != CARTESIAN
   LevelData<FArrayBox> m_dV;
  #endif

  // Patch integrator
  PatchPluto* m_patchPluto;

  // Number of ghost cells need locally for this level
  int m_numGhost;

  // Exchange copier
  Copier m_exchangeCopier;

 #if (TIME_STEPPING == RK2)
  // Flags to be passed from predictor to corrector (only RK2)
  LevelData< BaseFab<unsigned char> > m_Flags;
  // Temporary storage for conserved variables
  LevelData<FArrayBox> m_Utmp;
 #endif

  // Temporary storage space for conserved variables
  LevelData<FArrayBox> m_U;

  // Interpolator for filling in ghost cells from the next coarser level
   PiecewiseLinearFillPluto m_patcher;

  // Grid spacing
  Real m_dx;

  // Minimum cell size
  Real m_dl_min;

  // Level number
  int m_level;

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

  // Has this object been defined
  bool m_isDefined;

private:
  // Disallowed for all the usual reasons
  void operator=(const LevelPluto& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  LevelPluto(const LevelPluto& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
