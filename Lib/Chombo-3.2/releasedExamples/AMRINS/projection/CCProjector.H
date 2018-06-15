#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CCPROJECTOR_H_
#define _CCPROJECTOR_H_

#include <cstdlib>
#include <cmath>
#include <cassert>

#include <iostream>
using std::istream;
using std::cout;
using std::cerr;

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRMultiGrid.H"
#include "QuadCFInterp.H"
#include "PhysBCUtil.H"

#include "FluxBox.H"

#include "UsingNamespace.H"

/// this class manages the various forms of the CC projection
class CCProjector
/** this class performs the various cell-centered projections required
by the IAMR algorithm. These include the single-level operators levelMac
and levelProject, as well as the multilevel sync projection and volume
discrepancy solves. */
{
public:
  ///
  CCProjector();

  /// full constructor
  CCProjector(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_crseGridsPtr,
              const Box& a_domain,
              const Real  a_dx,
              CCProjector* a_finerProj,
              CCProjector* a_crseProj,
              int a_nRefCrse,
              int a_level,
              const PhysBCUtil& a_physBC);

  /// full constructor
  CCProjector(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_crseGridsPtr,
              const ProblemDomain& a_domain,
              const Real  a_dx,
              CCProjector* a_finerProj,
              CCProjector* a_crseProj,
              int a_nRefCrse,
              int a_level,
              const PhysBCUtil& a_physBC);

  /// destructor
  ~CCProjector();

  /// define function
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_crseGridsPtr,
              const Box& a_domain,
              const Real  a_dx,
              CCProjector* a_finerProj,
              CCProjector* a_crseProj,
              int a_nRefCrse,
              int a_level,
              const PhysBCUtil& a_physBC);

  /// define function
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_crseGridsPtr,
              const ProblemDomain& a_domain,
              const Real  a_dx,
              CCProjector* a_finerProj,
              CCProjector* a_crseProj,
              int a_nRefCrse,
              int a_level,
              const PhysBCUtil& a_physBC);

  /// initialize new projection with data from old projection
  void init(const CCProjector& a_oldProj);

  /// define static parts
  void variableSetUp();

  ///
  void setCrseProj(CCProjector* a_crseProj, int nRefCrse);

  ///
  void setFineProj(CCProjector* a_fineProj);

  /// set solver parameter
  void limitSolverCoarsening(bool a_limitSolverCoarsening);

  /// write checkpoint header
  void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// write this class to a checkpoint file for later restart
  void  writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// read the checkpoint header
  void readCheckpointHeader(HDF5Handle& a_handle);

  /// read this class from a checkpoint file
  void readCheckpointLevel(HDF5Handle& a_handle);

  /// do level MAC projection, correct uEdge.
  /*** a_oldTime is
       old-time on AMRLevel which is calling the projection.
       assumes BC's already set on uEdge values */
  void levelMacProject(LevelData<FluxBox>& uEdge,
                       Real a_oldTime, Real a_dt);

  /// do level projection and correct (cell-centered) velocities;
  /*** assumes physical and exchange (copy) boundary
       conditions already set.
  */
  void LevelProject(LevelData<FArrayBox>& a_velocity,
                    LevelData<FArrayBox>* a_crseVelPtr,
                    const Real a_newTime, const Real a_dt);

  /** defines AMRSolver and calls sync projection
      and freestream preservation solve. assumes all physical
      and copy BC's already set (does C/F BC's however) */
  void doSyncOperations(Vector<LevelData<FArrayBox>* >& a_velocity,
                        Vector<LevelData<FArrayBox>* >& a_lambda,
                        const Real a_newTime, const Real a_dtSync);

  /// Initialization functions
  /** performs initial level projection -- phys and copy BC's
      should be preset */
  void initialLevelProject(LevelData<FArrayBox>& a_velocity,
                           LevelData<FArrayBox>* a_crseVelPtr,
                           const Real a_oldTime, const Real a_newTime);

  /** defines AMRSolver and calls initial sync projection and
      freestream preservation solves */
  void doInitialSyncOperations(Vector<LevelData<FArrayBox>* >& a_vel,
                               Vector<LevelData<FArrayBox>* >& a_lambda,
                               const Real a_newTime, const Real a_dtSync);

  /// performs multilevel projection on velocity, modifies velocity in place
  /**  dtSync is used as a scaling factor for coarse-level BC's,
       not used in actual projection.  if homogeneousCFBC is false, C/F bc
       for l_base is taken from coarse-level eSync, o/w use l_base-1
       correction field = 0.  if applyCorrection == false, then correction
       is computed, but not applied (can be useful as a diagnostic).
       physical and copy BC's should already be set before this is called.
  */

  void initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_velocity,
                              bool a_homogeneousCFBC = true);

  /// same as other initialVelocityProject, but uses a pre-defined AMRSolver
  void initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_velocity,
                              AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
                              bool a_homogeneousCFBC = true);

  ////  performs the necessary re-initialization after regridding
  /*** does composite projection on new velocity field
       and also recomputes e_lambda and grad_eLambda */
  void doPostRegridOps(Vector<LevelData<FArrayBox>* >& a_velocity,
                       Vector<LevelData<FArrayBox>* >& a_lambda,
                       const Real a_dt, const Real a_time);

  /// access functions

  ///
  int getLevel() const;

  ///
  const ProblemDomain& dProblem() const;

  ///
  const DisjointBoxLayout& getBoxes() const;

  ///
  int nRefCrse() const;

  ///
  Real dx() const;

  /// returns coefficient for volume-discrepancy solve.
  Real etaLambda() const;

  ///
  CCProjector* fineProjPtr() const;

  ///
  CCProjector* crseProjPtr() const;

  /// returns MAC correction
  LevelData<FArrayBox>& phi();

  /// const version of accessor
  const LevelData<FArrayBox>& phi() const;

  /// returns edge-centered grad(phi) in direction dir
  void gradPhi(LevelData<FArrayBox>& a_gradPhi, int a_dir) const;

  /// returns all components of grad(phi) (gradPhi should be correct size)
  void gradPhi(LevelData<FluxBox>& a_gradPhi) const;

  /// returns level-projection pressure
  LevelData<FArrayBox>& Pi();

  /// returns grad(pi) in direction dir
  void gradPi(LevelData<FArrayBox>& a_gradPi, int a_dir) const;

  /// returns grad(pi) in all directions into SpaceDim-dimensioned gradPi
  void gradPi(LevelData<FArrayBox>& a_gradPi) const;

  /// returns synchronization correction
  LevelData<FArrayBox>& eSync();

  /// returns synchronization correction (const version)
  const LevelData<FArrayBox>& eSync() const;

  /// returns cell-centered grad(eSync) (composite gradient)
  void grad_eSync(LevelData<FArrayBox>& a_grad_eSync, int a_dir) const;

  /// returns cell-centered G^{comp}(eSync)
  void grad_eSync(LevelData<FArrayBox>& a_grad_eSync) const;

  /// returns volume-discrepancy correction
  LevelData<FArrayBox>& eLambda();

  /// returns volume-discrepancy correction (const version)
  const LevelData<FArrayBox>& eLambda() const;

  /// returns edge-centered grad(eLambda) in direction dir
  void grad_eLambda(LevelData<FArrayBox>& a_grad_eLambda, int a_dir) const;

  /// returns edge-centered grad(eLambda) in all directions
  /** non-constant reference returned because of need to rescale. */
  LevelData<FluxBox>&  grad_eLambda();

  ///
  const LevelData<FluxBox>& grad_eLambda() const;

  /// do sync projection?
  bool doSyncProjection() const;

  /// do volume discrepancy correction?
  bool doMacSync() const;

  /// has this object been completely initialized?
  bool isInitialized() const;

  /// use quadratic interpolation (instead of extrap) for c/f bc
  bool doQuadInterp() const;

  /// returns predefined quadratic coarse-fine interpolation object
  QuadCFInterp& quadCFInterpolator();

  /// returns predefined intvectset which shows coverage
  const LayoutData<IntVectSet>& getGridIVS();

  /// is this the finest level?
  bool isFinestLevel() const;

  /// set whether this is the finest level
  void isFinestLevel(bool a_finest_level);

  /// set verbosity
  void verbosity(int a_verbosity);

  /// returns verbosity
  int verbosity() const;

  /// sets physBCs
  void setPhysBC(const PhysBCUtil& a_bc);

  ///
  PhysBCUtil* getPhysBCPtr() const;

protected:
  /// returns uEdge - G(phi_mac); C/F BC is CFscale*Pi(crse)
  void applyMacCorrection(LevelData<FluxBox>& a_uEdge, Real CFscale);

  /// vel := vel - scale*G_CC(pi)
  void correctCCVelocities(LevelData<FArrayBox>& a_velocity,
                           const Real scale) const;

  /// do sync projection with this level as lBase; corrects velocity in place
  void doSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
                        const Real a_newTime, const Real a_dtSync,
                        AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                        );

  /// perform initial sync projection
  void initialSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
                             const Real a_newTime, const Real a_dtSync,
                             AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                             );

  /// do volume discrepancy solve, with this level as lBase
  void computeVDCorrection(Vector<LevelData<FArrayBox>* >& a_lambda,
                           const Real a_newTime,
                           const Real a_dtSync,
                           AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                           );

  ///  apply sync correction vel = vel - scale*G^{comp} eSync
  /** (recursively applies correction to all finer levels) */
  void applySyncCorrection(Vector<LevelData<FArrayBox>* >& a_velocity,
                           const Real scale,
                           LevelData<FArrayBox>* crseCorr);

  /// compute composite gradient of eLambda -- does this recursively
  void computeGrad_eLambda();

  /// rescales composite gradient of eLambda, recursively all finer levels
  /** rescales grad(eLambda) by (a_dx_lbase)/m_dx for this level and
      recursively calls this function for finer levels.  This is done
      for stability reasons.
  */
  void rescaleGrad_eLambda(Real a_dx_lbase);

  /// defines Multilevel AMRSolver used for composite projections
  /** uses grids in vel, lBase is this level.
      a_freestream solver is true if solver will be used for
      freestream preservation solve. */
  void defineMultiGrid(AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
                       const Vector<LevelData<FArrayBox>* >& a_vel,
                       bool a_freestreamSolver);

  /// defines m_solverMGlevel
  void defineSolverMGlevel(const DisjointBoxLayout& a_grids,
                           const DisjointBoxLayout* a_crseGridsPtr);

  /// solves using m_solverMGlevel
  void solveMGlevel(LevelData<FArrayBox>&         a_phi,
                    const LevelData<FArrayBox>*   a_phiCoarsePtr,
                    const LevelData<FArrayBox>&   a_rhs);

  /// takes care of setting ghost cell values after restarting
  void postRestart();

private:

  ///
  int m_level;

  ///
  Real m_dx;

  ///
  ProblemDomain m_domain;

  /// Quadratic C/F interpolation object
  QuadCFInterp m_cfInterp;

  /// cache intvectsets for gradient here (if relevant)
  LayoutData<IntVectSet> m_gradIVS;

  /// ptr to finer-level projector object (if it exists -- o/w == NULL)
  CCProjector* m_fineProjPtr;
  /// ptr to coarser-level projector (if level == 0, ptr == NULL)
  CCProjector* m_crseProjPtr;

  /// if coarse level exists, refinement ratio down to it
  int m_nRefCrse;

  /// cell centered pressure data at time n+1/2
  LevelData<FArrayBox> m_Pi;

  /// MAC correction field
  LevelData<FArrayBox> m_phi;

  /// synchronization correction
  LevelData<FArrayBox> m_eSync;
  /// volume-discrepancy correction
  LevelData<FArrayBox> m_eLambda;

  /// edge-centered gradient of eLambda
  LevelData<FluxBox> m_grad_eLambda;

  ///
  AMRMultiGrid<LevelData<FArrayBox> > m_solverMGlevel;

  ///
  static bool m_doSyncProjection;

  ///
  static bool m_applySyncCorrection;

  ///
  bool m_isInitialized;

  ///
  static bool m_doQuadInterp;

  /// solver parameter -- passed on through when defining solvers
  bool m_limitSolverCoarsening;

  ///
  static Real m_etaLambda;

  /// multigrid solver tolerance
  static Real s_solver_tolerance;

  /// multigrid solver parameter
  static int s_num_smooth_up;

  /// multigrid solver parameter
  static int s_num_smooth_down;

  /// number of smoothings in preconditioner for projection
  static int s_num_precond_smooth;

  /// use consistent scaling for lamdba solve?
  static bool s_constantLambdaScaling;

  /// timestep used for all subcycled lamdba solves
  static Real s_lambda_timestep;

  /// have the static variables been initialized (ParmParse)
  static bool pp_init;

  /// verbosity of output
  static int s_verbosity;

  /// if false, VD correction is set to zero once it's computed
  static bool s_applyVDCorrection;

  /// is this the finest extant level?
  bool m_finest_level;

  ///
  PhysBCUtil* m_physBCPtr;
};

#endif
