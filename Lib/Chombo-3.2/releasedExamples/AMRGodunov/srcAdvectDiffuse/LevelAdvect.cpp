#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "AMRIO.H"
#include "SPMD.H"
#include "PhysIBC.H"
#include "LoHiSide.H"
#include "CH_Timer.H"
#include "AdvectPhysicsF_F.H"
#include "LevelAdvect.H"

#include "NamespaceHeader.H"

/***/
void
LevelAdvect::averageVelToCC(FArrayBox&           a_normalVel,
                            const FluxBox&       a_advectionVel,
                            const Box&           a_box) const
{
  CH_TIME("leveladvect::aveVeltoCC");
  for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      //treat every cell as regular
      FORT_AVEFACETOCELL(CHF_FRA(a_normalVel),
                         CHF_CONST_FRA1(a_advectionVel[faceDir],0),
                         CHF_CONST_INT(faceDir),
                         CHF_BOX(a_box));
    }
}

// Define the object so that time stepping can begin
void
LevelAdvect::define(const AdvectPhysics&        a_gphys,
                    const DisjointBoxLayout&    a_thisDisjointBoxLayout,
                    const DisjointBoxLayout&    a_coarserDisjointBoxLayout,
                    const ProblemDomain&        a_domain,
                    const int&                  a_refineCoarse,
                    const bool&                 a_useLimiting,
                    const Real&                 a_dx,
                    const bool&                 a_hasCoarser,
                    const bool&                 a_hasFiner)
{
  CH_TIME("LevelAdvect::define");

  // Sanity checks
  CH_assert(a_refineCoarse > 0);
  CH_assert(a_dx > 0.0);

  // Make a copy of the current grids
  m_grids  = a_thisDisjointBoxLayout;

  // Order of the normal predictor (1 -> PLM, 2-> PPM)
  m_normalPredOrder = 1;

  m_useFourthOrderSlopes = false;
  m_usePrimLimiting      = a_useLimiting;
  m_useCharLimiting      = false;
  m_useFlattening        = false;

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = false;
  m_artificialViscosity    = -1.0;

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;

  AdvectPhysics* godPhys = (AdvectPhysics*)(a_gphys.new_godunovPhysics());
  godPhys->define(m_domain, m_dx);
  m_patchGodunov.define( m_domain,  m_dx,
                         godPhys,
                         m_normalPredOrder,
                         m_useFourthOrderSlopes,
                         m_usePrimLimiting,
                         m_useCharLimiting,
                         m_useFlattening,
                         m_useArtificialViscosity,
                         m_artificialViscosity);

  //not used again--was just used as a factory
  //this makes me nervous because if someone changes the internals
  //of PatchGodunov to just store the pointer, this will break.
  delete godPhys;

  //need to store the one that lives in patchGodunov
  //so we can set the velocities
  m_advectPhysics = (AdvectPhysics*)(m_patchGodunov.getGodunovPhysicsPtr());
  m_advectPhysics->define(m_domain, m_dx);
  PhysIBC* physIBCPtr = m_advectPhysics->getPhysIBC();
  physIBCPtr->define(m_domain, m_dx);

  m_numGhost  = 4;
  m_numCons   = 1;

  m_exchangeCopier.exchangeDefine(a_thisDisjointBoxLayout,
                                  m_numGhost*IntVect::Unit);

  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

  // Set up the interpolator if there is a coarser level
  if (m_hasCoarser)
    {
      m_patcher.define(a_thisDisjointBoxLayout,
                       a_coarserDisjointBoxLayout,
                       m_numCons,
                       coarsen(a_domain,a_refineCoarse),
                       a_refineCoarse,
                       m_numGhost);
    }

  // Everything is defined
  m_isDefined = true;
}

void
LevelAdvect::fillGhost(LevelData<FArrayBox>&       a_U,
                       const LevelData<FArrayBox>& a_UCoarseOld,
                       const Real&                 a_TCoarseOld,
                       const LevelData<FArrayBox>& a_UCoarseNew,
                       const Real&                 a_TCoarseNew,
                       const Real&                 a_dt,
                       const Real&                 a_time)
{
  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

  a_U.exchange(m_exchangeCopier);

  // Fill a_U's ghost cells using fillInterp
  if (m_hasCoarser)
    {
      // Fraction "a_time" falls between the old and the new coarse times
      Real alpha = (a_time - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);

      // Truncate the fraction to the range [0,1] to remove floating-point
      // subtraction roundoff effects
      Real eps = 0.04 * a_dt / m_refineCoarse;

      if (Abs(alpha) < eps)
        {
          alpha = 0.0;
        }

      if (Abs(1.0-alpha) < eps)
        {
          alpha = 1.0;
        }

      // Current time before old coarse time
      if (alpha < 0.0)
        {
          MayDay::Error( "LevelAdvect::step: alpha < 0.0");
        }

      // Current time after new coarse time
      if (alpha > 1.0)
        {
          MayDay::Error( "LevelAdvect::step: alpha > 1.0");
        }

      // Interpolate ghost cells from next coarser level using both space
      // and time interpolation
      m_patcher.fillInterp(a_U,
                           a_UCoarseOld,
                           a_UCoarseNew,
                           alpha,
                           0,0,m_numCons);
    }
}

Real
LevelAdvect::step(LevelData<FArrayBox>&       a_U,
                  LevelFluxRegister&          a_finerFluxRegister,
                  LevelFluxRegister&          a_coarserFluxRegister,
                  LevelData<FluxBox>&         a_advectionVelocity,
                  const LevelData<FArrayBox>& a_S,
                  const LevelData<FArrayBox>& a_UCoarseOld,
                  const Real&                 a_TCoarseOld,
                  const LevelData<FArrayBox>& a_UCoarseNew,
                  const Real&                 a_TCoarseNew,
                  const Real&                 a_time,
                  const Real&                 a_dt)
{
  CH_TIME("LevelAdvect::step");

  CH_assert(m_isDefined);
  Interval UInterval(0, m_numCons-1); //(0,0)
  // Clear flux registers with next finer level
  if (m_hasFiner)
    {
      a_finerFluxRegister.setToZero();
    }

  //fill ghost cells of a_U
  fillGhost(a_U, a_UCoarseOld,  a_TCoarseOld,  a_UCoarseNew,  a_TCoarseNew, a_dt, a_time);



  // Potentially used in boundary conditions
  m_patchGodunov.setCurrentTime(a_time);

  // Dummy source used if source term passed in is empty
  FArrayBox zeroSource;

  // Use to restrict maximum wave speed away from zero
  Real maxWaveSpeed = 1.0e-12;

  // Beginning of loop through patches/grids.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      // The current box
      Box curBox = m_grids.get(dit());

      // The current grid of conserved variables
      FArrayBox& curU = a_U[dit];

      // The current source terms if they exist
      const FArrayBox* source = &zeroSource;
      if (a_S.isDefined())
        {
          source = &a_S[dit];
        }

      FluxBox& advVel = a_advectionVelocity[dit];
      const Box& cellBox = advVel.box();
      FArrayBox cellVel(cellBox, SpaceDim);
      averageVelToCC(cellVel, advVel, cellBox);
      m_advectPhysics->setVelocities(&cellVel, &advVel);
      // Set the current box for the patch integrator
      m_patchGodunov.setCurrentBox(curBox);

      Real maxWaveSpeedGrid;

      // The fluxes computed for this grid - used for refluxing and returning
      // other face centered quantities
      FluxBox flux;

      // Update the current grid's conserved variables, return the final
      // fluxes used for this, and the maximum wave speed for this grid
      m_patchGodunov.updateState(curU,
                                 flux,
                                 maxWaveSpeedGrid,
                                 *source,
                                 a_dt,
                                 curBox);

      // Clamp away from zero
      maxWaveSpeed = Max(maxWaveSpeed,maxWaveSpeedGrid);

      // Do flux register updates
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          // Increment coarse flux register between this level and the next
          // finer level - this level is the next coarser level with respect
          // to the next finer level
          if (m_hasFiner)
            {
              a_finerFluxRegister.incrementCoarse(flux[idir],a_dt,dit(),
                                                  UInterval,
                                                  UInterval,idir);
            }

          // Increment fine flux registers between this level and the next
          // coarser level - this level is the next finer level with respect
          // to the next coarser level
          if (m_hasCoarser)
            {
              a_coarserFluxRegister.incrementFine(flux[idir],a_dt,dit(),
                                                  UInterval,
                                                  UInterval,idir);
            }
        }

    }

  // Find the minimum of dt's over this level
  Real local_dtNew = m_dx / maxWaveSpeed;
  Real dtNew;

#ifdef CH_MPI
  int result = MPI_Allreduce(&local_dtNew, &dtNew, 1, MPI_CH_REAL,
                             MPI_MIN, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    {
      MayDay::Error("LevelAdvect::step: MPI communcation error");
    }
#else
  dtNew = local_dtNew;
#endif

  // Return the maximum stable time step
  return dtNew;
}

// Find the maximum wave speed on the current level
Real
LevelAdvect::getMaxWaveSpeed(const LevelData<FArrayBox>& a_U,
                             LevelData<FluxBox>&         a_advectionVelocity)
{
  CH_TIME("LevelAdvect::getMaxWaveSpeed");

  const DisjointBoxLayout& disjointBoxLayout = a_U.disjointBoxLayout();
  DataIterator dit = disjointBoxLayout.dataIterator();

  // Initial maximum wave speed
  Real speed = 0.0;

  // This computation doesn't need to involve a time, but the time being set
  // is checked by PatchGodunov::getMaxWaveSpeed(), so we have to set it
  // to something...
  m_patchGodunov.setCurrentTime(0.0);

  // Loop over all grids to get the maximum wave speed
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& curBox = disjointBoxLayout.get(dit());

      FluxBox& advVel = a_advectionVelocity[dit];
      const Box& cellBox = advVel.box();
      FArrayBox cellVel(cellBox, SpaceDim);
      averageVelToCC(cellVel, advVel, cellBox);
      m_advectPhysics->setVelocities(&cellVel, &advVel);

      // Set the current box and get the maximum wave speed on the current grid
      m_patchGodunov.setCurrentBox(curBox);

      // Get maximum wave speed for this grid
      Real speedOverBox = m_advectPhysics->getMaxWaveSpeed(a_U[dit], curBox);

      // Compute a running maximum
      speed = Max(speed,speedOverBox);
    }

  // Gather maximum wave speeds and broadcast the maximum over these
  Vector<Real> allSpeeds;

  gather(allSpeeds,speed,uniqueProc(SerialTask::compute));

  if (procID() == uniqueProc(SerialTask::compute))
    {
      speed = allSpeeds[0];
      for (int i = 1; i < allSpeeds.size (); ++i)
        {
          speed = Max(speed,allSpeeds[i]);
        }
    }

  broadcast(speed,uniqueProc(SerialTask::compute));

  // Return the maximum wave speed
  return speed;
}

#include "NamespaceFooter.H"
