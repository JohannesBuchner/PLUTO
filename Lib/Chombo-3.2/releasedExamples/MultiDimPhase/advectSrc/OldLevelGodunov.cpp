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

#include "OldLevelGodunov.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LevelFluxRegister.H"
#include "PiecewiseLinearFillPatch.H"
#include "ExtrapFillPatch.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "OldPhysIBC.H"
#include "LoHiSide.H"
#include "UsingNamespace.H"

// Constructor - set up some defaults
OldLevelGodunov::OldLevelGodunov()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchGodunov = NULL;
}

// Destructor - free up storage
OldLevelGodunov::~OldLevelGodunov()
{
  if (m_patchGodunov != NULL)
    {
      delete m_patchGodunov;
    }
}

// Define the object so that time stepping can begin
void OldLevelGodunov::define(const DisjointBoxLayout&     a_thisDisjointBoxLayout,
                             const DisjointBoxLayout&     a_coarserDisjointBoxLayout,
                             const ProblemDomain&         a_domain,
                             const int&                   a_refineCoarse,
                             const Real&                  a_dx,
                             const OldPatchGodunov* const a_patchGodunovFactory,
                             const bool&                  a_hasCoarser,
                             const bool&                  a_hasFiner)
{
  // Sanity checks
  CH_assert(a_refineCoarse > 0);
  CH_assert(a_dx > 0.0);

  // Make a copy of the current grids
  m_grids  = a_thisDisjointBoxLayout;

  // Cache data
  m_dx = a_dx;
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;

  // Remove old patch integrator (if any), create a new one, and initialize
  if (m_patchGodunov != NULL)
    {
      delete m_patchGodunov;
    }

  m_patchGodunov = a_patchGodunovFactory->new_patchGodunov();
  m_patchGodunov->define(m_domain,m_dx);

  // Determing the number of ghost cells necessary here
  if (m_patchGodunov->useFourthOrderSlopes())
    {
      m_numGhost = 4;
    }
  else
    {
      CH_assert(!(m_patchGodunov->useFlattening()));
      m_numGhost = 2;
    }

  // Get the number of conserved variable and face centered fluxes
  m_numCons   = m_patchGodunov->numConserved();
  m_numFluxes = m_patchGodunov->numFluxes();

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
  m_defined = true;
}

// Advance the solution by "a_dt" by using an unsplit method.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// If source terms do not exist, "a_S" should be null constructed and not
// defined (i.e. its define() should not be called).
Real OldLevelGodunov::step(LevelData<FArrayBox>&       a_U,
                           LevelData<FArrayBox>        a_flux[CH_SPACEDIM],
                           LevelFluxRegister&          a_finerFluxRegister,
                           LevelFluxRegister&          a_coarserFluxRegister,
                           const LevelData<FArrayBox>& a_S,
                           const LevelData<FArrayBox>& a_UCoarseOld,
                           const Real&                 a_TCoarseOld,
                           const LevelData<FArrayBox>& a_UCoarseNew,
                           const Real&                 a_TCoarseNew,
                           const Real&                 a_time,
                           const Real&                 a_dt)
{
  // Make sure everything is defined
  CH_assert(m_defined);

  // Clear flux registers with next finer level
  if (m_hasFiner)
    {
      a_finerFluxRegister.setToZero();
    }

  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

  // Create temporary storage with a layer of "m_numGhost" ghost cells
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  LevelData<FArrayBox> U(m_grids,m_numCons,ivGhost);

  for (DataIterator dit = U.dataIterator(); dit.ok(); ++dit)
    {
      U[dit()].setVal(0.);
    }

  // Copy the current conserved variables into the temporary storage
  a_U.copyTo(UInterval,U,UInterval);

  // Fill U's ghost cells using fillInterp
  if (m_hasCoarser)
    {
      // Truncate the fraction to the range [0,1] to remove floating-point
      // subtraction roundoff effects
      Real eps = 0.04 * a_dt / m_refineCoarse;

      // check for current time too far outside the coarse time range
      if ( a_time+eps < a_TCoarseOld || a_time-eps > a_TCoarseNew )
      {
        pout() << "error: OldLevelGodunov::step: a_time [" << a_time
               << "] is outside the old,new range of coarse level times ["
               << a_TCoarseOld << "," << a_TCoarseNew << "]" << endl ;
        MayDay::Error( "OldLevelGodunov::step: new time is outside acceptable range" );
      }
      // if just a little outside the range (e.g. roundoff errors), just fix it
      Real curtime = a_time ;
      if ( a_time < a_TCoarseOld ) curtime = a_TCoarseOld ;
      if ( a_time > a_TCoarseNew ) curtime = a_TCoarseNew ;

      // "time" falls in the range of the old and the new coarse times
      Real alpha = (curtime - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);

      // Interpolate ghost cells from next coarser level using both space
      // and time interpolation
      m_patcher.fillInterp(U,
                           a_UCoarseOld,
                           a_UCoarseNew,
                           alpha,
                           0,0,m_numCons);
    }

  // Exchange all the data between grids at this level
  // I don't think this is necessary
  //U.exchange(UInterval);

  // Potentially used in boundary conditions
  m_patchGodunov->setCurrentTime(a_time);

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
      FArrayBox& curU = U[dit()];

      // The current source terms if they exist
      const FArrayBox* source = &zeroSource;
      if (a_S.isDefined())
        {
          source = &a_S[dit()];
        }

      // The fluxes computed for this grid - used for refluxing and returning
      // other face centered quantities
      FArrayBox flux[SpaceDim];

    // Set the current box for the patch integrator
      m_patchGodunov->setCurrentBox(curBox);

      Real maxWaveSpeedGrid;

    // Update the current grid's conserved variables, return the final
    // fluxes used for this, and the maximum wave speed for this grid
      m_patchGodunov->updateState(curU,
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
                                                  UInterval,idir,Side::Lo);
              a_coarserFluxRegister.incrementFine(flux[idir],a_dt,dit(),
                                                  UInterval,
                                                  UInterval,idir,Side::Hi);
            }
        }
    }

  // Now that we have completed the updates of all the patches, we copy the
  // contents of temporary storage, U, into the permanent storage, a_U.
  // U.copyTo(UInterval,a_U,UInterval);

  for (DataIterator dit = U.dataIterator(); dit.ok(); ++dit)
    {
      a_U[dit()].copy(U[dit()]);
    }

  // Find the minimum of dt's over this level
  Real local_dtNew = m_dx / maxWaveSpeed;
  Real dtNew;

#ifdef CH_MPI
  int result = MPI_Allreduce(&local_dtNew, &dtNew, 1, MPI_CH_REAL,
                                 MPI_MIN, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
  {
    MayDay::Error("sorry, but I had a communcation error on new dt");
  }
#else
  dtNew = local_dtNew;
#endif

  // Return the maximum stable time step
  return dtNew;
}

// Find the maximum wave speed on the current level
Real OldLevelGodunov::getMaxWaveSpeed(const LevelData<FArrayBox>& a_U)
{
  const DisjointBoxLayout& disjointBoxLayout = a_U.disjointBoxLayout();
  DataIterator dit = disjointBoxLayout.dataIterator();

  // Initial maximum wave speed
  Real speed = 0.0;

  // This computation doesn't need involve a time but the time being set
  // is checked by OldPatchGodunov::getMaxWaveSpeed so we have to set it
  // to something...
  m_patchGodunov->setCurrentTime(0.0);

  // Loop over all grids to get the maximum wave speed
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& currentBox = disjointBoxLayout.get(dit());

      // Set the current box and get the maximum wave speed on the current grid
      m_patchGodunov->setCurrentBox(currentBox);
      Real speedOverBox = m_patchGodunov->getMaxWaveSpeed(a_U[dit()],
                                                          currentBox);

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
