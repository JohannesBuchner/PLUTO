#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*****************/
/*****************/

#include "AMRNavierStokes.H"
#include "SetValLevel.H"

// ---------------------------------------------------------------
Real
AMRNavierStokes::computeInitialDt()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::computeInitialDt " << m_level << endl;
    }

  return s_init_shrink*computeDt();
}

// ---------------------------------------------------------------
Real
AMRNavierStokes::computeDt()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::computeDt " << m_level << endl;
    }

  Real dt;

  // if prescribedDt > 0 then use that
  // (don't worry about not being on level 0 -- Amr class
  // will take care of that
  if (s_prescribedDt > 0.0)
    {
      dt = s_prescribedDt;
    }
  else
    {
      Real MaxVel = 0;

      Interval velComps(0,SpaceDim-1);
      // compute maximum velocity  (MaxNorm)
      MaxVel = norm(newVel(), velComps, 0);

      if (MaxVel != 0)
        {
          dt = m_cfl*m_dx/MaxVel;
        }
      else
        {
          dt = s_max_dt;
        }

      // check to see if boundary conditions restrict dt
      m_physBCPtr->computeBoundaryDt(dt, m_cfl, m_dx);
    }

  return dt;
}

// ---------------------------------------------------------------
void
AMRNavierStokes::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

// ---------------------------------------------------------------
void
AMRNavierStokes::refinementThreshold(Real a_thresh)
{
  m_refine_threshold = a_thresh;
}

// ---------------------------------------------------------------
void
AMRNavierStokes::limitSolverCoarsening(bool a_limitSolverCoarsening)
{
  m_limitSolverCoarsening = a_limitSolverCoarsening;
  // also pass through to projection
  m_projection.limitSolverCoarsening(a_limitSolverCoarsening);
}

// ---------------------------------------------------------------
bool
AMRNavierStokes::finestLevel() const
{
  return m_finest_level;
}

// ---------------------------------------------------------------
bool
AMRNavierStokes::isEmpty() const
{
  return m_is_empty;
}

// ---------------------------------------------------------------
void
AMRNavierStokes::finestLevel(bool a_finest_level)
{
  m_finest_level = a_finest_level;
  m_projection.isFinestLevel(a_finest_level);
}

// ---------------------------------------------------------------
Real
AMRNavierStokes::Dx() const
{
  return m_dx;
}

// ---------------------------------------------------------------
AMRNavierStokes*
AMRNavierStokes::finerNSPtr() const
{
  AMRNavierStokes* fine_ns_ptr;
  if (m_finer_level_ptr != NULL)
    {
      fine_ns_ptr = dynamic_cast <AMRNavierStokes*> (m_finer_level_ptr);
      if (fine_ns_ptr == NULL)
        {
          MayDay::Error ("AMRNavierStokes::fineNSPtr -- fineptr not castable to AMRNavierStokes*");
        }
    }
  else
    {
      fine_ns_ptr = NULL;
    }

  return fine_ns_ptr;
}

// ---------------------------------------------------------------
AMRNavierStokes*
AMRNavierStokes::crseNSPtr() const
{
  AMRNavierStokes* crse_ns_ptr;
  if (m_coarser_level_ptr != NULL)
    {
      crse_ns_ptr = dynamic_cast <AMRNavierStokes*> (m_coarser_level_ptr);
      if (crse_ns_ptr == NULL)
        {
          MayDay::Error ("AMRNavierStokes::crseNSPtr -- crseptr not castable to AMRNavierStokes*");
        }
    }
  else
    {
      crse_ns_ptr = NULL;
    }

  return crse_ns_ptr;
}

// ---------------------------------------------------------------
LevelData<FArrayBox>&
AMRNavierStokes::newVel()
{
  CH_assert(m_vel_new_ptr != NULL);

  return *(m_vel_new_ptr);
}

// ---------------------------------------------------------------
LevelData<FArrayBox>&
AMRNavierStokes::oldVel()
{
  CH_assert(m_vel_old_ptr != NULL);

  return *(m_vel_old_ptr);
}

// ---------------------------------------------------------------
LevelData<FArrayBox>&
AMRNavierStokes::newLambda()
{
  CH_assert(m_lambda_new_ptr != NULL);

  return *(m_lambda_new_ptr);
}

// ---------------------------------------------------------------
LevelData<FArrayBox>&
AMRNavierStokes::oldLambda()
{
  CH_assert(m_lambda_old_ptr != NULL);

  return *(m_lambda_old_ptr);
}

// ---------------------------------------------------------------
LevelData<FArrayBox>&
AMRNavierStokes::newScal(const int a_comp)
{
  CH_assert(m_scal_new[a_comp] != NULL);

  return *(m_scal_new[a_comp]);
}

// ---------------------------------------------------------------
LevelData<FArrayBox>&
AMRNavierStokes::oldScal(const int a_comp)
{
  CH_assert(m_scal_old[a_comp] != NULL);

  return *(m_scal_old[a_comp]);
}

// ---------------------------------------------------------------
const LevelData<FArrayBox>&
AMRNavierStokes::newVel() const
{
  CH_assert(m_vel_new_ptr != NULL);

  return *(m_vel_new_ptr);
}

// ---------------------------------------------------------------
const LevelData<FArrayBox>&
AMRNavierStokes::oldVel() const
{
  CH_assert(m_vel_old_ptr != NULL);

  return *(m_vel_old_ptr);
}

// ---------------------------------------------------------------
const LevelData<FArrayBox>&
AMRNavierStokes::newLambda() const
{
  CH_assert(m_lambda_new_ptr != NULL);

  return *(m_lambda_new_ptr);
}

// ---------------------------------------------------------------
const LevelData<FArrayBox>&
AMRNavierStokes::oldLambda() const
{
  CH_assert(m_lambda_old_ptr != NULL);

  return *(m_lambda_old_ptr);
}

// ---------------------------------------------------------------
const LevelData<FArrayBox>&
AMRNavierStokes::newScal(const int a_comp) const
{
  CH_assert(m_scal_new[a_comp] != NULL);

  return *(m_scal_new[a_comp]);
}

// ---------------------------------------------------------------
const LevelData<FArrayBox>&
AMRNavierStokes::oldScal(const int a_comp) const
{
  CH_assert(m_scal_old[a_comp] != NULL);

  return *(m_scal_old[a_comp]);
}

// ---------------------------------------------------------------
void
AMRNavierStokes::setPhysBC(const PhysBCUtil& a_BC)
{
  m_physBCPtr = a_BC.newPhysBCUtil();
}

// ---------------------------------------------------------------
PhysBCUtil*
AMRNavierStokes::getPhysBCPtr() const
{
  return m_physBCPtr;
}

// ---------------------------------------------------------------
// this also changes time levels
void
AMRNavierStokes::swapOldAndNewStates()
{
  LevelData<FArrayBox>* tempPtr;

  // swap both velocity and scalar data
  tempPtr = m_vel_new_ptr;
  m_vel_new_ptr = m_vel_old_ptr;
  m_vel_old_ptr = tempPtr;

  tempPtr = m_lambda_new_ptr;
  m_lambda_new_ptr = m_lambda_old_ptr;
  m_lambda_old_ptr = tempPtr;

  // loop over scalar components
  CH_assert (m_scal_new.size() == m_scal_old.size());

  for (int comp=0; comp<m_scal_new.size(); comp++)
    {
      tempPtr = m_scal_new[comp];
      m_scal_new[comp] = m_scal_old[comp];
      m_scal_old[comp] = tempPtr;
    }

  time(m_time + m_dt);
}

// ---------------------------------------------------------------
// kinda like swapOldAndNewStates, but resets time to a_time
void
AMRNavierStokes::resetStates(const Real a_time)
{
  LevelData<FArrayBox>* tempPtr;
  tempPtr = m_vel_new_ptr;
  m_vel_new_ptr = m_vel_old_ptr;
  m_vel_old_ptr = tempPtr;

  tempPtr = m_lambda_new_ptr;
  m_lambda_new_ptr = m_lambda_old_ptr;
  m_lambda_old_ptr = tempPtr;

  // loop over scalar components
  CH_assert (m_scal_new.size() == m_scal_old.size());

  for (int comp=0; comp<m_scal_new.size(); comp++)
    {
      tempPtr = m_scal_new[comp];
      m_scal_new[comp] = m_scal_old[comp];
      m_scal_old[comp] = tempPtr;
    }

  time(a_time);
}


// ---------------------------------------------------------------
void
AMRNavierStokes::setAllBogus()
{
  if (s_set_bogus_values)
    {
      setValLevel(*m_vel_old_ptr, s_bogus_value);
      setValLevel(*m_vel_new_ptr, s_bogus_value);
      setValLevel(*m_lambda_old_ptr, s_bogus_value);
      setValLevel(*m_lambda_new_ptr, s_bogus_value);

      setValLevels(m_scal_old, 0, s_num_scal_comps-1, s_bogus_value);
      setValLevels(m_scal_new, 0, s_num_scal_comps-1, s_bogus_value);
    }
}
