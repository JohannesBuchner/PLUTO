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
#include "timeInterp.H"
#include "VelBCHolder.H"
#include "SetValLevel.H"

// ---------------------------------------------------------------
void
AMRNavierStokes::fillLambda(LevelData<FArrayBox>& a_lambda, Real a_time) const
{
  CH_assert (!a_lambda.isDefined());

  IntVect ghostVect(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));
  const DisjointBoxLayout& levelGrids = m_vel_new_ptr->getBoxes();
  int ncomp = 1;

  a_lambda.define(levelGrids, ncomp, ghostVect);

  if (s_set_bogus_values)
    {
      setValLevel(a_lambda, s_bogus_value);
    }

  Interval lambdaComps(0, ncomp-1);

  Real old_time = m_time - m_dt;
  // first copy lambda to a_lambda
  if (abs(a_time - old_time) < TIME_EPS)
    {
      // copy from old values
      m_lambda_old_ptr->copyTo(lambdaComps, a_lambda, lambdaComps);
    }
  else if (abs(a_time - m_time) < TIME_EPS)
    {
      // copy from new values
      m_lambda_new_ptr->copyTo(lambdaComps, a_lambda, lambdaComps);
    }
  else
    {
      // interp in time
      timeInterp(a_lambda, a_time, oldLambda(), old_time,
                 newLambda(), m_time, lambdaComps);
    }

  if (m_level > 0)
    {
      // interpolate boundary values from coarse grids
      AMRNavierStokes* crseLevelPtr = crseNSPtr();
      CH_assert (crseLevelPtr != NULL);
      LevelData<FArrayBox>& crse_old_lambda = crseLevelPtr->oldLambda();
      LevelData<FArrayBox>& crse_new_lambda = crseLevelPtr->newLambda();
      const DisjointBoxLayout& crseGrids = crse_new_lambda.getBoxes();
      const ProblemDomain& crseDomain = crseLevelPtr->problemDomain();
      int nRefCrse = crseLevelPtr->refRatio();

      Real crse_dt = crseLevelPtr->dt();
      // crse time should be the new time
      Real crse_old_time = crseLevelPtr->time() - crse_dt;
      Real time_interp_coeff;

      // check for "essentially 0 or 1"
      if (abs(a_time - crse_old_time) < TIME_EPS)
        {
          time_interp_coeff = 0.0;
        }
      else if (abs(a_time - crseLevelPtr->time()) < TIME_EPS)
        {
          time_interp_coeff = 1.0;
        }
      else
        {
          time_interp_coeff = (a_time - crse_old_time)/crse_dt;
        }

      PiecewiseLinearFillPatch filpatcher(levelGrids, crseGrids,
                                          ncomp, crseDomain, nRefCrse,
                                          ADVECT_GROW);

      filpatcher.fillInterp(a_lambda, crse_old_lambda, crse_new_lambda,
                            time_interp_coeff, 0, 0, ncomp);
    }

  const ProblemDomain& physDomain = problemDomain();

  // now do physical boundary conditions
  BCHolder lambdaBC = m_physBCPtr->lambdaFuncBC();
  // loop over boxes
  DataIterator dit = a_lambda.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      lambdaBC(a_lambda[dit], levelGrids[dit],
               physDomain, m_dx,
               false); // inhomogeneous
    }
}

// ---------------------------------------------------------------
void
AMRNavierStokes::fillVelocity(LevelData<FArrayBox>& a_vel, Real a_time)
{
  Interval velComps(0,SpaceDim-1);

  if (s_set_bogus_values)
  {
    setValLevel(a_vel, s_bogus_value);
  }

  Real old_time = m_time - m_dt;

  if (abs(a_time - old_time) < TIME_EPS)
    {
      m_vel_old_ptr->copyTo(velComps, a_vel, velComps);
    }
  else if (abs(a_time - m_time) < TIME_EPS)
    {
      m_vel_new_ptr->copyTo(velComps, a_vel, velComps);
    }
  else
    {
      // do linear interpolation in time
      timeInterp(a_vel, a_time, oldVel(), old_time,
                 newVel(), m_time, velComps);
    }

  // if necessary, do interpolation from coarser levels
  if (m_level > 0)
    {
      const DisjointBoxLayout& levelGrids = a_vel.getBoxes();
      const DisjointBoxLayout& thisLevelsGrids = newVel().getBoxes();

      const IntVect& velGrowVect = a_vel.ghostVect();
      int velGrow = velGrowVect[0];

      // if grids for a_vel are the same as those for this level, and
      // there are no ghost cells, we don't need to do this at all
      if (!( (velGrow == 0) && (levelGrids == thisLevelsGrids)))
        {
          // fill in coarse-fine BC data by conservative linear interp
          AMRNavierStokes& crseLevel = *crseNSPtr();

          LevelData<FArrayBox>& oldCrseVel = crseLevel.oldVel();
          LevelData<FArrayBox>& newCrseVel = crseLevel.newVel();
          const DisjointBoxLayout& crseGrids = oldCrseVel.getBoxes();
          const ProblemDomain& crseDomain = crseLevel.problemDomain();
          int nRefCrse = crseLevel.refRatio();

          Real crse_new_time = crseLevel.m_time;
          Real crse_dt = crseLevel.dt();
          Real crse_old_time = crse_new_time - crse_dt;
          Real crse_time_interp_coeff;

          // check for "essentially 0 or 1"
          if (abs(a_time - crse_old_time) < TIME_EPS)
            {
              crse_time_interp_coeff = 0.0;
            }
          else if (abs(a_time - crse_new_time) < TIME_EPS)
            {
              crse_time_interp_coeff = 1.0;
            }
          else
            {
              crse_time_interp_coeff = (a_time - crse_old_time)/crse_dt;
            }

          PiecewiseLinearFillPatch filpatcher(levelGrids, crseGrids,
                                              SpaceDim, crseDomain,
                                              nRefCrse, velGrow);

          filpatcher.fillInterp(a_vel, oldCrseVel, newCrseVel,
                                crse_time_interp_coeff, 0,0,SpaceDim);
        }
    }

  // need to set physical boundary conditions here
  const DisjointBoxLayout& levelGrids = a_vel.getBoxes();
  if (s_nu > 0.0)
    {
      VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());
      velBC.applyBCs(a_vel, levelGrids,
                     m_problem_domain, m_dx,
                     false); // inhomogeneous
    }
  else
    {
      VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
      velBC.applyBCs(a_vel, levelGrids,
                     m_problem_domain, m_dx,
                     false); // inhomogeneous
    }
  // note that there is no exchange here -- exchanges are the client's
  // problem (want to avoid extra exchanges)
}

// ---------------------------------------------------------------
void
AMRNavierStokes::fillVelocity(LevelData<FArrayBox>& a_vel, Real a_time,
                              int vel_comp, int dest_comp, int num_comp)
{
  Interval srcComps(vel_comp,vel_comp+num_comp-1);
  Interval destComps(dest_comp, dest_comp+num_comp-1);

  Real old_time = m_time - m_dt;
  // first do copy on this level
  if (abs(a_time - old_time) < TIME_EPS)
    {
      m_vel_old_ptr->copyTo(srcComps, a_vel, destComps);
    }
  else if (abs(a_time - m_time) < TIME_EPS)
    {
      m_vel_new_ptr->copyTo(srcComps, a_vel, destComps);
    }
  else
    {
      // do linear interpolation in time
      timeInterp(a_vel, a_time, oldVel(), old_time,
                 newVel(), m_time, srcComps,destComps);
    }

  // if necessary, do interpolation from coarser levels
  if (m_level > 0)
    {
      const DisjointBoxLayout& levelGrids = a_vel.getBoxes();
      const DisjointBoxLayout& thisLevelsGrids = newVel().getBoxes();

      const IntVect& velGrowVect = a_vel.ghostVect();
      int velGrow = velGrowVect[0];

      // if grids for a_vel are the same as those for this level, and
      // there are no ghost cells, we don't need to do this at all
      if (!( (velGrow == 0) && (levelGrids == thisLevelsGrids)))
        {
          // fill in coarse-fine BC data by conservative linear interp
          AMRNavierStokes& crseLevel = *crseNSPtr();

          LevelData<FArrayBox>& oldCrseVel = crseLevel.oldVel();
          LevelData<FArrayBox>& newCrseVel = crseLevel.newVel();
          const DisjointBoxLayout& crseGrids = oldCrseVel.getBoxes();
          const ProblemDomain& crseDomain = crseLevel.problemDomain();
          int nRefCrse = crseLevel.refRatio();

          Real crse_new_time = crseLevel.m_time;
          Real crse_dt = crseLevel.dt();
          Real crse_old_time = crse_new_time - crse_dt;
          Real crse_time_interp_coeff;

          // check for "essentially 0 or 1"
          if (abs(a_time - crse_old_time) < TIME_EPS)
            {
              crse_time_interp_coeff = 0.0;
            }
          else if (abs(a_time - crse_new_time) < TIME_EPS)
            {
              crse_time_interp_coeff = 1.0;
            }
          else
            {
              crse_time_interp_coeff = (a_time - crse_old_time)/crse_dt;
            }

          PiecewiseLinearFillPatch filpatcher(levelGrids, crseGrids,
                                              num_comp, crseDomain,
                                              nRefCrse, velGrow);

          filpatcher.fillInterp(a_vel, oldCrseVel, newCrseVel,
                                crse_time_interp_coeff, vel_comp, dest_comp,
                                num_comp);
        }
    }

  // need to set physical boundary conditions here
  const DisjointBoxLayout& levelGrids = a_vel.getBoxes();
  for (int comp = vel_comp; comp < vel_comp + num_comp; comp++)
    {
      int var = comp - vel_comp;
      Interval dataInterval(var, var);
      LevelData<FArrayBox> velComp;
      aliasLevelData(velComp, &a_vel, dataInterval);
      BCHolder velBC;
      if (s_nu > 0.0)
        {
          velBC = BCHolder(m_physBCPtr->viscousSolveFuncBC(comp));
        }
      else
        {
          velBC = BCHolder(m_physBCPtr->tracingSolveFuncBC(comp));
        }
      DataIterator dit = a_vel.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          velBC(velComp[dit], levelGrids[dit],
                m_problem_domain, m_dx,
                false); // inhomogeneous
        }
    }
}

// ------------------------------------------------------------------
void
AMRNavierStokes::velocity(LevelData<FArrayBox>& a_vel, Real a_time) const
{
  Interval velComps(0,SpaceDim-1);

  if (abs(a_time-m_time) < TIME_EPS)
    {
      m_vel_new_ptr->copyTo(velComps, a_vel, velComps);
    }
  else if (abs(a_time-(m_time-m_dt)) < TIME_EPS)
    {
      m_vel_old_ptr->copyTo(velComps, a_vel, velComps);
    }
  else
    {
      Real old_time = m_time - m_dt;
      timeInterp(a_vel, a_time, *m_vel_old_ptr, old_time,
                 *m_vel_new_ptr, m_time, velComps);
    }
}

// ------------------------------------------------------------------
void
AMRNavierStokes::lambda(LevelData<FArrayBox>& a_lambda, Real a_time) const
{
  Interval lambdaComps(0,0);

  if (abs(a_time-m_time) < TIME_EPS)
    {
      m_lambda_new_ptr->copyTo(lambdaComps, a_lambda, lambdaComps);
    }
  else if (abs(a_time-(m_time-m_dt)) < TIME_EPS)
    {
      m_lambda_old_ptr->copyTo(lambdaComps, a_lambda, lambdaComps);
    }
  else
    {
      Real old_time = m_time - m_dt;
      timeInterp(a_lambda, a_time, *m_lambda_old_ptr, old_time,
                 *m_lambda_new_ptr, m_time, lambdaComps);
    }
}

// ---------------------------------------------------------------
void
AMRNavierStokes::fillScalars(LevelData<FArrayBox>& a_scal,
                             Real a_time, const int a_comp) const
{
  // this will handle allocation as well...
  CH_assert (!a_scal.isDefined());

  IntVect ghostVect(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));
  const DisjointBoxLayout& levelGrids = m_vel_new_ptr->getBoxes();

  a_scal.define(levelGrids, m_scal_new[a_comp]->nComp(), ghostVect);

  if (s_set_bogus_values)
    {
      setValLevel(a_scal, s_bogus_value);
    }

  Interval scalComps = a_scal.interval();

  Real old_time = m_time - m_dt;

  // first do copy on this level
  if (abs(a_time - old_time) < TIME_EPS)
    {
    m_scal_old[a_comp]->copyTo(scalComps, a_scal, scalComps);
    }
  else if (abs(a_time - m_time) < TIME_EPS)
    {
      m_scal_new[a_comp]->copyTo(scalComps, a_scal, scalComps);
    }
  else
    {
      // do linear interpolation in time
      timeInterp(a_scal, a_time, oldScal(a_comp), old_time,
                 newScal(a_comp), m_time, scalComps);
    }

  if (m_level > 0)
    {
      // fill in coarse-fine BC data by conservative linear interp
      AMRNavierStokes& crseLevel = *crseNSPtr();

      LevelData<FArrayBox>& oldCrseScal = crseLevel.oldScal(a_comp);
      LevelData<FArrayBox>& newCrseScal = crseLevel.newScal(a_comp);
      const DisjointBoxLayout& crseGrids = oldCrseScal.getBoxes();
      const ProblemDomain& crseDomain = crseLevel.problemDomain();
      int nRefCrse = crseLevel.refRatio();

      Real crse_new_time = crseLevel.m_time;
      Real crse_dt = crseLevel.dt();
      Real crse_old_time = crse_new_time - crse_dt;
      Real crse_time_interp_coeff;

      // check for "essentially 0 or 1"
      if (abs(a_time - crse_old_time) < TIME_EPS)
        {
          crse_time_interp_coeff = 0.0;
        }
      else if (abs(a_time - crse_new_time) < TIME_EPS)
        {
          crse_time_interp_coeff = 1.0;
        }
      else
        {
          crse_time_interp_coeff = (a_time - crse_old_time)/crse_dt;
        }

      const IntVect& scalGrowVect = a_scal.ghostVect();
      int scalGrow = scalGrowVect[0];

      PiecewiseLinearFillPatch filpatcher(levelGrids, crseGrids,
                                          a_scal.nComp(), crseDomain,
                                          nRefCrse, scalGrow);

      filpatcher.fillInterp(a_scal, oldCrseScal, newCrseScal,
                            crse_time_interp_coeff, 0,0,a_scal.nComp());
    }

  const ProblemDomain& physDomain = problemDomain();

  BCHolder thisBC = m_physBCPtr->scalarTraceFuncBC(a_comp);
  // loop over boxes
  DataIterator dit = a_scal.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      thisBC(a_scal[dit], levelGrids[dit],
             physDomain, m_dx,
             false); // inhomogeneous
    }
}
