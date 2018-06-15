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
#include "AdvectUtil.H"
#include "CellToEdge.H"
#include "CH_Timer.H"

#include "AdvectionPhysics.H"
#include "VelBCHolder.H"
#include "SetValLevel.H"

// ---------------------------------------------------------------
// advance solution by one timestep
// return max safe timestep
Real
AMRNavierStokes::advance()
{
  if (s_verbosity >= 2)
  {
    pout () << "AMRNavierStokes::advance " << m_level
            << ", starting time = "
            << setiosflags(ios::fixed) << setprecision(6)
            << setw(12) << m_time
            << ", dt = "
            << setiosflags(ios::fixed) << setprecision(6) << setw(12) << dt()
            << endl;
  }

#ifdef DEBUG
  // if we're at the beginning of a composite timestep,
  // save current multilevel state
  if (m_level == 0)
    {
      AMRNavierStokes* levelPtr = this;
      // figure out finest level
      while (!levelPtr->finestLevel())
        {
          levelPtr = levelPtr->finerNSPtr();
        }
    int finest_level = levelPtr->m_level;

    levelPtr = this;
    for (int lev=0; lev<= finest_level; lev++)
      {
        levelPtr->m_saved_time = m_time;

        levelPtr->newVel().copyTo(levelPtr->newVel().interval(),
                                  *levelPtr->m_vel_save_ptr,
                                  levelPtr->m_vel_save_ptr->interval());

        levelPtr->newLambda().copyTo(levelPtr->newLambda().interval(),
                                     *levelPtr->m_lambda_save_ptr,
                                     levelPtr->m_lambda_save_ptr->interval());

        for (int comp=0; comp<s_num_scal_comps; comp++)
          {
            levelPtr->newScal(comp).copyTo(levelPtr->newScal(comp).interval(),
                                           *levelPtr->m_scal_save[comp],
                                           levelPtr->m_scal_save[comp]->interval());
          }

        levelPtr = levelPtr->finerNSPtr();
      }
    }
#endif

  Real old_time = m_time;
  Real new_time = m_time + m_dt;
  m_dt_save = m_dt;
  swapOldAndNewStates();

  const DisjointBoxLayout levelGrids = newVel().getBoxes();

  // initialize flux registers
  if (!finestLevel())
    {
      m_flux_register.setToZero();
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          m_scal_fluxreg_ptrs[comp]->setToZero();
        }
      m_lambda_flux_reg.setToZero();
    }

  // compute advection velocities -- if we're using the
  // patchGodunov approach, these will have ghost cells.
  IntVect advVelGhost(IntVect::Unit);
  LevelData<FluxBox> advectionVel(levelGrids, 1, advVelGhost);

  if (s_set_bogus_values)
    {
      setValLevel(advectionVel, s_bogus_value);
    }

  computeAdvectionVelocities(advectionVel);

  // now that advection velocities have been computed,
  // update advected scalars

  // do scalar boundary conditions on lambda
  // (also grow to appropriate size)

  LevelData<FArrayBox> lambdaTemp;

  fillLambda(lambdaTemp, old_time);

  LevelFluxRegister* crseLambdaFluxRegPtr = NULL;
  if (m_level > 0)
    {
      crseLambdaFluxRegPtr = &(crseNSPtr()->m_lambda_flux_reg);
    }

  advectScalar(newLambda(), lambdaTemp, advectionVel,
               crseLambdaFluxRegPtr, m_lambda_flux_reg,
               m_patchGodLambda,
               m_dt);

  // now do advected-diffused scalars
  if (s_num_scal_comps > 0)
    {
      // loop over scalar components
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          BCHolder scalPhysBCs = m_physBCPtr->scalarTraceFuncBC(comp);

          LevelData<FArrayBox> scalarTemp;
          fillScalars(scalarTemp, old_time,comp);

          // now get coarse-level scalars for BC's if necessary
          LevelData<FArrayBox>* newCrseScalPtr = NULL;
          LevelData<FArrayBox>* oldCrseScalPtr = NULL;
          Real oldCrseTime = -s_bogus_value;
          Real newCrseTime = s_bogus_value;
          LevelFluxRegister* crseScalFluxRegPtr = NULL;

          if (m_level > 0)
            {
              newCrseScalPtr = &(crseNSPtr()->newScal(comp));
              oldCrseScalPtr = &(crseNSPtr()->oldScal(comp));
              newCrseTime = crseNSPtr()->time();
              oldCrseTime = newCrseTime - crseNSPtr()->dt();
              crseScalFluxRegPtr = (crseNSPtr()->m_scal_fluxreg_ptrs[comp]);
            }

          advectDiffuseScalar(newScal(comp), scalarTemp,
                              advectionVel,
                              s_scal_coeffs[comp],
                              oldCrseScalPtr, newCrseScalPtr,
                              oldCrseTime, newCrseTime,
                              crseScalFluxRegPtr,
                              *(m_scal_fluxreg_ptrs[comp]),
                              *(m_patchGodScalars[comp]),
                              scalPhysBCs,
                              m_dt, comp);
        } // end loop over components
    } // end advected-diffused scalars

  // now predict velocities -- this returns the advection term
  // put this in "new_vel"
  LevelData<FArrayBox>& uStar = newVel();
  predictVelocities(uStar, advectionVel);

  computeUStar(uStar);

  // if a coarser level exists, will need coarse-level data for proj
  LevelData<FArrayBox>* crseVelPtr = NULL;

  if (m_level > 0)
    {
      const DisjointBoxLayout& crseGrids = crseNSPtr()->newVel().getBoxes();
      crseVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
      // coarse velocity BC data is interpolated in time
      crseNSPtr()->fillVelocity(*crseVelPtr, new_time);
    }

  // need to do physical boundary conditions and exchanges
  bool isViscous = (s_nu > 0);
  VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));
  velBC.applyBCs(uStar, levelGrids,
                 m_problem_domain, m_dx,
                 false); // inhomogeneous

  // noel -- all time in level project
  m_projection.LevelProject(uStar, crseVelPtr, new_time, m_dt);

  // as things stand now, physical BC's are re-set in LevelProjection

  // compute maximum safe timestep for next iteration
  Real newDt = computeDt();

  // clean up temp storage
  if (crseVelPtr != NULL)
    {
      delete crseVelPtr;
      crseVelPtr = NULL;
    }

  ++m_level_steps;

  if (s_verbosity >= 2)
    {

      pout () << "AMRNavierStokes::advance " << m_level
              << ",      end time = "
              << setiosflags(ios::fixed) << setprecision(6) << setw(12) << m_time
              << ", dt = "
              << setiosflags(ios::fixed) << setprecision(6)
              << setw(12) << dt()
              << endl;
    }

  return newDt;
}

// ---------------------------------------------------------------
void
AMRNavierStokes::computeAdvectionVelocities(LevelData<FluxBox>& a_advVel)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::computeAdvectionVelocities: "
             << m_level << endl;
    }

  bool isViscous = (s_nu > 0.0);

  const DisjointBoxLayout& levelGrids = newVel().getBoxes();

  /// need to build grown grids to get be able to do all of
  /// tracing properly
  IntVect advect_grow(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));

  LevelData<FArrayBox> old_vel(levelGrids, SpaceDim, advect_grow);
  LevelData<FArrayBox> viscousSource(levelGrids, SpaceDim, advect_grow);
  LevelData<FArrayBox>* crseVelPtr = NULL;

  if (s_set_bogus_values)
    {
      setValLevel(old_vel, s_bogus_value);
      setValLevel(viscousSource, s_bogus_value);
    }

  // m_time contains the time at which the new state is centered
  Real old_time = m_time - m_dt;
  fillVelocity(old_vel, old_time);

  // set physical boundary conditions here
  // set physical boundary conditions on velocity

  if (isViscous)
    {
      LevelData<FArrayBox> viscousVel(levelGrids, SpaceDim, advect_grow);
      DataIterator dit = viscousVel.dataIterator();
      // rather than resetting BC's on old_vel, and then setting them
      // back, just copy old_vel to a temporary holder, and set
      // BCs on that one.
      for (dit.begin(); dit.ok(); ++dit)
        {
          viscousVel[dit].copy(old_vel[dit]);
        }
      VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());
      velBC.applyBCs(viscousVel, levelGrids,
                     m_problem_domain, m_dx,
                     false); // inhomogeneous

      // if crse level exists, fill coarse velocity BC
      if (m_level > 0)
        {
          const DisjointBoxLayout& crseGrids = crseNSPtr()->newVel().getBoxes();
          crseVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
          crseNSPtr()->fillVelocity(*crseVelPtr, old_time);
        }

      computeLapVel(viscousSource, viscousVel, crseVelPtr);

      for (dit.reset(); dit.ok(); ++dit)
        {
          viscousSource[dit].mult(s_nu);
        }
    }
  else
    {
      setValLevel(viscousSource, 0.0);
    }

  // tracing will use inviscid BC's
  {
    VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
    velBC.applyBCs(old_vel, levelGrids,
                   m_problem_domain, m_dx,
                   false); // inhomogeneous
  }

  // call utility function to do tracing
  traceAdvectionVel(a_advVel, old_vel, viscousSource,
                    m_patchGodVelocity, old_time, m_dt);

  EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectionVelFuncBC(isViscous));
  edgeVelBC.applyBCs(a_advVel, levelGrids,
                     m_problem_domain, m_dx,
                     false); // inhomogeneous

  // noel  levelMacProject is big guy
  // now MAC project
  m_projection.levelMacProject(a_advVel, old_time, m_dt);

  // finally, add volume discrepancy correction
  if (m_projection.etaLambda() > 0 && s_applyFreestreamCorrection)
    {
      LevelData<FluxBox>& grad_e_Lambda = m_projection.grad_eLambda();

      DataIterator dit = levelGrids.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          FluxBox& thisGrad_eLambda = grad_e_Lambda[dit];
          FluxBox& thisAdvVel = a_advVel[dit];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              thisAdvVel[dir] += thisGrad_eLambda[dir];
            }
        }
    }

  edgeVelBC.applyBCs(a_advVel, levelGrids,
                     m_problem_domain, m_dx,
                     false); // inhomogeneous

  // clean up storage
  if (crseVelPtr != NULL)
    {
      delete crseVelPtr;
      crseVelPtr = NULL;
    }
}
