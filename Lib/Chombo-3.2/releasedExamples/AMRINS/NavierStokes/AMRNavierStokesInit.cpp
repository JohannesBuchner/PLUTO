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
#include "FABView.H"
#include "Gradient.H"

#include "probF_F.H"
#include "AMRNSF_F.H"

#include "AdvectionPhysics.H"
#include "VelBCHolder.H"

// ---------------------------------------------------------------
// initialize grid
void
AMRNavierStokes::initialGrid(const Vector<Box>& a_new_grids)
{
  if (s_verbosity >= 3)
  {
    pout () << "AMRNavierStokes::initialGrid " << m_level << endl;
  }
  m_level_grids = a_new_grids;
  if (a_new_grids.size() > 0)
  {
    m_is_empty = false;
  }

  // note -- this assumes levels are defined from coarse->fine
  // (which is a reasonable assumption, and just happens to be
  // the way things _are_ implemented)
  if (!m_is_empty)
    {
      finestLevel(true);
      if (m_coarser_level_ptr != NULL)
        {
          AMRNavierStokes* crse_amrns_ptr = dynamic_cast<AMRNavierStokes*> (m_coarser_level_ptr);
          if (crse_amrns_ptr != NULL)
            {
              crse_amrns_ptr->finestLevel(false);
            }
          else
            {
              MayDay::Error("in AMRNavierStokes::initialGrid: m_coarser_level is not castable to AMRNavierStokes*");
            }
        }
    }

  const DisjointBoxLayout level_domain = loadBalance(a_new_grids);

  if (s_verbosity >= 4)
    {
      pout () << "new grids: " << endl;
      for (LayoutIterator lit = level_domain.layoutIterator(); lit.ok(); ++lit)
        {
          pout() << level_domain[lit] << endl;
        }
    }

  IntVect ghostVect(D_DECL(1,1,1));

  // this is necessary because of the way initialization
  // is done -- AMRNavierStokes::initialGrid is called each time
  // through level-by-level building of the initial grid
  // hierarchy

  if (m_vel_new_ptr != NULL)
    {
      delete m_vel_new_ptr;
      m_vel_new_ptr = NULL;
    }
  m_vel_new_ptr = new LevelData<FArrayBox>(level_domain,
                                           s_num_vel_comps,
                                           ghostVect);

  if (m_vel_old_ptr != NULL)
    {
      delete m_vel_old_ptr;
      m_vel_old_ptr = NULL;
    }
  m_vel_old_ptr = new LevelData<FArrayBox>(level_domain,
                                           s_num_vel_comps,
                                           ghostVect);

#ifdef DEBUG
  if (m_vel_save_ptr != NULL)
    {
      delete m_vel_save_ptr;
      m_vel_save_ptr = NULL;
    }
  m_vel_save_ptr = new LevelData<FArrayBox>(level_domain,
                                            s_num_vel_comps,
                                            ghostVect);
#endif

  if (m_lambda_new_ptr != NULL)
    {
      delete m_lambda_new_ptr;
      m_lambda_new_ptr = NULL;
    }
  m_lambda_new_ptr = new LevelData<FArrayBox>(level_domain,
                                              1, ghostVect);

  if (m_lambda_old_ptr != NULL)
    {
      delete m_lambda_old_ptr;
      m_lambda_old_ptr = NULL;
    }
  m_lambda_old_ptr = new LevelData<FArrayBox>(level_domain,
                                              1, ghostVect);

#ifdef DEBUG
  if (m_lambda_save_ptr != NULL)
    {
      delete m_lambda_save_ptr;
      m_lambda_save_ptr = NULL;
    }
  m_lambda_save_ptr = new LevelData<FArrayBox>(level_domain,
                                               1,
                                               ghostVect);
#endif

  // now do scalars
  if (m_scal_new.size() != s_num_scal_comps)
    {
      //  initialize pointers to null
      Vector<LevelData<FArrayBox>* > tempVect(s_num_scal_comps,NULL);
      m_scal_new = tempVect;
    }
  if (m_scal_old.size() != s_num_scal_comps)
    {
      Vector<LevelData<FArrayBox>* > tempVect(s_num_scal_comps,NULL);
      m_scal_old = tempVect;
    }

#ifdef DEBUG
  if (m_scal_save.size() != s_num_scal_comps)
    {
      Vector<LevelData<FArrayBox>* > tempVect(s_num_scal_comps,NULL);
      m_scal_save = tempVect;
    }
#endif

  // do flux registers as well
  if (m_scal_fluxreg_ptrs.size() < s_num_scal_comps)
    {
      m_scal_fluxreg_ptrs.resize(s_num_scal_comps,NULL);
    }

  for (int comp=0; comp<s_num_scal_comps; comp++)
    {
      if (m_scal_new[comp] != NULL)
        {
          delete m_scal_new[comp];
          m_scal_new[comp] = NULL;
        }
      m_scal_new[comp] = new LevelData<FArrayBox>(level_domain,
                                                  s_compsPerScalar,ghostVect);

      if (m_scal_old[comp] != NULL)
        {
          delete m_scal_old[comp];
          m_scal_old[comp] = NULL;
        }
      m_scal_old[comp] = new LevelData<FArrayBox>(level_domain,
                                                  s_compsPerScalar, ghostVect);

#ifdef DEBUG
      if (m_scal_save[comp] != NULL)
        {
          delete m_scal_save[comp];
          m_scal_save[comp] = NULL;
        }
      m_scal_save[comp] = new LevelData<FArrayBox>(level_domain,
                                                   s_compsPerScalar, ghostVect);
#endif

      if (m_scal_fluxreg_ptrs[comp] != NULL)
        {
          delete m_scal_fluxreg_ptrs[comp];
          m_scal_fluxreg_ptrs[comp] = NULL;
        }
      m_scal_fluxreg_ptrs[comp] = new LevelFluxRegister;
    } // end loop over scalar components

  setAllBogus();

  if (SpaceDim == 1)
    {
      m_dx = s_domLength[0]/problemDomain().domainBox().size(0);
    }
  else
    {
      m_dx = s_domLength[1]/problemDomain().domainBox().size(1);
    }

  // set up data structures -- actually, defer this until postInitialize
  //levelSetup(level_domain);
}

// ---------------------------------------------------------------
// set up data structures associated with this level
void
AMRNavierStokes::levelSetup(const DisjointBoxLayout& a_grids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::levelSetup " << m_level << endl;
    }

  CCProjector* crseProjPtr = NULL;
  CCProjector* fineProjPtr = NULL;
  int nRefCrse = -1;

  const DisjointBoxLayout* crseGridsPtr = NULL;

  if (m_coarser_level_ptr != NULL)
    {
      AMRNavierStokes* amrns_ptr = dynamic_cast<AMRNavierStokes*> (m_coarser_level_ptr);
      if (amrns_ptr != NULL)
        {
          crseGridsPtr = &(amrns_ptr->newVel().getBoxes());

          crseProjPtr = &(amrns_ptr->m_projection);
          nRefCrse = amrns_ptr->m_ref_ratio;

          m_coarse_average.define(a_grids, s_num_vel_comps, nRefCrse);
          m_coarse_average_lambda.define(a_grids, 1, nRefCrse);

          // flux registers
          amrns_ptr->m_flux_register.undefine();
          amrns_ptr->m_flux_register.define(a_grids,
                                            *crseGridsPtr,
                                            m_problem_domain,
                                            nRefCrse,
                                            s_num_vel_comps);

          amrns_ptr->m_lambda_flux_reg.undefine();
          amrns_ptr->m_lambda_flux_reg.define(a_grids,
                                              *crseGridsPtr,
                                              m_problem_domain,
                                              nRefCrse,
                                              1);

          for (int comp=0; comp<s_num_scal_comps; comp++)
            {
              amrns_ptr->m_scal_fluxreg_ptrs[comp]->undefine();
              amrns_ptr->m_scal_fluxreg_ptrs[comp]->define(a_grids,
                                                           *crseGridsPtr,
                                                           m_problem_domain,
                                                           nRefCrse,
                                                           1);
            }

          amrns_ptr->m_flux_register.setToZero();
          amrns_ptr->m_lambda_flux_reg.setToZero();
          for (int comp=0; comp<s_num_scal_comps; comp++)
            {
              amrns_ptr->m_scal_fluxreg_ptrs[comp]->setToZero();
            }

          // define coarse-fine interpolation operator
          m_velCFInterp.define(a_grids, crseGridsPtr,
                               m_dx, nRefCrse, SpaceDim,
                               m_problem_domain);

          // check to see if scalars are diffusive
          bool isDiffusive = false;
          for (int comp=0; comp<s_num_scal_comps; comp++)
            {
              if (s_scal_coeffs[comp] > 0) isDiffusive = true;
            }

          if (isDiffusive)
            {
              BCHolder bcHolder = m_physBCPtr->scalarTraceFuncBC(0);
              m_scalarsAMRPoissonOp.define(a_grids, crseGridsPtr,
                                           m_dx, nRefCrse, m_problem_domain,
                                           bcHolder); // will be replaced
            }

          // define viscous solver here as well
          if ((s_nu > 0.0) || isDiffusive)
            {
              BCHolder bcHolder = m_physBCPtr->viscousFuncBC();
              m_velocityAMRPoissonOp.define(a_grids, crseGridsPtr,
                                            m_dx, nRefCrse, m_problem_domain,
                                            bcHolder); // will be replaced
              defineViscousMGSolver(a_grids, crseGridsPtr, nRefCrse);
            }

        }
      else
        {
          MayDay::Error("in AMRNavierStokes::levelSetup: m_coarser_level_ptr is not castable to AMRNavierStokes*");
        }

    }
  else
    {
      // if no coarser level, still need to define viscous solver
      // check to see if scalars are diffusive
      bool isDiffusive = false;
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          if (s_scal_coeffs[comp] > 0) isDiffusive = true;
        }

      if (isDiffusive)
        {
          BCHolder bcHolder = m_physBCPtr->scalarTraceFuncBC(0);
          m_scalarsAMRPoissonOp.define(a_grids, crseGridsPtr,
                                       m_dx, nRefCrse, m_problem_domain,
                                       bcHolder); // will be replaced
        }

      if ((s_nu > 0.0) || isDiffusive)
        {
          BCHolder bcHolder = m_physBCPtr->viscousFuncBC();
          m_velocityAMRPoissonOp.define(a_grids, crseGridsPtr,
                                        m_dx, nRefCrse, m_problem_domain,
                                        bcHolder); // will be replaced

          defineViscousMGSolver(a_grids, crseGridsPtr, nRefCrse);
        } // end if viscous


    } // end if no coarser level exists

  if (m_finer_level_ptr != NULL)
    {
      AMRNavierStokes* amrns_ptr = dynamic_cast<AMRNavierStokes*> (m_finer_level_ptr);
      if (amrns_ptr != NULL)
        {
          fineProjPtr = &(amrns_ptr->m_projection);
        }
      else
        {
          MayDay::Error("in AMRNavierStokes::levelSetup: m_finer_level_ptr is not castable to AMRNavierStokes*");
        }
    }
  else
    {
      m_flux_register.undefine();
    }

  // pass on verbosity to projection
  m_projection.verbosity(s_verbosity);
  // pass on whether to limit coarsenings
  m_projection.limitSolverCoarsening(m_limitSolverCoarsening);

  m_projection.define(a_grids, crseGridsPtr, m_problem_domain,
                      m_dx, fineProjPtr,
                      crseProjPtr, nRefCrse, m_level, *m_physBCPtr);

  // 1 -> PLM, 2 -> PPM
  int normalPredOrder = 1;

  // Use 4th order slope computations
  bool useFourthOrderSlopes = true;

  // Don't do slope limiting
  bool usePrimLimiting = false;
  bool useCharLimiting = false;

  // Don't do flattening
  bool useFlattening = false;

  // No artificial viscosity
  bool useArtVisc = false;
  Real artVisc = -1.0;

  // Define objects for the advection of lambda
  AdvectionPhysics advectionPhysicsLambda;

  advectionPhysicsLambda.define(m_problem_domain,m_dx);
  advectionPhysicsLambda.setNComp(1);

  PhysIBC* lambdaIBC = m_physBCPtr->lambdaTraceIBC();
  advectionPhysicsLambda.setPhysIBC(lambdaIBC);

  m_patchGodLambda.define(m_problem_domain,
                          m_dx,
                          &advectionPhysicsLambda,
                          normalPredOrder,
                          useFourthOrderSlopes,
                          usePrimLimiting,
                          useCharLimiting,
                          useFlattening,
                          useArtVisc,
                          artVisc);

  // Define objects for the advection of the velocity
  AdvectionPhysics advectionPhysicsVelocity;

  advectionPhysicsVelocity.define(m_problem_domain,m_dx);
  advectionPhysicsVelocity.setNComp(SpaceDim);

  PhysIBC* velocityIBC = m_physBCPtr->advectionVelIBC();
  advectionPhysicsVelocity.setPhysIBC(velocityIBC);

  m_patchGodVelocity.define(m_problem_domain,
                            m_dx,
                            &advectionPhysicsVelocity,
                            normalPredOrder,
                            useFourthOrderSlopes,
                            usePrimLimiting,
                            useCharLimiting,
                            useFlattening,
                            useArtVisc,
                            artVisc);

  // Define objects for the advection of the scalars

  m_patchGodScalars.resize(s_num_scal_comps);

  for (int comp = 0; comp < s_num_scal_comps; comp++)
  {
    AdvectionPhysics advectionPhysicsScalars;

    advectionPhysicsScalars.define(m_problem_domain,m_dx);
    advectionPhysicsScalars.setNComp(s_compsPerScalar);

    PhysIBC* scalarIBC = m_physBCPtr->scalarTraceIBC(comp);
    advectionPhysicsScalars.setPhysIBC(scalarIBC);

    m_patchGodScalars[comp] = new PatchGodunov;
    m_patchGodScalars[comp]->define(m_problem_domain,
                                    m_dx,
                                    &advectionPhysicsScalars,
                                    normalPredOrder,
                                    useFourthOrderSlopes,
                                    usePrimLimiting,
                                    useCharLimiting,
                                    useFlattening,
                                    useArtVisc,
                                    artVisc);
  }
}

// ---------------------------------------------------------------
void
AMRNavierStokes::readParameters()
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::readParameters " << endl;
    }

  ParmParse ppNS("ns");

  Real defaultDomainLength = 1.0;
  vector<Real> tempVect(SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      tempVect[dir] = defaultDomainLength;
    }
  ppNS.queryarr("domainLength", tempVect, 0, SpaceDim);
  s_domLength = tempVect;

  ppNS.query("init_shrink", s_init_shrink);

  ppNS.query("max_dt", s_max_dt);

  ppNS.query("setDt", s_prescribedDt);

  int temp = s_project_initial_vel;
  ppNS.query("project_initial_vel", temp);
  s_project_initial_vel = (temp == 1);

  temp = s_initialize_pressures;
  ppNS.query("init_pressures", temp);
  s_initialize_pressures = (temp == 1);

  temp = s_smooth_after_regrid;
  ppNS.query("smooth_after_regrid", temp);
  s_smooth_after_regrid = (temp == 1);

  ppNS.query("postRegrid_smoothing_coeff", s_regrid_smoothing_coeff);

  if (s_initialize_pressures)
    {
      ppNS.query("num_init_passes", s_num_init_passes);
    }

  temp = s_reflux_momentum;
  ppNS.query("reflux_momentum", temp);
  s_reflux_momentum = (temp == 1);

  temp = s_reflux_normal_momentum;
  ppNS.query("reflux_normal_momentum", temp);
  s_reflux_normal_momentum = (temp == 1);

  temp = s_tag_vorticity;
  ppNS.query("tag_vorticity", temp);
  s_tag_vorticity = (temp == 1);

  temp = s_tag_theta;
  ppNS.query("tag_theta", temp);
  s_tag_theta = (temp == 1);

  ppNS.query("vorticity_tagging_factor", s_vort_factor);

  ppNS.query("tags_grow", s_tags_grow);

  ppNS.query("viscosity", s_nu);

  temp = s_implicit_reflux;
  ppNS.query("implicit_reflux", temp);
  s_implicit_reflux = (temp == 1);

  temp = s_applyFreestreamCorrection;
  ppNS.query("applyFreestreamCorrection",temp);
  s_applyFreestreamCorrection = (temp == 1);

  temp = s_reflux_scal;
  ppNS.query("reflux_scalars", temp);
  s_reflux_scal = (temp == 1);

  temp = s_implicit_scal_reflux;
  ppNS.query("implicit_scal_reflux", temp);
  s_implicit_scal_reflux = (temp == 1);

  ppNS.query("viscous_solver_type", s_viscous_solver_type);

  ppNS.query("viscous_solver_tolerance", s_viscous_solver_tol);

  ppNS.query("viscous_num_smooth_up", s_viscous_num_smooth_up);
  ppNS.query("viscous_num_smooth_down", s_viscous_num_smooth_down);

  ppNS.query("num_scalars", s_num_scal_comps);

  if (s_num_scal_comps > 0)
    {
      s_scal_coeffs.resize(s_num_scal_comps);
      ppNS.getarr("scal_diffusion_coeffs", tempVect,
                  0,s_num_scal_comps);
      s_scal_coeffs = tempVect;
    }

  temp = s_specifyInitialGrids;
  ppNS.query("specifyInitialGrids", temp);
  s_specifyInitialGrids = (temp == 1);

  if (s_specifyInitialGrids)
    {
      ppNS.get("initialGridFile", s_initialGridFile);
    }

  temp = s_init_vel_from_vorticity;
  ppNS.query("initVelFromVorticity", temp);
  s_init_vel_from_vorticity = (temp == 1);
  // if 1-D, can't really intialize from vorticity
  if (SpaceDim ==1)
  {
    s_init_vel_from_vorticity = 0;
  }

  ppNS.query("backgroundVelocity", s_backgroundVel);

  temp = s_write_divergence;
  ppNS.query("writeDivergence", temp);
  s_write_divergence = (temp == 1);

  temp =  s_write_lambda;
  ppNS.query("writeLambda", temp);
  s_write_lambda = (temp == 1);

  temp = s_write_time_derivatives;
  ppNS.query("writeTimeDerivatives", temp);
  s_write_time_derivatives = (temp == 1);

  temp = s_write_vorticity;
  ppNS.query("writeVorticity", temp);
  s_write_vorticity = (temp == 1);

  temp = s_write_scalars;
  ppNS.query("writeScalars", temp);
  s_write_scalars = (temp == 1);

  /// include d(scalars)/dt in plotfile?
  temp = s_write_dScalar_dt;
  ppNS.query("writeDscalDt", temp);
  s_write_dScalar_dt = (temp == 1);

  /// include strain-rates in plotfile? (2D only)
  temp = s_write_strains;
  ppNS.query("writeStrains", temp);
  s_write_strains = (temp == 1);

  /// include freestream preservation correction in plotfile?
  temp = s_write_grad_eLambda;
  ppNS.query("writeGradELambda", temp);
  s_write_grad_eLambda = (temp == 1);

  /// include solution error (if there is an exact solution)
  temp = s_write_error;
  ppNS.query("writeError", temp);
  s_write_error = (temp == 1);

  // include processor id information in plotfile?
  temp = s_write_proc_ids;
  ppNS.query("writeProcIDs", temp);
  s_write_proc_ids = (temp == 1);

  // compute and report scalar error after composite timestep?
  temp = s_compute_scal_err;
  ppNS.query("computeScalErr", temp);
  s_compute_scal_err = (temp ==1);

  temp = s_write_grids;
  ppNS.query("writeGrids", temp);
  s_write_grids = (temp == 1);

  // now dump this all back out...
  if (s_verbosity > 2)
    {
      pout() << "NavierStokes inputs:" << endl;

      pout() << "  domainLength = " << s_domLength << endl;

      pout() << "  init_shrink = " <<  s_init_shrink << endl;

      pout() << "  max_dt = " << s_max_dt << endl;

      pout() << "  setDt = " << s_prescribedDt << endl;

      pout() << "  project_initial_vel"  << " = " << s_project_initial_vel
             << endl;

      pout() << "  init_pressures = " << s_initialize_pressures << endl;

      pout() << "  smooth_after_regrid = " << s_smooth_after_regrid
             << endl;

      pout() << "  postRegrid_smoothing_coeff = "
             << s_regrid_smoothing_coeff  << endl;

      if (s_initialize_pressures)
        {
          pout() << "  num_init_passes = " << s_num_init_passes << endl;
        }

      pout() << "  reflux_momentum = " << s_reflux_momentum << endl;

      pout() << "  reflux_normal_momentum = "
             << s_reflux_normal_momentum << endl;

      pout() << "  tag_vorticity = " << s_tag_vorticity << endl;

      pout() << "  tag_theta = " << s_tag_theta << endl;

      pout() << "  vorticity_tagging_factor = " << s_vort_factor << endl;

      pout() << "  tags_grow = " << s_tags_grow << endl;

      pout() << "  viscosity = " << s_nu << endl;

      pout() << "  implicit_reflux = " << s_implicit_reflux << endl;

      pout() << "  applyFreestreamCorrection = "
             << s_applyFreestreamCorrection << endl;

      pout() << "  reflux_scalars = " << s_reflux_scal << endl;

      pout() << "  implicit_scal_reflux = " << s_implicit_scal_reflux
             << endl;

      pout() << "  viscous_solver_type = " << s_viscous_solver_type
             << endl;

      pout() << "  viscous_solver_tolerance = " << s_viscous_solver_tol
             << endl;

      pout() << "  viscous_solver_num_smooth_up = "
             << s_viscous_num_smooth_up << endl;

      pout() << "  viscous_solver_num_smooth_down = "
             << s_viscous_num_smooth_down << endl;

      pout() << "  num_scalars = " << s_num_scal_comps << endl;

      if (s_num_scal_comps > 0)
        {
          pout () << "  scal_diffusion_coeffs = " << s_scal_coeffs << endl;
        }

      pout() << "  specifyInitialGrids = " << s_specifyInitialGrids
             << endl;

      if (s_specifyInitialGrids)
        {
          pout() << "  initialGridFile = " << s_initialGridFile << endl;
        }

      pout() << "  initVelFromVorticity = " << s_init_vel_from_vorticity
             << endl;

      pout() << "  backgroundVelocity = " << s_backgroundVel << endl;

      pout() << "  writeDivergence = " << s_write_divergence << endl;

      pout() << "  writeLambda = " << s_write_lambda << endl;

      pout() << "  writeTimeDerivatives = " << s_write_time_derivatives
             << endl;

      pout() << "  writeVorticity = " << s_write_vorticity << endl;

      pout() << "  writeScalars = " << s_write_scalars << endl;

      /// include d(scalars)/dt in plotfile?
      pout() << "  writeDscalDt = " << s_write_dScalar_dt << endl;

      /// include strain-rates in plotfile? (2D only)
      pout() << "  writeStrains = " << s_write_strains << endl;

      /// include freestream preservation correcion in plotfile?
      pout() << "  writeGradELambda = " << s_write_grad_eLambda << endl;

      /// include solution error (if there is an exact solution)
      pout() << "  writeError = " << s_write_error  << endl;

      // include processor id information in plotfile?
      pout() << "  writeProcIDs = " << s_write_proc_ids << endl;

      pout() << "  compute scalar error = " << s_compute_scal_err << endl;
    }

  // set flag to indicate that we've done this
  s_ppInit  =true;

}

// ---------------------------------------------------------------
void
AMRNavierStokes::initialData()
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::initialData " << m_level << endl;
    }

  DataIterator dit = newVel().dataIterator();

  if (s_init_vel_from_vorticity)
    {
      // first initialize vorticity field -- do this in an AMR way
      if (finestLevel())
        {
          Vector<LevelData<FArrayBox>*> initialVort(m_level+1, NULL);
          Vector<LevelData<FArrayBox>*> streamFunc(m_level+1, NULL);
          Vector<DisjointBoxLayout> amrGrids(m_level+1);
          Vector<ProblemDomain> amrDomains(m_level+1);
          Vector<Real> amrDx(m_level+1);
          Vector<int> refRatios(m_level+1);

          int nVortComp;
          if (SpaceDim == 2)
            {
              nVortComp = 1;
            }
          else
            {
              nVortComp = SpaceDim;
            }

          // start at finest (this) level and progress to coarser levels
          AMRNavierStokes* thisLevelPtr = this;
          for (int lev=m_level; lev>=0; lev--)
            {
              const DisjointBoxLayout& levelGrids
                = thisLevelPtr->newVel().getBoxes();

              amrGrids[lev] = levelGrids;
              initialVort[lev] = new LevelData<FArrayBox>(levelGrids,
                                                          nVortComp,
                                                          IntVect::Zero);

              streamFunc[lev] = new LevelData<FArrayBox>(levelGrids, nVortComp,
                                                         IntVect::Unit);

              amrDomains[lev] = thisLevelPtr->problemDomain();
              Real levelDx = thisLevelPtr->Dx();
              amrDx[lev] = levelDx;
              refRatios[lev] = thisLevelPtr->refRatio();

              DataIterator levelDit = levelGrids.dataIterator();
              for (levelDit.reset(); levelDit.ok(); ++levelDit)
                {
                  // since there is only one component of vorticity
                  // in this case, put it in the 0th component
                  for (int comp=0; comp<nVortComp; comp++)
                    {
                      FORT_INITVORT(CHF_FRA1((*initialVort[lev])[levelDit],
                                             comp),
                                    CHF_INT(comp),
                                    CHF_REAL(levelDx));
                    }
                }

              // index to next coarser level
              if (lev > 0)
                {
                  thisLevelPtr = thisLevelPtr->crseNSPtr();
                }
            }

          int numLevels = m_level+1;

          AMRPoissonOpFactory localPoissonOpFactory;
          localPoissonOpFactory.define(amrDomains[0],
                                       amrGrids,
                                       refRatios,
                                       amrDx[0],
                                       m_physBCPtr->streamBC());

          RelaxSolver<LevelData<FArrayBox> > bottomSolver;
          bottomSolver.m_verbosity = s_verbosity;

          AMRMultiGrid<LevelData<FArrayBox> > streamSolver;
          AMRLevelOpFactory<LevelData<FArrayBox> >& downCastFact = (AMRLevelOpFactory<LevelData<FArrayBox> >& ) (localPoissonOpFactory);
          streamSolver.define(amrDomains[0],
                              downCastFact,
                              &bottomSolver,
                              numLevels);
          streamSolver.m_verbosity = s_verbosity;
          streamSolver.m_eps = 1e-10;

          if (SpaceDim == 2)
            {
              streamSolver.solve(streamFunc, initialVort, numLevels-1, 0,
                                 true); // initialize streamFunc (all levels) to zero

              thisLevelPtr = this;
              // reapply BCs and compute velocity field
              CH_assert (thisLevelPtr->m_level == this->m_level);

              Interval psiComps(0,0);
              BCFunc psiBCFunc = m_physBCPtr->streamBC();

              for (int lev=m_level; lev>=0; lev--)
                {
                  LevelData<FArrayBox> gradPsi(amrGrids[lev],SpaceDim,
                                               IntVect::Zero);

                  // do averaging down if necessary
                  if (lev < m_level)
                    {
                      LevelData<FArrayBox>& finePsi = *streamFunc[lev+1];
                      const DisjointBoxLayout& fineGrids = finePsi.getBoxes();
                      CoarseAverage PsiAvgDown(fineGrids, finePsi.nComp(),
                                               refRatios[lev]);
                      PsiAvgDown.averageToCoarse(*streamFunc[lev],
                                                 finePsi);
                    }

                  // apply physical BC's to stream function - inhomogeneous BC
                  bool isHomogeneous = false;
                  DataIterator levelDit = streamFunc[lev]->dataIterator();
                  for (levelDit.reset(); levelDit.ok(); ++levelDit)
                    {
                      const Box& streamFuncBox = (amrGrids[lev])[levelDit];
                      FArrayBox& streamFuncFAB = (*streamFunc[lev])[levelDit];

                      psiBCFunc(streamFuncFAB,
                                streamFuncBox,
                                amrDomains[lev],
                                amrDx[lev],
                                isHomogeneous);
                    }

                  // now compute derivatives of psi
                  LevelData<FArrayBox>* crsePsiPtr = NULL;
                  LevelData<FArrayBox>* finePsiPtr = NULL;
                  int nRefCrse=-1;
                  int nRefFine=-1;

                  if (lev > 0)
                    {
                      crsePsiPtr = streamFunc[lev-1];
                      nRefCrse = refRatios[lev-1];
                    }

                  if (lev < m_level)
                    {
                      finePsiPtr = streamFunc[lev+1];
                      nRefFine = refRatios[lev];
                    }

                  Gradient::compGradientCC(gradPsi, *streamFunc[lev],
                                           crsePsiPtr, finePsiPtr,
                                           amrDx[lev], nRefCrse, nRefFine,
                                           amrDomains[lev]);

                  // now copy derivatives of psi into the proper
                  // velocity components
                  Interval xVelComp(0,0);
                  Interval yVelComp(1,1);
                  Interval xDerivativeComp(0,0);
                  Interval yDerivativeComp(1,1);

                  LevelData<FArrayBox>& levelVel = thisLevelPtr->newVel();

                  gradPsi.copyTo(xDerivativeComp,levelVel,yVelComp);
                  gradPsi.copyTo(yDerivativeComp,levelVel,xVelComp);

                  // multiply x-component of vel by -1
                  for (levelDit.reset(); levelDit.ok(); ++levelDit)
                    {
                      levelVel[levelDit].mult(-1.0,0,1);
                    }

                  // average down if necessary
                  if (lev < m_level)
                    {
                      LevelData<FArrayBox>& fineVel = thisLevelPtr->finerNSPtr()->newVel();

                      const DisjointBoxLayout& fineGrids = fineVel.getBoxes();
                      CoarseAverage VelAvgDown(fineGrids, fineVel.nComp(),
                                               refRatios[lev]);
                      VelAvgDown.averageToCoarse(levelVel,
                                                 fineVel);
                    }

                  thisLevelPtr = thisLevelPtr->crseNSPtr();
                } // end loop over levels
            } // end if SpaceDim == 2
          else if (SpaceDim == 3)
            {
              // need to loop over directions in this case and do each one
              // independently.

              // first define temporary storage
              Vector<LevelData<FArrayBox>*> compVort(m_level+1,NULL);
              Vector<LevelData<FArrayBox>*> compPsi(m_level+1,NULL);

              // start at finest (this) level and progress to coarse levels
              for (int lev=0; lev<compVort.size(); lev++)
                {
                  compVort[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
                                                           IntVect::Zero);

                  compPsi[lev] = new LevelData<FArrayBox>(amrGrids[lev],1,
                                                          IntVect::Unit);
                }

              // now compute stream function one direction at a time
              for (int dim=0; dim<SpaceDim; ++dim)
                {
                  // copy vorticity into temp storage
                  Interval dirComps(dim,dim);
                  Interval destComps(0,0);

                  for (int lev=0; lev < initialVort.size(); lev++)
                    {
                      initialVort[lev]->copyTo(dirComps, *compVort[lev],
                                               destComps);
                    }

                  streamSolver.solve(compPsi, compVort, numLevels-1, 0,
                                     true); // initialize compPsi (all levels) to zero

                  // now copy to composite psi
                  for (int lev=0; lev<compPsi.size(); lev++)
                    {
                      compPsi[lev]->copyTo(destComps, *streamFunc[lev],
                                           dirComps);
                    }

                } // end loop over dimensions

              // now compute vel = curl(streamFunc)
              thisLevelPtr = this;

              Interval psiComps(0,SpaceDim-1);
              BCFunc psiBCFunc = m_physBCPtr->streamBC();

              for (int lev=m_level; lev >= 0; lev--)
                {
                  const DisjointBoxLayout& levelGrids = streamFunc[lev]->getBoxes();


                  LevelData<FArrayBox>& levelVel = thisLevelPtr->newVel();

                  // do averaging down if necessary
                  if (lev < m_level)
                    {
                      LevelData<FArrayBox>& finePsi = *streamFunc[lev+1];
                      const DisjointBoxLayout& fineGrids = finePsi.getBoxes();
                      CoarseAverage PsiAvgDown(fineGrids, finePsi.nComp(),
                                               refRatios[lev]);
                      PsiAvgDown.averageToCoarse(*streamFunc[lev],
                                                 finePsi);
                    }

                  // apply physical BC's to stream function - inhomogeneous BC
                  bool isHomogeneous = false;
                  DataIterator levelDit = streamFunc[lev]->dataIterator();
                  for (levelDit.reset(); levelDit.ok(); ++levelDit)
                    {
                      const Box& streamFuncBox = (amrGrids[lev])[levelDit];
                      FArrayBox& streamFuncFAB = (*streamFunc[lev])[levelDit];

                      psiBCFunc(streamFuncFAB,
                                streamFuncBox,
                                amrDomains[lev],
                                amrDx[lev],
                                isHomogeneous);
                    }

                  // do coarse-fine BC's
                  if (lev > 0)
                    {
                      QuadCFInterp interpolator(amrGrids[lev],
                                                &amrGrids[lev-1],
                                                amrDx[lev],
                                                refRatios[lev-1],
                                                SpaceDim, amrDomains[lev]);

                      interpolator.coarseFineInterp(*streamFunc[lev],
                                                    *streamFunc[lev-1]);

                    }

                  streamFunc[lev]->exchange(psiComps);

                  for (levelDit.reset(); levelDit.ok(); ++levelDit)
                    {
                      for (int dir=0; dir<SpaceDim; dir++)
                        {
                          FORT_COMPUTEVORT(CHF_FRA1(levelVel[levelDit],dir),
                                           CHF_CONST_FRA((*streamFunc[lev])[levelDit]),
                                           CHF_BOX(levelGrids[levelDit]),
                                           CHF_CONST_REAL(amrDx[lev]),
                                           CHF_CONST_INT(dir));
                        } // end loop over directions

                      // get signs right!
                      levelVel[levelDit].mult(-1.0);
                      // add in background velocity
                      levelVel[levelDit].plus(s_backgroundVel, SpaceDim-1,1);
                    } // end loop over grids

                  // average down if necessary
                  if (lev < m_level)
                    {
                      LevelData<FArrayBox>& fineVel = thisLevelPtr->finerNSPtr()->newVel();

                      const DisjointBoxLayout& fineGrids = fineVel.getBoxes();
                      CoarseAverage VelAvgDown(fineGrids, fineVel.nComp(),
                                               refRatios[lev]);
                      VelAvgDown.averageToCoarse(levelVel,
                                                 fineVel);
                    }

                  thisLevelPtr = thisLevelPtr->crseNSPtr();
                } // end loop over levels

              // clean up temporary storage!
              for (int lev=0; lev<compVort.size(); lev++)
                {
                  if (compVort[lev] != NULL)
                    {
                      delete compVort[lev];
                      compVort[lev] = NULL;
                    }

                  if (compPsi[lev] != NULL)
                    {
                      delete compPsi[lev];
                      compPsi[lev] = NULL;
                    }
                }
            } // end if 3d
          else
            {
              MayDay::Error("HEY! What's with the oddly dimensional space???");
            }

          // clean up memory
          for (int lev=0; lev<amrGrids.size(); lev++)
            {
              if (initialVort[lev] != NULL)
                {
                  delete initialVort[lev];
                  initialVort[lev] = NULL;
                }
              if (streamFunc[lev] != NULL)
                {
                  delete streamFunc[lev];
                  streamFunc[lev] = NULL;
                }
            }
        } // end if finest level
    } // end if initialize velocity from vorticity field
  else
    {
      for (dit.reset(); dit.ok(); ++dit)
        {
          FORT_INITVEL(CHF_FRA(newVel()[dit]),
                       CHF_REAL(m_dx));
        }
    }  // end setting velocity field

  for (dit.reset(); dit.ok(); ++dit)
    {
      // now initialize Lambda
      newLambda()[dit].setVal(1.0, 0);

      // now init (cell-centered) scalars
      int isEdgeCentered = 0;
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          FORT_INITSCAL(CHF_FRA(newScal(comp)[dit]),
                        CHF_FRA(newVel()[dit]),
                        CHF_REAL(m_dx),
                        CHF_REAL(m_time),
                        CHF_INT(comp),
                        CHF_BOX(m_problem_domain.domainBox()),
                        CHF_INT(isEdgeCentered));

          // may also want old scalars initialized
          FORT_INITSCAL(CHF_FRA(oldScal(comp)[dit]),
                        CHF_FRA(newVel()[dit]),
                        CHF_REAL(m_dx),
                        CHF_REAL(m_time),
                        CHF_INT(comp),
                        CHF_BOX(m_problem_domain.domainBox()),
                        CHF_INT(isEdgeCentered));
        }
    }
}

// ---------------------------------------------------------------
void
AMRNavierStokes::postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::postInitialize " << m_level << endl;
    }

  // initialize data structures which haven't yet been initialized
  const DisjointBoxLayout& grids = newVel().getBoxes();
  levelSetup(grids);

  // set up for initial velocity projection -- do this from level 0,
  // since post_initialize is called from fine->coarse
  if (m_level == 0)
    {
      AMRNavierStokes* thisLevelData = this;
      while (!thisLevelData->finestLevel())
        {
          thisLevelData = thisLevelData->finerNSPtr();
        }
      CH_assert (thisLevelData->finestLevel());
      int numLevels = thisLevelData->m_level + 1;

      Vector<LevelData<FArrayBox>* > amrVel(numLevels);

      thisLevelData = this;
      for (int lev=m_level; lev <numLevels; lev++)
        {
          amrVel[lev] = &(thisLevelData->newVel());

          if (!thisLevelData->finestLevel())
            thisLevelData =
              dynamic_cast <AMRNavierStokes*> (thisLevelData->m_finer_level_ptr);
        }

      CCProjector& level0Proj = m_projection;

      // set physical boundary conditions on velocity
      bool isViscous = (s_nu > 0.0);
      thisLevelData = this;

      VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));

      for (int lev=0; lev< numLevels; lev++)
        {
          const ProblemDomain& levelDomain = thisLevelData->problemDomain();
          Real levelDx = thisLevelData->Dx();
          LevelData<FArrayBox>& thisAmrVel = *amrVel[lev];
          const DisjointBoxLayout& thisLevelGrids = thisAmrVel.getBoxes();
          velBC.applyBCs(thisAmrVel, thisLevelGrids,
                         levelDomain, levelDx,
                         false); // inhomogeneous
          if (thisLevelData->m_finer_level_ptr != NULL)
            {
              thisLevelData = dynamic_cast<AMRNavierStokes*> (thisLevelData->m_finer_level_ptr);
            }
        }

      bool homoBC = false;
      if (s_project_initial_vel)
        {
          level0Proj.initialVelocityProject(amrVel, homoBC);
        }

      thisLevelData = this;

      // need to reset boundary conditions here
      for (int lev=0; lev< numLevels; lev++)
        {
          const ProblemDomain& levelDomain = thisLevelData->problemDomain();
          Real levelDx = thisLevelData->Dx();
          LevelData<FArrayBox>& thisAmrVel = *amrVel[lev];
          const DisjointBoxLayout& thisLevelGrids = thisAmrVel.getBoxes();
          velBC.applyBCs(thisAmrVel, thisLevelGrids,
                         levelDomain, levelDx,
                         false); // inhomogeneous
          if (thisLevelData->m_finer_level_ptr != NULL)
            {
              thisLevelData = dynamic_cast<AMRNavierStokes*>
                (thisLevelData->m_finer_level_ptr);
            }
        }

      if (s_initialize_pressures)
        {
          initializeGlobalPressure();
        }

      // finally, if desired, dump out grids and processor mappings
      if (s_write_grids)
        {
#ifdef CH_MPI
          if (procID() == 0)
          {
#endif
            ofstream os("grids.out", ios::out);
            if (os.fail())
              {
                pout() << "cannot open grid output file " << os << endl;
                MayDay::Error();
              }

            thisLevelData = this;
            os << "NumLevels = " << numLevels << endl;
            for (int lev=0; lev<numLevels; lev++)
              {
                const DisjointBoxLayout& thisLevelGrids
                  = thisLevelData->m_vel_new_ptr->getBoxes();
                os << "Number of Grids on Level " << lev << ": "
                   << thisLevelGrids.size() << endl;
                // note that we don't need an endl here because
                // BoxLayout::operator<< provides a '\n' at the end
                os << thisLevelGrids;
                thisLevelData = thisLevelData->finerNSPtr();
              }
            os.close();
#ifdef CH_MPI
          } // end if procID == 0
#endif
        } // end if writing grids to file
    } // end if level 0
}
