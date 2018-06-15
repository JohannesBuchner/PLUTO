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

// initialize static variables here
const int AMRNavierStokes::s_num_vel_comps;
int AMRNavierStokes::s_num_scal_comps = 0;
Vector<Real> AMRNavierStokes::s_scal_coeffs;
Vector<Real> AMRNavierStokes::s_domLength(SpaceDim);
bool AMRNavierStokes::s_ppInit = false;
Real AMRNavierStokes::s_init_shrink = 1.0;
Real AMRNavierStokes::s_max_dt = 1.0e8;
Real AMRNavierStokes::s_prescribedDt = -1.0;
bool AMRNavierStokes::s_project_initial_vel = true;
bool AMRNavierStokes::s_initialize_pressures = true;
bool AMRNavierStokes::s_smooth_after_regrid = false;
Real AMRNavierStokes::s_regrid_smoothing_coeff = 4.0;
int  AMRNavierStokes::s_num_init_passes = 1;
bool AMRNavierStokes::s_reflux_momentum = true;
bool AMRNavierStokes::s_reflux_normal_momentum = true;
bool AMRNavierStokes::s_set_bogus_values = true;
bool AMRNavierStokes::s_tag_vorticity = false;
bool AMRNavierStokes::s_tag_theta = false;
Real AMRNavierStokes::s_vort_factor = 0.30;
int  AMRNavierStokes::s_tags_grow = 1;
Real AMRNavierStokes::s_bogus_value = 1.0e8;
Real AMRNavierStokes::s_nu = 0.0;
bool AMRNavierStokes::s_reflux_scal = true;
bool AMRNavierStokes::s_implicit_reflux = true;
bool AMRNavierStokes::s_implicit_scal_reflux = true;
bool AMRNavierStokes::s_applyFreestreamCorrection = true;
int  AMRNavierStokes::s_viscous_solver_type = 2;
Real AMRNavierStokes::s_viscous_solver_tol = 1e-7;
int AMRNavierStokes::s_viscous_num_smooth_up = 4;
int AMRNavierStokes::s_viscous_num_smooth_down = 4;
bool AMRNavierStokes::s_specifyInitialGrids = false;
string AMRNavierStokes::s_initialGridFile;
bool AMRNavierStokes::s_init_vel_from_vorticity = false;
Real AMRNavierStokes::s_backgroundVel = 0.0;

// initialize plotfile options
bool AMRNavierStokes::s_write_divergence = true;
bool AMRNavierStokes::s_write_lambda = true;
bool AMRNavierStokes::s_write_time_derivatives = false;
bool AMRNavierStokes::s_write_vorticity = true;
bool AMRNavierStokes::s_write_scalars = false;
bool AMRNavierStokes::s_write_dScalar_dt = false;
bool AMRNavierStokes::s_write_strains = false;
bool AMRNavierStokes::s_write_grad_eLambda = false;
bool AMRNavierStokes::s_write_error = false;
bool AMRNavierStokes::s_write_proc_ids = false;
bool AMRNavierStokes::s_compute_scal_err = false;
bool AMRNavierStokes::s_write_grids = false;

// ---------------------------------------------------------------
// names of state variables
const char* AMRNavierStokes::s_vel_names[s_num_vel_comps] =
{
  "x-velocity",
#if CH_SPACEDIM >= 2
  "y-velocity",
#endif
#if CH_SPACEDIM >= 3
  "z-velocity"
#endif
};

// ---------------------------------------------------------------
const char* AMRNavierStokes::s_scal_names[s_max_scal_comps] =
{
  "xMarker",
#if CH_SPACEDIM >= 2
  "yMarker",
#endif
#if CH_SPACEDIM >= 3
  "zMarker"
#endif
};

// ---------------------------------------------------------------
// default constructor
AMRNavierStokes::AMRNavierStokes() : m_cfl(0.5),
                                     m_level_steps(0),
                                     m_refine_threshold(0.5),
                                     m_finest_level(false),
                                     m_is_empty(true)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRNavierStokes default constructor" << endl;
  }
  m_vel_new_ptr = NULL;
  m_vel_old_ptr = NULL;
  m_lambda_new_ptr = NULL;
  m_lambda_old_ptr = NULL;

#ifdef DEBUG
  m_vel_save_ptr = NULL;
#endif

  m_regrid_smoothing_done = false;

  m_limitSolverCoarsening = false;
}

// ---------------------------------------------------------------
// full constructor
AMRNavierStokes::AMRNavierStokes(AMRLevel* a_coarser_level_ptr,
                                 const Box& a_prob_domain,
                                 int a_level,
                                 int a_ref_ratio)
  : m_cfl(0.5),
    m_limitSolverCoarsening(false),
    m_level_steps(0),
    m_refine_threshold(0.5),
    m_finest_level(false),
    m_is_empty(true)
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes full constructor" << endl;
    }

  m_vel_new_ptr = NULL;
  m_vel_old_ptr = NULL;
  m_lambda_new_ptr = NULL;
  m_lambda_old_ptr = NULL;

#ifdef DEBUG
  m_vel_save_ptr = NULL;
#endif

  ProblemDomain physdomain(a_prob_domain);
  define(a_coarser_level_ptr, physdomain, a_level, a_ref_ratio);
}

// ---------------------------------------------------------------
// full constructor
AMRNavierStokes::AMRNavierStokes(AMRLevel* a_coarser_level_ptr,
                                 const ProblemDomain& a_prob_domain,
                                 int a_level,
                                 int a_ref_ratio)
  : m_cfl(0.5),
    m_limitSolverCoarsening(false),
    m_level_steps(0),
    m_refine_threshold(0.5),
    m_finest_level(false),
    m_is_empty(true)
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes full constructor" << endl;
    }

  m_vel_new_ptr = NULL;
  m_vel_old_ptr = NULL;
  m_lambda_new_ptr = NULL;
  m_lambda_old_ptr = NULL;

#ifdef DEBUG
  m_vel_save_ptr = NULL;
#endif

  define(a_coarser_level_ptr, a_prob_domain, a_level, a_ref_ratio);
}

// ---------------------------------------------------------------
// destructor
AMRNavierStokes::~AMRNavierStokes()
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes destructor" << endl;
    }

  if (m_vel_new_ptr != NULL)
    {
      delete m_vel_new_ptr;
      m_vel_new_ptr = NULL;
    }

  if (m_vel_old_ptr != NULL)
    {
      delete m_vel_old_ptr;
      m_vel_old_ptr = NULL;
    }

  if (m_lambda_new_ptr != NULL)
    {
      delete m_lambda_new_ptr;
      m_lambda_new_ptr = NULL;
    }

  if (m_lambda_old_ptr != NULL)
    {
      delete m_lambda_old_ptr;
      m_lambda_old_ptr = NULL;
    }

  // loop over scalars and delete
  int nScalComp = m_scal_new.size();
  for (int comp=0; comp<nScalComp; comp++)
    {
      if (m_scal_new[comp] != NULL)
        {
          delete m_scal_new[comp];
          m_scal_new[comp] = NULL;
        }
    }

  nScalComp = m_scal_old.size();
  for (int comp=0; comp<nScalComp; comp++)
    {
      if (m_scal_old[comp] != NULL)
        {
          delete m_scal_old[comp];
          m_scal_old[comp] = NULL;
        }
    }

#ifdef DEBUG
  nScalComp = m_scal_save.size();
  for (int comp=0; comp<nScalComp; comp++)
    {
      if (m_scal_save[comp] != NULL)
        {
          delete m_scal_save[comp];
          m_scal_save[comp] = NULL;
        }
    }

  if (m_vel_save_ptr != NULL)
    {
      delete m_vel_save_ptr;
      m_vel_save_ptr = NULL;
    }

  if (m_lambda_save_ptr != NULL)
    {
      delete m_lambda_save_ptr;
      m_lambda_save_ptr = NULL;
    }
#endif

  if (m_physBCPtr != NULL)
    {
      delete m_physBCPtr;
    }

  nScalComp = m_scal_fluxreg_ptrs.size();
  for (int comp=0; comp<nScalComp; comp++)
    {
      if (m_scal_fluxreg_ptrs[comp] != NULL)
        {
          delete m_scal_fluxreg_ptrs[comp];
          m_scal_fluxreg_ptrs[comp] = NULL;
        }
    }

  nScalComp = m_patchGodScalars.size();
  for (int comp = 0; comp < nScalComp; comp++)
  {
    if (m_patchGodScalars[comp] != NULL)
    {
      delete m_patchGodScalars[comp];
      m_patchGodScalars[comp] = NULL;
    }
  }
}

// ---------------------------------------------------------------
// define for starting from scratch
void
AMRNavierStokes::define(AMRLevel* a_coarser_level_ptr,
                        const Box& a_problem_domain,
                        int a_level,
                        int a_ref_ratio)
{
  ProblemDomain physdomain(a_problem_domain);
  define(a_coarser_level_ptr, physdomain, a_level, a_ref_ratio);
}

// ---------------------------------------------------------------
// define for starting from scratch
void
AMRNavierStokes::define(AMRLevel* a_coarser_level_ptr,
                        const ProblemDomain& a_problem_domain,
                        int a_level,
                        int a_ref_ratio)
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::define " << a_level << endl;
    }

  AMRLevel::define(a_coarser_level_ptr, a_problem_domain,
                   a_level, a_ref_ratio);

  // read stuff from parmParse database
  if (!s_ppInit)
    {
      readParameters();
    }

  if (a_coarser_level_ptr != NULL)
    {
      AMRNavierStokes* amrns_ptr = dynamic_cast <AMRNavierStokes*> (a_coarser_level_ptr);
      if (amrns_ptr != NULL)
        {
          m_cfl = amrns_ptr->m_cfl;
          m_refine_threshold = amrns_ptr->m_refine_threshold;
        }
      else
        {
          MayDay::Error ("AMRNavierStokes::define: a_coarser_level_ptr is not castable to AMRNavierStokes*");
        }
      // set coarser levels finer_level_ptr to point to this
      a_coarser_level_ptr->finerLevelPtr(static_cast<AMRLevel*> (this));
    }

  m_regrid_smoothing_done = false;
}

// ---------------------------------------------------------------
// "virtual" full constructor
AMRLevel*
AMRNavierStokes::makeAMRLevel (AMRLevel* a_coarser_level_ptr,
                               const Box& a_prob_domain,
                               int a_level,
                               int a_ref_ratio) const
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::makeAMRLevel " << a_level << endl;
    }

  ProblemDomain physdomain(a_prob_domain);
  AMRNavierStokes* amrns_ptr = new AMRNavierStokes(a_coarser_level_ptr,
                                                   physdomain,
                                                   a_level, a_ref_ratio);
  return (static_cast <AMRLevel*> (amrns_ptr));
}

// ---------------------------------------------------------------
// "virtual" full constructor
AMRLevel*
AMRNavierStokes::makeAMRLevel (AMRLevel* a_coarser_level_ptr,
                               const ProblemDomain& a_prob_domain,
                               int a_level,
                               int a_ref_ratio) const
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::makeAMRLevel " << a_level << endl;
    }

  AMRNavierStokes* amrns_ptr = new AMRNavierStokes(a_coarser_level_ptr,
                                                   a_prob_domain,
                                                   a_level,
                                                   a_ref_ratio);
  return (static_cast <AMRLevel*> (amrns_ptr));
}
