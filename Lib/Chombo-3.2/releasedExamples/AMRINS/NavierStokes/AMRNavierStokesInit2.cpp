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
#include "VelBCHolder.H"
#include "SetValLevel.H"

// -------------------------------------------------------------
// this function manages the pressure initialization after
// initialization and regridding
void
AMRNavierStokes::initializeGlobalPressure()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::initializeGlobalPressure "
             << m_level << endl;
    }

  // index through levels to find out what finest level is
  AMRNavierStokes* thisNSPtr = this;
  while (!thisNSPtr->finestLevel())
    {
      thisNSPtr = thisNSPtr->finerNSPtr();
    }
  CH_assert(thisNSPtr->finestLevel());
  int finest_level = thisNSPtr->m_level;

  // first compute dtInit
  Real cur_time = m_time;
  // also save dt for each level so it can be reset later
  Vector<Real> dtSave(finest_level+1, 1.0e8);
  Real dtLevel;
  Real dtInit = 10000000.0;
  thisNSPtr = this;
  for (int lev=m_level; lev<=finest_level; lev++)
    {
      dtSave[lev] = thisNSPtr->m_dt;
      dtLevel = thisNSPtr->computeDt();
      if (dtLevel < dtInit) dtInit = dtLevel;
      thisNSPtr = thisNSPtr->finerNSPtr();
    }
  dtInit *= 0.5;

  // pressures should already be initialized to first
  // guess (most likely 0)

  // now loop through levels and compute estimate of pressure
  for (int iter=0; iter<s_num_init_passes; iter++)
    {
      int lbase = m_level;

      AMRNavierStokes* thisNSPtr = this;

      for (int lev=lbase; lev<=finest_level; lev++)
        {
          thisNSPtr->initializeLevelPressure(cur_time, dtInit);
          thisNSPtr = thisNSPtr->finerNSPtr();
        }

      // now reset times and dt
      thisNSPtr = this;
      for (int lev=lbase; lev<=finest_level; lev++)
        {
          thisNSPtr->resetStates(cur_time);
          thisNSPtr->dt(dtSave[lev]);
          thisNSPtr = thisNSPtr->finerNSPtr();
        }
    } // end loop over init passes
}

// -------------------------------------------------------------
// this function does all the level-based operations for
// initializing the pressure
void
AMRNavierStokes::initializeLevelPressure(Real a_currentTime, Real a_dtInit)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::initializeLevelPressure "
             << m_level << endl;
    }

  Real new_time = a_currentTime + a_dtInit;
  m_dt = a_dtInit;
  swapOldAndNewStates();

  const DisjointBoxLayout levelGrids = newVel().getBoxes();

  // just to be sure, initialize flux registers
  if (!finestLevel())
    {
      m_flux_register.setToZero();
      m_lambda_flux_reg.setToZero();
      for (int comp=0; comp<m_scal_fluxreg_ptrs.size(); comp++)
        {
          m_scal_fluxreg_ptrs[comp]->setToZero();
        }
    }

  // compute advection velocities
  LevelData<FluxBox> advectionVel(levelGrids, 1, IntVect::Unit);

  if (s_set_bogus_values)
    {
      setValLevel(advectionVel, s_bogus_value);
    }

  computeAdvectionVelocities(advectionVel);

  // since this is just for initialization, no need to update
  // advected scalars.

  // now predict velocities -- this returns the advection
  // term.  put this in new_vel.
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
 {
   VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));
   velBC.applyBCs(uStar, levelGrids,
                  problemDomain(), m_dx,
                  false); // inhomogeneous
 }

  // (Level projection calls exchange)
  //uStar.exchange(velComps);

  m_projection.LevelProject(uStar, crseVelPtr, new_time, m_dt);

  // as things stand now, physical BC's are reset in Level Projection

  // clean up temp storage
  if (crseVelPtr != NULL)
    {
      delete crseVelPtr;
      crseVelPtr = NULL;
    }
}

// -------------------------------------------------------------
void
AMRNavierStokes::smoothVelocityField(int a_lbase)
{
  // in this function, we take the filpatched fields stored in
  // the old-time storage (which have been modified to contain
  // s - mu*lap(s) in the regrid() function and perform an
  // elliptic solve which should result in a smoothed s

  // first need to loop over levels and put together amr storage
  // this should be called on the lbase level
  CH_assert (m_level == a_lbase);
  CH_assert (m_regrid_smoothing_done);

  AMRNavierStokes* thisNSPtr = this;
  while (!thisNSPtr->finestLevel())
    {
      thisNSPtr = thisNSPtr->finerNSPtr();
    }
  CH_assert (thisNSPtr->finestLevel());
  int finest_level = thisNSPtr->m_level;

  thisNSPtr = this;
  int startLev = m_level;
  if (m_level > 0)
    {
      startLev = m_level-1;
      thisNSPtr = thisNSPtr->crseNSPtr();
    }

  // now set up multilevel stuff
  Vector<LevelData<FArrayBox>*> oldS(finest_level+1, NULL);
  Vector<LevelData<FArrayBox>*> newS(finest_level+1, NULL);
  Vector<DisjointBoxLayout> amrGrids(finest_level+1);
  Vector<int> amrRefRatios(finest_level+1,0);
  Vector<Real> amrDx(finest_level+1,0);
  Vector<ProblemDomain> amrDomains(finest_level+1);
  // also will need to avg down new stuff
  Vector<CoarseAverage*> amrAvgDown(finest_level+1,NULL);

  // loop over levels, allocate temp storage for velocities,
  // set up for amrsolves
  for (int lev=startLev; lev<=finest_level; lev++)
    {
      const DisjointBoxLayout& levelGrids = thisNSPtr->newVel().getBoxes();
      // since AMRSolver can only handle one component at a
      // time, need to allocate temp space to copy stuff
      // into, and then back out of to compute Laplacian
      IntVect ghostVect(D_DECL(1,1,1));
      newS[lev] = new LevelData<FArrayBox>(levelGrids, 1, ghostVect);
      oldS[lev] = new LevelData<FArrayBox>(levelGrids,1,ghostVect);

      amrGrids[lev] = levelGrids;
      amrRefRatios[lev] = thisNSPtr->refRatio();
      amrDx[lev] = thisNSPtr->Dx();
      amrDomains[lev] = thisNSPtr->problemDomain();
      thisNSPtr = thisNSPtr->finerNSPtr();
      if (lev>startLev)
        {
          amrAvgDown[lev] = new CoarseAverage(levelGrids,1,
                                              amrRefRatios[lev-1]);
        }
    }

  AMRPoissonOpFactory localPoissonOpFactory;

  defineRegridAMROp(localPoissonOpFactory, amrGrids,
                    amrDomains, amrDx,
                    amrRefRatios, m_level);

  RelaxSolver<LevelData<FArrayBox> > bottomSolver;
  bottomSolver.m_verbosity = s_verbosity;

  AMRMultiGrid<LevelData<FArrayBox> > streamSolver;
  streamSolver.define(amrDomains[0],
                      localPoissonOpFactory,
                      &bottomSolver,
                      finest_level+1);
  streamSolver.m_verbosity = s_verbosity;
  streamSolver.m_eps = 1e-10;

  // now loop over velocity components
  if (s_nu > 0)
    {
      for (int dir=0; dir<SpaceDim; ++dir)
        {
          Interval velComps(dir,dir);
          Interval tempComps(0,0);
          thisNSPtr = this;
          if (startLev<m_level) thisNSPtr = thisNSPtr->crseNSPtr();
          for (int lev=startLev; lev<=finest_level; lev++)
            {
              thisNSPtr->oldVel().copyTo(velComps, *oldS[lev], tempComps);
              // do this as initial guess?
              thisNSPtr->oldVel().copyTo(velComps, *newS[lev], tempComps);
              thisNSPtr = thisNSPtr->finerNSPtr();
            }

          // now do elliptic solve
          // last two args are max level and base level
          streamSolver.solve(newS, oldS, finest_level, m_level,
                             true); // initialize newS to zero (converted AMRSolver)

          // average new s down to invalid regions
          for (int lev=finest_level; lev>startLev; lev--)
            {
              amrAvgDown[lev]->averageToCoarse(*newS[lev-1], *newS[lev]);
            }

          // now copy back to new velocity field
          thisNSPtr = this;
          for (int lev=m_level; lev <= finest_level; lev++)
            {
              newS[lev]->copyTo(tempComps, thisNSPtr->newVel(), velComps);

              thisNSPtr = thisNSPtr->finerNSPtr();
            }
        } // end loop over velocity components

    }// end if viscous

  // since scalars won't need temp storge, clean it up here
  for (int lev=startLev; lev<= finest_level; lev++)
    {
      if (newS[lev] != NULL)
        {
          delete newS[lev];
          newS[lev] = NULL;
        }
      if (oldS[lev] != NULL)
        {
          delete oldS[lev];
          oldS[lev] = NULL;
        }
    }

  // now do scalars
  for (int scalComp=0; scalComp < s_num_scal_comps; scalComp++)
    {
      // only do this if scalar is diffused
      if (s_scal_coeffs[scalComp] > 0)
        {
          thisNSPtr =this;
          for (int lev=startLev; lev <= finest_level; lev++)
            {
              newS[lev] = thisNSPtr->m_scal_new[scalComp];
              oldS[lev] = thisNSPtr->m_scal_old[scalComp];
              thisNSPtr = thisNSPtr->finerNSPtr();
            }

          // last two args are max level and base level
          streamSolver.solve(newS, oldS, finest_level, m_level,
                             true); // initialize newS to zero (converted AMRSolver)

          // now do averaging down
          for (int lev=finest_level; lev>m_level; lev--)
            {
              amrAvgDown[lev]->averageToCoarse(*newS[lev-1], *newS[lev]);
            }
        }  // end if this scalar is diffused
    }  // end loop over scalars

  // finally, loop over levels, clean up storage, and reset boolean
  thisNSPtr = this;
  for (int lev=m_level; lev<= finest_level; lev++)
    {
      if (amrAvgDown[lev] != NULL)
        {
          delete amrAvgDown[lev];
          amrAvgDown[lev] = NULL;
        }
      thisNSPtr->m_regrid_smoothing_done = false;
      thisNSPtr = thisNSPtr->finerNSPtr();
    }
}


// -------------------------------------------------------------
void
AMRNavierStokes::defineRegridAMROp(AMRPoissonOpFactory& a_factory,
                                   const Vector<DisjointBoxLayout>& a_grids,
                                   const Vector<ProblemDomain>& a_domains,
                                   const Vector<Real>& a_amrDx,
                                   const Vector<int>& a_refRatios,
                                   const int& a_lBase)
{
  // want to use dt from lbase
  Real dtLBase = m_dt;
  if (a_lBase > m_level)
    {
      AMRNavierStokes* thisNSPtr = finerNSPtr();
      while (a_lBase > thisNSPtr->m_level)
        {
          thisNSPtr = thisNSPtr->finerNSPtr();
        }
      dtLBase = thisNSPtr->m_dt;
    }
  else if (a_lBase < m_level)
    {
      AMRNavierStokes* thisNSPtr = crseNSPtr();
      while (a_lBase < thisNSPtr->m_level)
        {
          thisNSPtr = thisNSPtr->crseNSPtr();
        }
      dtLBase = thisNSPtr->m_dt;
    }

  // define coefficient
  Real mu = -s_regrid_smoothing_coeff*dtLBase*s_nu;

  // Would like to use extrap BC's, since they're probably the safest

  a_factory.define(a_domains[0],
                   a_grids,
                   a_refRatios,
                   a_amrDx[0],
                   m_physBCPtr->streamBC(), // BCHolder, reads from input file
                   1., // alpha
                   mu); // beta
}


// -------------------------------------------------------------
void AMRNavierStokes::defineViscousMGSolver(const DisjointBoxLayout& a_grids,
                                            const DisjointBoxLayout* a_crseGridsPtr,
                                            int                      a_refCrse)
{
  Tuple< RefCountedPtr<AMRLevelOpFactory< LevelData<FArrayBox> > >, SpaceDim>
    velTGAOpFactoryPtrs;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      // separate velTGAOpFactory just to hold s_nu to be sent to LevelTGA
      velTGAOpFactoryPtrs[idir] =
        RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
        ((AMRLevelOpFactory<LevelData<FArrayBox> >*)
         (new AMRPoissonOpFactory()));
    }

  Vector< RefCountedPtr<AMRLevelOpFactory< LevelData<FArrayBox> > > >
    scalTGAOpFactoryPtrs;
  if (s_num_scal_comps > 0)
    {
      scalTGAOpFactoryPtrs.resize(s_num_scal_comps);
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          scalTGAOpFactoryPtrs[comp] =
            RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
            ((AMRLevelOpFactory<LevelData<FArrayBox> >*)
             (new AMRPoissonOpFactory()));
        }
    }

  ProblemDomain baseDomain(m_problem_domain); // on this level
  // operator = alpha*I + beta*laplacian
  Real alpha = 1.; // default is 0.
  // Real beta = 1.; // default is 1.
  int numSolverLevels = (a_crseGridsPtr == NULL) ? 1 : 2;
  Vector<DisjointBoxLayout> allGrids(numSolverLevels);
  // not sure what becomes of this when numSolverLevels == 1
  Vector<int> refRatios(1, a_refCrse);
  if (a_crseGridsPtr != NULL)
    { // coarser level exists:  define solver on two levels
      baseDomain.coarsen(a_refCrse);
      allGrids[0] = *a_crseGridsPtr;
      allGrids[1] = a_grids;
      Real dxCrse = a_refCrse * m_dx;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          AMRPoissonOpFactory& velTGAOpFactoryDir =
            (AMRPoissonOpFactory&)(*velTGAOpFactoryPtrs[idir]);
          velTGAOpFactoryDir.define(baseDomain,
                                    allGrids,
                                    refRatios,
                                    dxCrse,
                                    m_physBCPtr->viscousSolveFuncBC(idir),
                                    alpha, s_nu);
        }
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          AMRPoissonOpFactory& scalTGAOpFactoryComp =
            (AMRPoissonOpFactory&)(*scalTGAOpFactoryPtrs[comp]);
          scalTGAOpFactoryComp.define(baseDomain,
                                      allGrids,
                                      refRatios,
                                      dxCrse,
                                      m_physBCPtr->extrapFuncBC(2), // order 2
                                      alpha, s_scal_coeffs[comp]);
        }
    }
  else
    { // no coarser level:  define solver on only one level
      allGrids[0] = a_grids;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          AMRPoissonOpFactory& velTGAOpFactoryDir =
            (AMRPoissonOpFactory&)(*velTGAOpFactoryPtrs[idir]);
          velTGAOpFactoryDir.define(baseDomain,
                                    a_grids,
                                    m_dx,
                                    m_physBCPtr->viscousSolveFuncBC(idir),
                                    -1, // max depth, unused
                                    alpha, s_nu);
        }
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          AMRPoissonOpFactory& scalTGAOpFactoryComp =
            (AMRPoissonOpFactory&)(*scalTGAOpFactoryPtrs[comp]);
          scalTGAOpFactoryComp.define(baseDomain,
                                      a_grids,
                                      m_dx,
                                      m_physBCPtr->extrapFuncBC(2), // order 2
                                      -1, // max depth, unused
                                      alpha, s_scal_coeffs[comp]);
        }
    }

  // You really need to delete this when you're done with the solvers.
  // But m_bottomSolver is a protected field of AMRMultiGrid.
  RelaxSolver<LevelData<FArrayBox> >* bottomSolverPtr = new
    RelaxSolver<LevelData<FArrayBox> >;
  bottomSolverPtr->m_verbosity = s_verbosity;

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_velMGsolverPtrs[idir] =
        RefCountedPtr< AMRMultiGrid<LevelData<FArrayBox> > >
        (new AMRMultiGrid<LevelData<FArrayBox> >() );
      m_velMGsolverPtrs[idir]->define(baseDomain, // on either this level or coarser level
                                      *velTGAOpFactoryPtrs[idir],
                                      bottomSolverPtr,
                                      numSolverLevels);
      m_velMGsolverPtrs[idir]->m_verbosity = s_verbosity;
      m_velMGsolverPtrs[idir]->m_eps = s_viscous_solver_tol;
      m_velMGsolverPtrs[idir]->m_pre = s_viscous_num_smooth_down;
      m_velMGsolverPtrs[idir]->m_post = s_viscous_num_smooth_up;

      m_TGAsolverPtrs[idir] = RefCountedPtr<LevelTGA>
        (new LevelTGA(allGrids, refRatios, baseDomain,
                      velTGAOpFactoryPtrs[idir], m_velMGsolverPtrs[idir]));
    }

  if (s_num_scal_comps > 0)
    {
      m_scalMGsolverPtrs.resize(s_num_scal_comps);
      m_TGAsolverScalPtrs.resize(s_num_scal_comps);
      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          m_scalMGsolverPtrs[comp] =
            RefCountedPtr< AMRMultiGrid<LevelData<FArrayBox> > >
            (new AMRMultiGrid<LevelData<FArrayBox> >() );
          m_scalMGsolverPtrs[comp]->define(baseDomain, // on either this level or coarser level
                                           *scalTGAOpFactoryPtrs[comp],
                                           bottomSolverPtr,
                                           numSolverLevels);
          m_scalMGsolverPtrs[comp]->m_verbosity = s_verbosity;
          m_scalMGsolverPtrs[comp]->m_eps = s_viscous_solver_tol;
          m_scalMGsolverPtrs[comp]->m_pre = s_viscous_num_smooth_down;
          m_scalMGsolverPtrs[comp]->m_post = s_viscous_num_smooth_up;

          m_TGAsolverScalPtrs[comp] = RefCountedPtr<LevelTGA>
            (new LevelTGA(allGrids, refRatios, baseDomain,
                          scalTGAOpFactoryPtrs[comp], m_scalMGsolverPtrs[comp]));
        }
    }
}
