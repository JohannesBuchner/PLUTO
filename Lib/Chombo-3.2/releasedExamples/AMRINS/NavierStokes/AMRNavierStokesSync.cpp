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
#include "computeNorm.H"
#include "computeSum.H"
#include "Divergence.H"
#include "probF_F.H"

// ---------------------------------------------------------------
// things do do after a basic timestep
void
AMRNavierStokes::postTimeStep()
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::postTimeStep " << m_level << endl;
    }

  // none of these things are done if this is the finest level

  if (!finestLevel())
    {
      AMRNavierStokes* fineAMRNavierStokesPtr = finerNSPtr();
      // create boundary condition object
      bool isViscous = (s_nu > 0.0);
      VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));

      // first do refluxing and avgDown for conservation
      // do momentum refluxing first
      Real refluxScale = -1.0/m_dx; // petermc made negative, 7 Dec 07
      LevelData<FArrayBox>& newVel = *m_vel_new_ptr;
      if (s_reflux_momentum)
        {
          if (s_implicit_reflux)
            {
              // defer this until we're doing sync projection
            }
          else
            {
              m_flux_register.reflux(newVel, refluxScale);
              CoarseAverage& velAvgDown = fineAMRNavierStokesPtr->m_coarse_average;
              velAvgDown.averageToCoarse(newVel,*fineAMRNavierStokesPtr->m_vel_new_ptr);
            }
        }

      // now do lamdba
      LevelData<FArrayBox>& newLambda = *m_lambda_new_ptr;
      m_lambda_flux_reg.reflux(newLambda,refluxScale);
      CoarseAverage& lambdaAvgDown
        = fineAMRNavierStokesPtr->m_coarse_average_lambda;
      lambdaAvgDown.averageToCoarse(newLambda,*fineAMRNavierStokesPtr->m_lambda_new_ptr);

      // finally, do scalars -- loop over components
      if (s_reflux_scal)
        {
          if (s_implicit_scal_reflux)
            {
              // if implicit refluxing, defer until we're doing multilevel stuff
            }
          else
            {
              for (int comp=0; comp<s_num_scal_comps; comp++)
                {
                  LevelData<FArrayBox>& newScalars = *m_scal_new[comp];
                  m_scal_fluxreg_ptrs[comp]->reflux(newScalars, refluxScale);
                  CoarseAverage& scalAvgDown =
                    fineAMRNavierStokesPtr->m_coarse_average_lambda;

                  scalAvgDown.averageToCoarse(newScalars,
                                              *fineAMRNavierStokesPtr->m_scal_new[comp]);
                }
            }
        } // end if reflux scalars

      // now call sync projection if necessary
      Real crseTime = 0.0;
      if (m_level > 0) crseTime = m_coarser_level_ptr->time();

      // do multilevel operations if this is the coarsest level or if
      // coarser level is not at the same time as this level
      if (m_level == 0 || (abs(crseTime-m_time)>TIME_EPS))
        {
          int finest_level = m_level;
          // index through levels to find out what finest level is
          AMRNavierStokes* thisNSPtr = this;
          while (!thisNSPtr->finestLevel())
            {
              thisNSPtr = thisNSPtr->finerNSPtr();
            }
          CH_assert(thisNSPtr->finestLevel());
          finest_level = thisNSPtr->m_level;

          // first, do implicit refluxing for scalars
          if (s_implicit_scal_reflux && s_reflux_scal)
            {
              Vector<LevelData<FArrayBox>* > scalRefluxCorr(finest_level+1,NULL);
              Vector<LevelData<FArrayBox>* > scalRefluxRHS(finest_level+1,NULL);
              Vector<DisjointBoxLayout> AmrGrids(finest_level+1);
              Vector<int> AmrRefRatios(finest_level+1);
              Vector<Real> AmrDx(finest_level+1);
              Vector<ProblemDomain> AmrDomains(finest_level+1);

              // loop over levels, allocate storage, set up for AMRMultiGrid
              // solve
              thisNSPtr = this;
              int startLev=m_level;
              // if crser level exists, define it as well for BC's
              if (startLev > 0)
                {
                  startLev = startLev-1;
                  thisNSPtr = thisNSPtr->crseNSPtr();
                }
              for (int lev=startLev; lev<=finest_level; lev++)
                {
                  const DisjointBoxLayout& levelGrids = thisNSPtr->newVel().getBoxes();
                  // recall that AMRMultiGrid can only do one component.
                  // rhs has no ghost cells
                  scalRefluxRHS[lev] = new LevelData<FArrayBox>(levelGrids, 1);
                  //soln has one layer of ghost cells
                  IntVect ghostVect(D_DECL(1,1,1));
                  scalRefluxCorr[lev] = new LevelData<FArrayBox>(levelGrids,
                                                                 1,ghostVect);
                  AmrGrids[lev] = levelGrids;
                  AmrRefRatios[lev] = thisNSPtr->refRatio();
                  AmrDx[lev] = thisNSPtr->Dx();
                  AmrDomains[lev] = thisNSPtr->problemDomain();

                  // initialize corr, RHS to 0
                  DataIterator levelDit = scalRefluxRHS[lev]->dataIterator();
                  LevelData<FArrayBox>& levelRefluxRHS = *(scalRefluxRHS[lev]);
                  LevelData<FArrayBox>& levelRefluxCorr = *(scalRefluxCorr[lev]);
                  for (levelDit.reset(); levelDit.ok(); ++levelDit)
                    {
                      levelRefluxRHS[levelDit()].setVal(0.0);
                      levelRefluxCorr[levelDit()].setVal(0.0);
                    }
                  thisNSPtr = thisNSPtr->finerNSPtr();
                }

              Interval solverComps(0,0);

              // now do each scalar component
              for (int scalComp=0; scalComp<s_num_scal_comps; scalComp++)
                {
                  // loop over levels and establish RHS, corr
                  thisNSPtr = this;

                  for (int lev=m_level; lev<=finest_level; lev++)
                    {
                      if (lev < finest_level)
                        {
                          LevelData<FArrayBox>& levelRefluxRHS = *(scalRefluxRHS[lev]);
                          // first set to 0
                          DataIterator levelDit = levelRefluxRHS.dataIterator();
                          for (levelDit.reset(); levelDit.ok(); ++levelDit)
                            {
                              levelRefluxRHS[levelDit()].setVal(0.0);
                            }
                          refluxScale = -1.0/AmrDx[lev]; // petermc made negative, 7 Dec 07
                          thisNSPtr->m_scal_fluxreg_ptrs[scalComp]->reflux(levelRefluxRHS,
                                                                           refluxScale);
                          // initial guess for correction is RHS
                          levelRefluxRHS.copyTo(solverComps,
                                                *(scalRefluxCorr[lev]),
                                                solverComps);
                        }

                      thisNSPtr = thisNSPtr->finerNSPtr();
                    } // end loop over levels to compute RHS

                  // define levelOp and solver for diffusive solves;
                  // need to define
                  // new AMRMultiGrid for each component because they
                  // may have different coefficients

                  Real nuComp = s_scal_coeffs[scalComp];
                  // only do all of this if  actually diffusive
                  if (nuComp > 0)
                    {
                      int numLevels = finest_level+1;

                      // This is a Helmholtz operator
                      Real alpha = 1.0;
                      Real beta = -nuComp*m_dt;

                      AMRPoissonOpFactory diffusiveOpFactory;
                      diffusiveOpFactory.define(AmrDomains[0],
                                                AmrGrids,
                                                AmrRefRatios,
                                                AmrDx[0],
                                                m_physBCPtr->scalarRefluxSolveBC(scalComp),
                                                alpha,
                                                beta);

                      RelaxSolver<LevelData<FArrayBox> > bottomSolver;
                      bottomSolver.m_verbosity = s_verbosity;

                      AMRMultiGrid<LevelData<FArrayBox> > diffusionSolver;
                      AMRLevelOpFactory<LevelData<FArrayBox> >& castFact = ( AMRLevelOpFactory<LevelData<FArrayBox> >&) diffusiveOpFactory;
                      diffusionSolver.define(AmrDomains[0],
                                             castFact,
                                             &bottomSolver,
                                             numLevels);
                      diffusionSolver.m_verbosity = s_verbosity;
                      diffusionSolver.m_eps = s_viscous_solver_tol;

                      // set convergence metric to be the norm of the scalar
                      Interval comps(scalComp, scalComp);
                      Vector<LevelData<FArrayBox>* > vectScal(finest_level+1, NULL);
                      thisNSPtr = this;
                      for (int lev=m_level; lev<=finest_level; lev++)
                        {
                          vectScal[lev] = &thisNSPtr->newScal(scalComp);
                          if (lev < finest_level)
                            {
                              thisNSPtr = thisNSPtr->finerNSPtr();
                            }
                        }
                      // This needs to be fixed when AMRMultiGrid has this
                      // capability.
                      //
                      int normType = 0;
                      Real scalNorm = computeNorm(vectScal, AmrRefRatios,
                                                  m_dx, comps, normType,
                                                  m_level);
                      // if scalNorm = 0, not sure what we're doing here,
                      // but let solver default to whatever it normally
                      // does in that case
                      if (scalNorm != 0.)
                        {
                          // diffusionSolver.setConvergenceMetric(scalNorm,0);
                          diffusionSolver.m_convergenceMetric = scalNorm;
                        }

                      // quick check -- want sum(RHS) to be 0
#ifdef DEBUG
                      Interval solverComps(0,0);
                      Real sumRHS;
                      sumRHS = computeSum(scalRefluxRHS, AmrRefRatios, m_dx,
                                          solverComps, m_level);
                      pout() << "Sum(RHS) for implicit reflux for "
                             << s_scal_names[scalComp] << " = "
                             << setiosflags(ios::scientific) << sumRHS << endl;
#endif

                      // now solve
                      diffusionSolver.solve(scalRefluxCorr, scalRefluxRHS,
                                            finest_level, m_level,
                                            false); // don't initialize to zero

                      // now increment scalars with reflux correction
                      // go from finest->coarsest so that we can also avgDown
                      // increment NS pointer to finest level
                      thisNSPtr = this;
                      while (!thisNSPtr->finestLevel())
                        {
                          thisNSPtr = thisNSPtr->finerNSPtr();
                        }
                      CH_assert(thisNSPtr->m_level == finest_level);

                      for (int lev=finest_level; lev>=m_level; lev--)
                        {
                          LevelData<FArrayBox>& levelScal = thisNSPtr->newScal(scalComp);
                          LevelData<FArrayBox>& levelCorr = *(scalRefluxCorr[lev]);
                          DataIterator levelDit = levelCorr.dataIterator();
                          for (levelDit.reset(); levelDit.ok(); ++levelDit)
                            {
                              levelScal[levelDit()] += levelCorr[levelDit()];
                            }

                          // if a finer level exists, do avgDown as well
                          if (lev < finest_level)
                            {
                              AMRNavierStokes& fineNS = *(thisNSPtr->finerNSPtr());
                              // quick sanity check
                              CH_assert (fineNS.m_level = lev+1);
                              CoarseAverage& scalAvgDown = fineNS.m_coarse_average_lambda;
                              LevelData<FArrayBox>& fineScal = fineNS.newScal(scalComp);
                              scalAvgDown.averageToCoarse(levelScal, fineScal);
                            }
                          thisNSPtr = thisNSPtr->crseNSPtr();
                        } // end loop over levels

                    } // end if diffusive
                  else  // not diffusive
                    {
                      // if this component is not diffusive, just add RHS
                      // to scalar now increment scalars with reflux
                      // correction
                      // go from finest->coarsest so that we can also avgDown
                      // increment NS pointer to finest level
                      thisNSPtr = this;
                      while (!thisNSPtr->finestLevel())
                        {
                          thisNSPtr = thisNSPtr->finerNSPtr();
                        }
                      CH_assert(thisNSPtr->m_level == finest_level);
                      for (int lev=finest_level; lev>=m_level; lev--)
                        {
                          LevelData<FArrayBox>& levelScal = thisNSPtr->newScal(scalComp);
                          LevelData<FArrayBox>& levelCorr = *(scalRefluxRHS[lev]);
                          DataIterator levelDit = levelCorr.dataIterator();
                          for (levelDit.reset(); levelDit.ok(); ++levelDit)
                            {
                              levelScal[levelDit()] += levelCorr[levelDit()];
                            }

                          // if a finer level exists, do avgDown as well
                          if (lev < finest_level)
                            {
                              AMRNavierStokes& fineNS = *(thisNSPtr->finerNSPtr());
                              // quick sanity check
                              CH_assert (fineNS.m_level == lev+1);
                              CoarseAverage& scalAvgDown = fineNS.m_coarse_average_lambda;
                              LevelData<FArrayBox>& fineScal = fineNS.newScal(scalComp);
                              scalAvgDown.averageToCoarse(levelScal, fineScal);
                            }
                          thisNSPtr = thisNSPtr->crseNSPtr();
                        } // end loop over levels
                    } // end if not diffusive
                } // end loop over scalar components

              // clean up temporary scalar storage
              for (int lev=0; lev<=finest_level; lev++)
                {
                  if (scalRefluxRHS[lev] != NULL)
                    {
                      delete scalRefluxRHS[lev];
                      scalRefluxRHS[lev] = NULL;
                    }

                  if (scalRefluxCorr[lev] != NULL)
                    {
                      delete scalRefluxCorr[lev];
                      scalRefluxCorr[lev] = NULL;
                    }
                }
            } // end implicit scalar refluxing

          // now allocate container for composite velocity and lambda
          Vector< LevelData<FArrayBox>* > compVel(finest_level+1,NULL);
          Vector< LevelData<FArrayBox>* > compLambda(finest_level+1,NULL);

          if (m_level > 0)
            {
              // coarser level will be coarse-fine BC
              AMRNavierStokes* crseAMRNavierStokesPtr = crseNSPtr();
              DisjointBoxLayout crseGrids = crseAMRNavierStokesPtr->newVel().getBoxes();

              LevelData<FArrayBox>* crseVelPtr =
                new LevelData<FArrayBox>(crseGrids, SpaceDim);

              crseAMRNavierStokesPtr->velocity(*crseVelPtr, m_time);
              compVel[m_level-1] = crseVelPtr;

              LevelData<FArrayBox>* crseLambdaPtr =
                new LevelData<FArrayBox>(crseGrids, 1);

              crseAMRNavierStokesPtr->lambda(*crseLambdaPtr, m_time);
              compLambda[m_level-1] = crseLambdaPtr;
            }

          thisNSPtr = this;
          // now loop over levels and construct composite vel, lambda
          for (int lev=m_level; lev <= finest_level; lev++)
            {
              compVel[lev] = thisNSPtr->m_vel_new_ptr;
              compLambda[lev] = thisNSPtr->m_lambda_new_ptr;
              const DisjointBoxLayout& levelGrids = compVel[lev]->getBoxes();
              velBC.applyBCs(*compVel[lev], levelGrids,
                             thisNSPtr->problemDomain(), thisNSPtr->m_dx,
                             false); // not homogeneous
              thisNSPtr = thisNSPtr->finerNSPtr();
            }

          // before we do sync projection, do implicit refluxing,
          // if appropriate
          if (s_reflux_momentum && s_implicit_reflux)
            {
              // loop over levels and compute RHS
              Vector<LevelData<FArrayBox>* > refluxRHS(finest_level+1,NULL);
              Vector<LevelData<FArrayBox>* > refluxCorr(finest_level+1,NULL);
              // this is necessary because while solver only can do
              // one component, levelfluxRegister::reflux can only
              // do all of them at once.
              Vector<LevelData<FArrayBox>* > tempRefluxData(finest_level+1,NULL);
              Vector<DisjointBoxLayout> AmrGrids(finest_level+1);
              Vector<int> AmrRefRatios(finest_level+1);
              Vector<Real> AmrDx(finest_level+1);
              Vector<ProblemDomain> AmrDomains(finest_level+1);

              // loop over levels, allocate storage, set up for AMRMultiGrid
              // solve
              thisNSPtr = this;
              int startLev=m_level;
              // if crser level exists, define it as well for BC's
              if (startLev > 0)
                {
                  startLev = startLev-1;
                  thisNSPtr = thisNSPtr->crseNSPtr();
                }
              for (int lev=startLev; lev<=finest_level; lev++)
                {
                  const DisjointBoxLayout& levelGrids = thisNSPtr->newVel().getBoxes();
                  // recall that AMRMultiGrid can only do one component
                  // rhs has no ghost cells
                  refluxRHS[lev] = new LevelData<FArrayBox>(levelGrids, 1);
                  tempRefluxData[lev] = new LevelData<FArrayBox>(levelGrids,
                                                                 SpaceDim);
                  //soln has one layer of ghost cells
                  IntVect ghostVect(D_DECL(1,1,1));
                  refluxCorr[lev] = new LevelData<FArrayBox>(levelGrids,1,
                                                             ghostVect);
                  AmrGrids[lev] = levelGrids;
                  AmrRefRatios[lev] = thisNSPtr->refRatio();
                  AmrDx[lev] = thisNSPtr->Dx();
                  AmrDomains[lev] = thisNSPtr->problemDomain();

                  // initialize rhs to 0
                  DataIterator levelDit = tempRefluxData[lev]->dataIterator();
                  LevelData<FArrayBox>& levelRefluxData = *(tempRefluxData[lev]);
                  for (levelDit.reset(); levelDit.ok(); ++levelDit)
                    {
                      levelRefluxData[levelDit()].setVal(0.0);
                    }

                  // while we're here, do refluxing.
                  // (recall that startLev may be coarser than m_level
                  // for BC's,
                  // however, don't want to do anything to that level)
                  if  ((lev >= m_level) && (lev < finest_level))
                    {
                      refluxScale = -1.0/AmrDx[lev]; // petermc made negative, 7 Dec 07
                      thisNSPtr->m_flux_register.reflux(levelRefluxData,
                                                        refluxScale);
                    }

                  thisNSPtr = thisNSPtr->finerNSPtr();
                } // end loop over levels

              // coarse BC is 0
              if (m_level > 0)
                {
                  LevelData<FArrayBox>& crseBCData = *refluxCorr[m_level-1];
                  DataIterator crseDit = crseBCData.dataIterator();
                  for (crseDit.reset(); crseDit.ok(); ++crseDit)
                    {
                      crseBCData[crseDit()].setVal(0.0);
                    }
                }

              // set convergence metric to be norm of velocity
              // for now just do component-wise maxnorm
              // compute norm over all directions and then use for
              // all velocity components (in order to be consistent)
              Vector<LevelData<FArrayBox>* > vectVel(finest_level+1, NULL);
              thisNSPtr = this;
              for (int lev=m_level; lev<=finest_level; lev++)
                {
                  vectVel[lev] = thisNSPtr->m_vel_new_ptr;
                  if (lev < finest_level)
                    {
                      thisNSPtr = thisNSPtr->finerNSPtr();
                    }
                }

              int normType = 0;
              Interval allVelComps(0, SpaceDim-1);

              Real velNorm = computeNorm(vectVel, AmrRefRatios,
                                         m_dx, allVelComps, normType, m_level);

              // now loop over directions
              Interval solverComps(0,0);
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  Interval velComps(dir,dir);
                  for (int lev=m_level; lev<=finest_level; ++lev)
                    {
                      // copy rhs to single-component holder
                      tempRefluxData[lev]->copyTo(velComps,*(refluxRHS[lev]),
                                                  solverComps);
                      // initial guess for correction is RHS.
                      LevelData<FArrayBox>& levelCorr = *(refluxCorr[lev]);
                      DataIterator levelDit = levelCorr.dataIterator();
                      for (levelDit.reset(); levelDit.ok(); ++levelDit)
                        {
                          levelCorr[levelDit()].setVal(0.0);
                        }
                      refluxRHS[lev]->copyTo(solverComps,*(refluxCorr[lev]),
                                             solverComps);
                    } // end set initial guess

                  // now set up solver
                  int numLevels = finest_level+1;

                  // This is a Helmholtz operator
                  Real alpha = 1.0;
                  Real beta = -s_nu*m_dt;

                  AMRPoissonOpFactory viscousOpFactory;
                  viscousOpFactory.define(AmrDomains[0],
                                          AmrGrids,
                                          AmrRefRatios,
                                          AmrDx[0],
                                          m_physBCPtr->viscousRefluxBC(dir),
                                          alpha,
                                          beta);

                  RelaxSolver<LevelData<FArrayBox> > bottomSolver;
                  bottomSolver.m_verbosity = s_verbosity;

                  AMRMultiGrid<LevelData<FArrayBox> > viscousSolver;
                  AMRLevelOpFactory<LevelData<FArrayBox> >& viscCastFact = (AMRLevelOpFactory<LevelData<FArrayBox> >&) viscousOpFactory;
                  viscousSolver.define(AmrDomains[0],
                                       viscCastFact,
                                       &bottomSolver,
                                       numLevels);

                  viscousSolver.m_verbosity = s_verbosity;

                  viscousSolver.m_eps = s_viscous_solver_tol;

                  viscousSolver.m_pre  = s_viscous_num_smooth_down;
                  viscousSolver.m_post = s_viscous_num_smooth_up;

                  // This needs to be fixed when AMRMultiGrid has this
                  // capability.
                  //
                  // viscousSolver.setConvergenceMetric(velNorm, 0);
                  viscousSolver.m_convergenceMetric = velNorm;

                  // now solve
                  viscousSolver.solve(refluxCorr, refluxRHS,
                                      finest_level, m_level,
                                      false); // don't initialize to zero

                  // now increment velocity with reflux correction
                  for (int lev=m_level; lev<=finest_level; ++lev)
                    {
                      LevelData<FArrayBox>& levelVel = *(compVel[lev]);
                      LevelData<FArrayBox>& levelCorr = *(refluxCorr[lev]);
                      DataIterator levelDit = levelCorr.dataIterator();
                      for (levelDit.reset(); levelDit.ok(); ++levelDit)
                        {
                          levelVel[levelDit()].plus(levelCorr[levelDit()],0,dir,1);
                        }
                    }
                } // end loop over directions

              // clean up storage
              for (int lev=startLev; lev<=finest_level; lev++)
                {
                  if (refluxRHS[lev] != NULL)
                    {
                      delete refluxRHS[lev];
                      refluxRHS[lev] = NULL;
                    }

                  if (refluxCorr[lev] != NULL)
                    {
                      delete refluxCorr[lev];
                      refluxCorr[lev] = NULL;
                    }

                  if (tempRefluxData[lev] != NULL)
                    {
                      delete tempRefluxData[lev];
                      tempRefluxData[lev] = NULL;
                    }
                }
            } // end implicit reflux

          // re-apply physical boundary conditions
          thisNSPtr = this;
          for (int lev = m_level; lev < finest_level; lev++)
            {
              LevelData<FArrayBox>& levelVel = thisNSPtr->newVel();
              const DisjointBoxLayout& levelGrids = levelVel.getBoxes();
              velBC.applyBCs(levelVel, levelGrids,
                             thisNSPtr->problemDomain(), thisNSPtr->m_dx,
                             false); // not homogeneous
              thisNSPtr = thisNSPtr->finerNSPtr();
            }

          // call projection-based sync operations
          m_projection.doSyncOperations(compVel, compLambda, m_time, m_dt);

          // re-apply physical boundary conditions
          thisNSPtr = this;
          for (int lev = m_level; lev < finest_level; lev++)
            {
              LevelData<FArrayBox>& levelVel = thisNSPtr->newVel();
              const DisjointBoxLayout& levelGrids = levelVel.getBoxes();
              velBC.applyBCs(levelVel, levelGrids,
                             thisNSPtr->problemDomain(), thisNSPtr->m_dx,
                             false); // not homogeneous
              thisNSPtr = thisNSPtr->finerNSPtr();
            }

          // this is a diagnostic -- compute max(divergence)
          // also compute max(vorticity)
          // and integrals of kinetic energy and vorticity.
          if (m_level == 0)
            {
              LevelData<FArrayBox>* crseVelPtr = NULL;
              LevelData<FArrayBox>* fineVelPtr = NULL;
              int nRefCrse = -1;
              int nRefFine = -1;

              AMRNavierStokes* thisLevelPtr = this;

              Real maxDiv = 0.0;
              Real maxVort = 0.0;
              Real maxLambdaErr = 0.0;
              Real integralLambda = 0.0;
              Real integralKE = 0.0;

              for (int lev=m_level; lev <= finest_level; lev++)
                {
                  LevelData<FArrayBox>& levelVel = *compVel[lev];
                  const DisjointBoxLayout& levelGrids = levelVel.getBoxes();
                  LevelData<FArrayBox> levelDiv(levelGrids, 1);
                  LevelData<FArrayBox> levelKE(levelGrids,1);

                  LevelData<FArrayBox> levelLambda(levelGrids,1);
                  thisLevelPtr->newLambda().copyTo(levelLambda.interval(),
                                                   levelLambda,
                                                   levelLambda.interval());

                  // since we want max(error(lambda)), subtract 1...
                  DataIterator levelDit = levelLambda.dataIterator();
                  for (levelDit.begin(); levelDit.ok(); ++levelDit)
                    {
                      levelLambda[levelDit()] -= 1.0;
                    }

                  int numVortComps = 1;
                  if (SpaceDim == 3) numVortComps = SpaceDim;
                  LevelData<FArrayBox> levelVort(levelGrids, numVortComps);

                  const ProblemDomain& levelDomain = thisLevelPtr->problemDomain();

                  if (lev < finest_level)
                    {
                      fineVelPtr = compVel[lev+1];
                      nRefFine = thisLevelPtr->m_ref_ratio;
                    }
                  else
                    {
                      fineVelPtr = NULL;
                      nRefFine = -1;
                    }

                  Divergence::compDivergenceCC(levelDiv, levelVel,
                                               crseVelPtr, fineVelPtr,
                                               thisLevelPtr->m_dx,
                                               nRefCrse, nRefFine,
                                               levelDomain, true);

                  thisLevelPtr->computeVorticity(levelVort);
                  thisLevelPtr->computeKineticEnergy(levelKE);

                  DataIterator dit = levelDiv.dataIterator();
                  Real levelMax = 0.0;
                  Real levelMaxVort = 0.0;
                  Real levelMaxLambda = 0.0;
                  Real levelIntegralLambda = 0.0;
                  Real levelSumKE = 0.0;
                  Interval divComps(0,0);
                  Interval vortComps(0,numVortComps-1);

                  if (fineVelPtr != NULL)
                    {
                      const DisjointBoxLayout* finerGridsPtr = &fineVelPtr->getBoxes();
                      levelMax = computeNorm(levelDiv, finerGridsPtr,
                                             nRefFine, thisLevelPtr->m_dx,
                                             divComps, 0);

                      levelMaxVort = computeNorm(levelVort, finerGridsPtr,
                                                 nRefFine, thisLevelPtr->m_dx,
                                                 vortComps,0);

                      levelMaxLambda = computeNorm(levelLambda, finerGridsPtr,
                                                   nRefFine,
                                                   thisLevelPtr->m_dx,
                                                   levelLambda.interval(), 0);

                      levelIntegralLambda = computeSum(levelLambda,
                                                       finerGridsPtr,
                                                       nRefFine,
                                                       thisLevelPtr->m_dx,
                                                       levelKE.interval());

                      levelSumKE = computeSum(levelKE, finerGridsPtr, nRefFine,
                                              thisLevelPtr->m_dx,
                                              levelKE.interval());
                    } // end if there is a finer level
                  else
                    {
                      const DisjointBoxLayout* finerGridsPtr = NULL;
                      levelMax = computeNorm(levelDiv, finerGridsPtr,
                                             nRefFine,
                                             thisLevelPtr->m_dx,
                                             divComps, 0);

                      levelMaxVort = computeNorm(levelVort, finerGridsPtr,
                                                 nRefFine,
                                                 thisLevelPtr->m_dx,
                                                 vortComps, 0);

                      levelMaxLambda = computeNorm(levelLambda,
                                                   finerGridsPtr,
                                                   nRefFine,
                                                   thisLevelPtr->m_dx,
                                                   levelLambda.interval(),0);

                      levelIntegralLambda = computeSum(levelLambda,
                                                       finerGridsPtr,
                                                       nRefFine,
                                                       thisLevelPtr->m_dx,
                                                       levelKE.interval());

                      levelSumKE = computeSum(levelKE, finerGridsPtr,
                                              nRefFine,
                                              thisLevelPtr->m_dx,
                                              levelKE.interval());
                    } // end if there is no finer level
                  if (levelMax > maxDiv) maxDiv = levelMax;
                  if (levelMaxVort > maxVort) maxVort = levelMaxVort;
                  if (levelMaxLambda > maxLambdaErr) maxLambdaErr = levelMaxLambda;
                  integralLambda += levelIntegralLambda;
                  integralKE += levelSumKE;

                  // now index to next level
                  nRefCrse = nRefFine;
                  crseVelPtr = &levelVel;
                  thisLevelPtr = thisLevelPtr->finerNSPtr();
                } // end loop  over levels

              pout() << setiosflags(ios::scientific) << setprecision(15);
              pout() << "Time = " << setw(15) << m_time
                     << setw(30) << " Max Div(u) = "
                     << setw(23)  << maxDiv << endl;
              pout() << "Time = " << setw(15) << m_time
                     << setw(30) << " Max(Vorticity) = "
                     << setw(23) << maxVort  << endl;
              pout() << "Time = " << setw(15) << m_time
                     << setw(30) << " Max(Err(Lambda) = "
                     << setw(23) << maxLambdaErr  << endl;
              pout() << "Time = " << setw(15) << m_time
                     << setw(30) << " Integral(Lambda) = "
                     << setw(23) << integralLambda << endl;
              pout() << "Time = " << setw(15) << m_time
                     << setw(30) << " Integral(Kinetic Energy) = "
                     << setw(23) << integralKE << endl;

              // this is something else we only do after full composite timestep
              if (s_compute_scal_err)
                {
                  thisLevelPtr = this;
                  Vector<LevelData<FArrayBox>* > scalErrVect(compVel.size(),
                                                             NULL);
                  Vector<int> AmrRefRatios(finest_level+1);

                  for (int scalComp = 0; scalComp<s_num_scal_comps; scalComp++)
                    {
                      // do this once for each scalar component
                      for (int lev=m_level; lev<=finest_level; lev++)
                        {
                          const DisjointBoxLayout& levelGrids =
                            compVel[lev]->getBoxes();
                          Real dxLevel = thisLevelPtr->m_dx;
                          const ProblemDomain& levelDomain
                            = thisLevelPtr->m_problem_domain;
                          AmrRefRatios[lev] = thisLevelPtr->refRatio();
                          const LevelData<FArrayBox>& levelScal
                            = thisLevelPtr->newScal(scalComp);
                          int ncomp = levelScal.nComp();
                          scalErrVect[lev]
                            = new LevelData<FArrayBox>(levelGrids,
                                                       ncomp);

                          // now compute exact scalar solution, then subtract
                          // existing solution to get error
                          LevelData<FArrayBox>& levelErr = *scalErrVect[lev];
                          int isEdgeCentered = 0;
                          DataIterator dit = levelScal.dataIterator();
                          for (dit.begin(); dit.ok(); ++dit)
                            {
                              FArrayBox& thisErr = levelErr[dit()];
                              FORT_INITSCAL(CHF_FRA(thisErr),
                                            CHF_FRA((*compVel[lev])[dit]),
                                            CHF_REAL(dxLevel),
                                            CHF_REAL(m_time),
                                            CHF_INT(scalComp),
                                            CHF_BOX(levelDomain.domainBox()),
                                            CHF_INT(isEdgeCentered));

                              thisErr.minus(levelScal[dit()],
                                            levelGrids[dit()],
                                            0, 0, ncomp);
                            } // end loop over grids on this level

                          // increment level pointer
                          thisLevelPtr = thisLevelPtr->finerNSPtr();
                        } // end loop over levels for error computation

                      // now compute norms of errors
                      // do this only for the 0th component of each scalar
                      Interval errInterval(0,0);
                      Real L1Err, L2Err, MaxErr;
                      int normType = 0;
                      int lBase = m_level;

                      MaxErr = computeNorm(scalErrVect, AmrRefRatios, m_dx,
                                           errInterval, normType, lBase);

                      normType = 1;
                      L1Err = computeNorm(scalErrVect, AmrRefRatios, m_dx,
                                          errInterval, normType, lBase);

                      normType = 2;
                      L2Err = computeNorm(scalErrVect, AmrRefRatios, m_dx,
                                          errInterval, normType, lBase);

                      // now report norms of errors
                      pout() << "scalar error, comp " << scalComp << ":"
                             << endl;
                      pout() << "Time = " << m_time << ",    L1 = "
                             << setw(11) << setprecision(4)
                             << setiosflags(ios::showpoint)
                             << setiosflags(ios::scientific)
                             << L1Err << ", L2 = "
                             << L2Err << ", Max(Err) = " << MaxErr << endl;

                      pout() << endl;

                      // clean up storage for error
                      for (int lev=m_level; lev<=finest_level; lev++)
                        {
                          if (scalErrVect[lev] != NULL)
                            {
                              delete scalErrVect[lev];
                              scalErrVect[lev] = NULL;
                            }
                        }
                    } // end loop over scalar components
                } // end if we compute scalar errors
            } // end if m_level is 0

          if (m_level > 0)
            {
              // need to clean up storage
              if (compLambda[m_level-1] != NULL)
                {
                  delete compLambda[m_level-1];
                  compLambda[m_level-1] = NULL;
                }

              if (compVel[m_level-1] != NULL)
                {
                  delete compVel[m_level-1];
                  compVel[m_level-1] = NULL;
                }
            }
        } // end sync projection operations

      // end if not finest level
    }
  else
    {
      // in case this is the only level
      if (m_level == 0)
        {
          // create boundary condition object and set physical BC's
          bool isViscous = (s_nu > 0.0);
          VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));

          const DisjointBoxLayout& grids = newVel().getBoxes();
          velBC.applyBCs(newVel(), grids,
                         m_problem_domain, m_dx,
                         false); // not homogeneous

          LevelData<FArrayBox>* crseVelPtr = NULL;
          LevelData<FArrayBox>* fineVelPtr = NULL;
          int nRefCrse = -1;
          int nRefFine = -1;

          Real maxDiv = 0.0;

          LevelData<FArrayBox> thisDiv(grids, 1);

          Interval velComps = newVel().interval();
          newVel().exchange(velComps);

          Divergence::compDivergenceCC(thisDiv, newVel(),
                                       crseVelPtr, fineVelPtr,
                                       m_dx, nRefCrse, nRefFine,
                                       m_problem_domain,
                                       true);

          DataIterator dit=thisDiv.dataIterator();
          Interval divComps(0,0);

          const DisjointBoxLayout* finerGridsPtr = NULL;
          maxDiv = computeNorm(thisDiv, finerGridsPtr, nRefFine,
                               m_dx, divComps, 0);

          pout() << setiosflags(ios::scientific) << setprecision(15);
          pout() << "Time = " << setw(15) << m_time
                 << setw(30) << " Max Div(u) = "
                 << setw(23)  << maxDiv << endl;

          // this is something else we only do after full composite timestep
          if (s_compute_scal_err)
            {
              for (int scalComp = 0; scalComp<s_num_scal_comps; scalComp++)
                {
                  // do this once for each scalar component
                  const LevelData<FArrayBox>& levelScal = newScal(scalComp);
                  const DisjointBoxLayout& levelGrids = levelScal.getBoxes();

                  int ncomp = levelScal.nComp();
                  LevelData<FArrayBox> levelErr(levelGrids, ncomp);

                  // now compute exact scalar solution, then subtract
                  // existing solution to get error
                  int isEdgeCentered = 0;
                  DataIterator dit = levelScal.dataIterator();
                  for (dit.begin(); dit.ok(); ++dit)
                    {
                      FArrayBox& thisErr = levelErr[dit()];
                      FORT_INITSCAL(CHF_FRA(thisErr),
                                    CHF_FRA(newVel()[dit]),
                                    CHF_REAL(m_dx),
                                    CHF_REAL(m_time),
                                    CHF_INT(scalComp),
                                    CHF_BOX(m_problem_domain.domainBox()),
                                    CHF_INT(isEdgeCentered));

                      thisErr.minus(levelScal[dit()], levelGrids[dit()],
                                    0, 0, ncomp);
                    } // end loop over grids on this level

                  // now compute norms of errors
                  // do this only for the 0th component of each scalar
                  Interval errInterval(0,0);
                  DisjointBoxLayout* finerGridsPtr = NULL;
                  int nRef = 0;
                  Real L1Err, L2Err, MaxErr;
                  int normType = 0;

                  MaxErr = computeNorm(levelErr, finerGridsPtr, nRef, m_dx,
                                       errInterval, normType);

                  normType = 1;
                  L1Err = computeNorm(levelErr, finerGridsPtr, nRef, m_dx,
                                      errInterval, normType);

                  normType = 2;
                  L2Err = computeNorm(levelErr, finerGridsPtr, nRef, m_dx,
                                      errInterval, normType);

                  // now report norms of errors
                  pout() << "scalar error, comp " << scalComp << ":" << endl;
                  pout() << "Time = " << m_time << ",    L1 = "
                         << setw(11) << setprecision(4)
                         << setiosflags(ios::showpoint)
                         << setiosflags(ios::scientific)
                         << L1Err << ", L2 = " << L2Err << ", Max(Err) = "
                         << MaxErr << endl;

                  pout() << endl;
                } // end loop over scalar components
            } // end if we compute scalar errors

        } // end if single-grid case...
    } // end if finest level
}
