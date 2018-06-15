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
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "FineInterp.H"
#include "SetValLevel.H"

#include "AMRNSF_F.H"

// ---------------------------------------------------------------
// tag cells for regridding
void
AMRNavierStokes::tagCells(IntVectSet & a_tags)
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::tagCells " << m_level << endl;
    }

  IntVectSet local_tags;

  // create tags based on something or other
  // for now, don't do anything
  const DisjointBoxLayout& level_domain = newVel().getBoxes();

  if (s_tag_vorticity)
    {
      LevelData<FArrayBox> vorticity;
      LevelData<FArrayBox> mag_vorticity;
      if (SpaceDim == 2)
        {
          vorticity.define(level_domain,1);
        }
      else if (SpaceDim == 3)
        {
          vorticity.define(level_domain,SpaceDim);
          mag_vorticity.define(level_domain,1);
        }
      computeVorticity(vorticity);

      Interval vortInterval(0,0);
      if (SpaceDim == 3)
        {
          vortInterval = Interval(0,SpaceDim-1);
        }
      Real tagLevel = norm(vorticity, vortInterval, 0);
      // Real tagLevel = 1.0;
      //tagLevel *= s_vort_factor/m_dx;
      // actually tag where vort*dx > s_vort_factor*max(vort)
      tagLevel = s_vort_factor/m_dx;

      if (tagLevel > 0)
        {
          // now tag where vorticity magnitude is greater than or equal
          // to tagLevel
          DataIterator dit = vorticity.dataIterator();
          for (dit.reset(); dit.ok(); ++dit)
            {
              FArrayBox& vortFab = vorticity[dit];

              // this only needs to be done in 3d...
              if (SpaceDim==3)
                {
                  FArrayBox& magVortFab = mag_vorticity[dit];
                  FORT_MAGVECT(CHF_FRA1(magVortFab,0),
                               CHF_CONST_FRA(vortFab),
                               CHF_BOX(level_domain[dit]));
                }

              BoxIterator bit(vortFab.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  if (SpaceDim == 2)
                    {
                      if (abs(vortFab(iv)) >= tagLevel)
                        {
                          local_tags |= iv;
                        }
                    }
                  else if (SpaceDim == 3)
                    {
                      FArrayBox& magVortFab = mag_vorticity[dit];
                      if (abs(magVortFab(iv)) >= tagLevel)
                        {
                          local_tags |= iv;
                        }
                    } // end if DIM=3
                } // end loop over interior of box
            } // loop over grids
        } // if taglevel > 0
    } // if tagging on vorticity

  a_tags = local_tags;
}

// ---------------------------------------------------------------
// create tags at initial time
void
AMRNavierStokes::tagCellsInit(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::tagCellsInit " << m_level << endl;
    }

  if (s_specifyInitialGrids)
    {
      Vector<Vector<Box> > amrGrids;

#ifdef CH_MPI
      MPI_Barrier(Chombo_MPI::comm);
      if (procID() == 0)
        {
#endif
          // read in predefined grids
          ifstream is(s_initialGridFile.c_str(), ios::in);
          if (is.fail())
            {
              pout() << "cannot open grids file " << s_initialGridFile << endl;
              MayDay::Error();
            }
          // format of file -- same as other grid files -- number of
          // levels, then for each level starting with level 1,
          // number of grids on level, list of boxes.
          int in_numLevels;
          is >> in_numLevels;
          CH_assert (in_numLevels >= m_level+1);
          if (s_verbosity >= 3)
            {
              pout() << "numLevels = " << in_numLevels << endl;
            }
          while (is.get() != '\n');
          amrGrids.resize(in_numLevels);
          // now loop over levels, starting with level 1
          int ngrid;
          for (int lev=1; lev< in_numLevels; lev++)
            {
              is >> ngrid;
              if (s_verbosity > 3)
                {
                  pout() << "level " << lev << " numGrids = " << ngrid << endl;
                  pout() << "Grids: ";
                }
              while (is.get() != '\n');
              amrGrids[lev].resize(ngrid);
              for (int i=0; i<ngrid; i++)
                {
                  Box bx;
                  is >> bx;
                  while (is.get() != '\n');
                  if (s_verbosity >= 3)
                    {
                      pout() << bx << endl;
                    }
                  amrGrids[lev][i] = bx;
                } // end loop over boxes on this level
            } // end loop over levels
#ifdef CH_MPI
        } // end if procID = 0
      MPI_Barrier(Chombo_MPI::comm);
      broadcast(amrGrids, 0);
      MPI_Barrier(Chombo_MPI::comm);
#endif

      // now coarsen next levels boxes down to this level and add them
      // to the tagged IntvectSet
      int nRefFine = refRatio();
      for (int n=0; n<amrGrids.size(); n++)
        {
          Box coarseFineBox(amrGrids[m_level+1][n]);
          coarseFineBox.coarsen(nRefFine);
          a_tags |= coarseFineBox;
        }
    }
  else
    {
      // use same function at initial time as in the rest of the time
      tagCells(a_tags);
    }
}

// ---------------------------------------------------------------
// create a load-balanced DisjointBoxLayout from a collection of boxes
DisjointBoxLayout
AMRNavierStokes::loadBalance(const Vector<Box>& a_grids)
{
  if (s_verbosity >= 3)
  {
    pout () << "AMRNavierStokes::loadBalance " << m_level
            << "  Nboxes=" << a_grids.size() << endl;
  }

  Vector<int> proc_map;
  LoadBalance(proc_map, a_grids);

  if (s_verbosity >= 4)
  {
    pout () << "AMRNavierStokes::loadBalance: processor map: " << endl;
    for (int igrid=0; igrid< a_grids.size(); ++igrid)
    {
      pout() << igrid << ": " << proc_map[igrid] << "  " << endl;
    }
    pout() << endl;
  }

  return (DisjointBoxLayout(a_grids, proc_map, m_problem_domain) );
}

// ---------------------------------------------------------------
// set up data on this level after regridding
void
AMRNavierStokes::regrid(const Vector<Box>& a_new_grids)
{
  if (s_verbosity >= 3)
  {
    pout () << "AMRNavierStokes::regrid " << m_level
            << " Nboxes=" << a_new_grids.size() << endl;
  }

  m_level_grids = a_new_grids;

  mortonOrdering(m_level_grids);
  const DisjointBoxLayout level_domain = loadBalance(a_new_grids);

  if (s_verbosity >= 4)
    {
      pout () << "new grids: " << endl;
      for (LayoutIterator lit = level_domain.layoutIterator(); lit.ok(); ++lit)
        {
          pout () << level_domain[lit] << endl;
        }
    }

  // save data for later copy -- only do this if level is not empty
  if (a_new_grids.size() > 0)
    {
      // also set flag to indicate that grids exist on this level
      m_is_empty = false;

      // may be necessary to set finestLevel flag
      // the way this works is that if the finer level is smpty, then
      // we define this level as the finest level; if it turns out that we
      // haven't gotten to that level yet, but it will have grids, then we
      // correct this in the finer level's regrid call
      if (finerNSPtr() != NULL)
        {
          if (finerNSPtr()->isEmpty())
            {
              finestLevel(true);
              crseNSPtr()->finestLevel(false);
            }
        } // end if finerNSPtr != Null
      else
        {
          // if finerNSPtr == NULL, then this is obviously the finest level
          finestLevel(true);
          crseNSPtr()->finestLevel(false);
        }

      if (s_smooth_after_regrid)
        {
          // first compute Laplacians of diffused quantities
          // for later use in smoothing ops.
          // make use of the fact that all levels should be
          // at the same time, so that we can store this info
          // in the old-data registers for later use

          if (!m_regrid_smoothing_done)
            {
              // since we need a composite laplacian here _and_
              // we need to do averaging down, do this all at
              // once at beginning of regrid cycle for velocity
              // and for all diffused scalars (anything that's hit
              // with a diffusion operator).
              // for now, do this the easy way using AMRSolvers (ouch!)
              // and doing composite ops.  In all cases, we store
              // the Laplacian in the old-data holder (since we
              // shouldn't need it until the next timestep -- all
              // levels should be at the same time)

              AMRNavierStokes* thisNSPtr = this;
              while (!thisNSPtr->finestLevel())
                {
                  thisNSPtr = thisNSPtr->finerNSPtr();
                }
              CH_assert (thisNSPtr->finestLevel());
              int finest_level = thisNSPtr->m_level;

              // this looks a lot like what we do in post-timestep
              thisNSPtr = this;
              // if this level > 0, we will actually need up to two
              // more coarser levels for C/F BCs (recall that we need to
              // apply the operator to the level from which we will be
              // interpolating data, so we need one level coarser than
              // that one for a C/F bc.
              int startLev = m_level;
              // this moves startLev to the level on which we need to apply
              // the operator
              if (startLev > 0)
                {
                  startLev = startLev -1;
                  thisNSPtr = thisNSPtr->crseNSPtr();
                }
              // this moves startLev to the level from which we get C/F BCs for
              // that level
              if (startLev > 0)
                {
                  startLev = startLev - 1;
                  thisNSPtr = thisNSPtr->crseNSPtr();
                }

              // now set up multilevel stuff
              Vector<LevelData<FArrayBox>* > amrS(finest_level+1,NULL);
              Vector<LevelData<FArrayBox>* > amrLapS(finest_level+1);
              Vector<DisjointBoxLayout> amrGrids(finest_level+1);
              Vector<int> amrRefRatios(finest_level+1,0);
              Vector<Real> amrDx(finest_level+1,0);
              Vector<ProblemDomain> amrDomains(finest_level+1);
              // will use this to do averaging down
              Vector<CoarseAverage*> amrAvgDown(finest_level+1, NULL);

              // loop over levels, allocate temp storage for velocities,
              // set up for amrsolves

              for (int lev=startLev; lev<=finest_level; lev++)
                {
                  const DisjointBoxLayout& levelGrids = thisNSPtr->newVel().getBoxes();
                  // since AMRSolver can only handle one component at a
                  // time, need to allocate temp space to copy stuff
                  // into, and then back out of to compute laplacian
                  IntVect ghostVect(D_DECL(1,1,1));
                  amrS[lev] = new LevelData<FArrayBox>(levelGrids, 1,
                                                       ghostVect);
                  amrLapS[lev] = new LevelData<FArrayBox>(levelGrids,1,
                                                          ghostVect);

                  amrGrids[lev] = levelGrids;
                  amrRefRatios[lev] = thisNSPtr->refRatio();
                  amrDx[lev] = thisNSPtr->Dx();
                  amrDomains[lev] = thisNSPtr->problemDomain();
                  thisNSPtr = thisNSPtr->finerNSPtr();
                  if (lev>= m_level)
                    {
                      amrAvgDown[lev] = new CoarseAverage(levelGrids,1,
                                                          amrRefRatios[lev-1]);
                    }
                } // end loop over levels

              AMRPoissonOpFactory localPoissonOpFactory;
              // allocate solver -- do this in a separate function in
              // order to make sure the two steps which involve
              // elliptic operators are consistent.

              // reset startLevel if necessary -- at this point, startLev
              // will be the coarsest level on which we apply the operator.
              if (m_level > 0) startLev=m_level-1;

              // note that we use startLev instead of m_level
              // as lBase for the solver because we will need
              // to compute operator on startLev
              defineRegridAMROp(localPoissonOpFactory, amrGrids,
                                amrDomains, amrDx,
                                amrRefRatios, startLev);

              int numLevs = finest_level - startLev + 1;
              Vector<AMRPoissonOp*> localPoissonOpPtrVec(numLevs);
              for (int lev=startLev; lev<=finest_level; lev++)
                {
                  // virtual  AMRLevelOp<LevelData< FArrayBox> >*
                  localPoissonOpPtrVec[lev - startLev] =
                    (AMRPoissonOp*)
                    localPoissonOpFactory.AMRnewOp(amrDomains[lev]);
                }

              // now loop over velocity components
              for (int dir=0; dir<SpaceDim; ++dir)
                {
                  Interval velComps(dir,dir);
                  Interval tempComps(0,0);
                  // for each velocity component, copy all levels
                  // of data into temp storage.  then, loop over
                  // all levels and compute laplacian of this comonent.
                  // finally, copy lap(vel) into old_vel storage space.
                  thisNSPtr = this;
                  if (startLev < m_level) thisNSPtr = thisNSPtr->crseNSPtr();

                  for (int lev=startLev; lev<=finest_level; lev++)
                    {
                      thisNSPtr->newVel().copyTo(velComps, *amrS[lev],
                                                 tempComps);
                      thisNSPtr = thisNSPtr->finerNSPtr();
                    }

                  // now loop over levels and apply multilevel operator.
                  // also do averaging down here if necessary

                  for (int lev=startLev; lev<=finest_level; lev++)
                    {
                      LevelData<FArrayBox>& levelS = *amrS[lev];
                      LevelData<FArrayBox>& levelLapS = *amrLapS[lev];
                      int indVec = lev - startLev;
                      if (lev == 0)
                        { // no coarser level
                          if (lev == finest_level)
                            { // no finer level:  this is all there is
                              localPoissonOpPtrVec[indVec]->
                                applyOp(levelLapS, levelS);
                            }
                          else
                            { // finer level exists
                              localPoissonOpPtrVec[indVec]->
                                AMROperatorNC(levelLapS,
                                              *amrS[lev + 1],
                                              levelS, false,
                                              localPoissonOpPtrVec[indVec+1]);
                            }
                        }
                      else if (lev < finest_level - 1)
                        { // three-level operator
                          localPoissonOpPtrVec[indVec]->
                            AMROperator(levelLapS,
                                        *amrS[lev + 1],
                                        levelS,
                                        *amrS[lev - 1], false,
                                        localPoissonOpPtrVec[indVec+1]);
                        }
                      else
                        { // no finer level
                          localPoissonOpPtrVec[indVec]->
                            AMROperatorNF(levelLapS,
                                          levelS,
                                          *amrS[lev - 1], false);
                        }
                    }

                  // need to do average down from finest level on down
                  // to speed this up eventually, may want to do this
                  // _after_ we copy into multicomponent velocity LDF
                  for (int lev= finest_level; lev>startLev; lev--)
                    {
                      amrAvgDown[lev]->averageToCoarse(*amrLapS[lev-1],
                                                       *amrLapS[lev]);
                    }

                  // finally, loop over levels and copy to old-vel storage
                  thisNSPtr = this;
                  if (m_level > 0) thisNSPtr = thisNSPtr->crseNSPtr();

                  for (int lev=startLev; lev<=finest_level; lev++)
                    {
                      amrLapS[lev]->copyTo(tempComps,thisNSPtr->oldVel(),
                                           velComps);
                      thisNSPtr = thisNSPtr->finerNSPtr();
                    }

                } // end loop over velocity components

              // since we're done with velocity, can clear up storage space
              for (int lev=0; lev<=finest_level; lev++)
                {
                  if (amrS[lev] != NULL)
                  {
                    delete amrS[lev];
                    amrS[lev] = NULL;
                  }
                  if (amrLapS[lev] != NULL)
                  {
                    delete amrLapS[lev];
                    amrLapS[lev] = NULL;
                  }
                }

              // now do scalars -- since (at the moment),
              // we only have single-component scalars, can
              // avoid copies and work directly with the scalars
              // themselves
              for (int scalComp=0; scalComp < s_num_scal_comps; scalComp++)
                {
                  // only do all this if scalar is diffused
                  if (s_scal_coeffs[scalComp] > 0)
                    {
                      thisNSPtr = this;
                      if (startLev < m_level)
                        {
                          thisNSPtr = thisNSPtr->crseNSPtr();
                          // may need even coarser level for C/F BCs
                          if (startLev > 0)
                            {
                              startLev = startLev -1;
                              thisNSPtr = thisNSPtr->crseNSPtr();
                            }
                        }

                      for (int lev=startLev; lev<=finest_level; lev++)
                        {
                          amrS[lev] = thisNSPtr->m_scal_new[scalComp];
                          amrLapS[lev] = thisNSPtr->m_scal_old[scalComp];
                          thisNSPtr = thisNSPtr->finerNSPtr();
                        }

                      if (m_level > 0) startLev = m_level-1;

                      // now compute laplacian and average down if necessary
                      for (int lev=startLev; lev<=finest_level; lev++)
                        {
                          LevelData<FArrayBox>& levelS = *amrS[lev];
                          LevelData<FArrayBox>& levelLapS = *amrLapS[lev];
                          int indVec = lev - startLev;
                          if (lev == 0)
                            { // no coarser level
                              if (lev == finest_level)
                                { // no finer level:  this is all there is
                                  localPoissonOpPtrVec[indVec]->
                                    applyOp(levelLapS, levelS);
                                }
                              else
                                { // finer level exists
                                  localPoissonOpPtrVec[indVec]->
                                    AMROperatorNC(levelLapS,
                                                  *amrS[lev + 1],
                                                  levelS, false,
                                                  localPoissonOpPtrVec[indVec+1]);
                                }
                            }
                          else if (lev < finest_level - 1)
                            { // three-level operator
                              localPoissonOpPtrVec[indVec]->
                                AMROperator(levelLapS,
                                            *amrS[lev + 1],
                                            levelS,
                                            *amrS[lev - 1], false,
                                            localPoissonOpPtrVec[indVec+1]);
                            }
                          else
                            { // no finer level
                              localPoissonOpPtrVec[indVec]->
                                AMROperatorNF(levelLapS,
                                              levelS,
                                              *amrS[lev - 1], false);
                            }
                          if (lev > m_level)
                            {
                              amrAvgDown[lev]->averageToCoarse(*amrLapS[lev-1],
                                                               *amrLapS[lev]);
                            }
                        } // end loop over levels

                      for (int lev=startLev; lev<=finest_level; lev++)
                        {
                          delete localPoissonOpPtrVec[lev - startLev];
                        }

                      // finally, multiply by mu, then add old scal
                      Real mu = s_regrid_smoothing_coeff*m_dt*s_scal_coeffs[scalComp];
                      for (int lev=startLev; lev <= finest_level; ++lev)
                        {
                          LevelData<FArrayBox>& levelS = *amrS[lev];
                          LevelData<FArrayBox>& levelLapS = *amrLapS[lev];
                          DataIterator levelDit = levelS.dataIterator();
                          for (levelDit.reset(); levelDit.ok(); ++levelDit)
                            {
                              levelLapS[levelDit] *= -mu;
                              levelLapS[levelDit] += levelS[levelDit];
                            }
                        } // end loop over levels
                    } // end if scalar is diffused
                } // end loop over scalar components

              // finally, clean up storage for coarseAverages
              for (int lev=0; lev<=finest_level; ++lev)
                {
                  if (amrAvgDown[lev] != NULL)
                    {
                      delete amrAvgDown[lev];
                      amrAvgDown[lev] = NULL;
                    }
                }
              thisNSPtr = this;
              if (startLev < m_level) thisNSPtr = thisNSPtr->crseNSPtr();
              for (int lev=startLev; lev <= finest_level; ++lev)
                {
                  thisNSPtr->m_regrid_smoothing_done = true;
                  thisNSPtr = thisNSPtr->finerNSPtr();
                }
            } // end if post-regridding smoothing hasn't been done yet
        } // end if we're doing post-regrid smoothing

      // since we're using pointers here, it's easy to save old
      // data -- then clean up afterwards
      LevelData<FArrayBox>* old_newVelPtr = m_vel_new_ptr;
      LevelData<FArrayBox>* old_oldVelPtr = m_vel_old_ptr;
      LevelData<FArrayBox>* old_newLambdaPtr = m_lambda_new_ptr;
      LevelData<FArrayBox>* old_oldLambdaPtr = m_lambda_old_ptr;
      Vector<LevelData<FArrayBox>*> old_oldScals(m_scal_old.size(),NULL);
      Vector<LevelData<FArrayBox>*> old_newScals(m_scal_new.size(),NULL);
#ifdef DEBUG
      LevelData<FArrayBox>* old_velSavePtr = m_vel_save_ptr;
      Vector<LevelData<FArrayBox>*> old_scalSave(m_scal_save.size(),NULL);
      LevelData<FArratBox>* old_lambdaSavePtr = m_lambda_save_ptr;
#endif
      // now loop over scalar components
      for (int comp=0; comp<m_scal_old.size(); comp++)
        {
          old_oldScals[comp] = m_scal_old[comp];
          old_newScals[comp] = m_scal_new[comp];
#ifdef DEBUG
          old_scalSave[comp] = m_scal_save[comp];
#endif
        }

      // reshape state with new grids
      IntVect ghostVect(D_DECL(1,1,1));
      m_vel_new_ptr = new LevelData<FArrayBox>(level_domain,
                                               s_num_vel_comps, ghostVect);
      m_vel_old_ptr = new LevelData<FArrayBox>(level_domain,
                                               s_num_vel_comps, ghostVect);

      m_lambda_new_ptr = new LevelData<FArrayBox>(level_domain,
                                                  1, ghostVect);
      m_lambda_old_ptr = new LevelData<FArrayBox>(level_domain,
                                                  1, ghostVect);

#ifdef DEBUG
      m_vel_save_ptr = new LevelData<FArrayBox>(level_domain,
                                                s_num_vel_comps,
                                                ghostVect);
      m_lambda_save_ptr = new LevelData<FArrayBox>(level_domain, 1,
                                                   ghostVect);
#endif

      m_scal_new.resize(s_num_scal_comps);
      m_scal_old.resize(s_num_scal_comps);
#ifdef DEBUG
      m_scal_save.resize(s_num_scal_comps);
#endif

      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          m_scal_new[comp] = new LevelData<FArrayBox>(level_domain,
                                                      s_compsPerScalar,
                                                      ghostVect);
          m_scal_old[comp] = new LevelData<FArrayBox>(level_domain,
                                                      s_compsPerScalar,
                                                      ghostVect);
#ifdef DEBUG
          m_scal_save[comp] = new LevelData<FArrayBox>(level_domain,
                                                       s_compsPerScalar,
                                                       ghostVect);
#endif
        }

      setAllBogus();

#ifdef DEBUG
      if (s_set_bogus_values)
        {
          setValLevel(*m_vel_save_ptr, s_bogus_value);
          setValLevel(*m_lambda_save_ptr, s_bogus_value);
          setValLevels(m_scal_save, 0, s_num_scal_comps-1, s_bogus_value);
        }
#endif

      // set up data structures
      levelSetup(level_domain);

      // interpolate from coarser level
      if (m_coarser_level_ptr != NULL)
        {
          AMRNavierStokes* amrns_ptr = dynamic_cast<AMRNavierStokes*> (m_coarser_level_ptr);
          if (amrns_ptr != NULL)
            {
              int nRefCrse = amrns_ptr->m_ref_ratio;
              FineInterp fine_interp(level_domain,s_num_vel_comps,nRefCrse,
                                     m_problem_domain);

              fine_interp.interpToFine (newVel(), amrns_ptr->newVel());

              fine_interp.interpToFine (oldVel(), amrns_ptr->oldVel());

#ifdef DEBUG
              fine_interp.interpToFine (*m_vel_save_ptr,
                                        *(amrns_ptr->m_vel_save_ptr));

#endif

              FineInterp fine_interp_lambda(level_domain, 1,
                                            nRefCrse, m_problem_domain);

              fine_interp_lambda.interpToFine(newLambda(), amrns_ptr->newLambda());

#ifdef DEBUG
              fine_interp_lambda.interpToFine(*m_lambda_save_ptr,
                                              *(amrns_ptr->m_lambda_save_ptr));
#endif

              for (int comp=0; comp<s_num_scal_comps; comp++)
                {
                  fine_interp_lambda.interpToFine(newScal(comp),
                                                  amrns_ptr->newScal(comp));
                }

            }
          else
            {
              MayDay::Error("in AMRNavierStokes::regrid: m_coarser_level is not castable to AMRNavierStokes*");
            }
        } // end if there is a coarser level

      // copy from old state
      if (old_newVelPtr != NULL)
        {
          old_newVelPtr->copyTo(old_newVelPtr->interval(),
                                newVel(),
                                newVel().interval());
        }

      if (old_oldVelPtr != NULL)
        {
          old_oldVelPtr->copyTo(old_oldVelPtr->interval(),
                                oldVel(),
                                oldVel().interval());
        }

      if (old_newLambdaPtr != NULL)
        {
          old_newLambdaPtr->copyTo(old_newLambdaPtr->interval(),
                                   newLambda(), newLambda().interval() );
        }

#ifdef DEBUG
      if (old_velSavePtr != NULL)
        {
          old_velSavePtr->copyTo(old_velSavePtr->interval(),
                                 *m_vel_save_ptr, m_vel_save_ptr->interval());
        }

      if (old_lambdaSavePtr != NULL)
        {
          old_lambdaSavePtr->copyTo(old_lambdaSavePtr->interval(),
                                    *m_lambda_save_ptr, m_lambda_save_ptr->interval());
        }
#endif

      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          if (old_oldScals[comp] != NULL)
            {
              old_oldScals[comp]->copyTo(old_oldScals[comp]->interval(),
                                         *m_scal_old[comp],
                                         m_scal_old[comp]->interval());
            }

          if (old_newScals[comp] != NULL)
            {
              old_newScals[comp]->copyTo(old_newScals[comp]->interval(),
                                         *m_scal_new[comp],
                                         m_scal_new[comp]->interval());
            }

#ifdef DEBUG
          if (old_scalSave[comp] != NULL)
            {
              old_scalSave[comp]->copyTo(old_scalSave[comp]->interval(),
                                         *m_scal_save[comp],
                                         m_scal_save[comp]->interval());
            }
#endif
        } // end loop over scalar components

      // clean up before these pointers go out of scope...
      if (old_newVelPtr != 0)
        {
          delete old_newVelPtr;
          old_newVelPtr = NULL;
        }

      if (old_oldVelPtr != NULL)
        {
          delete old_oldVelPtr;
          old_oldVelPtr = NULL;
        }

      if (old_newLambdaPtr != NULL)
        {
          delete old_newLambdaPtr;
          old_newLambdaPtr = NULL;
        }

      if (old_oldLambdaPtr != NULL)
        {
          delete old_oldLambdaPtr;
          old_oldLambdaPtr = NULL;
        }

#ifdef DEBUG
      if (old_velSavePtr != NULL)
        {
          delete old_velSavePtr;
          old_velSavePtr = NULL;
        }

      if (old_lambdaSavePtr != NULL)
        {
          delete old_lambdaSavePtr;
          old_lambdaSavePtr = NULL;
        }
#endif

      for (int comp=0; comp<s_num_scal_comps; comp++)
        {
          if (old_oldScals[comp] != NULL)
            {
              delete old_oldScals[comp];
              old_oldScals[comp]  = NULL;
            }
          if (old_newScals[comp] != NULL)
            {
              delete old_newScals[comp];
              old_newScals[comp] = NULL;
            }
#ifdef DEBUG
          if (old_scalSave[comp] != NULL)
            {
              delete old_scalSave[comp];
              old_scalSave[comp] = NULL;
            }
#endif
        } // end loop over components for scalars
    }
  else // new level empty
    {
      // if new level is empty, then just clear everything

      // should just be able to delete pointers here
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

#ifdef DEBUG
      if (m_vel_save_ptr != NULL)
        {
          delete m_vel_save_ptr;
          m_vel_save_ptr = NULL;
        }
#endif

      for (int comp=0; comp<m_scal_new.size(); comp++)
        {
          if (m_scal_new[comp] != NULL)
            {
              delete m_scal_new[comp];
              m_scal_new[comp] = NULL;
            }
          if (m_scal_old[comp] != NULL)
            {
              delete m_scal_old[comp];
              m_scal_old[comp] = NULL;
            }

#ifdef DEBUG
          if (m_scal_save[comp] != NULL)
            {
              delete m_scal_save[comp];
              m_scal_save[comp] = NULL;
            }
#endif
        } // end loop over scalar components

      // this also implies that this is no longer the finest level
      finestLevel(false);
      // if a coarser level exists, let it know that it's the finest level
      if (m_coarser_level_ptr != NULL)
        {
          AMRNavierStokes* crse_amrns_ptr = dynamic_cast<AMRNavierStokes*> (m_coarser_level_ptr);
          if (crse_amrns_ptr != NULL)
            {
              if (!crse_amrns_ptr->isEmpty()) crse_amrns_ptr->finestLevel(true);
            }
        }
      else
        {
          MayDay::Error("in AMRNavierStokes::regrid: m_coarser_level is not castable to AMRNavierStokes*");
        }
    } // new level empty
}

// ---------------------------------------------------------------
void
AMRNavierStokes::postRegrid(int a_lBase)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::postRegrid " << m_level << endl;
    }

  // this looks remarkably similar to postInitialize
  // need to re-project velocity and recompute e_lambda;
  // do this from the finest level

  if (finestLevel() && (m_level != a_lBase))
    {

      int numLevels = m_level+1;
      Vector<LevelData<FArrayBox>* > amrVel(numLevels);
      Vector<LevelData<FArrayBox>* > amrLambda(numLevels);

      AMRNavierStokes* thisLevelData = this;
      for (int lev=m_level; lev >=a_lBase; --lev)
        {
          amrVel[lev] = &(thisLevelData->newVel());
          amrLambda[lev] = &(thisLevelData->newLambda());

          if (lev > a_lBase)
            {
              thisLevelData = dynamic_cast<AMRNavierStokes*> (thisLevelData->m_coarser_level_ptr);
            }
        }

      CH_assert(thisLevelData->m_level == a_lBase);

      AMRNavierStokes* lBaseAMRNavierStokes = thisLevelData;

      CCProjector& lBaseProj = thisLevelData->m_projection;
      Real dt_lBase = thisLevelData->m_dt;

      if (s_smooth_after_regrid)
        {
          lBaseAMRNavierStokes->smoothVelocityField(a_lBase);
        }

      // set physical boundary conditions on velocity
      bool isViscous = (s_nu > 0);
      VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));

      for (int lev=a_lBase; lev <= m_level; lev++)
        {
          const ProblemDomain& levelDomain = thisLevelData->problemDomain();
          Real levelDx = thisLevelData->Dx();
          const DisjointBoxLayout& levelGrids = amrVel[lev]->getBoxes();
          velBC.applyBCs(*amrVel[lev], levelGrids,
                         levelDomain, levelDx,
                         false); // inhomogeneous

          if (thisLevelData->m_finer_level_ptr != NULL)
            {
              thisLevelData = dynamic_cast<AMRNavierStokes*> (thisLevelData->m_finer_level_ptr);
            }
        } // end loop over levels

      // may need coarser CF-BC
      LevelData<FArrayBox> crseVel;
      LevelData<FArrayBox> crseLambda;
      if (a_lBase > 0)
        {
          AMRNavierStokes* crseLevelPtr = dynamic_cast<AMRNavierStokes*> (lBaseAMRNavierStokes->m_coarser_level_ptr);
          const DisjointBoxLayout& crseGrids = crseLevelPtr->newVel().getBoxes();
          // this shouldn't need any ghost cells (proper nesting)
          crseVel.define(crseGrids, SpaceDim);

          crseLevelPtr->fillVelocity(crseVel, m_time);
          crseLevelPtr->fillLambda(crseLambda, m_time);

          amrVel[a_lBase-1] = &crseVel;
          amrLambda[a_lBase-1] = &crseLambda;
        }

      // note that initial velocity projection is done holding
      // lBase velocity constant (so we call the finer-level projection)
      bool homogeneousCFBC = true;
      lBaseProj.fineProjPtr()->initialVelocityProject(amrVel,homogeneousCFBC);

      // compute freestream preservation correction on lbase as well
      // (since grad(eLambda) on lbase accounts for the grids on
      // lbase+1
      lBaseProj.doPostRegridOps(amrVel, amrLambda,
                                dt_lBase, m_time);
      thisLevelData = lBaseAMRNavierStokes;

      // need to reset boundary conditions here
      for (int lev=a_lBase; lev<=m_level; lev++)
        {
          const ProblemDomain& levelDomain = thisLevelData->problemDomain();
          Real levelDx = thisLevelData->Dx();
          const DisjointBoxLayout& levelGrids = amrVel[lev]->getBoxes();
          velBC.applyBCs(*amrVel[lev], levelGrids,
                         levelDomain, levelDx,
                         false); // inhomogeneous

          if (thisLevelData->m_finer_level_ptr != NULL)
            {
              thisLevelData = dynamic_cast<AMRNavierStokes*> (thisLevelData->m_finer_level_ptr);
            }
        } // end loop over levels to reset physBC's

      if (s_initialize_pressures)
        {
          lBaseAMRNavierStokes->initializeGlobalPressure();
        }
    } // end if finest level
}
