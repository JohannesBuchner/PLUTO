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
#include "EdgeToCell.H"
#include "Gradient.H"
#include "AMRPoissonOpF_F.H"
#include "AMRNSF_F.H"
#include "AdvectUtil.H"
#include "AdvectUtilF_F.H"
#include "AdvectionPhysics.H"
#include "VelBCHolder.H"
#include "SetValLevel.H"

// ---------------------------------------------------------------
void
AMRNavierStokes::predictVelocities(LevelData<FArrayBox>& a_uDelU,
                                   LevelData<FluxBox>& a_advVel)
{
  if (s_verbosity >= 3)
    {
      pout () << "AMRNavierStokes::predictVelocities: " << m_level << endl;
    }

  bool isViscous = (s_nu > 0.0);

  const DisjointBoxLayout& levelGrids = a_advVel.getBoxes();

  // for tracing, will need to fill in boundary values for
  // grown copy of velocity
  IntVect ghostVect(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));
  LevelData<FArrayBox> viscousSource(levelGrids, SpaceDim, ghostVect);

  if (s_set_bogus_values)
    {
      setValLevel(viscousSource, s_bogus_value);
    }

  LevelData<FArrayBox> traceVel(levelGrids, SpaceDim, ghostVect);
  if (s_set_bogus_values)
    {
      setValLevel(traceVel, s_bogus_value);
    }

  // now fill in velocity field
  // m_time is new time
  Real old_time = m_time - m_dt;
  fillVelocity(traceVel, old_time);

  // do physical boundary conditions on traceVel

  if (isViscous)
    {
      // rather than resetting BC's on old_vel, and then setting them
      // back, just copy old_vel to a temporary holder, and set
      // BCs on that one.
      LevelData<FArrayBox> viscousVel(levelGrids, SpaceDim, ghostVect);
      DataIterator dit = viscousVel.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          viscousVel[dit].copy(traceVel[dit]);
        }

      {
        VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());
        velBC.applyBCs(viscousVel, levelGrids,
                       m_problem_domain, m_dx,
                       false); // not homogeneous
      }

      // if crse level exists, fill coarse velocity BC for laplacian
      LevelData<FArrayBox>* crseVelPtr = NULL;
      if (m_level > 0)
        {
          const DisjointBoxLayout& crseGrids = crseNSPtr()->newVel().getBoxes();
          crseVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
          crseNSPtr()->fillVelocity(*crseVelPtr, old_time);
        }

      computeLapVel(viscousSource, viscousVel, crseVelPtr);

      if (crseVelPtr != NULL)
        {
          delete crseVelPtr;
          crseVelPtr = NULL;
        }

      for (dit.reset(); dit.ok(); ++dit)
        {
          viscousSource[dit].mult(s_nu);
        }

    } // end if viscous
  else
    {
      // inviscid case
      setValLevel(viscousSource, 0.0);
    } // end if inviscid

  // tracing uses inviscid BC's
  {
    VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
    velBC.applyBCs(traceVel, levelGrids,
                   m_problem_domain, m_dx,
                   false); // not homogeneous
  }

  // will need edge-centered storage for all velocity components
  LevelData<FluxBox> uHalf(levelGrids, SpaceDim);

  { //this is to allow this fluxbox to go out of scope
    // will also need grad(eLambda) and grad(phi_mac)
    LevelData<FluxBox>& grad_eLambda = m_projection.grad_eLambda();
    LevelData<FluxBox> gradPhi(levelGrids, SpaceDim);
    m_projection.gradPhi(gradPhi);

    if (s_set_bogus_values)
      {
        setValLevel(uHalf, s_bogus_value);
      }

    computePredictedVelocities(uHalf, traceVel, a_advVel,
                               viscousSource,
                               m_patchGodVelocity,
                               grad_eLambda,
                               gradPhi,
                               s_applyFreestreamCorrection,
                               old_time, m_dt);


  } // FluxBox goes out of scope, memory reclaimed.

  // now loop over grids to compute advection and update flux registers
  DataIterator dit = a_uDelU.dataIterator();

  // for nonconservative update, also need a cell-centered
  // advection velocity -- get this by averaging edge-centered
  // advection velocity to cell centers.  this one shouldn't need
  // any ghost cells.

  // reduce, reuse, recycle...
  //LevelData<FArrayBox> cellAdvVel(levelGrids, SpaceDim);
  LevelData<FArrayBox>&  cellAdvVel = viscousSource;

  EdgeToCell(a_advVel, cellAdvVel);

  for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& thisUHalf = uHalf[dit];
      FArrayBox& thisCellAdvVel = cellAdvVel[dit];
      const Box& thisBox = levelGrids[dit];
      FArrayBox& this_uDelU = a_uDelU[dit];
      this_uDelU.setVal(0.0);

      // to do this in a dimensionality independent way,
      // loop over directions.
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& uHalfDir = thisUHalf[dir];

          // this will increment this_uDelU with
          // cellVel*d(uHalfDir)/d[xyz] for all velocity components
          FORT_UDELS(CHF_FRA(this_uDelU),
                     CHF_FRA(uHalfDir),
                     CHF_FRA(thisCellAdvVel),
                     CHF_BOX(thisBox),
                     CHF_REAL(m_dx),
                     CHF_INT(dir));
        }

      // now need to do flux register stuff
      // first compute fluxes (u_edge*u_adv)
      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& this_uEdge = uHalf[dit][dir];
          FArrayBox& thisAdvVel = a_advVel[dit][dir];

          // now need to loop over traced velocity components
          for (int velComp=0; velComp<SpaceDim; ++velComp)
            {
              this_uEdge.mult(thisAdvVel,0,velComp,1);
            }
        }

      Interval velComps(0,SpaceDim-1);
      if (!finestLevel())
        {
          CH_assert (m_flux_register.isDefined());
          for (int dir=0; dir < SpaceDim; dir++)
            {
              // because we eventually want -D_R(U-del-U), signs
              // on fluxMult are opposite of what you might expect...
              Real fluxMult = m_dt; // petermc made positive, 7 Dec 07
              if (s_reflux_normal_momentum)
                {
                  m_flux_register.incrementCoarse(uHalf[dit][dir],
                                                  fluxMult, dit(), velComps,
                                                  velComps, dir);
                }
              else
                {
                  // if we're not refluxing normal momentum component,
                  // then do this component by component
                  for (int velComp=0; velComp<SpaceDim; velComp++)
                    {
                      if (velComp != dir)
                        {
                          Interval velInt(velComp, velComp);
                          m_flux_register.incrementCoarse(uHalf[dit][dir],
                                                          fluxMult, dit(),
                                                          velInt,
                                                          velInt, dir);
                        } // end if a tangential velocity component
                    } // end loop over velocity components
                } // end if we're not refluxing normal momentum
            } // end loop over directions
        } // end if a finer level exists

      if (m_level > 0)
        {
          LevelFluxRegister& crseFR = crseNSPtr()->m_flux_register;
          CH_assert(crseFR.isDefined());
          // because we will eventually want -D_R(U-del-U) signs
          // on fluxMult are opposite of what you might expect...
          Real fluxMult = m_dt; // petermc made positive, 7 Dec 07
          for (int dir=0; dir<SpaceDim; dir++)
            {
              if (s_reflux_normal_momentum)
                {
                  crseFR.incrementFine(uHalf[dit][dir],
                                       fluxMult, dit(), velComps,
                                       velComps, dir, Side::Lo);

                  crseFR.incrementFine(uHalf[dit][dir],
                                       fluxMult, dit(), velComps,
                                       velComps, dir, Side::Hi);
                }
              else
                {
                  // for case where we only want to reflux tangential
                  // momentum components, do this component by component
                  // so that we can isolate the normal component.
                  for (int velComp = 0; velComp<SpaceDim; velComp++)
                    {
                      if (velComp != dir)
                        {
                          Interval velInt(velComp, velComp);

                          crseFR.incrementFine(uHalf[dit][dir],
                                               fluxMult, dit(), velInt,
                                               velInt, dir, Side::Lo);

                          crseFR.incrementFine(uHalf[dit][dir],
                                               fluxMult, dit(), velInt,
                                               velInt, dir, Side::Hi);
                        } // end if this is a tangential component
                    } // end loop over velocity components
                } // end if we're not refluxing normal momentum
            } // end loop over directions
        } // end if level > 0
    } // end loop over grids
  // at this point, we've computed the advective term, and updated
  // the flux registers with the momentum flux
}


// ---------------------------------------------------------------
// the approximation to u-dot-del-u is passed into this
// function in a_uStar -- at end of function, a_uStar
// will contain approximation to new velocity:
// uStar + dt*gradPi
// CALLS LEVELTGA
void
AMRNavierStokes::computeUStar(LevelData<FArrayBox>& a_uStar)
{
  bool isViscous = (s_nu > 0.0);

  // -> input:  phiOld
  // -> output:  phiNew

  // -> need to call updateSoln once for each dimension.

  DataIterator dit = a_uStar.dataIterator();

  if (!isViscous)
    {
      /// compute uStar = oldVel - dt*(u-dot-del-u)
      for (dit.reset(); dit.ok(); ++dit)
        {
          Real mult = -1.0*m_dt;
          a_uStar[dit] *= mult;
          a_uStar[dit] += oldVel()[dit];
        }
    }  // end if inviscid
  else
    {
      // viscous case is more elaborate -- need to do viscous solves.
      const DisjointBoxLayout& levelGrids = a_uStar.getBoxes();

      // m_time is new time
      Real old_time = m_time - m_dt;

      // define temporary storage
      LevelData<FArrayBox> gradPi(levelGrids, SpaceDim);

      // first get gradPi
      m_projection.gradPi(gradPi);

      // entering computeUStar(), a_uStar is u-dot-del-u;
      // leaving computeUStar(), a_uStar should be new velocity.
      // "f" is -u-del-u -gradPi +source
      LevelData<FArrayBox> src(levelGrids, SpaceDim, a_uStar.ghostVect());
      for (dit.reset(); dit.ok(); ++dit)
        {
          //          a_uStar[dit] += gradPi[dit];
          //          a_uStar[dit] *= -1.0;
          src[dit].copy(a_uStar[dit]);
          src[dit] += gradPi[dit];
          src[dit] *= -1.0;
        }

      int numberMGlevels = (m_level == 0) ? 0 : 1;
      pout() << "computeUStar at level " << numberMGlevels << endl;
      // -> updateSoln():  a_src = a_uStar here
      // update fine flux registers if and only if NOT finest level
      LevelFluxRegister* fineFluxRegisterPtr = NULL;
      if (!finestLevel()) fineFluxRegisterPtr = &m_flux_register;
      Real bogusTime = 0.;
      Real crseNewTime = bogusTime;
      Real crseOldTime = bogusTime;
      if (numberMGlevels > 0)
        {
          crseNewTime = crseNSPtr()->m_time;
          Real crseDt = crseNSPtr()->m_dt;
          crseOldTime = crseNewTime - crseDt;
        }

      LevelData<FArrayBox> oldVelocity(levelGrids, SpaceDim, IntVect::Unit);
      fillVelocity(oldVelocity, old_time);
      {
        // VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());
        VelBCHolder velBC(m_physBCPtr->uDelUFuncBC(true));
        velBC.applyBCs(oldVelocity, levelGrids,
                       m_problem_domain, m_dx,
                       false); // not homogeneous
      }

      for (int comp = 0; comp < SpaceDim; comp++)
        {
          Interval intvl(comp, comp);

          // set convergence metric to be norm of velocity field
          int normType = 0;
          Real velNorm = norm(oldVel(), intvl, normType);
          // m_velMGsolverPtr should be passed into LevelTGA
          m_velMGsolverPtrs[comp]->m_convergenceMetric = velNorm;

          LevelData<FArrayBox> compNewVelocity;
          aliasLevelData(compNewVelocity, &a_uStar, intvl);

          LevelData<FArrayBox> compOldVelocity;
          aliasLevelData(compOldVelocity, &oldVelocity, intvl);

          LevelData<FArrayBox> compSrc;
          aliasLevelData(compSrc, &src, intvl);

          pout() << "calling updateSoln on component " << comp << endl;
          if (numberMGlevels == 0)
            {
              m_TGAsolverPtrs[comp]->updateSoln(compNewVelocity,
                                                compOldVelocity,
                                                compSrc,
                                                fineFluxRegisterPtr, NULL,
                                                NULL, NULL,
                                                old_time, crseOldTime, crseNewTime, m_dt,
                                                numberMGlevels, false, comp); // do not initialize to zero
            }
          else // numberMGlevels == 1
            {
              LevelData<FArrayBox> compCrseOldVel;
              LevelData<FArrayBox> compCrseNewVel;
              aliasLevelData(compCrseOldVel,
                             crseNSPtr()->m_vel_old_ptr, intvl);

              aliasLevelData(compCrseNewVel,
                             crseNSPtr()->m_vel_new_ptr, intvl);
              m_TGAsolverPtrs[comp]->updateSoln(compNewVelocity,
                                                compOldVelocity,
                                                compSrc,
                                                fineFluxRegisterPtr,
                                                &(crseNSPtr()->m_flux_register),
                                                &compCrseOldVel,
                                                &compCrseNewVel,
                                                old_time, crseOldTime, crseNewTime, m_dt,
                                                numberMGlevels, false, comp); // do not initialize to zero
            }
        }

      // final thing to do -- subtract off dt*gradPi term
      for (dit.reset(); dit.ok(); ++dit)
        {
          gradPi[dit] *= m_dt;
          a_uStar[dit] += gradPi[dit];
        }
    } // end if viscous
}
