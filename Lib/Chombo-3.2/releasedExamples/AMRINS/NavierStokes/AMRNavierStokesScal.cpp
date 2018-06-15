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
#include "Divergence.H"
#include "Gradient.H"
#include "CellToEdge.H"
#include "timeInterp.H"
#include "AdvectUtil.H"

// this is necessary because we call library fortran directly
#include "AMRPoissonOpF_F.H"

#include "AdvectionPhysics.H"

// ---------------------------------------------------------------
void
AMRNavierStokes::advectScalar(LevelData<FArrayBox>& a_new_scal,
                              LevelData<FArrayBox>& a_old_scal,
                              LevelData<FluxBox>&   a_advVel,
                              LevelFluxRegister*    a_crseFluxRegPtr,
                              LevelFluxRegister&    a_fineFluxReg,
                              PatchGodunov&         a_patchGodScalar,
                              Real                  a_dt)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::advectScalar: " << m_level << endl;
    }

  const DisjointBoxLayout& levelGrids = a_new_scal.getBoxes();

  /// need to build grown grids to get be able to do all of
  /// tracing properly -- this is necessary to sync with old
  /// code
  IntVect advect_grow(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));

  LevelData<FArrayBox> local_old_vel(levelGrids, SpaceDim, advect_grow);
  // m_time contains the time at which the new state is centered
  Real old_time = m_time - m_dt;
  fillVelocity(local_old_vel, old_time);

  DataIterator dit = a_new_scal.dataIterator();
  dit.reset();
  int nScal = a_new_scal.nComp();
  Interval scalComps(0,nScal-1);

  IntVect ghostVect = IntVect::Unit;
  LevelData<FluxBox> edgeScal(levelGrids, nScal, ghostVect);
  LevelData<FArrayBox> src(levelGrids, nScal, a_old_scal.ghostVect());
  for (dit.begin(); dit.ok(); ++dit)
    {
      src[dit].setVal(0.0);
    }

  // use PatchGodunov::updateState function to do scalar advection

  // get AdvectionPhysics object w/in the PatchGodunov object
  AdvectionPhysics* advectionPhysics =
    dynamic_cast<AdvectionPhysics*>(a_patchGodScalar.getGodunovPhysicsPtr());
  if (advectionPhysics == NULL)
    {
      MayDay::Error("AMRNavierStokes::advectDiffuseScalar - unable to upcast GodunovPhysics to AdvectionPhysics");
    }

  // set up patchGodunov for this problem
  a_patchGodScalar.setCurrentTime(old_time);

  // also need to build a grown advection velocity
  IntVect advectVelGrow(IntVect::Unit);
  LevelData<FluxBox> grownAdvVel(levelGrids, 1, advectVelGrow);
  // this is inefficient, but we'll see how bad it really is
  // for now, fill extra cells with cell-to-face averaged velocities
  CellToEdge(local_old_vel, grownAdvVel);
  // now overwrite with advection velocities wherever possible
  a_advVel.copyTo(a_advVel.interval(), grownAdvVel,
                  grownAdvVel.interval());

  // copy old-time scalar to new-time holder
  a_old_scal.copyTo(a_new_scal);

  // now compute update

  Real maxWaveSpeed = 0.0;
  for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& thisEdgeScal = edgeScal[dit()];
      FluxBox& thisAdvVel = grownAdvVel[dit()];
      FArrayBox& thisCellVel = local_old_vel[dit()];
      FArrayBox& thisOldScal = a_old_scal[dit()];
      FArrayBox& thisNewScal = a_new_scal[dit()];
      FArrayBox& thisSrc = src[dit()];

      // PatchGodunov-based approach to advection
      const Box& curBox = levelGrids[dit()];
      a_patchGodScalar.setCurrentBox(curBox);
      // set cell-centered velocity field
      advectionPhysics->setCellVelPtr(&thisCellVel);
      // set advection velocity field
      advectionPhysics->setAdvVelPtr(&thisAdvVel);

      // compute updated scalar -- note that this updates
      // thisOldScal in place
      a_patchGodScalar.updateState(thisOldScal, thisEdgeScal,
                                   maxWaveSpeed,
                                   thisSrc, a_dt, curBox);
      // now copy oldScal->newScal (need to do it this way because
      // oldScal has the number of ghost cells needed for this
      // computation, while newScal most likely doesn't.
      thisNewScal.copy(thisOldScal, thisNewScal.box());
    }

  // now do flux register stuff
  for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& thisFlux = edgeScal[dit];

      Real fluxMult = a_dt; // petermc made positive, 7 Dec 07; still 13 Dec??
      if (!finestLevel())
        {
          CH_assert(a_fineFluxReg.isDefined());
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // since we will want the increment to be -U-del-U,
              // make this mult negative
              a_fineFluxReg.incrementCoarse(thisFlux[dir], fluxMult, dit(),
                                            scalComps, scalComps, dir);
            }
        } // end if a finer level exists

      if (m_level > 0)
        {
          CH_assert(a_crseFluxRegPtr != NULL);
          CH_assert(a_crseFluxRegPtr->isDefined());
          for (int dir=0; dir<SpaceDim; dir++)
            {
              a_crseFluxRegPtr->incrementFine(thisFlux[dir], fluxMult, dit(),
                                              scalComps, scalComps, dir,
                                              Side::Lo);

              a_crseFluxRegPtr->incrementFine(thisFlux[dir], fluxMult, dit(),
                                              scalComps, scalComps, dir,
                                              Side::Hi);
            }
        } // end if coarser level exists
    } // end loop over grids

  // that should be it!
}


// ---------------------------------------------------------------
void
AMRNavierStokes::
advectDiffuseScalar(
                    LevelData<FArrayBox>&       a_new_scal,
                    LevelData<FArrayBox>&       a_old_scal,
                    LevelData<FluxBox>&         a_adv_vel,
                    const Real&                 a_diffusive_coeff,
                    const LevelData<FArrayBox>* a_crse_scal_old,
                    const LevelData<FArrayBox>* a_crse_scal_new,
                    const Real                  a_old_crse_time,
                    const Real                  a_new_crse_time,
                    LevelFluxRegister*          a_crseFluxRegPtr,
                    LevelFluxRegister&          a_fineFluxReg,
                    PatchGodunov&               a_patchGodScalar,
                    BCHolder&                   a_scalPhysBC,
                    Real                        a_dt,
                    int                         a_comp)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRNavierStokes::advectDiffuseScalar: " << m_level << endl;
    }

  //check to see if we need diffusive stuff or not...
  bool isDiffusive = (a_diffusive_coeff > 0);

  if (!isDiffusive)
    { // not diffusive: just turn around and call nondiffusive scalar advection

      advectScalar(a_new_scal, a_old_scal, a_adv_vel, a_crseFluxRegPtr,
                   a_fineFluxReg, a_patchGodScalar,
                   a_dt);
    }
  else
    { // otherwise, do diffusive algorithm

      int numScal = a_new_scal.nComp();
      Interval scalComps(0, numScal-1);
      Real old_time = m_time - m_dt;

      // at the moment, can only do one component at a time for diffusion
      CH_assert (numScal == 1);

      const DisjointBoxLayout& levelGrids = a_new_scal.getBoxes();

      // define temporary storage
      LevelData<FArrayBox> diffusiveSrc(levelGrids,numScal,
                                        IntVect::Unit);
      LevelData<FArrayBox> tempStorage(levelGrids,numScal,
                                       IntVect::Unit);
      LevelData<FArrayBox>* crseDataPtr = NULL;

      if (m_level > 0)
        {
          // allocate crseBC info
          const DisjointBoxLayout& crseGrids = crseNSPtr()->newVel().getBoxes();
          crseDataPtr = new LevelData<FArrayBox>(crseGrids,1);
        }

      DataIterator dit = a_new_scal.dataIterator();

      if (s_set_bogus_values)
        {
          for (dit.reset(); dit.ok(); ++dit)
            {
              diffusiveSrc[dit].setVal(s_bogus_value);
              tempStorage[dit].setVal(s_bogus_value);
            }
        } // end bogus-value initialization

      // now compute viscous source term
      // copy scalar into temp storage so we don't overwrite
      // C/F BC's

      // define our own copier here to prevent ghost cells
      // being filled (avoiding a potentially expensive exchange)
      // since they will be filled a bit later when we compute the Laplacian
      Copier thisCopier(a_old_scal.getBoxes(), tempStorage.getBoxes(),
                        IntVect::Zero);
      a_old_scal.copyTo(a_old_scal.interval(), tempStorage,
                        tempStorage.interval(), thisCopier);
      // now set up crse level BC
      if (crseDataPtr != NULL)
        {
          CH_assert (a_crse_scal_old != NULL);
          CH_assert (a_crse_scal_new != NULL);
          // need to interpolate in time
          timeInterp(*crseDataPtr, old_time,
                     *a_crse_scal_old, a_old_crse_time,
                     *a_crse_scal_new, a_new_crse_time,
                     scalComps);
        }
      // exchange for tempStorage (scal) not necessary since it's
      // done in PoissonOp::ApplyOpI

      computeLapScal(diffusiveSrc, tempStorage, a_scalPhysBC, crseDataPtr);

      // exchange for diffusiveSrc not necessary because it's done
      // in computeLapScal

      for (dit.reset(); dit.ok(); ++dit)
        {
          FArrayBox& thisSrc = diffusiveSrc[dit];
          thisSrc *= a_diffusive_coeff;
        }

      // this is pretty much identical to advectScalar()
      /// need to build grown grids to get be able to do all of
      /// tracing properly -- this is necessary to sync with old
      /// code
      IntVect advect_grow(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));

      LevelData<FArrayBox> local_old_vel(levelGrids, SpaceDim, advect_grow);
      // m_time contains the time at which the new state is centered
      fillVelocity(local_old_vel, old_time);

      IntVect ghostVect = IntVect::Unit;
      LevelData<FluxBox> edgeScal(levelGrids, numScal, ghostVect);

      computeScalarAdvectiveFlux(edgeScal, a_old_scal, a_adv_vel,
                                 local_old_vel,
                                 diffusiveSrc, a_patchGodScalar,
                                 old_time, m_dt);

      LevelData<FArrayBox> src(levelGrids, numScal, a_new_scal.ghostVect());

      // compute div(us)
      Divergence::levelDivergenceMAC(src, edgeScal, m_dx);

      // now begin to set up viscous solves...
      //"f" term is -u-del-s, for scalars
      for (dit.reset(); dit.ok(); ++dit)
        {
          src[dit] *= -1.0;
        }
      // now compute diffused version of this

      // now add advective flux to diffusive flux
      // (edgeScal already contains fluxes)
      // advective flux has opposite sign of diffusive flux
      // now actually increment flux registers with entire flux
      Real fluxMult = m_dt; // petermc made negative, 7 Dec 07; then positive again because using edgeScal = -diffusiveFlux
      for (dit.reset(); dit.ok(); ++dit)
        {
          FluxBox& localFlux = edgeScal[dit];

          if (m_level > 0)
            {
              LevelFluxRegister& crseFR = *a_crseFluxRegPtr;
              CH_assert (crseFR.isDefined());
              // use same multipliers and sign convention as in computeUStar
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  crseFR.incrementFine(localFlux[dir],
                                       fluxMult, dit(), scalComps,
                                       scalComps, dir, Side::Lo);

                  crseFR.incrementFine(localFlux[dir],
                                       fluxMult, dit(), scalComps,
                                       scalComps, dir, Side::Hi);
                } // end loop over directions
            } // end if coarser level exists

          if (!finestLevel())
            {
              CH_assert(a_fineFluxReg.isDefined());
              for (int dir=0; dir<SpaceDim; ++dir)
                {
                  a_fineFluxReg.incrementCoarse(localFlux[dir],
                                                fluxMult, dit(),
                                                scalComps, scalComps,
                                                dir);
                } // end loop over directions
            } // end if finer level exists
        } // end loop over grids for flux register updates

      Interval thisComp(0, 0);
      int normType = 0;
      Real scalNorm = norm(a_old_scal, thisComp, normType);
      m_scalMGsolverPtrs[a_comp]->m_convergenceMetric = scalNorm;

      int numberMGlevels = (m_level == 0) ? 0 : 1;
      pout() << "calling updateSoln on scalar component " << a_comp
             << " on level " << m_level << endl;
      LevelFluxRegister* fineFluxRegPtr = NULL;
      if (!finestLevel())
        {
          CH_assert(a_fineFluxReg.isDefined());
          fineFluxRegPtr = &a_fineFluxReg;
        }
      m_TGAsolverScalPtrs[a_comp]->updateSoln(a_new_scal,
                                              a_old_scal,
                                              src,
                                              fineFluxRegPtr,
                                              a_crseFluxRegPtr,
                                              a_crse_scal_old,
                                              a_crse_scal_new,
                                              old_time,
                                              a_old_crse_time,
                                              a_new_crse_time, m_dt,
                                              numberMGlevels,
                                              false); // do not initialize to zero
      if (m_level > 0)
        {
          if (crseDataPtr != NULL)
            {
              delete crseDataPtr;
            }
        }
    } // end if is diffusive
}
