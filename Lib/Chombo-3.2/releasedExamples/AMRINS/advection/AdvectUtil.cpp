#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectUtil.H"
#include "CellToEdge.H"
#include "AdvectionPhysics.H"

void
traceAdvectionVel(LevelData<FluxBox>& a_advVel,
                  LevelData<FArrayBox>& a_old_vel,
                  const LevelData<FArrayBox>& a_viscousSource,
                  PatchGodunov& a_patchGodVelocity,
                  Real a_old_time, Real a_dt)
{

  // Get AdvectionPhysics object within the PatchGodunov object
  AdvectionPhysics* advectionPhysics =
    dynamic_cast<AdvectionPhysics*>(a_patchGodVelocity.getGodunovPhysicsPtr());
  if (advectionPhysics == NULL)
    {
      MayDay::Error("AMRNavierStokes::computeAdvectionVelocities - unable to upcast GodunovPhysics to AdvectionPhysics");
    }

  a_patchGodVelocity.setCurrentTime(a_old_time);

  // this should now be unnecessary
  // also need to fill grown advection velocity
  // by cell-to-edge averaging old-time velocity
  CellToEdge(a_old_vel, a_advVel);

  // loop over grids and predict face-centered velocities at the half
  // time using the patchGodunov infrastructure
  const DisjointBoxLayout& levelGrids = a_advVel.getBoxes();
  DataIterator dit = a_advVel.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box gridBox = levelGrids[dit()];
      FArrayBox& thisOldVel = a_old_vel[dit()];
      FluxBox& thisAdvVel = a_advVel[dit()];
      const FArrayBox& srcFab = a_viscousSource[dit()];

      // temporary storage for the advected velocities
      FluxBox tempAdvVel(gridBox,SpaceDim);

      // set up PatchGodunov and AdvectionPhysics for this patch
      a_patchGodVelocity.setCurrentBox(gridBox);
      advectionPhysics->setCellVelPtr(&thisOldVel);
      advectionPhysics->setAdvVelPtr(&thisAdvVel);

      // predict face-centered velocities
      a_patchGodVelocity.computeWHalf(tempAdvVel, thisOldVel, srcFab,
                                      a_dt, gridBox);

      // copy the normal component of the advected velocity into "thisAdvVel"
      // (the rest is discarded)
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        int srcComp = idir;
        int destComp = 0;
        int numComp = 1;

        thisAdvVel[idir].copy(tempAdvVel[idir],srcComp,destComp,numComp);
      }
    } // end loop over grids

}


void
computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal,
                           LevelData<FluxBox>& a_edgeScalHalf,
                           LevelData<FArrayBox>& a_old_scal,
                           LevelData<FluxBox>& a_adv_vel,
                           LevelData<FArrayBox>& a_old_vel,
                           LevelData<FArrayBox>& a_diffusiveSrc,
                           PatchGodunov& a_patchGodScalar,
                           Real a_old_time, Real a_dt)
{

  int numScal = a_edgeScal.nComp();

  // Get AdvectionPhysics object within the PatchGodunov object
    AdvectionPhysics* advectionPhysics =
      dynamic_cast<AdvectionPhysics*>(a_patchGodScalar.getGodunovPhysicsPtr());
    if (advectionPhysics == NULL)
    {
      MayDay::Error("AMRNavierStokes::advectDiffuseScalar - unable to upcast GodunovPhysics to AdvectionPhysics");
    }

    // set up patchGodunov for this problem
    a_patchGodScalar.setCurrentTime(a_old_time);

    // also need to build a grown adveaction velocity
    IntVect advectVelGrow(IntVect::Unit);
    const DisjointBoxLayout& levelGrids = a_old_vel.getBoxes();
    LevelData<FluxBox> grownAdvVel(levelGrids, 1, advectVelGrow);
    // this is inefficient, but we'll see how bad it really is
    // for now, fill extra cells with cell-to-face averaged velocities
    CellToEdge(a_old_vel, grownAdvVel);
    // now overwrite with advection velocities wherever possible
    a_adv_vel.copyTo(a_adv_vel.interval(), grownAdvVel,
                     grownAdvVel.interval());


    // now trace scalars to edges at time n+1/2
    DataIterator dit = levelGrids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
      {
        FluxBox& thisEdgeScal = a_edgeScal[dit()];
        FluxBox& thisEdgeScalHalf = a_edgeScalHalf[dit()];
        FluxBox& thisAdvVel = grownAdvVel[dit()];
        FArrayBox& thisCellVel = a_old_vel[dit()];
        FArrayBox& thisOldScal = a_old_scal[dit()];
        FArrayBox& thisSrc = a_diffusiveSrc[dit()];

        // PatchGodunov-based approach to advection
        const Box& curBox = levelGrids[dit()];
        a_patchGodScalar.setCurrentBox(curBox);
        // set cell-centered velocity field
        advectionPhysics->setCellVelPtr(&thisCellVel);
        // set advection velocity field
        advectionPhysics->setAdvVelPtr(&thisAdvVel);

        // compute face-centered, predicted scalars
        a_patchGodScalar.computeWHalf(thisEdgeScal, thisOldScal,
                                      thisSrc, a_dt, curBox);

        for (int dir=0; dir<SpaceDim; dir++)
          {
            (thisEdgeScalHalf[dir]).copy(thisEdgeScal[dir]) ;

            // multiply by edge velocity to get flux
            // do componentwise
            for (int comp=0; comp<numScal; comp++)
              {
                thisEdgeScal[dir].mult(thisAdvVel[dir],0,comp,1);
              }
          } // end loop over tracing directions
      } // end loop over grids for tracing

}


void
computeScalarAdvectiveFlux(LevelData<FluxBox>& a_edgeScal,
                           LevelData<FArrayBox>& a_old_scal,
                           LevelData<FluxBox>& a_adv_vel,
                           LevelData<FArrayBox>& a_old_vel,
                           LevelData<FArrayBox>& a_diffusiveSrc,
                           PatchGodunov& a_patchGodScalar,
                           Real a_old_time, Real a_dt)
{

  int numScal = a_edgeScal.nComp();

  // Get AdvectionPhysics object within the PatchGodunov object
    AdvectionPhysics* advectionPhysics =
      dynamic_cast<AdvectionPhysics*>(a_patchGodScalar.getGodunovPhysicsPtr());
    if (advectionPhysics == NULL)
    {
      MayDay::Error("AMRNavierStokes::advectDiffuseScalar - unable to upcast GodunovPhysics to AdvectionPhysics");
    }

    // set up patchGodunov for this problem
    a_patchGodScalar.setCurrentTime(a_old_time);

    // also need to build a grown adveaction velocity
    IntVect advectVelGrow(IntVect::Unit);
    const DisjointBoxLayout& levelGrids = a_old_vel.getBoxes();
    LevelData<FluxBox> grownAdvVel(levelGrids, 1, advectVelGrow);
    // this is inefficient, but we'll see how bad it really is
    // for now, fill extra cells with cell-to-face averaged velocities
    CellToEdge(a_old_vel, grownAdvVel);
    // now overwrite with advection velocities wherever possible
    a_adv_vel.copyTo(a_adv_vel.interval(), grownAdvVel,
                     grownAdvVel.interval());


    // now trace scalars to edges at time n+1/2
    DataIterator dit = levelGrids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
      {
        FluxBox& thisEdgeScal = a_edgeScal[dit()];
        FluxBox& thisAdvVel = grownAdvVel[dit()];
        FArrayBox& thisCellVel = a_old_vel[dit()];
        FArrayBox& thisOldScal = a_old_scal[dit()];
        FArrayBox& thisSrc = a_diffusiveSrc[dit()];

        // PatchGodunov-based approach to advection
        const Box& curBox = levelGrids[dit()];
        a_patchGodScalar.setCurrentBox(curBox);
        // set cell-centered velocity field
        advectionPhysics->setCellVelPtr(&thisCellVel);
        // set advection velocity field
        advectionPhysics->setAdvVelPtr(&thisAdvVel);

        // compute face-centered, predicted scalars
        a_patchGodScalar.computeWHalf(thisEdgeScal, thisOldScal,
                                      thisSrc, a_dt, curBox);

        for (int dir=0; dir<SpaceDim; dir++)
          {
            // multiply by edge velocity to get flux
            // do componentwise
            for (int comp=0; comp<numScal; comp++)
              {
                thisEdgeScal[dir].mult(thisAdvVel[dir],0,comp,1);
              }
          } // end loop over tracing directions
      } // end loop over grids for tracing

}


void
computePredictedVelocities(LevelData<FluxBox>& a_uHalf,
                           LevelData<FArrayBox>& a_traceVel,
                           LevelData<FluxBox>& a_advVel,
                           LevelData<FArrayBox>& a_viscousSource,
                           PatchGodunov& a_patchGodVelocity,
                           LevelData<FluxBox>& a_grad_eLambda,
                           LevelData<FluxBox>& a_gradPhi,
                           bool a_applyFreestreamCorrection,
                           Real a_old_time, Real a_dt)
{

  // Get AdvectionPhysics object within the PatchGodunov object
  AdvectionPhysics* advectionPhysics =
    dynamic_cast<AdvectionPhysics*>(a_patchGodVelocity.getGodunovPhysicsPtr());
  if (advectionPhysics == NULL)
    {
      MayDay::Error("AMRNavierStokes::predictVelocities - unable to upcast GodunovPhysics to AdvectionPhysics");
    }

  a_patchGodVelocity.setCurrentTime(a_old_time);

  const DisjointBoxLayout& levelGrids = a_traceVel.getBoxes();
  // loop over grids and do prediction
    DataIterator dit = levelGrids.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        Box gridBox = levelGrids[dit()];
        FArrayBox& thisOldVel = a_traceVel[dit()];
        FluxBox& thisAdvVel = a_advVel[dit()];
        FluxBox& thisPredictedVel = a_uHalf[dit()];
        FArrayBox& srcFab = a_viscousSource[dit()];

        // set up patchGodunov for this patch
        a_patchGodVelocity.setCurrentBox(gridBox);
        advectionPhysics->setCellVelPtr(&thisOldVel);
        advectionPhysics->setAdvVelPtr(&thisAdvVel);

        // compute face-centered variables using "Godunov Box"
        a_patchGodVelocity.computeWHalf(thisPredictedVel, thisOldVel,
                                        srcFab, a_dt, gridBox);

        // now loop over directions -- normal direction velocities
        // are copied from advection velocities; for tangential direction
        // we use face-centered predicted velocities (flux divided by
        // advection velocity)
        for (int velComp = 0; velComp<SpaceDim; velComp++)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& thisVelDir = thisPredictedVel[dir];
                FArrayBox& thisAdvVelDir = thisAdvVel[dir];

                if (dir == velComp)
                  {
                    // normal direction -- copy from advVel->uHalf
                    // srcComp is always 0, since advVel only has
                    // one component
                    int srcComp=0;
                    int destComp = velComp;

                    thisVelDir.copy(thisAdvVelDir, srcComp, destComp, 1);

                    if (a_applyFreestreamCorrection)
                      {
                        // now need to subtract off grad(eLambda)
                        // this is the due to the difference between
                        // _advecting_ and _advected_ velocities
                        FluxBox& thisCorr = a_grad_eLambda[dit()];
                        FArrayBox& thisCorrDir = thisCorr[dir];
                        thisVelDir.minus(thisCorrDir, 0, dir, 1);
                      }
                  } // end if dir is velcomp
                else
                  {
                    // need to add MAC correction to traced velocities
                    thisVelDir.minus(a_gradPhi[dit()][dir],
                                     velComp, velComp, 1);
                  } // end if tangential direction
              } // end loop over face directions
          } // end loop over velocity components
      } // end loop over grids
}
