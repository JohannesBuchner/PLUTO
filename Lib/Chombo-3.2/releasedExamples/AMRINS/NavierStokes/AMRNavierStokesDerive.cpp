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
#include "AMRNSF_F.H"

// ---------------------------------------------------------
void
AMRNavierStokes::computeVorticity(LevelData<FArrayBox>& a_vorticity) const
{
  // this breaks the "const"-ness of the function,
  // but is necessary to ensure that boundary
  // conditions are properly set.
  LevelData<FArrayBox>& vel =
    *(const_cast<LevelData<FArrayBox>*>(m_vel_new_ptr));

  const DisjointBoxLayout& grids = vel.getBoxes();

  if (m_level > 0)
    {
      // do quadratic C/F BC's
      // for now, assume that BC's are with new velocity
      // (may need to be changed for fine-level regridding)
      LevelData<FArrayBox>& crseVel = crseNSPtr()->newVel();
      const DisjointBoxLayout& crseGrids = crseVel.getBoxes();
      int nRefCrse = crseNSPtr()->refRatio();

      QuadCFInterp interpolator(grids, &crseGrids, m_dx, nRefCrse,
                                SpaceDim, m_problem_domain);

      interpolator.coarseFineInterp(vel, crseVel);
    }

  Interval velComps(0,SpaceDim-1);
  vel.exchange(velComps);

  DataIterator dit = vel.dataIterator();

  // at the moment, do extrap at physical boundaries
  // (same effect as one-sided differencing)
  BCHolder vortVelBC = m_physBCPtr->extrapFuncBC(1);

  for (dit.reset(); dit.ok(); ++dit)
    {
      // apply physical BC's
      vortVelBC(vel[dit], grids[dit],
                m_problem_domain, m_dx,
                true); // homogeneous
      if (SpaceDim == 3)
        {
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // and then compute vorticity
              FORT_COMPUTEVORT(CHF_FRA1(a_vorticity[dit],dir),
                               CHF_CONST_FRA(vel[dit]),
                               CHF_BOX(grids[dit]),
                               CHF_CONST_REAL(m_dx),
                               CHF_CONST_INT(dir));
            }
        }
      else if (SpaceDim == 2)
        {
          int dir = 2;
          FORT_COMPUTEVORT(CHF_FRA1(a_vorticity[dit],0),
                           CHF_CONST_FRA(vel[dit]),
                           CHF_BOX(grids[dit]),
                           CHF_CONST_REAL(m_dx),
                           CHF_CONST_INT(dir));
        } // end dimensionality choice
    } // end loop over grids
}

// ---------------------------------------------------------
void
AMRNavierStokes::computeLapVel(LevelData<FArrayBox>& a_lapVel,
                               LevelData<FArrayBox>& a_vel,
                               const LevelData<FArrayBox>* a_crseVelPtr)
{
  // set BC's
  VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());
  bool isHomogeneous = false;
  m_velocityAMRPoissonOp.applyOp(a_lapVel, a_vel, a_crseVelPtr,
                                 isHomogeneous, velBC);

  // may need to extend lapVel to cover ghost cells as well
  {
    BCHolder viscBC = m_physBCPtr->viscousFuncBC();
    const DisjointBoxLayout& grids = a_lapVel.getBoxes();
    DataIterator dit = a_lapVel.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
      {
        viscBC(a_lapVel[dit], grids[dit],
               m_problem_domain, m_dx,
               false); // not homogeneous
      }
  }

  // finally, do exchange
  a_lapVel.exchange(a_lapVel.interval());
}

// ---------------------------------------------------------
void
AMRNavierStokes::computeLapScal(LevelData<FArrayBox>& a_lapScal,
                                LevelData<FArrayBox>& a_scal,
                                const BCHolder& a_physBC,
                                const LevelData<FArrayBox>* a_crseScalPtr)
{
  m_scalarsAMRPoissonOp.setBC(a_physBC);
  bool isHomogeneous = false;
  if (a_crseScalPtr != NULL)
    {
      m_scalarsAMRPoissonOp.AMROperatorNF(a_lapScal, a_scal, *a_crseScalPtr,
                                          isHomogeneous);
    }
  else
    {
      m_scalarsAMRPoissonOp.applyOpI(a_lapScal, a_scal,
                                     isHomogeneous);
    }

  BCHolder viscBC = m_physBCPtr->viscousFuncBC();

  DataIterator dit = a_lapScal.dataIterator();
  const DisjointBoxLayout& grids = a_lapScal.getBoxes();
  for (dit.reset(); dit.ok(); ++dit)
    {
      viscBC(a_lapScal[dit],
             grids[dit],
             m_problem_domain, m_dx,
             false); // not homogeneous
    }

  // finally, do exchange
  a_lapScal.exchange(a_lapScal.interval());
}

// ---------------------------------------------------------
void
AMRNavierStokes::computeKineticEnergy(LevelData<FArrayBox>& a_energy) const
{
  // this is simple since no BC's need to be set.
  const LevelData<FArrayBox>& levelVel = *m_vel_new_ptr;
  const DisjointBoxLayout& levelGrids = levelVel.getBoxes();

  DataIterator dit = a_energy.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FORT_KINETICENERGY(CHF_FRA1(a_energy[dit],0),
                         CHF_CONST_FRA(levelVel[dit]),
                         CHF_BOX(levelGrids[dit]));
    }
}
