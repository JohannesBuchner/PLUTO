#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBCrankNicholsonStrategy.H"
#include "EBLevelDataOps.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
EBCrankNicholsonStrategy::
EBCrankNicholsonStrategy(Real a_implicitness):
  EBImplicitIntegrationStrategy(),
  m_implicitness(a_implicitness)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBCrankNicholsonStrategy::
~EBCrankNicholsonStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBCrankNicholsonStrategy::
updateSolution(LevelData<EBCellFAB>& a_phiNew,
               LevelData<EBCellFAB>& a_phiOld,
               LevelData<EBCellFAB>& a_src,
               Real a_oldTime,
               Real a_dt)
{
  int ncomp = a_phiNew.nComp();
  Interval intervBase(0, ncomp-1);
  HelmholtzOpType& op = this->helmholtzOp();
  if (!op.isTimeDependent())
  {
    // Scale the source term by kappa (k).
    LevelData<EBCellFAB> kSrc;
    op.create(kSrc, a_src);
    op.assign(kSrc, a_src);
    EBLevelDataOps::kappaWeight(kSrc);

    //                                                        n
    // Assign [k * a + (1 - implicitness) * dt * k * L ] * phi   + dt * k * src -> rhs.
    LevelData<EBCellFAB> rhs;
    op.create(rhs, a_phiOld);
    applyHelmholtz(rhs, a_phiOld, 1.0 - m_implicitness, a_dt, true, a_oldTime);
    EBLevelDataOps::axby(rhs, rhs, kSrc, 1.0, a_dt);

    //                                                 n+1              n+1
    // Solve [k * a - implicitness * dt * k * L   ] phi   = rhs  for phi.
    solveHelmholtz(a_phiNew, rhs, m_implicitness, a_dt, a_oldTime + a_dt);
  }
  else
  {
    // Beware! Time-dependent operator.

    // Scale the source term by kappa (k).
    LevelData<EBCellFAB> kSrc;
    op.create(kSrc, a_src);
    op.assign(kSrc, a_src);
    EBLevelDataOps::kappaWeight(kSrc);

    //              n                                  n       n
    // Assign [k * a  + (1 - implicitness) * dt * k * L ] * phi   + dt * k * src -> rhs.
    LevelData<EBCellFAB> rhs;
    op.create(rhs, a_phiOld);
    op.setTime(a_oldTime, 0.0, a_dt);
    applyHelmholtz(rhs, a_phiOld, 1.0 - m_implicitness, a_dt, true, a_oldTime);
    EBLevelDataOps::axby(rhs, rhs, kSrc, 1.0, a_dt);

    //             n+1                            n+1     n+1              n+1
    // Solve [k * a    - implicitness * dt * k * L   ] phi   = rhs  for phi.
    op.setTime(a_oldTime, 1.0, a_dt);
    solveHelmholtz(a_phiNew, rhs, m_implicitness, a_dt, a_oldTime + a_dt);
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBCrankNicholsonStrategy::
computeDiffusion(LevelData<EBCellFAB>& a_diffusiveTerm,
                 LevelData<EBCellFAB>& a_phiOld,
                 LevelData<EBCellFAB>& a_src,
                 Real a_oldTime,
                 Real a_dt)
{
  HelmholtzOpType& op = this->helmholtzOp();
  if (!op.isTimeDependent())
  {
    // The operator has no time-dependent parameters. Life is easier.

    // First compute the updated solution.
    LevelData<EBCellFAB> phiNew;
    op.create(phiNew, a_phiOld);
    op.setToZero(phiNew);
    updateSolution(phiNew, a_phiOld, a_src, a_oldTime, a_dt);

    // Now compute the diffusive term.
    EBLevelDataOps::axby(a_diffusiveTerm, phiNew, a_phiOld, 1.0/a_dt, -1.0/a_dt);

    // Scale the term by the a coefficient (no kappa weighting).
    op.diagonalScale(a_diffusiveTerm, false);

    // Finally, subtract off the source.
    for (DataIterator dit = phiNew.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      a_diffusiveTerm[dit()] -= a_src[dit()];
  }
  else
  {
    // The operator has time-dependent coefficients. We must be more careful!

    // Compute the a coefficient at time n, n+1/2, and n+1.
    LevelData<EBCellFAB> aOld, aHalf, aNew;
    op.create(aOld, a_phiOld);
    op.create(aHalf, a_phiOld);
    op.create(aNew, a_phiOld);
    for (DataIterator dit = a_phiOld.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
    {
      aOld[dit()].setVal(1.);
      aNew[dit()].setVal(1.);
      aHalf[dit()].setVal(1.);
    }
    op.setTime(a_oldTime, 0.0, a_dt);
    op.diagonalScale(aOld, false); // No kappa weighting!
    op.setTime(a_oldTime, 1.0, a_dt);
    op.diagonalScale(aNew, false); // No kappa weighting!
    op.setTime(a_oldTime, 0.5, a_dt);
    op.diagonalScale(aHalf, false); // No kappa weighting!

    // Compute the updated solution for phi.
    LevelData<EBCellFAB> phiNew;
    op.create(phiNew, a_phiOld);
    op.setToZero(phiNew);
    updateSolution(phiNew, a_phiOld, a_src, a_oldTime, a_dt);

    //              n+1     n                 n+1   n
    //  n+1/2   (phi   - phi )      n+1/2   (a   - a )
    // a      * -------------- + phi      * ---------- - a_src  -> a_diffusiveTerm.
    //                dt                       dt
    LevelData<EBCellFAB> phidadt, dadt;
    op.create(phidadt, a_phiOld);
    op.create(dadt, a_phiOld);
    for (DataIterator dit = a_phiOld.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
    {
      // Evaluate the first term above.
      a_diffusiveTerm[dit()].axby(phiNew[dit()], a_phiOld[dit()], 1.0/a_dt, -1.0/a_dt);
      a_diffusiveTerm[dit()] *= aHalf[dit()];

      // Compute a second order estimate of phidadt.
      dadt[dit()].axby(aNew[dit()], aOld[dit()], 1.0/a_dt, -1.0/a_dt);
      phidadt[dit()].axby(phiNew[dit()], a_phiOld[dit()], 0.5, 0.5);
      phidadt[dit()] *= dadt[dit()];

      // Add phidadt and subtract the original source term.
      a_diffusiveTerm[dit()] += phidadt[dit()];
      a_diffusiveTerm[dit()] -= a_src[dit()];
    }
  }
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
