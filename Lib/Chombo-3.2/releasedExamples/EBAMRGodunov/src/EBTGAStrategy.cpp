#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBTGAStrategy.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
EBTGAStrategy::
EBTGAStrategy():
  EBImplicitIntegrationStrategy(),
  m_mu1(0.0), m_mu2(0.0), m_mu3(0.0), m_mu4(0.0), m_r1(0.0)
{
  Real tgaEpsilon = 1.e-12;
  Real a = 2.0 - std::sqrt(2.0) - tgaEpsilon;
  m_mu1 = (a - std::sqrt(pow(a,2) - 4.0*a + 2.0))/2.0 ;
  m_mu2 = (a + std::sqrt(pow(a,2) - 4.0*a + 2.0))/2.0 ;
  m_mu3 = (1.0 - a);
  m_mu4 = 0.5 - a;

  Real discr = std::sqrt(a*a - 4.0*a + 2.0);
  m_r1 = (2.0*a - 1.0)/(a + discr);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBTGAStrategy::
~EBTGAStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBTGAStrategy::
updateSolution(LevelData<EBCellFAB>& a_phiNew,
               LevelData<EBCellFAB>& a_phiOld,
               LevelData<EBCellFAB>& a_src,
               Real a_oldTime,
               Real a_dt)
{
  int ncomp = a_phiNew.nComp();
  Interval intervBase(0, ncomp-1);
  HelmholtzOpType& op = this->helmholtzOp();

  LevelDataType rhst, srct, phis;
  op.create(rhst, a_src);
  op.create(srct, a_phiNew);
  op.create(phis, a_phiNew);

  op.setToZero(srct);
  op.setToZero(rhst);
  op.incr(srct, a_src, 1.0);

  // Divide the source S by the identity coefficient a and scale it by dt.
  op.divideByIdentityCoef(srct);
  op.scale(srct, a_dt);

  // from here on k is kappa and L is kappa L
  // this makes rhs hold dt * (k*a I + mu4 L) (S/a)
  // 'false', ie, use the extrap data
  applyHelmholtz(rhst, srct, m_mu4, a_dt, false, a_oldTime);

  // This makes a_phiNew hold (k*a I + mu3 L) phi^n
  // true -> apply CF and domain BC
  applyHelmholtz(a_phiNew, a_phiOld, m_mu3, a_dt, true, a_oldTime);

  // This makes rhs hold (k*a I + mu3 L) phi^n + dt(k*a I +mu4 L) S/a
  op.incr(rhst, a_phiNew, 1.0);

  // This makes phinew = (k*a I - mu2 L)^-1 (rhs) at the intermediate time.
  solveHelmholtz(a_phiNew, rhst, m_mu2, a_dt, a_oldTime + (1.0 - m_r1)*a_dt);

  // This puts the answer into rhst so we can do the final solve.
  op.assign(rhst, a_phiNew);

  // This makes rhs hold k*a [(k*a I - mu2 L)^-1 (rhs)]
  op.diagonalScale(rhst);

  // This makes phinew = (k*a I - mu1 L)^-1 [ka ((k*a I - mu2 L)^-1 (orig rhs))]
  solveHelmholtz(a_phiNew, rhst, m_mu1, a_dt, a_oldTime + a_dt);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBTGAStrategy::
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

    // first compute updated solution
    LevelData<EBCellFAB> phiNew;
    op.create(phiNew, a_phiOld);
    op.setToZero(phiNew);
    updateSolution(phiNew, a_phiOld, a_src, a_oldTime, a_dt);

    // Now compute the diffusive term.
    for (DataIterator dit = a_diffusiveTerm.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      a_diffusiveTerm[dit()].axby(phiNew[dit()], a_phiOld[dit()], 1.0/a_dt, -1.0/a_dt);

    // Scale the term by the a coefficient (no kappa weighting).
    op.diagonalScale(a_diffusiveTerm, false);

    // Finally, subtract off the source.
    for (DataIterator dit = a_diffusiveTerm.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
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

    // (We leave the operator at the half-step time for the remainder of
    //  this calculation.)

    // Subtract a first-order estimate of phi * da/dt from the source.
    LevelData<EBCellFAB> phidadt, dadt, rhs;
    op.create(phidadt, a_phiOld);
    op.create(dadt, a_phiOld);
    op.create(rhs, a_phiOld);
    op.assign(phidadt, a_phiOld);
    for (DataIterator dit = a_phiOld.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
    {
      dadt[dit()].axby(aNew[dit()], aOld[dit()], 1.0/a_dt, -1.0/a_dt);
      phidadt[dit()] *= dadt[dit()];
      rhs[dit()].axby(a_src[dit()], phidadt[dit()], 1.0, -1.0);
    }

    // Compute the updated solution.
    LevelData<EBCellFAB> phiNew;
    op.create(phiNew, a_phiOld);
    op.setToZero(phiNew);
    updateSolution(phiNew, a_phiOld, rhs, a_oldTime, a_dt);

    //              n+1     n                 n+1   n
    //  n+1/2   (phi   - phi )      n+1/2   (a   - a )
    // a      * -------------- - phi      * ---------- - a_src  -> a_diffusiveTerm.
    //                dt                       dt
    for (DataIterator dit = a_phiOld.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
    {
      // Evaluate the first term above.
      a_diffusiveTerm[dit()].axby(phiNew[dit()], a_phiOld[dit()], 1.0/a_dt, -1.0/a_dt);
      a_diffusiveTerm[dit()] *= aHalf[dit()];

      // Correct our estimate of phidadt to second order.
      phidadt[dit()].axby(phiNew[dit()], a_phiOld[dit()], 0.5, 0.5);
      phidadt[dit()] *= dadt[dit()];

      // Subtract phidadt and the original source term.
      a_diffusiveTerm[dit()] -= phidadt[dit()];
      a_diffusiveTerm[dit()] -= a_src[dit()];
    }
  }
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
