#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBConsBackwardEulerIntegrator.H"
#include "AMRTGA.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
EBConsBackwardEulerIntegrator::
EBConsBackwardEulerIntegrator(RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >& a_solver):
  m_solver(a_solver),
  m_op(0),
  m_eblg()
{
  Vector<MGLevelOp<LevelData<EBCellFAB> >*> ops = m_solver->getAllOperators();
  m_op = dynamic_cast<EBConductivityOp*>(ops[0]);
  CH_assert(m_op != 0);
  m_eblg = m_op->getEBLG();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBConsBackwardEulerIntegrator::
~EBConsBackwardEulerIntegrator()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBConsBackwardEulerIntegrator::
updateSolution(LevelData<EBCellFAB>& a_newTemperature,
               const LevelData<EBCellFAB>& a_oldTemperature,
               Real a_specificHeat,
               const LevelData<EBCellFAB>& a_oldDensity,
               const LevelData<EBCellFAB>& a_newDensity,
               const LevelData<EBCellFAB>& a_source,
               Real a_time,
               Real a_timeStep,
               bool a_zeroTemp)
{
  // Create new storage for the RHS.
  LevelData<EBCellFAB> RHS, oldRhoT;
  m_op->create(RHS, a_source);
  m_op->setToZero(RHS);
  m_op->create(oldRhoT, a_oldTemperature);
  m_op->setToZero(oldRhoT);

  // Assign (a_oldDensity * a_oldTemperature + dt * a_source/Cv) -> RHS.
  m_op->assign(RHS, a_source);
  m_op->assign(oldRhoT, a_oldDensity);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
  {
    RHS[dit()] *= a_timeStep / a_specificHeat;
    oldRhoT[dit()] *= a_oldTemperature[dit()];
  }
  m_op->incr(RHS, oldRhoT, 1.0);

  // Set alpha = 1 and beta = -dt/Cv. Make sure to set all of the multigrid operators, not just
  // the top level one.
  m_op->setAlphaAndBeta(1.0, -a_timeStep/a_specificHeat);
  Vector<MGLevelOp<LevelData<EBCellFAB> >*> solverOps = m_solver->getAllOperators();
  for (int iop = 0; iop < solverOps.size(); ++iop)
  {
    TGAHelmOp<LevelData<EBCellFAB> >* helmOp = dynamic_cast<TGAHelmOp<LevelData<EBCellFAB> >*>(solverOps[iop]);
    CH_assert(helmOp != 0);
    helmOp->setAlphaAndBeta(1.0, -a_timeStep/a_specificHeat);

    // Make sure that the operator is set to use the coefficients at the new time.
    helmOp->setTime(a_time, 1.0, a_timeStep);
  }

  // Solve [a_newDensity - dt * L / Cv] newT = RHS
  // for newT.
  if (a_zeroTemp)
    m_op->setToZero(a_newTemperature);
  Vector<LevelData<EBCellFAB>*> newTemps(1, &a_newTemperature),
                                RHSes(1, &RHS);
  m_solver->solve(newTemps, RHSes, 0, 0, a_zeroTemp);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBConsBackwardEulerIntegrator::
computeDiffusion(LevelData<EBCellFAB>& a_diffusionTerm,
                 const LevelData<EBCellFAB>& a_oldTemperature,
                 Real a_specificHeat,
                 const LevelData<EBCellFAB>& a_oldDensity,
                 const LevelData<EBCellFAB>& a_newDensity,
                 const LevelData<EBCellFAB>& a_source,
                 Real a_time,
                 Real a_timeStep,
                 bool a_zeroTemp)
{
  // Allocate storage for the new temperature.
  LevelData<EBCellFAB> newTemp;
  EBCellFactory fact(m_eblg.getEBISL());
  newTemp.define(m_eblg.getDBL(), 1, 4*IntVect::Unit, fact);
  EBLevelDataOps::setVal(newTemp, 0.0);

  //          n+1
  // Compute T.
  updateSolution(newTemp, a_oldTemperature, a_specificHeat, a_oldDensity,
                 a_newDensity, a_source, a_time, a_timeStep, a_zeroTemp);

  //            n+1             n
  // ([rho Cv T]    - [rho Cv T] )
  // ----------------------------- - a_source -> a_diffusionTerm.
  //              dt
  a_oldTemperature.copyTo(a_diffusionTerm);
  for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
  {
    EBCellFAB& diff = a_diffusionTerm[dit()];

    //           n
    // [rho Cv T]   -> diff
    diff *= a_oldDensity[dit()];
    diff *= a_specificHeat;

    //           n+1
    // [rho Cv T]   -> newRhoCvT
    EBCellFAB& newRhoCvT = newTemp[dit()];
    newRhoCvT *= a_newDensity[dit()];
    newRhoCvT *= a_specificHeat;

    // Subtract newRhoCvT from diff, and divide by -dt.
    diff -= newRhoCvT;
    diff /= -a_timeStep;

    // Subtract off the source.
    diff -= a_source[dit()];
  }
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
