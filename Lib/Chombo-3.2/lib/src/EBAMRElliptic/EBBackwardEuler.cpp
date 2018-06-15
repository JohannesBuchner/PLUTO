#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AMRMultiGrid.H"
#include "EBBackwardEuler.H"
#include "EBLevelDataOps.H"
#include "NamespaceHeader.H"

void EBBackwardEuler::
resetAlphaAndBeta(const Real& a_alpha,
                  const Real& a_beta)
{
  CH_TIME("EBBackwardEuler::resetAlphaAndBeta");
  Vector<MGLevelOp< LevelData<EBCellFAB> >* > ops = m_solver->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
    {

      TGAHelmOp< LevelData<EBCellFAB> >* helmop = (TGAHelmOp< LevelData<EBCellFAB> >*) ops[iop];
      helmop->setAlphaAndBeta(a_alpha, a_beta);
    }
}

/*****/
TGAHelmOp<LevelData<EBCellFAB> >*
EBBackwardEuler::
newOp(const ProblemDomain&                             a_indexSpace,
      const AMRLevelOpFactory<LevelData<EBCellFAB> >&  a_opFact)
{
  AMRLevelOpFactory<LevelData<EBCellFAB> >& opFact = (AMRLevelOpFactory<LevelData<EBCellFAB> >&) a_opFact;
  TGAHelmOp<LevelData<EBCellFAB> >* retval = (TGAHelmOp<LevelData<EBCellFAB> >*) opFact.AMRnewOp(a_indexSpace);
  return retval;
}
/*****/
EBBackwardEuler::
EBBackwardEuler(const RefCountedPtr<AMRMultiGrid< LevelData<EBCellFAB> > >& a_solver,
                const AMRLevelOpFactory<LevelData<EBCellFAB> > &            a_opFact,
                const ProblemDomain&                                        a_level0Domain,
                const Vector<int>&                                          a_refRat,
                int a_numLevels, int a_verbosity)
{
  CH_TIME("EBBackwardEuler::EBBackwardEuler");
  m_verbosity = a_verbosity;
  m_level0Domain = a_level0Domain;
  m_refRat = a_refRat;
  m_solver  = a_solver;
  m_numLevels = a_numLevels;
  if (m_numLevels < 0)
    {
      m_numLevels = a_refRat.size();
    }

  m_ops.resize(m_numLevels);
  Vector< AMRLevelOp< LevelData<EBCellFAB> > * >& amrops =  m_solver->getAMROperators();
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      m_ops[ilev] = dynamic_cast<TGAHelmOp< LevelData< EBCellFAB > >* >(amrops[ilev]);
      if (m_ops[ilev]==NULL)
        {
          MayDay::Error("dynamic cast failed---is that operator really a TGAHelmOp?");
        }
    }

  m_dataCreated = false;
}
/*****/
void
EBBackwardEuler::
computeDiffusion(Vector< LevelData<EBCellFAB>* >&       a_diffusiveTerm,
                 Vector< LevelData<EBCellFAB>* >&       a_phiOld,
                 Vector< LevelData<EBCellFAB>* >&       a_src,
                 Real a_oldTime,
                 Real a_dt,
                 int a_lbase, int a_lmax, bool a_zeroPhi)
{
  CH_TIME("EBBackwardEuler::computeDiffusion");
  // first compute updated solution
  Vector<LevelData<EBCellFAB>* > tempSoln(a_lmax+1, NULL);

  int minLev = Max(a_lbase-1, 0);

  for (int ilev = minLev; ilev <= a_lmax; ilev++)
    {
      tempSoln[ilev] = new LevelData<EBCellFAB>();
    }
  for (int ilev = minLev; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->create(*a_diffusiveTerm[ilev], *a_phiOld[ilev]);
      m_ops[ilev]->create(       *tempSoln[ilev], *a_phiOld[ilev]);
      m_ops[ilev]->setToZero(    *tempSoln[ilev]);
      if (!a_zeroPhi)
        {
          m_ops[ilev]->assign(*tempSoln[ilev], *a_phiOld[ilev]);
        }
    }
  oneStep(tempSoln, a_phiOld, a_src,
          a_dt, a_lbase, a_lmax, a_zeroPhi);

  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      // now subtract everything off to leave us with diffusive term
      m_ops[ilev]->incr( *tempSoln[ilev], *a_phiOld[ilev], -1.0);
      m_ops[ilev]->scale(*tempSoln[ilev], 1.0/a_dt);

      //now multiply by a if there is an a
      m_ops[ilev]->diagonalScale(*tempSoln[ilev]);

      // and finally, subtract off a_src
      m_ops[ilev]->incr(*tempSoln[ilev], *a_src[ilev], -1.0);

      // what's left should be the time-centered diffusive part of the update
      m_ops[ilev]->assign(*a_diffusiveTerm[ilev], *tempSoln[ilev]);
    }
  for (int ilev = minLev; ilev <= a_lmax; ilev++)
    {
      delete tempSoln[ilev];
    }
}
/*****/
EBBackwardEuler::~EBBackwardEuler()
{
  CH_TIME("EBBackwardEuler::destructor");
  for (int ilev = 0; ilev < m_rhst.size(); ilev++)
    {
      if (m_rhst[ilev] != NULL)
        {
          delete m_rhst[ilev];
          m_rhst[ilev] = NULL;
        }
    }
}
/***/
void
EBBackwardEuler::
createData(Vector<LevelData<EBCellFAB>* >&       a_source,
           int a_lbase,int a_lmax)
{
  CH_TIME("EBBackwardEuler::createData");
  m_dataCreated = true;
  m_rhst.resize(a_source.size(), NULL);
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_rhst[ilev] = new LevelData<EBCellFAB>();
      m_ops[ilev]->create(*m_rhst[ilev], *a_source[ilev]);
    }
}
/***/
void EBBackwardEuler::
oneStep(Vector<LevelData<EBCellFAB>* >&       a_phiNew,
        Vector<LevelData<EBCellFAB>* >&       a_phiOld,
        Vector<LevelData<EBCellFAB>* >&       a_source,
        const Real&        a_dt,
        int                a_lbase,
        int                a_lmax,
        bool               a_zeroPhi,
        const bool         a_kappaWeighted)
{
  if (!m_dataCreated)
    {
      createData(a_source, a_lbase, a_lmax);
    }

  if (m_verbosity >= 3)
    {
      pout() << "  EBBackwardEuler:: making rhs" << std::endl;
    }
  createEulerRHS(m_rhst, a_source, a_phiOld, a_lbase, a_lmax, a_dt, a_kappaWeighted);

  if (m_verbosity >= 3)
    {
      pout() << "  EBBackwardEuler:: solving" << std::endl;
    }
  //have to set one level below lbase in case we are using
  //a coarser level for boundary conditions
  int minLev = Max(a_lbase-1, 0);
  for (int ilev = minLev; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->assign( *a_phiNew[ilev], *a_phiOld[ilev]);
    }

  solveHelm(a_phiNew, m_rhst, a_lbase, a_lmax,a_dt, a_zeroPhi);
}
/*******/
//fills a_ans = dt*kappa*a_source + kappa*acoef*phiOld
void EBBackwardEuler::
createEulerRHS(Vector<LevelData<EBCellFAB>* >&   a_ans,
               Vector<LevelData<EBCellFAB>* >&   a_rho,
               Vector<LevelData<EBCellFAB>* >&   a_phiOld,
               int                               a_lbase,
               int                               a_lmax,
               Real                              a_dt,
               const bool                        a_kappaWeighted)

{
  CH_TIME("EBBackwardEuler::createEulerRHS");
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      DisjointBoxLayout grids = a_phiOld[ilev]->disjointBoxLayout();
      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          //this makes rhs  = 0
          (*a_ans[ilev])[dit()].setVal(0.);
          //this makes rhs = phiOld
          (*a_ans[ilev])[dit()] +=  ((*a_phiOld[ilev])[dit()]);
        }

      // //this makes rhs = kappa*acoef*(phi^n)
      // //do not kappa weight if already coming in kappa weighted
      m_ops[ilev]->diagonalScale(*a_ans[ilev], !a_kappaWeighted);

      for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
        {
          const EBCellFAB& rho  = (*a_rho[ilev])[dit()];
          EBCellFAB scaleRho(rho.getEBISBox(), rho.box(), rho.nComp());
          scaleRho.setVal(0.);

          scaleRho += rho;
          //this makes scalerho = kappa*a_rho
          if (!a_kappaWeighted) EBLevelDataOps::kappaWeight(scaleRho);

          //this makes scalerho = dt*kappa*a_rho
          scaleRho *= a_dt;

          //this makes rhs = kappa*acoef*(phi^n) + dt*kappa*a_rho
          (*a_ans[ilev])[dit()] += scaleRho;
        }
    }
}
/*******/
void EBBackwardEuler::
solveHelm(Vector<LevelData<EBCellFAB>*>&       a_ans,
          Vector<LevelData<EBCellFAB>*>&       a_rhs,
          int               a_lbase,
          int               a_lmax,
          Real              a_dt,
          bool              a_zeroPhi)
{
  CH_TIME("EBBackwardEuler::solveHelm");
  Real factor  = -a_dt;
  resetAlphaAndBeta(1.0, factor);

  m_solver->solve(a_ans, a_rhs, a_lmax, a_lbase, a_zeroPhi);

  if (m_solver->m_exitStatus==1 && m_verbosity>3)
    {
      pout() << "EBBackwardEuler:: WARNING: solver exitStatus == 1" << std::endl;
    }
}

#include "NamespaceFooter.H"

