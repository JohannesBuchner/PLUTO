#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iomanip>

#include "parstream.H"

#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "AMRIO.H"
#include "computeSum.H"
#include "computeNorm.H"

#include "AMRLevelAdvectDiffuse.H"
#include "AdvectPhysicsF_F.H"
#include "AdvectTestIBC.H"
#include "AdvectTestF_F.H"
#include "GodunovUtilitiesF_F.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

RefCountedPtr<LevelTGA>                                   AMRLevelAdvectDiffuse::s_diffuseLevTGA = RefCountedPtr<LevelTGA>();
RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >       AMRLevelAdvectDiffuse::s_diffuseAMRMG  = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >();
RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >  AMRLevelAdvectDiffuse::s_diffuseOpFact = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >();
BiCGStabSolver<LevelData< FArrayBox> >                    AMRLevelAdvectDiffuse::s_botSolver;

/*******/
AMRLevelAdvectDiffuse::
~AMRLevelAdvectDiffuse()
{
  s_diffuseLevTGA  = RefCountedPtr<LevelTGA>();
  s_diffuseAMRMG   = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >();
  s_diffuseOpFact  = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >();
}

/*******/
void
AMRLevelAdvectDiffuse::
getHierarchyAndGrids(Vector<AMRLevelAdvectDiffuse*>&        a_hierarchy,
                     Vector<DisjointBoxLayout>&             a_grids,
                     Vector<int>&                           a_refRat,
                     ProblemDomain&                         a_lev0Dom,
                     Real&                                  a_lev0Dx)
{
  Vector<AMRLevel*> hierarchy = AMRLevel::getAMRLevelHierarchy();
  int nlevels = hierarchy.size();

  a_hierarchy.resize(nlevels);
  a_refRat.resize(   nlevels);
  a_grids.resize(    nlevels);

  AMRLevelAdvectDiffuse* coarsestLevel = (AMRLevelAdvectDiffuse*)(hierarchy[0]);
  a_lev0Dx       = coarsestLevel->m_dx;
  a_lev0Dom      = coarsestLevel->m_problem_domain;

  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      AMRLevelAdvectDiffuse* adLevel = (AMRLevelAdvectDiffuse*)(hierarchy[ilev]);

      a_hierarchy[ilev] = adLevel;
      a_grids [ilev] = adLevel->m_grids;
      a_refRat[ilev] = adLevel->m_ref_ratio;
    }
}

/*******/
void
AMRLevelAdvectDiffuse::
defineSolvers()
{
  if (m_hasDiffusion)
    {
      CH_TIME("AMRLevelMHD::defineSolvers");

      int numSmooth, numMG, maxIter, mgverb;
      Real tolerance, hang, normThresh;

      ParmParse pp("amrmultigrid");
      pp.get("num_smooth", numSmooth);
      pp.get("num_mg",     numMG);
      pp.get("hang_eps",   hang);
      pp.get("norm_thresh",normThresh);
      pp.get("tolerance",  tolerance);
      pp.get("max_iter",   maxIter);
      pp.get("verbosity",  mgverb);

      Vector<AMRLevelAdvectDiffuse*>  hierarchy;
      Vector<DisjointBoxLayout>       grids;
      Vector<int>                     refRat;
      ProblemDomain                   lev0Dom;
      Real                            lev0Dx;
      getHierarchyAndGrids(hierarchy, grids, refRat, lev0Dom, lev0Dx);

      s_botSolver.m_verbosity = mgverb-3;

      AMRPoissonOpFactory* amrpop =  new AMRPoissonOpFactory();
      amrpop->define(lev0Dom, grids, refRat, lev0Dx, m_bcFunc, 1.0, m_nu);
      s_diffuseOpFact  = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >(amrpop);

      s_diffuseAMRMG = RefCountedPtr<AMRMultiGrid<LevelData<FArrayBox> > >
        (new AMRMultiGrid<LevelData<FArrayBox> >());

      s_diffuseAMRMG->define(lev0Dom, *s_diffuseOpFact, &s_botSolver, hierarchy.size());
      s_diffuseAMRMG->setSolverParameters(numSmooth, numSmooth, numSmooth, numMG,
                                          maxIter, tolerance, hang, normThresh);
      s_diffuseAMRMG->m_verbosity = mgverb;

      s_diffuseLevTGA = RefCountedPtr<LevelTGA>
        (new LevelTGA(grids, refRat, lev0Dom, s_diffuseOpFact, s_diffuseAMRMG));
    }
}

/********/
void
AMRLevelAdvectDiffuse::
define(const AdvectPhysics&        a_gphys,
       AdvectionVelocityFunction   a_advFunc,
       BCHolder                    a_bcFunc,
       const Real&                 a_cfl,
       const Real&                 a_domainLength,
       const Real&                 a_refineThresh,
       const int&                  a_tagBufferSize,
       const Real&                 a_initialDtMultiplier,
       const bool&                 a_useLimiting,
       const Real&                 a_nu)
{
  m_isDefined = true;
  m_cfl = a_cfl;
  m_domainLength = a_domainLength;
  m_refineThresh = a_refineThresh;
  m_tagBufferSize = a_tagBufferSize;
  m_initialDtMultiplier = a_initialDtMultiplier;
  m_useLimiting = a_useLimiting;
  m_nu = a_nu;
  m_doImplicitReflux = (m_nu > 0);
  m_hasDiffusion     = (m_nu > 0);

  m_advFunc = a_advFunc;
  m_bcFunc  = a_bcFunc;
  GodunovPhysics* gphysPtr = a_gphys.new_godunovPhysics();
  m_advPhys = RefCountedPtr<AdvectPhysics>((AdvectPhysics*)gphysPtr);
}

/********/
void AMRLevelAdvectDiffuse::define(AMRLevel*            a_coarserLevelPtr,
                                   const ProblemDomain& a_problemDomain,
                                   int                  a_level,
                                   int                  a_refRatio)
{
  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  if (a_coarserLevelPtr != NULL)
    {
      AMRLevelAdvectDiffuse* amrGodPtr = dynamic_cast<AMRLevelAdvectDiffuse*>(a_coarserLevelPtr);

      if (amrGodPtr != NULL)
        {
          define(*amrGodPtr->m_advPhys,
                 amrGodPtr->m_advFunc,
                 amrGodPtr->m_bcFunc,
                 amrGodPtr->m_cfl,
                 amrGodPtr->m_domainLength,
                 amrGodPtr->m_refineThresh,
                 amrGodPtr->m_tagBufferSize,
                 amrGodPtr->m_initialDtMultiplier,
                 amrGodPtr->m_useLimiting,
                 amrGodPtr->m_nu);
        }
      else
        {
          MayDay::Error("AMRLevelAdvectDiffuse::define: a_coarserLevelPtr is not castable to AMRLevelAdvectDiffuse*");
        }
    }

  // Compute the grid spacing
  m_dx = m_domainLength/a_problemDomain.domainBox().longside();

  m_numGhost = 4;
  m_stateNames  = Vector<string>(1, string("scalar"));
  m_advPhys->define(m_problem_domain, m_dx);
  PhysIBC* physIBCPtr = m_advPhys->getPhysIBC();
  physIBCPtr->define(m_problem_domain, m_dx);
}

/********/
void
AMRLevelAdvectDiffuse::
getCoarseDataPointers(LevelData<FArrayBox>** a_coarserDataOldPtr,
                      LevelData<FArrayBox>** a_coarserDataNewPtr,
                      LevelFluxRegister**    a_coarserFRPtr,
                      LevelFluxRegister**    a_finerFRPtr,
                      Real& a_tCoarserOld,
                      Real& a_tCoarserNew)
{
  *a_coarserDataOldPtr = NULL;
  *a_coarserDataNewPtr = NULL;
  *a_coarserFRPtr = NULL;
  *a_finerFRPtr   = NULL;

  a_tCoarserOld = 0.0;
  a_tCoarserNew = 0.0;

  // A coarser level exists
  if (m_hasCoarser)
    {
      AMRLevelAdvectDiffuse* coarserPtr = getCoarserLevel();

      // Recall that my flux register goes between my level and the next
      // finer level
      *a_coarserFRPtr = &coarserPtr->m_fluxRegister;

      *a_coarserDataOldPtr = &coarserPtr->m_UOld;
      *a_coarserDataNewPtr = &coarserPtr->m_UNew;

      a_tCoarserNew = coarserPtr->m_time;
      a_tCoarserOld = a_tCoarserNew - coarserPtr->m_dt;
    }

  // A finer level exists
  if (m_hasFiner)
    {
      // Recall that my flux register goes between my level and the next
      // finer level
      *a_finerFRPtr = &m_fluxRegister;
    }
}

/*******/
Real
AMRLevelAdvectDiffuse::
advance()
{
  if (s_verbosity >= 2)
    {
      pout() << "AMRLevelAdvectDiffuse::advance " << m_level << endl;
    }

  // Copy the new to the old
  m_UNew.copyTo(m_UNew.interval(),
                m_UOld,
                m_UOld.interval());

  LevelData<FArrayBox> diffusiveSrc(m_grids, 1, IntVect::Unit);
  makeDiffusiveSource(diffusiveSrc);

  Real newDt = diffusiveAdvance(diffusiveSrc);

  // Update the time and store the new timestep
  m_time += m_dt;
  Real returnDt = m_cfl * newDt;

  m_dtNew = returnDt;

  return returnDt;
}

/*********/
Real
AMRLevelAdvectDiffuse::
diffusiveAdvance(LevelData<FArrayBox>& a_diffusiveSrc)
{
  LevelFluxRegister* coarserFRPtr=NULL;
  LevelFluxRegister* finerFRPtr  =NULL;
  LevelData<FArrayBox>* coarserDataOldPtr = NULL;
  LevelData<FArrayBox>* coarserDataNewPtr = NULL;
  Real tCoarserOld, tCoarserNew;

  getCoarseDataPointers(&coarserDataOldPtr,
                        &coarserDataNewPtr,
                        &coarserFRPtr,
                        &finerFRPtr,
                        tCoarserOld, tCoarserNew);

  // Advance the solve one timestep
  Real newDt = m_levelGodunov.step(m_UNew,
                                   *finerFRPtr,
                                   *coarserFRPtr,
                                   m_advVel,
                                   a_diffusiveSrc,
                                   *coarserDataOldPtr,
                                   tCoarserOld,
                                   *coarserDataNewPtr,
                                   tCoarserNew,
                                   m_time,
                                   m_dt);
  if (m_hasDiffusion)
    {
      //compute Du = unew-uold and put uold back into unew
      //so we can advance using leveltga
      Interval interv(0, 0);
      m_UNew.copyTo(interv, m_dU,   interv);
      m_UOld.copyTo(interv, m_UNew, interv);
      for (DataIterator dit=m_dU.dataIterator(); dit.ok(); ++dit)
        {
          m_dU[dit()] -= m_UOld[dit()];
          m_dU[dit()] /= m_dt;
        }
      s_diffuseLevTGA->updateSoln(m_UNew, m_UOld, m_dU,
                                  finerFRPtr, coarserFRPtr,
                                  coarserDataOldPtr, coarserDataNewPtr,
                                  m_time, tCoarserOld, tCoarserNew,
                                  m_dt, m_level);
    }
  return newDt;
}

/*********/
void
AMRLevelAdvectDiffuse::
makeDiffusiveSource(LevelData<FArrayBox>& a_diffusiveSrc)
{
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      a_diffusiveSrc[dit()].setVal(0);
    }
  if (m_hasDiffusion)
    {

      RefCountedPtr<AMRPoissonOp> amrpop = RefCountedPtr<AMRPoissonOp>((AMRPoissonOp*)s_diffuseOpFact->AMRnewOp(m_problem_domain));
      LevelData<FArrayBox> zero;
      amrpop->create(zero, a_diffusiveSrc);
      amrpop->setToZero(zero);
      amrpop->setAlphaAndBeta(0., 1.0);
      if (m_level == 0)
        {
          amrpop->residual(a_diffusiveSrc, m_UOld, zero, false);
        }
      else
        {
          LevelFluxRegister* coarserFRPtr=NULL;
          LevelFluxRegister* finerFRPtr  =NULL;
          LevelData<FArrayBox>* coarserDataOldPtr = NULL;
          LevelData<FArrayBox>* coarserDataNewPtr = NULL;
          Real tCoarserOld, tCoarserNew;

          getCoarseDataPointers(&coarserDataOldPtr,
                                &coarserDataNewPtr,
                                &coarserFRPtr,
                                &finerFRPtr,
                                tCoarserOld, tCoarserNew);

          LevelData<FArrayBox> UCoarse(coarserDataOldPtr->disjointBoxLayout(),
                                       coarserDataOldPtr->nComp(),
                                       coarserDataOldPtr->ghostVect());

          Real interpTime = m_time + 0.5*m_dt;
          interpolateInTime(UCoarse, *coarserDataOldPtr, *coarserDataNewPtr,
                            interpTime, tCoarserOld, tCoarserNew);

          amrpop->AMRResidualNF(a_diffusiveSrc, m_UOld, UCoarse, zero, false);
        }

      amrpop->scale(a_diffusiveSrc, -1.0);

      /// Over the coarse-fine interface, the diffusive source is set
      /// to zero. At fine-fine interfaces, it is filled in by
      /// neighboring boxes.
      a_diffusiveSrc.exchange();
    }
}

/*******/
void
AMRLevelAdvectDiffuse::
interpolateInTime(LevelData<FArrayBox>&          a_interp,
                  const LevelData<FArrayBox>&    a_old,
                  const LevelData<FArrayBox>&    a_new,
                  Real a_time, Real a_tOld, Real a_tNew)
{
  CH_assert(a_tNew >= a_tOld);
  //interp = alpha* unew + (1-alpha) uold
  Real alpha = 0;
  Real diff = a_tNew-a_tOld;
  if (diff > 0)
    {
      CH_assert(a_time >= a_tOld);
      alpha = (a_time - a_tOld)/diff;
    }
  for (DataIterator dit = a_interp.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox tempInte(a_interp[dit()].box(), a_interp.nComp());
      tempInte.copy(       a_new[dit()]);
      a_interp[dit()].copy(a_old[dit()]);
      tempInte        *= alpha;     //temp has alpha*unew
      a_interp[dit()] *= 1.0-alpha; //interp has (1-alpha)uold
      a_interp[dit()] += tempInte;  //interp now has (1-alpha)uold + alpha*unew
    }
}
void
AMRLevelAdvectDiffuse::
setSolverCoef(Real a_alpha, Real a_beta)
{
  // now set alpha and beta of all the operators in the solver.
  // this includes resetting the relaxation coefficient
  Vector<MGLevelOp<LevelData<FArrayBox> >* > ops = s_diffuseAMRMG->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
    {
      LevelTGAHelmOp<FArrayBox,FluxBox>* helmop = (LevelTGAHelmOp<FArrayBox,FluxBox>*) ops[iop];
      helmop->setAlphaAndBeta(a_alpha, a_beta);
    }
}

/*******/
void
AMRLevelAdvectDiffuse::
postTimeStep()
{
  if (s_verbosity >= 2)
    {
      pout() << "AMRLevelAdvectDiffuse::postTimeStep " << m_level << endl;
    }

  if (m_hasFiner)
    {
      if (m_doImplicitReflux)
        {
          doImplicitReflux();
        }
      else
        {
          // explicit Reflux
          Real scale = -1.0/m_dx;
          m_fluxRegister.reflux(m_UNew,scale);
        } // end if we're doing explicit refluxing

      // Average from finer level data
      AMRLevelAdvectDiffuse* amrGodFinerPtr = getFinerLevel();

      amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                      amrGodFinerPtr->m_UNew);
    } // end if there is a finer level

}

/*******/
void
AMRLevelAdvectDiffuse::
doImplicitReflux()
{
  // first find out what time coarser level is (if it exists)
  Real time_eps = 1.0e-10;
  // do multilevel operations if this is the coarsest level or if
  // coarser level is not at same time as this level. otherwise,
  // defer this until we get down to the coarsest level at this time.
  if (m_level == 0 || (abs(m_coarser_level_ptr->time() - m_time) > time_eps))
    {
      Vector<AMRLevelAdvectDiffuse*>         hierarchy;
      Vector<int>                            refRat;
      Vector<DisjointBoxLayout>              grids;
      Real                                   lev0Dx;
      ProblemDomain                          lev0Domain;
      getHierarchyAndGrids(hierarchy, grids, refRat, lev0Domain, lev0Dx);

      int finest_level = hierarchy.size()-1;

      // now do implicit refluxing
      // Vector of pointers to LevelData of FABS
      Vector<LevelData<FArrayBox>* > refluxCor(finest_level+1, NULL);
      Vector<LevelData<FArrayBox>* > refluxRHS(finest_level+1, NULL);
      // collectUN: AMR vector containing soln at new time
      Vector<LevelData<FArrayBox>* > collectUN(finest_level+1, NULL);

      // loop over levels, allocate storage, set up for AMRSolve
      // if coarser level exists, define it as well for BCs.
      int startLev = Max(m_level-1, 0);

      for (int lev = startLev; lev<= finest_level; lev++)
        {
          // rhs has no ghost cells, correction does
          refluxRHS[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Zero);
          refluxCor[lev]  = new LevelData<FArrayBox>(grids[lev], 1, IntVect::Unit);
          collectUN[lev]  = &(hierarchy[lev]->m_UNew);
          for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit)
            {
              (*refluxRHS[lev])[dit()].setVal(0.0);
              (*refluxCor[lev])[dit()].setVal(0.0);
            }
        } // end loop over levels for setup.

      // now loop over levels and set up RHS
      // note that we don't need to look at finest level here,
      // since it will have no finer level to get a reflux correction
      // from.   Also this starts at m_level instead of startLev since,
      // in the case m_level > 0, m_level-1 is only used for boundary conditions
      for (int lev= m_level; lev < finest_level; lev++)
        {
          Real dxLev = hierarchy[lev]->m_dx;
          Real refluxScale = 1.0/dxLev;

          hierarchy[lev]->m_fluxRegister.reflux(*refluxRHS[lev], refluxScale);
        }


      int lbase = m_level;
      int lmax  = finest_level;
      // this resets the coeffients including eta, alpha, beta
      Real alpha = 1.0;
      Real beta  = -m_dt;
      setSolverCoef(alpha, beta);
      s_diffuseAMRMG->solve(refluxCor, refluxRHS, lmax, lbase);

      for (int lev= m_level; lev <= finest_level; lev++)
        {
          for (DataIterator dit = grids[lev].dataIterator(); dit.ok(); ++dit)
            {
              (*hierarchy[lev]).m_UNew[dit()] += (*refluxCor[lev])[dit()];
            }
        }

      //remember that startLev can be different from m_level
      for (int lev = startLev; lev<= finest_level; lev++)
        {
          delete refluxRHS[lev];
          delete refluxCor[lev];
          refluxRHS[lev] = NULL;
          refluxCor[lev] = NULL;
        }
    } //end if times are in alignment
}

/*******/
void
AMRLevelAdvectDiffuse::
tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::tagCells " << m_level << endl;
    }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::tagCellsInit " << m_level << endl;
    }

  // Create tags based on undivided gradient of density
  IntVectSet localTags;
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
    {
      const AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();

      PiecewiseLinearFillPatch pwl(levelDomain,
                                   amrGodCoarserPtr->m_UNew.disjointBoxLayout(),
                                   1,
                                   amrGodCoarserPtr->m_problem_domain,
                                   amrGodCoarserPtr->m_ref_ratio,
                                   1);

      pwl.fillInterp(m_UNew,
                     amrGodCoarserPtr->m_UNew,
                     amrGodCoarserPtr->m_UNew,
                     1.0,
                     0,
                     0,
                     1);
    }
  m_UNew.exchange(Interval(0,1-1));

  // Compute undivided gradient
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& b = levelDomain[dit()];
      FArrayBox gradFab(b,SpaceDim);
      const FArrayBox& UFab = m_UNew[dit()];

      for (int dir = 0; dir < SpaceDim; ++dir)
        {
          const Box bCenter = b & grow(m_problem_domain,-BASISV(dir));
          const Box bLo     = b & adjCellLo(bCenter,dir);
          const int hasLo = ! bLo.isEmpty();
          const Box bHi     = b & adjCellHi(bCenter,dir);
          const int hasHi = ! bHi.isEmpty();
          FORT_GETGRADF(CHF_FRA1(gradFab,dir),
                        CHF_CONST_FRA1(UFab,0),
                        CHF_CONST_INT(dir),
                        CHF_BOX(bLo),
                        CHF_CONST_INT(hasLo),
                        CHF_BOX(bHi),
                        CHF_CONST_INT(hasHi),
                        CHF_BOX(bCenter));
        }

      FArrayBox gradMagFab(b,1);
      FORT_MAGNITUDEF(CHF_FRA1(gradMagFab,0),
                      CHF_CONST_FRA(gradFab),
                      CHF_BOX(b));

      // Tag where gradient exceeds threshold
      BoxIterator bit(b);
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();

          if (gradMagFab(iv) >= m_refineThresh)
            {
              localTags |= iv;
            }
        }
    }

  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;

  a_tags = localTags;
}

/*******/
void
AMRLevelAdvectDiffuse::
tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}

/*******/
void
AMRLevelAdvectDiffuse::
regrid(const Vector<Box>& a_newGrids)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::regrid " << m_level << endl;
    }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);

  // Save data for later
  LevelData<FArrayBox> UOld;
  UOld.define(m_UNew);
  m_UNew.copyTo(m_UNew.interval(),
                UOld,
                UOld.interval());

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,1,ivGhost);
  m_UOld.define(m_grids,1,ivGhost);
  m_dU.define(m_grids,1,ivGhost);
  m_advVel.define(m_grids,1,ivGhost);
  fillAdvectionVelocity();

  // Set up data structures
  levelSetup();

  // Interpolate from coarser level
  if (m_hasCoarser)
    {
      AMRLevelAdvectDiffuse* amrGodCoarserPtr = getCoarserLevel();
      m_fineInterp.interpToFine(m_UNew,amrGodCoarserPtr->m_UNew);
    }

  // Copy from old state
  UOld.copyTo(UOld.interval(),
              m_UNew,
              m_UNew.interval());
}

/*******/
void
AMRLevelAdvectDiffuse::
initialGrid(const Vector<Box>& a_newGrids)
{

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  Vector<int> procs;
  LoadBalance(procs, a_newGrids);
  m_grids = DisjointBoxLayout(a_newGrids, procs, m_problem_domain);

  // Define old and new state data structures
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,1,ivGhost);
  m_UOld.define(m_grids,1,ivGhost);
  m_dU.define(m_grids,1,ivGhost);
  m_advVel.define(m_grids,1,ivGhost);
  fillAdvectionVelocity();

  // Set up data structures
  levelSetup();
}

/*******/
void
AMRLevelAdvectDiffuse::
initialData()
{

  PhysIBC* physIBCPtr = m_advPhys->getPhysIBC();
  physIBCPtr->initialize(m_UNew);
  physIBCPtr->initialize(m_UOld);
}

/*******/
void
AMRLevelAdvectDiffuse::
postInitialize()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::postInitialize " << m_level << endl;
    }

  // solve for post-regrid smoothing may eventually go here.
}

#ifdef CH_USE_HDF5

/*******/
void
AMRLevelAdvectDiffuse::
writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::writeCheckpointHeader" << endl;
    }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = 1;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < 1; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = m_stateNames[comp];
    }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }
}

/*******/
void
AMRLevelAdvectDiffuse::
writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::writeCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_int ["tag_buffer_size"] = m_tagBufferSize;
  header.m_real["dx"]              = m_dx;
  header.m_real["dt"]              = m_dt;
  header.m_real["time"]            = m_time;
  header.m_real["nu"]              = m_nu;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();


  // Setup the periodicity info
  D_TERM(
         if (m_problem_domain.isPeriodic(0))
           header.m_int ["is_periodic_0"] = 1;
         else
           header.m_int ["is_periodic_0"] = 0; ,

         if (m_problem_domain.isPeriodic(1))
           header.m_int ["is_periodic_1"] = 1;
         else
           header.m_int ["is_periodic_1"] = 0; ,

         if (m_problem_domain.isPeriodic(2))
           header.m_int ["is_periodic_2"] = 1;
         else
           header.m_int ["is_periodic_2"] = 0; );

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_UNew,"data");
}

/*******/
void
AMRLevelAdvectDiffuse::
readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::readCheckpointHeader" << endl;
    }

  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << "hdf5 header data:" << endl;
      pout() << header << endl;
    }

  // Get the number of components
  if (header.m_int.find("num_components") == header.m_int.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointHeader: checkpoint file does not have num_components");
    }

  int numStates = header.m_int["num_components"];
  if (numStates != 1)
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointHeader: num_components in checkpoint file does not match solver");
    }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < 1; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      if (header.m_string.find(compStr) == header.m_string.end())
        {
          MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointHeader: checkpoint file does not have enough component names");
        }

      stateName = header.m_string[compStr];
      if (stateName != m_stateNames[comp])
        {
          MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointHeader: state_name in checkpoint does not match solver");
        }
    }
}

/*******/
void
AMRLevelAdvectDiffuse::
readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::readCheckpointLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << "hdf5 header data:" << endl;
      pout() << header << endl;
    }

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
    {
      pout() << "read ref_ratio = " << m_ref_ratio << endl;
    }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain tag_buffer_size");
    }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
    {
      pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
    }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain dx");
    }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
    {
      pout() << "read dx = " << m_dx << endl;
    }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
    {
      pout() << "read dt = " << m_dt << endl;
    }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain time");
    }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
    {
      pout() << "read time = " << m_time << endl;
    }

  // Get nu
  if (header.m_real.find("nu") == header.m_real.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain nu");
    }
  m_nu = header.m_real["nu"];

  if (s_verbosity >= 2)
    {
      pout() << "read nu = " << m_nu << endl;
    }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain prob_domain");
    }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
           isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
         else
           isPeriodic[0] = false; ,

         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
           isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
         else
           isPeriodic[1] = false; ,

         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
           isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
         else
           isPeriodic[2] = false;);

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);

  if (gridStatus != 0)
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain a Vector<Box>");
    }

  // Create level domain
  Vector<int> procs;
  LoadBalance(procs, grids);
  m_grids = DisjointBoxLayout(grids,procs, m_problem_domain);

  // Indicate/guarantee that the indexing below is only for reading
  // otherwise an error/assertion failure occurs
  const DisjointBoxLayout& constGrids = m_grids;

  LayoutIterator lit = constGrids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = constGrids[lit()];
      m_level_grids.push_back(b);
    }

  if (s_verbosity >= 4)
    {
      pout() << "read level domain: " << endl;
      LayoutIterator lit = m_grids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& b = m_grids[lit()];
          pout() << lit().intCode() << ": " << b << endl;
        }
      pout() << endl;
    }

  // Reshape state with new grids
  m_UNew.define(m_grids,1);
  m_dU.define(  m_grids,1);
  const int dataStatus = read<FArrayBox>(a_handle,
                                         m_UNew,
                                         "data",
                                         m_grids);

  if (dataStatus != 0)
    {
      MayDay::Error("AMRLevelAdvectDiffuse::readCheckpointLevel: file does not contain state data");
    }
  m_UOld.define(m_grids,1);
  m_advVel.define(m_grids,1, m_numGhost*IntVect::Unit);
  fillAdvectionVelocity();

  // Set up data structures
  levelSetup();
}

/*******/
void
AMRLevelAdvectDiffuse::
writePlotHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::writePlotHeader" << endl;
    }

  // Setup the number of components -- include space for error
  HDF5HeaderData header;

  int numPlotComp = 1;
  header.m_int["num_components"] = numPlotComp;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < 1; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = m_stateNames[comp];
    }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }
}

/*******/
void
AMRLevelAdvectDiffuse::
writePlotLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::writePlotLevel" << endl;
    }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx;
  header.m_real["dt"]          = m_dt;
  header.m_real["time"]        = m_time;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // Write the data for this level
  int numPlotVar = 1;
  LevelData<FArrayBox> plotData(m_UNew.getBoxes(), numPlotVar);

  // first copy data to plot data holder
  m_UNew.copyTo(m_UNew.interval(), plotData, m_UNew.interval());

  write(a_handle,m_UNew.boxLayout());
  write(a_handle,plotData,"data");

}
#endif

/*******/
Real
AMRLevelAdvectDiffuse::
computeDt()
{
  Real newDT = m_cfl * m_dx
    / m_levelGodunov.getMaxWaveSpeed(m_UNew, m_advVel);
  return newDT;
}

/*******/
Real
AMRLevelAdvectDiffuse::
computeInitialDt()
{
  Real newDT = m_initial_dt_multiplier * m_dx
    / m_levelGodunov.getMaxWaveSpeed(m_UNew, m_advVel);

  return newDT;
}

/*******/
void
AMRLevelAdvectDiffuse::
levelSetup()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelAdvectDiffuse::levelSetup " << m_level << endl;
    }

  AMRLevelAdvectDiffuse* amrADCoarserPtr = getCoarserLevel();
  AMRLevelAdvectDiffuse* amrADFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrADCoarserPtr != NULL);
  m_hasFiner   = (amrADFinerPtr   != NULL);

  int nRefCrse = -1;
  DisjointBoxLayout* crseGridsPtr = NULL;

  if (m_hasCoarser)
    {
      nRefCrse = m_coarser_level_ptr->refRatio();
      crseGridsPtr = &amrADCoarserPtr->m_grids;

      const DisjointBoxLayout& coarserLevelDomain = amrADCoarserPtr->m_grids;



      m_coarseAverage.define(m_grids,
                             1,
                             nRefCrse);

      m_fineInterp.define(m_grids,
                          1,
                          nRefCrse,
                          m_problem_domain);


      // Maintain levelGodunov
      m_levelGodunov.define(*m_advPhys,
                            m_grids,
                            coarserLevelDomain,
                            m_problem_domain,
                            nRefCrse,
                            m_useLimiting,
                            m_dx,
                            m_hasCoarser,
                            m_hasFiner);

      // This may look twisted but you have to do this this way because the
      // coarser levels get setup before the finer levels so, since a flux
      // register lives between this level and the next FINER level, the finer
      // level has to do the setup because it is the only one with the
      // information at the time of construction.

      // Maintain flux registers
      amrADCoarserPtr->m_fluxRegister.define(m_grids,
                                             amrADCoarserPtr->m_grids,
                                             m_problem_domain,
                                             amrADCoarserPtr->m_ref_ratio,
                                             1);
      amrADCoarserPtr->m_fluxRegister.setToZero();
    }
  else
    {
      m_levelGodunov.define(*m_advPhys,
                            m_grids,
                            DisjointBoxLayout(),
                            m_problem_domain,
                            m_ref_ratio,
                            m_useLimiting,
                            m_dx,
                            m_hasCoarser,
                            m_hasFiner);
    }

  defineSolvers();
}

/*******/
AMRLevelAdvectDiffuse*
AMRLevelAdvectDiffuse::
getCoarserLevel() const
{
  AMRLevelAdvectDiffuse* amrADCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
    {
      amrADCoarserPtr = dynamic_cast<AMRLevelAdvectDiffuse*>(m_coarser_level_ptr);

      if (amrADCoarserPtr == NULL)
        {
          MayDay::Error("AMRLevelAdvectDiffuse::getCoarserLevel: dynamic cast failed");
        }
    }

  return amrADCoarserPtr;
}

/*******/
AMRLevelAdvectDiffuse*
AMRLevelAdvectDiffuse::
getFinerLevel() const
{
  AMRLevelAdvectDiffuse* amrADFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
    {
      amrADFinerPtr = dynamic_cast<AMRLevelAdvectDiffuse*>(m_finer_level_ptr);

      if (amrADFinerPtr == NULL)
        {
          MayDay::Error("AMRLevelAdvectDiffuse::getFinerLevel: dynamic cast failed");
        }
    }

  return amrADFinerPtr;
}

/*******/
void
AMRLevelAdvectDiffuse::
fillAdvectionVelocity()
{
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FluxBox& advVelFlux = m_advVel[dit()];
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          FArrayBox& velDir = advVelFlux[faceDir];
          const Box& fabBox = velDir.box();
          for (BoxIterator bit(fabBox); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              RealVect loc;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  if (idir == faceDir)
                    {
                      loc[idir] = m_dx*Real(iv[idir]);
                    }
                  else
                    {
                      loc[idir] = m_dx*(0.5 + Real(iv[idir]));
                    }
                }
              velDir(iv, 0) = m_advFunc(loc, faceDir);
            }
        }
    }
}
/*******/

#include "NamespaceFooter.H"
