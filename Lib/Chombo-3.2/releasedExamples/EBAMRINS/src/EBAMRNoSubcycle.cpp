#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ParmParse.H"
#include "EBAMRNoSubcycle.H"
#include "EBAMRNoSubcycleF_F.H"
#include "EBAMRPoissonOpF_F.H"
#include "MeshRefine.H"
#include "BRMeshRefine.H"
#include "EBEllipticLoadBalance.H"
#include "EBArith.H"
#include "EBPWLFineInterp.H"
#include "EBCoarseAverage.H"
#include "EBFluxFactory.H"
#include "EBCellFactory.H"
#include "EBLevelAdvect.H"
#include "EBGradDivFilter.H"
#include "EBPatchAdvect.H"
#include "REAL.H"
#include "EBPhysIBCFactory.H"
#include "EBAMRIO.H"
#include "BaseIFFactory.H"
#include "EBLevelRedist.H"
#include "BaseIVFactory.H"
#include "EBConstantCFInterp.H"
#include "EBArith.H"
#include "EBAMRDataOps.H"
#include "NeumannPoissonEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "InflowOutflowIBC.H"
#include "EBNormalizeByVolumeFraction.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include "memusage.H"
#include "memtrack.H"

#include "NamespaceHeader.H"

extern Real g_simulationTime;
#define debugIV IntVect(D_DECL(16, 3, 0))

/**********************/
EBAMRNoSubcycle::
EBAMRNoSubcycle(const AMRParameters&      a_params,
                const EBIBCFactory&       a_ibcfact,
                const ProblemDomain&      a_coarsestDomain,
                Real                      a_viscosity,
                const EBIndexSpace* const a_ebisPtr):
  m_ebisPtr(a_ebisPtr)
{
  CH_TIME("EBAMRNoSubcycle::EBAMRNoSubcycle");
  if (a_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::EBAMRNoSubcycle" << endl;
    }

  //set parameters of the run
  m_params    = a_params;
  m_viscosity = a_viscosity;
  m_viscousCalc = (m_viscosity > 0);

  //create initial and boundary condition object
  m_ibc    =   a_ibcfact.create();

  //resize vectors and set them where we can
  Real coarsestDx = m_params.m_domainLength/Real(a_coarsestDomain.size(0));
  int nlevels = m_params.m_maxLevel + 1;
  m_domain.resize(nlevels);
  m_dx.resize(nlevels);
  m_grids.resize(nlevels);
  m_ebisl.resize(nlevels);
  m_eblg.resize(nlevels);

  m_quadCFI.resize(nlevels);
  m_aveOper.resize(nlevels);
  m_aveSpac.resize(nlevels);
  m_ebLevAd.resize(nlevels);
  m_fluxReg.resize(nlevels);
  m_velo.resize(nlevels, NULL);
  m_pres.resize(nlevels, NULL);
  m_gphi.resize(nlevels, NULL);
  m_advVel.resize(nlevels, NULL);
  m_coveredAdvVelLo.resize(nlevels, NULL);
  m_coveredAdvVelHi.resize(nlevels, NULL);
  m_coveredFaceLitLo.resize(nlevels, NULL);
  m_coveredFaceLitHi.resize(nlevels, NULL);
  m_coveredSetsLitLo.resize(nlevels, NULL);
  m_coveredSetsLitHi.resize(nlevels, NULL);

  allocateDataHolders();

  m_domain[0] = a_coarsestDomain;
  m_dx[0]     =   coarsestDx;
  for (int ilev = 1; ilev < nlevels; ilev++)
    {
      CH_assert(m_params.m_refRatio[ilev-1] > 0);
      m_domain[ilev] = refine(m_domain[ilev-1], m_params.m_refRatio[ilev-1]);
      m_dx[ilev] = m_dx[ilev-1]/Real(m_params.m_refRatio[ilev-1]);
    }
  m_prescribedDt = -1.0;
  m_useFixedDt = false;
  m_steadyState = false;
  m_stopAdvance = false;

  m_ccProjector  = NULL;
  m_macProjector = NULL;
  m_time = 0.0;
  m_curStep = 0;
  m_dt = -1.0;
  //setup still needs to get called
  m_isSetup  = false;

  m_pointsUpdated = 0;
}
/**********/
void
EBAMRNoSubcycle::
allocateDataHolders()
{
  CH_TIME("EBAMRNoSubcycle::allocateDataHolders");
  for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      m_velo[ilev]   = new LevelData<EBCellFAB>();
      m_pres[ilev]   = new LevelData<EBCellFAB>();
      m_gphi[ilev]   = new LevelData<EBCellFAB>();
      m_advVel[ilev] = new LevelData<EBFluxFAB>();

      m_coveredFaceLitLo[ilev] = new LayoutData< Vector< Vector<VolIndex> > >();
      m_coveredFaceLitHi[ilev] = new LayoutData< Vector< Vector<VolIndex> > >();
      m_coveredSetsLitLo[ilev] = new LayoutData< Vector< IntVectSet > >      ();
      m_coveredSetsLitHi[ilev] = new LayoutData< Vector< IntVectSet > >      ();

      m_coveredAdvVelLo[ilev]  = new LayoutData<Vector<BaseIVFAB<Real> * > >();
      m_coveredAdvVelHi[ilev]  = new LayoutData<Vector<BaseIVFAB<Real> * > >();
    }
}
/**********/
EBAMRNoSubcycle::
~EBAMRNoSubcycle()
{
  CH_TIME("EBAMRNoSubcycle::~EBAMRNoSubcycle");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::~EBAMRNoSubcycle" << endl;
    }
  delete m_ibc;
  for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      delete m_velo[ilev];
      delete m_pres[ilev];
      delete m_gphi[ilev];
      delete m_advVel[ilev];

      delete m_coveredFaceLitLo[ilev];
      delete m_coveredFaceLitHi[ilev];
      delete m_coveredSetsLitLo[ilev];
      delete m_coveredSetsLitHi[ilev];

      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*m_coveredAdvVelLo[ilev])[dit()][idir];
              delete (*m_coveredAdvVelHi[ilev])[dit()][idir];
            }
        }
      delete m_coveredAdvVelLo[ilev];
      delete m_coveredAdvVelHi[ilev];
    }
  for (int ilev=0; ilev<m_params.m_maxLevel; ilev++)
    {
      m_advVel[ilev]           = NULL;
      m_coveredAdvVelLo[ilev]  = NULL;
      m_coveredAdvVelHi[ilev]  = NULL;
    }
  if (m_ccProjector !=  NULL)
    {
      delete m_ccProjector;
    }
  if (m_macProjector !=  NULL)
    {
      delete m_macProjector;
    }
}
/**********/
void
EBAMRNoSubcycle::
setupForAMRRun()
{
  CH_TIME("EBAMRNoSubcycle::setupForAMRMRun");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::setupForAMRRun" << endl;
    }
  m_isSetup= true;
  m_doRestart    = false;
  //we're generating the hiearachy dyamically
  //modelled on AMR::initialGrid()

  Vector<Vector<Box> > old_boxes(1);
  Vector<Vector<Box> > new_boxes;

  //also keep old tags around
  Vector<IntVectSet> oldTags;

  //define base mesh
  //chop up base level into grids to satisfy box size requirements
  domainSplit(m_domain[0], old_boxes[0], m_params.m_maxBoxSize,
              m_params.m_blockFactor);

  m_finestLevel = 0;

  //now generate more levels if necessary
  int top_level = 0;

  bool moreLevels = (m_params.m_maxLevel > 0);

  //create grid generation object
  BRMeshRefine meshrefine;

  if (moreLevels)
    {
      meshrefine.define(m_domain[0],              m_params.m_refRatio,
                        m_params.m_fillRatio,     m_params.m_blockFactor,
                        m_params.m_nestingRadius, m_params.m_maxBoxSize);
    }

  //now intialize data for existing hierarchy
  initialGrid(old_boxes);

  initialData();

  while (moreLevels)
    {
      //default is moreLevels = false
      //(only repeat loop in the case where a new level
      //is generated which is still coarser than maxLevel)
      moreLevels = false;

      int base_level = 0;
      int old_top_level = top_level;

      Vector<IntVectSet> tagsVect(top_level+1);
      tagCells(tagsVect);

      int new_finest = meshrefine.regrid(new_boxes, tagsVect,
                                         base_level, top_level,
                                         old_boxes);
      if (new_finest > top_level) top_level++;

      old_boxes = new_boxes;

      //now see if we need another pass through grid generation
      if ((top_level<m_params.m_maxLevel) && (top_level > old_top_level))
        moreLevels = true;


      //if we added another level, reinintialize everything again
      if (top_level > old_top_level)
        {
          initialGrid(new_boxes);

          initialData();
        }
    } //end loop over regridding passes

  defineIrregularData();
  //finally, call post-initialization
  postInitialize();

}
void
EBAMRNoSubcycle::
filter(Vector<LevelData<EBCellFAB>* >&   a_veloc)
{
  //finestlevel vs. maxlevel fix
  CH_TIME("EBAMRNoSubcycle::filter");
  for (int ifilter = 0; ifilter < m_params.m_numFilterIterations; ifilter++)
    {
      averageDown(a_veloc);
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          const DisjointBoxLayout*           coarGridsPtr = NULL;
          const EBISLayout*                  coarEBISLPtr = NULL;
          const LevelData<EBCellFAB>*        coarVelocPtr = NULL;
          int refRat = -1;
          if (ilev > 0)
            {
              coarGridsPtr = &m_grids[ilev-1];
              coarEBISLPtr = &m_ebisl[ilev-1];
              coarVelocPtr =  a_veloc[ilev-1];
              refRat = m_params.m_refRatio[ilev-1];
            }
          Real lambdaScale = 0.125; // default value in EBGradDivFilter.H
          EBGradDivFilter gdFilt(m_grids[ilev],
                                 coarGridsPtr,
                                 m_ebisl[ilev],
                                 coarEBISLPtr,
                                 m_domain[ilev],
                                 m_dx[ilev]*RealVect::Unit,
                                 refRat,
                                 lambdaScale,
                                 m_ebisPtr);
          //flux velocity used for boundary conditions.  this assumes (correctly) that the
          //projection has just happened.
          Vector< LevelData<EBFluxFAB>* >& fluxVel = m_ccProjector->getMacVelocity();
          gdFilt.filter(*a_veloc[ilev], *fluxVel[ilev], coarVelocPtr);

        }
    }
}
/**************************/
Real
EBAMRNoSubcycle::
run(Real a_maxTime, int a_maxStep)
{
  CH_TIME("EBAMRNoSubcycle::run");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::run" << endl;
    }
  CH_assert(m_isSetup);

  ParmParse pp;

  pp.query("stokesFlow", m_params.m_stokesFlow);
  if (m_params.m_stokesFlow)
    {
      pout() << "Doing Stokes flow" << endl;
    }
  else
    {
      pout() << "Doing advective flow" << endl;
    }

  pp.query("use_stokes_dt", m_params.m_useStokesDt);
  if (m_params.m_useStokesDt)
    {
      pout() << "Using Stokes dt" << endl;
    }

  pp.query("do_steady_state", m_params.m_doSteadyState);
  if (m_params.m_doSteadyState)
    {
      pout() << "Stopping velocity predictor-corrector at steady-state" << endl;
    }

  //only call computeInitialDt if we're not doing a restart
  if (!m_doRestart)
    {
      // do initial velocity projection here instead of
      // in postInitialize to save time if we are not using the object for
      // a run
      if (m_params.m_verbosity >= 2)
        {
          pout() << "EBAMRNoSubcycle: cc projecting initial velocity" << endl;
        }
      Interval interv(0, SpaceDim-1);
      Vector<LevelData<EBCellFAB>* > tempLDPtr2;
      tempLDPtr2.resize(m_finestLevel+1);
      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          tempLDPtr2[ilev] = new LevelData<EBCellFAB>();
          tempLDPtr2[ilev]->define(m_grids[ilev], SpaceDim,  IntVect::Zero, ebcellfact);
          m_gphi[ilev]->copyTo(interv, *(tempLDPtr2[ilev]), interv);
        }

      // set time used in setting BCs ahead of projection
      Real time = 0.0;
      m_ccProjector->setTime(time);
      m_ccProjector->project(m_velo, m_gphi);
      filter(m_velo);

      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          tempLDPtr2[ilev]->copyTo(interv, *m_gphi[ilev], interv);
        }

      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          delete tempLDPtr2[ilev];
        }

      m_advanceGphiOnly = true;

      m_curStep = 0;
      m_time = time;
      m_dt = computeInitialDt();
      if ((m_time + m_dt) > a_maxTime)
      {
        m_dt = a_maxTime - m_time;
      }

      if (m_params.m_initIterations>0)
      {
        pout() << "EBAMRNoSubcycle: iterating to get initial pressure gradient" << endl;
      }
      for (int iter = 0; iter < m_params.m_initIterations; iter++)
        {
          pout() << "##############################################################################" << endl;
          pout() << "EBAMRNoSubcycle: gphi iteration " << iter  << ", dt = " <<  m_dt  << endl;
          advance();
          postTimeStep();
        }

      m_advanceGphiOnly = false;

      //dump plotfile and checkpointfile before regridding
#ifdef CH_USE_HDF5
      if ((m_curStep%m_params.m_plotInterval == 0) && (m_params.m_plotInterval > 0))
        {
          if (m_params.m_verbosity > 1)
            {
              pout() << "EBAMRNoSubcycle::writePlotFile" << endl;
            }
          writePlotFile();
        }

      if ((m_curStep%m_params.m_checkpointInterval == 0) && (m_params.m_checkpointInterval  > 0))
        {
          if (m_params.m_verbosity > 1)
            {
              pout() << "EBAMRNoSubcycle::writeCheckpointFile" << endl;
            }
          writeCheckpointFile();
        }
#endif
    }
  else
    {
      m_dt = computeDt();
      if ((m_time + m_dt) > a_maxTime)
      {
        m_dt = a_maxTime - m_time;
      }
    }
  m_advanceGphiOnly = false;

  if (m_params.m_verbosity > 0)
    {
      pout() << "EBAMRNoSubcycle: starting run " << endl;
    }

  //skip computeDt if in the first step since it was already calculated above
  bool skipDtComputation = true;

  Real initDt = m_dt;

  pout() << "EBAMRNoSubcycle: starting time advance  " << endl;

  //advance solution until done
  while ((a_maxTime > m_time) && (m_curStep < a_maxStep) && (!m_stopAdvance))
    {
      g_simulationTime = m_time;

      //advance step number
      ++m_curStep;

      pout() << "##############################################################################" << endl;
      pout() << "EBAMRNoSubcycle: step " << m_curStep << endl;

      //do regridding if appropriate
      if ( (m_curStep%m_params.m_regridInterval == 0) &&
           (m_params.m_regridInterval > 0) &&
           (m_params.m_maxLevel > 0) )
        {
          regrid();
        }

      //compute new dt
      if (skipDtComputation)//do not skip computation in future iterations
        {
          skipDtComputation = false;
        }
      else
        {
          m_dt = computeDt();
        }

      if (m_dt < 1e-5 * initDt)
      {
        pout() << "EBAMRNoSubcycle::run -- Time step too small" << endl;
        break;
      }

      if ((m_time + m_dt) > a_maxTime)
      {
        m_dt = a_maxTime - m_time;
      }

      pout() << "Beginning of time step " << m_curStep
             << ", start time = " << m_time
             << ", dt = " << m_dt << endl;

      //do timestep
      advance();

      postTimeStep();

      //advance time
      m_time += m_dt;

      //dump plotfile and checkpointfile before regridding
#ifdef CH_USE_HDF5
      if ((m_curStep%m_params.m_plotInterval == 0) && (m_params.m_plotInterval > 0))
        {
          if (m_params.m_verbosity > 0)
            {
              pout() << "EBAMRNoSubcycle: writing plot file" << endl;
            }
          writePlotFile();
        }

      if ((m_curStep%m_params.m_checkpointInterval == 0) && (m_params.m_checkpointInterval > 0))
        {
          if (m_params.m_verbosity > 0)
            {
              pout() << "EBAMRNoSubcycle: writing checkpoint file" << endl;
            }
          writeCheckpointFile();
        }
#endif

      pout() << "End of time step " << m_curStep << ", end time = " << m_time << endl;

    }//end while loop over timesteps

  pout() << "##############################################################################" << endl;
  pout() << "total number of points updated = " << m_pointsUpdated << endl;
  pout() << "" << endl;

  conclude();

  return m_time;
}
/******************/
void
EBAMRNoSubcycle::
conclude()
{
  CH_TIME("EBAMRNoSubcycle::conclude");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::conclude" << endl;
    }
#ifdef CH_USE_HDF5
  if (m_params.m_plotInterval >= 0)
    {
      writePlotFile();
    }

  if (m_params.m_checkpointInterval >= 0)
    {
      writeCheckpointFile();
    }
#endif
}
/*****************/
void
EBAMRNoSubcycle::tagCells(Vector<IntVectSet>& a_tags)
{
  CH_TIME("EBAMRNoSubcycle::tagCells");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::tagCells" << endl;
    }
  int numlevels = a_tags.size();
  if (numlevels > m_finestLevel+1) numlevels = m_finestLevel+1;

  for (int lev=0; lev<numlevels; lev++)
    {
      IntVectSet& levelTags = a_tags[lev];

      tagCellsLevel(levelTags, lev);
#ifdef CH_MPI
      int gatherint = 0;
      if (!levelTags.isEmpty()) gatherint = 1;
      int itags;
      MPI_Allreduce(&gatherint, &itags, 1, MPI_INT,
                    MPI_MAX, Chombo_MPI::comm);
      bool thereAreTags = (itags==1);
      if (!thereAreTags)
        {
          MayDay::Error("EBAMRNoSubcycle::tagCells -- numlevels > 1 and no cells tagged (probably all reg geometry with no initial vorticity)");
        }
#endif
    }
}
/*****************/
void
EBAMRNoSubcycle::regrid()
{
  CH_TIME("EBAMRNoSubcycle::regrid");
  //smoothing to make regridding less heinous
  preRegrid();

  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::regrid" << endl;
    }
  //don't regrid base grid
  int lbase = 0;

  if (m_params.m_maxLevel > 0)
    {
      //first, construct tags
      int top_level = Min(m_finestLevel, m_params.m_maxLevel-1);
      Vector<IntVectSet> tagsVect(top_level+1);
      Vector<Vector<Box> > new_grids;
      Vector<Vector<Box> > vectBoxes(top_level+1);
      for (int ilev = 0; ilev <= top_level; ilev++)
        {
//           vectBoxes[ilev] = Vector<Box>(1, m_domain[ilev].domainBox());
          domainSplit(m_domain[ilev], vectBoxes[ilev], m_params.m_maxBoxSize,
                      m_params.m_blockFactor);
        }

      // do tagging
      tagCells(tagsVect);

      int new_finest_level;

      BRMeshRefine meshrefine(m_domain[0], m_params.m_refRatio,
                              m_params.m_fillRatio, m_params.m_blockFactor,
                              m_params.m_nestingRadius, m_params.m_maxBoxSize);




      new_finest_level = meshrefine.regrid(new_grids,
                                           tagsVect,
                                           lbase,
                                           top_level,
                                           vectBoxes);



      //can only add one level at a time
      new_finest_level = Min(m_finestLevel+1, new_finest_level);

      if ((m_finestLevel != new_finest_level) && (m_params.m_verbosity >= 2))
        {
          pout() << "finest level changes here from "
                 << m_finestLevel << " to "
                 << new_finest_level << endl;
        }

      //allow for levels to change
      m_finestLevel = Min(new_finest_level, m_params.m_maxLevel);


      //now redefine grid hierarchy
      regrid(new_grids);

      defineIrregularData();

      //finish up
      postRegrid();
    } //end if max level > 0
}
/*********************/
void
EBAMRNoSubcycle::preRegrid()
{
  CH_TIME("EBAMRNoSubcycle::preRegrid");
  if (m_viscousCalc && m_params.m_doRegridSmoothing)
    {
      if (m_params.m_verbosity > 0)
        {
          pout() << " smoothing velocity before regridding " << endl;
        }
      //allocate bunch of scratch stuff
      allocateTemporaries();
      allocateExtraTemporaries();
      //make cellscratch2 == 0 for residual calc
      EBAMRDataOps::setVal(m_cellScratc2, 0.0);
      EBAMRDataOps::setVal(m_cellScratc1, 0.0);

      //make vel = (I-mu lapl)vel
      Real alpha = 1.0;
      Real beta = -4.0*m_dt;
      int coarsestLevel = 0;
      for (int velComp = 0; velComp < SpaceDim; velComp++)
        {
          DirichletPoissonEBBC::s_velComp = velComp;
          //copy velo into cellscratc1, make cellscratch2 == 0
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              Interval srcInterv(velComp, velComp);
              Interval dstInterv(0 , 0);
              m_velo[ilev]->copyTo(srcInterv, *m_cellScratch[ilev], dstInterv);
            }

          if (m_params.m_orderTimeIntegration == 2)
            {
              m_tgaSolver[velComp]->setTime(m_time);
              m_tgaSolver[velComp]->resetAlphaAndBeta(alpha, beta);
            }
          else if (m_params.m_orderTimeIntegration == 1)
            {
              m_backwardSolver[velComp]->resetAlphaAndBeta(alpha, beta);
            }
          else
            {
              MayDay::Error("EBAMRNoSubcycle::preRegrid -- bad order time integration");
            }

          //apply the operator (by computing the residual with rhs = 0, and * -1
          Real junkNorm = 0;
          junkNorm = m_solver[velComp]->computeAMRResidual(m_cellScratc1, //comes out holding -(alpha + beta lapl) vel
                                                           m_cellScratch, //holds velocity component
                                                           m_cellScratc2, //holds zero
                                                           m_finestLevel,
                                                           coarsestLevel,
                                                           false);        //not homogeneous bcs


          //make cellscratch hold (alpha + beta lapl) vel
          EBAMRDataOps::scale(m_cellScratc1,-1.0);
          averageDown(m_cellScratc1);

          //copy into velocity containers (makes vel := (alpha + beta lapl)vel)
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              Interval dstInterv(velComp, velComp);
              Interval srcInterv(0 , 0);
              m_cellScratc1[ilev]->copyTo(srcInterv, *m_velo[ilev], dstInterv);
            }
        }//end loop over velocity components
      //remove all that scratch space
      deleteTemporaries();
      deleteExtraTemporaries();
    }
}
/*********************/
void
EBAMRNoSubcycle::postRegrid()
{
  CH_TIME("EBAMRNoSubcycle::postRegrid");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::postRegrid" << endl;
    }
  if (m_viscousCalc && m_params.m_doRegridSmoothing)
    {
      if (m_params.m_verbosity > 0)
        {
          pout() << " do smoothing velocity after regridding " << endl;
        }
      //allocate bunch of scratch stuff
      allocateTemporaries();
      allocateExtraTemporaries();

      //make vel = (I-mu lapl)^-1 vel
      Real alpha = 1.0;
      //beta != -4.0*m_dt*m_viscosity (no viscosity)
      Real beta = -4.0*m_dt;
      for (int velComp = 0; velComp < SpaceDim; velComp++)
        {
          DirichletPoissonEBBC::s_velComp = velComp;
          //copy velo into cellscratc2 for rhs of tga solve
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              Interval srcInterv(velComp, velComp);
              Interval dstInterv(0 , 0);
              m_velo[ilev]->copyTo(srcInterv, *m_cellScratc2[ilev], dstInterv);
            }

          if (m_params.m_orderTimeIntegration == 2)
            {
              m_tgaSolver[velComp]->setTime(m_time);
              m_tgaSolver[velComp]->resetAlphaAndBeta(alpha, beta);
            }
          else if (m_params.m_orderTimeIntegration == 1)
            {
              m_backwardSolver[velComp]->resetAlphaAndBeta(alpha, beta);
            }
          else
            {
              MayDay::Error("EBAMRNoSubcycle::postRegrid -- bad order time integration");
            }

          EBAMRDataOps::setToZero(m_cellScratch);
          //solve (I - mu lapl)velnew = velo
          m_solver[velComp]->solve(m_cellScratch, //comes out holding desmoothed vel
                                   m_cellScratc2, //rhs = vel
                                   m_finestLevel, 0);

          averageDown(m_cellScratch);

          //copy into velocity containers (makes vel := (alpha + beta lapl)^-1 vel)
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              Interval dstInterv(velComp, velComp);
              Interval srcInterv(0 , 0);
              m_cellScratch[ilev]->copyTo(srcInterv, *m_velo[ilev], dstInterv);
            }
        }//end loop over velocity components

      //remove all that scratch space
      deleteTemporaries();
      deleteExtraTemporaries();
    }

  if (m_params.m_verbosity > 0)
    {
      pout() << "EBAMRNoSubcycle:cc projecting after regrid" << endl;
    }

  //in case gphi is not zero, need to save it and copy it back after projection
  Interval interv(0, SpaceDim-1);
  Vector<LevelData<EBCellFAB>* > tempLDPtr2;
  tempLDPtr2.resize(m_finestLevel+1);
  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      tempLDPtr2[ilev] = new LevelData<EBCellFAB>();
      tempLDPtr2[ilev]->define(m_grids[ilev], SpaceDim,  IntVect::Zero, ebcellfact);
      m_gphi[ilev]->copyTo(interv, *(tempLDPtr2[ilev]), interv);
    }

  m_ccProjector->project(m_velo, m_gphi);

  filter(m_velo);

  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      tempLDPtr2[ilev]->copyTo(interv, *m_gphi[ilev], interv);
    }

  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      delete tempLDPtr2[ilev];
    }

  //re-initialize pressure here after projecting velocity (as done in init iterations)
  //do this at thy peril (use, instead, pre-regrid data as initial guess)
  if (m_params.m_verbosity > 1)
    {
      pout() << "EBAMRNoSubcycle: iterating to get gradp after regrid" << endl;
    }
  m_advanceGphiOnly = true;
  for (int iter = 0; iter < m_params.m_gphiIterations; iter++)
    {
      pout() << "EBAMRNoSubcycle: gphi iteration " << iter  << endl;
      advance();
      postTimeStep();
    }
  m_advanceGphiOnly = false;
}
/*********************/
void EBAMRNoSubcycle::
averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data)
{
  CH_TIME("EBAMRNoSubcycle::averageDown(fluxData)");
  //do average down here
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if (ncomp == SpaceDim)
        {
          avePtr = m_aveSpac[ilev];
        }
      EBCoarseAverage& ebaverage = *avePtr;

      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } //end loop over levels
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      a_data[ilev]->exchange(interval);
    }
}
/**************************/
void
EBAMRNoSubcycle::
averageDown(Vector<LevelData<EBCellFAB>* >&  a_data)
{
  CH_TIME("EBAMRNoSubcycle::averageDown(celldata)");
  //do average down here
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if (ncomp == SpaceDim)
        {
          avePtr = m_aveSpac[ilev];
        }
      EBCoarseAverage& ebaverage = *avePtr;

      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } //end loop over levels
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      a_data[ilev]->exchange(interval);
    }
}
/*********************/
void
EBAMRNoSubcycle::defineGrids(const Vector<Vector<Box> >& a_vectBoxes)
{
  CH_TIME("EBAMRNoSubcycle::defineGrids");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::defineGrids" << endl;
    }
  m_finestLevel = 0;
  //now define all of the storage we need
  int start = 0;
  int end = a_vectBoxes.size()-1;
  for (int ilev = start; ilev<=end; ilev++)
    {
      if (a_vectBoxes[ilev].size() > 0)
        {
          m_finestLevel = ilev;
          //first do load balance
          Vector<int> procAssign;
          mortonOrdering((Vector<Box>&)(a_vectBoxes[ilev]));
          EBEllipticLoadBalance(procAssign,  a_vectBoxes[ilev], m_domain[ilev],
            false, m_ebisPtr );
          m_grids[ilev] = DisjointBoxLayout();
          m_grids[ilev].define(a_vectBoxes[ilev], procAssign);
        }
    }
}
/*********************/
void
EBAMRNoSubcycle::defineEBISLs()
{
  CH_TIME("EBAMRNoSubcycle::defineEBISLs");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::defineEBISLs" << endl;
    }
  int numEBGhost = 6;//number of ghost cells in EBISL
  m_eblg.resize(m_finestLevel+1);
  //now define all of the storage we need
  RefCountedPtr<EBPhysIBCFactory> advectBC = m_ibc->getVelAdvectBC(0); //this gets reset
  RefCountedPtr<EBPatchAdvectFactory> fact = RefCountedPtr<EBPatchAdvectFactory> (new EBPatchAdvectFactory(advectBC, m_params.m_useLimiting));

  for (int ilev = 0; ilev<= m_finestLevel; ilev++)
    {
      m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_domain[ilev], numEBGhost, m_ebisPtr);
      m_ebisl[ilev] = m_eblg[ilev].getEBISL();
      DisjointBoxLayout coarDBL;
      EBISLayout        coarEBISL;
      int refRat = 2;
       if (ilev > 0)
         {
           coarDBL =  m_grids[ilev-1];
           coarEBISL= m_ebisl[ilev-1];
           refRat = m_params.m_refRatio[ilev-1];
         }
       bool hasCoarser = (ilev > 0);
       bool hasFiner   = (ilev < m_finestLevel);
       m_ebLevAd[ilev]  = RefCountedPtr<EBLevelAdvect>(new EBLevelAdvect(m_grids[ilev],
                                                                         coarDBL,
                                                                         m_ebisl[ilev],
                                                                         coarEBISL,
                                                                         ProblemDomain(m_domain[ilev]),
                                                                         refRat,
                                                                         m_dx[ilev]*RealVect::Unit,
                                                                         hasCoarser,
                                                                         hasFiner,
                                                                         &(*fact),
                                                                         m_ebisPtr));

       /**/
       /**/
      if (ilev > 0)
        {
          //always one component for quadcfi---only way to get reuse
          int nvarquad = 1;
          m_quadCFI[ilev]  = RefCountedPtr<EBQuadCFInterp>(new  EBQuadCFInterp(m_grids[ilev], m_grids[ilev-1],
                                                                               m_ebisl[ilev], m_ebisl[ilev-1],
                                                                               m_domain[ilev-1],
                                                                               m_params.m_refRatio[ilev-1], nvarquad,
                                                                               *m_eblg[ilev].getCFIVS(),
                                                                               m_ebisPtr ));

          m_aveOper[ilev]  = RefCountedPtr<EBCoarseAverage>(new EBCoarseAverage(m_grids[ilev], m_grids[ilev-1],
                                                                                m_ebisl[ilev], m_ebisl[ilev-1],
                                                                                m_domain[ilev-1],
                                                                                m_params.m_refRatio[ilev-1], nvarquad,
                                                                                m_ebisPtr));

          m_aveSpac[ilev]  = RefCountedPtr<EBCoarseAverage>(new EBCoarseAverage(m_grids[ilev], m_grids[ilev-1],
                                                                                m_ebisl[ilev], m_ebisl[ilev-1],
                                                                                m_domain[ilev-1],
                                                                                m_params.m_refRatio[ilev-1], SpaceDim,
                                                                                m_ebisPtr));
        }
    }
  for (int ilev = 0; ilev<= m_finestLevel; ilev++)
    {
      if (ilev < m_finestLevel)
        {
          m_fluxReg[ilev]  = RefCountedPtr<EBFastFR>(new EBFastFR(m_eblg[ilev+1], m_eblg[ilev], m_params.m_refRatio[ilev], SpaceDim));
        }
    }
  long long totalPoints = 0;
  long long totalBoxes  = 0;
  int numLevels = m_finestLevel + 1;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      long long pointsThisLevel = 0;
      for (LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
        {
          pointsThisLevel += m_grids[ilev][lit()].numPts();
        }
      totalPoints += pointsThisLevel;
      totalBoxes += m_grids[ilev].size();
      pout() << "getAllIrregRefineLayouts:level[" << ilev
             << "], number of boxes = " << m_grids[ilev].size()
             << ", number of points = " << pointsThisLevel << endl;
    }
  pout() << "getAllIrregRefineLayouts:"
         <<  "   total boxes = " << totalBoxes
         <<  ", total points = " << totalPoints <<  endl;
}
/*********************/
void
EBAMRNoSubcycle::initialGrid(const Vector<Vector<Box> >& a_vectBoxes)
{
  CH_TIME("EBAMRNoSubcycle::initialGrid");
  if (m_params.m_verbosity >= 3)
    {
      pout () << "EBAMRNoSubcycle::initialGrid "  << endl;
    }
  defineGrids(a_vectBoxes);
  defineEBISLs();
  defineExtraEBISLs();
  defineNewVel();
  definePressure();
  defineExtraTerms();
  defineProjections();
}
/*********************/
void
EBAMRNoSubcycle::
defineNewVel(const int a_startLevel)
{
  CH_TIME("EBAMRNoSubcycle::defineNewVel");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::defineNewVel" << endl;
    }
  int startLevel = 0;
  if (a_startLevel > startLevel)
    {
      startLevel=a_startLevel;
    }
  for (int ilev = startLevel; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_velo[ilev]->define(m_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact);
      EBFluxFactory ebfluxfact(m_ebisl[ilev]);
      m_advVel[ilev]->define(m_grids[ilev], 1,  3*IntVect::Unit, ebfluxfact);
    }
  EBAMRDataOps::setToZero(m_velo);
  EBAMRDataOps::setToZero(m_advVel);
}
/*********************/
void
EBAMRNoSubcycle::
definePressure(const int a_startLevel)
{
  CH_TIME("EBAMRNoSubcycle::definePressure");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::definePressure" << endl;
    }
  int startLevel = 0;
  if (a_startLevel > startLevel)
    {
      startLevel=a_startLevel;
    }
  for (int ilev = startLevel; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_gphi[ilev]->define(m_grids[ilev], SpaceDim,    IntVect::Zero, ebcellfact);
      m_pres[ilev]->define(m_grids[ilev],        1,    IntVect::Zero, ebcellfact);
    }
  EBAMRDataOps::setToZero(m_gphi);
  EBAMRDataOps::setToZero(m_pres);
}
/*********************/
void
EBAMRNoSubcycle::
defineProjections()
{
  CH_TIME("EBAMRNoSubcycle::defineProjections");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::defineProjections" << endl;
    }
  if (m_ccProjector != NULL)
    {
      delete m_ccProjector;
      m_ccProjector = NULL;
    }
  if (m_macProjector != NULL)
    {
      delete m_macProjector;
      m_macProjector = NULL;
    }
  RefCountedPtr<BaseDomainBCFactory> macBCVel = m_ibc->getMACVelBC();
  RefCountedPtr<BaseDomainBCFactory> celBCPhi = m_ibc->getPressBC();
  RefCountedPtr<BaseEBBCFactory>     ebbcVelo = m_ibc->getVelocityEBBC(0);
  RefCountedPtr<BaseEBBCFactory>     ebbcPhi  = m_ibc->getPressureEBBC();

  int numLevels = m_finestLevel + 1;
  RealVect coarDxVect = m_dx[0]*RealVect::Unit;

  ParmParse pp;

  pp.query("mg_num_smooths", m_params.m_numSmooth);

  pp.query("mg_relax_type", m_params.m_relaxType);
  bool lazy = false;
  pp.query("mg_relax_lazy", lazy);

  int bottomSolverType = 0;
  pp.query("mg_bottom_solver", bottomSolverType);

  int numPreCond = 4;
  pp.query("mg_num_precond", numPreCond);

  pp.query("mg_hang", m_params.m_hang);
  pp.query("mg_tolerance", m_params.m_tolerance);

  EBAMRPoissonOp::doLazyRelax(lazy);

  int mgCoarsenLimit = 2;
  pp.query("mg_coarsen_limit", mgCoarsenLimit);
  EBAMRPoissonOpFactory::setTestRef(mgCoarsenLimit);
  int whichReflux = 0;
  EBAMRPoissonOpFactory::setWhichReflux(whichReflux);

  m_macProjector =   new EBCompositeMACProjector(m_eblg, m_params.m_refRatio, m_quadCFI,
                                                 coarDxVect, RealVect::Zero,
                                                 macBCVel, celBCPhi, ebbcPhi,
                                                 m_params.m_subtractOffMean, numLevels,
                                                 m_params.m_verbosity,numPreCond,m_time,m_params.m_relaxType,bottomSolverType);


  m_macProjector->setSolverParams(m_params.m_numSmooth, m_params.m_iterMax, m_params.m_mgCycle,
                                  m_params.m_hang, m_params.m_tolerance,
                                  m_params.m_verbosity);

  m_ccProjector  =   new EBCompositeCCProjector( m_eblg, m_params.m_refRatio, m_quadCFI,
                                                 coarDxVect, RealVect::Zero,
                                                 macBCVel, celBCPhi, ebbcPhi,
                                                 m_params.m_subtractOffMean, numLevels,
                                                 m_params.m_verbosity, numPreCond, m_time, m_params.m_relaxType,
                                                 bottomSolverType,m_macProjector);
  if (m_viscousCalc)
    {
      if (m_params.m_orderTimeIntegration == 2)
        {
          pout() << "using TGA for viscous solver" << endl;
        }
      else if (m_params.m_orderTimeIntegration == 1)
        {
          pout() << "using backward Euler for viscous solver" << endl;
        }
      else
        {
          MayDay::Error("EBAMRNoSubcycle::defineProjections -- bad order time integration");
        }

      int lbase = 0;
      int lmax = m_finestLevel;

      LinearSolver<LevelData<EBCellFAB> >* bottomSolverPtr = NULL;
      if (bottomSolverType == 0)//BiCGStab
        {
          bottomSolverPtr = &m_bottomSolver;
        }
      else if (bottomSolverType == 1)//simple relaxation
        {
          bottomSolverPtr = &m_bottomSolverSimp;
        }
      else
        {
          MayDay::Error("EBAMRNoSubcycle::defineProjections (1) -- bad bottom solver type");
        }

      Vector<LevelData<EBCellFAB>* > phi(numLevels);
      Vector<LevelData<EBCellFAB>* > rhs(numLevels);
      for (int ilev = 0; ilev < numLevels; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          phi[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1,   3*IntVect::Unit, ebcellfact);
          rhs[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1,   3*IntVect::Unit, ebcellfact);
        }

      m_bottomSolver.m_verbosity = m_params.m_verbosity;
      int numBotSmooth = 64;
      pp.query("mg_num_bottom_smooths", numBotSmooth);
      m_params.m_numBotSmooth = numBotSmooth;
      m_bottomSolverSimp.setNumSmooths(m_params.m_numBotSmooth*m_params.m_numSmooth);
      pp.query("mg_num_cycles", m_params.m_mgCycle);
      pp.query("mg_iter_max", m_params.m_iterMax);
      pp.query("mg_norm_thresh", m_params.m_normThresh);
      int numLevels = m_finestLevel + 1;

      Real alpha = 1.0;
      Real beta = m_viscosity;

      ProblemDomain level0Dom = m_eblg[0].getDomain();

      for (int idir = 0;  idir < SpaceDim; idir++)
        {
          DirichletPoissonEBBC::s_velComp = idir;
          RefCountedPtr<BaseDomainBCFactory> celBCVel = m_ibc->getVelBC(idir);

          EBAMRPoissonOpFactory opFactory(m_eblg,
                                          m_params.m_refRatio,
                                          m_quadCFI,
                                          coarDxVect,
                                          RealVect::Zero,
                                          numPreCond,
                                          m_params.m_relaxType,
                                          celBCVel,
                                          ebbcVelo,
                                          alpha,
                                          beta,
                                          m_time,
                                          3*IntVect::Unit,
                                          3*IntVect::Unit);

          m_solver[idir] = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >);
          m_solver[idir]->define(level0Dom,
                                 opFactory,
                                 bottomSolverPtr,
                                 numLevels);

          m_solver[idir]->m_pre        =  m_params.m_numSmooth;
          m_solver[idir]->m_post       =  m_params.m_numSmooth;
          m_solver[idir]->m_bottom     =  m_params.m_numSmooth;
          m_solver[idir]->m_hang       =  m_params.m_hang;
          m_solver[idir]->m_eps        =  m_params.m_tolerance;
          m_solver[idir]->m_verbosity  =  m_params.m_verbosity;
          m_solver[idir]->m_iterMax    =  m_params.m_iterMax;
          m_solver[idir]->m_normThresh =  m_params.m_normThresh;
          m_solver[idir]->setMGCycle(m_params.m_mgCycle);

          m_solver[idir]->init(phi,
                               rhs,
                               lmax,
                               lbase);

          if (m_params.m_orderTimeIntegration == 2)
            {
              m_tgaSolver[idir] = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
                (new AMRTGA<LevelData<EBCellFAB> >(m_solver[idir],
                                                   opFactory,
                                                   level0Dom,
                                                   m_params.m_refRatio,
                                                   numLevels,
                                                   m_params.m_verbosity));
            }
          else if (m_params.m_orderTimeIntegration == 1)
            {
              m_backwardSolver[idir] = RefCountedPtr<EBBackwardEuler>
                (new EBBackwardEuler(m_solver[idir],
                                     opFactory,
                                     level0Dom,
                                     m_params.m_refRatio,
                                     numLevels,
                                     m_params.m_verbosity));
            }
          else
            {
              MayDay::Error("EBAMRNoSubcycle::defineProjections (2) -- bad order time integration");
            }

        }
      for (int ilev = 0; ilev < numLevels; ilev++)
        {
          delete phi[ilev];
          delete rhs[ilev];
        }
    }
  defineExtraSolvers();
}
/*********************/
void
EBAMRNoSubcycle::initialData()
{
  CH_TIME("EBAMRNoSubcycle::initialData");
  if (m_params.m_verbosity >= 3)
    {
      pout () << "EBAMRNoSubcycle::initialData "  << endl;
    }

  for (int ilev = 0; ilev<= m_finestLevel; ilev++)
    {
      RealVect dxLev = m_dx[ilev]*RealVect::Unit;
      m_ibc->initializeVelocity(*m_velo[ilev], m_grids[ilev], m_ebisl[ilev], m_domain[ilev], RealVect::Zero, m_time, dxLev);
      EBLevelDataOps::setVal(   *m_gphi[ilev], 0.0);
      EBLevelDataOps::setVal(   *m_pres[ilev], 0.0);
    }

  initialExtraData();
}
/*********************/
void
EBAMRNoSubcycle::postInitialize()
{
  CH_TIME("EBAMRNoSubcycle::postInitialize");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::postInitialize" << endl;
    }
}
/*********************/
void
EBAMRNoSubcycle::useFixedDt(Real a_dt)
{
  CH_TIME("EBAMRNoSubcycle::useFixedDt");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::useFixedDt" << endl;
    }
  m_dt = a_dt;
  m_prescribedDt = a_dt;
  m_useFixedDt = true;
}
/*********************/
Real
EBAMRNoSubcycle::computeDt()
{

  CH_TIMERS("EBAMRNoSubcycle::computeDt");
  CH_TIMER("computeDt_fillMask",t1);
  CH_TIMER("computeDt_allReduce",t2);
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::computeDt " << endl;
    }

  Real dt;
  if (m_useFixedDt)
    {
      dt = m_prescribedDt;
    }
  else
    {
      int numLevels = m_finestLevel+1;
      Real maxVel = 0.0;
      for (int ilev=0; ilev < numLevels; ilev++)
        {
          LevelData<EBCellFAB>& velNewLD = *m_velo[ilev];

          for (DataIterator dit = velNewLD.dataIterator(); dit.ok(); ++dit)
            {
              BaseFab<Real>& velFAB = (velNewLD[dit()]).getSingleValuedFAB();
              int iRegIrregCovered;
              CH_START(t1);
              const BaseFab<int>& maskFAB = m_ebisl[ilev][dit()].getEBGraph().getMask(iRegIrregCovered);
              CH_STOP(t1);
              if (iRegIrregCovered != -1)//not all covered
                {
                  if (iRegIrregCovered == 0)//has irreg
                    {
                      Box box = m_grids[ilev].get(dit());
                      for (int idir = 0; idir < SpaceDim; idir++)
                        {
                          FORT_MAXNORMMASK(CHF_REAL(maxVel),
                                           CHF_CONST_FRA1(velFAB,idir),
                                           CHF_BOX(box),
                                           CHF_CONST_FIA1(maskFAB,0));
                        }
                      IntVectSet ivs = m_ebisl[ilev][dit()].getMultiCells(box);
                      for (VoFIterator vofit(ivs, m_ebisl[ilev][dit()].getEBGraph()); vofit.ok(); ++vofit)
                        {
                          for (int idir = 0; idir < SpaceDim; idir++)
                            {
                              if (Abs(velNewLD[dit()](vofit(), idir)) > maxVel)
                              {
                                maxVel = Abs(velNewLD[dit()](vofit(), idir));
                              }
                            }
                        }
                    }
                  else//all reg
                    {
                      Box box = m_grids[ilev].get(dit());
                      for (int idir = 0; idir < SpaceDim; idir++)
                        {
                          FORT_MAXNORM(CHF_REAL(maxVel),
                                       CHF_CONST_FRA1(velFAB,idir),
                                       CHF_BOX(box));
                        }
                    }
                }
            }
        }
      CH_START(t2);
#ifdef CH_MPI
      Real tmp = 1.;
      int result = MPI_Allreduce(&maxVel, &tmp, 1, MPI_CH_REAL,
                                 MPI_MAX, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
        {
          MayDay::Error("communication error on norm");
        }
      maxVel = tmp;
#endif
      CH_STOP(t2);

      if (maxVel > 1.0e-10)
        {
          dt = m_params.m_cfl*m_dx[m_finestLevel]/maxVel;
        }
      else
        {
          dt = m_params.m_maxDt;
        }
    } //end if no prescribed dt

  //finally, enforce max dt grow factor
  if ( (m_dt > 0.) && (m_params.m_maxDtGrow > 0) &&
       (dt > m_params.m_maxDtGrow*m_dt) )
    {
      dt = m_params.m_maxDtGrow*m_dt;
    }

  return dt;
}
/*********************/
Real
EBAMRNoSubcycle::
computeInitialDt()
{
  CH_TIME("EBAMRNoSubcycle::computeInitialDt");
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRNoSubcycle::computeInitialDt " << endl;
    }
  Real retval;
  if (m_useFixedDt)
    {
      retval = m_prescribedDt;
    }
  else
    {
      retval =  (m_params.m_initCFL/m_params.m_cfl)*computeDt();
    }
  return retval;
}
/*********************/
void
EBAMRNoSubcycle::
regrid(const Vector<Vector<Box> >& a_newGrids)
{
  CH_TIME("EBAMRNoSubcycle::regrid");

  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRNoSubcycle::regrid " << endl;
    }

  Interval interv(0, SpaceDim-1);
  Vector<LevelData<EBCellFAB>* > tempLDPtr, tempLDPtr2;
  tempLDPtr.resize(m_finestLevel+1);
  tempLDPtr2.resize(m_finestLevel+1);
  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      tempLDPtr[ilev] = new LevelData<EBCellFAB>();
      tempLDPtr[ilev]->define(m_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact);
      m_velo[ilev]->copyTo(interv, *(tempLDPtr[ilev]), interv);
      tempLDPtr2[ilev] = new LevelData<EBCellFAB>();
      tempLDPtr2[ilev]->define(m_grids[ilev], SpaceDim,  IntVect::Zero, ebcellfact);
      m_gphi[ilev]->copyTo(interv, *(tempLDPtr2[ilev]), interv);
    }

  cacheExtraTerms();

  //this changes m_grids and m_ebisl
  defineGrids(a_newGrids);
  defineEBISLs();
  defineExtraEBISLs();

  //this redefines new data with new set of grids
  //this can also change m_finestLevel
  defineNewVel();
  definePressure();
  defineProjections();
  defineExtraTerms();
  //now fill new data holders, copying from old over old grids and interpolating where there is new grid
  tempLDPtr[0]->copyTo(interv, *m_velo[0], interv);
  tempLDPtr2[0]->copyTo(interv, *m_gphi[0], interv);
  for (int ilev=1; ilev<= m_finestLevel; ilev++)
    {
      //interpolate everywhere
      EBPWLFineInterp ebInterpVec(m_grids[ ilev  ],
                                  m_grids[ ilev-1],
                                  m_ebisl[ ilev  ],
                                  m_ebisl[ ilev-1],
                                  m_domain[ilev-1],
                                  m_params.m_refRatio[ilev-1],
                                  SpaceDim,
                                  m_ebisPtr);

      ebInterpVec.interpolate(*m_velo[ilev  ],
                              *m_velo[ilev-1],
                              interv);
      tempLDPtr[ilev]->copyTo(interv, *m_velo[ilev], interv);

      ebInterpVec.interpolate(*m_gphi[ilev  ],
                              *m_gphi[ilev-1],
                              interv);
      tempLDPtr2[ilev]->copyTo(interv, *m_gphi[ilev], interv);
    } //end loop over levels

  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      delete tempLDPtr[ilev];
      delete tempLDPtr2[ilev];
    }

  interpolateExtraTerms();
}
/*****************/
void
EBAMRNoSubcycle::postTimeStep()
{
  CH_TIME("EBAMRNoSubcycle::postTimeStep");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::postTimeStep" << endl;
    }
  //do average down here
  averageDown(m_velo);
  averageDown(m_pres);
  averageDown(m_gphi);
  averageDownExtraTerms();
}
/*****************/
void
EBAMRNoSubcycle::
computeVorticity(LevelData<EBCellFAB>& a_vort, int a_level)
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::computeVorticity" << endl;
    }
  CH_assert(m_isSetup);
  EBCellFactory ebcellfact(m_ebisl[a_level]);
  //define the data holder
  //vorticity is a vector in 3d, scalar in 2d
  int ncomp = 1;
  if (SpaceDim==3)
    {
      ncomp = SpaceDim;
    }
  a_vort.define(m_grids[a_level], ncomp, IntVect::Zero, ebcellfact);

  //interpolate velocity at coarse-fine interfaces
  //need some scratch space to do quadcfi
  for (int ilev=0; ilev < m_cellScratch.size(); ilev++)
    {
      if (m_cellScratch[ilev] != NULL)
        {
          delete m_cellScratch[ilev];
          m_cellScratch[ilev] = NULL;
        }
    }

  int numLevels = m_finestLevel+1;
  m_cellScratch.resize(numLevels, NULL);
  for (int ilev=0; ilev < numLevels; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_cellScratch[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 3*IntVect::Unit, ebcellfact);
    }

  if (a_level > 0)
    {

      LevelData<EBCellFAB>& phiCoar = *m_cellScratch[a_level-1];
      LevelData<EBCellFAB>& phiFine = *m_cellScratch[a_level  ];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Interval phiInterv(0, 0);
          Interval velInterv(idir, idir);
          m_velo[a_level-1]->copyTo(velInterv, phiCoar, phiInterv);
          m_velo[a_level  ]->copyTo(velInterv, phiFine, phiInterv);
          m_quadCFI[a_level]->interpolate(phiFine,
                                          phiCoar,
                                          phiInterv);

          //on copy back, we need ghost cells, so do the data iterator loop
          for (DataIterator dit = phiFine.dataIterator(); dit.ok(); ++dit)
            {
              Box region = phiFine[dit()].getRegion(); //includes ghost cells
              (*m_velo[a_level])[dit()].copy(region, velInterv, region, phiFine[dit()], phiInterv);
            }

          EBLevelDataOps::setVal(phiFine, 0.0);
          EBLevelDataOps::setVal(phiCoar, 0.0);
        }

    }
  m_velo[a_level]->exchange(Interval(0, SpaceDim-1));

  //clean up temp space
  for (int ilev=0; ilev < m_cellScratch.size(); ilev++)
    {
      if (m_cellScratch[ilev] != NULL)
        {
          delete m_cellScratch[ilev];
          m_cellScratch[ilev] = NULL;
        }
    }

  Real dxLev = m_dx[a_level];
  //do actual computation
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&  vortFAB =                a_vort[dit()];
      EBCellFAB&  veloFAB = (*m_velo[a_level])[dit()];
      const EBGraph& ebgraph =   m_ebisl[a_level][dit()].getEBGraph();
      vortFAB.setVal(0.);

      BaseFab<Real>& regVort = vortFAB.getSingleValuedFAB();
      BaseFab<Real>& regVelo = veloFAB.getSingleValuedFAB();
      Box interiorBox = m_grids[a_level].get(dit());
      interiorBox.grow(1);
      interiorBox &= m_domain[a_level];
      interiorBox.grow(-1);
      int vortIndex;
#if CH_SPACEDIM==3
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          vortIndex = idir;
#else
          vortIndex = 0;
          int idir = 3;
#endif

          //compute on all cells as regular
          FORT_COMPUTEVORT(CHF_FRA1(regVort, vortIndex),
                           CHF_CONST_FRA(regVelo),
                           CHF_BOX(interiorBox),
                           CHF_CONST_REAL(dxLev),
                           CHF_CONST_INT(idir));;

          //do cells on domain boundary as if they were irregular
          IntVectSet ivsIrreg(m_grids[a_level].get(dit()));
          ivsIrreg -= interiorBox;
          ivsIrreg |= ebgraph.getIrregCells(interiorBox);
          int diffDirVec[2];
          if (SpaceDim==2 || idir ==2)
            {
              diffDirVec[0] = 0;
              diffDirVec[1] = 1;
            }
          else if (idir == 0)
            {
              diffDirVec[0] = 1;
              diffDirVec[1] = 2;
            }
          else if (idir==1)
            {
              diffDirVec[0] = 2;
              diffDirVec[1] = 0;
            }
          else
            {
              MayDay::Error("missed a case");
            }

          for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
            {
              const VolIndex vof = vofit();
              Real vortValue = 0.0;

              //taking the derivative d(udvdir)/d(xdiffdir)
              for (int idiff = 0; idiff < 2; idiff++)
                {

                  Real signDiff = 0;
                  int diffDir= -1;
                  int velComp = -1;
                  if (idiff == 0)
                    {
                      signDiff = 1.0;
                      diffDir  = diffDirVec[0];
                      velComp  = diffDirVec[1];
                    }
                  else if (idiff == 1)
                    {
                      signDiff = -1.0;
                      diffDir  = diffDirVec[1];
                      velComp  = diffDirVec[0];
                    }
                  else
                    {
                      MayDay::Error("missed a case");
                    }

                  Vector<FaceIndex> hiFaces =  ebgraph.getFaces(vof, diffDir, Side::Hi);
                  Vector<FaceIndex> loFaces =  ebgraph.getFaces(vof, diffDir, Side::Lo);

                  bool hasHi = (hiFaces.size() == 1) && (!hiFaces[0].isBoundary());
                  bool hasLo = (loFaces.size() == 1) && (!loFaces[0].isBoundary());
                  Real diffValue = 0.0;
                  if (hasHi && hasLo)
                    {
                      const VolIndex& vofHi = hiFaces[0].getVoF(Side::Hi);
                      const VolIndex& vofLo = loFaces[0].getVoF(Side::Lo);
                      Real hiValue = veloFAB(vofHi, velComp);
                      Real loValue = veloFAB(vofLo, velComp);
                      diffValue =0.5*(hiValue - loValue)/dxLev;
                    }
                  else if (hasHi)
                    {
                      const VolIndex& vofHi = hiFaces[0].getVoF(Side::Hi);
                      Real hiValue = veloFAB(vofHi, velComp);
                      Real loValue = veloFAB(vof,   velComp);
                      diffValue =(hiValue - loValue)/dxLev;
                    }
                  else if (hasLo)
                    {
                      const VolIndex& vofLo = loFaces[0].getVoF(Side::Lo);
                      Real hiValue = veloFAB(vof  , velComp);
                      Real loValue = veloFAB(vofLo, velComp);
                      diffValue =  (hiValue - loValue)/dxLev;
                    }
                  else
                    {
                      diffValue = 0.0;
                    }
                  vortValue += signDiff*diffValue;
                } //end loop over idiff

              vortFAB(vof, vortIndex) = vortValue;

            } //end loop over irregular vofs
#if CH_SPACEDIM==3
        } //end loop over vort components in 3d
#endif

    }
}
/**********/
void
EBAMRNoSubcycle::
tagCellsLevel(IntVectSet& a_tags, int a_level)
{
  CH_TIME("EBAMRNoSubcycle::tagCellsLevel");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::tagCellsLevel" << endl;
    }
  CH_assert(m_isSetup);
  LevelData<EBCellFAB> vort;
  computeVorticity(vort, a_level);
  a_tags.makeEmpty();

  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit)
    {
      const EBCellFAB& vortFAB = vort[dit()];
      const Box& grid =        m_grids[a_level].get(dit());
      const EBGraph& ebgraph = m_ebisl[a_level][dit()].getEBGraph();
      IntVectSet ivsTot(grid);
      Box shrunkDom = m_domain[a_level].domainBox();
      int shrinkNumCells = m_params.m_tagShrinkDomain;
      int shrinkRefRatio = 1;
      for (int ilev = 1; ilev <= a_level; ilev++)
        {
          shrinkRefRatio *= m_params.m_refRatio[ilev-1];
        }
      for (int idir = 0;  idir < SpaceDim; idir++)
        {
          if (idir != m_params.m_flowDir)
            {
              shrunkDom.grow(idir, -shrinkNumCells*shrinkRefRatio);
            }
        }
      ivsTot &= shrunkDom;

      for (VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          Real vortmag = 0.0;

          for (int icomp = 0; icomp < vortFAB.nComp(); icomp++)
            {
              Real vortDirVal = vortFAB(vof, icomp);
              vortmag += vortDirVal*vortDirVal;
            }
          vortmag = sqrt(vortmag);
          if (vortmag >= m_params.m_refineThreshold)
            {
              a_tags |= iv;
            }
        } //end loop over vofs

      //refine all irregular cells
      IntVectSet irregIVS = ebgraph.getIrregCells(grid);
      a_tags |= irregIVS;
    }

  a_tags.grow(m_params.m_tagBuffer);

}
/*****************/
void
EBAMRNoSubcycle::
defineIrregularData()
{
  CH_TIME("EBAMRNoSubcycle::defineIrregularData");
  for (int ilev = 0; ilev <=  m_finestLevel; ilev++)
    {
      m_coveredFaceLitLo[ilev]->define(m_grids[ilev]);
      m_coveredFaceLitHi[ilev]->define(m_grids[ilev]);
      m_coveredSetsLitLo[ilev]->define(m_grids[ilev]);
      m_coveredSetsLitHi[ilev]->define(m_grids[ilev]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          Box litBox = m_grids[ilev].get(dit());
          litBox.grow(1);
          litBox &= m_domain[ilev];
          (*m_coveredFaceLitLo[ilev])[dit()].resize(SpaceDim);
          (*m_coveredFaceLitHi[ilev])[dit()].resize(SpaceDim);
          (*m_coveredSetsLitLo[ilev])[dit()].resize(SpaceDim);
          (*m_coveredSetsLitHi[ilev])[dit()].resize(SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {

              IntVectSet irregIVSPlus, irregIVSMinu;
              //get the covered sets and faces
              EBArith::computeCoveredFaces((*m_coveredFaceLitHi[ilev])[dit()][idir],
                                           (*m_coveredSetsLitHi[ilev])[dit()][idir],
                                           irregIVSPlus,idir, Side::Hi, m_ebisl[ilev][dit()], litBox);
              EBArith::computeCoveredFaces((*m_coveredFaceLitLo[ilev])[dit()][idir],
                                           (*m_coveredSetsLitLo[ilev])[dit()][idir],
                                           irregIVSMinu, idir, Side::Lo,  m_ebisl[ilev][dit()], litBox);

            }
        }

    }//end loop over levels

  for (int ilev = 0; ilev <=  m_finestLevel; ilev++)
    {
      m_coveredAdvVelLo[ilev]->define(m_grids[ilev]);
      m_coveredAdvVelHi[ilev]->define(m_grids[ilev]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*m_coveredAdvVelLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*m_coveredAdvVelHi[ilev])[dit()].resize(SpaceDim, NULL);
          const EBGraph& ebgraph = m_ebisl[ilev][dit()].getEBGraph();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              (*m_coveredAdvVelLo[ilev])[dit()][idir]  = new BaseIVFAB<Real>((*m_coveredSetsLitLo[ilev])[dit()][idir], ebgraph, 1);
              (*m_coveredAdvVelHi[ilev])[dit()][idir]  = new BaseIVFAB<Real>((*m_coveredSetsLitHi[ilev])[dit()][idir], ebgraph, 1);
            }
        }
    }
}
/*****************/
void
EBAMRNoSubcycle::
allocateTemporaries()
{
  CH_TIME("EBAMRNoSubcycle::allocateTemporaries");
  //number of active levels
  int numLevels = m_finestLevel+1;
  m_uDotDelU.resize(numLevels, NULL);

  m_cellScratch.resize(numLevels, NULL);
  m_cellScratc1.resize(numLevels, NULL);
  m_cellScratc2.resize(numLevels, NULL);
  m_macGradient.resize(numLevels, NULL);
  m_macScratch1.resize(numLevels, NULL);
  m_macScratch2.resize(numLevels, NULL);
  m_coveredScratchLo.resize(numLevels, NULL);
  m_coveredScratchHi.resize(numLevels, NULL);

  //allocate storage
  for (int ilev=0; ilev < numLevels; ilev++)
    {
      EBFluxFactory ebfluxfact(m_ebisl[ilev]);
      EBCellFactory ebcellfact(m_ebisl[ilev]);

      m_uDotDelU[ilev]    = new LevelData<EBCellFAB>(m_grids[ilev], SpaceDim, IntVect::Unit, ebcellfact);
      m_cellScratch[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],        1, 3*IntVect::Unit, ebcellfact);
      m_cellScratc1[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],        1, 3*IntVect::Unit, ebcellfact);
      m_cellScratc2[ilev] = new LevelData<EBCellFAB>(m_grids[ilev],        1, 3*IntVect::Unit, ebcellfact);


      m_macGradient[ilev] = new LevelData<EBFluxFAB>(m_grids[ilev],        1, 3*IntVect::Unit, ebfluxfact);
      m_macScratch1[ilev] = new LevelData<EBFluxFAB>(m_grids[ilev],        1, 3*IntVect::Unit, ebfluxfact);
      m_macScratch2[ilev] = new LevelData<EBFluxFAB>(m_grids[ilev],        1, 3*IntVect::Unit, ebfluxfact);

      m_coveredScratchLo[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      m_coveredScratchHi[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);

      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*m_coveredScratchLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*m_coveredScratchHi[ilev])[dit()].resize(SpaceDim, NULL);

          const EBGraph& ebgraph = m_ebisl[ilev][dit()].getEBGraph();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              (*m_coveredScratchLo[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitLo[ilev])[dit()][idir], ebgraph, 1);
              (*m_coveredScratchHi[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitHi[ilev])[dit()][idir], ebgraph, 1);
            }
        }
    }
  EBAMRDataOps::setToZero(m_uDotDelU  );
  EBAMRDataOps::setToZero(m_cellScratch);
  EBAMRDataOps::setToZero(m_cellScratc1);
  EBAMRDataOps::setToZero(m_cellScratc2);
}
/*****************/
void
EBAMRNoSubcycle::
deleteTemporaries()
{
  CH_TIME("EBAMRNoSubcycle::deleteTemporaries");
  //number of active levels
  int numLevels = m_finestLevel+1;
  //clean up storage
  for (int ilev=0; ilev<numLevels; ilev++)
    {
      delete m_uDotDelU[ilev];
      delete m_cellScratch[ilev];
      delete m_cellScratc1[ilev];
      delete m_cellScratc2[ilev];
      delete m_macGradient[ilev];
      delete m_macScratch1[ilev];
      delete m_macScratch2[ilev];
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*m_coveredScratchLo[ilev])[dit()][idir];
              delete (*m_coveredScratchHi[ilev])[dit()][idir];
            }
        }
      delete m_coveredScratchLo[ilev];
      delete m_coveredScratchHi[ilev];
    }

  for (int ilev=0; ilev<numLevels; ilev++)
    {
      m_uDotDelU[ilev]         = NULL;
      m_cellScratch[ilev]      = NULL;
      m_cellScratc1[ilev]      = NULL;
      m_cellScratc2[ilev]      = NULL;
      m_macGradient[ilev]      = NULL;
      m_macScratch1[ilev]      = NULL;
      m_macScratch2[ilev]      = NULL;
      m_coveredScratchLo[ilev] = NULL;
      m_coveredScratchHi[ilev] = NULL;
    }
}
/*****************/
void
EBAMRNoSubcycle::
advance()
{

  CH_TIME("EBAMRNoSubcycle::advance");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::advance, nstep " << m_curStep
             << ", starting time = "  << m_time
             << ", dt = " << m_dt << endl;
    }

  //allocate space for integration
  allocateTemporaries();
  allocateExtraTemporaries();

  predictor();
  corrector();

  pointsUpdated();

  //remove all unnecessary data
  deleteTemporaries();
  deleteExtraTemporaries();

}
/******************/
void
EBAMRNoSubcycle::
predictor()
{
  predictExtraFieldsPreVelocity();
  if (!m_steadyState)
    {
      predictVelocity();
    }
  predictExtraFieldsPostVelocity();
}
/******************/
void
EBAMRNoSubcycle::
predictVelocity()
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::predictVelocity" << endl;
    }
  CH_TIMERS("EBAMRNoSubcycle::predictVelocity");

  if (m_params.m_verbosity > 1)
    {
      pout() << "EBAMRNoSubcycle::normalVelocityPredictor" << endl;
    }
  normalVelocityPredictor(m_advVel, m_coveredAdvVelLo, m_coveredAdvVelHi, m_velo);

  if (m_params.m_verbosity > 1)
    {
      pout() << "EBAMRNoSubcycle::transverseVelocityPredictor" << endl;
    }
  transverseVelocityPredictor(m_uDotDelU, m_velo, true);

  if (m_params.m_stokesFlow)
    {
      if (m_params.m_verbosity > 1)
        {
          pout() << "EBAMRNoSubcycle::predictVelocity: <<< Stokes flow ---> udelu==0 >>>" << endl;
        }
      EBAMRDataOps::setToZero(m_uDotDelU);
    }
}

/*****************/
void
EBAMRNoSubcycle::
normalVelocityPredictor(Vector<LevelData<EBFluxFAB> *>&                       a_advVel,
                        Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
                        Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
                        const Vector<LevelData<EBCellFAB> *>&                 a_velo)
{
  CH_TIMERS("EBAMRNoSubcycle::normalVelocityPredictor");
  CH_TIMER("extrapolation_to_faces", t1);
  CH_TIMER("mac_projection", t2);
  CH_TIMER("post_projection", t3);

  CH_START(t1);
  //extrapolate velocities to faces
  //normal velocities at edges
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      //initially fill advective velocity with average to faces of
      //cell-centered solution
      RealVect dxLev = m_dx[ilev]*RealVect::Unit;
      LayoutData<IntVectSet> cfivs;
      ccpAverageVelocityToFaces(*m_macScratch1[ilev], *m_velo[ilev],
                                m_grids[ilev], m_ebisl[ilev], m_domain[ilev], dxLev,
                                cfivs);
    }
  //extrapolate to get advective velocity at covered faces
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          //copy cell centered velocity into scratch space for each direction
          Interval srcInterv(idir, idir);
          Interval dstInterv(0 , 0);
          m_velo[ilev]->copyTo(srcInterv, *m_cellScratch[ilev], dstInterv);
        }
      extrapolateToCoveredFaces(m_coveredScratchLo,
                                m_coveredScratchHi,
                                m_macScratch1,
                                m_cellScratch, idir);
    }

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          //copy cell centered velocity into scratch space for each direction
          Interval srcInterv(idir, idir);
          Interval dstInterv(0 , 0);
          EBLevelDataOps::setVal(*m_cellScratch[ilev], 0.0);
          a_velo[ilev]->copyTo(srcInterv, *m_cellScratch[ilev], dstInterv);
        }

      //fill in viscous source term where necessary
      Vector<LevelData<EBCellFAB>*>* source = NULL;
      if (m_viscousCalc)
        {
          //cellscratch already holds velocity component

          DirichletPoissonEBBC::s_velComp = idir;
          viscousSourceForAdvect(m_cellScratc2,   //returns holding source term = nu*lapl
                                 m_cellScratch,   //holds cell centered vel comp n
                                 m_cellScratc1,   //will hold the zero for the residual calc (zeroed inside routine)
                                 idir,            //velocity component
                                 m_time);         //time, for BC setting

          source = &m_cellScratc2;
        }

      Vector<LevelData<EBCellFAB>*> extraSource;
      extraSource.resize(m_finestLevel+1, NULL);
      for (int ilev=0; ilev <= m_finestLevel; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          extraSource[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 3*IntVect::Unit, ebcellfact);
        }
      EBAMRDataOps::setToZero(extraSource);

      computeExtraSourceForPredictor(extraSource, idir);

      if (source != NULL)
        {
          for (int ilev=0; ilev<= m_finestLevel; ilev++)
            {
              for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
                {
                  (*(*source)[ilev])[dit()] += (*extraSource[ilev])[dit()];
                }
            }
        }
      else
        {
          source = &extraSource;
        }

      if (m_params.m_verbosity > 3)
        {
          pout() << "EBAMRNoSubcycle::computing advection velocities, component " << idir << endl;
        }
      EBPatchGodunov::setCurComp(idir);
      EBPatchGodunov::setDoingVel(1);
      EBPatchGodunov::setDoingAdvVel(1);

      RefCountedPtr<EBPhysIBCFactory> advectBC = m_ibc->getVelAdvectBC(idir);
      //extrapolate velocity component to faces
      extrapolateScalarCol(m_macScratch2,
                           m_coveredScratchLo,
                           m_coveredScratchHi,
                           advectBC,
                           m_macScratch1,//contains initial advective vel
                           source,
                           m_cellScratch,
                           a_velo);

      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          delete extraSource[ilev];
        }

      //now copy the result to the appropriate faces in advVel. (advVel only has normal velocities)
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              Interval interv(0, 0);
              EBFaceFAB&  extrapFAB = (*m_macScratch2[ilev])[dit()][idir];
              EBFaceFAB&  advVelFAB =  (*a_advVel[ilev])[dit()][idir];
              Box region = extrapFAB.getCellRegion();
              advVelFAB.copy(region, interv, region, extrapFAB, interv);

              //same goes for covered velocites.  only normal ones get into covered adv vel.
              (*a_coveredAdvVelLo[ilev])[dit()][idir]->copy(region, interv, region, *(*m_coveredScratchLo[ilev])[dit()][idir], interv);
              (*a_coveredAdvVelHi[ilev])[dit()][idir]->copy(region, interv, region, *(*m_coveredScratchHi[ilev])[dit()][idir], interv);
            }
        }
    } //end loop over velocity directions (idir)

  CH_STOP(t1);

  //MAC project the velocity
  if (m_params.m_verbosity >= 2)
    {
      pout() << "EBAMRNoSubcycle: mac projecting advection velocity" << endl;
    }

  CH_START(t2);

  Real time = m_time + 0.5*m_dt;
  m_macProjector->setTime(time);
  EBAMRDataOps::scale(m_pres, m_dt/2.);
  m_macProjector->setInitialPhi(m_pres);
  m_macProjector->project(a_advVel, m_macGradient);

  CH_STOP(t2);

  // EBLevelMACProjector::setVerbose(false);
  if (m_params.m_verbosity > 3)
    {

      m_macProjector->kappaDivergence(m_cellScratch, a_advVel);
      averageDown(m_cellScratch);

      Real norm[3];
      for (int inorm = 0; inorm < 3; inorm++)
        {
          norm[inorm] = EBArith::norm(m_cellScratch, m_grids, m_ebisl,
                                      m_params.m_refRatio, 0, inorm, EBNormType::OverBoth);
        }
      pout() << setprecision(8)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific);
      pout() << "div(adv vel) after mac projection: " <<
        "L_inf = " << norm[0]  <<
        ", L_1 = " << norm[1] <<
        ", L_2 = " << norm[2] << endl;

    }

  CH_START(t3);
  //average down advection velocity so that it makes sense at coarse-fine interfaces
  averageDown(a_advVel     );
  averageDown(m_macGradient);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      //faceDir, velcomp are the same thing for advection velocities
      int faceDir = idir;
      int velComp = idir;
      //correct covered velocity with extrapolation of pressure gradient
      m_macProjector->correctVelocityComponent(a_coveredAdvVelLo,
                                               a_coveredAdvVelHi,
                                               m_coveredFaceLitLo,
                                               m_coveredFaceLitHi,
                                               m_coveredSetsLitLo,
                                               m_coveredSetsLitHi,
                                               m_macGradient, faceDir, velComp);
    }
  CH_STOP(t3);
}
/*****************/
void
EBAMRNoSubcycle::
transverseVelocityPredictor(Vector<LevelData<EBCellFAB>* >&    a_uDotDelU,
                            Vector<LevelData<EBCellFAB>* >&    a_scalOld,
                            bool                               a_reallyVelocity)
{
  CH_TIME("EBAMRNoSubcycle::computeAdvectiveDerivative");
  int ncomp = a_uDotDelU[0]->nComp();
  int nlevels = m_finestLevel+1;

  //make temporaries with right number of variables
  Vector<LevelData<EBFluxFAB>* >                       macScratchVec(nlevels, NULL);
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScratchVecLo(nlevels, NULL);
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScratchVecHi(nlevels, NULL);

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      EBFluxFactory ebfluxfact(m_ebisl[ilev]);
      macScratchVec[ilev]       = new LevelData<EBFluxFAB>(m_grids[ilev], ncomp, 3*IntVect::Unit, ebfluxfact);

      coveredScratchVecLo[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      coveredScratchVecHi[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {

          (*coveredScratchVecLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*coveredScratchVecHi[ilev])[dit()].resize(SpaceDim, NULL);

          const EBGraph& ebgraph = m_ebisl[ilev][dit()].getEBGraph();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              (*coveredScratchVecLo[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitLo[ilev])[dit()][idir], ebgraph, ncomp);
              (*coveredScratchVecHi[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitHi[ilev])[dit()][idir], ebgraph, ncomp);
            }
        }
    }
  for (int icomp = 0; icomp < ncomp; icomp++)
    {
      EBPatchGodunov::setCurComp(icomp);
      EBPatchGodunov::setDoingAdvVel(0);

      RefCountedPtr<EBPhysIBCFactory> advectBC;
      if (a_reallyVelocity)
        {
          advectBC  = m_ibc->getVelAdvectBC(icomp);
          EBPatchGodunov::setDoingVel(1);
        }
      else
        {
          advectBC  = m_ibc->getScalarAdvectBC(icomp);
          EBPatchGodunov::setDoingVel(0);
        }

      for (int ilev=0; ilev <= m_finestLevel; ilev++)
        {
          Interval srcInterv(icomp, icomp);
          Interval dstInterv(0, 0);
          a_scalOld[ilev]->copyTo(srcInterv, *m_cellScratch[ilev], dstInterv);
        }

      //fill in viscous source term where necessary
      Vector<LevelData<EBCellFAB>* > * source = NULL;
      if (a_reallyVelocity && m_viscousCalc)
        {
          //cellscratch already holds velocity component

          DirichletPoissonEBBC::s_velComp = icomp;
          viscousSourceForAdvect(m_cellScratc2,   //returns holding source term = nu*lapl
                                 m_cellScratch,   //holds cell centered vel comp n
                                 m_cellScratc1,   //will hold the zero for the residual calc (zeroed inside routine)
                                 icomp,           //velocity component
                                 m_time);         //time, for BC setting
          source = &m_cellScratc2;
        }

      Vector<LevelData<EBCellFAB>*> extraSource;
      extraSource.resize(m_finestLevel+1, NULL);
      for (int ilev=0; ilev <= m_finestLevel; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          extraSource[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 3*IntVect::Unit, ebcellfact);
        }
      EBAMRDataOps::setToZero(extraSource);

      computeExtraSourceForPredictor(extraSource, icomp);

      if (source != NULL)
        {
          for (int ilev=0; ilev<= m_finestLevel; ilev++)
            {
              for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
                {
                  (*(*source)[ilev])[dit()] += (*extraSource[ilev])[dit()];
                }
            }
        }
      else
        {
          source = &extraSource;
        }

      //cellscratch used for consstate
      extrapolateScalarCol(m_macScratch1,
                           m_coveredScratchLo,
                           m_coveredScratchHi,
                           advectBC,
                           m_advVel,
                           source,
                           m_cellScratch,
                           m_velo);

      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          delete extraSource[ilev];
        }

      if (a_reallyVelocity)
        {
          //correct with previously stored pressure gradient
          Real time = m_time + 0.5*m_dt;
          m_macProjector->setTime(time);
          m_macProjector->correctVelocityComponent(m_macScratch1, m_macGradient, icomp);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              int faceDir = idir;
              int velComp = icomp;
              //correct covered velocity with extrapolation of pressure gradient
              m_macProjector->correctVelocityComponent(m_coveredScratchLo,
                                                       m_coveredScratchHi,
                                                       m_coveredFaceLitLo,
                                                       m_coveredFaceLitHi,
                                                       m_coveredSetsLitLo,
                                                       m_coveredSetsLitHi,
                                                       m_macGradient, faceDir, velComp);
            }

          //overwrite extrapolated with advective velocity if appropriate
          //(the normal velocity here had the wrong boundary conditions)
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
                {
                  Box box = grow(m_grids[ilev].get(dit()), 1);
                  Interval interv(0,0);
                  //copy over covered advective velocity
                  (*m_coveredScratchLo[ilev])[dit()][icomp]->copy(box, interv, box, *(*m_coveredAdvVelLo[ilev])[dit()][icomp], interv);
                  (*m_coveredScratchHi[ilev])[dit()][icomp]->copy(box, interv, box, *(*m_coveredAdvVelHi[ilev])[dit()][icomp], interv);

                  EBFluxFAB& macExtrapFAB    = (*m_macScratch1[ilev])[dit()];
                  const EBFluxFAB& advVelFAB = (*m_advVel[ilev])[dit()];
                  //icomp is the the component of the velocity and
                  //therefore also the face for which it is the normal component
                  macExtrapFAB[icomp].copy(advVelFAB[icomp]);
                }
            }
        }

      //copy extrapolated values into vector holders
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          Interval srcInterv(0, 0);
          Interval dstInterv(icomp, icomp);
          m_macScratch1[ilev]->copyTo(srcInterv, *macScratchVec[ilev],  dstInterv);
          int ibox = 0;
          for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              Box box = grow(m_grids[ilev].get(dit()), 1);
              box &= m_domain[ilev];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  (*coveredScratchVecLo[ilev])[dit()][idir]->copy(box, dstInterv, box, *(*m_coveredScratchLo[ilev])[dit()][idir], srcInterv);
                  (*coveredScratchVecHi[ilev])[dit()][idir]->copy(box, dstInterv, box, *(*m_coveredScratchHi[ilev])[dit()][idir], srcInterv);
                }

              ibox++;
            }
        }

      //May need predicted velocity for EBAMRNoSubcycle extensions
      storePredictedVelocity(m_macScratch1,
                             m_coveredScratchLo,
                             m_coveredScratchHi,
                             icomp);

    } //end loop over velocity components  (icomp)

  //average down face centered stuff so that it makes sense at coarse-fine interfaces
  averageDown(macScratchVec);

  //compute the actual advective derivative
  //split this out to make it separately testable
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle:computing advective derivative" << endl;
    }
  computeAdvectiveDerivative(a_uDotDelU, m_advVel, macScratchVec,
                             m_coveredAdvVelLo, m_coveredAdvVelHi,
                             coveredScratchVecLo, coveredScratchVecHi);

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      delete macScratchVec[ilev];
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*coveredScratchVecLo[ilev])[dit()][idir];
              delete (*coveredScratchVecHi[ilev])[dit()][idir];
            }
        }
      delete coveredScratchVecLo[ilev];
      delete coveredScratchVecHi[ilev];
    }
}
/*******************/
/*******************/
void
EBAMRNoSubcycle::
computeAdvectiveDerivative(Vector<LevelData<EBCellFAB>* >                     &  a_uDotDelU,
                           Vector<LevelData<EBFluxFAB>* >                     &  a_macAdvVel,
                           Vector<LevelData<EBFluxFAB>* >                     &  a_macScalar,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelLo,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredAdvVelHi,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarLo,
                           Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredScalarHi,
                           bool                                                  a_nonConsOnly,
                           bool                                                  a_consOnly,
                           Vector<RefCountedPtr<EBLevelAdvect> >              *  a_ebLevAd)
{
  CH_TIME("EBAMRNoSubcycle::computeAdvectiveDerivative2");
  int ncomp = a_uDotDelU[0]->nComp();
  //compute the actual advective derivative
  if (a_ebLevAd == NULL)
  {
    a_ebLevAd = &m_ebLevAd;
  }
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      a_macScalar[ilev]->exchange(Interval(0,ncomp-1));
      a_macAdvVel[ilev]->exchange(Interval(0,0));
      //compute scalar advective derivative
      //then copy into the appropriate
      //component of the input data holder
      RealVect dxLev = m_dx[ilev]*RealVect::Unit;
      ProblemDomain domLev = m_domain[ilev];
      IntVectSet cfivs;
      int jbox = 0;
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          EBPatchAdvect& patcher = (*a_ebLevAd)[ilev]->getPatchAdvect(dit());
          //EBCellFAB& udeluFAB = (*a_uDotDelU[ilev])[dit()];
          //this takes the non-conservative advective derivative = udotdelu
          patcher.advectiveDerivative((*a_uDotDelU[ilev])[dit()],
                                      (*a_macScalar[ilev])[dit()],        //contains extrapolated vel component (facerho)
                                      (*a_macAdvVel[ilev])[dit()],        //facevel
                                      (*a_coveredScalarLo[ilev])[dit()],  //coveredRhoMinu
                                      (*a_coveredScalarHi[ilev])[dit()],  //coveredRhoPlus
                                      (*a_coveredAdvVelLo[ilev])[dit()],  //coveredVelMinu
                                      (*a_coveredAdvVelHi[ilev])[dit()],  //coveredVelPlus
                                      (*m_coveredFaceLitLo[ilev])[dit()], //coveredFaceMinu
                                      (*m_coveredFaceLitHi[ilev])[dit()], //coveredFacePlus
                                      m_grids[ilev].get(dit()) );

          jbox++;
        }

      //these are not grown by one.
      LayoutData<IntVectSet> irregSetsSmall;
      //these are grown by one in the directions != idir
      LayoutData<IntVectSet> irregSetsGrown[SpaceDim];
      LevelData<BaseIFFAB<Real> > fluxInterpolants[SpaceDim];

      irregSetsSmall.define(m_grids[ilev]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = m_ebisl[ilev][dit()];
          if (!ebisBox.isAllCovered())
            {
              const Box&  thisBox = m_grids[ilev].get(dit());
              irregSetsSmall[dit()] = ebisBox.getIrregIVS(thisBox);
            }
        }

      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          EBArith::defineFluxInterpolant(fluxInterpolants[faceDir],
                                         irregSetsGrown  [faceDir],
                                         m_grids[ilev], m_ebisl[ilev],
                                         m_domain[ilev], ncomp, faceDir);
        }


      //set up flux = u*s on irregular sets
      //first put cell-face centered u*s into fluxInterpolant
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          //first put cell-face centered u*s into fluxInterpolant
          const EBISBox& ebisBox =  m_ebisl[ilev][dit()];
          const Box& cellBox = m_grids[ilev].get(dit());
          //put the interpolant = vel*s  into interpolantGrid
          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              IntVectSet ivsIrregGrown = irregSetsGrown[faceDir][dit()];
              ivsIrregGrown &= cellBox;
              FaceStop::WhichFaces stopCrit = FaceStop::SurroundingWithBoundary;

              BaseIFFAB<Real>& interpol = fluxInterpolants[faceDir][dit()];
              interpol.setVal(7.7777e7);
              EBFaceFAB& velDir = (*a_macAdvVel[ilev])[dit()][faceDir];
              EBFaceFAB& scaDir = (*a_macScalar[ilev])[dit()][faceDir];
              for (FaceIterator faceit(ivsIrregGrown, ebisBox.getEBGraph(),
                                      faceDir, stopCrit);
                  faceit.ok(); ++faceit)
                {
                  for (int ivar = 0; ivar < ncomp; ivar++)
                    {
                      Real vel = velDir(faceit(), 0);
                      Real sca = scaDir(faceit(), ivar);
                      interpol(faceit(), ivar) = vel*sca;
                    }
                }
            }
        }

      //exchange ghost cell data for flux interpolant
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          fluxInterpolants[faceDir].exchange(Interval(0, ncomp-1));
        }

      //just doing redistribution over a level here since the EB and CF do not intersect
      BaseIVFactory<Real> ivfact(m_ebisl[ilev], irregSetsSmall);
      LevelData<BaseIVFAB<Real> > massDiffLD(m_grids[ilev], ncomp, 2*IntVect::Unit, ivfact);
      LevelData<BaseIVFAB<Real> > consDivLD(m_grids[ilev], ncomp, IntVect::Unit, ivfact);//used when advectingScalar
      Interval consInterv(0, ncomp-1);
      EBLevelRedist levelRedist(m_grids[ilev], m_ebisl[ilev], m_domain[ilev], ncomp);
      levelRedist.setToZero();

      int ibox = 0;
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          const EBISBox& ebisBox = m_ebisl[ilev][dit()];

          //define centroid flux and do interpolation
          BaseIFFAB<Real> centroidFlux[SpaceDim];
          const BaseIFFAB<Real>* interpolantGrid[SpaceDim];
          const IntVectSet& ivsIrregSmall = irregSetsSmall[dit()];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              const BaseIFFAB<Real>& interpol = fluxInterpolants[idir][dit()];
              interpolantGrid[idir] = &interpol;
              BaseIFFAB<Real>& fluxDir= centroidFlux[idir];
              fluxDir.define(ivsIrregSmall, ebisBox.getEBGraph(), idir, ncomp);
            }

          EBPatchAdvect& patcher  = (*a_ebLevAd)[ilev]->getPatchAdvect(dit());

          EBPatchGodunov& patchGod = patcher;
          patchGod.interpolateFluxToCentroids(centroidFlux,
                                              interpolantGrid,
                                              ivsIrregSmall);

          BaseIVFAB<Real>  ebFlux(ivsIrregSmall, ebisBox.getEBGraph(), ncomp);
          BaseIVFAB<Real> consDiv(ivsIrregSmall, ebisBox.getEBGraph(), ncomp);
          ebFlux.setVal(0.0);
          consDiv.setVal(0.0);

          //this fills  consDiv with kappa*div(US)
          patchGod.consUndividedDivergence(consDiv, centroidFlux, ebFlux, ivsIrregSmall);

          //a_udotdelu holds udels(NC).  make it hold (kappa*div(us) + (1-kappa)udels(NC))
          EBCellFAB&        udelsFAB = (*a_uDotDelU[ilev])[dit()];
          BaseIVFAB<Real>&  massDiff = massDiffLD[dit()];
          massDiff.setVal(0.);
          for (VoFIterator vofit(ivsIrregSmall, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              Real kappa = ebisBox.volFrac(vof);
              for (int ivar = 0; ivar < ncomp; ivar++)
                {
                  Real nonConsDiv  = udelsFAB(vof, ivar);
                  Real kappaConsDiv = consDiv(vof, ivar);
                  const IntVect iv = vof.gridIndex();
                  CH_assert(!a_nonConsOnly || !a_consOnly);
                  if (!a_consOnly && !a_nonConsOnly)
                    {
                      udelsFAB(vof, ivar) = (1.0-kappa)*nonConsDiv + kappaConsDiv;
                      if (ilev == m_finestLevel)
                        {
                          massDiff(vof, ivar) = -(1.0-kappa)*(kappa*nonConsDiv - kappaConsDiv);
                        }
                    }
                  else if (a_nonConsOnly)
                    {
                      udelsFAB(vof, ivar) = nonConsDiv;
                    }
                  else if (a_consOnly)
                    {
                      udelsFAB(vof, ivar) = kappaConsDiv;
                    }
                  else
                    {
                      MayDay::Error("computeAdvectiveDerivative: should not get here");
                    }
                }//ivar
            }//vofit

          if (!a_consOnly && !a_nonConsOnly && ilev == m_finestLevel)
            {
              levelRedist.increment(massDiff, dit(), consInterv);
            }

          ibox++;

        }//dit

      /**/
      //smush the mass difference back in
      if (!a_consOnly && !a_nonConsOnly && ilev== m_finestLevel)
        {
          levelRedist.redistribute(*a_uDotDelU[ilev], consInterv);
        }

    }//ilev
  //for udotdelu, reflux at coarse fine interfaces
  /**/
  if ((ncomp == SpaceDim) && (!a_consOnly) && (!a_nonConsOnly)) //those cases are for convergence tests
    {
      refluxUDotDelU(a_uDotDelU, a_macAdvVel, a_macScalar);
    }
  averageDown(a_uDotDelU);
}
/*******************/
void
EBAMRNoSubcycle::
refluxUDotDelU(Vector<LevelData<EBCellFAB>* >                     &  a_uDotDelU,
               Vector<LevelData<EBFluxFAB>* >                     &  a_macAdvVel,
               Vector<LevelData<EBFluxFAB>* >                     &  a_vecVel)
{
  int velcomp = (*a_vecVel[0]).nComp();
  int uducomp = (*a_uDotDelU[0]).nComp();
  int maccomp = (*a_macAdvVel[0]).nComp();
  if ((velcomp != SpaceDim) || (uducomp != SpaceDim) || (maccomp != 1))
    {
      MayDay::Error("refluxUdotDelu assumptions about components violated");
    }
  pout() << "doing flux matching for u dot del u"  << endl;
  Interval interv(0, SpaceDim-1);
  //flux register for i to i+1 lives at level i
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      if (ilev < m_finestLevel) 
        {
          m_fluxReg[ilev]->setToZero();
        }
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          const Box&     grid    = m_eblg[ilev].getDBL()[dit()];
          const EBISBox& ebisBox = m_eblg[ilev].getEBISL()[dit()];
          //flux is u u  = advection vel at
          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
            {
              EBFaceFAB flux(ebisBox, grid, faceDir, SpaceDim);
              flux.setVal(0.);
              //set flux = vec vel* advection velocity
              flux.plus((*a_vecVel[ilev])[dit()][faceDir], 0, 0, SpaceDim);
              for (int comp = 0; comp < SpaceDim; comp++)
                {
                  flux.mult((*a_macAdvVel[ilev])[dit()][faceDir], 0, comp, 1);
                }
              Real scale = 1; 
              //now increment flux registers with the flux
              if (ilev < m_finestLevel)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      //both here means regular and irregular.   the name is a
                      //relic of a more civilized age
                      m_fluxReg[ilev]->incrementCoarseBoth(flux, scale, dit(), interv, faceDir, sit());
                    }
                }
              if (ilev > 0)
                {
                  for (SideIterator sit; sit.ok(); ++sit)
                    {
                      //both here means regular and irregular.   the name is a
                      //relic of a more civilized age
                      m_fluxReg[ilev-1]->incrementFineBoth(flux, scale, dit(), interv, faceDir, sit());
                    }
                }
            }
        }
    }
  //now change the advective deriviate based upon the flux difference
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      if (ilev < m_finestLevel) 
        {
          Real scale = 1.0/m_dx[ilev];
          m_fluxReg[ilev]->reflux((*a_uDotDelU[ilev]),  interv,  interv, scale);
          m_fluxReg[ilev]->setToZero();
        }
    }
  
}

/*******************/
 void
EBAMRNoSubcycle::
extrapolateScalarCol(Vector<LevelData<EBFluxFAB>* >                     &  a_macScalar,
                     Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                     Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                     const RefCountedPtr<EBPhysIBCFactory>              &  a_advectBC,
                     const Vector<LevelData<EBFluxFAB>* >               &  a_advectiveVel,
                     const Vector<LevelData<EBCellFAB>*>*                  a_sourceTerm,
                     const Vector<LevelData<EBCellFAB>* >               &  a_cellScalar,
                     const Vector<LevelData<EBCellFAB>* >               &  a_cellVelocity,
                     Vector<RefCountedPtr<EBLevelAdvect> >              *  a_ebLevAd)
{
  CH_TIME("EBAMRNoSubcycle::extrapolateScalarCol");

  if (a_ebLevAd == NULL)
  {
    a_ebLevAd = &m_ebLevAd;
  }
  //extrapolate velocity component to faces
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {

      EBPatchGodunov::setCurLevel(ilev);
      LevelData<EBCellFAB>* coarDataOld = NULL;
      LevelData<EBCellFAB>* coarDataNew = NULL;

      LevelData<EBCellFAB>* coarVeloOld = NULL;
      LevelData<EBCellFAB>* coarVeloNew = NULL;

      DisjointBoxLayout coarDBL;
      EBISLayout        coarEBISL;
      int refRat = 2;
      if (ilev > 0)
        {
          coarDBL =  m_grids[ilev-1];
          coarEBISL= m_ebisl[ilev-1];
          coarDataOld = a_cellScalar[ilev-1];
          coarDataNew = a_cellScalar[ilev-1];
          coarVeloOld = a_cellVelocity[ilev-1];
          coarVeloNew = a_cellVelocity[ilev-1];
          refRat = m_params.m_refRatio[ilev-1];
        }

      RefCountedPtr<EBLevelAdvect> ebLevelAdvect = (*a_ebLevAd)[ilev];
      //need to reset boundary conditions because this object was defined with dummy bcs
      //because it is reused over several different variables

      ebLevelAdvect->resetBCs(a_advectBC);

      //this advects all SpaceDim components of velocity to
      //each face.  put into macscratch.
      LevelData<EBCellFAB>* source = NULL;
      LevelData<EBCellFAB>* sourceCoarOld = NULL;
      LevelData<EBCellFAB>* sourceCoarNew = NULL;
      if (a_sourceTerm != NULL)
        {
          source = (*a_sourceTerm)[ilev];
          if (ilev > 0)
            {
              sourceCoarOld = (*a_sourceTerm)[ilev-1];
              sourceCoarNew = (*a_sourceTerm)[ilev-1];
            }
        }
      ebLevelAdvect->advectToFacesCol(*a_macScalar[ilev],     //extrapolated component idir
                                      *a_coveredMacLo[ilev],
                                      *a_coveredMacHi[ilev],
                                      *m_coveredFaceLitLo[ilev],
                                      *m_coveredFaceLitHi[ilev],
                                      *m_coveredSetsLitLo[ilev],
                                      *m_coveredSetsLitHi[ilev],
                                      *a_cellScalar[ilev],    //consstate
                                      *a_cellVelocity[ilev],  //use velo as normal velocity
                                      *a_advectiveVel[ilev],  //contains initial advective velocity
                                      coarDataOld,
                                      coarDataNew,
                                      coarVeloOld,
                                      coarVeloNew,
                                      m_time, m_time, m_time, m_dt,
                                      source, sourceCoarOld, sourceCoarNew);

      a_macScalar[ilev]->exchange(Interval(0,0));
    }

}
/*******************/
/*******************/
void
EBAMRNoSubcycle::
extrapolateToCoveredFaces(Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                          Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                          Vector<LevelData<EBFluxFAB>* >&                       a_macOpen,
                          Vector<LevelData<EBCellFAB>* >&                       a_cellOpen,
                          int                                                   a_idir,
                          Vector<RefCountedPtr<EBLevelAdvect> >              *  a_ebLevAd)
{
  CH_TIME("EBAMRNoSubcycle::extrapolateToCoveredFaces");
  if (a_ebLevAd == NULL)
  {
    a_ebLevAd = &m_ebLevAd;
  }
  // averageDown(a_macOpen);
  //extrapolate the projected velocity to covered faces
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      ProblemDomain curDomain(m_domain[ilev]);

      LevelData<EBFluxFAB>& faceValue = (LevelData<EBFluxFAB>&) *a_macOpen[ilev];
      faceValue.exchange(Interval(0,0));

      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {

          EBPatchAdvect& patcher  = (*a_ebLevAd)[ilev]->getPatchAdvect(dit());
          patcher.extrapToCoveredFaces(*(*a_coveredMacLo[ilev])[dit()][a_idir],
                                       (*a_macOpen[ilev])[dit()][a_idir],
                                       (*a_cellOpen[ilev])[dit()],
                                       (*m_coveredFaceLitLo[ilev])[dit()][a_idir],
                                       a_idir, Side::Lo, m_grids[ilev].get(dit()));

          patcher.extrapToCoveredFaces(*(*a_coveredMacHi[ilev])[dit()][a_idir],
                                       (*a_macOpen[ilev])[dit()][a_idir],
                                       (*a_cellOpen[ilev])[dit()],
                                       (*m_coveredFaceLitHi[ilev])[dit()][a_idir],
                                       a_idir, Side::Hi, m_grids[ilev].get(dit()));
        }
    }
}
/*****************/
/*****************/
void
EBAMRNoSubcycle::
viscousSourceForAdvect(Vector<LevelData<EBCellFAB>* >&       a_source,
                       Vector<LevelData<EBCellFAB>* >&       a_velComp,
                       Vector<LevelData<EBCellFAB>* >&       a_zero,
                       int                                   a_icomp,
                       Real                                  a_time)
{
  CH_TIME("EBAMRNoSubcycle::viscousSourceForAdvect");
  EBAMRDataOps::setToZero(a_source);
  EBAMRDataOps::setToZero(a_zero);
  CH_assert(a_velComp[0]->nComp() == 1);

  applyEBAMROp(a_source,  //returns holding kappa*nu*lapl(vel comp)
               a_velComp, //holds cell centered vel comp
               a_zero,    //holds the zero for the residual calc
               a_icomp,
               a_time);  //velocity component

  //now fill the ghost cells of the laplacian over coarse-fine interfaces with
  //constant extrapolation from neighboring valid values
  //this includes an exchange
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      Interval interv(0, 0);
      EBNormalizeByVolumeFraction normalized_source(m_eblg[ilev]);
      normalized_source(*a_source[ilev],interv);
      if (ilev==0)
        {
          //fill fine-fine ghost cells only
          EBLevelDataOps::exchangeAll(*a_source[ilev]);
        }
      else//exchange happens in here
        {
          IntVect ivGhost = a_source[ilev]->ghostVect();
          EBConstantCFInterp interpolator(m_grids[ilev], m_ebisl[ilev], m_domain[ilev], ivGhost);
          interpolator.interpolate(*a_source[ilev]);
        }
    }
}
/*****************/
/*****************/
void
EBAMRNoSubcycle::
applyEBAMROp(Vector<LevelData<EBCellFAB>* >&       a_lap,
             Vector<LevelData<EBCellFAB>* >&       a_phi,
             Vector<LevelData<EBCellFAB>* >&       a_zero,
             int                                   a_velComp,
             Real                                  a_time)
{
  CH_TIME("EBAMRNoSubcycle::applyEBAMROp");
  Real alpha = 0.0;
  Real beta = 1.0;
  if (m_params.m_orderTimeIntegration == 2)
    {
      m_tgaSolver[a_velComp]->setTime(a_time);
      m_tgaSolver[a_velComp]->resetAlphaAndBeta(alpha, beta);
    }
  else if (m_params.m_orderTimeIntegration == 1)
    {
      m_backwardSolver[a_velComp]->resetAlphaAndBeta(alpha, beta);
    }
  else
    {
      MayDay::Error("EBAMRNoSubcycle::applyEBAMROp -- bad order time integration");
    }

  int coarsestLevel = 0;
  bool homogeneousBC = false;
  //apply the operator (by computing the residual with rhs = 0, and * -1
  m_solver[a_velComp]->computeAMRResidual(a_lap,
                                          a_phi,
                                          a_zero,
                                          m_finestLevel,
                                          coarsestLevel,
                                          homogeneousBC);

  EBAMRDataOps::scale(a_lap,-1.0);
  EBAMRDataOps::setCoveredVal(a_lap,0.0);
  EBAMRDataOps::setCoveredAMRVal(a_lap,m_ebisl,m_params.m_refRatio,0.0);

}
/*****************/
/******************/
void
EBAMRNoSubcycle::
corrector()
{
  correctExtraFieldsPreVelocity();
  if (!m_steadyState)
    {
      correctVelocity();
    }
  correctExtraFieldsPostVelocity();
}
/*****************/
/*****************/
void
EBAMRNoSubcycle::
correctVelocity()
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::correctVelocity" << endl;
    }
  CH_TIMERS("EBAMRNoSubcycle::correctVelocity");
  CH_TIMER("inviscid_advance", t1);
  CH_TIMER("viscous_advance", t2);
  CH_TIMER("cell-centered_projection", t3);

  Interval interv(0, SpaceDim-1);
  Vector<LevelData<EBCellFAB>* > tempLDPtr;
  tempLDPtr.resize(m_finestLevel+1);
  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      tempLDPtr[ilev] = new LevelData<EBCellFAB>();
      tempLDPtr[ilev]->define(m_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact);
      m_velo[ilev]->copyTo(interv, *(tempLDPtr[ilev]), interv);
    }

  Vector<LevelData<EBCellFAB>*> extraSource;
  extraSource.resize(m_finestLevel+1, NULL);
  for (int ilev=0; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      extraSource[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);
    }
  EBAMRDataOps::setToZero(extraSource);

  //computes extra source term vector to be added below
  computeExtraSourceForCorrector(extraSource);

  if (!m_viscousCalc)
    {
      CH_START(t1);
      //for inviscid calc, do not add pressure gradient into
      //the update so we don't have to iterate for the pressure gradient after regrid.
      //If you use the pressure grad in the update and do not iterate after regrid,
      //the solution can ring.
      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              EBCellFAB& newVel = (*m_velo[ilev])[dit()];
              EBCellFAB& uDotDelU = (*m_uDotDelU[ ilev])[dit()];
              EBCellFAB& gradPres = (*m_gphi [ ilev])[dit()];
              EBCellFAB& source = (*extraSource[ ilev])[dit()];
              //make udotdelu = udotdelu + gradp
              uDotDelU += gradPres;

              //make udotdelu = -udotdelu - gradp
              uDotDelU *= -1.0;

              //adding extra source
              uDotDelU += source;

              //now udotdelu = dt*(-udotdelu - gradp);
              uDotDelU *=  m_dt;

              //put change of velocity into newvel
              newVel += uDotDelU;
            }
        }
      CH_STOP(t1);
    }
  else
    {
      CH_START(t2);
      //add gradient of pressure into udotdelu
      EBAMRDataOps::incr(m_uDotDelU, m_gphi, 1.0);
      //make udelu = -udelu-gradp == the source term of heat eqn
      EBAMRDataOps::scale(m_uDotDelU, -1.0);

      EBAMRDataOps::incr(m_uDotDelU, extraSource, 1.0);

      if (m_params.m_verbosity >= 2)
        {
          pout() << "EBAMRNoSubcycle::solving implicitly for viscous and any extra source terms" << endl;
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //make cellscratch = udotdelu
          EBAMRDataOps::setToZero(m_cellScratc2);
          EBAMRDataOps::setToZero(m_cellScratc1);
          EBAMRDataOps::setToZero(m_cellScratch);
          //put component of rhs into cellscratc2
          //put component of velo into cellscratch
          //new velocity comes out in cellscratc1
          Interval vecInterv(idir, idir);
          Interval scaInterv(0, 0);
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              m_uDotDelU[ilev]->copyTo(vecInterv,  *m_cellScratc2[ilev], scaInterv);
              m_velo[ ilev]->copyTo(vecInterv,  *m_cellScratch[ilev], scaInterv);
              m_velo[ ilev]->copyTo(vecInterv,  *m_cellScratc1[ilev], scaInterv);
            }

          //tell EBBC which velocity component we are solving for
          DirichletPoissonEBBC::s_velComp = idir;

          //solve the stinking equation
          int lbase = 0;
          int lmax = m_finestLevel;
          if (m_params.m_orderTimeIntegration == 2)
            {
              m_tgaSolver[idir]->oneStep(m_cellScratc1, //vel new
                                         m_cellScratch, //vel old
                                         m_cellScratc2, //source
                                         m_dt,
                                         lbase,
                                         lmax,
                                         m_time);
            }
          else if (m_params.m_orderTimeIntegration == 1)
            {
              m_backwardSolver[idir]->oneStep(m_cellScratc1, //vel new
                                              m_cellScratch, //vel old
                                              m_cellScratc2, //source
                                              m_dt,
                                              lbase,
                                              lmax,
                                              false);//do not zero phi
            }
          else
            {
              MayDay::Error("EBAMRNoSubcycle::correctVelocity -- bad order time integration");
            }

          //now copy the answer back from scratch into velo
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              m_cellScratc1[ilev]->copyTo(scaInterv, *m_velo[ilev],  vecInterv);
            }
        }
      CH_STOP(t2);
    }

  if (m_params.m_verbosity >= 2)
    {
      pout() << "EBAMRNoSubcycle::cc projecting velocity" << endl;
    }

  EBAMRDataOps::setCoveredVal(m_velo,0.0);
  EBAMRDataOps::setCoveredAMRVal(m_velo,m_ebisl,m_params.m_refRatio,0.0);
  EBAMRDataOps::setCoveredVal(m_gphi,0.0);
  EBAMRDataOps::setCoveredAMRVal(m_gphi,m_ebisl,m_params.m_refRatio,0.0);

  //make u* := u* + dt*gphi
  //this makes the output pressure (I-P)u*= gphi*dt
  EBAMRDataOps::incr(m_velo, m_gphi, m_dt);
  //this puts  into gphi (new pressure gradient)*dt
  //since we have added only a pure gradient, P(u*) = u^n+1 still
  CH_START(t3);
  Real time = m_time;
  if (!m_advanceGphiOnly)//rob thinks we want to use the new time for BCs only if not iterating grad(phi)
    {
      time += m_dt;
    }
  m_ccProjector->setTime(time);
  //m_pres already scaled with dt/2 above in MAC projection
  EBAMRDataOps::scale(m_pres, 2.);
  m_ccProjector->setInitialPhi(m_pres);
  m_ccProjector->project(m_velo, m_gphi);

  const Vector<LevelData<EBCellFAB>* >& projPhi = m_ccProjector->getPhi();
  EBAMRDataOps::assign(m_pres, projPhi);
  CH_STOP(t3);

  filter(m_velo);

  //pres now holds (new pressure)*dt
  //gphi now holds (new pressure gradient)*dt
  //so divide out the dt
  EBAMRDataOps::scale(m_pres, 1.0/m_dt);
  EBAMRDataOps::scale(m_gphi, 1.0/m_dt);

  //post-processing
  if (m_params.m_verbosity > 3)
    {
      m_ccProjector->kappaDivergence(m_cellScratch, m_velo);
      Real norm[3];
      for (int inorm = 0; inorm < 3; inorm++)
        {
          norm[inorm] = EBArith::norm(m_cellScratch, m_grids, m_ebisl,
                                      m_params.m_refRatio, 0, inorm, EBNormType::OverBoth);
        }
      pout() << "div(vel) after cell project and filter: " <<
        "L_inf = " << norm[0]  <<
        ", L_1 = " << norm[1] <<
        ", L_2 = " << norm[2] << endl;
    }

  //check steady-state
  //don't change tempLDPtr if using it to cache velocity during priming
  if (!m_advanceGphiOnly && m_params.m_doSteadyState)
    {
      EBAMRDataOps::incr(tempLDPtr, m_velo, -1.0);
      EBAMRDataOps::scale(tempLDPtr, 1./m_dt);
      Interval scaInterv(0, 0);
      m_steadyState = true;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Interval vecInterv(idir, idir);
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              tempLDPtr[ ilev]->copyTo(vecInterv,  *m_cellScratch[ilev], scaInterv);
            }
          Real norm[3];
          for (int inorm = 0; inorm < 3; inorm++)
            {
              norm[inorm] = EBArith::norm(m_cellScratch, m_grids, m_ebisl,
                                          m_params.m_refRatio, 0, inorm, EBNormType::OverBoth);
            }

          if (norm[0] < 1.e-8)
            {
              pout() << setprecision(8)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific);
              pout() << "Steady-state velocity" << idir << ": " <<
                "L_inf = " << norm[0] <<
                ", L_1 = " << norm[1] <<
                ", L_2 = " << norm[2] << endl;
            }
          else
            {
              m_steadyState = false;
            }
        }
      if (m_steadyState)
        {
          m_stopAdvance = true;
        }
    }

  //overwrite velocity with previous during initialization or after regridding
  if (m_advanceGphiOnly)
    {
      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          tempLDPtr[ilev]->copyTo(interv, *m_velo[ilev], interv);
        }
    }

  for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
      delete tempLDPtr[ilev];
      delete extraSource[ilev];
    }

  EBAMRDataOps::setCoveredVal(m_velo, 0.0);
  EBAMRDataOps::setCoveredVal(m_gphi, 0.0);
  EBAMRDataOps::setCoveredVal(m_pres, 0.0);
}
/*****************/
/*****************/
void
EBAMRNoSubcycle::
pointsUpdated()
{
  CH_TIME("EBAMRNoSubcycle::pointsUpdated");
  long long totalPoints = 0;
  long long totalBoxes  = 0;
  int numLevels = m_finestLevel + 1;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      long long pointsThisLevel = 0;
      for (LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
        {
          pointsThisLevel += m_grids[ilev][lit()].numPts();
        }
      totalPoints += pointsThisLevel;
      totalBoxes += m_grids[ilev].size();
    }

  m_pointsUpdated += totalPoints;
}
/*****************/
#ifdef CH_USE_HDF5
/*****************/
void
EBAMRNoSubcycle::writePlotFile()
{
  CH_TIME("EBAMRNoSubcycle::writePlotFile");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::writePlotFile" << endl;
    }
  int curNumLevels = m_finestLevel + 1;

  Vector<string> presnames(1, string("pressure"));
#if CH_SPACEDIM == 2
  Vector<string> vortnames(1, string("vorticity"));
#else
  Vector<string> vortnames(SpaceDim);
#endif

  Vector<string> velonames(SpaceDim);
  Vector<string> gphinames(SpaceDim);
  Vector<string> names;
  for (int idir = 0; idir < SpaceDim; idir++)
    {

      char velochar[100];
      char gphichar[100];
      sprintf(velochar, "velocity%d", idir);
      sprintf(gphichar, "gradPres%d", idir);
      velonames[idir] = string(velochar);
      gphinames[idir] = string(gphichar);

#if CH_SPACEDIM==3
      char vortchar[100];
      sprintf(vortchar, "vorticity%d", idir);
      vortnames[idir] = string(vortchar);
#endif

    }

  names = velonames;
  names.append(gphinames);
  names.append(presnames);
  names.append(vortnames);
  //For adding extra data
  names.append(extraNames());

  char fileChar[1000];
  int ncells = m_domain[0].size(0);
  sprintf(fileChar, "plot.nx%d.step.%07d.%dd.hdf5", ncells, m_curStep, SpaceDim);

  bool replaceCovered = false;
  Vector<Real> coveredValues;

  int nlev = m_finestLevel + 1;
  Vector<LevelData<EBCellFAB>* > vorticity(nlev, NULL);
  Vector<LevelData<EBCellFAB>* > outputData(nlev, NULL);

  for (int ilev = 0; ilev < nlev; ilev++)
    {
      vorticity[ilev] = new LevelData<EBCellFAB>();
      computeVorticity(*vorticity[ilev], ilev);

#if CH_SPACEDIM==2
      int nvar = 2*SpaceDim + 2;
#else
      int nvar = 3*SpaceDim + 1;
#endif
      //For adding extra data
      int numExtraVars = getNumExtraVars();
      nvar += numExtraVars;
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      outputData[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], nvar, IntVect::Zero, ebcellfact);

      Interval srcInterv, dstInterv;
      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(0, SpaceDim-1);
      m_velo[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(SpaceDim, 2*SpaceDim-1);
      m_gphi[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

      srcInterv = Interval(0, 0);
      dstInterv = Interval(2*SpaceDim, 2*SpaceDim);
      m_pres[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

#if CH_SPACEDIM==2
      srcInterv = Interval(0, 0);
      dstInterv = Interval(2*SpaceDim+1, 2*SpaceDim+1);
      vorticity[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

      //Add extra data
      int startDstInterv = 2*SpaceDim+2;
      addExtraOutputData(*outputData[ilev], startDstInterv, ilev);
#else
      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(2*SpaceDim+1, 3*SpaceDim);
      vorticity[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

      //Add extra data
      int startDstInterv = 3*SpaceDim+1;
      addExtraOutputData(*outputData[ilev],  startDstInterv, ilev);
#endif
      setCoveredStuffToZero(*outputData[ilev]);
    }

  string filename(fileChar);
  writeEBHDF5(filename,
              m_grids,
              outputData,
              names,
              m_domain[0].domainBox(),
              m_dx[0],
              m_dt,
              m_time,
              m_params.m_refRatio,
              curNumLevels,
              replaceCovered,
              coveredValues);

  addExtraDatatoPlotFile(filename);

  for (int ilev = 0; ilev <nlev; ilev++)
    {
      delete vorticity[ilev];
      delete outputData[ilev];
    }
}
/*****************/
/*****************/
void
EBAMRNoSubcycle::writeCheckpointFile()
{
  CH_TIME("EBAMRNoSubcycle::writeCheckpointFile");
  CH_assert(m_isSetup);
  //Setup the level header information
  HDF5HeaderData header;

  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::writeCheckpointFile" << endl;
    }

  //all the stuff in m_params had to come in at define time
  //this is also true of viscosity, the domains, and the boundary conditions
  //viscous calc gets set there too

  //bool conversion to int
  int iuseFixedDt       = 0;
  if (m_useFixedDt) iuseFixedDt = 1;

  header.m_real["time"]                   = m_time;
  header.m_real["dt"]                     = m_dt;
  header.m_real["prescribed_dt"]          = m_prescribedDt;
  header.m_int ["cur_step"]               = m_curStep;
  header.m_int ["finest_level"]           = m_finestLevel;
  header.m_int ["use_fixed_dt"]           = iuseFixedDt;
  header.m_int ["steady_state"]           = m_steadyState;

  char iter_str[100];

  int ncells = m_domain[0].size(0);
  sprintf(iter_str, "check%d.nx%d.%dd.hdf5", m_curStep, ncells, SpaceDim);

  HDF5Handle handleOut(iter_str, HDF5Handle::CREATE);
  //Write the header for this level
  header.writeToFile(handleOut);
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      handleOut.setGroupToLevel(ilev);
      write(handleOut,m_grids[ilev]);
      write(handleOut,*m_velo[ilev],"velo");
      write(handleOut,*m_gphi[ilev],"gphi");
      write(handleOut,*m_pres[ilev],"pres");
      write(handleOut,*m_advVel[ilev],"advVel");
      writeExtraDataToCheckpoint(handleOut, ilev);
    }
  handleOut.close();
}
/*****************/
/*****************/
void
EBAMRNoSubcycle::readCheckpointFile(const string& a_restartFile)
{
  CH_TIME("EBAMRNoSubcycle::readCheckpointFile");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::readCheckpointFile" << endl;
    }

  HDF5Handle handleIn(a_restartFile, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handleIn);

  //all the stuff in m_params had to come in at define time
  //this is also true of viscosity, the domains, and the boundary conditions
  m_time          =   header.m_real["time"]                   ;
  m_dt            =   header.m_real["dt"]                     ;
  m_prescribedDt  =   header.m_real["prescribed_dt"]          ;
  m_curStep       =   header.m_int ["cur_step"]               ;
  m_finestLevel   =   header.m_int ["finest_level"]           ;
  int iuseFixedDt =   header.m_int["use_fixed_dt"];
  m_steadyState   =   header.m_int["steady_state"];
  m_useFixedDt =  (iuseFixedDt == 1);

  int finestLevelFromParmParse   =   m_params.m_maxLevel;
  if (m_finestLevel > finestLevelFromParmParse)
    {
      m_finestLevel = finestLevelFromParmParse;
    }

  //get all the grids
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      handleIn.setGroupToLevel(ilev);
      //Get the grids
      Vector<Box> vboxGrids;
      const int gridStatus = read(handleIn, vboxGrids);
      if (gridStatus != 0)
        {
          MayDay::Error("readCheckpointLevel: file has no grids");
        }

      Vector<int> proc_map;
      EBEllipticLoadBalance(proc_map,vboxGrids, m_domain[ilev], false, m_ebisPtr );

      m_grids[ilev]= DisjointBoxLayout(vboxGrids,proc_map);
    }

  //define stuff using grids
  defineEBISLs();
  defineExtraEBISLs();
  defineNewVel();
  definePressure();
  defineExtraTerms();
  //now input the actual data
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      handleIn.setGroupToLevel(ilev);
      read<EBCellFAB>(handleIn, *m_velo[ilev], "velo", m_grids[ilev], Interval(), false);
      read<EBCellFAB>(handleIn, *m_gphi[ilev], "gphi", m_grids[ilev], Interval(), false);
      read<EBCellFAB>(handleIn, *m_pres[ilev], "pres", m_grids[ilev], Interval(), false);
      read<EBFluxFAB>(handleIn, *m_advVel[ilev], "advVel", m_grids[ilev], Interval(), false);
      readExtraDataFromCheckpoint(handleIn, ilev);
    }
  handleIn.close();

}
/*****************/
#endif //CH_USE_HDF5
/*****************/
void
EBAMRNoSubcycle::
setupForRestart(const string& a_restartFile)
{
  CH_TIME("EBAMRNoSubcycle::setupForRestart");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::setupForRestart" << endl;
    }
  m_isSetup = true;
  m_doRestart  = true;
#ifdef CH_USE_HDF5
  readCheckpointFile(a_restartFile);
#else
  MayDay::Error("cannot restart from checkpoint without hdf5");
#endif
  int finestLevelFromParmParse = m_params.m_maxLevel;
  int oldFinestLevel = m_finestLevel;
  if (oldFinestLevel < finestLevelFromParmParse)
    {
      //cache data on old levels
      Interval interv(0, SpaceDim-1);
      Vector<LevelData<EBCellFAB>* > tempLDPtr, tempLDPtr2;
      tempLDPtr.resize(oldFinestLevel+1);
      tempLDPtr2.resize(oldFinestLevel+1);
      for (int ilev=0; ilev<= oldFinestLevel; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          tempLDPtr[ilev] = new LevelData<EBCellFAB>();
          tempLDPtr[ilev]->define(m_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact);
          m_velo[ilev]->copyTo(interv, *(tempLDPtr[ilev]), interv);
          tempLDPtr2[ilev] = new LevelData<EBCellFAB>();
          tempLDPtr2[ilev]->define(m_grids[ilev], SpaceDim,  IntVect::Zero, ebcellfact);
          m_gphi[ilev]->copyTo(interv, *(tempLDPtr2[ilev]), interv);
        }

      cacheExtraTerms(oldFinestLevel);

      m_finestLevel = finestLevelFromParmParse;

      bool moreLevels = (m_finestLevel > oldFinestLevel);

      BRMeshRefine meshrefine;
      if (moreLevels)
        {
          meshrefine.define(m_domain[0],              m_params.m_refRatio,
                            m_params.m_fillRatio,     m_params.m_blockFactor,
                            m_params.m_nestingRadius, m_params.m_maxBoxSize);
        }

      Vector<Vector<Box> > new_boxes;
      Vector<Vector<Box> > old_boxes(m_finestLevel+1);

      for (int ilev=0; ilev<m_finestLevel+1; ilev++)
        {
          if (ilev>oldFinestLevel)
            {
              old_boxes[ilev] = Vector<Box>(1, m_domain[ilev].domainBox());
            }
          else
            {
              old_boxes[ilev] = m_grids[ilev].boxArray();
            }
        }

      int base_level = oldFinestLevel;
      int top_level = m_finestLevel-1;

      IntVectSet tagsVect;
      tagCellsLevel(tagsVect, 0);

      int new_finest = meshrefine.regrid(new_boxes, tagsVect,
                                         base_level, top_level,
                                         old_boxes);
      if (new_finest != m_finestLevel)
        {
          MayDay::Error("EBAMRNoSubcycle::setupForRestart -- new_finest not equal m_finestLevel");
        }

      //from defineGrids
      for (int ilev=oldFinestLevel+1; ilev<=m_finestLevel; ilev++)
        {
          Vector<int> procAssign;
          mortonOrdering((Vector<Box>&)(new_boxes[ilev]));
          EBEllipticLoadBalance(procAssign,  new_boxes[ilev], m_domain[ilev], false, m_ebisPtr );
          m_grids[ilev] = DisjointBoxLayout();
          m_grids[ilev].define(new_boxes[ilev], procAssign);
        }

      int startLevel = oldFinestLevel+1;
      int endLevel = oldFinestLevel;
      defineEBISLs();
      defineExtraEBISLs();
      defineNewVel(startLevel);
      definePressure(startLevel);
      defineExtraTerms(startLevel);
      defineProjections();

      for (int ilev=0; ilev<=oldFinestLevel; ilev++)
        {
          tempLDPtr[ilev]->copyTo(interv, *m_velo[ilev], interv);
          tempLDPtr2[ilev]->copyTo(interv, *m_gphi[ilev], interv);
        }

      for (int ilev=oldFinestLevel+1; ilev<=m_finestLevel; ilev++)
        {
          //interpolate everywhere
          EBPWLFineInterp ebInterpVec(m_grids[ ilev  ],
                                      m_grids[ ilev-1],
                                      m_ebisl[ ilev  ],
                                      m_ebisl[ ilev-1],
                                      m_domain[ilev-1],
                                      m_params.m_refRatio[ilev-1],
                                      SpaceDim,
                                      m_ebisPtr);
          ebInterpVec.interpolate(*m_velo[ilev  ],
                                  *m_velo[ilev-1],
                                  interv);
          ebInterpVec.interpolate(*m_gphi[ilev  ],
                                  *m_gphi[ilev-1],
                                  interv);
        }

      interpolateExtraTermsForRestart(startLevel, endLevel);

      m_ccProjector->project(m_velo, m_gphi);
      filter(m_velo);
      defineIrregularData();
      postInitialize();
    }
  else  if (oldFinestLevel > finestLevelFromParmParse)
    {
      m_finestLevel = finestLevelFromParmParse;
      defineProjections();
      defineIrregularData();
      postInitialize();
    }
  else//no change in finest level from inputs
    {
      defineProjections();
      defineIrregularData();
      postInitialize();
    }
}
/*****************/
/*****************/
void
EBAMRNoSubcycle::
setupForFixedHierarchyRun(const Vector<Vector<Box> >& a_grids)
{
  CH_TIME("EBAMRNoSubcycle::setupForFixedHierarchyRun");
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRNoSubcycle::setupForFixedHierarchyRun" << endl;
    }
  //turn off regridding
  m_params.m_regridInterval = -1;
  m_isSetup = true;

  m_finestLevel = a_grids.size() - 1;
  for (int ilev = 0; ilev < m_finestLevel; ilev++)
    {
      CH_assert(a_grids[ilev].size() > 0);
    }

  initialGrid(a_grids);
  initialData();
  defineIrregularData();

  //finally, call post-initialization
  postInitialize();
  m_doRestart  = false;
}
/*****************/
void
EBAMRNoSubcycle::
setCoveredStuffToZero(LevelData<EBCellFAB>& a_vort)
{
  CH_TIME("setCoveredStuffToZero");
  for (DataIterator dit = a_vort.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&  vortFAB =a_vort[dit()];
      Real covVal = 0.0;
      for (int icomp = 0; icomp < vortFAB.nComp(); icomp++)
        {
          vortFAB.setCoveredCellVal(covVal, icomp);
        }
    }
}

#include "NamespaceFooter.H"
