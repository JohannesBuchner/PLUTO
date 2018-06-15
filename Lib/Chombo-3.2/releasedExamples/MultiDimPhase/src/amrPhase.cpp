#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include "CH_Timer.H"
#include "FABView.H"
#include "GradientF_F.H"
#include "ReductionOps.H"
#include "1DPotentialSolve.H"

using std::ifstream;
using std::ios;

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::min;

#ifndef CH_MULTIDIM
#include "AMRSolver.H"
#include "PoissonBC.H"
#include "GhostBC.H"
#endif


#include "Box.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "ParmParse.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "parstream.H"
//#include "CoarseAverage.H"
//#include "CoarseAverageEdge.H"
//#include "FineInterp.H"
#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "MayDay.H"
#include "amrPhase.H"
#include "computeNorm.H"
#include "Divergence.H"
#include "EdgeToCell.H"
//#include "PiecewiseLinearFillPatch.H"
#include "AdvectIBC.H"
#include "CellToEdge.H"

#include "CH_HDF5.H"

#include "phaseF_F.H"
#include "PROB_F_F.H"

#include "UsingNamespace.H"


///
enum BCType
{
  ///
  bogusBC = -1,
  ///
  DIRICHLETBC,
  ///
  NEUMANNBC,
  ///
  SOLIDWALLBC,
  ///
  NUM_BC_TYPES
};

bool
amrPHASE::isDefined() const
{
  return m_is_defined;
}


amrPHASE::amrPHASE()
{
  setDefaults();
}

void
amrPHASE::setDefaults()
{
  // set some bogus values as defaults
  m_is_defined = false;
  m_max_level = -1;
  m_finest_level = -1;
  m_block_factor = -1;
  m_fill_ratio = -1;
  m_num_comps = -1;
  m_do_restart = false;
  m_restart_step = -1;
  m_dt = -10000000.0;
  m_scalIBC = NULL;


  // set the rest of these to reasonable defaults
  m_nesting_radius = 1;
  m_tagging_val = 0.1;
  m_tags_grow = 1;
  m_cfl = 0.25;
  m_max_dt_grow = 1.5;
  m_dt = 1.0e20;
  m_max_box_size = 64;
  m_refine_initial_domain = false;
  // domain size default is unit square/cube
  m_domainSize = RealVect::Unit;
  m_G = 1.0;

  // origin default is zero vector
  m_origin = RealVect::Zero;

  m_plot_prefix = "plot";
  m_plot_interval = 10000000;
  m_check_prefix = "chk";
  m_check_interval = -1;

}

amrPHASE::~amrPHASE()
{
  // clean up memory
  for (int lev = 0; lev < m_old_phi.size(); lev++)
    {
      if (m_old_phi[lev] != NULL)
        {
          delete m_old_phi[lev];
          m_old_phi[lev] = NULL;
        }
    }

  for (int lev=0; lev< m_new_phi.size(); lev++)
    {
      if (m_new_phi[lev] != NULL)
        {
          delete m_new_phi[lev];
          m_new_phi[lev] = NULL;
        }
    }

  for (int lev=0; lev< m_CCvelocity.size(); lev++)
    {
      if (m_CCvelocity[lev] != NULL)
        {
          delete m_CCvelocity[lev];
          m_CCvelocity[lev] = NULL;
        }
    }

  for (int lev=0; lev< m_potential.size(); lev++)
    {
      if (m_potential[lev] != NULL)
        {
          delete m_potential[lev];
          m_potential[lev] = NULL;
        }
    }

  for (int lev=0; lev<m_reduceCopiers.size(); lev++)
    {
      if (m_reduceCopiers[lev] != NULL)
        {
          delete m_reduceCopiers[lev];
          m_reduceCopiers[lev] = NULL;
        }

      if (m_spreadCopiers[lev] != NULL)
        {
          delete m_spreadCopiers[lev];
          m_spreadCopiers[lev] = NULL;
        }
    }

  for (int lev=0; lev<m_patchGod.size(); lev++)
    {
      if (m_patchGod[lev] != NULL)
        {
          delete m_patchGod[lev];
          m_patchGod[lev] = NULL;
        }
    }

  if (m_scalIBC != NULL)
    {
      delete m_scalIBC;
      m_scalIBC = NULL;
    }
  // that should be it!

}

void
amrPHASE::initialize()
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::initialize" << endl;
    }

  // first, read in info from parmParse file
  ParmParse ppPhase("phase");
  Vector<int> ancells(SpaceDim);
  Vector<Real> domsize(SpaceDim);
  Vector<Real> origin(SpaceDim);

  // assumption is that domains are not periodic
  bool is_periodic[SpaceDim];
  for (int dir=0; dir<SpaceDim; dir++)
    is_periodic[dir] = false;
  Vector<int> is_periodic_int(SpaceDim, 0);

  int tempBool = m_refine_initial_domain;
  ppPhase.query("refineInitialDomain", tempBool);
  m_refine_initial_domain = (tempBool == 1);

  ppPhase.get("maxLevel", m_max_level);

  ppPhase.get("num_comps", m_num_comps);

  ppPhase.getarr("num_cells", ancells, 0, SpaceDim);

  if (ppPhase.contains("domain_size"))
    {
      ppPhase.getarr("domain_size", domsize, 0, SpaceDim);
      m_domainSize = RealVect(D_DECL(domsize[0], domsize[1], domsize[2]));
    }


  if (ppPhase.contains("origin"))
    {
      ppPhase.getarr("origin", origin, 0, SpaceDim);
      m_origin = RealVect(D_DECL(origin[0], origin[1], origin[2]));
    }


  ppPhase.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      is_periodic[dir] = (is_periodic_int[dir] == 1);
    }

  ppPhase.query("cfl", m_cfl);

  m_initial_cfl = m_cfl;
  ppPhase.query("initial_cfl", m_initial_cfl);

  ppPhase.query("max_dt_grow_factor", m_max_dt_grow);

  ppPhase.query("plot_interval", m_plot_interval);

  ppPhase.query("plot_prefix", m_plot_prefix);

  ppPhase.query("check_interval", m_check_interval);

  ppPhase.query("check_prefix", m_check_prefix);

  if (m_max_level > 0)
    {
      ppPhase.getarr("ref_ratio", m_refinement_ratios, 0, m_max_level);
    }
  else
    {
      m_refinement_ratios.resize(1);
      m_refinement_ratios[0] = -1;
    }

  ppPhase.query("verbosity", m_verbosity);

  ppPhase.get("regrid_interval", m_regrid_interval);

  ppPhase.get("blockFactor", m_block_factor);

  ppPhase.get("fill_ratio", m_fill_ratio);

  ppPhase.query("nestingRadius", m_nesting_radius);

  ppPhase.query("tagging_val", m_tagging_val);

  ppPhase.query("tags_grow", m_tags_grow);

  ppPhase.query("max_box_size", m_max_box_size);

  // get gravitational constant, if present
  ppPhase.query("gravitational_constant", m_G);

  // get boundary conditions for potential solve
  Vector<int> loBCvect(SpaceDim), hiBCvect(SpaceDim);
  ppPhase.getarr("lo_bc", loBCvect, 0, SpaceDim);
  ppPhase.getarr("hi_bc", hiBCvect, 0, SpaceDim);

  Real inflowVel = 1.0;
  ppPhase.query("inflowVel", inflowVel);

  // this is code designed to emulate the lower-dimension potential
  // solve -- unused in the true multidim case
#ifndef CH_MULTIDIM
  int bcInt;
  SideIterator sit;
  // note that the only BC's which are relevant for the potential
  // solve are those in the x-direction.  The y-direction ones are
  // hardwired to Neumann BCs
  int dir = 0;
  for (sit.begin(); sit.ok(); ++sit)
    {
      if (sit() == Side::Lo)
        {
          bcInt = loBCvect[dir];
        }
      else
        {
          bcInt = hiBCvect[dir];
        }

      if (bcInt == DIRICHLETBC)
        {
          DirichletBC thisBC(dir,sit());
          m_potentialBC.setBoxGhostBC(thisBC);
        }
      else if (bcInt == NEUMANNBC)
        {
          NeumannBC thisBC(dir,sit());
          m_potentialBC.setBoxGhostBC(thisBC);
        }
      else if (bcInt == SOLIDWALLBC)
        {
          // homogeneous Neumann
          NeumannBC thisBC(dir,sit());
          m_potentialBC.setBoxGhostBC(thisBC);
        }
      else
        {
          // invalid bc type
          MayDay::Error("Invalid BC Type!");
        }
    } // end loop over sides


  int reducedSpacedim = SpaceDim -1;
  for (dir=1; dir<reducedSpacedim; dir++)
    {
      for (sit.begin(); sit.ok(); ++sit)
        {
          NeumannBC thisBC(dir,sit());
          //D1::NeumannBC thisBC(dir,sit());
          m_potentialBC.setBoxGhostBC(thisBC);
        } // end loop over sides
    } // end loop over directions

  // now add bc's to LevelOps
  m_poissOp.setDomainGhostBC(m_potentialBC);
#endif // end if not true multidim

  // now set up problem domains
  {
    IntVect loVect = IntVect::Zero;
    IntVect hiVect(D_DECL(ancells[0]-1, ancells[1]-1, ancells[2]-1));
    ProblemDomain baseDomain(loVect, hiVect);
    IntVect reducedHiVect(hiVect);
    for (int dir=1; dir<SpaceDim; dir++)
      {
        reducedHiVect[dir] = 0;
      }
    // this is the nx1 domain used for the 1-D solves
    ProblemDomain reducedBaseDomain(loVect, reducedHiVect);

    // now set periodicity
    for (int dir=0; dir<SpaceDim; dir++)
      {
        baseDomain.setPeriodic(dir, is_periodic[dir]);
        reducedBaseDomain.setPeriodic(dir, is_periodic[dir]);
      }

    // now set up vector of domains
    m_amrDomains.resize(m_max_level+1);
    m_amrReducedDomains.resize(m_max_level+1);
    m_amrDx.resize(m_max_level+1);

    m_amrDomains[0] = baseDomain;
    m_amrReducedDomains[0] = reducedBaseDomain;
    m_amrDx[0] = m_domainSize[0]/baseDomain.domainBox().size(0);

    for (int lev=1; lev<= m_max_level; lev++)
      {
        m_amrDomains[lev] = refine(m_amrDomains[lev-1],
                                   m_refinement_ratios[lev-1]);
        m_amrReducedDomains[lev] = refine(m_amrReducedDomains[lev-1],
                                          m_refinement_ratios[lev-1]);
        m_amrDx[lev] = m_amrDx[lev-1]/m_refinement_ratios[lev-1];
      }

  } // leaving problem domain setup scope

  // check to see if we're using predefined grids
  bool usePredefinedGrids = false;
  std::string gridFile;
  if (ppPhase.contains("gridsFile"))
    {
      usePredefinedGrids = true;
      ppPhase.get("gridsFile",gridFile);
    }



  // check to see if we're restarting from a checkpoint file
  if (!ppPhase.contains("restart_file"))
    {
      // if we're not restarting

      // now set up data holders
      m_old_phi.resize(m_max_level+1, NULL);
      m_new_phi.resize(m_max_level+1, NULL);
      m_CCvelocity.resize(m_max_level+1, NULL);
      m_potential.resize(m_max_level+1, NULL);

      int finest_level = -1;
      if (usePredefinedGrids)
        {
          setupFixedGrids(gridFile);
        }
      else
        {
          // now create initial grids
          initGrids(finest_level);
        }

      // set time to be 0
      m_time = 0.0;
      m_cur_step = 0;

      // that should be it
    }
  else
    {
      // we're restarting from a checkpoint file
      string restart_file;
      ppPhase.get("restart_file", restart_file);
      m_do_restart = true;

      this->restart(restart_file);
    }

  // set up counter of number of cells
  m_num_cells.resize(m_max_level+1, 0);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LayoutIterator lit = levelGrids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box& thisBox = levelGrids.get(lit());
          m_num_cells[lev] += thisBox.numPts();
        }
    }


  // set up advection stuff
  Vector<Real> loPhiBCVal(SpaceDim, 0.0);
  Vector<Real> hiPhiBCVal(SpaceDim,0.0);
  ppPhase.queryarr("loPhiBCVal", loPhiBCVal, 0, SpaceDim);
  ppPhase.queryarr("hiPhiBCVal", hiPhiBCVal, 0, SpaceDim);

  AdvectIBC* phiBCptr = new AdvectIBC;

  for (int dir=0; dir<SpaceDim; dir++)
    {
      phiBCptr->setBoundaryValue(loPhiBCVal[dir], dir, Side::Lo);
      phiBCptr->setBoundaryValue(hiPhiBCVal[dir], dir, Side::Hi);
    }

  if (m_scalIBC != NULL) delete m_scalIBC;
  m_scalIBC = phiBCptr;

  m_patchGod.resize(m_max_level+1, NULL);
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      m_patchGod[lev] = new PatchAdvection;

      m_patchGod[lev]->setPhysIBC(m_scalIBC);

      m_patchGod[lev]->define(m_amrDomains[lev], m_amrDx[lev]);
      // don't use artificial viscosity
      bool useArtVisc = false;
      Real artVisc = -1.0;
      m_patchGod[lev]->setArtificialViscosity(useArtVisc,
                                              artVisc);

      // slope computation parameters -- may want to make these
      // user-controlled
      bool fourthOrderSlopes = true;
      //bool fourthOrderSlopes = false;
      bool useFlattening = false;
      bool limitSlopes = false;
      m_patchGod[lev]->setSlopeParameters(fourthOrderSlopes,
                                          useFlattening,
                                          limitSlopes);


    }

}


void
amrPHASE::run(Real a_max_time, int a_max_step)
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::run -- max_time= " << a_max_time
             << ", max_step = " << a_max_step << endl;
    }


  //Real dt;
  // only call computeInitialDt if we're not doing restart
  if (!m_do_restart)
    {
      m_dt = computeInitialDt();
    }
  else
    {
      m_dt = computeDt();
    }


  // advance solution until done
  while ( (a_max_time > m_time) && (m_cur_step < a_max_step)
          && (m_dt > TIME_EPS))
    {

      // dump plotfile before regridding
      if ((m_cur_step%m_plot_interval == 0) && m_plot_interval > 0)
        {
          writePlotFile();
        }

      if ((m_cur_step != 0) && (m_cur_step%m_regrid_interval ==0))
        {
          regrid();
        }

      if ((m_cur_step%m_check_interval == 0) && (m_check_interval > 0)
          && (m_cur_step != m_restart_step))
        {
          writeCheckpointFile();
        }

      timeStep(m_dt);
      m_dt = computeDt();
      if (a_max_time - m_time < m_dt) m_dt = a_max_time - m_time;
    }

  // dump out final plotfile, if appropriate
  if (m_plot_interval >= 0)
    {
      writePlotFile();
    }

}

void
amrPHASE::timeStep(Real a_dt)
{

  if (m_verbosity >=2)
    {
      pout() << "Timestep " << m_cur_step
             << " Advancing solution from time "
             << m_time << " with dt = " << a_dt << endl;
    }

  // first swap old and new data; since these
  // are pointers anyway, we can do this w/o doing the copy
  CH_assert (m_old_phi.size() == m_new_phi.size());
  LevelData<FArrayBox>* tempDataPtr;

  // swap old and new states
  for (int lev=0; lev< m_old_phi.size(); lev++)
    {
      tempDataPtr = m_old_phi[lev];
      m_old_phi[lev] = m_new_phi[lev];
      m_new_phi[lev] = tempDataPtr;
    }

  // compute velocity from potential
  computePhaseVelocity(m_CCvelocity, m_potential,
                       m_old_phi,m_time);


  // advect scalar using velocity field
  advectScalar(m_new_phi, m_old_phi,
               m_CCvelocity, m_scalIBC, m_dt);

#if 0
  // average solutions down to coarser levels
  if (m_finest_level != 0)
    {
      for (int lev=m_finest_level; lev>0; lev--)
        {
          // probably will make sense to predefine these and
          // keep them around rather than re-allocating them
          // every timestep
          CoarseAverage avgDown(m_amrGrids[lev],
                                m_num_comps, m_refinement_ratios[lev-1]);
          avgDown.averageToCoarse(*m_new_phi[lev-1], *m_new_phi[lev]);
        }
    }
#endif

  // finally, update to new time and increment current step
  m_dt = a_dt;
  m_time += a_dt;
  m_cur_step += 1;

  pout () << "amrPHASE::timestep "
          << ",      end time = "
          << setiosflags(ios::fixed) << setprecision(6) << setw(12) << m_time
          << ", dt = "
          << setiosflags(ios::fixed) << setprecision(6) << setw(12) << a_dt
          << endl;


  // report max and min of phi
  if (m_verbosity >= 2)
    {
      Real phiMax, phiMin;
      for (int comp=0; comp<m_num_comps; comp++)
        {
          Interval comps(comp,comp);
          phiMax = computeMax(m_new_phi, m_refinement_ratios,
                               comps, 0);

          phiMin = computeMin(m_new_phi, m_refinement_ratios,
                               comps, 0);

          pout() << "   Time = " << m_time << " phi"
                 << comp << ": Min= " << phiMin
                 << ", Max= " << phiMax << endl;

        }
    }




  int totalCellsAdvanced = 0;
  for (int lev=0; lev<m_num_cells.size(); lev++)
    {
      totalCellsAdvanced += m_num_cells[lev];
    }

  pout() << "Time = " << m_time << " cells advanced = "
         << totalCellsAdvanced << endl;
  for (int lev=0; lev<m_num_cells.size(); lev++)
    {
      pout () << "Time = " << m_time
              << "  level " << lev << " cells advanced = "
              << m_num_cells[lev] << endl;
    }


}


void
amrPHASE::computePhaseVelocity(Vector<LevelData<FArrayBox>* > & a_velocity,
                               Vector<LevelData<FArrayBox>* > & a_potential,
                               Vector<LevelData<FArrayBox>* >& a_phi,
                               Real a_dt)
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::computePhaseVelocity" << endl;
    }

  int lBase = 0;
  int numPotentialComps = 1;

  // temp storage for rhs for Laplacian solve
  Vector<LevelData<FArrayBox>* > rhs(m_finest_level+1, NULL);
  Vector<LevelData<FArrayBox>* > gradient(m_finest_level+1, NULL);


  // loop over levels to set up for solve
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      //Real dxLevel = m_amrDx[lev];
      const DisjointBoxLayout& levelGrids = m_amrReducedGrids[lev];

      // allocate temp storage for this level
      rhs[lev] = new LevelData<FArrayBox>(levelGrids,
                                          numPotentialComps);

      // allocate temp storage for this level
      gradient[lev] = new LevelData<FArrayBox>(levelGrids,
                                               numPotentialComps);


    } // end loop over levels

  // compute RHS (if maxLevel > 0, need to rethink this)
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      if (lev > 0) MayDay::Error("Potential solve not implemented for AMR");

      LevelData<FArrayBox>& rhsLev = *(rhs[lev]);
      LevelData<FArrayBox>& phiLev = *(a_phi[lev]);

      // set rhs to 0
      DataIterator dit = rhsLev.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          rhsLev[dit()].setVal(0.0);
        }

      // RHS is phi integrated in the y- (v-) direction
      int transverseDir = 1;
      SumOp op(transverseDir);
      Real Pi = 4.0*atan(1.0);
      Real integrationScale = 4.0*Pi*m_G*m_amrDx[lev];
      op.scale = integrationScale;

      Interval srcInterval(0,0);
      Interval destInterval(0,0);

      if (m_verbosity > 5)
        {
          pout() << "reduction Copier: " << endl;
          m_reduceCopiers[lev]->print();
        }

      phiLev.copyTo(srcInterval, rhsLev, destInterval,
                    *m_reduceCopiers[lev], op);

    }

  if (m_verbosity >= 3)
    {
      pout() << "entering 1d potential solve" << endl;
    }

  // this is where things can go to 1D
  // put this in a separate function
  int status = doOneDimensionalPotentialSolve(rhs,
                                              gradient,
                                              m_amrDomains,
                                              m_amrDx,
                                              m_refinement_ratios,
                                              m_finest_level,
                                              lBase,
                                              numPotentialComps,
                                              m_verbosity);
  if (m_verbosity >= 3)
    {
      pout() << "leaving 1d potential solve" << endl;
    }

  if (status != 0)
    {
      pout() << "1D potential solve returned nonzero error code "
             << status << endl;
    }
  // this is where we can bring things back from 1D

  // so now rhs contains the cell-centered gradient of
  // the potential in the x-direction, which will be the
  // y-direction phase velocity. Spread this back out over
  // the rest of the domain

  if (m_verbosity >= 3)
    {
      pout() << "beginning spreading of gradient..." << endl;
    }

  for (int lev=0; lev<=m_finest_level; lev++)
    {
  if (m_verbosity >= 3)
    {
      pout() << " ... level " << lev << "... " << endl;
    }

      LevelData<FArrayBox>& levelGradient = *gradient[lev];

      int spreadingDir = 1;
      SpreadingOp spreadOp(spreadingDir);
      spreadOp.scale = 1.0;

      LevelData<FArrayBox>& levelVel = *a_velocity[lev];
      Interval srcInterval(0,0);
      // spread into the "v" component of velocity
      Interval destInterval(1,1);

      if (m_verbosity > 5)
        {
          pout() << "spreading Copier: " << endl;

          m_spreadCopiers[lev]->print();
        }

      levelGradient.copyTo(srcInterval, levelVel,
                           destInterval, *m_spreadCopiers[lev],
                           spreadOp);

      levelVel.exchange(destInterval);
    }

  if (m_verbosity >= 3)
    {
      pout() << "end spreading of gradient" << endl;
    }

  // finally, clean up local storage
  for (int lev=0; lev<rhs.size(); lev++)
    {
      if (rhs[lev] != NULL)
        {
          delete rhs[lev];
          rhs[lev] = NULL;
        }

      if (gradient[lev] != NULL)
        {
          delete gradient[lev];
          gradient[lev] = NULL;
        }
    }

}



void
amrPHASE::advectScalar(Vector<LevelData<FArrayBox>* >& a_new_scal,
                       Vector<LevelData<FArrayBox>* >& a_old_scal,
                       Vector<LevelData<FArrayBox>* >&a_CC_vel,
                       OldPhysIBC* a_scalIBC,
                       Real a_dt)
{

  if (m_verbosity >= 3)
    {
      pout() << "amrPHASE::advectScalar " << endl;
    }



  /// need to build grown grids to get be able to do all of
  /// tracing properly
  int numlevels = m_finest_level+1;
  IntVect advect_grow(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));
  IntVect ghostVect = IntVect::Unit;
  // cell-centered velocity
  Vector<LevelData<FArrayBox>* > local_old_vel(numlevels, NULL);
  Vector<LevelData<FluxBox>* > edgeScal(numlevels, NULL);
  Vector<LevelData<FluxBox>* > adv_vel(numlevels, NULL);
  int nScal = a_new_scal[0]->nComp();

  for (int lev=0; lev<numlevels; lev++)
    {

      const DisjointBoxLayout& levelGrids = a_new_scal[lev]->getBoxes();
      local_old_vel[lev] = new LevelData<FArrayBox>(levelGrids,
                                                    SpaceDim, advect_grow);

      edgeScal[lev] = new LevelData<FluxBox>(levelGrids, nScal,
                                             ghostVect);


      adv_vel[lev] = new LevelData<FluxBox>(levelGrids, 1,
                                            ghostVect);

    }

  // m_time contains the time at which the new state is centered
  Real old_time = m_time - m_dt;
  fillCCvelocity(local_old_vel, a_CC_vel, old_time);

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      // set up patchGodunov for this problem
      m_patchGod[lev]->setPhysIBC(a_scalIBC);
      m_patchGod[lev]->setNumVar(nScal);
      m_patchGod[lev]->setCurrentTime(old_time);


      IntVect ghostVect = IntVect::Unit;
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      LevelData<FArrayBox>& newScal = *a_new_scal[lev];
      LevelData<FArrayBox>& oldScal = *a_old_scal[lev];
      LevelData<FArrayBox>& oldVel = *local_old_vel[lev];
      LevelData<FluxBox>& advVel = *adv_vel[lev];
      LevelData<FluxBox>& levelEdgeScal = *edgeScal[lev];
      DataIterator dit = newScal.dataIterator();

      // average cc velocity to faces
      CellToEdge(oldVel, advVel);

      // if level > 0, need to fill in ghost cells for phi
      if (lev > 0)
        {
          MayDay::Error("AMR not yet implemented");
#if 0
          int nPhiComp = newScal.nComp();
          int nRefCrse = m_refinement_ratios[lev-1];
          int interpRadius = oldScal.ghostVect()[0];
          const ProblemDomain& crseDomain = m_amrDomains[lev-1];
          LevelData<FArrayBox>& oldCrseScal = *a_old_scal[lev-1];
          PiecewiseLinearFillPatch patcher(oldScal.getBoxes(),
                                           oldCrseScal.getBoxes(),
                                           nPhiComp,
                                           crseDomain,
                                           nRefCrse,
                                           interpRadius);

          // set time coefficient to zero
          Real timeInterpCoeff = 0.0;
          patcher.fillInterp(oldScal, oldCrseScal, oldCrseScal,
                             timeInterpCoeff, 0, 0, nPhiComp);
#endif
        }

      // call exchange to make sure that ghost cells are in agreement
      oldScal.exchange();

#ifdef CH_MPI
      if (m_verbosity >= 3)
        {
          //pout() << "procID: " << procID() << "tracing scalars, level = "
          // << lev << endl;
        }
#endif

      FArrayBox srcFab;
      // now trace scalars to edges at time n+1/2
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box& thisGrownBox = oldScal[dit()].box();
          FluxBox& thisEdgeScal = levelEdgeScal[dit()];
          FluxBox& thisAdvVel = advVel[dit()];

#ifdef CH_MPI
          if (m_verbosity >= 3)
            {
              //pout() << "procID: " << procID()
              // << "tracing grid " << thisGrownBox << endl;
            }
#endif


          // this is necessary to synchronize with old-version code
          FArrayBox&  thisCellVel = oldVel[dit()];

          FArrayBox& thisOldScal = oldScal[dit()];

          srcFab.resize(thisGrownBox, nScal);
          srcFab.setVal(0.0);

          // new OldPatchGodunov-based approach to advection
          const Box& curBox = levelGrids[dit()];
          m_patchGod[lev]->setCurrentBox(curBox);
          // set cell-centered velocity field
          m_patchGod[lev]->setCellVelPtr(&thisCellVel);
          // set advection velocity field
          m_patchGod[lev]->setAdvVelPtr(&thisAdvVel);

          // set up fluxBox->FArrayBox[SpaceDim] conversion
          FArrayBox faceScalArray[SpaceDim];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              // use aliasing to handle this
              Interval scalInterval = thisEdgeScal[dir].interval();
              faceScalArray[dir].define(scalInterval, thisEdgeScal[dir]);
            }

          // now compute fluxes
          m_patchGod[lev]->computeFluxes(thisOldScal,
                                         faceScalArray, srcFab,
                                         m_dt, curBox);

          for (int dir=0; dir<SpaceDim; dir++)
            {
              //multiply by edge velocity to get flux
              // do componentwise
              for (int comp=0; comp<nScal; comp++)
              {
                thisEdgeScal[dir].mult(thisAdvVel[dir],0,comp,1);
              }
            }

        } // end loop over grids on this level

    } // end loop over levels


  // average fluxes down to coarser levels
  for (int lev=m_finest_level; lev>0; lev--)
    {
      MayDay::Error("flux averaging not implemented");
#if 0
      CoarseAverageEdge averager(m_amrGrids[lev], nScal,
                                 m_refinement_ratios[lev-1]);
      averager.averageToCoarse(*edgeScal[lev-1],
                               *edgeScal[lev]);
#endif
    }


  // compute div(us)
#ifdef CH_MPI
  if (m_verbosity >= 3)
    {
      //pout() << "procID: " << procID() << "taking divergence" << endl;
    }
#endif

  for (int lev=0; lev<numlevels; lev++)
    {
      LevelData<FArrayBox>& newScal = *a_new_scal[lev];
      LevelData<FArrayBox>& oldScal = *a_old_scal[lev];
      levelDivergenceMAC(newScal, *edgeScal[lev], m_amrDx[lev]);

      // now multiply div(us) by -dt and add to old_scal
      DataIterator dit = newScal.dataIterator();

      for (dit.reset(); dit.ok(); ++dit)
        {
          Real mult = -1.0*a_dt;
          newScal[dit()] *= mult;
          newScal[dit()] += oldScal[dit()];
        }


    } // end loop over levels

  // memory leaks bad
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      delete local_old_vel[lev];
      local_old_vel[lev] = NULL;

      delete adv_vel[lev];
      adv_vel[lev] = NULL;

      delete edgeScal[lev];
      edgeScal[lev] = NULL;
    }
  // that should be it!

}



// do regridding
void
amrPHASE::regrid()
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::regrid" << endl;
    }

  // defer implementation of this guy
  MayDay::Error("amrPHASE::regrid not implemented yet");
  // only do any of this if the max level > 0
  if (m_max_level > 0)
    {

      // in this code, lbase is always 0
      int lbase =0;

      // first generate tags
      Vector<IntVectSet> tagVect(m_max_level);
      Vector<IntVectSet> oneDTagVect(m_max_level);
      tagCells(tagVect, oneDTagVect);

      // now generate new boxes; first do regular grids, then
      // do reduced grids
      int top_level = min(m_finest_level, m_max_level-1);
      Vector<Vector<Box> > old_grids(m_finest_level+1);
      Vector<Vector<Box> > new_grids;

      // this is clunky, but i don't know of a better way to turn
      // a DisjointBoxLayout into a Vector<Box>
      for (int lev=0; lev<= m_finest_level; lev++)
        {
          const DisjointBoxLayout& levelDBL = m_amrGrids[lev];
          old_grids[lev].resize(levelDBL.size());
          LayoutIterator lit = levelDBL.layoutIterator();
          int boxIndex = 0;
          for (lit.begin(); lit.ok(); ++lit, ++boxIndex)
            {
              old_grids[lev][boxIndex] = levelDBL[lit()];
            }
        }

      int new_finest_level;
      if (procID() == uniqueProc(SerialTask::compute))
        {
          BRMeshRefine meshrefine(m_amrDomains[0], m_refinement_ratios,
                                  m_fill_ratio, m_block_factor,
                                  m_nesting_radius, m_max_box_size);

          new_finest_level = meshrefine.regrid(new_grids, tagVect,
                                               lbase, top_level,
                                               old_grids);
        }

      broadcast(new_finest_level, uniqueProc(SerialTask::compute));

      //int numLevels = Min(new_finest_level, m_max_level)+1;

      broadcast(new_grids, uniqueProc(SerialTask::compute));

      const IntVect& phiGhostVect = m_new_phi[0]->ghostVect();

      // now loop through levels and redefine if necessary
      for (int lev=lbase+1; lev<= new_finest_level; ++lev)
        {
          int numGridsNew = new_grids[lev].size();
          Vector<int> procIDs(numGridsNew);
          LoadBalance(procIDs, new_grids[lev]);
          const DisjointBoxLayout newDBL(new_grids[lev], procIDs,
                                         m_amrDomains[lev]);
          m_amrGrids[lev] = newDBL;

          // build new storage
          LevelData<FArrayBox>* old_oldDataPtr = m_old_phi[lev];
          LevelData<FArrayBox>* old_newDataPtr = m_new_phi[lev];

          LevelData<FArrayBox>* new_oldDataPtr =
            new LevelData<FArrayBox>(newDBL, m_num_comps, phiGhostVect);
          LevelData<FArrayBox>* new_newDataPtr =
            new LevelData<FArrayBox>(newDBL, m_num_comps, phiGhostVect);

          int potentialComps = m_potential[0]->nComp();
          IntVect potentialGhost = m_potential[0]->ghostVect();

          LevelData<FArrayBox>* new_potentialPtr =
            new LevelData<FArrayBox>(newDBL, potentialComps, potentialGhost);


          int CCvelocityComps = m_CCvelocity[0]->nComp();
          IntVect CCvelocityGhost = m_CCvelocity[0]->ghostVect();

          LevelData<FArrayBox>* new_CCvelocityPtr =
            new LevelData<FArrayBox>(newDBL, CCvelocityComps,
                                     CCvelocityGhost);

          // first fill with interpolated data from coarser level
          // no need to fill old_data, just new_data

#if 0
          // may eventually want to do post-regrid smoothing on this!
          FineInterp interpolator(newDBL, m_num_comps,
                                  m_refinement_ratios[lev-1],
                                  m_amrDomains[lev]);

          interpolator.interpToFine(*new_newDataPtr, *m_new_phi[lev-1]);
#endif
          // now copy old-grid data into new holder
          if (old_newDataPtr != NULL)
            {
              Interval dataComps = new_newDataPtr->interval();
              old_newDataPtr->copyTo(dataComps, *new_newDataPtr, dataComps);

              // can now delete old data
              delete old_oldDataPtr;
              delete old_newDataPtr;
              // no need to copy old potential and velocity data
              delete m_potential[lev];
              delete m_CCvelocity[lev];
            }


          // now copy new holders into multilevel arrays
          m_old_phi[lev] = new_oldDataPtr;
          m_new_phi[lev] = new_newDataPtr;
          m_potential[lev] = new_potentialPtr;
          m_CCvelocity[lev] = new_CCvelocityPtr;

        } // end loop over currently defined levels

      // now ensure that any remaining levels are null pointers
      // (in case of de-refinement)
      for (int lev=new_finest_level+1; lev<m_old_phi.size(); lev++)
        {
          if (m_old_phi[lev] != NULL)
            {
              delete m_old_phi[lev];
              m_old_phi[lev] = NULL;
            }

          if (m_new_phi[lev] != NULL)
            {
              delete m_new_phi[lev];
              m_new_phi[lev] = NULL;
            }

          if (m_potential[lev] != NULL)
            {
              delete m_potential[lev];
              m_potential[lev] = NULL;
            }

          if (m_CCvelocity[lev] != NULL)
            {
              delete m_CCvelocity[lev];
              m_CCvelocity[lev] = NULL;
            }

          DisjointBoxLayout emptyDBL;
          m_amrGrids[lev] = emptyDBL;
        }

      m_finest_level = new_finest_level;



      // set up counter of number of cells
      for (int lev=0; lev<=m_max_level; lev++)
        {
          m_num_cells[lev] = 0;
          if (lev <= m_finest_level)
            {
              const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
              LayoutIterator lit = levelGrids.layoutIterator();
              for (lit.begin(); lit.ok(); ++lit)
                {
                  const Box& thisBox = levelGrids.get(lit());
                  m_num_cells[lev] += thisBox.numPts();
                }
            }
        }


    } // end if max level > 0 in the first place
}



void
amrPHASE::tagCells(Vector<IntVectSet>& a_tags,
                   Vector<IntVectSet>& a_oneDTags)
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::tagCells" << endl;
    }


  int top_level = a_tags.size();
  top_level = min(top_level-1, m_finest_level);
  // loop over levels
  for (int lev=0; lev<=top_level; lev++)
    {
      IntVectSet& levelTags = a_tags[lev];
      IntVectSet& leveloneDTags = a_oneDTags[lev];
      tagCellsLevel(levelTags, leveloneDTags, lev);
    }

}

void
amrPHASE::tagCellsLevel(IntVectSet& a_tags,
                        IntVectSet& a_oneDTags, int a_level)
{

  if (m_verbosity > 4)
    {
      pout() << "amrPHASE::tagCellsLevel " << a_level << endl;
    }



  // base tags on undivided gradient of phi
  // first stab -- don't do BC's; just do one-sided
  // stencils at box edges (hopefully good enough),
  // since doing BC's properly is somewhat expensive.

  DataIterator dit = m_new_phi[a_level]->dataIterator();
  LevelData<FArrayBox>& levelPhi = *m_new_phi[a_level];
  const DisjointBoxLayout& levelGrids = m_amrGrids[a_level];

  // for now, do taggine based on first component
  //int startComp = 0;
  //int numComp = 1;

  // need to ensure that ghost cells are set properly
  levelPhi.exchange(levelPhi.interval());

  IntVectSet local_tags;
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox gradPhi(levelGrids[dit()], m_num_comps);

      for (int dir=0; dir<SpaceDim; dir++)
        {
          const Box b = levelGrids[dit()];
          const Box bcenter = b & grow ( m_amrDomains[a_level], -BASISV(dir) );
          const Box blo = b & adjCellLo( bcenter, dir );
          const Box bhi = b & adjCellHi( bcenter, dir );
          const int haslo = ! blo.isEmpty();
          const int hashi = ! bhi.isEmpty();
          FORT_UNDIVIDEDGRAD ( CHF_FRA(gradPhi),
                               CHF_CONST_FRA(levelPhi[dit()]),
                               CHF_BOX(bcenter),
                               CHF_BOX(blo),
                               CHF_BOX(bhi),
                               CHF_CONST_INT(dir),
                               CHF_CONST_INT(haslo),
                               CHF_CONST_INT(hashi));


          // now tag cells based on values
          BoxIterator bit(levelGrids[dit()]);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (abs(gradPhi(iv)) > m_tagging_val)
                local_tags |= iv;
            } // end loop over cells
        } // end loop over directions

    } // end loop over grids

  // now buffer tags
  local_tags.grow(m_tags_grow);
  local_tags &= m_amrDomains[a_level];

  // gather all tags into one place
  Vector<IntVectSet> all_tags;

  const int dest_proc = uniqueProc (SerialTask::compute);

    gather(all_tags, local_tags, dest_proc);
    if (procID()==uniqueProc (SerialTask::compute))
    {
      for (int i=0; i!=all_tags.size(); ++i)
        a_tags |= all_tags[i];
    }
}

void
amrPHASE::tagCellsInit(Vector<IntVectSet>& a_tags,
                       Vector<IntVectSet>& a_oneDtags)
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::tagCellsInit" << endl;
    }


  // for initial time, refine entire domain for all levels

  if (m_refine_initial_domain)
    {
      for (int lev=0; lev<m_max_level; lev++)
        {
          a_tags[lev].define(m_amrDomains[lev].domainBox());
        }
    }
  else
    {
      tagCells(a_tags, a_oneDtags);
    }

}


void
amrPHASE::initGrids(int a_finest_level)
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::initGrids" << endl;
    }


  m_finest_level = 0;
  // first create base level
  Vector<Box> baseBoxes;
  domainSplit(m_amrDomains[0], baseBoxes, m_max_box_size,
              m_block_factor);

  Vector<int> procAssign(baseBoxes.size());
  LoadBalance(procAssign,baseBoxes);

  DisjointBoxLayout baseGrids(baseBoxes, procAssign, m_amrDomains[0]);

  Vector<Box> reducedBaseBoxes;
#if 0
  // this needs to be a 1D domainSplit
  domainSplit(m_amrReducedDomains[0], reducedBaseBoxes,
              m_max_box_size,
              m_block_factor);
#endif

  // do this by modifying boxes from the baseGrids
  for (int i=0; i<baseBoxes.size(); i++)
    {
      Box thisBox(baseBoxes[i]);
      // if this box touches the lower y bounday, then
      // project it down to the one-D space
      if (thisBox.smallEnd(1) == 0)
        {
          for (int dir=1; dir<SpaceDim; dir++)
            {
              thisBox.setSmall(dir, 0);
              thisBox.setBig(dir, 0);
            }
          reducedBaseBoxes.push_back(thisBox);
        }
    }

  Vector<int> reducedProcAssign(reducedBaseBoxes.size());
  LoadBalance(reducedProcAssign,reducedBaseBoxes);

  DisjointBoxLayout reducedBaseGrids(reducedBaseBoxes, reducedProcAssign,
                                     m_amrReducedDomains[0]);

  m_amrGrids.resize(m_max_level+1);
  m_amrReducedGrids.resize(m_max_level+1);
  m_amrGrids[0] = baseGrids;
  m_amrReducedGrids[0] = reducedBaseGrids;

  int transverseDir = 1;
  m_reduceCopiers.resize(m_max_level+1);
  m_reduceCopiers[0] = new ReductionCopier(baseGrids,
                                           reducedBaseGrids,
                                           m_amrDomains[0],
                                           transverseDir);

  m_spreadCopiers.resize(m_max_level+1);
  m_spreadCopiers[0] = new SpreadingCopier(reducedBaseGrids,
                                           baseGrids,
                                           m_amrDomains[0],
                                           transverseDir);
  IntVect phiGhostVect = ADVECT_GROW*IntVect::Unit;
  m_new_phi[0] = new LevelData<FArrayBox>(baseGrids, m_num_comps,
                                          phiGhostVect);

  m_old_phi[0] = new LevelData<FArrayBox>(baseGrids, m_num_comps,
                                          phiGhostVect);


  IntVect ghostVect = IntVect::Unit;

  int numCCVelComps = SpaceDim;
  m_CCvelocity[0] = new LevelData<FArrayBox>(baseGrids, numCCVelComps,
                                             ghostVect);


  int numPotentialComps = 1;
  m_potential[0] = new LevelData<FArrayBox>(reducedBaseGrids,
                                            numPotentialComps,
                                            ghostVect);

  // initialize base level data
  initData(m_new_phi, m_CCvelocity);

  int numLevels = 1;
  bool moreLevels = (m_max_level > 0);

  int baseLevel = 0;
  //int topLevel = m_finest_level;


  BRMeshRefine meshrefine;
  if (moreLevels)
    {
      meshrefine.define(m_amrDomains[0], m_refinement_ratios,
                        m_fill_ratio, m_block_factor,
                        m_nesting_radius, m_max_box_size);
    }

  Vector<IntVectSet> tagVect(m_max_level);
  Vector<IntVectSet> oneDTagVect(m_max_level);

  Vector<Vector<Box> > oldBoxes(1);
  Vector<Vector<Box> > newBoxes;
  oldBoxes[0] = baseBoxes;
  newBoxes = oldBoxes;
  int new_finest_level = 0;

  while (moreLevels)
    {
      // haven't quite figured out what to do for reduced grids
      // for AMR
      MayDay::Error("AMR not implemented for amrPHASE::initGrids");

      // default is moreLevels = false
      // (only repeat loop in the case where a new level is generated
      // which is still coarser than maxLevel)
      moreLevels = false;
      tagCellsInit(tagVect, oneDTagVect);

      // two possibilities -- need to generate grids
      // level-by-level, or we are refining all the
      // way up for the initial time.  check to
      // see which it is by seeing if the finest-level
      // tags are empty
      if (tagVect[m_max_level-1].isEmpty())
        {
          int top_level = m_finest_level;
          int old_top_level = top_level;
          new_finest_level = meshrefine.regrid(newBoxes,
                                               tagVect, baseLevel,
                                               top_level,
                                               oldBoxes);

          if (new_finest_level > top_level) top_level++;
          oldBoxes = newBoxes;

          // now see if we need another pass through grid generation
          if ((top_level < m_max_level) && (top_level > old_top_level))
            {
              moreLevels = true;
            }

        }
      else
        {

          // for now, define old_grids as just domains
          oldBoxes.resize(m_max_level+1);
          for (int lev=1; lev<=m_max_level; lev++)
            {
              oldBoxes[lev].push_back(m_amrDomains[lev].domainBox());
            }

          int top_level = m_max_level -1;
          new_finest_level = meshrefine.regrid(newBoxes,
                                               tagVect, baseLevel,
                                               top_level,
                                               oldBoxes);
        }


      broadcast(new_finest_level, uniqueProc(SerialTask::compute));

      numLevels = Min(new_finest_level, m_max_level)+1;

      broadcast(newBoxes, uniqueProc(SerialTask::compute));


      // now loop through levels and define
      for (int lev=baseLevel+1; lev<= new_finest_level; ++lev)
        {
          int numGridsNew = newBoxes[lev].size();
          Vector<int> procIDs(numGridsNew);
          LoadBalance(procIDs, newBoxes[lev]);
          const DisjointBoxLayout newDBL(newBoxes[lev], procIDs,
                                         m_amrDomains[lev]);
          m_amrGrids[lev] = newDBL;

          // build new storage, clearing previous storage if necessary
          if (m_new_phi[lev] != NULL)
            {
              delete m_new_phi[lev];
              m_new_phi[lev] = NULL;
            }
          m_new_phi[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                                    m_num_comps,
                                                    phiGhostVect);

          if (m_old_phi[lev] != NULL)
            {
              delete m_old_phi[lev];
              m_old_phi[lev] = NULL;
            }
          m_old_phi[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                                    m_num_comps,
                                                    phiGhostVect);

          if (m_potential[lev] != NULL)
            {
              delete m_potential[lev];
              m_potential[lev] = NULL;
            }
          m_potential[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                                     numPotentialComps,
                                                     ghostVect);


          if (m_CCvelocity[lev] != NULL)
            {
              delete m_CCvelocity[lev];
              m_CCvelocity[lev] = NULL;
            }
          m_CCvelocity[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                                       numCCVelComps,
                                                       ghostVect);
        }

      m_finest_level = new_finest_level;

      // finally, initialize data on final hierarchy
      // only do this if we've created new levels
      if (m_finest_level > 0)
        {
          initData(m_new_phi, m_CCvelocity);
        }
    } // end while more levels to do
}


void
amrPHASE::setupFixedGrids(const std::string& a_gridFile)
{
  Vector<Vector<Box> > gridvect;

  MayDay::Error("amrPHASE::setupFixedGrids not implemented");

  if (procID() == uniqueProc(SerialTask::compute))
    {
      gridvect.push_back(Vector<Box>(1,m_amrDomains[0].domainBox()));

      // read in predefined grids
      ifstream is(a_gridFile.c_str(), ios::in);

      if (is.fail())
        {
          MayDay::Error("Cannot open grids file");
        }

      // format of file:
      //   number of levels, then for each level (starting with level 1):
      //   number of grids on level, list of boxes
      int inNumLevels;
      is >> inNumLevels;

      CH_assert (inNumLevels <= m_max_level+1);

      if (m_verbosity > 3)
        {
          pout() << "numLevels = " << inNumLevels << endl;
        }

      while (is.get() != '\n');

      gridvect.resize(inNumLevels);

      // check to see if coarsest level needs to be broken up
      domainSplit(m_amrDomains[0],gridvect[0], m_max_box_size,
                  m_block_factor);

      if (m_verbosity >= 3)
        {
          pout() << "level 0: ";
          for (int n=0; n < gridvect[0].size(); n++)
            {
              pout() << gridvect[0][n] << endl;
            }
        }

      // now loop over levels, starting with level 1
      int numGrids = 0;
      for (int lev=1; lev<inNumLevels; lev++)
        {
          is >> numGrids;

          if (m_verbosity >= 3)
            {
              pout() << "level " << lev << " numGrids = "
                     << numGrids <<  endl;
              pout() << "Grids: ";
            }

          while (is.get() != '\n');

          gridvect[lev].resize(numGrids);

          for (int i=0; i<numGrids; i++)
            {
              Box bx;
              is >> bx;

              while (is.get() != '\n');

              // quick check on box size
              Box bxRef(bx);

              if (bxRef.longside() > m_max_box_size)
                {
                  pout() << "Grid " << bx << " too large" << endl;
                  MayDay::Error();
                }

              if (m_verbosity >= 3)
                {
                  pout() << bx << endl;
                }

              gridvect[lev][i] = bx;
            } // end loop over boxes on this level
        } // end loop over levels
    } // end if serial proc

  // broadcast results
  broadcast(gridvect, uniqueProc(SerialTask::compute));

  // now create disjointBoxLayouts and allocate grids

  m_amrGrids.resize(m_max_level+1);
  for (int lev=0; lev<gridvect.size(); lev++)
    {
      int numGridsLev = gridvect[lev].size();
      Vector<int> procIDs(numGridsLev);
      LoadBalance(procIDs, gridvect[lev]);
      const DisjointBoxLayout newDBL(gridvect[lev],
                                     procIDs,
                                     m_amrDomains[lev]);

      m_amrGrids[lev] = newDBL;

      // build storage for this level
      IntVect ghostVect(IntVect::Unit);
      m_new_phi[lev] = new LevelData<FArrayBox>(newDBL,
                                                m_num_comps,
                                                ghostVect);

      m_old_phi[lev] = new LevelData<FArrayBox>(newDBL,
                                                m_num_comps,
                                                ghostVect);

      int numPotentialComps = 1;
      m_potential[lev] = new LevelData<FArrayBox>(newDBL,
                                                 numPotentialComps,
                                                 ghostVect);


      int numCCVelComps = SpaceDim;
      m_CCvelocity[lev] = new LevelData<FArrayBox>(newDBL,
                                                   numCCVelComps,
                                                   ghostVect);


    }

  // finally set finest level and initialize data on hierarchy
  m_finest_level = gridvect.size() -1;

  initData(m_new_phi, m_CCvelocity);

}






void
amrPHASE::initData(Vector<LevelData<FArrayBox>* >& a_phi,
                   Vector<LevelData<FArrayBox>* >& a_ccVelocity)
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::initData" << endl;
    }

  for (int lev=0; lev<=m_finest_level; lev++)
    {
      LevelData<FArrayBox>& levelPhi = *(a_phi[lev]);
      LevelData<FArrayBox>& levelVel = *(a_ccVelocity[lev]);

      DataIterator levelDit = levelPhi.dataIterator();
      for (levelDit.begin(); levelDit.ok(); ++levelDit)
        {
          FArrayBox& thisPhi = levelPhi[levelDit()];
          FArrayBox& thisVel = levelVel[levelDit()];
          FORT_SETPHI(CHF_FRA(thisPhi),
                      CHF_BOX(thisPhi.box()),
                      CHF_CONST_REAL(m_amrDx[lev]),
                      CHF_CONST_REALVECT(m_domainSize),
                      CHF_CONST_REALVECT(m_origin));

          FORT_SETVEL(CHF_FRA(thisVel),
                      CHF_BOX(thisVel.box()),
                      CHF_CONST_REAL(m_amrDx[lev]),
                      CHF_CONST_REALVECT(m_domainSize),
                      CHF_CONST_REALVECT(m_origin));
        }
    }

#if 0
  // may be necessary to average down here
  for (int lev=m_finest_level; lev>0; lev--)
    {
      CoarseAverage avgDown(m_amrGrids[lev],
                            m_num_comps, m_refinement_ratios[lev-1]);
      avgDown.averageToCoarse(*m_new_phi[lev-1], *m_new_phi[lev]);
    }
#endif

  // finally, need to initialize potential and velocity
  computePhaseVelocity(a_ccVelocity, m_potential, a_phi, m_time);

}

// compute timestep
Real
amrPHASE::computeDt()
{
  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::computeDt" << endl;
    }


  Real dt = 1.0e8;
  for (int lev=0; lev<=m_finest_level; lev++)
    {
      Real dtLev = dt;
      LevelData<FArrayBox>& levelVel = *m_CCvelocity[lev];
      Real maxVelLev = -1.0;

#if 0
      const DisjointBoxLayout& levelGrids = m_amrGrids[lev];
      DataIterator levelDit = levelVel.dataIterator();
      for (levelDit.reset(); levelDit.ok(); ++levelDit)
        {
          const Box& thisBox = levelGrids[levelDit()];
          int normType = 0;
          Real maxVelGrid = levelVel[levelDit].norm(thisBox,
                                                    normType,
                                                    0,SpaceDim);

          if (maxVelGrid > maxVelLev) maxVelLev = maxVelGrid;
        }
#endif
      Interval velComps(0,SpaceDim-1);
      maxVelLev = norm(levelVel, velComps, 0);
      dtLev = m_amrDx[lev]/maxVelLev;
      dt = min(dt, dtLev);
    }

  if (m_cur_step == 0)
    {
      dt *= m_initial_cfl;
    }
  else
    {
      dt *= m_cfl;
    }

  // also check to see if max grow rate applies
  // (m_dt > 0 test screens out initial time, when we set m_dt to a negative
  // number by default)
  if (dt > m_max_dt_grow*m_dt && (m_dt > 0) )
    dt = m_max_dt_grow*m_dt;

  return dt;

}

Real
amrPHASE::computeInitialDt()
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::computeInitialDt" << endl;
    }


  // for now, just call computeDt;
  Real dt = computeDt();
  return dt;
}


#ifdef CH_USE_HDF5

  /// write hdf5 plotfile to the standard location
void
amrPHASE::writePlotFile()
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::writePlotFile" << endl;
    }

  // write scalar + velocity
  int numPlotComps = m_num_comps + SpaceDim;



  // generate data names
  string phiName("phi");

  Vector<string> vectName(numPlotComps);
  int plotComp = 0;
  for (int comp=0; comp<m_num_comps; comp++)
    {
      char idx[4]; sprintf(idx, "%d", comp);
      vectName[plotComp] = phiName+string(idx);
      plotComp++;
    }

  D_TERM(vectName[plotComp] = "xVel";,
        vectName[plotComp+1] = "yVel";,
        vectName[plotComp+2] = "zVel";)

  Box domain = m_amrDomains[0].domainBox();
  Real dt = 1.;
  int numLevels = m_finest_level +1;

  // compute plot data
  Vector<LevelData<FArrayBox>* > plotData(m_new_phi.size(), NULL);
  for (int lev=0; lev<numLevels; lev++)
    {
      // first allocate storage
      // ghost vect makes things simpler
      IntVect ghostVect(IntVect::Unit);
      plotData[lev] = new LevelData<FArrayBox>(m_amrGrids[lev],
                                               numPlotComps,
                                               ghostVect);
      // now copy new-time solution into plotData
      Interval srcComps = m_new_phi[lev]->interval();
      m_new_phi[lev]->copyTo(srcComps, *(plotData[lev]), srcComps);

      int numPhi = m_new_phi[lev]->nComp();
      Interval destInterval(numPhi, numPhi+SpaceDim-1);
      m_CCvelocity[lev]->copyTo(m_CCvelocity[lev]->interval(),
                                *(plotData[lev]), destInterval);

    } // end loop over levels for computing plot data

  // generate plotfile name
  char iter_str[100];
  if (m_cur_step < 10)
    {
      sprintf(iter_str, "%s0000%d.%dd.hdf5", m_plot_prefix.c_str(),
              m_cur_step, SpaceDim);
    }
  else if (m_cur_step < 100)
    {
      sprintf(iter_str, "%s000%d.%dd.hdf5", m_plot_prefix.c_str(),
              m_cur_step, SpaceDim);
    }
  else if (m_cur_step < 1000)
    {
      sprintf(iter_str, "%s00%d.%dd.hdf5", m_plot_prefix.c_str(),
              m_cur_step, SpaceDim);
    }
  else if (m_cur_step < 10000)
    {
      sprintf(iter_str, "%s0%d.%dd.hdf5", m_plot_prefix.c_str(),
              m_cur_step, SpaceDim);
    }
  else
    {
      sprintf(iter_str, "%s%d.%dd.hdf5", m_plot_prefix.c_str(),
              m_cur_step, SpaceDim);
    }

  string filename(iter_str);
  WriteAMRHierarchyHDF5(filename, m_amrGrids, plotData, vectName,
                        domain, m_amrDx[0], dt, m_time, m_refinement_ratios,
                        numLevels);

  // need to delete plotData
  for (int lev=0; lev<numLevels; lev++)
    {
      if (plotData[lev] != NULL)
        {
          delete plotData[lev];
          plotData[lev] = NULL;
        }
    }
}

  /// write checkpoint file out for later restarting
void
amrPHASE::writeCheckpointFile() const
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::writeCheckpointfile" << endl;
    }

#ifdef CH_USE_HDF5

  string phiName("phi");
  Vector<string> vectName(m_num_comps);
  for (int comp=0; comp<m_num_comps; comp++)
    {
      char idx[4]; sprintf(idx, "%d", comp);
      vectName[comp] = phiName+string(idx);
    }
  Box domain = m_amrDomains[0].domainBox();
  //int numLevels = m_finest_level +1;

  // generate checkpointfile name
  char (iter_str[100]);
#if 0
  if (m_cur_step < 10)
    {
      sprintf(iter_str, "%s000%d.%dd.hdf5", m_check_prefix.c_str(),
              m_cur_step, SpaceDim);
    }
  else if (m_cur_step < 100)
    {
      sprintf(iter_str, "%s00%d.%dd.hdf5", m_check_prefix.c_str(),
              m_cur_step, SpaceDim);
    }
  else if (m_cur_step < 1000)
    {
      sprintf(iter_str, "%s0%d.%dd.hdf5", m_check_prefix.c_str(),
              m_cur_step, SpaceDim);
    }
  else
    {
      sprintf(iter_str, "%s%d.%dd.hdf5", m_check_prefix.c_str(),
              m_cur_step, SpaceDim);
    }
  #endif
  // overwrite the same checkpoint file, rather than re-writing them
  sprintf(iter_str, "%s.%d.hdf5", m_check_prefix.c_str(), SpaceDim);

  if (m_verbosity > 3)
    {
      pout() << "checkpoint file name = " << iter_str << endl;
    }

  HDF5Handle handle(iter_str, HDF5Handle::CREATE);

  // write amr data -- only dump out things which are essential
  // to restarting the computation (i.e. max_level, finest_level,
  // time, refinement ratios, etc.).  Other paramters (regrid
  // intervals, block-factor, etc can be changed by the inputs
  // file of the new run.
  // At the moment, the maximum level is not allowed to change,
  // although in principle, there is no real reason why it couldn't
  //
  HDF5HeaderData header;
  header.m_int["max_level"] = m_max_level;
  header.m_int["finest_level"] = m_finest_level;
  header.m_int["current_step"] = m_cur_step;
  header.m_real["time"] = m_time;
  header.m_real["dt"] = m_dt;
  header.m_int["num_comps"] = m_num_comps;
  // at the moment, save cfl, but it can be changed by the inputs
  // file if desired.
  header.m_real["cfl"] = m_cfl;

  // periodicity info
  D_TERM(
         if (m_amrDomains[0].isPeriodic(0))
         header.m_int["is_periodic_0"] = 1;
         else
         header.m_int["is_periodic_0"] = 0; ,

         if (m_amrDomains[0].isPeriodic(1))
         header.m_int["is_periodic_1"] = 1;
         else
         header.m_int["is_periodic_1"] = 0; ,

         if (m_amrDomains[0].isPeriodic(2))
         header.m_int["is_periodic_2"] = 1;
         else
         header.m_int["is_periodic_2"] = 0;
         );


  // set up component names
  char compStr[30];
  //string phiName("phi");
  string compName;
  for (int comp=0; comp < m_num_comps; comp++)
    {
      // first generate component name
      char idx[4]; sprintf(idx, "%d", comp);
      compName = phiName + string(idx);
      sprintf(compStr, "component_%d", comp);
      header.m_string[compStr] = compName;
    }

  header.writeToFile(handle);

  // now loop over levels and write out each level's data
  // note that we loop over all allowed levels, even if they
  // are not defined at the moment.
  for (int lev=0; lev<= m_max_level; lev++)
    {
      // set up the level string
      char levelStr[20];
      sprintf(levelStr, "%d", lev);
      const std::string label = std::string("level_") + levelStr;

      handle.setGroup(label);

      // set up the header info
      HDF5HeaderData levelHeader;
      if (lev < m_max_level)
        {
          levelHeader.m_int["ref_ratio"] = m_refinement_ratios[lev];
        }
      levelHeader.m_real["dx"] = m_amrDx[lev];
      levelHeader.m_box["prob_domain"] = m_amrDomains[lev].domainBox();

      levelHeader.writeToFile(handle);

      // now write the data for this level
      // only try to write data if level is defined.
      if (lev <= m_finest_level)
        {
          write(handle, m_amrGrids[lev]);
          write(handle, *(m_new_phi[lev]), "phiData");
        }
    }// end loop over levels

  handle.close();
#endif
}


/// read checkpoint file for restart
void
amrPHASE::readCheckpointFile(HDF5Handle& a_handle)
{

  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::readCheckpointFile" << endl;
    }

#ifndef CH_USE_HDF5
  MayDay::Error("code must be compiled with HDF5 to read checkpoint files");
#endif

#ifdef CH_USE_HDF5
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (m_verbosity >= 3)
    {
      pout() << "hdf5 header data: " << endl;
      pout() << header << endl;
    }

  // read max level
  if (header.m_int.find("max_level") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain max_level");
    }
  // we can change max level upon restart
  int max_level_check = header.m_int["max_level"];
  if (max_level_check != m_max_level)
    {
      if (m_verbosity > 0)
        {
          pout() << "Restart file has a different max level than inputs file"
                 << endl;
          pout() << "     max level from inputs file = "
                 << m_max_level << endl;
          pout() << "     max level in checkpoint file = "
                 << max_level_check << endl;
          pout() << "Using max level from inputs file" << endl;
        }
    }
  // read finest level
  if (header.m_int.find("finest_level") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain finest_level");
    }

  m_finest_level = header.m_int["finest_level"];
  if (m_finest_level > m_max_level)
    {
      MayDay::Error("finest level in restart file > max allowable level!");
    }

  // read current step
  if (header.m_int.find("current_step") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain current_step");
    }

  m_cur_step = header.m_int["current_step"];
  m_restart_step = m_cur_step;

  // read time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("checkpoint file does not contain time");
    }

  m_time = header.m_real["time"];

  // read timestep
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("checkpoint file does not contain dt");
    }

  m_dt = header.m_real["dt"];

  // read num comps
  if (header.m_int.find("num_comps") == header.m_int.end())
    {
      MayDay::Error("checkpoint file does not contain num_comps");
    }

  m_num_comps = header.m_int["num_comps"];

  // read cfl
  if (header.m_real.find("cfl") == header.m_real.end())
    {
      MayDay::Error("checkpoint file does not contain cfl");
    }

  Real check_cfl = header.m_real["cfl"];
  ParmParse ppCheck("phase");

    if (ppCheck.contains("cfl"))
      {
        // check for consistency and warn if different
        if (check_cfl != m_cfl)
          {
            if (m_verbosity > 0)
              {
                pout() << "CFL in checkpoint file different from inputs file"
                       << endl;
                pout() << "     cfl in inputs file = " << m_cfl << endl;
                pout() << "     cfl in checkpoint file = " << check_cfl
                       << endl;
                pout() << "Using cfl from inputs file" << endl;
              }
          }  // end if cfl numbers differ
      } // end if cfl present in inputs file
    else
      {
        m_cfl = check_cfl;
      }

  // read periodicity info
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

  // now resize stuff
  m_amrDomains.resize(m_max_level+1);
  m_amrGrids.resize(m_max_level+1);
  m_amrDx.resize(m_max_level+1);
  m_old_phi.resize(m_max_level+1, NULL);
  m_new_phi.resize(m_max_level+1, NULL);

  // now read in level-by-level data
  for (int lev=0; lev<= m_max_level; lev++)
    {
      // set up the level string
      char levelStr[20];
      sprintf(levelStr, "%d", lev);
      const std::string label = std::string("level_") + levelStr;

      a_handle.setGroup(label);

      // read header info
      HDF5HeaderData header;
      header.readFromFile(a_handle);

      if (m_verbosity >= 3)
        {
          pout() << "level " << lev << " header data" << endl;
          pout() << header << endl;
        }

  // Get the refinement ratio
      if (lev < m_max_level)
        {
          int checkRefRatio;
          if (header.m_int.find("ref_ratio") == header.m_int.end())
            {
              MayDay::Error("checkpoint file does not contain ref_ratio");
            }
          checkRefRatio = header.m_int["ref_ratio"];

          // check for consistency
          if (checkRefRatio != m_refinement_ratios[lev])
            {
              MayDay::Error("inputs file and checkpoint file ref ratios inconsistent");
            }
        }

      // read dx
      if (header.m_real.find("dx") == header.m_real.end())
        {
          MayDay::Error("checkpoint file does not contain dx");
        }

      m_amrDx[lev] = header.m_real["dx"];

      // read problem domain box
      if (header.m_box.find("prob_domain") == header.m_box.end())
        {
          MayDay::Error("checkpoint file does not contain prob_domain");
        }
      Box domainBox = header.m_box["prob_domain"];

      m_amrDomains[lev] = ProblemDomain(domainBox, isPeriodic);


      // the rest is only applicable if this level is defined
      if (lev <= m_finest_level)
        {
          // read grids
          Vector<Box> grids;
          const int grid_status = read(a_handle, grids);
          if (grid_status != 0)
            {
              MayDay::Error("checkpoint file does not contain a Vector<Box>");
            }
          // do load balancing
          int numGrids = grids.size();
          Vector<int> procIDs(numGrids);
          LoadBalance(procIDs, grids);
          DisjointBoxLayout levelDBL(grids, procIDs, m_amrDomains[lev]);
          m_amrGrids[lev] = levelDBL;

          // allocate this level's storage
          IntVect ghostVect(IntVect::Unit);
          m_old_phi[lev] = new LevelData<FArrayBox>(levelDBL, m_num_comps,
                                                    ghostVect);
          m_new_phi[lev] = new LevelData<FArrayBox>(levelDBL, m_num_comps,
                                                    ghostVect);

          // read this level's data
          LevelData<FArrayBox>& new_phi = *m_new_phi[lev];
          const int dataStatus = read<FArrayBox>(a_handle,
                                                 new_phi,
                                                 "phiData",
                                                 levelDBL);

          if (dataStatus != 0)
            {
              MayDay::Error("checkpoint file does not contain phi data");
            }
        } // end if this level is defined
    } // end loop over levels


  // do we need to close the handle?
#endif

}

/// set up for restart
void
amrPHASE::restart(string& a_restart_file)
{
  if (m_verbosity > 3)
    {
      pout() << "amrPHASE::restart" << endl;
    }

  HDF5Handle handle(a_restart_file, HDF5Handle::OPEN_RDONLY);
  // first read in data from checkpoint file
  readCheckpointFile(handle);

  // don't think I need to do anything else, do I?


}

#endif

/// fill cell-centered velocity field
void
amrPHASE::fillCCvelocity(Vector<LevelData<FArrayBox>* >& a_fillVel,
                        Vector<LevelData<FArrayBox>* >& a_originalVel,
                        Real a_oldTime)
{
  // fill velocity by copying from valid regions, and by
  // interpolating from coarse level where this level's data is
  // not available
  for (int lev=0; lev<= m_finest_level; lev++)
    {
      LevelData<FArrayBox>& destVel = *(a_fillVel[lev]);
      LevelData<FArrayBox>& srcVel = *(a_originalVel[lev]);

      srcVel.copyTo(srcVel.interval(), destVel, destVel.interval());


      if (lev > 0)
        {
          MayDay::Error("AMR not yet implemented");
#if 0
          const LevelData<FArrayBox>& crseVel = *(a_originalVel[lev-1]);
          const ProblemDomain& crseDomain = m_amrDomains[lev-1];
          int nRefCrse = m_refinement_ratios[lev-1];
          // fill in ghost cells by interpolation from coarse data
          int interpRadius = destVel.ghostVect()[0];
          int nComp = destVel.nComp();
          PiecewiseLinearFillPatch patcher(destVel.getBoxes(),
                                           crseVel.getBoxes(),
                                           nComp,
                                           crseDomain,
                                           nRefCrse,
                                           interpRadius);

          // set time coefficient to one
          Real timeInterpCoeff = 1.0;
          patcher.fillInterp(destVel, crseVel, crseVel,
                             timeInterpCoeff, 0, 0, nComp);
#endif
        }



    }

}

