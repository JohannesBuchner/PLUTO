#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CCProjector.H"
#include "Gradient.H"
#include "GradientF_F.H"
#include "Divergence.H"
#include "AMRPoissonOp.H"
#include "RelaxSolver.H"
#include "CoarseAverageEdge.H"
#include "CornerCopier.H"
#include "AMRIO.H"
#include "computeSum.H"
#include "SetValLevel.H"
#include "CoarseAverage.H"

#ifdef CH_USE_HDF5
#include "CH_HDF5.H"
#endif

/// define static variables here
bool CCProjector::m_doSyncProjection = true;
bool CCProjector::m_applySyncCorrection = true;
bool CCProjector::s_applyVDCorrection = true;
bool CCProjector::m_doQuadInterp = true;
Real CCProjector::m_etaLambda = 0.9;

#if defined(CH_USE_DOUBLE)
Real CCProjector::s_solver_tolerance = 1.0e-10; // pmc, 18 oct 2007:  was -10
#elif defined(CH_USE_FLOAT)
Real CCProjector::s_solver_tolerance = 1.0e-5; // pmc, 18 oct 2007:  was -10
#endif

int  CCProjector::s_num_smooth_up = 4;
int  CCProjector::s_num_smooth_down = 4;
int  CCProjector::s_num_precond_smooth = 4;
bool CCProjector::s_constantLambdaScaling = false;
Real CCProjector::s_lambda_timestep = 0.0;
bool CCProjector::pp_init = false;
int  CCProjector::s_verbosity = 0;

/// first define quick-n-easy access functions

// ------------------------------------------------------
int CCProjector::getLevel() const
{
  return m_level;
}

// ------------------------------------------------------
const ProblemDomain& CCProjector::dProblem() const
{
  return m_domain;
}

// ------------------------------------------------------
const DisjointBoxLayout& CCProjector::getBoxes() const
{
  return m_Pi.getBoxes();
}

// ------------------------------------------------------
int CCProjector::nRefCrse() const
{
  return m_nRefCrse;
}

// ------------------------------------------------------
Real CCProjector::dx() const
{
  return m_dx;
}

// ------------------------------------------------------
Real CCProjector::etaLambda() const
{
  return m_etaLambda;
}

// ------------------------------------------------------
CCProjector* CCProjector::fineProjPtr() const
{
  return m_fineProjPtr;
}

// ------------------------------------------------------
CCProjector* CCProjector::crseProjPtr() const
{
  return m_crseProjPtr;
}

// ------------------------------------------------------
LevelData<FArrayBox>& CCProjector::phi()
{
  return m_phi;
}

// ------------------------------------------------------
const LevelData<FArrayBox>& CCProjector::phi() const
{
  return m_phi;
}

// ------------------------------------------------------
LevelData<FArrayBox>& CCProjector::eSync()
{
  return m_eSync;
}

// ------------------------------------------------------
const LevelData<FArrayBox>& CCProjector::eSync() const
{
  return m_eSync;
}

// ------------------------------------------------------
LevelData<FArrayBox>& CCProjector::Pi()
{
  return m_Pi;
}

// ------------------------------------------------------
LevelData<FArrayBox>& CCProjector::eLambda()
{
  return m_eLambda;
}

// ------------------------------------------------------
const LevelData<FArrayBox>& CCProjector::eLambda() const
{
  return m_eLambda;
}

// ------------------------------------------------------
LevelData<FluxBox>& CCProjector::grad_eLambda()
{
  return m_grad_eLambda;
}

// ------------------------------------------------------
const LevelData<FluxBox>& CCProjector::grad_eLambda() const
{
  return m_grad_eLambda;
}

// ------------------------------------------------------
bool CCProjector::doSyncProjection() const
{
  return m_doSyncProjection;
}

// ------------------------------------------------------
bool CCProjector::doMacSync() const
{
  return (m_etaLambda > 0.0);
}

// ------------------------------------------------------
bool CCProjector::isInitialized() const
{
  return m_isInitialized;
}

// ------------------------------------------------------
bool CCProjector::doQuadInterp() const
{
  return m_doQuadInterp;
}

// ------------------------------------------------------
QuadCFInterp& CCProjector::quadCFInterpolator()
{
  return m_cfInterp;
}

/// returns predefined intvectset which shows coverage
const LayoutData<IntVectSet>& CCProjector::getGridIVS()
{
  return m_gradIVS;
}

// ------------------------------------------------------
bool CCProjector::isFinestLevel() const
{
  return m_finest_level;
}

// ------------------------------------------------------
void CCProjector::isFinestLevel(bool a_finest_level)
{
  m_finest_level = a_finest_level;
}

// ------------------------------------------------------
void CCProjector::verbosity(int a_verbosity)
{
  s_verbosity = a_verbosity;
}

// ------------------------------------------------------
int CCProjector::verbosity() const
{
  return s_verbosity;
}

// ------------------------------------------------------
void CCProjector::setPhysBC(const PhysBCUtil& a_bc)
{
  if (m_physBCPtr != NULL)
    {
      delete m_physBCPtr;
      m_physBCPtr = NULL;
    }
  m_physBCPtr = a_bc.newPhysBCUtil();
}

// ------------------------------------------------------
PhysBCUtil* CCProjector::getPhysBCPtr() const
{
  return m_physBCPtr;
}

// ------------------------------------------------------
// default constructor
CCProjector::CCProjector()
{
  m_level = -1;
  m_nRefCrse = -1;
  m_isInitialized = false;
  // have to initialize this to _something_!
  m_finest_level = false;
  m_physBCPtr = NULL;
  m_limitSolverCoarsening = false;
}

// ------------------------------------------------------
// destructor
CCProjector::~CCProjector()
{
  if (m_physBCPtr != NULL)
    {
      delete m_physBCPtr;
      m_physBCPtr = NULL;
    }

  // everything else should be automatic here
}

// ------------------------------------------------------
CCProjector::CCProjector(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_crseGridsPtr,
                         const Box& a_domain,
                         const Real a_dx,
                         CCProjector* a_finerProjPtr,
                         CCProjector* a_crseProjPtr,
                         int a_nRefCrse,
                         int a_level,
                         const PhysBCUtil& a_physBC)
{
  m_physBCPtr = NULL;
  ProblemDomain physdomain(a_domain);
  define(a_grids, a_crseGridsPtr, physdomain, a_dx, a_finerProjPtr,
         a_crseProjPtr, a_nRefCrse, a_level, a_physBC);
}

// ------------------------------------------------------
CCProjector::CCProjector(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_crseGridsPtr,
                         const ProblemDomain& a_domain,
                         const Real a_dx,
                         CCProjector* a_finerProjPtr,
                         CCProjector* a_crseProjPtr,
                         int a_nRefCrse,
                         int a_level,
                         const PhysBCUtil& a_physBC)
{
  m_physBCPtr = NULL;
  define(a_grids, a_crseGridsPtr, a_domain, a_dx, a_finerProjPtr,
         a_crseProjPtr, a_nRefCrse, a_level, a_physBC);
}

// ------------------------------------------------------
void CCProjector::define(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_crseGridsPtr,
                         const Box& a_domain,
                         const Real a_dx,
                         CCProjector* a_finerProjPtr,
                         CCProjector* a_crseProjPtr,
                         int a_nRefCrse,
                         int a_level,
                         const PhysBCUtil& a_physBC)
{
  ProblemDomain physdomain(a_domain);
  define(a_grids, a_crseGridsPtr, physdomain, a_dx, a_finerProjPtr,
         a_crseProjPtr, a_nRefCrse, a_level, a_physBC);
}

// ------------------------------------------------------
void CCProjector::define(const DisjointBoxLayout& a_grids,
                         const DisjointBoxLayout* a_crseGridsPtr,
                         const ProblemDomain& a_domain,
                         const Real a_dx,
                         CCProjector* a_finerProjPtr,
                         CCProjector* a_crseProjPtr,
                         int a_nRefCrse,
                         int a_level,
                         const PhysBCUtil& a_physBC)
{
  // set physical boundary condition object
  setPhysBC(a_physBC);

  if (!pp_init)
    {
      variableSetUp();
    }

  CH_assert (a_level > -1);
  m_level = a_level;
  m_dx = a_dx;
  m_domain = a_domain;

  m_fineProjPtr = a_finerProjPtr;

  m_crseProjPtr = a_crseProjPtr;
  if (m_crseProjPtr != NULL)
    {
      m_nRefCrse = a_nRefCrse;
      m_crseProjPtr->setFineProj(this);
      m_crseProjPtr->isFinestLevel(false);
      CH_assert (a_crseGridsPtr != NULL);

      m_cfInterp.define(a_grids, a_crseGridsPtr, a_dx, a_nRefCrse,
                        1,a_domain);
    }

  int nGrow = 1;
  Gradient::createGridIVS(m_gradIVS, a_grids, nGrow);

  IntVect ghostVect(IntVect::Unit);

  // define persistent storage
  m_Pi.define(a_grids,1,ghostVect);
  m_phi.define(a_grids,1,ghostVect);
  m_eLambda.define(a_grids,1,ghostVect);
  m_eSync.define(a_grids,1,ghostVect);

  m_grad_eLambda.define(a_grids,1);

  // set grad_e_lambda to be 0
  setValLevel(m_grad_eLambda, 0.0);
  // also initialize all other fields to a bogus value
  setValLevel(m_Pi, 1.0e5);
  setValLevel(m_phi, 1.0e5);
  // these need to start as 0 for boundary conditions
  // for first intermediate solve if more than one level
  // of refinement
  setValLevel(m_eLambda, 0.0);
  setValLevel(m_eSync, 0.0);

  defineSolverMGlevel(a_grids, a_crseGridsPtr);

  m_isInitialized = false;
}

// -------------------------------------------------------
void CCProjector::init(const CCProjector& a_oldProj)
{
  // this function performs a partial initialization
  // of this object using the oldProj...

  // at the current state of this algorithm,
  // this doesn't need to do anything (everything will
  // need to be re-initialized anyway)
}

void CCProjector::variableSetUp()
{
  // first set up parm parse object
  ParmParse pp("projection");

  int tempBool;

  tempBool = (int) m_doSyncProjection;
  pp.query("doSyncProjection",tempBool);
  m_doSyncProjection = (tempBool == 1);

  tempBool = (int) m_applySyncCorrection;
  pp.query("applySyncCorrection", tempBool);
  m_applySyncCorrection = (tempBool == 1);

  tempBool = (int) s_applyVDCorrection;
  pp.query("applyFreestreamCorrection", tempBool);
  s_applyVDCorrection = (tempBool == 1);

  tempBool = (int) m_doQuadInterp;
  pp.query("doQuadInterp", tempBool);
  m_doQuadInterp = (tempBool == 1);

  pp.query("eta", m_etaLambda);

  pp.query("solverTol", s_solver_tolerance);
  pp.query("numSmoothUp", s_num_smooth_up);
  pp.query("numSmoothDown", s_num_smooth_down);
  pp.query("numPrecondSmooth", s_num_precond_smooth);

  tempBool = (int) s_constantLambdaScaling;
  pp.query("constantLambdaScaling", tempBool);
  s_constantLambdaScaling = (tempBool == 1);

  // now dump this all back out
  if (s_verbosity > 2)
    {
      pout () << "Projection inputs: " << endl;

      pout() << "  doSyncProjection = "  << m_doSyncProjection << endl;

      pout() << "  applySyncCorrection = "  << m_applySyncCorrection << endl;

      pout() << "  applyFreestreamCorrection = "  << s_applyVDCorrection
             << endl;

      pout() << "  doQuadInterp = "  << m_doQuadInterp << endl;

      pout() << "  eta = "  <<  m_etaLambda << endl;

      pout() << "  constantLambdaScaling = "  << s_constantLambdaScaling
             << endl;
    }

  // set flag so we don't have to call this function again.
  pp_init = true;

}

void CCProjector::setCrseProj(CCProjector* a_crseProjPtr, int a_nRefCrse)
{
  if (m_crseProjPtr != NULL)
    {
      // if a coarse level is defined, should never undefine it
      CH_assert (a_crseProjPtr != NULL);
      // nRefCrse should never change either.
      CH_assert (m_nRefCrse == a_nRefCrse);
      // check to see if coarse-grid info will need to be redone
      if (!(m_crseProjPtr->getBoxes()==a_crseProjPtr->getBoxes()))
        {
          // new coarse boxes -- need to redefine levelSolver
          const DisjointBoxLayout& crseGrids = a_crseProjPtr->getBoxes();
          defineSolverMGlevel(getBoxes(), &crseGrids);

        } // end if old and new grids are different
      else
        {
          // if old and new coarse grids are the same, no need to
          // do this
        }
    } // end if no coarser level
  else
    {
      // if old crse level didn't exist and new one does, then need
      // to initialize levelsolver
      if (a_crseProjPtr != NULL)
        {
          m_nRefCrse = a_nRefCrse;
          const DisjointBoxLayout& crseGrids = a_crseProjPtr->getBoxes();
          defineSolverMGlevel(getBoxes(), &crseGrids);
        } // end if there's a crse level

    } // end if new level

  // don't forget to actually _set_ the pointer!
  m_crseProjPtr = a_crseProjPtr;
}

// -----------------------------------------------------------
void CCProjector::setFineProj(CCProjector* a_fineProjPtr)
{
  // this is much simpler than resetting crseProj, since
  // there is no need to rebuild levelsolvers or anything like
  // that
  m_fineProjPtr = a_fineProjPtr;
}

void CCProjector::limitSolverCoarsening(bool a_limitSolverCoarsening)
{
  m_limitSolverCoarsening = a_limitSolverCoarsening;
}

#ifdef CH_USE_HDF5
// --------------------------------------------------------------
void CCProjector::writeCheckpointHeader(HDF5Handle& a_handle) const
{

  // nothing needs to be done here -- just
  // included for completeness...

}

// --------------------------------------------------------------
void CCProjector::writeCheckpointLevel(HDF5Handle& a_handle) const
{

  char level_str[20];
  sprintf (level_str, "%d", m_level);
  const std::string label = std::string("level_") + level_str
    + std::string("/projection");

  a_handle.setGroup(label);

  HDF5HeaderData header;

  header.m_real["dx"] = m_dx;
  header.m_box["domain"] = m_domain.domainBox();
  header.m_int["n_ref_crse"] = m_nRefCrse;
  if (m_finest_level)
  {
    header.m_int["finest_level"] = 1;
  }
  else
  {
    header.m_int["finest_level"] = 0;
  }

  // as in AMREuler, don't write static data (comes from ParmParse)

  header.writeToFile(a_handle);

  // now dump out this level's data
  write (a_handle, m_Pi.boxLayout());
  // also want to dump out ghost-cell vales
  write (a_handle, m_Pi, "pi", m_Pi.ghostVect());
  write (a_handle, m_phi, "phi", m_phi.ghostVect());
  write (a_handle, m_eSync, "e_sync", m_eSync.ghostVect());
  write (a_handle, m_eLambda, "e_lambda", m_eLambda.ghostVect());
  write (a_handle, m_grad_eLambda, "grad_e_lambda", m_grad_eLambda.ghostVect());
}

// --------------------------------------------------------------
void CCProjector::readCheckpointHeader(HDF5Handle& a_handle)
{
  //  this doesn't actually do anything...
}

// --------------------------------------------------------------
void CCProjector::readCheckpointLevel(HDF5Handle& a_handle)
{
  char level_str[20];
  sprintf (level_str, "%d", m_level);
  const std::string label = std::string("level_") + level_str
    + std::string("/projection");

  a_handle.setGroup(label);

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout () << "hdf5 header data: " << endl;
      pout () << header << endl;
    }

  // read dx
  if (header.m_real.find("dx") == header.m_real.end())
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain dx");
    }
  if (m_dx != header.m_real["dx"])
    {
      MayDay::Error("CCProjector::readCheckpointLevel: dx in checkfile does not match that in solver");
    }

  // read domain
  if (header.m_box.find("domain") == header.m_box.end())
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain domain");
      if (m_domain.domainBox() != header.m_box["domain"])
        {
          MayDay::Error("CCProjector::readCheckpointLevel: domain in checkfile does not match that in solver!");
        }
    }

  // read nRefCrse
  if (header.m_int.find("n_ref_crse") == header.m_int.end())
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain n_ref_crse");
      if (m_nRefCrse != header.m_int["n_ref_crse"])
      {
        MayDay::Error("CCProjector::readCheckpointLevel: n_ref_crse in checkfile does not match that in solver!");
      }
    }

  // read nRefCrse
  if (header.m_int.find("finest_level") == header.m_int.end())
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain finest_level");
    }
  m_finest_level = (header.m_int["finest_level"] == 1);

  // read grids
  Vector<Box> grids;
  const int grid_status = read(a_handle, grids);
  if (grid_status != 0)
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain a Vector<Box>");
    }
  // can't really check for equality, but at least check for same size...
  if (grids.size() != m_Pi.getBoxes().size())
    {
      MayDay::Error("CCProjector::readCheckpointLevel: grids in checkfile different size than in solver!");
    }

  const DisjointBoxLayout& levelGrids = m_Pi.getBoxes();
  const int piData_status = read<FArrayBox> (a_handle, m_Pi, "pi", levelGrids);

  if (piData_status != 0)
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain pi data");
    }

  const int phiData_status = read<FArrayBox> (a_handle, m_phi, "phi", levelGrids);

  if (phiData_status != 0)
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain phi data");
    }

  const int eSyncData_status = read<FArrayBox> (a_handle, m_eSync, "e_sync", levelGrids);

  if (eSyncData_status != 0)
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain e_sync data");
    }

  const int eLambdaData_status = read<FArrayBox> (a_handle, m_eLambda, "e_lambda", levelGrids);

  if (eLambdaData_status != 0)
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain eLambda data");
    }

  const int grad_eLambdaData_status = read<FluxBox> (a_handle, m_grad_eLambda, "grad_e_lambda", levelGrids);

  if (grad_eLambdaData_status != 0)
    {
      MayDay::Error("CCProjector::readCheckpointLevel: file does not contain grad_e_lambda data");
    }

  // now need to set boundary conditions on all these things
  postRestart();
}

#endif

// --------------------------------------------------------------
// note that this function assumes that all BC's on phi have been set
void CCProjector::gradPhi(LevelData<FArrayBox>& a_gradPhi, int a_dir) const
{
  // first compute edge-centering of gradPhi
  DataIterator dit = a_gradPhi.dataIterator();
  dit.reset();

  int edgeDir = -1;
  const Box& edgeBox = a_gradPhi[dit].box();
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (edgeBox.type(dir) == IndexType::NODE)
        {
          CH_assert (edgeDir == -1);
          edgeDir = dir;
        }
    }
  CH_assert(edgeDir != -1);

  const DisjointBoxLayout& grids = getBoxes();

  for (dit.reset(); dit.ok(); ++dit)
    {
      // call gradient subroutine directly, since this is one-directional

      // I don't think this function is called by anybody ever
      // I think this is correct, but just in case, flag it
      // (dfmartin@lbl.gov)
      CH_assert (false);

      Gradient::singleBoxMacGrad(a_gradPhi[dit],
                                 m_phi[dit],
                                 0, 0, 1, grids[dit], m_dx,
                                 a_dir, edgeDir, m_gradIVS[dit]);
    } // end loop over grids
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on phi have been set
void CCProjector::gradPhi(LevelData<FluxBox>& a_gradPhi) const
{
  Gradient::levelGradientMAC(a_gradPhi, m_phi, m_dx);
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on Pi have been set
void CCProjector::gradPi(LevelData<FArrayBox>& a_gradPi, int a_dir) const
{
  DataIterator dit = a_gradPi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
    {
      // call FORTRAN subroutine directly...
      FORT_GRADCC(CHF_FRA1(a_gradPi[dit],0),
                  CHF_CONST_FRA1(m_Pi[dit],0),
                  CHF_BOX(getBoxes()[dit]),
                  CHF_CONST_REAL(m_dx),
                  CHF_INT(a_dir));
    }
}

// --------------------------------------------------------------
// note that this function assumes that all BC's on Pi have been set
void CCProjector::gradPi(LevelData<FArrayBox>& a_gradPi) const
{
  Gradient::levelGradientCC(a_gradPi, m_Pi, m_dx);
}

// -------------------------------------------------------------
// note that this function assumes all BC's on eSync have been set
// also only returns valid grad(eSync) in valid regions of grids
void CCProjector::grad_eSync(LevelData<FArrayBox>& a_grad_eSync) const
{
  if (m_fineProjPtr != NULL)
    {
      const LevelData<FArrayBox>& fine_eSync = m_fineProjPtr->eSync();
      int nRefFine = m_fineProjPtr->nRefCrse();

      Gradient::compGradientCC(a_grad_eSync, m_eSync,  &fine_eSync,
                               m_dx, nRefFine, m_domain);
    }
  else
    {
      // no finer level -- can do level-operator gradient
      Gradient::levelGradientCC(a_grad_eSync, m_eSync, m_dx);
    }
}

// --------------------------------------------------------------
// note that this function assumes all BC's on eSync have been set
// also only returns valid grad(eSync) in valid regions of grids
void CCProjector::grad_eSync(LevelData<FArrayBox>& a_grad_eSync, int a_dir) const
{
  // this is not yet implemented (not sure if i'll really need it!)
  MayDay::Warning("directional grad_eSync not yet implemented!");
}

// --------------------------------------------------------------
void CCProjector::levelMacProject(LevelData<FluxBox>& a_uEdge,
                                  Real a_oldTime, Real a_dt)
{
  // MAC rhs should have no ghost values
  LevelData<FArrayBox> MacRHS(getBoxes(),1);

  Divergence::levelDivergenceMAC(MacRHS, a_uEdge, m_dx);

  LevelData<FArrayBox>* crseBCPtr = NULL;
  Real CFscale = 0.5*a_dt;

  if (m_crseProjPtr != NULL)
    {
      // coarse-fine BC is 0.5*dt*(coarse Pi)
      const DisjointBoxLayout crseGrids = m_crseProjPtr->getBoxes();
      crseBCPtr = new LevelData<FArrayBox>(crseGrids,1);
      LevelData<FArrayBox>& crseBC = *crseBCPtr;

      Interval comps(0,0);
      const LevelData<FArrayBox>& crsePi = m_crseProjPtr->Pi();
      crsePi.copyTo(comps,crseBC,comps);

      DataIterator dit = crseBC.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
        {
          crseBC[dit] *= CFscale;
        }
    }

  // report sum(rhs)
  if (s_verbosity >= 2)
    {
      DisjointBoxLayout* finerGridsPtr = NULL;
      int nRefFine = -1;
      Real sumRHS = computeSum(MacRHS, finerGridsPtr,
                               nRefFine, m_dx, MacRHS.interval());
      pout() << "MAC projection -- sum(RHS) = " << sumRHS << endl;
    }

  // now solve for phi
  solveMGlevel(m_phi, crseBCPtr, MacRHS);

  // now correct

  // first re-apply physical and copy BC's
  Interval phiComps(0,0);
  m_phi.exchange(phiComps);
  BCHolder bcHolder = m_physBCPtr->gradMacPressureFuncBC();
  const DisjointBoxLayout& levelGrids = getBoxes();

  DataIterator dit = m_phi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      bcHolder.operator()(m_phi[dit],
                          levelGrids[dit],
                          m_domain,
                          m_dx,
                          false); // not homogeneous
    }

#define CORNERCOPIER
#ifdef CORNERCOPIER
  CornerCopier cornerExchangeCopier(m_phi.getBoxes(),
                                    m_phi.getBoxes(),
                                    m_domain,
                                    m_phi.ghostVect(),
                                    true);

  m_phi.exchange(m_phi.interval(), cornerExchangeCopier);
#endif

  applyMacCorrection(a_uEdge, CFscale);

  if (crseBCPtr != NULL)
    {
      delete crseBCPtr;
      crseBCPtr = NULL;
    }
  // that should be it!
}

// --------------------------------------------------------------
void CCProjector::applyMacCorrection(LevelData<FluxBox>& a_uEdge,
                                     Real CFscale)
{
  // this function computes a_uEdge - grad(phi)
  // (assumes all BC's already set)

  const DisjointBoxLayout& levelGrids = a_uEdge.getBoxes();
  LevelData<FluxBox> gradPhi(levelGrids,1);

  LevelData<FArrayBox>* crseBCDataPtr=NULL;
  // should C/F BC's already be set? assume not and reset them here
  if (m_crseProjPtr != NULL)
    {
      const DisjointBoxLayout crseGrids = m_crseProjPtr->getBoxes();
      crseBCDataPtr = new LevelData<FArrayBox>(crseGrids,1);
      LevelData<FArrayBox>& crseBCData = *crseBCDataPtr;

      Interval comps(0,0);
      m_crseProjPtr->Pi().copyTo(comps,crseBCData,comps);

      // now scale the coarse-level BC
      DataIterator ditCrse = crseBCData.dataIterator();

      for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
        {
          crseBCData[ditCrse] *= CFscale;
        }
    }

  // compute gradient
  Gradient::levelGradientMAC(gradPhi,m_phi, crseBCDataPtr, m_dx,
                             m_gradIVS,
                             m_cfInterp);

  // now rescale and subtract from uEdge
  DataIterator dit = a_uEdge.dataIterator();

  // assumes that uEdge and phi can use same dataIterator
  for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& thisGrad = gradPhi[dit];
      FluxBox& thisEdgeVel = a_uEdge[dit];

      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisGradDir = thisGrad[dir];
          FArrayBox& thisVelDir = thisEdgeVel[dir];

          thisVelDir -= thisGradDir;
        }
    }

  if (crseBCDataPtr!=NULL)
    {
      delete crseBCDataPtr;
      crseBCDataPtr = NULL;
    }
}

// ---------------------------------------------------------------
void CCProjector::LevelProject(LevelData<FArrayBox>& a_velocity,
                               LevelData<FArrayBox>* a_crseVelPtr,
                               const Real a_newTime, const Real a_dt)
{
  LevelData<FArrayBox>* crseDataPtr=NULL;
  LevelData<FArrayBox> levelProjRHS(getBoxes(),1);

  // just to be safe.  proper place for this may be outside this function
  Interval velComps(0,SpaceDim-1);
  a_velocity.exchange(velComps);

  if (m_level > 0)
    {
      if (doQuadInterp())
        {
          CH_assert(a_crseVelPtr != NULL);
          CH_assert(m_crseProjPtr != NULL);

          // will need to add dt*gradPi from crse vel.
          const DisjointBoxLayout& crseLevelGrids = m_crseProjPtr->getBoxes();
          crseDataPtr= new LevelData<FArrayBox>(crseLevelGrids, SpaceDim);
          LevelData<FArrayBox>& crseData = *crseDataPtr;
          m_crseProjPtr->gradPi(crseData);

          // now multiply by dt and add crseVel
          DataIterator ditCrse = crseData.dataIterator();
          LevelData<FArrayBox>& crseVel = *a_crseVelPtr;

          // inherent assumption that crseVel and crseGradPi share
          // the same layout here
          for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
            {
              FArrayBox& thisCrseVel = crseVel[ditCrse];
              FArrayBox& thisCrseGradPi = crseData[ditCrse];

              // note that we multiply by dt on _this_ level
              thisCrseGradPi *= a_dt;

              // now add crseVel
              thisCrseGradPi.plus(thisCrseVel);
            }
        }
      else
        {
          // alternative is not implemented at this point
          MayDay::Error("non-quadratic interpolation BC's not implemented");
        }
    } // end if coarser level exists

  QuadCFInterp velCFInterp;
  if (doQuadInterp() && m_crseProjPtr != NULL)
    {
      // define two-component CF-interp object for velocities
      velCFInterp.define(a_velocity.getBoxes(),
                         &(crseDataPtr->getBoxes()),
                         m_dx, m_nRefCrse, SpaceDim, m_domain);
    }

  // now compute RHS for projection
  Divergence::levelDivergenceCC(levelProjRHS, a_velocity, crseDataPtr,
                                m_dx, doQuadInterp(), velCFInterp);

  // for proper scaling of pi, divide this by dt
  DataIterator dit = levelProjRHS.dataIterator();
  Real dtScale = 1.0/a_dt;
  for (dit.reset(); dit.ok(); ++dit)
    {
      levelProjRHS[dit] *= dtScale;
    }

  // set up coarse BC's for solve, then solve
  const LevelData<FArrayBox>* crsePiPtr = NULL;
  if (m_crseProjPtr != NULL) crsePiPtr = &(m_crseProjPtr->Pi());

  // report sum(rhs)
  if (s_verbosity >= 2)
    {
      DisjointBoxLayout* finerGridsPtr = NULL;
      int nRefFine = -1;
      Real sumRHS = computeSum(levelProjRHS, finerGridsPtr,
                               nRefFine, m_dx, levelProjRHS.interval());
      pout() << "Level projection -- sum(RHS) = " << sumRHS << endl;
    }

  // now solve for pi

  solveMGlevel(m_Pi, crsePiPtr, levelProjRHS);

  // apply appropriate physical BC's
  BCHolder bcHolder = m_physBCPtr->gradPiFuncBC();
  const DisjointBoxLayout& levelGrids = getBoxes();

  // loop over grids...
  for (dit.reset(); dit.ok(); ++dit)
    {
      bcHolder.operator()(m_Pi[dit],
                          levelGrids[dit],
                          m_domain,
                          m_dx,
                          false); // not homogeneous
    }

  if (m_crseProjPtr != NULL)
    {
      // reapply coarse-fine BC's here if necessary
      m_cfInterp.coarseFineInterp(m_Pi, *crsePiPtr);
      // also clean up after ourselves!
      delete crseDataPtr;
      crseDataPtr = NULL;

    }

  dtScale = a_dt;
  correctCCVelocities(a_velocity, dtScale);

  // check resulting velocity field
  // to do this, need to reset physical BC's
}

// ---------------------------------------------------------------
void CCProjector::correctCCVelocities(LevelData<FArrayBox>& a_velocity,
                                      const Real scale) const
{
  const DisjointBoxLayout& grids = getBoxes();
  LevelData<FArrayBox> gradPi(grids, SpaceDim);

  // assume that all relevant BC's have already been set
  Gradient::levelGradientCC(gradPi, m_Pi, m_dx);

  // vel = vel - scale*gradPi
  DataIterator dit = gradPi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& thisGradPi = gradPi[dit];
      FArrayBox& thisVel = a_velocity[dit];

      thisGradPi *= scale;
      thisVel -= thisGradPi;
    }
}

// ---------------------------------------------------------------
void CCProjector::doSyncOperations(Vector<LevelData<FArrayBox>* >& a_velocity,
                                   Vector<LevelData<FArrayBox>* >& a_lambda,
                                   const Real a_newTime, const Real a_dtLevel)
{
  AMRMultiGrid<LevelData<FArrayBox> >* dspSolver = new
    AMRMultiGrid<LevelData<FArrayBox> >;

  defineMultiGrid(*dspSolver, a_velocity, false);

  // now call sync projection
  doSyncProjection(a_velocity, a_newTime, a_dtLevel, *dspSolver);

  // Need to delete dspSolver->m_bottomSolver from AMRMultiGrid.
  delete dspSolver;

  AMRMultiGrid<LevelData<FArrayBox> >* cvdcSolver = new
    AMRMultiGrid<LevelData<FArrayBox> >;
  defineMultiGrid(*cvdcSolver, a_velocity, false);

  // now do freestream preservation solve
  computeVDCorrection(a_lambda, a_newTime, a_dtLevel, *cvdcSolver);

  // Need to delete cvdcSolver->m_bottomSolver from AMRMultiGrid.
  delete cvdcSolver;
}

// ---------------------------------------------------------------
void CCProjector::doPostRegridOps(Vector<LevelData<FArrayBox>* >& a_velocity,
                                  Vector<LevelData<FArrayBox>* >& a_lambda,
                                  const Real a_dt, const Real a_time)
{
  AMRMultiGrid<LevelData<FArrayBox> >* bigSolverPtr = new
    AMRMultiGrid<LevelData<FArrayBox> >;

  defineMultiGrid(*bigSolverPtr, a_lambda, false);

  // for inviscid flow, only do this
  // now do freestream preservation solve
  computeVDCorrection(a_lambda, a_time, a_dt, *bigSolverPtr);

  // Need to delete bigSolverPtr->m_bottomSolver from AMRMultiGrid.
  delete bigSolverPtr;
}

// ---------------------------------------------------------------
void CCProjector::doSyncProjection(Vector<LevelData<FArrayBox>* >& a_velocity,
                                   const Real a_newTime, const Real a_dtSync,
                                   AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                                   )
{
  if (doSyncProjection())
  {
    // set up temp storage
    int vectorSize = a_velocity.size();

    Vector<LevelData<FArrayBox>* > syncRHS(vectorSize,NULL);
    Vector<LevelData<FArrayBox>* > syncCorr(vectorSize,NULL);

    // initially
    CCProjector* levelProjPtr = this;
    while (!levelProjPtr->isFinestLevel())
      {
        levelProjPtr = levelProjPtr->fineProjPtr();
      }

    int finestLevel = levelProjPtr->getLevel();

    // reset projection ptr to this level
    levelProjPtr = this;

    // this is a bogus leveldata pointer for the composite divergences
    LevelData<FArrayBox>* crseVelPtr=NULL;
    LevelData<FArrayBox>* fineVelPtr = NULL;

    // loop over levels to allocate storage and compute RHS
    for (int lev = m_level; lev<= finestLevel; lev++)
      {
        const DisjointBoxLayout& levelGrids = a_velocity[lev]->getBoxes();

        syncRHS[lev]= new LevelData<FArrayBox>(levelGrids,1);
        syncCorr[lev] = &(levelProjPtr->eSync());

        Real dx = levelProjPtr->dx();
        const ProblemDomain& levelDomain = levelProjPtr->dProblem();
        int nRefCrse = levelProjPtr->nRefCrse();
        int nRefFine = -1;

        if (lev > 0) crseVelPtr = a_velocity[lev-1];
        if (lev < finestLevel)
          {
            fineVelPtr = a_velocity[lev+1];
            nRefFine = levelProjPtr->fineProjPtr()->nRefCrse();
          }
        else
          {
            fineVelPtr = NULL;
          }

        // do exchange here
        a_velocity[lev]->exchange(a_velocity[lev]->interval());

        // now compute divergence
        Divergence::compDivergenceCC(*syncRHS[lev],*a_velocity[lev],
                                     crseVelPtr, fineVelPtr, dx,
                                     nRefCrse, nRefFine, levelDomain,
                                     doQuadInterp());

        // increment projection level
        levelProjPtr = levelProjPtr->fineProjPtr();

      } // end loop over levels for computation of RHS

    // take care of C/F Boundary condition for lBase
    if (m_level > 0)
      {
        syncCorr[m_level-1] = &(m_crseProjPtr->eSync());
        // need to rescale crse-level BC
        LevelData<FArrayBox>& crseSyncCorr = *(syncCorr[m_level-1]);
        DataIterator ditCrse = crseSyncCorr.dataIterator();
        for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
          {
            crseSyncCorr[ditCrse] *= a_dtSync;
          }
      }

    // debugging tool -- want sum(rhs) = 0
    Vector<int> nRefFineVect(syncRHS.size());
    levelProjPtr = this;
    for (int lev=m_level; lev<=finestLevel; lev++)
      {
        if (lev < finestLevel)
          {
            nRefFineVect[lev] = levelProjPtr->fineProjPtr()->nRefCrse();
          }
        else
          {
            nRefFineVect[lev] = 1;
          }
        levelProjPtr = levelProjPtr->fineProjPtr();
      }

    Interval sumComps(0,0);
    Real sumRHS;
    sumRHS =  computeSum(syncRHS, nRefFineVect, m_dx, sumComps, m_level);
    pout() << "Sum(RHS) for sync solve = "
           << setiosflags(ios::scientific) << sumRHS << endl;

    // now solve
    a_solver.solve(syncCorr, syncRHS, vectorSize-1, m_level,
                   true); // initialize syncCorr to zero

    // apply physical boundary conditions before taking gradients
    levelProjPtr = this;

    BCHolder bcHolder = m_physBCPtr->gradESyncFuncBC();

    for (int lev=m_level; lev<=finestLevel; lev++)
      {

        // loop over grids
        LevelData<FArrayBox>& level_eSync = *syncCorr[lev];
        const ProblemDomain& levelDomain = levelProjPtr->dProblem();
        Real levelDx = levelProjPtr->dx();
        const DisjointBoxLayout& levelGrids = levelProjPtr->getBoxes();
        DataIterator ditLev = level_eSync.dataIterator();

        for (ditLev.reset(); ditLev.ok(); ++ditLev)
          {
            bcHolder.operator()(level_eSync[ditLev],
                                levelGrids[ditLev],
                                levelDomain,
                                levelDx,
                                false); // not homogeneous
          }
        // exchange here should be unnecessary
        //level_eSync.exchange(corrComps);
        levelProjPtr = levelProjPtr->fineProjPtr();
      }

    // now apply sync correction
    Real scale = 1.0;
    LevelData<FArrayBox>* crseCorrPtr = NULL;
    if (m_level > 0) crseCorrPtr = syncCorr[m_level-1];

    applySyncCorrection(a_velocity, scale, crseCorrPtr);

    // cleaning up after we're done --
    // now need to rescale eSync to look like other pressures
    // start with level-1 because we rescaled it for the BC
    // also need to delete syncRHS

    Real corrScale = 1.0/a_dtSync;
    for (int lev=m_level-1; lev<=finestLevel; lev++)
      {
        if (lev >=0)
          {
            LevelData<FArrayBox>& levelCorr = *syncCorr[lev];
            DataIterator dit = levelCorr.dataIterator();
            for (dit.reset(); dit.ok(); ++dit)
              {
                levelCorr[dit] *= corrScale;
              }
          } // end if level > 0
      } // end loop over levels for rescaling eSync

    // loop over levels to delete syncRHS
    for (int lev= m_level; lev<=finestLevel; lev++)
      {
        delete syncRHS[lev];
        syncRHS[lev] = NULL;
      }
  } // end if do sync in the first place
}

// ---------------------------------------------------------------
void CCProjector::applySyncCorrection(Vector<LevelData<FArrayBox>* >& a_velocities,
                                      const Real a_scale,
                                      LevelData<FArrayBox>* crseCorrPtr)
{
  if (m_applySyncCorrection)
    {
      const DisjointBoxLayout& levelGrids = getBoxes();

      LevelData<FArrayBox> levelGradCorr(levelGrids, SpaceDim);
      LevelData<FArrayBox>* fineCorrPtr = NULL;
      int nRefFine = -1;

      if (!isFinestLevel())
        {
          fineCorrPtr = &(m_fineProjPtr->eSync());
          nRefFine = m_fineProjPtr->nRefCrse();
        }

      Gradient::compGradientCC(levelGradCorr, m_eSync,
                               crseCorrPtr, fineCorrPtr,
                               m_dx, nRefFine, m_domain,
                               m_gradIVS,
                               m_cfInterp);

      // now loop over grids and correct
      LevelData<FArrayBox>& levelVel = *a_velocities[m_level];

      DataIterator dit = levelVel.dataIterator();
       for (dit.reset(); dit.ok(); ++dit)
        {
          levelGradCorr[dit] *= a_scale;
          levelVel[dit] -= levelGradCorr[dit];
        }

      // now recursively call finer level
      if (!isFinestLevel())
        {
          m_fineProjPtr->applySyncCorrection(a_velocities, a_scale,
                                             &m_eSync);

          // average down finer level to this level
          const DisjointBoxLayout& fineGrids = m_fineProjPtr->getBoxes();
          int nRefFine = m_fineProjPtr->nRefCrse();
          CoarseAverage avgDownThingy(fineGrids, SpaceDim, nRefFine);

          avgDownThingy.averageToCoarse(*a_velocities[m_level],
                                        *a_velocities[m_level+1]);
        }
    }
}

// ---------------------------------------------------------------
void CCProjector::computeVDCorrection(Vector<LevelData<FArrayBox>* >& a_lambda,
                                      const Real a_newTime, const Real a_dtSync,
                                      AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                                   )
{
  // only do this if eta > 0
  if (m_etaLambda > 0)
    {

      if (m_level == 0) s_lambda_timestep = a_dtSync;
      // also need to handle initial case where level 0
      // hasn't been hit yet. in this case, back out
      // what the level 0 timestep is.
      if (s_lambda_timestep ==0.0)
        {
          int timestepMult = 1;
          CCProjector* levelProjPtr = this;
          for (int lev=m_level; lev>0; lev--)
            {
              timestepMult *= levelProjPtr->nRefCrse();
              levelProjPtr = levelProjPtr->crseProjPtr();
            }
          s_lambda_timestep = a_dtSync*timestepMult;
        }

      // set up temp storage
      int vectorSize = a_lambda.size();

      Vector<LevelData<FArrayBox>* > VDCorr(vectorSize);

      CCProjector* levelProjPtr = this;
      while (!levelProjPtr->isFinestLevel())
        {
          levelProjPtr = levelProjPtr->fineProjPtr();
        }

      int finestLevel = levelProjPtr->getLevel();

      // reset projection ptr to this level
      levelProjPtr = this;

      CH_assert (a_dtSync != 0);
      Real lambdaMult;
      if (s_constantLambdaScaling)
        {
          lambdaMult = m_etaLambda/s_lambda_timestep;
        }
      else
        {
          lambdaMult = m_etaLambda/a_dtSync;
        }

      for (int lev = m_level; lev<=finestLevel; lev++)
        {
          // stuff eLambda into VDCorr array
          VDCorr[lev] = &(levelProjPtr->eLambda());

          // use lambda array as RHS for solve
          LevelData<FArrayBox>& thisLambda = *(a_lambda[lev]);

          DataIterator dit = thisLambda.dataIterator();

          // RHS = (lambda-1)*eta/dtSync
          for (dit.reset(); dit.ok(); ++dit)
            {
              thisLambda[dit] -= 1.0;
              thisLambda[dit] *= lambdaMult;
              levelProjPtr->eLambda()[dit].setVal(0.0);
            }

          // advance projection pointer to next finer level
          levelProjPtr = levelProjPtr->fineProjPtr();
        } // end loop over levels

      // if necessary, define coarse BC here
      if (m_level > 0)
        {
          VDCorr[m_level-1] = &(m_crseProjPtr->eLambda());
          // rescale crse-level BC
          // try setting coarse-level BC to 0
          setValLevel(m_crseProjPtr->eLambda(), 0.0);
        }

      // debugging tool -- want sum(rhs) = 0
      Vector<int> nRefFineVect(a_lambda.size());
      levelProjPtr = this;
      for (int lev=m_level; lev<=finestLevel; lev++)
        {
          if (lev < finestLevel)
            {
              nRefFineVect[lev] = levelProjPtr->fineProjPtr()->nRefCrse();
            }
          else
            {
              nRefFineVect[lev] = 1;
            }
          levelProjPtr = levelProjPtr->fineProjPtr();
        }

      Interval sumComps(0,0);
      Real sumRHS;
      sumRHS =  computeSum(a_lambda, nRefFineVect, m_dx, sumComps, m_level);
      pout() << "Time = " << a_newTime << " Sum(RHS) for VD solve = "
             << setiosflags(ios::scientific) << sumRHS << endl;

      // now solve
      a_solver.solve(VDCorr, a_lambda, vectorSize-1, m_level,
                     true); // initialize VDCorr to zero

      // apply all bc's here

      BCHolder bcHolder = m_physBCPtr->gradELambdaFuncBC();

      Interval corrComps(0,0);
      levelProjPtr = this;
      for (int lev=m_level; lev<=finestLevel; lev++)
      {
        LevelData<FArrayBox>& levelCorr = *(VDCorr[lev]);
        const ProblemDomain& levelDomain = levelProjPtr->dProblem();
        Real levelDx = levelProjPtr->dx();
        const DisjointBoxLayout& levelGrids = levelProjPtr->getBoxes();
        DataIterator ditLev = levelCorr.dataIterator();
        for (ditLev.reset(); ditLev.ok(); ++ditLev)
          {
            bcHolder.operator()(levelCorr[ditLev],
                                levelGrids[ditLev],
                                levelDomain,
                                levelDx,
                                false); // not homogeneous
          }
        levelCorr.exchange(corrComps);
        levelProjPtr = levelProjPtr->fineProjPtr();
      }

      // compute correction field
      computeGrad_eLambda();

      lambdaMult = 1.0/lambdaMult;
      // now return lambda to previous state
      for (int lev = m_level; lev<=finestLevel; lev++)
        {
          LevelData<FArrayBox>& levelLambda = *(a_lambda[lev]);
          DataIterator dit = levelLambda.dataIterator();

          for (dit.reset(); dit.ok(); ++dit)
            {
              levelLambda[dit] *= lambdaMult;
              levelLambda[dit] += 1.0;
            }
        }
    } // end if eta > 0
}

// ---------------------------------------------------------------
void CCProjector::computeGrad_eLambda()
{
  // do this recursively starting at the finest level

  if (!isFinestLevel())
    {
      m_fineProjPtr->computeGrad_eLambda();
    }

  // now compute composite gradient on this level
  LevelData<FArrayBox>* crseBCPtr = NULL;

  if (m_crseProjPtr != NULL) crseBCPtr = &(m_crseProjPtr->eLambda());

  // recall that composite MAC gradient is the same as the level
  // gradient, since finer level is not considered to be part
  // of this level (take care of covered regions with avgDown)

  Gradient::levelGradientMAC(m_grad_eLambda, m_eLambda,
                             crseBCPtr, m_dx,
                             m_gradIVS,
                             m_cfInterp);

  if (!s_applyVDCorrection)
    {
      setValLevel(m_grad_eLambda, 0.0);
    }

  // now average down finer level grad_eLambda if necessary
  if (!isFinestLevel())
    {
      const DisjointBoxLayout& fineGrids = m_fineProjPtr->getBoxes();
      int nRefFine = m_fineProjPtr->nRefCrse();
      int nComp = 1;

      CoarseAverageEdge avgDownObject(fineGrids, nComp, nRefFine);

      const LevelData<FluxBox> & fineEdgeGrad = m_fineProjPtr->grad_eLambda();

      avgDownObject.averageToCoarse(m_grad_eLambda, fineEdgeGrad);
    }
}

// ---------------------------------------------------------------
void CCProjector::rescaleGrad_eLambda(Real a_dx_lbase)
{
  // do this recursively starting at the finest level

  if (!isFinestLevel())
    {
      m_fineProjPtr->rescaleGrad_eLambda(a_dx_lbase);
    }

  // now rescale correction on this level.  we rescale correction
  // by (dt(lbase)/dt(this_level) = (dx(lbase)/m_dx for consistency
  Real scale = a_dx_lbase/m_dx;
  DataIterator dit = m_grad_eLambda.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int dir=0; dir<SpaceDim; dir++)
        {
          m_grad_eLambda[dit][dir] *= scale;
        }
    }
}

// ---------------------------------------------------------------
void CCProjector::initialLevelProject(LevelData<FArrayBox>& a_velocity,
                                      LevelData<FArrayBox>* a_crseVelPtr,
                                      const Real a_newTime, const Real a_dt)
{
  LevelData<FArrayBox>* crseDataPtr=NULL;
  LevelData<FArrayBox> levelProjRHS(getBoxes(),1);

  // do coarse-fine BC
  if (m_level > 0)
    {
      if (doQuadInterp())
        {
          CH_assert(a_crseVelPtr != NULL);
          CH_assert(m_crseProjPtr != NULL);

          // will need to add dt*gradPi to crse vel
          const DisjointBoxLayout& crseLevelGrids = m_crseProjPtr->getBoxes();
          crseDataPtr = new LevelData<FArrayBox>(crseLevelGrids, SpaceDim);
          LevelData<FArrayBox>& crseData = *crseDataPtr;
          m_crseProjPtr->gradPi(crseData);

          // now multiply by dt and add crseVel
          DataIterator ditCrse = crseData.dataIterator();
          LevelData<FArrayBox>& crseVel = *a_crseVelPtr;

          // inherent assumption that crseVel & crseGradPi share same layout
          for (ditCrse.reset(); ditCrse.ok(); ++ditCrse)
            {
              FArrayBox& thisCrseVel = crseVel[ditCrse];
              FArrayBox& thisCrseGradPi = crseData[ditCrse];

              // note that we multiply by this level's dt
              thisCrseGradPi *= a_dt;
              // now add crseVel
              thisCrseGradPi.plus(thisCrseVel);
            }

        }
      else
        {
          // alternative is not implemented at this point
          MayDay::Error("non-quadratic interpolation BC's not implemented");
        }
    } // end if coarser level exists

  // just to be safe...
  Interval velComps(0,SpaceDim-1);
  a_velocity.exchange(velComps);

  // now compute RHS for projection
  Divergence::levelDivergenceCC(levelProjRHS, a_velocity, crseDataPtr,
                                m_dx, doQuadInterp(), m_cfInterp);

  // for proper scaling of pi, divide this by dt
  DataIterator dit = levelProjRHS.dataIterator();
  Real dtScale = 1.0/a_dt;

  for (dit.reset(); dit.ok(); ++dit)
    {
      levelProjRHS[dit] *= dtScale;
    }

  // set up coarse BC's for solve, then solve
  const LevelData<FArrayBox>* crsePiPtr = NULL;
  if (m_crseProjPtr != NULL) crsePiPtr = &(m_crseProjPtr->Pi());

  // now solve for Pi
  solveMGlevel(m_Pi, crsePiPtr, levelProjRHS);

  if (m_crseProjPtr != NULL)
    {
      // reapply coarse-fine BC's here if necessary
      m_cfInterp.coarseFineInterp(m_Pi, *crsePiPtr);
    }

  // apply appropriate physical BC's
  BCHolder bcHolder = m_physBCPtr->gradPiFuncBC();
  const DisjointBoxLayout& levelGrids = getBoxes();

  // loop over grids...
  for (dit.reset(); dit.ok(); ++dit)
    {
      bcHolder.operator()(m_Pi[dit],
                          levelGrids[dit],
                          m_domain,
                          m_dx,
                          false); // not homogeneous
    }

  Interval piComps(0,0);
  m_Pi.exchange(piComps);

  dtScale = a_dt;
  correctCCVelocities(a_velocity,dtScale);

  // clean up storage
  if (crseDataPtr!=NULL)
    {
      delete crseDataPtr;
      crseDataPtr=NULL;
    }
}

// ---------------------------------------------------------------
void CCProjector::initialSyncProjection(Vector<LevelData<FArrayBox>* >& a_vel,
                                        const Real a_newTime, const Real a_dtSync,
                                        AMRMultiGrid<LevelData<FArrayBox> >& a_solver
                                        )
{
  // actually, in the current algorithm, i _think_ this is the same
  // as the normal sync projection
  doSyncProjection(a_vel, a_newTime, a_dtSync, a_solver);
}

// ---------------------------------------------------------------
void CCProjector::doInitialSyncOperations(Vector<LevelData<FArrayBox>* >& a_vel,
                                          Vector<LevelData<FArrayBox>* >& a_lambda,
                                          const Real a_newTime, const Real a_dtSync)
{
  // now can define multilevel solver
  AMRMultiGrid<LevelData<FArrayBox> >* ispSolverPtr = new
    AMRMultiGrid<LevelData<FArrayBox> >;

  defineMultiGrid(*ispSolverPtr, a_vel, false);

  // now call sync projection
  initialSyncProjection(a_vel, a_newTime, a_dtSync, *ispSolverPtr);

  // Need to delete ispSolverPtr->m_bottomSolver from AMRMultiGrid.
  delete ispSolverPtr;
  AMRMultiGrid<LevelData<FArrayBox> >* cvdcSolverPtr = new
    AMRMultiGrid<LevelData<FArrayBox> >;
  defineMultiGrid(*cvdcSolverPtr, a_lambda, true);

  // now do freestream preservation solve
  computeVDCorrection(a_lambda, a_newTime, a_dtSync, *cvdcSolverPtr);
}

// ---------------------------------------------------------------
void CCProjector::initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_vel,
                                         bool a_homogeneousCFBC)
{
  // first need to define solver
  AMRMultiGrid<LevelData<FArrayBox> > bigSolver;
  defineMultiGrid(bigSolver, a_vel, false);

  initialVelocityProject(a_vel, bigSolver, a_homogeneousCFBC);
}

// ---------------------------------------------------------------
void CCProjector::initialVelocityProject(Vector<LevelData<FArrayBox>* >& a_vel,
                                         AMRMultiGrid<LevelData<FArrayBox> >& a_bigSolver,
                                         bool a_homogeneousCFBC)
{
  // set up temp storage
  int vectorSize = a_vel.size();

  Vector<LevelData<FArrayBox>* > projRHS(vectorSize,NULL);
  Vector<LevelData<FArrayBox>* > projCorr(vectorSize,NULL);

  CCProjector* levelProjPtr = this;

  int finestLevel = vectorSize - 1;

  // this is a bogus levelData pointer for the composite divergence
  LevelData<FArrayBox>* crseVelPtr = NULL;
  LevelData<FArrayBox>* fineVelPtr = NULL;
  Vector<int> nRefFineVect(vectorSize, -1);

  // loop over levels to allocate storage and compute RHS
  for (int lev=m_level; lev <= finestLevel; lev++)
    {
      const DisjointBoxLayout& levelGrids = a_vel[lev]->getBoxes();

      projRHS[lev] = new LevelData<FArrayBox>(levelGrids,1);
      projCorr[lev] = &(levelProjPtr->eSync());

      Real dx = levelProjPtr->dx();
      const ProblemDomain& levelDomain = levelProjPtr->dProblem();
      int nRefCrse = levelProjPtr->nRefCrse();
      int nRefFine = -1;

      if (lev > 0)
        {
          crseVelPtr = a_vel[lev-1];
        }
      if (lev < finestLevel)
        {
          fineVelPtr = a_vel[lev+1];
          nRefFine = levelProjPtr->fineProjPtr()->nRefCrse();
        }
      else
        {
          fineVelPtr = NULL;
        }
      nRefFineVect[lev] = nRefFine;

      // just in case...
      Interval velComps(0,SpaceDim-1);
      a_vel[lev]->exchange(velComps);

      // now compute divergence
      Divergence::compDivergenceCC(*projRHS[lev],*a_vel[lev],crseVelPtr,
                                   fineVelPtr,dx,nRefCrse,nRefFine,
                                   levelDomain,doQuadInterp());

      // increment levelProjPtr
      levelProjPtr= levelProjPtr->fineProjPtr();

    } // end loop over levels for computation of RHS

  // take care of C/F Boundary condition for lBase

  LevelData<FArrayBox>* crseBCPtr = NULL;
  if (m_level > 0)
    {
      CH_assert(m_crseProjPtr != NULL);
      if (a_homogeneousCFBC)
        {
          // need to define coarse BC and set to 0
          const DisjointBoxLayout& crseGrids = m_crseProjPtr->getBoxes();

          crseBCPtr = new LevelData<FArrayBox>(crseGrids,1);

          setValLevel(*crseBCPtr, 0.0);

          projCorr[m_level-1] = crseBCPtr;
        }
      else
        {
          // not sure if this should ever be called
          MayDay::Warning("inhomogeneous BC's on initialVelProject");
          projCorr[m_level-1] = &m_crseProjPtr->eSync();
        }

    } // end if level > 0

  Vector<DisjointBoxLayout> allGrids(vectorSize);
  for (int lev=m_level; lev<=finestLevel; lev++)
    {
      allGrids[lev] = projCorr[lev]->getBoxes();
    }

  Interval sumComps(0,0);
  Real sumRHS;
  sumRHS =  computeSum(projRHS, nRefFineVect, m_dx, sumComps, m_level);
  pout() << "Sum(RHS) for IVP solve = "
         << setiosflags(ios::scientific) << sumRHS << endl;

  // now solve!
  a_bigSolver.solve(projCorr, projRHS, finestLevel, m_level,
                    true); // initialize projCorr to zero

  // apply physical boundary conditions before
  // taking gradients
  Interval corrComps(0,0);
  levelProjPtr = this;
  BCHolder bcHolder = m_physBCPtr->gradESyncFuncBC();

  for (int lev=m_level; lev<=finestLevel; lev++)
    {

      // loop over grids
      LevelData<FArrayBox>& levelCorr = *projCorr[lev];
      const ProblemDomain& levelDomain = levelProjPtr->dProblem();
      Real levelDx = levelProjPtr->dx();
      const DisjointBoxLayout& levelGrids = allGrids[lev];
      // levelProjPtr->getBoxes();
      DataIterator ditLev = levelCorr.dataIterator();

      for (ditLev.reset(); ditLev.ok(); ++ditLev)
        {
          bcHolder.operator()(levelCorr[ditLev],
                              levelGrids[ditLev],
                              levelDomain,
                              levelDx,
                              false); // not homogeneous
        }

      levelCorr.exchange(corrComps);

      levelProjPtr = levelProjPtr->fineProjPtr();
    }

  // now apply correction
  Real scale = 1.0;

  applySyncCorrection(a_vel, scale, crseBCPtr);

  // delete temp storage -- RHS and crse BC
  for (int lev=m_level; lev<=finestLevel; lev++)
    {
      delete projRHS[lev];
      projRHS[lev] = NULL;
    }

  if (crseBCPtr != NULL)
    {
      delete crseBCPtr;
      crseBCPtr = NULL;
    }
}

// ---------------------------------------------------------------
void CCProjector::defineMultiGrid(AMRMultiGrid<LevelData<FArrayBox> >& a_solver,
                                  const Vector<LevelData<FArrayBox>* >& a_vel,
                                  bool a_freestreamSolve)
{
  // assume that this level is lBase...
  int vectorSize = a_vel.size();

  // determine finest existing level
  CCProjector* levelProj = this;

  int finestLevel = vectorSize - 1;

  // put level grids togther into a vector
  Vector<DisjointBoxLayout> allGrids(vectorSize);
  Vector<ProblemDomain> allDomains(vectorSize);
  Vector<Real> allDx(vectorSize);
  Vector<int> allRefRatios(vectorSize);

  levelProj = this;
  // Now levelProj is at level m_level.

  for (int lev = m_level; lev<=finestLevel; lev++)
    {
      // for now, get grids from velocity arrays
      const DisjointBoxLayout& levelGrids = a_vel[lev]->getBoxes();
      allGrids[lev] = levelGrids;
      // allDomains[lev] = levelProj->dProblem();
      // allDx[lev] = levelProj->dx();

      // since nRef should be refinement ratio to next finer level,
      // advance the levelProj pointer to next finer level here.
      // Now levelProj is at level lev.
      if (!levelProj->isFinestLevel())
        {
          // Make levelProj be at level lev+1, instead of at level lev.
          levelProj = levelProj->fineProjPtr();
          // Now levelProj is at level lev+1.
          // refinement ratio from level lev to level lev+1
          allRefRatios[lev] = levelProj->nRefCrse();
        }
    }
  // Now levelProj is at level finestLevel.

  // if coarser level exists. also need to include that for boundary
  // conditions

  //DisjointBoxLayout* crseGridsPtr = NULL;  //not used?  (ndk)

  ProblemDomain baseDomain;
  Real baseDx;
  if (m_level > 0)
    {
      CH_assert(m_crseProjPtr != NULL);
      const DisjointBoxLayout& crseLevelGrids
        = a_vel[m_level-1]->getBoxes();
      allGrids[m_level-1] = crseLevelGrids;
      // allDomains[m_level-1] = m_crseProjPtr->dProblem();
      // allDx[m_level-1] = m_crseProjPtr->dx();
      // refinement ratio from level m_level-1 to level m_level
      allRefRatios[m_level-1] = nRefCrse();
      if (m_level == 1)
        {
          levelProj = this->crseProjPtr();
          // Now levelProj is at level 0.
          baseDomain = levelProj->dProblem();
          baseDx = levelProj->dx();
        }
    }
  else // m_level == 0
    {
      baseDomain = dProblem();
      baseDx = dx();
    }

  // We have filled in allRefRatios[m_level-1:finestLevel-1].
  // Need to fill in allRefRatios[0:m_level-2].
  // Need to fill in allGrids[0:m_level-2].
  // Also need baseDomain == allDomains[0] and baseDx == allDx[0].

  if (m_level > 1) // m_level >= 2
    {
      levelProj = this->crseProjPtr();
      // Now levelProj is at level m_level-1 > 0.
      // refinement ratio from level m_level-2 to level m_level-1
      allRefRatios[m_level-2] = levelProj->nRefCrse();
      levelProj = levelProj->crseProjPtr();
      // Now levelProj is at level m_level-2 >= 0.
      allGrids[m_level-2] = a_vel[m_level-2]->getBoxes();
      for (int lev = m_level-3; lev >= 0; lev--) // m_level >= 3
        {
          // Now levelProj is at level lev+1.
          // refinement ratio from level lev to level lev+1
          allRefRatios[lev] = levelProj->nRefCrse();
          allGrids[lev] = a_vel[lev]->getBoxes();
          if (lev > 0)
            {
              // Now levelProj is at level lev+1.
              levelProj = levelProj->crseProjPtr();
              // Now levelProj is at level lev.
            }
        }
      // Now levelProj is at level 0.
      baseDomain = levelProj->dProblem();
      baseDx = levelProj->dx();
    }

  AMRPoissonOpFactory localPoissonOpFactory;
  // physBCPtr = m_physBCPtr->FreestreamCorrBC();
  // physBCPtr = m_physBCPtr->SyncProjBC();
  localPoissonOpFactory.define(baseDomain,
                               allGrids,
                               allRefRatios,
                               baseDx,
                               ((a_freestreamSolve) ?
                                m_physBCPtr->FreestreamCorrFuncBC() :
                                m_physBCPtr->SyncProjFuncBC() ));

  // You really need to delete this when you're done with a_solver.
  // But m_bottomSolver is a protected field of AMRMultiGrid.
  RelaxSolver<LevelData<FArrayBox> >* bottomSolverPtr = new
    RelaxSolver<LevelData<FArrayBox> >;
  bottomSolverPtr->m_verbosity = s_verbosity;

  // AMRMultiGrid<LevelData<FArrayBox> >& a_solver
  a_solver.define(baseDomain,
                  localPoissonOpFactory,
                  bottomSolverPtr,
                  finestLevel+1);
  // We also want to use m_limitCoarsening:
  // if true, only do multigrid coarsening down to next coarser
  // AMR level (only coarsen by a_nRefCrse).
  // If false, coarsen as far as possible. Only relevant when lBase > 0.
  a_solver.m_verbosity = s_verbosity;
  a_solver.m_eps = 1e-10;
  if (s_solver_tolerance > 0)
    {
      a_solver.m_eps = s_solver_tolerance;
    }
  a_solver.m_pre = s_num_smooth_down; // smoothings before avging
  a_solver.m_post = s_num_smooth_up; // smoothings after avging
}


// -------------------------------------------------------------
void CCProjector::postRestart()
{
  // set C/F BC's -- this depends on the fact that the coarser level
  // has already been read in!
  if (m_level > 0)
    {
      // need to do Pi
      CH_assert (m_cfInterp.isDefined());
      LevelData<FArrayBox>& crsePi = m_crseProjPtr->Pi();
      m_cfInterp.coarseFineInterp(m_Pi, crsePi);
    }

  Interval piComps = m_Pi.interval();
  m_Pi.exchange(piComps);

  // also need to do physical BC's
  BCHolder bcHolder = m_physBCPtr->gradPiFuncBC();
  const DisjointBoxLayout& levelGrids = getBoxes();

  // loop over grids
  DataIterator dit = m_Pi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      bcHolder.operator()(m_Pi[dit],
                          levelGrids[dit],
                          m_domain,
                          m_dx,
                          false); // not homogeneous
    }
}


// -------------------------------------------------------------
void CCProjector::defineSolverMGlevel(const DisjointBoxLayout& a_grids,
                                      const DisjointBoxLayout* a_crseGridsPtr)
{
  AMRPoissonOpFactory localPoissonOpFactory;
  ProblemDomain baseDomain(m_domain); // on this level
  int numSolverLevels;
  if (a_crseGridsPtr != NULL)
    { // coarser level exists:  define solver on two levels
      numSolverLevels = 2;
      // this returns a null domain for me:
      // baseDomain = m_crseProjPtr->dProblem();
      baseDomain.coarsen(m_nRefCrse);
      // baseDomain = ProblemDomain(boundingBox(*a_crseGridsPtr));
      // LayoutIterator litFirst = a_crseGridsPtr->layoutIterator();
      // const Box& bxFirst = a_crseGridsPtr->operator[](litFirst);
      // Box bxBound = bxFirst;
      // for (LayoutIterator lit = a_crseGridsPtr->layoutIterator(); lit.ok(); ++lit)
      // {
      // const Box& thisBox = a_crseGridsPtr->operator[](lit);
      // bxBound.minBox(thisBox);
      // }
      // baseDomain = ProblemDomain(bxBound);
      // }
      Vector<DisjointBoxLayout> allGrids(2);
      allGrids[0] = *a_crseGridsPtr;
      allGrids[1] = a_grids;
      Vector<int> refRatios(1, m_nRefCrse);
      // this returns zero for me:
      // Real dxCrse = m_crseProjPtr->dx();
      Real dxCrse = m_nRefCrse * m_dx;
      localPoissonOpFactory.define(baseDomain,
                                   allGrids,
                                   refRatios,
                                   dxCrse,
                                   m_physBCPtr->LevelPressureFuncBC());
    }
  else
    { // no coarser level:  define solver on only one level
      numSolverLevels = 1;
      localPoissonOpFactory.define(m_domain,
                                   a_grids,
                                   m_dx,
                                   m_physBCPtr->LevelPressureFuncBC());
    }

  // You really need to delete this when you're done with a_solver.
  // But m_bottomSolver is a protected field of AMRMultiGrid.
  RelaxSolver<LevelData<FArrayBox> >* bottomSolverPtr = new
    RelaxSolver<LevelData<FArrayBox> >;
  bottomSolverPtr->m_verbosity = s_verbosity;

  // AMRMultiGrid<LevelData<FArrayBox> >& a_solver
  m_solverMGlevel.define(baseDomain, // on either this level or coarser level
                         localPoissonOpFactory,
                         bottomSolverPtr,
                         numSolverLevels);
  m_solverMGlevel.m_verbosity = s_verbosity;
  m_solverMGlevel.m_eps = 1e-10;
  if (s_solver_tolerance > 0)
    {
      m_solverMGlevel.m_eps = s_solver_tolerance;
    }
  m_solverMGlevel.m_pre = s_num_smooth_down; // smoothings before avging
  m_solverMGlevel.m_post = s_num_smooth_up; // smoothings after avging
}


// -------------------------------------------------------------
void CCProjector::solveMGlevel(LevelData<FArrayBox>&   a_phi,
                               const LevelData<FArrayBox>*   a_phiCoarsePtr,
                               const LevelData<FArrayBox>&   a_rhs)
{
  // Initialize a_phi to zero.
  setValLevel(a_phi, 0.0);
  Vector< LevelData<FArrayBox>* > phiVect;
  Vector< LevelData<FArrayBox>* > rhsVect;
  LevelData<FArrayBox>& rhsRef =
    const_cast< LevelData<FArrayBox>& >(a_rhs);
  int maxLevel;
  if (a_phiCoarsePtr != NULL)
    {
      maxLevel = 1;
      LevelData<FArrayBox>* phiCoarsePtrRef =
        const_cast< LevelData<FArrayBox>* >(a_phiCoarsePtr);
      phiVect.push_back(phiCoarsePtrRef);
      rhsVect.push_back(NULL); // I don't think this will be used
    }
  else
    {
      maxLevel = 0;
    }
  phiVect.push_back(&a_phi);
  rhsVect.push_back(&rhsRef);
  // l_max = maxLevel, l_base = maxLevel
  m_solverMGlevel.solve(phiVect, rhsVect, maxLevel, maxLevel,
                        false); // don't initialize to zero
}
