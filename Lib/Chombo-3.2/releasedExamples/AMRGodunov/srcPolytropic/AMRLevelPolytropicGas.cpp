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
#include "ParmParse.H"
#include "LayoutIterator.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "computeSum.H"
#include "CH_HDF5.H"
#include "AMRIO.H"

#include "AMRLevelPolytropicGas.H"

#include "PolytropicPhysicsF_F.H"

#include "NamespaceHeader.H"

// Constructor
AMRLevelPolytropicGas::AMRLevelPolytropicGas()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas default constructor" << endl;
    }

  m_gdnvPhysics = NULL;
  m_paramsDefined = false;
}

// Destructor
AMRLevelPolytropicGas::~AMRLevelPolytropicGas()
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas destructor" << endl;
    }

  if (m_gdnvPhysics != NULL)
    {
      delete m_gdnvPhysics;
      m_gdnvPhysics = NULL;
    }

  m_paramsDefined = false;
}

void AMRLevelPolytropicGas::defineParams(const Real&                 a_cfl,
                                         const Real&                 a_domainLength,
                                         const int&                  a_verbosity,
                                         const Real&                 a_refineThresh,
                                         const int&                  a_tagBufferSize,
                                         const Real&                 a_initialDtMultiplier,
                                         const GodunovPhysics* const a_godunovPhysics,
                                         const int&                  a_normalPredOrder,
                                         const bool&                 a_useFourthOrderSlopes,
                                         const bool&                 a_usePrimLimiting,
                                         const bool&                 a_useCharLimiting,
                                         const bool&                 a_useFlattening,
                                         const bool&                 a_useArtificialViscosity,
                                         const Real&                 a_artificialViscosity,
                                         const bool&                 a_useSourceTerm,
                                         const Real&                 a_sourceTermScaling,
                                         const bool&                 a_highOrderLimiter)
{
  // Set the CFL number
  m_cfl = a_cfl;

  // Set the physical dimension of the longest side of the domain
  m_domainLength = a_domainLength;

  verbosity(a_verbosity);

  // Set the refinement threshold
  m_refineThresh = a_refineThresh;

  // Set the tag buffer size
  m_tagBufferSize = a_tagBufferSize;

  initialDtMultiplier(a_initialDtMultiplier);

  if (m_gdnvPhysics != NULL)
    {
      delete m_gdnvPhysics;
      m_gdnvPhysics = NULL;
    }

  m_gdnvPhysics = a_godunovPhysics->new_godunovPhysics();

  m_normalPredOrder = a_normalPredOrder;

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_usePrimLimiting      = a_usePrimLimiting;
  m_useCharLimiting      = a_useCharLimiting;
  m_useFlattening        = a_useFlattening;

  // Artificial viscosity coefficient must be greater than zero
  CH_assert(!a_useArtificialViscosity || (a_artificialViscosity >= 0.0));

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = a_useArtificialViscosity;
  m_artificialViscosity    = a_artificialViscosity;

  // Supply a source term to the computation
  m_useSourceTerm     = a_useSourceTerm;
  m_sourceTermScaling = a_sourceTermScaling;

  // Use a high-order limiter?
  m_highOrderLimiter = a_highOrderLimiter;

  m_paramsDefined = true;
}

// This instance should never get called - historical
void AMRLevelPolytropicGas::define(AMRLevel*  a_coarserLevelPtr,
                                   const Box& a_problemDomain,
                                   int        a_level,
                                   int        a_refRatio)
{
  ProblemDomain physdomain(a_problemDomain);

  MayDay::Error("AMRLevelPolytropicGas::define:\n\tShould never be called with a Box for a problem domain");
}

// Define new AMR level
void AMRLevelPolytropicGas::define(AMRLevel*            a_coarserLevelPtr,
                                   const ProblemDomain& a_problemDomain,
                                   int                  a_level,
                                   int                  a_refRatio)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::define " << a_level << endl;
    }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
    {
      AMRLevelPolytropicGas* amrGodPtr = dynamic_cast<AMRLevelPolytropicGas*>(a_coarserLevelPtr);

      if (amrGodPtr != NULL)
        {
          m_cfl           = amrGodPtr->m_cfl;
          m_domainLength  = amrGodPtr->m_domainLength;
          m_refineThresh  = amrGodPtr->m_refineThresh;
          m_tagBufferSize = amrGodPtr->m_tagBufferSize;
        }
      else
        {
          MayDay::Error("AMRLevelPolytropicGas::define: a_coarserLevelPtr is not castable to AMRLevelPolytropicGas*");
        }
    }

  // Compute the grid spacing
  m_dx = m_domainLength / a_problemDomain.domainBox().size(0);

  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more
  m_numGhost = 1;

  CH_assert(m_gdnvPhysics != NULL);
  CH_assert(isDefined());
  m_gdnvPhysics->define(m_problem_domain,m_dx);

  // Number and names of conserved states
  m_numStates  = m_gdnvPhysics->numConserved();
  m_stateNames = m_gdnvPhysics->stateNames();
}

// Advance by one timestep
Real AMRLevelPolytropicGas::advance()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::advance level " << m_level << " to time " << m_time + m_dt << endl;
    }

  // Copy the new to the old
  for (DataIterator dit = m_UNew.dataIterator(); dit.ok(); ++dit)
    {
      m_UOld[dit()].copy(m_UNew[dit()]);
    }

  Real newDt = 0.0;

  // Set up arguments to LevelGodunov::step based on whether there are
  // coarser and finer levels

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  const LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  const LevelData<FArrayBox>* coarserDataOld = &dummyData;
  const LevelData<FArrayBox>* coarserDataNew = &dummyData;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  // A coarser level exists
  if (m_hasCoarser)
    {
      AMRLevelPolytropicGas* coarserPtr = getCoarserLevel();

      // Recall that my flux register goes between my level and the next
      // finer level
      coarserFR = &coarserPtr->m_fluxRegister;

      coarserDataOld = &coarserPtr->m_UOld;
      coarserDataNew = &coarserPtr->m_UNew;

      tCoarserNew = coarserPtr->m_time;
      tCoarserOld = tCoarserNew - coarserPtr->m_dt;
    }

  // A finer level exists
  if (m_hasFiner)
    {
      // Recall that my flux register goes between my level and the next
      // finer level
      finerFR = &m_fluxRegister;
    }

  // Source term leveldata
  LevelData<FArrayBox> sourceData;

  // Set up source term for hyperbolic update
  if (m_useSourceTerm)
    {
      // Define source term leveldata
      IntVect ivGhost = m_numGhost * IntVect::Unit;
      sourceData.define(m_grids,m_gdnvPhysics->numPrimitives(),ivGhost);

      for (DataIterator dit = sourceData.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& sourceFAB = sourceData[dit()];
          const FArrayBox& consFAB = m_UNew[dit()];

          FORT_SETSOURCEPRIM(CHF_FRA(sourceFAB),
                             CHF_CONST_FRA(consFAB),
                             CHF_CONST_REAL(m_sourceTermScaling),
                             CHF_BOX(sourceFAB.box()));
        }
    }

  // we don't need the flux in the simple hyperbolic case...
  LevelData<FArrayBox> flux[SpaceDim];

  // Advance the solve one timestep
  newDt = m_levelGodunov.step(m_UNew,
                              flux,
                              *finerFR,
                              *coarserFR,
                              sourceData,
                              *coarserDataOld,
                              tCoarserOld,
                              *coarserDataNew,
                              tCoarserNew,
                              m_time,
                              m_dt);

  // Update with source term (2nd order accurate)
  if (m_useSourceTerm)
    {
      for (DataIterator dit = sourceData.dataIterator(); dit.ok(); ++dit)
        {
          const FArrayBox& consOldFAB = m_UOld[dit()];
          FArrayBox&       consNewFAB = m_UNew[dit()];

          FArrayBox sourceFAB(consNewFAB.box(),consNewFAB.nComp());

          FORT_SETSOURCECONS(CHF_FRA(sourceFAB),
                             CHF_CONST_FRA(consOldFAB),
                             CHF_CONST_REAL(m_sourceTermScaling),
                             CHF_BOX(sourceFAB.box()));

          sourceFAB *= m_dt;
          consNewFAB += sourceFAB;

          FArrayBox sourceTildeFAB(consNewFAB.box(),consNewFAB.nComp());

          FORT_SETSOURCECONS(CHF_FRA(sourceTildeFAB),
                             CHF_CONST_FRA(consNewFAB),
                             CHF_CONST_REAL(m_sourceTermScaling),
                             CHF_BOX(sourceTildeFAB.box()));

          sourceTildeFAB *= m_dt;
          sourceTildeFAB -= sourceFAB;
          sourceTildeFAB *= 0.5;

          consNewFAB += sourceTildeFAB;
        }
    }

  // Update the time and store the new timestep
  m_time += m_dt;
  Real returnDt = m_cfl * newDt;

  m_dtNew = returnDt;

  return returnDt;
}

// Things to do after a timestep
void AMRLevelPolytropicGas::postTimeStep()
{
  CH_assert(allDefined());

  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::postTimeStep " << m_level << endl;
    }

  if (m_hasFiner)
    {
      // Reflux
      Real scale = -1.0/m_dx;
      m_fluxRegister.reflux(m_UNew,scale);

      // Average from finer level data
      AMRLevelPolytropicGas* amrGodFinerPtr = getFinerLevel();

      amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                      amrGodFinerPtr->m_UNew);
    }

  if (s_verbosity >= 2 && m_level == 0)
    {
      int nRefFine = 1;

      pout() << "AMRLevelPolyTropicGas::postTimeStep:" << endl;
      pout() << "  Sums:" << endl;
      for (int comp = 0; comp < m_numStates; comp++)
        {
          Interval curComp(comp,comp);
          Real integral = computeSum(m_UNew,NULL,nRefFine,m_dx,curComp);

          pout() << "    " << setw(23)
                 << setprecision(16)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << integral
                 << " --- " << m_stateNames[comp];

          if (comp == 0 && !first)
            {
              pout() << " (" << setw(23)
                     << setprecision(16)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << (integral-last_integral)
                     << " " << setw(23)
                     << setprecision(16)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << (integral-orig_integral)
                     << ")";
            }

          pout() << endl;

          if (comp == 0)
            {
              if (first)
                {
                  orig_integral = integral;
                  first = false;
                }

              last_integral = integral;
            }
        }
    }

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::postTimeStep " << m_level << " finished" << endl;
    }
}

// Create tags for regridding
void AMRLevelPolytropicGas::tagCells(IntVectSet& a_tags)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::tagCells " << m_level << endl;
    }

  // Create tags based on undivided gradient of density
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  IntVectSet localTags;

  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
    {
      const AMRLevelPolytropicGas* amrGodCoarserPtr = getCoarserLevel();

      PiecewiseLinearFillPatch pwl(levelDomain,
                                   amrGodCoarserPtr->m_UNew.disjointBoxLayout(),
                                   m_numStates,
                                   amrGodCoarserPtr->m_problem_domain,
                                   amrGodCoarserPtr->m_ref_ratio,
                                   1);

      pwl.fillInterp(m_UNew,
                     amrGodCoarserPtr->m_UNew,
                     amrGodCoarserPtr->m_UNew,
                     1.0,
                     0,
                     0,
                     m_numStates);
    }
  m_UNew.exchange(Interval(0,m_numStates-1));

  // Compute relative gradient
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

          FORT_GETRELGRADF(CHF_FRA1(gradFab,dir),
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

// Create tags at initialization
// Since tags are calculated using only current time step data, use
// the same tagging function for initialization and for regridding.
void AMRLevelPolytropicGas::tagCellsInit(IntVectSet& a_tags)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::tagCellsInit " << m_level << endl;
    }

  tagCells(a_tags);
}

// Set up data on this level after regridding
void AMRLevelPolytropicGas::regrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::regrid " << m_level << endl;
    }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  if (s_verbosity >= 4)
    {
      // Indicate/guarantee that the indexing below is only for reading
      // otherwise an error/assertion failure occurs
      const DisjointBoxLayout& constGrids = m_grids;

      pout() << "new grids: " << endl;

      for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
        {
          pout() << constGrids[lit()] << endl;
        }
    }

  // Save data for later
  DataIterator dit = m_UNew.dataIterator();
  for (; dit.ok(); ++dit)
    {
      m_UOld[dit()].copy(m_UNew[dit()]);
    }

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost * IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);

  // Set up data structures
  levelSetup();

  // Interpolate from coarser level
  if (m_hasCoarser)
    {
      AMRLevelPolytropicGas* amrGodCoarserPtr = getCoarserLevel();
      m_fineInterp.interpToFine(m_UNew,amrGodCoarserPtr->m_UNew);
    }

  // Copy from old state
  m_UOld.copyTo(m_UOld.interval(),
                m_UNew,
                m_UNew.interval());

  m_UOld.define(m_grids,m_numStates,ivGhost);
}

// Initialize grids
void AMRLevelPolytropicGas::initialGrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::initialGrid " << m_level << endl;
    }

  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);

  if (s_verbosity >= 4)
    {
      // Indicate/guarantee that the indexing below is only for reading
      // otherwise an error/assertion failure occurs
      const DisjointBoxLayout& constGrids = m_grids;

      pout() << "new grids: " << endl;
      for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
        {
          pout() << constGrids[lit()] << endl;
        }
    }

  // Define old and new state data structures
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  m_UOld.define(m_grids,m_numStates,ivGhost);

  // Set up data structures
  levelSetup();
}

// Initialize data
void AMRLevelPolytropicGas::initialData()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::initialData " << m_level << endl;
    }

  PhysIBC* physIBCPtr = m_gdnvPhysics->getPhysIBC();
  physIBCPtr->initialize(m_UNew);
}

// Things to do after initialization
void AMRLevelPolytropicGas::postInitialize()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::postInitialize " << m_level << endl;
    }

  if (m_hasFiner)
    {
      // Volume weighted average from finer level data
      AMRLevelPolytropicGas* amrGodFinerPtr = getFinerLevel();

      amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                      amrGodFinerPtr->m_UNew);
    }
}

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelPolytropicGas::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::writeCheckpointHeader" << endl;
    }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
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

// Write checkpoint data for this level
void AMRLevelPolytropicGas::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::writeCheckpointLevel" << endl;
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
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM6(if (m_problem_domain.isPeriodic(0))
            {
              header.m_int ["is_periodic_0"] = 1;
            }
          else
            {
              header.m_int ["is_periodic_0"] = 0;
            } ,

          if (m_problem_domain.isPeriodic(1))
            {
              header.m_int ["is_periodic_1"] = 1;
            }
          else
            {
              header.m_int ["is_periodic_1"] = 0;
            } ,

          if (m_problem_domain.isPeriodic(2))
            {
              header.m_int ["is_periodic_2"] = 1;
            }
          else
            {
              header.m_int ["is_periodic_2"] = 0;
            } ,

          if (m_problem_domain.isPeriodic(3))
            {
              header.m_int ["is_periodic_3"] = 1;
            }
          else
            {
              header.m_int ["is_periodic_3"] = 0;
            } ,

          if (m_problem_domain.isPeriodic(4))
            {
              header.m_int ["is_periodic_4"] = 1;
            }
          else
            {
              header.m_int ["is_periodic_4"] = 0;
            } ,

          if (m_problem_domain.isPeriodic(5))
            {
              header.m_int ["is_periodic_5"] = 1;
            }
          else
            {
              header.m_int ["is_periodic_5"] = 0;
            } );

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

// Read checkpoint header
void AMRLevelPolytropicGas::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::readCheckpointHeader" << endl;
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
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointHeader: checkpoint file does not have num_components");
    }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates)
    {
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointHeader: num_components in checkpoint file does not match solver");
    }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
      {
        MayDay::Error("AMRLevelPolytropicGas::readCheckpointHeader: checkpoint file does not have enough component names");
      }

    stateName = header.m_string[compStr];
    if (stateName != m_stateNames[comp])
      {
        MayDay::Error("AMRLevelPolytropicGas::readCheckpointHeader: state_name in checkpoint does not match solver");
      }
  }
}

// Read checkpoint data for this level
void AMRLevelPolytropicGas::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::readCheckpointLevel" << endl;
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
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain ref_ratio");
    }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
    {
      pout() << "read ref_ratio = " << m_ref_ratio << endl;
    }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
    {
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain tag_buffer_size");
    }
  m_tagBufferSize = header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
    {
      pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
    }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
    {
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain dx");
    }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
    {
      pout() << "read dx = " << m_dx << endl;
    }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
    {
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain dt");
    }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
    {
      pout() << "read dt = " << m_dt << endl;
    }

  // Get time
  if (header.m_real.find("time") == header.m_real.end())
    {
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain time");
    }
  m_time = header.m_real["time"];

  if (s_verbosity >= 2)
    {
      pout() << "read time = " << m_time << endl;
    }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
    {
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain prob_domain");
    }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM6(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
            {
              isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
            }
          else
            {
              isPeriodic[0] = false;
            } ,

          if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
            {
              isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
            }
          else
            {
              isPeriodic[1] = false;
            } ,

          if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
            {
              isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
            }
          else
            {
              isPeriodic[2] = false;
            } ,

          if (!(header.m_int.find("is_periodic_3") == header.m_int.end()))
            {
              isPeriodic[3] =  (header.m_int["is_periodic_3"] == 1);
            }
          else
            {
              isPeriodic[3] = false;
            } ,

          if (!(header.m_int.find("is_periodic_4") == header.m_int.end()))
            {
              isPeriodic[4] =  (header.m_int["is_periodic_4"] == 1);
            }
          else
            {
              isPeriodic[4] = false;
            } ,

          if (!(header.m_int.find("is_periodic_5") == header.m_int.end()))
            {
              isPeriodic[5] =  (header.m_int["is_periodic_5"] == 1);
            }
          else
            {
              isPeriodic[5] = false;
            } );

  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);

  if (gridStatus != 0)
    {
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain a Vector<Box>");
    }

  // Create level domain
  m_grids = loadBalance(grids);

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
  m_UNew.define(m_grids,m_numStates);
  const int dataStatus = read<FArrayBox>(a_handle,
                                         m_UNew,
                                         "data",
                                         m_grids);

  if (dataStatus != 0)
    {
      MayDay::Error("AMRLevelPolytropicGas::readCheckpointLevel: file does not contain state data");
    }
  m_UOld.define(m_grids,m_numStates);

  // Set up data structures
  levelSetup();
}

// Write plotfile header
void AMRLevelPolytropicGas::writePlotHeader(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::writePlotHeader" << endl;
    }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      header.m_string[compStr] = m_stateNames[comp];
    }

  // Write the header
  header.writeToFile(a_handle);
  //  a_handle.setGroup("/Expressions");
//  HDF5HeaderData expressions;
//  m_levelGodunov.getGodunovPhysicsPtrConst()->expressions(expressions);
//  expressions.writeToFile(a_handle);

  if (s_verbosity >= 3)
    {
      pout() << header << endl;
    }
}

// Write plotfile data for this level
void AMRLevelPolytropicGas::writePlotLevel(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::writePlotLevel" << endl;
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
  write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_UNew,"data", IntVect::Unit);
}

#endif

// Returns the dt computed earlier for this level
Real AMRLevelPolytropicGas::computeDt()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::computeDt " << m_level << endl;
    }

  Real newDt;
  newDt = m_dtNew;

  return newDt;
}

// Compute dt using initial data
Real AMRLevelPolytropicGas::computeInitialDt()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::computeInitialDt " << m_level << endl;
    }

  Real newDT = m_initial_dt_multiplier * m_dx / m_levelGodunov.getMaxWaveSpeed(m_UNew);

  return newDT;
}

const LevelData<FArrayBox>& AMRLevelPolytropicGas::getStateNew() const
{
  CH_assert(allDefined());

  return m_UNew;
}

const LevelData<FArrayBox>& AMRLevelPolytropicGas::getStateOld() const
{
  CH_assert(allDefined());

  return m_UOld;
}

bool AMRLevelPolytropicGas::allDefined() const
{
  return isDefined()     &&
         m_paramsDefined ;
}

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelPolytropicGas::loadBalance(const Vector<Box>& a_grids)
{
  CH_assert(allDefined());

  // Load balance and create boxlayout
  Vector<int> procMap;

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
    {
      pout() << "AMRLevelPolytropicGas::loadBalance: procesor map: " << endl;

      for (int igrid = 0; igrid < a_grids.size(); ++igrid)
        {
          pout() << igrid << ": " << procMap[igrid] << "  " << endl;
        }

      pout() << endl;
    }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}

// Setup menagerie of data structures
void AMRLevelPolytropicGas::levelSetup()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::levelSetup " << m_level << endl;
    }

  AMRLevelPolytropicGas* amrGodCoarserPtr = getCoarserLevel();
  AMRLevelPolytropicGas* amrGodFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrGodCoarserPtr != NULL);
  m_hasFiner   = (amrGodFinerPtr   != NULL);

  if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();

      m_coarseAverage.define(m_grids,
                             m_numStates,
                             nRefCrse);

      m_fineInterp.define(m_grids,
                          m_numStates,
                          nRefCrse,
                          m_problem_domain);

      const DisjointBoxLayout& coarserLevelDomain = amrGodCoarserPtr->m_grids;

      // Maintain levelGodunov
      m_levelGodunov.define(m_grids,
                            coarserLevelDomain,
                            m_problem_domain,
                            nRefCrse,
                            m_dx,
                            m_gdnvPhysics,
                            m_normalPredOrder,
                            m_useFourthOrderSlopes,
                            m_usePrimLimiting,
                            m_useCharLimiting,
                            m_useFlattening,
                            m_useArtificialViscosity,
                            m_artificialViscosity,
                            m_hasCoarser,
                            m_hasFiner);
      m_levelGodunov.highOrderLimiter(m_highOrderLimiter);

      // This may look twisted but you have to do this this way because the
      // coarser levels get setup before the finer levels so, since a flux
      // register lives between this level and the next FINER level, the finer
      // level has to do the setup because it is the only one with the
      // information at the time of construction.

      // Maintain flux registers
      amrGodCoarserPtr->m_fluxRegister.define(m_grids,
                                              amrGodCoarserPtr->m_grids,
                                              m_problem_domain,
                                              amrGodCoarserPtr->m_ref_ratio,
                                              m_numStates);
      amrGodCoarserPtr->m_fluxRegister.setToZero();
    }
  else
    {
      m_levelGodunov.define(m_grids,
                            DisjointBoxLayout(),
                            m_problem_domain,
                            m_ref_ratio,
                            m_dx,
                            m_gdnvPhysics,
                            m_normalPredOrder,
                            m_useFourthOrderSlopes,
                            m_usePrimLimiting,
                            m_useCharLimiting,
                            m_useFlattening,
                            m_useArtificialViscosity,
                            m_artificialViscosity,
                            m_hasCoarser,
                            m_hasFiner);
      m_levelGodunov.highOrderLimiter(m_highOrderLimiter);
    }
}

// Get the next coarser level
AMRLevelPolytropicGas* AMRLevelPolytropicGas::getCoarserLevel() const
{
  CH_assert(allDefined());

  AMRLevelPolytropicGas* amrGodCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
    {
      amrGodCoarserPtr = dynamic_cast<AMRLevelPolytropicGas*>(m_coarser_level_ptr);

      if (amrGodCoarserPtr == NULL)
        {
          MayDay::Error("AMRLevelPolytropicGas::getCoarserLevel: dynamic cast failed");
        }
    }

  return amrGodCoarserPtr;
}

// Get the next finer level
AMRLevelPolytropicGas* AMRLevelPolytropicGas::getFinerLevel() const
{
  CH_assert(allDefined());

  AMRLevelPolytropicGas* amrGodFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
    {
      amrGodFinerPtr = dynamic_cast<AMRLevelPolytropicGas*>(m_finer_level_ptr);

      if (amrGodFinerPtr == NULL)
        {
          MayDay::Error("AMRLevelPolytropicGas::getFinerLevel: dynamic cast failed");
        }
    }

  return amrGodFinerPtr;
}

#include "NamespaceFooter.H"
