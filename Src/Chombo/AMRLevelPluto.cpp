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
#include "CH_HDF5.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "computeSum.H"
#include "AMRIO.H"
#include "AMRLevel.H"

#include "AMRLevelPluto.H"

#include "NamespaceHeader.H"

// Constructor
AMRLevelPluto::AMRLevelPluto()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto default constructor" << endl;
  }

  m_patchPlutoFactory = NULL;
  m_patchPluto = NULL;
  m_paramsDefined = false;
}

// Destructor
AMRLevelPluto::~AMRLevelPluto()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto destructor" << endl;
  }

  // Get rid of the patch integrator and its factory
  if (m_patchPluto != NULL)
  {
    delete m_patchPluto;
  }
  
  if (m_patchPlutoFactory != NULL)
  {
    delete m_patchPlutoFactory;
  }

  m_paramsDefined = false;
}

void  AMRLevelPluto::defineParams(const Real&              a_cfl,
                                  const Real&              a_domainLength,
                                  const int&               a_verbosity,
                                  const Real&              a_refineThresh,
                                  const int&               a_tagBufferSize,
                                  const Real&              a_initialDtMultiplier,
                                  const PatchPluto*  const a_patchPluto)
{
  // Store the CFL number
  m_cfl = a_cfl;

  // Store the physical dimension of the longest side of the domain
  m_domainLength = a_domainLength;

  // Store the verbosity of the object
  verbosity(a_verbosity);

  // Store the refinement threshold for gradient
  m_refineThresh = a_refineThresh;

  // Store the tag buffer size
  m_tagBufferSize = a_tagBufferSize;

  // Store the initial dt multiplier
  initialDtMultiplier(a_initialDtMultiplier);

  // Delete any existing patch factory object
  // (Factory has to be kept for the Patch integrator in LevelPluto ...
  //  ... go figure ...)
  if (m_patchPlutoFactory != NULL)
  {
    delete m_patchPlutoFactory;
  }

  m_patchPlutoFactory = a_patchPluto->new_patchPluto();

  // Delete any existing patch object
  if (m_patchPlutoFactory != NULL)
  {
    delete m_patchPluto;
  }

  m_patchPluto = a_patchPluto->new_patchPluto();

  m_paramsDefined = true;
}

// This instance should never get called - historical
void AMRLevelPluto::define(AMRLevel*  a_coarserLevelPtr,
                           const Box& a_problemDomain,
                           int        a_level,
                           int        a_refRatio)
{
  ProblemDomain physdomain(a_problemDomain);

  MayDay::Error("AMRLevelPluto::define -\n\tShould never be called with a Box for a problem domain");
}

// Define new AMR level
void AMRLevelPluto::define(AMRLevel*            a_coarserLevelPtr,
                           const ProblemDomain& a_problemDomain,
                           int                  a_level,
                           int                  a_refRatio)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::define " << a_level << endl;
  }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelPluto* amrGodPtr = dynamic_cast<AMRLevelPluto*>(a_coarserLevelPtr);

    if (amrGodPtr != NULL)
    {
      m_cfl = amrGodPtr->m_cfl;
      m_domainLength = amrGodPtr->m_domainLength;
      m_refineThresh = amrGodPtr->m_refineThresh;
      m_tagBufferSize = amrGodPtr->m_tagBufferSize;
    }
    else
    {
      MayDay::Error("AMRLevelPluto::define: a_coarserLevelPtr is not castable to AMRLevelPluto*");
    }
  }

  // Compute the grid spacing
  m_dx = m_domainLength / (a_problemDomain.domainBox().bigEnd(0)-a_problemDomain.domainBox().smallEnd(0)+1.);

  // Nominally, one layer of ghost cells is maintained permanently and
  // individual computations may create local data with more
  m_numGhost = 1;

  CH_assert(m_patchPluto != NULL);
  CH_assert(isDefined());
  m_patchPluto->define(m_problem_domain,m_dx,m_level,m_numGhost);

  // Get additional information from the patch integrator
  m_numStates      = m_patchPluto->numConserved();
  m_ConsStateNames = m_patchPluto->ConsStateNames();
  m_PrimStateNames = m_patchPluto->PrimStateNames();
}

// Advance by one timestep
Real AMRLevelPluto::advance()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::advance level " << m_level << " to time " << m_time + m_dt << endl;
  }

  // Copy the new to the old

  DataIterator dit = m_UNew.dataIterator();
  for(;dit.ok(); ++dit){
	m_UOld[dit()].copy(m_UNew[dit()]);
  }

  Real newDt = 0.0;

  // Set up arguments to LevelPluto::step based on whether there are
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
    AMRLevelPluto* coarserPtr = getCoarserLevel();

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
  
  // we don't need the flux in the simple hyperbolic case...
  LevelData<FArrayBox> flux[SpaceDim];

  g_intStage = 1;

  // Advance the solve one timestep
  newDt = m_levelPluto.step(m_UNew,
                            flux,
                            *finerFR,
                            *coarserFR,
                            m_split_tags,
                            *coarserDataOld,
                            tCoarserOld,
                            *coarserDataNew,
                            tCoarserNew,
                            m_time,
                            m_dt,
                            m_cfl);

#if (TIME_STEPPING == RK2)
  g_intStage = 2;
  Real DtCool; // The predictor returns the advective/diffusive timestep
               // The corrector returns the cooling timestep

  DtCool = m_levelPluto.step(m_UNew,
                            flux,
                            *finerFR,
                            *coarserFR,
                            m_split_tags,
                            *coarserDataOld,
                            tCoarserOld,
                            *coarserDataNew,
                            tCoarserNew,
                            m_time+m_dt,
                            m_dt,
                            m_cfl);

  #if (COOLING != NO)
   newDt = Min(newDt,DtCool);
  #endif
#endif

  m_time += m_dt;                  // Update m_time                                       
  Real returnDt = m_cfl * newDt;   // Store the new timestep

  m_dtNew = returnDt;

  return returnDt;
}

Real AMRLevelPluto::getDlMin()
{

  Real dlMin = m_levelPluto.getDlMin();

  return dlMin;

}

// Things to do after a timestep
void AMRLevelPluto::postTimeStep()
{
  CH_assert(allDefined());

  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::postTimeStep " << m_level << endl;
  }

  if (m_hasFiner)
  {
    // Reflux
    Real scale = -1.0/m_dx; 
    m_fluxRegister.reflux(m_UNew,scale);

    // Average from finer level data
    AMRLevelPluto* amrGodFinerPtr = getFinerLevel();

    amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                    amrGodFinerPtr->m_UNew);
  }

  if (s_verbosity >= 2 && m_level == 0)
  {
    int nRefFine = 1;

    pout() << "AMRLevelPluto::postTimeStep:" << endl;
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
             << " --- " << m_ConsStateNames[comp];

      if (comp == 0 && !first) {
        pout() << " (" << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << (integral-last_integral)/last_integral
               << " " << setw(23)
                      << setprecision(16)
                      << setiosflags(ios::showpoint)
                      << setiosflags(ios::scientific)
                      << (integral-orig_integral)/orig_integral
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
    pout() << "AMRLevelPluto::postTimeStep " << m_level << " finished" << endl;
  }
}

// Create tags for regridding
void AMRLevelPluto::tagCells(IntVectSet& a_tags) 
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::tagCells " << m_level << endl;
  }

  // Create tags based on undivided gradient of density
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  IntVectSet localTags;
  // If there is a coarser level interpolate undefined ghost cells
  if (m_hasCoarser)
  {
    const AMRLevelPluto* amrGodCoarserPtr = getCoarserLevel();

    PiecewiseLinearFillPluto pwl;

    pwl.define(levelDomain,
               amrGodCoarserPtr->m_UNew.disjointBoxLayout(),
               m_numStates,
               amrGodCoarserPtr->m_problem_domain,
               amrGodCoarserPtr->m_ref_ratio,
               m_dx,
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

  #if GEOMETRY != CARTESIAN
   const LevelData<FArrayBox>& dV = m_levelPluto.getdV();
  #endif

  // Compute relative gradient
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit){
    const Box& b = levelDomain[dit()];
    FArrayBox gradFab(b,1);
    FArrayBox& UFab = m_UNew[dit()];

    #if GEOMETRY != CARTESIAN
     const FArrayBox& curdV = dV[dit()];
    #else
     const FArrayBox curdV;
    #endif

    m_patchPluto->computeRefGradient(gradFab, UFab, curdV, b); 

    // Tag where gradient exceeds threshold
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit){
      const IntVect& iv = bit();
      if (gradFab(iv) >= m_refineThresh) {
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
void AMRLevelPluto::tagCellsInit(IntVectSet& a_tags) 
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
    {
      pout() << "AMRLevelPolytropicGas::tagCellsInit " << m_level << endl;
    }

  tagCells(a_tags);
}

// Set up data on this level after regridding
void AMRLevelPluto::regrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::regrid " << m_level << endl;
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
  for(;dit.ok(); ++dit){
	m_UOld[dit()].copy(m_UNew[dit()]);
  }

  // Reshape sweep with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);


  // Set up data structures
  levelSetup();

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    AMRLevelPluto* amrGodCoarserPtr = getCoarserLevel();
    m_fineInterp.interpToFine(m_UNew,amrGodCoarserPtr->m_UNew);
  }

  // Copy from old sweep
  m_UOld.copyTo(m_UOld.interval(),
                m_UNew,
                m_UNew.interval());

  m_UOld.define(m_grids,m_numStates,ivGhost);
} 

// Initialize grids
void AMRLevelPluto::initialGrid(const Vector<Box>& a_newGrids)
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::initialGrid " << m_level << endl;
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

  // Define old and new sweep data structures
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  m_UOld.define(m_grids,m_numStates,ivGhost);

  // Set up data structures
  levelSetup();
}

// Initialize data
void AMRLevelPluto::initialData()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::initialData " << m_level << endl;
  }

  m_patchPluto->initialize(m_UNew);

}

// Things to do after initialization
void AMRLevelPluto::postInitialize()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::postInitialize " << m_level << endl;
  }

  if (m_hasFiner)
  {
    // Volume weighted average from finer level data
    AMRLevelPluto* amrGodFinerPtr = getFinerLevel();

    amrGodFinerPtr->m_coarseAverage.averageToCoarse(m_UNew,
                                                    amrGodFinerPtr->m_UNew);
  }

}

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelPluto::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::writeCheckpointHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_ConsStateNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write checkpoint data for this level
void AMRLevelPluto::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::writeCheckpointLevel" << endl;
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
  header.m_int ["geometry"]        = GEOMETRY;
  header.m_int ["logr"]            = CHOMBO_LOGR;
 #if DIMENSIONS >= 2
  header.m_real["g_x2stretch"]        = g_x2stretch;
 #endif 
 #if DIMENSIONS == 3
  header.m_real["g_x3stretch"]        = g_x3stretch;
 #endif
  header.m_real["domBeg1"]         = g_domBeg[IDIR];
 #if DIMENSIONS >= 2
  header.m_real["domBeg2"]         = g_domBeg[JDIR];
 #endif
 #if DIMENSIONS == 3
  header.m_real["domBeg3"]         = g_domBeg[KDIR];
 #endif
  header.m_real["dt"]              = m_dt;
 #ifdef GLM_MHD
  header.m_real["ch"]              = glm_ch_max;
 #endif
  header.m_real["time"]            = m_time;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM(if (m_problem_domain.isPeriodic(0))
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
         } );

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  #if GEOMETRY != CARTESIAN
   LevelData<FArrayBox> tmp_U;
   IntVect ivGhost = m_numGhost*IntVect::Unit;
   tmp_U.define(m_grids,m_numStates,ivGhost);
   
   const LevelData<FArrayBox>& dV = m_levelPluto.getdV();

   DataIterator dit = m_UNew.dataIterator();
   for(;dit.ok(); ++dit){
    tmp_U[dit()].copy(m_UNew[dit()]);
    FArrayBox& curU = tmp_U[dit()];
    const FArrayBox& curdV = dV[dit()];
    for (int ivar = 0; ivar < curU.nComp(); ivar++) curU.divide(curdV,0,ivar); 
    #if CHOMBO_CONS_AM == YES
     #if ROTATING_FRAME == YES
      Box curBox = curU.box();
      for(BoxIterator bit(curBox); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        curU(iv,iMPHI) /= curdV(iv,1);
        curU(iv,iMPHI) -= curU(iv,RHO)*curdV(iv,1)*g_OmegaZ;
      }
     #else
      curU.divide(curdV,1,iMPHI);
     #endif
    #endif
   }
   write(a_handle,tmp_U.boxLayout());
   write(a_handle,tmp_U,"data");
  #else
   if (g_stretch_fact != 1.) {
    LevelData<FArrayBox> tmp_U;
    IntVect ivGhost = m_numGhost*IntVect::Unit;
    tmp_U.define(m_grids,m_numStates,ivGhost); 

    DataIterator dit = m_UNew.dataIterator();
    for(;dit.ok(); ++dit){
     tmp_U[dit()].copy(m_UNew[dit()]);    
     tmp_U[dit()] /= g_stretch_fact;
    }
    write(a_handle,tmp_U.boxLayout());
    write(a_handle,tmp_U,"data");
   } else {
    write(a_handle,m_UNew.boxLayout());
    write(a_handle,m_UNew,"data");
   }
  #endif

}

// Read checkpoint header
void AMRLevelPluto::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::readCheckpointHeader" << endl;
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
    MayDay::Error("AMRLevelPluto::readCheckpointHeader: checkpoint file does not have num_components");
  }

  int numStates = header.m_int["num_components"];
  if (numStates != m_numStates)
  {
    MayDay::Error("AMRLevelPluto::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }

  // Get the component names
  std::string sweepName;
  char compStr[60];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
    {
      MayDay::Error("AMRLevelPluto::readCheckpointHeader: checkpoint file does not have enough component names");
    }

    sweepName = header.m_string[compStr];
    if (sweepName != m_ConsStateNames[comp])
    {
      MayDay::Error("AMRLevelPluto::readCheckpointHeader: sweep_name in checkpoint does not match solver");
    }
  }
}

// Read checkpoint data for this level
void AMRLevelPluto::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::readCheckpointLevel" << endl;
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
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  {
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain tag_buffer_size");
  }
  m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real["dx"];

  if (s_verbosity >= 2)
  {
    pout() << "read dx = " << m_dx << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }

  // Get g_x2stretch and g_x3stretch factors
#if CH_SPACEDIM >= 2
  if (header.m_real.find("g_x2stretch") == header.m_real.end())
  {
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain g_x2stretch");
  }
  g_x2stretch = header.m_real["g_x2stretch"];

  if (s_verbosity >= 2)
  {
    pout() << "read g_x2stretch = " << g_x2stretch << endl;
  }
#endif
#if CH_SPACEDIM == 3
  if (header.m_real.find("g_x3stretch") == header.m_real.end())
  {
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain g_x3stretch");
  }
  g_x3stretch = header.m_real["g_x3stretch"];

  if (s_verbosity >= 2)
  {
    pout() << "read g_x3stretch = " << g_x3stretch << endl;
  }
#endif

#if GEOMETRY == CARTESIAN
 g_stretch_fact = g_x2stretch*g_x3stretch;
#endif

 #ifdef GLM_MHD
  // Get glm_ch
  if (header.m_real.find("ch") == header.m_real.end())
  {
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain glm_ch");
  }
  glm_ch_max = header.m_real["ch"];

  if (s_verbosity >= 2)
  {
    pout() << "read ch max = " << glm_ch_max << endl;
  }
 #endif

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain prob_domain");
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
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain a Vector<Box>");
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

  // Reshape sweep with new grids
  m_UNew.define(m_grids,m_numStates);
  const int dataStatus = read<FArrayBox>(a_handle,
                                         m_UNew,
                                         "data",
                                         m_grids);

  if (dataStatus != 0)
  {
    MayDay::Error("AMRLevelPluto::readCheckpointLevel: file does not contain sweep data");
  }

  m_UOld.define(m_grids,m_numStates);

  // Set up data structures
  levelSetup();

  #if GEOMETRY != CARTESIAN
   const LevelData<FArrayBox>& dV = m_levelPluto.getdV();
   DataIterator dit = m_UNew.dataIterator();
   for(;dit.ok(); ++dit){
    FArrayBox& curU = m_UNew[dit()]; 
    const FArrayBox& curdV = dV[dit()];
    #if CHOMBO_CONS_AM == YES
     #if ROTATING_FRAME == YES
      Box curBox = curU.box();
      for(BoxIterator bit(curBox); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        curU(iv,iMPHI) += curU(iv,RHO)*curdV(iv,1)*g_OmegaZ;
        curU(iv,iMPHI) *= curdV(iv,1);
      }
     #else
      curU.mult(curdV,1,iMPHI);
     #endif 
    #endif
    for (int ivar = 0; ivar < curU.nComp(); ivar++) curU.mult(curdV,0,ivar);
   }
  #else
   if (g_stretch_fact != 1.) {
    DataIterator dit = m_UNew.dataIterator();
    for(;dit.ok(); ++dit) {
     FArrayBox& curU = m_UNew[dit()];
     curU *= g_stretch_fact;
    }
   }
  #endif
}

// Write plotfile header
void AMRLevelPluto::writePlotHeader(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::writePlotHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_PrimStateNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);
  a_handle.setGroup("/Expressions");
  HDF5HeaderData expressions;
  DefineExpressions(expressions);
  expressions.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write plotfile data for this level
void AMRLevelPluto::writePlotLevel(HDF5Handle& a_handle) const
{
  CH_assert(allDefined());

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = m_dx;
  header.m_int ["geometry"]    = GEOMETRY;
  header.m_int ["logr"]        = CHOMBO_LOGR;
 #if DIMENSIONS >= 2
  header.m_real["g_x2stretch"]    = g_x2stretch;
 #endif
 #if DIMENSIONS == 3
  header.m_real["g_x3stretch"]    = g_x3stretch;
 #endif
  header.m_real["domBeg1"]     = g_domBeg[IDIR];
 #if DIMENSIONS >= 2
  header.m_real["domBeg2"]     = g_domBeg[JDIR];
 #endif
 #if DIMENSIONS == 3
  header.m_real["domBeg3"]     = g_domBeg[KDIR];
 #endif
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
  LevelData<FArrayBox> tmp_U;
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  tmp_U.define(m_grids,m_numStates,ivGhost);

  #if GEOMETRY != CARTESIAN
   const LevelData<FArrayBox>& dV = m_levelPluto.getdV();
  #endif

  DataIterator dit = m_UNew.dataIterator(); 
  for(;dit.ok(); ++dit){
    tmp_U[dit()].copy(m_UNew[dit()]);
    FArrayBox& curU = tmp_U[dit()];
    #if GEOMETRY != CARTESIAN
     const FArrayBox& curdV = dV[dit()];
     for (int ivar = 0; ivar < curU.nComp(); ivar++) curU.divide(curdV,0,ivar); 
     #if CHOMBO_CONS_AM == YES
      #if ROTATING_FRAME == YES
       Box curBox = curU.box();
       for(BoxIterator bit(curBox); bit.ok(); ++bit) {
         const IntVect& iv = bit();
         curU(iv,iMPHI) /= curdV(iv,1);
         curU(iv,iMPHI) -= curU(iv,RHO)*curdV(iv,1)*g_OmegaZ;
       }
      #else
       curU.divide(curdV,1,iMPHI); 
      #endif
     #endif
    #else
     if (g_stretch_fact != 1.) curU /= g_stretch_fact;
    #endif

    m_patchPluto->convertFArrayBox(curU); /* -- convert data to primitive -- */
  }
  write(a_handle,tmp_U.boxLayout());
  write(a_handle,tmp_U,"data");

}

void AMRLevelPluto::DefineExpressions(HDF5HeaderData& a_expressions) const
{

 char str_expr[128];

 // Grid/geometry expressions

 #if GEOMETRY == CARTESIAN
   sprintf (str_expr,"coords(Mesh)[0]+nodal_constant(Mesh,%12.6e)",g_domBeg[IDIR]);
   a_expressions.m_string["scalar X"] = str_expr;
  #if DIMENSIONS == 2
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[1]+nodal_constant(Mesh,%12.6e)",g_x2stretch,g_domBeg[JDIR]);
   a_expressions.m_string["scalar Y"] = str_expr;
   a_expressions.m_string["vector Displacement"] = "{X,Y}-coords(Mesh)";
  #endif
  #if DIMENSIONS == 3
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[1]+nodal_constant(Mesh,%12.6e)",g_x2stretch,g_domBeg[JDIR]);
   a_expressions.m_string["scalar Y"] = str_expr;
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[2]+nodal_constant(Mesh,%12.6e)",g_x3stretch,g_domBeg[KDIR]);
   a_expressions.m_string["scalar Z"] = str_expr;
   a_expressions.m_string["vector Displacement"] = "{X,Y,Z}-coords(Mesh)";
  #endif
 #endif

 #if GEOMETRY == CYLINDRICAL
   sprintf (str_expr,"coords(Mesh)[0]+nodal_constant(Mesh,%12.6e)",g_domBeg[IDIR]);
   a_expressions.m_string["scalar R"] = str_expr;
  #if DIMENSIONS == 2
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[1]+nodal_constant(Mesh,%12.6e)",g_x2stretch,g_domBeg[JDIR]);
   a_expressions.m_string["scalar Z"] = str_expr;
   a_expressions.m_string["vector Displacement"] = "{R,Z}-coords(Mesh)";
  #endif  
 #endif 

 #if GEOMETRY == SPHERICAL
  #if CHOMBO_LOGR == NO
   sprintf (str_expr,"coords(Mesh)[0]+nodal_constant(Mesh,%12.6e)",g_domBeg[IDIR]);
   a_expressions.m_string["scalar R"] = str_expr;
  #else 
   sprintf (str_expr,"exp(coords(Mesh)[0])*nodal_constant(Mesh,%12.6e)",g_domBeg[IDIR]);
   a_expressions.m_string["scalar R"] = str_expr;
  #endif
  #if DIMENSIONS == 2
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[1]+nodal_constant(Mesh,%12.6e)",g_x2stretch,g_domBeg[JDIR]);
   a_expressions.m_string["scalar Theta"] = str_expr;
   a_expressions.m_string["scalar X"] = "R*sin(Theta)";
   a_expressions.m_string["scalar Z"] = "R*cos(Theta)";
   a_expressions.m_string["vector Displacement"] = "{X,Z}-coords(Mesh)"; 
  #endif
  #if DIMENSIONS == 3
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[1]+nodal_constant(Mesh,%12.6e)",g_x2stretch,g_domBeg[JDIR]);
   a_expressions.m_string["scalar Theta"] = str_expr;
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[2]+nodal_constant(Mesh,%12.6e)",g_x3stretch,g_domBeg[KDIR]);
   a_expressions.m_string["scalar Phi"] = str_expr;
   a_expressions.m_string["scalar X"] = "R*sin(Theta)*cos(Phi)";
   a_expressions.m_string["scalar Y"] = "R*sin(Theta)*sin(Phi)";
   a_expressions.m_string["scalar Z"] = "R*cos(Theta)";
   a_expressions.m_string["vector Displacement"] = "{X,Y,Z}-coords(Mesh)";
  #endif
 #endif

 #if GEOMETRY == POLAR
  #if CHOMBO_LOGR == NO
   sprintf (str_expr,"coords(Mesh)[0]+nodal_constant(Mesh,%12.6e)",g_domBeg[IDIR]);
   a_expressions.m_string["scalar R"] = str_expr;
  #else
   sprintf (str_expr,"exp(coords(Mesh)[0])*nodal_constant(Mesh,%12.6e)",g_domBeg[IDIR]);
   a_expressions.m_string["scalar R"] = str_expr;
  #endif
  #if DIMENSIONS == 2
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[1]+nodal_constant(Mesh,%12.6e)",g_x2stretch,g_domBeg[JDIR]);
   a_expressions.m_string["scalar Phi"] = str_expr;
   a_expressions.m_string["scalar X"] = "R*cos(Phi)";
   a_expressions.m_string["scalar Y"] = "R*sin(Phi)";
   a_expressions.m_string["vector Displacement"] = "{X,Y}-coords(Mesh)";
  #endif
  #if DIMENSIONS == 3
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[1]+nodal_constant(Mesh,%12.6e)",g_x2stretch,g_domBeg[JDIR]);
   a_expressions.m_string["scalar Phi"] = str_expr;
   a_expressions.m_string["scalar X"] = "R*cos(Phi)";
   a_expressions.m_string["scalar Y"] = "R*sin(Phi)";
   sprintf (str_expr,"nodal_constant(Mesh,%12.6e)*coords(Mesh)[2]+nodal_constant(Mesh,%12.6e)",g_x3stretch,g_domBeg[KDIR]);
   a_expressions.m_string["scalar Z"] = str_expr;
   a_expressions.m_string["vector Displacement"] = "{X,Y,Z}-coords(Mesh)";
  #endif
 #endif

}

#endif

// Returns the dt computed earlier for this level
Real AMRLevelPluto::computeDt()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::computeDt " << m_level << endl;
  }

  Real newDt;
  newDt = m_dtNew;

  return newDt;
}

// Compute dt using initial data
Real AMRLevelPluto::computeInitialDt()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::computeInitialDt " << m_level << endl;
  }

  Real newDT = m_initial_dt_multiplier*m_dx/m_domainLength;

  return newDT;
}

const LevelData<FArrayBox>& AMRLevelPluto::getStateNew() const
{
  CH_assert(allDefined());

  return m_UNew;
}

const LevelData<FArrayBox>& AMRLevelPluto::getStateOld() const
{
  CH_assert(allDefined());

  return m_UOld;
}

bool AMRLevelPluto::allDefined() const
{
  return isDefined()     &&
         m_paramsDefined ;
}

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelPluto::loadBalance(const Vector<Box>& a_grids)
{
  CH_assert(allDefined());

  // Load balance and create boxlayout
  Vector<int> procMap;

// use the four-argument call to LoadBalance to have more control over it.
// this may be useful for dynamic load balancing (e.g. when cooling is 
// employed.

  Vector<long long> computeLoads(a_grids.size());

  for (int i = 0; i < a_grids.size(); ++i) {
//  if (i < 7) computeLoads[i] = a_grids[i].numPts()*(long long)100;
    computeLoads[i] = a_grids[i].numPts();
  }

  LoadBalance(procMap, computeLoads, a_grids, numProc());

  // appears to be faster for all procs to do the loadbalance (ndk)
//  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelPluto::loadBalance (level "<<m_level<<"): procesor map: " << endl;

    for (int igrid = 0; igrid < a_grids.size(); ++igrid)
    {
      pout() << igrid << "(load = " << computeLoads[igrid] <<  
                         "): " << procMap[igrid] << "  " << endl;
    }
    pout() << endl;
  }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}

// Setup menagerie of data structures
void AMRLevelPluto::levelSetup()
{
  CH_assert(allDefined());

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelPluto::levelSetup " << m_level << endl;
  }
 
  AMRLevelPluto* amrGodCoarserPtr = getCoarserLevel();
  AMRLevelPluto* amrGodFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrGodCoarserPtr != NULL);
  m_hasFiner   = (amrGodFinerPtr   != NULL);


#ifdef SKIP_SPLIT_CELLS
  // Mark split/unsplit cells
  m_split_tags.define(m_grids,1);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit) {
       m_split_tags[dit()].setVal(1.);
  }
#endif

  if (m_hasCoarser)
  {
    int nRefCrse = m_coarser_level_ptr->refRatio();
  
    m_coarseAverage.define(m_grids,
                           m_numStates,
                           nRefCrse);

    m_fineInterp.define(m_grids,
                        m_numStates,
                        nRefCrse,
                        m_dx,
                        m_problem_domain);

    const DisjointBoxLayout& coarserLevelDomain = amrGodCoarserPtr->m_grids;

#ifdef SKIP_SPLIT_CELLS
    // Mark split/unsplit cells of the coarser level
    amrGodCoarserPtr->mark_split(m_grids);
#endif

    // Maintain levelPluto
    m_levelPluto.define(m_grids,
                        coarserLevelDomain,
                        m_problem_domain,
                        nRefCrse,
                        m_level,
                        m_dx,
                        m_patchPlutoFactory,
                        m_hasCoarser,
                        m_hasFiner);

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
    m_levelPluto.define(m_grids,
                        DisjointBoxLayout(),
                        m_problem_domain,
                        m_ref_ratio,
                        m_level,
                        m_dx,
                        m_patchPlutoFactory,
                        m_hasCoarser,
                        m_hasFiner);
  }
}

// Get the next coarser level
AMRLevelPluto* AMRLevelPluto::getCoarserLevel() const
{
  CH_assert(allDefined());

  AMRLevelPluto* amrGodCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrGodCoarserPtr = dynamic_cast<AMRLevelPluto*>(m_coarser_level_ptr);

    if (amrGodCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelPluto::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrGodCoarserPtr;
}

// Get the next finer level
AMRLevelPluto* AMRLevelPluto::getFinerLevel() const
{
  CH_assert(allDefined());

  AMRLevelPluto* amrGodFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrGodFinerPtr = dynamic_cast<AMRLevelPluto*>(m_finer_level_ptr);

    if (amrGodFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelPluto::getFinerLevel: dynamic cast failed");
    }
  }

  return amrGodFinerPtr;
}

// Mark Split/unsplit cells of this level
void AMRLevelPluto::mark_split(const DisjointBoxLayout& finerLevelDomain)
{
      CH_assert(allDefined());

      DisjointBoxLayout finerDomain;

      coarsen(finerDomain,finerLevelDomain,m_ref_ratio);
  
      DataIterator dit = m_grids.dataIterator(); 
      LayoutIterator lit = finerDomain.layoutIterator();
      for (dit.begin(); dit.ok(); ++dit)  
      {
       FArrayBox& splitU = m_split_tags[dit()];
       splitU.setVal(1.);
       Box splitBox = splitU.box();

         for (lit.begin(); lit.ok(); ++lit)
         {
          Box finBox = finerDomain.get(lit());
          const Box intBox = splitBox & finBox;

          BoxIterator bit(intBox);
          for (bit.begin(); bit.ok(); ++bit)
          {
           const IntVect& iv = bit();
           splitU(iv) = 0.0;
          }
         }
      }

}

#include "NamespaceFooter.H"
