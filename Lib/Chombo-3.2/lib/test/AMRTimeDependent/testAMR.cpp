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
#include <iostream>
using std::endl;

#include "BRMeshRefine.H"
#include "AMRLevel.H"
#include "AMR.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "FArrayBox.H"
#include "LevelData.H"
#include "LayoutIterator.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testAMR" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

class AMRDerivedClass : public AMRLevel
{
public:
  AMRDerivedClass();

  virtual
  ~AMRDerivedClass();

// define
  virtual
  void
  define(AMRLevel* a_coarser_level_ptr,
         const Box& a_problem_domain,
         int a_level,
         int a_ref_ratio);

// define
  virtual
  void
  define(AMRLevel* a_coarser_level_ptr,
         const ProblemDomain& a_problem_domain,
         int a_level,
         int a_ref_ratio);

// advance by one timestep
  virtual
  Real
  advance();

// things to do after a timestep
  virtual
  void
  postTimeStep();

// create tags
  virtual
  void
  tagCells(IntVectSet& a_tags) ;

// create tags at initialization
  virtual
  void
  tagCellsInit(IntVectSet& a_tags) ;

// regrid
  virtual
  void
  regrid(const Vector<Box>& a_new_grids);

// initialize grids
  virtual
  void
  initialGrid(const Vector<Box>& a_new_grids);

// initialize data
  virtual
  void
  initialData();

// things to do after initialization
  virtual
  void
  postInitialize();

// set data after regrid
// (regridData should be removed)
  virtual
  void
  regridData();

#ifdef CH_USE_HDF5
  virtual
  void
  writeCheckpointHeader(HDF5Handle& a_handle) const;

  virtual
  void
  writeCheckpointLevel(HDF5Handle& a_handle) const;

  virtual
  void
  readCheckpointHeader(HDF5Handle& a_handle);

  virtual
  void
  readCheckpointLevel(HDF5Handle& a_handle);

  virtual
  void
  writePlotHeader(HDF5Handle& a_handle) const;

  virtual
  void
  writePlotLevel(HDF5Handle& a_handle) const;
#endif

// compute dt
  virtual
  Real
  computeDt();

// compute dt with initial data
  virtual
  Real
  computeInitialDt();

protected:
  DisjointBoxLayout
  loadBalance(const Vector<Box>& a_grids);

protected:
// state vector at old time
  LevelData<FArrayBox> m_state_old;
// state vector at new time
  LevelData<FArrayBox> m_state_new;
// number of components of m_state
  static const int s_num_comps = 2;
// names of components
  static const char* s_state_names[s_num_comps];
// grid spacing
  Real m_dx;

  CoarseAverage m_coarse_average;
  FineInterp m_fine_interp;
};

const int
AMRDerivedClass::s_num_comps;

const char*
AMRDerivedClass::s_state_names[s_num_comps] =
{
  "blah",
  "blarg"
};

AMRDerivedClass::AMRDerivedClass ()
{
  if ( verbose ) pout () << "AMRDerivedClass default constructor" << endl;
}

AMRDerivedClass::~AMRDerivedClass ()
{
}

void
AMRDerivedClass::define (AMRLevel* a_coarser_level_ptr,
                        const Box& a_problem_domain,
                        int a_level,
                        int a_ref_ratio)
{
  ProblemDomain physdomain(a_problem_domain);
  define(a_coarser_level_ptr, a_problem_domain, a_level, a_ref_ratio);
}

void
AMRDerivedClass::define (AMRLevel* a_coarser_level_ptr,
                        const ProblemDomain& a_problem_domain,
                        int a_level,
                        int a_ref_ratio)
{
  if ( verbose ) pout () << "AMRDerivedClass::define " << a_level << endl;

  AMRLevel::define (a_coarser_level_ptr,
                    a_problem_domain,
                    a_level,
                    a_ref_ratio);
  m_dx = 1. / a_problem_domain.domainBox().longside ();
}

// advance by one timestep

Real
AMRDerivedClass::advance ()
{
  if ( verbose ) pout () << "AMRDerivedClass::advance " << m_level << " at " << m_time << endl;

  m_state_new.copyTo ( m_state_new.interval (),
                       m_state_old,
                       m_state_old.interval () );

// blahblahblah.step(...)

  m_time += m_dt;

  return (.99 * m_dt);
}

// things to do after a timestep

void
AMRDerivedClass::postTimeStep ()
{
  if ( verbose ) pout () << "AMRDerivedClass::postTimeStep " << m_level << endl;

  if (m_finer_level_ptr != NULL)
  {
    AMRDerivedClass* amrd_ptr =
      dynamic_cast<AMRDerivedClass*> (m_finer_level_ptr);
    if ( amrd_ptr != NULL)
    {
      amrd_ptr->m_coarse_average.averageToCoarse (m_state_new,
                                                  amrd_ptr->m_state_new);
    }
    else
    {
      MayDay::Error ("in AMRDerivedClass::postTimeStep: m_coarser_level_ptr is not castable to AMRDerivedClass*");
    }
  }
  if ( verbose ) pout () << "AMRDerivedClass::postTimeStep " << m_level << " finished" << endl;
}

// create tags

void
AMRDerivedClass::tagCells (IntVectSet& a_tags)
{
  if ( verbose ) pout () << "AMRDerivedClass::tagCells " << m_level << endl;

//  if (m_level > 1) return;
  int losize = int(4 * pow (Real(cos (m_time)), Real(2.)));
  int hisize = int(4 * (1.- pow (Real(cos (m_time)), Real(2.))));
  Box bhi (m_problem_domain.domainBox().bigEnd () - hisize * IntVect::Unit,
           m_problem_domain.domainBox().bigEnd () );
  a_tags |= bhi;
  Box blo (m_problem_domain.domainBox().smallEnd (),
           m_problem_domain.domainBox().smallEnd () + losize * IntVect::Unit );
  a_tags |= blo;
  IntVect ivc = (m_problem_domain.domainBox().smallEnd () + m_problem_domain.domainBox().bigEnd ()) / 2
    + 2 * (losize - 2) * IntVect::Unit;
  a_tags |= Box (ivc -  m_level    * IntVect::Unit,
                 ivc + (m_level+1) * IntVect::Unit);
}

// create tags at initialization

void
AMRDerivedClass::tagCellsInit (IntVectSet& a_tags)
{
  if ( verbose ) pout () << "AMRDerivedClass::tagCellsInit " << m_level << endl;

  IntVectSet local_tags;
//  if (m_level > 1) return;
  local_tags |= m_problem_domain.domainBox().bigEnd ();
  local_tags |= m_problem_domain.domainBox().smallEnd ();
  IntVect ivc = (m_problem_domain.domainBox().smallEnd () + m_problem_domain.domainBox().bigEnd ()) / 2;
  local_tags |= Box (ivc -  m_level    * IntVect::Unit,
                     ivc + (m_level+1) * IntVect::Unit);

  Vector<IntVectSet> all_tags;

  const int dest_proc = 0;
  gather(all_tags, local_tags, dest_proc);

  if (procID() == uniqueProc (SerialTask::compute) )
  {
      for (int i = 0; i < all_tags.size(); ++i)
      {
          a_tags |= all_tags[i];
      }
  }
}

DisjointBoxLayout
AMRDerivedClass::loadBalance(const Vector<Box>& a_grids)
{
// load balance and create boxlayout
  Vector<int> proc_map;
  if (procID () == uniqueProc (SerialTask::compute) )
  {
    LoadBalance (proc_map, a_grids);
  }
  broadcast (proc_map, uniqueProc (SerialTask::compute) );

  if ( verbose ) pout () << "AMRDerivedClass::loadBalance: procesor map: " << endl;
  for (int igrid = 0; igrid < a_grids.size (); ++igrid)
  {
    if ( verbose ) pout () << igrid << ": " << proc_map[igrid] << "  " << endl;
  }
  if ( verbose ) pout () << endl;
  return ( DisjointBoxLayout (a_grids, proc_map) );
}

// regrid

void
AMRDerivedClass::regrid (const Vector<Box>& a_new_grids)
{
  if ( verbose ) pout () << "AMRDerivedClass::regrid " << m_level << endl;

  m_level_grids = a_new_grids;

  const DisjointBoxLayout level_domain = loadBalance (a_new_grids);

  for (LayoutIterator lit = level_domain.layoutIterator (); lit.ok (); ++lit)
  {
    if ( verbose ) pout () << level_domain[lit ()] << endl;
  }

// save data for later copy
  LevelData<FArrayBox> old_state;
  old_state.define (m_state_new);

// reshape state with new grids
  m_state_new.define (level_domain, s_num_comps);
  m_state_old.define (level_domain, s_num_comps);

// maintain interlevel stuff
  m_coarse_average.define (level_domain,
                          s_num_comps,
                          m_ref_ratio);
  m_fine_interp.define (level_domain,
                       s_num_comps,
                       m_ref_ratio,
                       m_coarser_level_ptr->problemDomain());

// interpolate from coarser level
  if (m_coarser_level_ptr != NULL)
  {
    AMRDerivedClass* amrd_ptr =
      dynamic_cast<AMRDerivedClass*> (m_coarser_level_ptr);
    if (amrd_ptr != NULL)
    {
      m_fine_interp.interpToFine (m_state_new,
                                  amrd_ptr->m_state_new);
    }
    else
    {
      MayDay::Error ("in AMRDerivedClass::regrid: m_coarser_level_ptr is not castable to AMRDerivedClass*");
    }
  }
// copy from old state
  old_state.copyTo (old_state.interval (),
                    m_state_new,
                    m_state_new.interval () );
}

// set data after regridding
// (regridData should be removed)

void
AMRDerivedClass::regridData ()
{
  if ( verbose ) pout () << "AMRDerivedClass::regridData " << m_level << endl;

}

// initialize grid

void
AMRDerivedClass::initialGrid (const Vector<Box>& a_new_grids)
{
  if ( verbose ) pout () << "AMRDerivedClass::initialGrid " << m_level << endl;

  m_level_grids = a_new_grids;

  const DisjointBoxLayout level_domain = loadBalance (a_new_grids);

  for (LayoutIterator lit = level_domain.layoutIterator (); lit.ok (); ++lit)
  {
    if ( verbose ) pout () << level_domain[lit ()] << endl;
  }
  m_state_new.define (level_domain, s_num_comps);
  m_state_old.define (level_domain, s_num_comps);

  m_coarse_average.define (level_domain,
                          s_num_comps,
                          m_ref_ratio);
  m_fine_interp.define (level_domain,
                       s_num_comps,
                       m_ref_ratio,
                       m_problem_domain);
}

// initialize data

void
AMRDerivedClass::initialData ()
{

  if ( verbose ) pout () << "AMRDerivedClass::initialData " << m_level << endl;

  DataIterator dit = m_state_new.dataIterator();
  for (dit.begin (); dit.ok (); ++dit)
  {
    FArrayBox& state_fab = m_state_new[dit ()];
    Box b = state_fab.box ();
    BoxIterator bit (b);
    for (int comp = 0; comp < m_state_new.nComp (); ++comp)
    {
      for (bit.begin (); bit.ok (); ++bit)
      {
        IntVect iv = bit ();
/*
// quadratic data
        Real rsq = 0.0;
        for (int d = 0; d < SpaceDim; ++d)
        {
          Real x = m_dx * (iv[d]+0.5) - 0.5;
          rsq += 4.0*x*x;
        }
        state_fab (iv,comp) = rsq + comp;
*/
// linear data
        Real a = 0.0;
        for (int d = 0; d < SpaceDim; ++d)
        {
          Real x = m_dx * (iv[d]+0.5) - 0.5;
          a += x;
        }
        state_fab (iv,comp) = a + comp;
/**/
      }
    }
  }
}

// things to do after initialization

void
AMRDerivedClass::postInitialize ()
{
  if ( verbose ) pout () << "AMRDerivedClass::postInitialize " << m_level << endl;
}

// write checkpoint header

#ifdef CH_USE_HDF5
void
AMRDerivedClass::writeCheckpointHeader (HDF5Handle& a_handle) const
{
  if ( verbose ) pout () << "AMRDerivedClass::writeCheckpointHeader" << endl;

  HDF5HeaderData header;
  header.m_int ["num_components"] = s_num_comps;
  char comp_str[30];
  for (int comp = 0; comp < s_num_comps; ++comp)
  {
    sprintf (comp_str, "component_%d", comp);
    header.m_string[comp_str] = s_state_names[comp];
  }
  header.writeToFile(a_handle);

  if ( verbose ) pout () << header << endl;
}

void
AMRDerivedClass::writeCheckpointLevel (HDF5Handle& a_handle) const
{
  if ( verbose ) pout () << "AMRDerivedClass::writeCheckpointLevel" << endl;

//  char* level_str = new char[int (log10 (m_level+1))+1];
  char level_str[20];
  sprintf (level_str, "%d", m_level);
  const std::string label = std::string ("level_") + level_str;
//  delete[] level_str;

  a_handle.setGroup (label);

  HDF5HeaderData header;

  header.m_int  ["ref_ratio"]   = m_ref_ratio;
  header.m_real ["dx"]          = m_dx;
  header.m_real ["dt"]          = m_dt;
  header.m_real ["time"]        = m_time;
  header.m_box  ["prob_domain"] = m_problem_domain.domainBox();
  header.writeToFile(a_handle);

  if ( verbose ) pout () << header << endl;

  write (a_handle, m_state_new.boxLayout ());
  write (a_handle, m_state_new, "data");
}

void
AMRDerivedClass::readCheckpointHeader  (HDF5Handle& a_handle)
{
  if ( verbose ) pout () << "AMRDerivedClass::readCheckpointHeader" << endl;

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if ( verbose ) pout () << "hdf5 header data:" << endl;
  if ( verbose ) pout () << header << endl;

// read number of components
  if (header.m_int.find ("num_components") == header.m_int.end())
  {
    MayDay::Error ("AMRDerivedClass::readCheckpointHeader: checkpoint file does not have num_components");
  }
  int num_comps = header.m_int ["num_components"];
  if (num_comps != s_num_comps)
  {
    MayDay::Error ("AMRDerivedClass::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }

// read component names
  std::string state_name;
  char comp_str[60];
  for (int comp = 0; comp < s_num_comps; ++comp)
  {
    sprintf (comp_str, "component_%d", comp);
    if (header.m_string.find (comp_str) == header.m_string.end())
    {
      MayDay::Error ("AMRDerivedClass::readCheckpointHeader: checkpoint file does not have enough component names");
    }
    state_name = header.m_string [comp_str];
    if (state_name != s_state_names[comp])
    {
      MayDay::Error("AMRDerivedClass::readCheckpointHeader: state_name in checkpoint does not match solver");
    }
  }

}

void
AMRDerivedClass::readCheckpointLevel (HDF5Handle& a_handle)
{
  if ( verbose ) pout () << "AMRDerivedClass::readCheckpointLevel" << endl;

//  char* level_str = new char[int (log10 (m_level+1))+1];
  char level_str[20];
  sprintf (level_str, "%d", m_level);
  const std::string label = std::string ("level_") + level_str;
//  delete[] level_str;

  a_handle.setGroup (label);

  HDF5HeaderData header;
  header.readFromFile (a_handle);

  if ( verbose ) pout () << "hdf5 header data:" << endl;
  if ( verbose ) pout () << header << endl;

// read refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
  {
    MayDay::Error("AMRDerivedClass::readCheckpointLevel: file does not contain ref_ratio");
  }
  m_ref_ratio = header.m_int ["ref_ratio"];

  if ( verbose ) pout () << "read ref_ratio = " << m_ref_ratio << endl;

// read dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRDerivedClass::readCheckpointLevel: file does not contain dx");
  }
  m_dx = header.m_real ["dx"];

  if ( verbose ) pout () << "read dx = " << m_dx << endl;

// read dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRDerivedClass::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real ["dt"];

  if ( verbose ) pout () << "read dt = " << m_dt << endl;

// read time
  if (header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMRDerivedClass::readCheckpointLevel: file does not contain time");
  }
  m_time = header.m_real ["time"];

  if ( verbose ) pout () << "read time = " << m_time << endl;

// read problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRDerivedClass::readCheckpointLevel: file does not contain prob_domain");
  }
  Box domainBox = header.m_box ["prob_domain"];
  m_problem_domain = ProblemDomain(domainBox);
  // need to figure out how to get periodic info into this
 CH_assert (false);

  // read grids
  Vector<Box> grids;
  const int grid_status = read (a_handle, grids);
  if (grid_status != 0)
  {
    MayDay::Error("AMRDerivedClass::readCheckpointLevel: file does not contain a Vector<Box>");
  }
// create level domain

  const DisjointBoxLayout level_domain = loadBalance (grids);

  if ( verbose ) pout () << "read level domain: " << endl;
  LayoutIterator lit = level_domain.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
  {
    const Box& b = level_domain[lit()];
    if ( verbose ) pout () << lit().intCode() << ": " << b << endl;
    m_level_grids.push_back(b);
  }
  if ( verbose ) pout () << endl;

// maintain interlevel stuff
  m_coarse_average.define (level_domain,
                          s_num_comps,
                          m_ref_ratio);
  m_fine_interp.define (level_domain,
                        s_num_comps,
                        m_ref_ratio,
                        m_problem_domain);

// reshape state with new grids
  m_state_new.define (level_domain, s_num_comps);
  const int data_status = read<FArrayBox> (a_handle,
                                           m_state_new,
                                           "data",
                                           level_domain);
  if (data_status != 0)
  {
    MayDay::Error("AMRDerivedClass::readCheckpointLevel: file does not contain state data");
  }
  m_state_old.define (level_domain, s_num_comps);

}

void
AMRDerivedClass::writePlotLevel (HDF5Handle& a_handle) const
{
  if ( verbose ) pout () << "AMRDerivedClass::writePlotLevel" << endl;
}

void
AMRDerivedClass::writePlotHeader (HDF5Handle& a_handle) const
{
  if ( verbose ) pout () << "AMRDerivedClass::writePlotHeader" << endl;
}
#endif

// compute dt

Real
AMRDerivedClass::computeDt ()
{
  if ( verbose ) pout () << "AMRDerivedClass::computeDt " << m_level << endl;

  m_dt *= 0.99;

  return m_dt;
}

// compute dt with initial data

Real
AMRDerivedClass::computeInitialDt ()
{
  if ( verbose ) pout () << "AMRDerivedClass::computeInitialDt " << m_level << endl;

  m_dt = m_initial_dt_multiplier * pow ( Real(0.4), Real (m_level));

  return m_dt;
}

class AMRDerivedClassFactory : public AMRLevelFactory
{
public:
  AMRDerivedClassFactory();

  virtual AMRLevel* new_amrlevel() const;

  virtual
  ~AMRDerivedClassFactory();

};

AMRDerivedClassFactory::AMRDerivedClassFactory ()
{
}
AMRDerivedClassFactory::~AMRDerivedClassFactory ()
{
}

// "virtual constructor"
AMRLevel*
AMRDerivedClassFactory::new_amrlevel() const
{
  AMRDerivedClass*
    amrdc_ptr = new AMRDerivedClass ();
  return (static_cast <AMRLevel*> (amrdc_ptr));
}

/// Prototypes:
int
testAMR();

void
parseTestOptions(int argc ,char* argv[]) ;

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testAMR();

  if ( status == 0 )
    pout() << indent << pgmname << " passed." << endl ;
  else
    pout() << indent << pgmname << " failed with return code " << status << endl ;

#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

int
testAMR ()
{
  Box prob_domain (IntVect::Zero,
                   3 * IntVect::Unit );

  int max_level = 3;
  Vector<int> ref_ratioes (max_level+1, 2);
  Vector<int> regrid_intervals (max_level+1, 4);
  AMRDerivedClassFactory amrd_fact;

  AMR amr;
  amr.define(max_level, ref_ratioes, prob_domain, &amrd_fact);

  amr.checkpointInterval(2);

  amr.setupForNewAMRRun();

  amr.run(8., 8);

  amr.conclude();

  return 0 ;
}

///
// Parse the standard test options (-v -q) out of the command line.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if ( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              // argv[i] = "" ;
            }
        }
    }
  return ;
}
