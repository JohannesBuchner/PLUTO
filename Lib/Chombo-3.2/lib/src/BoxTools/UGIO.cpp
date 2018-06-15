#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
using std::fstream;
using std::string;

#include "CH_HDF5.H"
#include "UGIO.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5

void WriteUGHDF5(const string&               a_filename,
                 const DisjointBoxLayout&    a_grids,
                 const LevelData<FArrayBox>& a_data,
                 const Box&                  a_domain)
{
  HDF5Handle handle(a_filename.c_str(), HDF5Handle::CREATE);
  WriteUGHDF5(handle, a_grids, a_data, a_domain);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  handle.close();
}

void WriteUGHDF5(HDF5Handle&                 a_handle,
                 const DisjointBoxLayout&    a_grids,
                 const LevelData<FArrayBox>& a_data,
                 const Box&                  a_domain)
{
  HDF5HeaderData header;
  int nComp = a_data.nComp();

  int numlevels = 1;
  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       =  numlevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
  {
    char labelChSt[100];
    sprintf(labelChSt, "component_%d", ivar);
    string label(labelChSt);
    header.m_string[label] = label;
  }

  header.writeToFile(a_handle);

  Real dtLevel = 1.0;
  Real dxLevel = 1.0;
  Real time = 1.0;
  int refLevel = 1;
  int ilev = 0;

  int eek = writeLevel(a_handle, ilev, a_data, dxLevel, dtLevel, time,
                       a_domain, refLevel);

  if (eek != 0)
  {
    MayDay::Error("WriteUGHDF5: Error in writeLevel");
  }
}

int ReadUGHDF5(const string&         a_filename,
               DisjointBoxLayout&    a_grids,
               LevelData<FArrayBox>& a_data,
               Box&                  a_domain)
{
  HDF5Handle handle;
  int err = handle.open(a_filename.c_str(),  HDF5Handle::OPEN_RDONLY);

  if (err < 0)
  {
    return -4;
  }

  int eekflag = ReadUGHDF5(handle, a_grids, a_data, a_domain);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  handle.close();
  return eekflag;
}

int ReadUGHDF5(HDF5Handle&           a_handle,
               DisjointBoxLayout&    a_grids,
               LevelData<FArrayBox>& a_data,
               Box&                  a_domain)
{
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  int numLevels = header.m_int["num_levels"];
  if (numLevels != 1)
  {
    MayDay::Warning("Tried to read multilevel data using ReadUGHDF5");
    return -5;
  }

  int nComp = header.m_int["num_components"];
  if (nComp <= 0)
  {
    MayDay::Warning("ReadUGHDF5: Bogus number of Components");
    return -2;
  }

  int ilev = 0;
  Real dtLevel;
  Real dxLevel;
  Real time;
  int refLevel =1;

  int eek = readLevel(a_handle, ilev, a_data, dxLevel, dtLevel, time,
                      a_domain, refLevel, Interval(), true);

  if (eek != 0)
  {
    MayDay::Warning("ReadUGHDF5: readLevel failed");
    return -3;
  }

  a_grids = a_data.getBoxes();

  return 0;
}

#endif
#include "NamespaceFooter.H"
