#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Fri, Dec 3, 1999

#include <fstream>
#include <cstdio>
#include <string>
using std::fstream;
using std::string;

#include <cstdio>
#ifndef __IBMCPP__
using std::sprintf;
#endif
using std::tmpnam;

#include <cstdlib>
using std::system;

#include <cmath>

#ifdef CH_USE_HDF5
#include "CH_HDF5.H"
#endif

#include "AMRIO.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "VisItChomboDriver.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5
/*
\\ write out hierarchy of amr data in HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_vectRatio == refinement ratio at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output
*/
void
WriteAMRHierarchyHDF5(const string& filename,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels)
{
  CH_TIMERS("WriteAMRHierarchyHDF5");
  CH_TIMER("CreateFile",createFile);
  CH_TIMER("WriteFile",writeFile);
  CH_TIMER("CloseFile",closeFile);

  CH_START(createFile);
  HDF5Handle handle(filename.c_str(),  HDF5Handle::CREATE);
  CH_STOP(createFile);

  CH_START(writeFile);
  WriteAMRHierarchyHDF5(handle, a_vectGrids, a_vectData, a_vectNames,
                        a_domain, a_dx, a_dt, a_time, a_refRatio, a_numLevels);
  CH_STOP(writeFile);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  CH_START(closeFile);
  handle.close();
  CH_STOP(closeFile);
}

/*
\\ write out hierarchy of anisotropic amr data in HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing in each direction at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_vectRatio == refinement ratio in each direction at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output
*/
void
WriteAnisotropicAMRHierarchyHDF5(
    const string& filename,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_refRatios,
    const int& a_numLevels)
{
  CH_TIMERS("WriteAnisotropicAMRHierarchyHDF5");
  CH_TIMER("CreateFile",createFile);
  CH_TIMER("WriteFile",writeFile);
  CH_TIMER("CloseFile",closeFile);

  CH_START(createFile);
  HDF5Handle handle(filename.c_str(),  HDF5Handle::CREATE);
  CH_STOP(createFile);

  CH_START(writeFile);
  WriteAnisotropicAMRHierarchyHDF5(
      handle, a_vectGrids, a_vectData, a_vectNames,
      a_domain, a_dx, a_dt, a_time, a_refRatios, a_numLevels);
  CH_STOP(writeFile);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  CH_START(closeFile);
  handle.close();
  CH_STOP(closeFile);
}

void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels)
{
  CH_assert(a_numLevels > 0);
  CH_assert(a_vectData.size()  >= a_numLevels);
  CH_assert(a_refRatio.size() >= a_numLevels-1);

  HDF5HeaderData header;
  int nComp = a_vectNames.size();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = a_numLevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      header.m_string[label] = a_vectNames[ivar];
    }
  header.writeToFile(handle);

  Box domainLevel = a_domain;
  Real dtLevel = a_dt;
  Real dxLevel = a_dx;
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int refLevel = 1;
      if (ilev != a_numLevels -1)
        {
          refLevel = a_refRatio[ilev];
        }
      if (ilev != 0)
        {
          domainLevel.refine(a_refRatio[ilev-1]);
          dtLevel /= a_refRatio[ilev-1];
          dxLevel /= a_refRatio[ilev-1];
        }
      CH_assert(a_vectData[ilev] != NULL);
      const LevelData<FArrayBox>& dataLevel = *a_vectData[ilev];
      CH_assert(dataLevel.nComp() == nComp);
      Interval comps(0,nComp-1);
      IntVect ghostVect = a_vectData[0]->ghostVect();
      int eek = writeLevel(handle, ilev, dataLevel,
                           dxLevel, dtLevel, a_time,
                           domainLevel, refLevel, ghostVect, comps);
      if (eek != 0)
        {
          MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
        }
    }
}

void
WriteAnisotropicAMRHierarchyHDF5(
    HDF5Handle& handle,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx, // Grid spacing in each direction
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_refRatios, // for each level, in each direction
    const int& a_numLevels)
{
  CH_assert(a_numLevels > 0);
  CH_assert(a_vectData.size()  >= a_numLevels);
  CH_assert(a_refRatios.size() >= a_numLevels-1);

  HDF5HeaderData header;
  int nComp = a_vectNames.size();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = a_numLevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      header.m_string[label] = a_vectNames[ivar];
    }
  header.writeToFile(handle);

  Box domainLevel = a_domain;
  Real dtLevel = a_dt;
  RealVect dxLevel = a_dx;
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      IntVect refLevel = IntVect::Unit;
      if (ilev != a_numLevels -1)
        {
          refLevel = a_refRatios[ilev];
        }
      if (ilev != 0)
        {
          domainLevel.refine(a_refRatios[ilev-1]);
          dtLevel /= a_refRatios[ilev-1][0]; // HACK - just use 0 dir ref ratio
          dxLevel /= a_refRatios[ilev-1];
        }
      CH_assert(a_vectData[ilev] != NULL);
      const LevelData<FArrayBox>& dataLevel = *a_vectData[ilev];
      CH_assert(dataLevel.nComp() == nComp);
      Interval comps(0,nComp-1);
      IntVect ghostVect = a_vectData[0]->ghostVect();
      int eek = writeLevel(handle, ilev, dataLevel,
                           dxLevel, dtLevel, a_time,
                           domainLevel, refLevel, ghostVect, comps);
      if (eek != 0)
        {
          MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
        }
    }
}

void
WriteAMRHierarchyHDF5(const string& filename,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Box& a_domain,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels)
{

  HDF5Handle handle(filename.c_str(),  HDF5Handle::CREATE);
  WriteAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                        a_domain, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();
}

void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<FArrayBox>* > & a_vectData,
                      const Box& a_domain,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels)
{
  CH_assert(a_numLevels > 0);
  CH_assert(a_vectData.size()  >= a_numLevels);
  CH_assert(a_refRatio.size() >= a_numLevels-1);
  Real dxin = 1.0;
  Real dtin = 1.0;
  Real time = 1.0;

  HDF5HeaderData header;
  int nComp = a_vectData[0]->nComp();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = a_numLevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      header.m_string[label] = label;
    }
  header.writeToFile(handle);

  Box domainLevel = a_domain;
  Real dtLevel = dtin;
  Real dxLevel = dxin;
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int refLevel = 1;
      if (ilev != a_numLevels -1)
        {
          refLevel = a_refRatio[ilev];
        }
      if (ilev != 0)
        {
          domainLevel.refine(a_refRatio[ilev-1]);
          dtLevel /= a_refRatio[ilev-1];
          dxLevel /= a_refRatio[ilev-1];
        }
      CH_assert(a_vectData[ilev] != NULL);
      const LevelData<FArrayBox>& dataLevel = *a_vectData[ilev];
      CH_assert(dataLevel.nComp() == nComp);
      Interval comps(0,nComp-1);
      IntVect ghostVect = a_vectData[0]->ghostVect();
      int eek = writeLevel(handle, ilev, dataLevel,
                           dxLevel, dtLevel, time,
                           domainLevel, refLevel, ghostVect, comps);
      if (eek != 0)
        {
          MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
        }
    }
}

//
/*
\\ Read in hierarchy of amr data in HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_vectRatio == refinement ratio at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output

return values:
0: success
-1: bogus number of levels
-2: bogus number of components
-3: error in readlevel
-4: file open failed
*/
int
ReadAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<int>& a_refRatio,
                     int& a_numLevels)
{
  HDF5Handle handle;
  int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
  if ( err < 0)
    {
      return -4;
    }
  int eekflag = ReadAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                                     a_vectNames, a_domain, a_dx, a_dt,
                                     a_time, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();

  return (eekflag);
}

int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<int>& a_refRatio,
                     int& a_numLevels)
{

  HDF5HeaderData header;
  header.readFromFile(handle);

  a_numLevels = header.m_int["num_levels"];
  if (a_numLevels <= 0)
  {
    MayDay::Warning("ReadAMRHierarchyHDF5: Bogus number of levels");
    return (-1);
  }
  a_vectData.resize(a_numLevels);
  a_refRatio.resize(a_numLevels);
  a_vectGrids.resize(a_numLevels);

  int nComp = header.m_int["num_components"];
  if (nComp <= 0)
  {
    MayDay::Warning("ReadAMRHierarchyHDF5: Bogus number of Components");
    return (-2);
  }
  a_vectNames.resize(nComp);

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      a_vectNames[ivar] = header.m_string[label];
    }
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int refLevel = 0;
      Box domainLevel;
      Real dtLevel;
      Real dxLevel;
      a_vectData[ilev] = new LevelData<FArrayBox>();
      int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                          dxLevel, dtLevel,  a_time,
                          domainLevel, refLevel, Interval(), true);
      if (eek != 0)
      {
        MayDay::Warning("ReadAMRHierarchyHDF5: readLevel failed");
        return (-3);
      }

      const DisjointBoxLayout& dbl = a_vectData[ilev]->getBoxes();
      a_vectGrids[ilev]= dbl;

      if (ilev == 0)
        {
          a_domain = domainLevel;
          a_dt = dtLevel;
          a_dx = dxLevel;
        }
      a_refRatio[ilev] = refLevel;
    }
  return (0);
}

int
ReadAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<int>& a_refRatio,
                     int& a_numLevels)
{
  HDF5Handle handle;
  int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
  if ( err < 0)
  {
    return -4;
  }

  int eekflag = ReadAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                                     a_domain, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();
  return (eekflag);
}

int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<int>& a_refRatio,
                     int& a_numLevels)
{
  HDF5HeaderData header;
  header.readFromFile(handle);

  a_numLevels = header.m_int["num_levels"];
  if (a_numLevels <= 0)
  {
    MayDay::Warning("ReadAMRHierarchyHDF5: Bogus number of levels");
    return (-1);
  }
  a_vectData.resize(a_numLevels);
  a_refRatio.resize(a_numLevels);
  a_vectGrids.resize(a_numLevels);

  //  int nComp = header.m_int["num_components"];
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      int refLevel = 0;
      Box domainLevel;
      Real dtLevel;
      Real dxLevel;
      Real time;
      a_vectData[ilev] = new LevelData<FArrayBox>();
      int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                          dxLevel, dtLevel,  time,
                          domainLevel, refLevel, Interval(), true);
      if (eek != 0)
      {
        MayDay::Warning("ReadAMRHierarchyHDF5: readLevel failed");
        return (-3);
      }

      const DisjointBoxLayout& dbl = a_vectData[ilev]->getBoxes();
      a_vectGrids[ilev]= dbl;

      if (ilev == 0)
        {
          a_domain = domainLevel;
        }
      a_refRatio[ilev] = refLevel;
    }

  return (0);
}

//
/*
\\ Read in hierarchy of amr data in ANISOTROPIC HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_vectRatio == refinement ratio at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output

return values:
0: success
-1: bogus number of levels
-2: bogus number of components
-3: error in readlevel
-4: file open failed
*/
int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     RealVect& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels)
{
  HDF5Handle handle;
  int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
  if ( err < 0)
    {
      return -4;
    }
  int eekflag = ReadAnisotropicAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                                     a_vectNames, a_domain, a_dx, a_dt,
                                     a_time, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();

  return (eekflag);
}

int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     RealVect& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels)
{

  HDF5HeaderData header;
  header.readFromFile(handle);

  a_numLevels = header.m_int["num_levels"];
  if (a_numLevels <= 0)
  {
    MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: Bogus number of levels");
    return (-1);
  }
  a_vectData.resize(a_numLevels);
  a_refRatio.resize(a_numLevels);
  a_vectGrids.resize(a_numLevels);

  int nComp = header.m_int["num_components"];
  if (nComp <= 0)
  {
    MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: Bogus number of Components");
    return (-2);
  }
  a_vectNames.resize(nComp);

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      a_vectNames[ivar] = header.m_string[label];
    }
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      IntVect refLevel = IntVect::Zero;
      Box domainLevel;
      Real dtLevel;
      RealVect dxLevel;
      a_vectData[ilev] = new LevelData<FArrayBox>();
      int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                          dxLevel, dtLevel,  a_time,
                          domainLevel, refLevel, Interval(), true);
      if (eek != 0)
      {
        MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: readLevel failed");
        return (-3);
      }

      const DisjointBoxLayout& dbl = a_vectData[ilev]->getBoxes();
      a_vectGrids[ilev]= dbl;

      if (ilev == 0)
        {
          a_domain = domainLevel;
          a_dt = dtLevel;
          a_dx = dxLevel;
        }
      a_refRatio[ilev] = refLevel;
    }
  return (0);
}

int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels)
{
  HDF5Handle handle;
  int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
  if ( err < 0)
  {
    return -4;
  }

  int eekflag = ReadAnisotropicAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                                     a_domain, a_refRatio, a_numLevels);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();
  return (eekflag);
}

int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels)
{
  HDF5HeaderData header;
  header.readFromFile(handle);

  a_numLevels = header.m_int["num_levels"];
  if (a_numLevels <= 0)
  {
    MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: Bogus number of levels");
    return (-1);
  }
  a_vectData.resize(a_numLevels);
  a_refRatio.resize(a_numLevels);
  a_vectGrids.resize(a_numLevels);

  //  int nComp = header.m_int["num_components"];
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      IntVect refLevel = IntVect::Zero;
      Box domainLevel;
      Real dtLevel;
      RealVect dxLevel;
      Real time;
      a_vectData[ilev] = new LevelData<FArrayBox>();
      int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                          dxLevel, dtLevel,  time,
                          domainLevel, refLevel, Interval(), true);
      if (eek != 0)
      {
        MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: readLevel failed");
        return (-3);
      }

      const DisjointBoxLayout& dbl = a_vectData[ilev]->getBoxes();
      a_vectGrids[ilev]= dbl;

      if (ilev == 0)
        {
          a_domain = domainLevel;
        }
      a_refRatio[ilev] = refLevel;
    }

  return (0);
}

// put simple debugging functions here, at the end of the hdf5 stuff

// put simple debugging functions here, at the end of the hdf5 stuff
static void
ChomboVisVisualizeFile(const char *fname)
{
  char command[2000];
  sprintf(command,"$CHOMBOVIS_HOME/bin/chombovis debug_level=0 %s &",fname);
  int ret;
  ret = system(command);
}

static void
ChomboBrowserBrowseFile(const char *fname)
{
  char command[2000];
  sprintf(command,"$CHOMBOVIS_HOME/bin/chombobrowser debug_level=0 %s &",fname);
  int ret;
  ret = system(command);
}

static VisItChomboDriver visit;
static void
VisItVisualizeFile(const char *fname)
{
  visit.VisualizeFile(fname);
}

static void
VisItBrowseFile(const char *fname)
{
  visit.BrowseFile(fname);
}

static void
VisualizeFile(const char *fname)
{
  const char *use_chombovis = getenv("CHOMBO_USE_CHOMBOVIS");
  if (use_chombovis)
      ChomboVisVisualizeFile(fname);
  else
      VisItVisualizeFile(fname);
}

static void
BrowseFile(const char *fname)
{
  const char *use_chombovis = getenv("CHOMBO_USE_CHOMBOVIS");
  if (use_chombovis)
      ChomboBrowserBrowseFile(fname);
  else
      VisItBrowseFile(fname);
}

void
writeFAB(const FArrayBox* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  char fname[12];
  sprintf(fname,"fab.%dd.hdf5", SpaceDim);
  writeFABname(a_dataPtr, fname);
}

void
viewFAB(const FArrayBox* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tmpnam(NULL);
  writeFABname(a_dataPtr, fname);
  VisualizeFile(fname);
}

void
viewBFI(const BaseFab<int>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  FArrayBox fab(a_dataPtr->box(), a_dataPtr->nComp());
  BoxIterator bit(fab.box());
  for (bit.begin(); bit.ok(); bit.next())
    {
      const IntVect& iv = bit();
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
        {
          fab(iv, ivar) = a_dataPtr->operator()(iv, ivar);
        }
    }

  const char* fname = tmpnam(NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
viewBFIV(const BaseFab<IntVect>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  FArrayBox fab(a_dataPtr->box(), SpaceDim * a_dataPtr->nComp());
  BoxIterator bit(fab.box());
  for (bit.begin(); bit.ok(); bit.next())
    {
      const IntVect& iv = bit();
      for (int ivar = 0; ivar < a_dataPtr->nComp(); ivar++)
        {
          const IntVect& ivVal = a_dataPtr->operator()(iv, ivar);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fab(iv, SpaceDim*ivar + idir) = ivVal[idir];
            }
        }
    }

  const char* fname = tmpnam(NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
viewBFRV(const BaseFab<RealVect>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  FArrayBox fab(a_dataPtr->box(), SpaceDim * a_dataPtr->nComp());
  BoxIterator bit(fab.box());
  for (bit.begin(); bit.ok(); bit.next())
    {
      const IntVect& iv = bit();
      for (int ivar = 0; ivar < a_dataPtr->nComp(); ivar++)
        {
          const RealVect& rvVal = a_dataPtr->operator()(iv, ivar);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fab(iv, SpaceDim*ivar + idir) = rvVal[idir];
            }
        }
    }

  const char* fname = tmpnam(NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
viewIVSFAB(const IVSFAB<Real>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const IntVectSet& ivs = a_dataPtr->getIVS();
  const Box& bx = ivs.minBox();
  int ncomp = a_dataPtr->nComp();

  FArrayBox fab(bx, ncomp);
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      for (int ivar = 0; ivar < ncomp; ivar++)
        {
          fab(iv, ivar) = a_dataPtr->operator()(iv, ivar);
        }
    }

  const char* fname = tmpnam(NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
viewIVSFABI(const IVSFAB<int>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const IntVectSet& ivs = a_dataPtr->getIVS();
  const Box& bx = ivs.minBox();
  int ncomp = a_dataPtr->nComp();

  FArrayBox fab(bx, ncomp);
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      for (int ivar = 0; ivar < ncomp; ivar++)
        {
          fab(iv, ivar) = a_dataPtr->operator()(iv, ivar);
        }
    }

  const char* fname = tmpnam(NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
viewIVSFABIV(const IVSFAB<IntVect>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const IntVectSet& ivs = a_dataPtr->getIVS();
  const Box& bx = ivs.minBox();
  int ncomp = a_dataPtr->nComp();

  FArrayBox fab(bx, SpaceDim * ncomp);
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      for (int ivar = 0; ivar < a_dataPtr->nComp(); ivar++)
        {
          const IntVect& ivVal = a_dataPtr->operator()(iv, ivar);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fab(iv, SpaceDim*ivar + idir) = ivVal[idir];
            }
        }
    }

  const char* fname = tmpnam(NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
viewIVSFABRV(const IVSFAB<RealVect>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const IntVectSet& ivs = a_dataPtr->getIVS();
  const Box& bx = ivs.minBox();
  int ncomp = a_dataPtr->nComp();

  FArrayBox fab(bx, SpaceDim * ncomp);
  for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      for (int ivar = 0; ivar < a_dataPtr->nComp(); ivar++)
        {
          const RealVect& rvVal = a_dataPtr->operator()(iv, ivar);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              fab(iv, SpaceDim*ivar + idir) = rvVal[idir];
            }
        }
    }

  const char* fname = tmpnam(NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);
}

void
browseFAB(const FArrayBox* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tmpnam(NULL);
  writeFABname(a_dataPtr, fname);
  BrowseFile(fname);

}

void
writeBFR(const BaseFab<Real>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  //copy the BaseFab over to a FArrayBox
  Box fabBox = (*a_dataPtr).box();
  int  ncomp = (*a_dataPtr).nComp();
  FArrayBox fab(fabBox,ncomp);
  fab.copy(*a_dataPtr);

  char fname[12];
  sprintf(fname,"fab.%dd.hdf5", SpaceDim);
  writeFABname(&fab, fname);
}

void
viewBFR(const BaseFab<Real>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  //copy the BaseFab over to a FArrayBox
  Box fabBox = (*a_dataPtr).box();
  int  ncomp = (*a_dataPtr).nComp();
  FArrayBox fab(fabBox,ncomp);
  fab.copy(*a_dataPtr);

  const char* fname = tmpnam(NULL);
  writeFABname(&fab, fname);
  VisualizeFile(fname);

}

void
writeFABname(const FArrayBox      * a_dataPtr,
             const char           * a_filename,
             const Vector<string> & a_compNames)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const FArrayBox& data = *a_dataPtr;

  HDF5Handle handle(a_filename, HDF5Handle::CREATE_SERIAL);
  HDF5HeaderData header;

  int numlevels= 1;
  int nComp = data.nComp();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = numlevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      if ( a_compNames.size() > ivar )
        {
          header.m_string[label] = a_compNames[ivar] ;
        }
      else
        {
          header.m_string[label] = label;
        }
    }
  header.writeToFile(handle);

  Box domainLevel = data.box();
  // put bogus numbers here
  Real dtLevel = 1.0;
  Real dxLevel = 1.0;
  Real time = 1.0;

  int refLevel = 1;

  // build bogus DisjointBoxLayout here
  Vector<Box> boxes(1,domainLevel);
  unsigned int myprocID= procID();
  Vector<int> procAssign(1,myprocID);
  DisjointBoxLayout grids(boxes, procAssign);
  LevelData<FArrayBox> ldf(grids, nComp);
  // now copy fab into ldf
  DataIterator dit = ldf.dataIterator();
  ldf[dit()].copy(data);

  int eek = writeLevel(handle, 0, ldf, dxLevel, dtLevel, time,
                       domainLevel, refLevel);
  if (eek != 0)
  {
    MayDay::Error("writeFABname: error in writeLEvel");
  }

  handle.close();
}

void
viewLevelNoFine(const LevelData<FArrayBox>* a_dataPtr,
                const LevelData<FArrayBox>* a_dataFinePtr,
                int a_refRatio)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  if (a_dataFinePtr == NULL)
  {
    return;
  }

  int nVar = a_dataPtr->nComp();
  const DisjointBoxLayout& layoutOrig = a_dataPtr->disjointBoxLayout();
  const DisjointBoxLayout& layoutFine = a_dataFinePtr->disjointBoxLayout();
  DisjointBoxLayout layoutFineCoarsened;
  coarsen(layoutFineCoarsened, layoutFine, a_refRatio);
  LevelData<FArrayBox>* dataFineCoarsenedPtr =
    new LevelData<FArrayBox>(layoutFineCoarsened, nVar);
  for (DataIterator dit = a_dataFinePtr->dataIterator(); dit.ok(); ++dit)
    {
      dataFineCoarsenedPtr->operator[](dit).setVal(0.);
    }
  LevelData<FArrayBox>* dataNoFinePtr =
    new LevelData<FArrayBox>(layoutOrig, nVar);
  // *dataNoFinePtr := *a_dataPtr
  a_dataPtr->copyTo(*dataNoFinePtr);
  // Zero out *dataNoFinePtr on coarsened grids of *a_dataFinePtr.
  dataFineCoarsenedPtr->copyTo(*dataNoFinePtr);
  delete dataFineCoarsenedPtr;

  viewLevel(dataNoFinePtr);
  delete dataNoFinePtr;
}

void
writeLevel(const LevelData<FArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  char fname[12];
  sprintf(fname, "LDF.%dd.hdf5", SpaceDim);
  writeLevelname(a_dataPtr, fname);
}

void
viewLevel(const LevelData<FArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tmpnam(NULL);
  writeLevelname(a_dataPtr, fname);
  VisualizeFile(fname);

}

void
browseLevel(const LevelData<FArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tmpnam(NULL);
  writeLevelname(a_dataPtr, fname);
  BrowseFile(fname);

}

void
viewLevelNoGhost(const LevelData<FArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const DisjointBoxLayout& grids = a_dataPtr->getBoxes();
  int nVar = a_dataPtr->nComp();

  LevelData<FArrayBox> temp(grids, nVar);
  a_dataPtr->copyTo(a_dataPtr->interval(), temp, temp.interval());

  const char* fname = tmpnam(NULL);
  writeLevelname(&temp, fname);

  VisualizeFile(fname);

}

void
viewVectorLevel(const Vector<LevelData<FArrayBox>*>* a_dataPtr,
                const Vector<int>*                   a_refRatios)
{
  const char* fname = tmpnam(NULL);
  writeVectorLevelName(a_dataPtr, a_refRatios, fname);

  VisualizeFile(fname);
}

void
writeLevelname(const LevelData<FArrayBox>* a_dataPtr,
               const char*                 a_filename)
{
  Vector<LevelData<FArrayBox>*> data(1);
  data[0]=(LevelData<FArrayBox>*)a_dataPtr;
  Vector<int> refRatios(1,1);
  writeVectorLevelName(&data, &refRatios, a_filename);
}

void
writeVectorLevelName(const Vector<LevelData<FArrayBox>*>* a_dataPtr,
                     const Vector<int>*          a_refRatios,
                     const char*                 a_filename)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  HDF5Handle handle(a_filename, HDF5Handle::CREATE_SERIAL);
  HDF5HeaderData header;

  int numlevels = a_dataPtr->size();
  int nComp =(a_dataPtr->operator[](0))->nComp();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int ["num_levels"]       = numlevels;
  header.m_int ["num_components"]    = nComp;

  for (int ivar = 0; ivar < nComp; ivar++)
    {
      char labelChSt[100];
      sprintf(labelChSt, "component_%d", ivar);
      string label(labelChSt);
      header.m_string[label] = label;
    }
  header.writeToFile(handle);

  // put bogus numbers here
  Real dtLevel = 1.0;
  Real dxLevel = 1.0;
  Real time = 1.0;

  for (int level=0; level<a_dataPtr->size(); level++)
    {
      const LevelData<FArrayBox>& data = *(a_dataPtr->operator[](level));

      // need to figure out what domain will contain this LevelData
      // This must be LayoutIterator instead of DataIterator because
      // we need domain over boxes in ALL procs.
      const DisjointBoxLayout& levelBoxes = data.getBoxes();
      LayoutIterator lit = levelBoxes.layoutIterator();
      lit.reset();
      Box domain;
      // check to see if DisjointBoxLayout contains a valid domain;
      // if so, use it. If not, compute the smallest domain box which
      // contains this DBL.
      const ProblemDomain& pdomain = levelBoxes.physDomain();
      domain = pdomain.domainBox();
      if (domain.isEmpty())
        {
          domain = levelBoxes.get(lit());
          for (lit.reset(); lit.ok(); ++lit)
            {
              const Box thisBox = levelBoxes.get(lit());
              D_TERM6(
                      if (thisBox.smallEnd(0)<domain.smallEnd(0))
                        domain.setSmall(0,thisBox.smallEnd(0)); ,
                      if (thisBox.smallEnd(1)<domain.smallEnd(1))
                        domain.setSmall(1,thisBox.smallEnd(1)); ,
                      if (thisBox.smallEnd(2)<domain.smallEnd(2))
                        domain.setSmall(2, thisBox.smallEnd(2)); ,
                      if (thisBox.smallEnd(3)<domain.smallEnd(3))
                        domain.setSmall(3,thisBox.smallEnd(3)); ,
                      if (thisBox.smallEnd(4)<domain.smallEnd(4))
                        domain.setSmall(4,thisBox.smallEnd(4)); ,
                      if (thisBox.smallEnd(5)<domain.smallEnd(5))
                        domain.setSmall(5, thisBox.smallEnd(5)); );

              D_TERM6(
                      if (thisBox.bigEnd(0)>domain.bigEnd(0))
                        domain.setBig(0,thisBox.bigEnd(0)); ,
                      if (thisBox.bigEnd(1)>domain.bigEnd(1))
                        domain.setBig(1,thisBox.bigEnd(1)); ,
                      if (thisBox.bigEnd(2)>domain.bigEnd(2))
                        domain.setBig(2, thisBox.bigEnd(2)); ,
                      if (thisBox.bigEnd(3)>domain.bigEnd(3))
                        domain.setBig(3,thisBox.bigEnd(3)); ,
                      if (thisBox.bigEnd(4)>domain.bigEnd(4))
                        domain.setBig(4,thisBox.bigEnd(4)); ,
                      if (thisBox.bigEnd(5)>domain.bigEnd(5))
                        domain.setBig(5, thisBox.bigEnd(5)); );

            } // end loop over boxes on level to determine "domain"
        } //end if ProblemDomain from DisjointBoxLayout wasn't defined

      int refLevel = 1;
      if (level < a_dataPtr->size()-1)
        {
          refLevel = a_refRatios->operator[](level);
        }

      Interval comps(0,nComp-1);
      IntVect ghostVect = data.ghostVect();
      int eek = writeLevel(handle, level, data, dxLevel, dtLevel, time,
                           domain, refLevel, ghostVect, comps);

      dtLevel /= refLevel;
      dxLevel /= refLevel;

      if (eek != 0)
      {
        MayDay::Error("writeLDFname: error in writeLEvel");
      }
    }
  handle.close();
}

void
writeDBL(const DisjointBoxLayout* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  char fname[12];
  sprintf(fname, "DBL.%dd.hdf5", SpaceDim);
  writeDBLname(a_dataPtr, fname);
}

void
viewDBL(const DisjointBoxLayout* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const char* fname = tmpnam(NULL);
  writeDBLname(a_dataPtr, fname);

  VisualizeFile(fname);

}

void
viewLevelBFI(const LevelData<BaseFab<int> >* a_dataPtr)
{
  if (a_dataPtr == NULL)
    {
      return;
    }

  LevelData<FArrayBox> fabData(a_dataPtr->getBoxes(),
                               a_dataPtr->nComp(),
                               a_dataPtr->ghostVect());

  DataIterator dit = fabData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisFabData = fabData[dit];
      const BaseFab<int>& intData = (*a_dataPtr)[dit];
      BoxIterator bit(intData.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          for (int ivar =0; ivar<intData.nComp(); ivar++)
            {
              thisFabData(iv,ivar) = intData(iv,ivar);
            }
        }
    }

  viewLevel(&fabData);

}

void writeCopier(const Copier* a_copier)
{
  Vector<Vector<Box> > boxes;
  Vector<Vector<int> > procs;
  Vector<Box> level;
  Vector<int> p;
  Vector<int> refRatio;
  Box domain;
  int last = -1;
  for (CopyIterator it(*a_copier, CopyIterator::LOCAL); it.ok(); ++it)
    {
      const MotionItem& item = it();
      const DataIndex& to = item.toIndex;
      if (last == -1)
        {
          last = to.intCode();
        }
      if (last != to.intCode())
        {
          if (level.size()==0)
            {
              break;
            }
          boxes.push_back(level);
          procs.push_back(p);
          level.resize(0);
          p.resize(0);
          refRatio.push_back(1);
          last = to.intCode();
        }
      p.push_back(0);
      level.push_back(item.toRegion);
      domain.minBox(item.toRegion);
    }

  int numLevels = boxes.size();
  Vector<DisjointBoxLayout> layouts(numLevels);
  Vector<LevelData<FArrayBox>*> ldf(numLevels);
  IntVect ghostVect(IntVect::Zero);
  for (int i=0; i<numLevels; ++i)
    {
      layouts[i].define(boxes[i], procs[i], ProblemDomain(domain));
      ldf[i] = new LevelData<FArrayBox>(layouts[i], 1, ghostVect);
    }
  WriteAMRHierarchyHDF5("copier.hdf5", layouts, ldf, domain, refRatio, numLevels);
  for (int i=0; i<numLevels; ++i)
    {
      delete ldf[i];
    }
}

void viewCopier(const Copier* a_copier)
{
  if (a_copier == NULL)
  {
    return;
  }

  writeCopier(a_copier);

  const char* fname = "copier.hdf5";
  VisualizeFile(fname);
}

void viewIVS(const IntVectSet* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  Vector<Box> boxes = a_dataPtr->boxes();
  Vector<int> procs(boxes.size(), 0);
  DisjointBoxLayout dbl(boxes, procs);
  viewDBL(&dbl);
}

void viewVectorBox(const Vector<Box>* boxes)
{
  if (boxes == NULL)
  {
    return;
  }

  Vector<int> procs(boxes->size(), 0);
  DisjointBoxLayout dbl(*boxes, procs);
  viewDBL(&dbl);
}

void
writeDBLname(const DisjointBoxLayout* a_dataPtr,
             const char*                 a_filename)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  const DisjointBoxLayout& grids = *a_dataPtr;

  // now create a LevelData<FArrayBox> based on the grids
  IntVect ghostVect(IntVect::Zero);
  LevelData<FArrayBox> data(grids, 1, ghostVect);

  // initialize data to procID.
  Real dataVal = procID();
  DataIterator dit = data.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      data[dit()].setVal(dataVal);
    }

  writeLevelname(&data, a_filename);
}

void
WritePartialAMRHierarchyHDF5(const string& filename,
                             const Vector<DisjointBoxLayout>& a_vectGrids,
                             const Vector<LevelData<FArrayBox>* > & a_vectData,
                             const Vector<string>& a_vectNames,
                             const Box& a_baseDomain,
                             const Real& a_baseDx,
                             const Real& a_dt,
                             const Real& a_time,
                             const Vector<int>& a_vectRatio,
                             const Interval& a_levels)
{
  int numLevels = a_levels.size();

  // now make new dataholders which only have numLevels entries,
  // and which will move the baseLevel to level 0

  Vector<DisjointBoxLayout> newVectGrids(numLevels);
  Vector<LevelData<FArrayBox> * > newVectData(numLevels);
  Vector<int> newVectRatio(numLevels);

  int leveloffset = a_levels.begin();
  for (int srcLevel = a_levels.begin(); srcLevel <= a_levels.end(); srcLevel++)
    {
      int destLevel = srcLevel-leveloffset;
      newVectGrids[destLevel] = a_vectGrids[srcLevel];
      newVectData[destLevel] = a_vectData[srcLevel];
      newVectRatio[destLevel] = a_vectRatio[srcLevel];
    }

  WriteAMRHierarchyHDF5(filename, newVectGrids, newVectData, a_vectNames,
                        a_baseDomain, a_baseDx, a_dt, a_time,
                        newVectRatio, numLevels);
}

void
viewFluxLevel(const LevelData<FluxBox>* a_dataPtr)
{
  if (a_dataPtr == NULL)
  {
    return;
  }

  // Create a LevelData that will store the flux components.
  int numComp = a_dataPtr->nComp();
  LevelData<FArrayBox>* temp = new LevelData<FArrayBox>(a_dataPtr->disjointBoxLayout(),
                                                        SpaceDim*numComp,
                                                        IntVect::Unit);

  // Copy the data.
  for (DataIterator dit = a_dataPtr->dataIterator(); dit.ok(); ++dit)
  {
    const FluxBox& F = (*a_dataPtr)[dit()];
    FArrayBox& f = (*temp)[dit()];
    for (int d = 0; d < SpaceDim; ++d)
    {
      f.shiftHalf(d, 1);
      const FArrayBox& Fd = F[d];
      Box box = Fd.box();
      box &= f.box();
      f.copy(Fd, box, 0, box, numComp*d, numComp);
      f.shiftHalf(d, -1);
    }
  }

  // Run it through the viewer.
  viewLevel(temp);
  delete temp;
}

#endif // CH_USE_HDF5

#include "NamespaceFooter.H"
