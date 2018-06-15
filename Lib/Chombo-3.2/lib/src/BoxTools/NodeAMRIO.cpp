#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeAMRIO.cpp
// petermc, 4 March 2003
// adapted from AMRIO by DTGraves, Fri, Dec 3, 1999

#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "NodeFArrayBox.H"
#include "HDF5Portable.H"
#include "NamespaceHeader.H"

using std::system;
using std::fstream;
using std::string;
#ifndef __IBMCPP__
using std::sprintf;
#endif

template <>
inline void dataSize(const NodeFArrayBox& item, Vector<int>& a_sizes,
                     const Box& box, const Interval& comps)
{
  Box boxNodes = surroundingNodes(box);
  a_sizes[0] = boxNodes.numPts() * comps.size();
}

template <>
inline const char* name(const NodeFArrayBox& a_dummySpecializationArg)
{
  // Attempt to get rid of warnings on IBM...
  //static const char* name = "NodeFArrayBox";
   const char* name = "NodeFArrayBox";
  return name;
}

#include "NamespaceFooter.H"

#ifdef CH_USE_HDF5

#include "CH_HDF5.H"
#include "NodeAMRIO.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "VisItChomboDriver.H"
#include "NamespaceHeader.H"

template <>
inline void dataTypes(Vector<hid_t>& a_types, const NodeFArrayBox& dummy)
{
  a_types.resize(1);
  a_types[0] = H5T_NATIVE_REAL;
}

// ---------------------------------------------------------
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
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels,
                      const RealVect& a_origin,
                      const Interval& a_comps )
{
  HDF5Handle handle(filename.c_str(),  HDF5Handle::CREATE);

  WriteAMRHierarchyHDF5(handle, a_vectGrids, a_vectData, a_vectNames,
                        a_domain, a_dx, a_dt, a_time, a_refRatio, a_numLevels,
                        a_origin, a_comps);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  handle.close();
}

// ---------------------------------------------------------
/*
\\ write out hierarchy of amr data in HDF5 format
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
    const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_refRatios,
    const int& a_numLevels,
    const RealVect& a_origin,
    const Interval& a_comps )
{
  HDF5Handle handle(filename.c_str(),  HDF5Handle::CREATE);

  WriteAnisotropicAMRHierarchyHDF5(handle, a_vectGrids, a_vectData, a_vectNames,
                        a_domain, a_dx, a_dt, a_time, a_refRatios, a_numLevels,
                        a_origin, a_comps);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  handle.close();
}

// ---------------------------------------------------------
void
WriteAMRHierarchyHDF5(const string& filename,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<Vector<LevelData<NodeFArrayBox>* > >& a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_vectRatio,
                      const int& a_numLevels)
{
  // We have a_vectData[component][level].
  // Each component should have the same number of levels.
  CH_assert(a_numLevels > 0);
  CH_assert(a_vectRatio.size() >= a_numLevels-1);
  int numvarsOrig = a_vectData.size();
  int numvars = 0;
  Vector<int> numvarsEach(numvarsOrig);
  for (int ivarOrig = 0; ivarOrig < numvarsOrig; ivarOrig++)
    {
      CH_assert( a_vectData[ivarOrig].size() >= a_numLevels );
      numvarsEach[ivarOrig] = a_vectData[ivarOrig][0]->nComp();
      numvars += numvarsEach[ivarOrig];
    }

  CH_assert( a_vectNames.size() == numvars );

  Vector<LevelData<NodeFArrayBox>* > vectVarsNode(a_numLevels, NULL);
  for (int ilev = 0; ilev < a_numLevels; ilev++)
    {
      const DisjointBoxLayout& grids = a_vectGrids[ilev];
      // What about ghosts in the LevelData?
      // Different variables may have different ghost vectors.
      // Solution:  dispense with ghosts.
      vectVarsNode[ilev] = new LevelData<NodeFArrayBox>(grids, numvars);
      LevelData<NodeFArrayBox>& vectVarsNodeLevel = *vectVarsNode[ilev];
      int ivar = 0;
      for (int ivarOrig = 0; ivarOrig < numvarsOrig; ivarOrig++)
        for (int ivarEach = 0; ivarEach < numvarsEach[ivarOrig]; ivarEach++)
          {
            Interval newInterval(ivar, ivar);
            Interval oldInterval(ivarEach, ivarEach);

            const LevelData<NodeFArrayBox>& dataLevel = *a_vectData[ivarOrig][ilev];
            // dataLevel.copyTo(oldInterval, vectVarsNodeLevel, newInterval);

            // Use copy() instead of copyTo() because it is more efficient.
            for (DataIterator dit = dataLevel.dataIterator(); dit.ok(); ++dit)
              {
                const NodeFArrayBox& dataLevelNfab = dataLevel[dit()];
                NodeFArrayBox& vectVarsNodeLevelNfab = vectVarsNodeLevel[dit()];
                Box bx(vectVarsNodeLevelNfab.box()); // no ghosts
                vectVarsNodeLevelNfab.copy(bx,
                                           newInterval,
                                           bx,
                                           dataLevelNfab,
                                           oldInterval);
              }

            ivar++;
          }
      vectVarsNodeLevel.exchange(vectVarsNodeLevel.interval());
    }

  WriteAMRHierarchyHDF5(filename,
                        a_vectGrids,
                        vectVarsNode,
                        a_vectNames,
                        a_domain,
                        a_dx, a_dt, a_time,
                        a_vectRatio,
                        a_numLevels);

  for (int ilev = 0; ilev < a_numLevels; ilev++)
    delete vectVarsNode[ilev];
}

// ---------------------------------------------------------
void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                      const Vector<string>& a_vectNames,
                      const Box& a_domain,
                      const Real& a_dx,
                      const Real& a_dt,
                      const Real& a_time,
                      const Vector<int>& a_refRatio,
                      const int& a_numLevels,
                      const RealVect& a_origin,
                      const Interval& a_comps )
{
  Interval comps_interval( a_comps );
  if ( a_comps.size() <= 0 ) comps_interval.define( 0,a_vectData[0]->nComp() -1 );

  int nComps = comps_interval.size() ;

  CH_assert(a_numLevels > 0);
  CH_assert(a_vectData.size()  >= a_numLevels);
  CH_assert(a_refRatio.size() >= a_numLevels-1);
  CH_assert(a_vectNames.size() == nComps);

  HDF5HeaderData header;

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  // string centeringdescriptor("node");
  // header.m_string ["data_centering"] = centeringdescriptor;
  header.m_int ["data_centering"] = 7; // 7 for node-centered data
  header.m_int ["num_levels"]       = a_numLevels;
  header.m_int ["num_components"]    = nComps;
  // write the grid origin if it isn't zero
  if ( a_origin != RealVect::Zero )
  {
    header.m_realvect["origin"] = a_origin ;
  }

  for (int ivar = 0; ivar < nComps; ivar++)
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
          refLevel = a_refRatio[ilev];
      if (ilev != 0)
        {
          domainLevel.refine(a_refRatio[ilev-1]);
          dtLevel /= a_refRatio[ilev-1];
          dxLevel /= a_refRatio[ilev-1];
        }
      CH_assert(a_vectData[ilev] != NULL);
      const LevelData<NodeFArrayBox>& dataLevel = *a_vectData[ilev];
      const IntVect ghostVect(dataLevel.ghostVect());
      CH_assert(dataLevel.nComp() >= nComps);
      int eek = writeLevel(handle, ilev, dataLevel,
                           dxLevel, dtLevel,  a_time,
                           domainLevel, refLevel,
                           ghostVect, comps_interval);
      if (eek != 0)
        MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
    }
}

// ---------------------------------------------------------
void
WriteAnisotropicAMRHierarchyHDF5(
    HDF5Handle& handle,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx, // for each direction
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_refRatios,
    const int& a_numLevels,
    const RealVect& a_origin,
    const Interval& a_comps )
{
  Interval comps_interval( a_comps );
  if ( a_comps.size() <= 0 ) comps_interval.define( 0,a_vectData[0]->nComp() -1 );

  int nComps = comps_interval.size() ;

  CH_assert(a_numLevels > 0);
  CH_assert(a_vectData.size()  >= a_numLevels);
  CH_assert(a_refRatios.size() >= a_numLevels-1);
  CH_assert(a_vectNames.size() == nComps);

  HDF5HeaderData header;

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  // string centeringdescriptor("node");
  // header.m_string ["data_centering"] = centeringdescriptor;
  header.m_int ["data_centering"] = 7; // 7 for node-centered data
  header.m_int ["num_levels"]       = a_numLevels;
  header.m_int ["num_components"]    = nComps;
  // write the grid origin if it isn't zero
  if ( a_origin != RealVect::Zero )
  {
    header.m_realvect["origin"] = a_origin ;
  }

  for (int ivar = 0; ivar < nComps; ivar++)
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
          refLevel = a_refRatios[ilev];
      if (ilev != 0)
        {
          domainLevel.refine(a_refRatios[ilev-1]);
          dtLevel /= a_refRatios[ilev-1][0]; // HACK use the 0 dir ref ratio
          dxLevel /= a_refRatios[ilev-1];
        }
      CH_assert(a_vectData[ilev] != NULL);
      const LevelData<NodeFArrayBox>& dataLevel = *a_vectData[ilev];
      const IntVect ghostVect(dataLevel.ghostVect());
      CH_assert(dataLevel.nComp() >= nComps);
      int eek = writeLevel(handle, ilev, dataLevel,
                           dxLevel, dtLevel,  a_time,
                           domainLevel, refLevel,
                           ghostVect, comps_interval);
      if (eek != 0)
        MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
    }
}

// ---------------------------------------------------------
void
WriteAMRHierarchyHDF5(const string& filename,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
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

// ---------------------------------------------------------
void
WriteAMRHierarchyHDF5(HDF5Handle& handle,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
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
  // string centeringdescriptor("node");
  // header.m_string ["data_centering"] = centeringdescriptor;
  header.m_int ["data_centering"] = 7; // 7 for node-centered data
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
          refLevel = a_refRatio[ilev];
      if (ilev != 0)
        {
          domainLevel.refine(a_refRatio[ilev-1]);
          dtLevel /= a_refRatio[ilev-1];
          dxLevel /= a_refRatio[ilev-1];
        }
      CH_assert(a_vectData[ilev] != NULL);
      const LevelData<NodeFArrayBox>& dataLevel = *a_vectData[ilev];
      CH_assert(dataLevel.nComp() == nComp);
      int eek = writeLevel(handle, ilev, dataLevel,
                           dxLevel, dtLevel,  time,
                           domainLevel, refLevel);
      if (eek != 0)
        MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
    }
}

// ---------------------------------------------------------
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
                     Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<int>& a_refRatio,
                     int& a_numLevels,
                     bool a_setGhost)
{
  HDF5Handle handle;
  int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
  if ( err < 0) return -4;
  int eekflag = ReadAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                                     a_vectNames, a_domain, a_dx, a_dt,
                                     a_time, a_refRatio, a_numLevels,
                                     a_setGhost);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();

  return (eekflag);
}

// ---------------------------------------------------------
int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     Real& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<int>& a_refRatio,
                     int& a_numLevels,
                     bool a_setGhost)
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
      a_vectData[ilev] = new LevelData<NodeFArrayBox>();
      int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                          dxLevel, dtLevel,  a_time,
                          domainLevel, refLevel, Interval(), a_setGhost);
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

// ---------------------------------------------------------
int
ReadAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<int>& a_refRatio,
                     int& a_numLevels,
                     bool a_setGhost)
{
  HDF5Handle handle;
  int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
  if ( err < 0) return -4;

  int eekflag = ReadAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                                     a_domain, a_refRatio, a_numLevels,
                                     a_setGhost);
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();
  return (eekflag);
}

// ---------------------------------------------------------
int
ReadAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<NodeFArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<int>& a_refRatio,
                     int& a_numLevels,
                     bool a_setGhost)
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
      a_vectData[ilev] = new LevelData<NodeFArrayBox>();
      int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                          dxLevel, dtLevel,  time,
                          domainLevel, refLevel, Interval(), a_setGhost);
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

// put simple debugging functions here, at the end of the hdf5 stuff
static void
ChomboVisVisualizeFile(const char *fname)
{
  char command[2000];
  sprintf(command,"$CHOMBOVIS_HOME/bin/chombovis debug_level=0 %s &",fname);
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
VisualizeFile(const char *fname)
{
  const char *use_chombovis = getenv("CHOMBO_USE_CHOMBOVIS");
  if (use_chombovis)
      ChomboVisVisualizeFile(fname);
  else
      VisItVisualizeFile(fname);
}

// ---------------------------------------------------------
void
writeNFAB(const NodeFArrayBox* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = "nfab.hdf5";
  writeNFABname(a_dataPtr, fname);
}

// ---------------------------------------------------------
void
viewNFAB(const NodeFArrayBox* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeNFABname(a_dataPtr, fname);
  VisualizeFile(fname);
}

// ---------------------------------------------------------
void
writeNFABname(const NodeFArrayBox* a_dataPtr,
              const char*      a_filename)
{
  if (a_dataPtr == NULL) return;

  const NodeFArrayBox& data = *a_dataPtr;

  HDF5Handle handle(a_filename, HDF5Handle::CREATE);
  HDF5HeaderData header;

  int numlevels= 1;
  int nComp = data.getFab().nComp();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  // string centeringdescriptor("node");
  // header.m_string ["data_centering"] = centeringdescriptor;
  header.m_int ["data_centering"] = 7; // 7 for node-centered data
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
  LevelData<NodeFArrayBox> ldf(grids, nComp);
  // now copy nfab into ldf
  DataIterator dit = ldf.dataIterator();
  ldf[dit()].copy(data);

  int eek = writeLevel(handle, 0, ldf, dxLevel, dtLevel, time,
                       domainLevel, refLevel);
  if (eek != 0)
    MayDay::Error("writeFABname: error in writeLEvel");

  handle.close();
}

// ---------------------------------------------------------
void
writeNodeLevel(const LevelData<NodeFArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = "LDF.hdf5";
  writeNodeLevelname(a_dataPtr, fname);
}

// ---------------------------------------------------------
void
viewNodeLevel(const LevelData<NodeFArrayBox>* a_dataPtr)
{
  if (a_dataPtr == NULL) return;

  const char* fname = tempnam(NULL,NULL);
  writeNodeLevelname(a_dataPtr, fname);
  VisualizeFile(fname);
}

// ---------------------------------------------------------
void
writeNodeLevelname(const LevelData<NodeFArrayBox>* a_dataPtr,
                   const char*                 a_filename)
{
  if (a_dataPtr == NULL) return;

  const LevelData<NodeFArrayBox>& data = *a_dataPtr;

  HDF5Handle handle(a_filename, HDF5Handle::CREATE);
  HDF5HeaderData header;

  int numlevels = 1;
  int nComp = data.nComp();

  string filedescriptor("VanillaAMRFileType");
  header.m_string ["filetype"]      = filedescriptor;
  // string centeringdescriptor("node");
  // header.m_string ["data_centering"] = centeringdescriptor;
  header.m_int ["data_centering"] = 7; // 7 for node-centered data
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

  // need to figure out what domain will contain this LevelData
  // This must be LayoutIterator instead of DataIterator because
  // we want domain over ALL procs.
  const DisjointBoxLayout& levelBoxes = data.getBoxes();
  LayoutIterator lit = levelBoxes.layoutIterator();
  lit.reset();
  Box domain = levelBoxes.get(lit());
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

  // put bogus numbers here
  Real dtLevel = 1.0;
  Real dxLevel = 1.0;
  Real time = 1.0;

  int refLevel = 1;

  IntVect ghostVect = data.ghostVect();
  int eek = writeLevel(handle, 0, data, dxLevel, dtLevel, time,
                       domain, refLevel, ghostVect);

  if (eek != 0)
    MayDay::Error("WriteNodeLevelname: error in writeNodeLevel");

  handle.close();
}

// ---------------------------------------------------------
void
WritePartialAMRHierarchyHDF5(const string& filename,
                             const Vector<DisjointBoxLayout>& a_vectGrids,
                             const Vector<LevelData<NodeFArrayBox>* > & a_vectData,
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
  Vector<LevelData<NodeFArrayBox> * > newVectData(numLevels);
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

// ---------------------------------------------------------
void
WritePartialAMRHierarchyHDF5(const string& filename,
                             const Vector<DisjointBoxLayout>& a_vectGrids,
                             const Vector<Vector<LevelData<NodeFArrayBox>* > >& a_vectData,
                             const Vector<string>& a_vectNames,
                             const Box& a_baseDomain,
                             const Real& a_baseDx,
                             const Real& a_dt,
                             const Real& a_time,
                             const Vector<int>& a_vectRatio,
                             const Interval& a_levels)
{
  // We have a_vectData[component][level].
  // Each component should have the same number of levels.
  int numvarsOrig = a_vectData.size();
  int numvars = 0;
  Vector<int> numvarsEach(numvarsOrig);
  // need at least one set of data.
  // let numlevels by the number of levels in the first set.
  CH_assert(a_vectData.size() >= 1);
  int numlevels = a_vectData[0].size();
  for (int ivarOrig = 0; ivarOrig < numvarsOrig; ivarOrig++)
    {
      CH_assert( a_vectData[ivarOrig].size() >= numlevels );
      numvarsEach[ivarOrig] = a_vectData[ivarOrig][0]->nComp();
      numvars += numvarsEach[ivarOrig];
    }

  CH_assert( a_vectNames.size() == numvars );

  Vector<LevelData<NodeFArrayBox>* > vectVarsNode(numlevels, NULL);
  for (int ilev = 0; ilev < numlevels; ilev++)
    {
      const DisjointBoxLayout& grids = a_vectGrids[ilev];
      // What about ghosts in the LevelData?
      // Different variables may have different ghost vectors.
      // Solution:  dispense with ghosts.
      vectVarsNode[ilev] = new LevelData<NodeFArrayBox>(grids, numvars);
      LevelData<NodeFArrayBox>& vectVarsNodeLevel = *vectVarsNode[ilev];
      int ivar = 0;
      for (int ivarOrig = 0; ivarOrig < numvars; ivarOrig++)
        for (int ivarEach = 0; ivarEach < numvarsEach[ivarOrig]; ivarEach++)
          {
            Interval newInterval(ivar, ivar);
            Interval oldInterval(ivarEach, ivarEach);

            const LevelData<NodeFArrayBox>& dataLevel = *a_vectData[ivar][ilev];
            // dataLevel.copyTo(oldInterval, vectVarsNodeLevel, newInterval);

            // Use copy() instead of copyTo() because it is more efficient.
            for (DataIterator dit = dataLevel.dataIterator(); dit.ok(); ++dit)
              {
                const NodeFArrayBox& dataLevelNfab = dataLevel[dit()];
                NodeFArrayBox& vectVarsNodeLevelNfab = vectVarsNodeLevel[dit()];
                Box bx(vectVarsNodeLevelNfab.box()); // no ghosts
                vectVarsNodeLevelNfab.copy(bx,
                                           newInterval,
                                           bx,
                                           dataLevelNfab,
                                           oldInterval);
              }

            ivar++;
          }
    }

  WritePartialAMRHierarchyHDF5(filename,
                               a_vectGrids,
                               vectVarsNode,
                               a_vectNames,
                               a_baseDomain,
                               a_baseDx, a_dt, a_time,
                               a_vectRatio,
                               a_levels);

  for (int ilev = 0; ilev < numlevels; ilev++)
    delete vectVarsNode[ilev];
}

#include "NamespaceFooter.H"

#endif // CH_USE_HDF5
