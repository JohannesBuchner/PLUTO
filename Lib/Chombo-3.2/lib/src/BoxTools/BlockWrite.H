#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BLOCKWRITE_H_
#define _BLOCKWRITE_H_

#ifdef CH_USE_HDF5  // if you don't have CH_USE_HDF5, then this file is useless

#include "CH_HDF5.H"

#include <iostream>
using std::cout;
using std::endl;

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "LevelData.H"
#include "HDF5Portable.H"
#include <string>
#include <map>
#include "RealVect.H"
#include "CH_Timer.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "Vector.H"
#include "memtrack.H"
#include "FluxBox.H"
#include "CH_HDF5.H"

#include "NamespaceHeader.H"

///
template <class T> int
blockReadLevel(HDF5Handle& a_handle,
               const int& a_level,
               LevelData<T>& a_data,
               Real& a_dx,
               Real& a_dt,
               Real& a_time,
               Box& a_domain,
               int& a_refRatio,
               const Interval& a_comps,
               bool setGhost);

///
template <class T> int
blockWriteLevel(HDF5Handle& a_handle,
                const int& a_level,
                const LevelData<T>& a_data,
                const Real& a_dx,
                const Real& a_dt,
                const Real& a_time,
                const Box& a_domain,
                const int& a_refRatio,
                const IntVect& outputGhost,
                const Interval& comps);

///
template <class T> int
blockWrite(HDF5Handle&             a_handle,
           const BoxLayoutData<T>& a_data,
           const std::string&      a_name,
           const IntVect&          a_outputGhost,
           const Interval&         a_comps);

///
template <class T>
int
blockWrite(HDF5Handle&         a_handle,
           const LevelData<T>& a_data,
           const std::string&  a_name,
           const IntVect&      a_outputGhost,
           const Interval&     a_in_comps);

///
template <class T> void
blockLocalOffsets(Vector<long long>&      a_localOffsets,
                  long long        &      a_localTotalSize,
                  Vector<Box>      &      a_localBoxes,
                  const BoxLayoutData<T>& a_data,
                  const Interval&         a_comps,
                  const IntVect&          a_outputGhost);

///
template <class T> void
blockWriteToBuffer(void*                   a_buffer,
                   const BoxLayoutData<T>& a_data,
                   const Interval&         a_comps,
                   const IntVect&          a_outputGhost);

///
int
gatherBoxesAndOffsets(long long&               a_offsetThisProc,
                      long long&               a_allProcSize,
                      Vector<long long>&       a_globalOffsets,
                      Vector<Box>&             a_globalBoxes,
                      const Vector<long long>& a_localOffsets,
                      const Vector<Box>&       a_localBoxes,
                      const long long&         a_localSize);

///
int
blockWriteBufferToFile(HDF5Handle&         a_handle,
                       void*               a_buffer,
                       const std::string&  a_name,
                       Vector<Box>&        a_boxes,
                       Vector<long long>&  a_offsets,
                       const Vector<hid_t>& a_types,
                       const BoxLayout&    a_layout,
                       const long long&    a_thisprocsize);

/// read BoxLayoutData named a_name from location specified by a_handle.
/**
    Read BoxLayoutData named a_name from location specified by a_handle.  User must supply the correct BoxLayout for this function if redefineData == true.  \\
    returns: success:      0\\
    bad location: 1\\
    HDF5 error:   negative error code.\\
*/
template <class T>
int blockRead(HDF5Handle&        a_handle,
              BoxLayoutData<T>&  a_data,
              const std::string& a_name,
              const BoxLayout&   a_layout,
              const Interval&    a_comps = Interval(),
              bool redefineData = true);

/// read LevelData named a_name from location specified by a_handle.
/**

    Read LevelData named a_name from location specified by a_handle.
    User must supply the correct BoxLayout for this function.

    Arg a_name is significant: the HDF5 group to which a_handle is set is
    is assumed to contain a dataset called <a_name>:datatype=<some integer>,
    a dataset called <a_name>:offsets=<some integer>, and a subgroup named
    <a_name>_attributes.  You will have all these items if you dumped your
    LevelData out using the corresponding write() function defined here.

    If arg redefineData==false, then the user must pass in a valid  LevelData.
    Otherwise, this function figures out how many components and ghost cells there
    are, and allocates the correct amount of space.  The actual FArray data held
    by the LevelData gets filled in here, regardless of redefineData; "redefine"
    alludes to the family of define() functions.

    returns: success:      0\\
    bad location: 1\\
    HDF5 error:   negative error code.\\
*/
template <class T>
int blockRead(HDF5Handle& a_handle,
              LevelData<T>& a_data,
              const std::string& a_name,
              const DisjointBoxLayout& a_layout,
              const Interval& a_comps = Interval(),
              bool redefineData = true);

//implementation
#include "BlockWriteI.H"

#endif  //ifdef   CH_USE_HDF5
#endif  //ifndef __BlockWrite__
