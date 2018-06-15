#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BLOCKWRITEI_H_
#define _BLOCKWRITEI_H_

#include "LevelData.H"
#include "HDF5Portable.H"
#include "CH_HDF5.H"
#include <string>
#include <map>
#include "RealVect.H"
#include "CH_Timer.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "Vector.H"
#include "memtrack.H"
#include "FluxBox.H"
#include "SPMD.H"

//==================================================================//
template <class T> void
blockLocalOffsets(Vector<long long>&      a_localSizes, //=numpts*ncomp
                                                        //for each box
                  long long        &      a_localTotalSize, //=numpts*ncomp
                  Vector<Box>      &      a_localBoxes,
                  const BoxLayoutData<T>& a_data,
                  const Interval&         a_comps,
                  const IntVect&          a_outputGhost)
{
  if (!(T::preAllocatable()==0))
    {
      MayDay::Error("non static preallocatable data not covered yet.");
    }

  Vector<int> thisSize(1);
  const BoxLayout& layout = a_data.boxLayout();
  T dummy;
  unsigned int curIndex = 0;
  DataIterator it(layout.dataIterator());

  a_localTotalSize = 0;
  a_localSizes.resize(it.size());
  a_localBoxes.resize(it.size());
  for (it.reset(); it.ok(); ++it)
    {
      long long  curSize = 0;
      const Box& curBox = layout[it()];
      Box region = curBox;
      region.grow(a_outputGhost);
      dataSize(dummy, thisSize, region, a_comps);

      curSize =  thisSize[0];
      a_localSizes[curIndex] = curSize;
      a_localBoxes[curIndex] = curBox;

      a_localTotalSize += thisSize[0];
      curIndex    += 1;
    }
}
//==================================================================//
template <class T> void
blockWriteToBuffer(void*                    a_buffer,
                   const Vector<long long>& a_sizes,  //ncomp*npts
                   const BoxLayoutData<T>&  a_data,
                   const Interval&          a_comps,
                   const IntVect&           a_outputGhost)
{
  int* curPtr = (int*)(a_buffer);
  int curIndex = 0;
  for (DataIterator it = a_data.dataIterator(); it.ok(); ++it)
    {
      Vector<void*> buffers(1, (void*) curPtr);
      const T& data = a_data[it()];
      //unsigned int index = a_data.boxLayout().index(it());
      Box box = a_data.box(it());
      box.grow(a_outputGhost);
      write(data, buffers, box, a_comps); //write T to buffer

      curPtr   += a_sizes[curIndex];
      curIndex += 1;
    }

}
//==================================================================//
int
gatherBoxesAndOffsets(long long&               a_offsetThisProc,
                      long long&               a_allProcSize,
                      Vector<long long>&       a_globalOffsets,
                      Vector<Box>&             a_globalBoxes,
                      const Vector<long long>& a_localBoxSizes,
                      const Vector<Box>&       a_localBoxes,
                      const long long&         a_localAllBoxSize) //= numpts * ncomp

{
  Vector<Vector<Box      > > allBoxes;
  Vector<Vector<long long> > allBoxSizes;
  Vector<long long>          allSizes;
  int iprocdest = 0;
  gather(allBoxes,    a_localBoxes,        iprocdest);
  gather(allBoxSizes, a_localBoxSizes,     iprocdest);
  gather(allSizes,    a_localAllBoxSize,   iprocdest);
  Vector<long long> offsetsProcs(numProc());
  if (procID() == iprocdest)
    {
      a_allProcSize = 0;
      long long numBoxes = 0;
      for (int iproc = 0; iproc < numProc(); iproc++)
        {
          a_allProcSize += allSizes[iproc];
          numBoxes += allBoxes[iproc].size();
        }
      a_globalBoxes.  resize(numBoxes);
      a_globalOffsets.resize(numBoxes+1); //need the last one to hold
                                          //the size
      long long  curOffset = 0;
      int        curIndex  = 0;
      for (int iproc = 0; iproc < numProc(); iproc++)
        {
          for (int ivec = 0; ivec < allBoxSizes[iproc].size(); ivec++)
            {
              if (ivec == 0)
                {
                  offsetsProcs[iproc] = curOffset;
                }
              long long boxSize = allBoxSizes[iproc][ivec];
              Box       curBox  =    allBoxes[iproc][ivec];

              a_globalBoxes  [curIndex] = curBox;
              a_globalOffsets[curIndex] = curOffset;

              curIndex      += 1;
              curOffset     += boxSize;
            }
        }
      if (curOffset != a_allProcSize)
        {
          return -43;  //offsets should add up to size
        }
      a_globalOffsets[curIndex] = curOffset; //make last one the total size

    }
  a_offsetThisProc = offsetsProcs[procID()];
  broadcast(a_allProcSize,   iprocdest);
  broadcast(a_globalOffsets, iprocdest);
  broadcast(a_globalBoxes,   iprocdest);
  return 0;
}
//==================================================================//
int
blockWriteBufferToFile(HDF5Handle&         a_handle,
                       void*               a_buffer,
                       const std::string&  a_name,
                       Vector<Box>&        a_boxes,
                       Vector<long long>&  a_sizes,  //ncomp*npts
                       const Vector<hid_t>& a_types,
                       const BoxLayout&    a_layout,
                       const long long&    a_thisprocsize)
{

  herr_t err;
  //need the offset into the entire data set for my proc
  long long offsetthisproc, sumprocsize;
  Vector<Box>       globalBoxes;
  Vector<long long> globalOffsets;
  //write boxes and offsets
  //get all procs boxes
  int ret = gatherBoxesAndOffsets(offsetthisproc, sumprocsize, globalOffsets, globalBoxes,  a_sizes, a_boxes, a_thisprocsize);
  if (ret!=0)
    {
      return ret;
    }

  if (globalBoxes.size()   != a_layout.size())
    {
      return -12;
    }

  {//offset write
    char offsetname[1024];
    sprintf(offsetname, "%s:offsets=0",a_name.c_str());

    hsize_t flatdims[1];
    flatdims[0] = globalOffsets.size();
    hid_t offsetspace  =  H5Screate_simple(1, flatdims, NULL);
    hid_t offsetdata   =  H5Dcreate(a_handle.groupID(), offsetname,
                                    H5T_NATIVE_LLONG, offsetspace, H5P_DEFAULT);

    if (procID() == 0)
      {
        hid_t memdataspace =  H5Screate_simple(1, flatdims, NULL);
        err = H5Dwrite(offsetdata, H5T_NATIVE_LLONG, memdataspace, offsetspace,
                       H5P_DEFAULT, &(globalOffsets[0]));
        if (err < 0) return  err;
        H5Sclose(memdataspace);
      }

    H5Sclose(offsetspace);
    H5Dclose(offsetdata);
  }
  { //data write
    hsize_t hs_procsize[1];
    hsize_t  hs_allprocsize[1];
    ch_offset_t ch_offset[1];

    ch_offset[0]      = offsetthisproc;
    hs_procsize[0]    = a_thisprocsize; //size of buffer on this proc
    hs_allprocsize[0] =     sumprocsize; //size of buffer on all procs

    char dataname[1024];
    sprintf(dataname, "%s:datatype=0",a_name.c_str());

    hid_t dataspace    = H5Screate_simple(1, hs_allprocsize, NULL);
    hid_t memdataspace = H5Screate_simple(1, hs_procsize, NULL);

    hid_t dataset   = H5Dcreate(a_handle.groupID(), dataname,
                                a_types[0], dataspace, H5P_DEFAULT);

    //select where in the file it will be written
    err =  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                               ch_offset,     NULL,
                               hs_procsize,   NULL);
    if (err < 0) return  err;

    err = H5Dwrite(dataset, a_types[0], memdataspace, dataspace,
                   H5P_DEFAULT,  a_buffer);
    if (err < 0) return  err;

    H5Sclose(memdataspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
  return 0;
}

template <class T>
int blockWrite(HDF5Handle&         a_handle,
               const LevelData<T>& a_data,
               const std::string&  a_name,
               const IntVect&      a_outputGhost,
               const Interval&     a_in_comps)
{
  HDF5HeaderData info;
  info.m_intvect["ghost"] = a_data.ghostVect();
  IntVect og(a_outputGhost);
  og.min(a_data.ghostVect());
  info.m_intvect["outputGhost"] = og;
  std::string group = a_handle.getGroup();
  a_handle.setGroup(group+"/"+a_name+"_attributes");
  info.writeToFile(a_handle);
  a_handle.setGroup(group);
  return blockWrite(a_handle, (const BoxLayoutData<T>&)a_data, a_name, og, a_in_comps);
}

//==================================================================//
template <class T> int
blockWrite(HDF5Handle&             a_handle,
           const BoxLayoutData<T>& a_data,
           const std::string&      a_name,
           const IntVect&          a_outputGhost,
           const Interval&         a_comps)
{
  CH_TIME("block write Level");
  int ret = 0;

  const BoxLayout& layout = a_data.boxLayout();
  T dummy; // used for preallocatable methods for dumb compilers.

  Vector<hid_t> types;
  dataTypes(types, dummy);
  if (types.size() != 1)
    {
      MayDay::Error("premature generality, I only deal with types of size 1");
    }

  //all this stuff is local to this proc
  Vector<long long> sizes;     //=numpts*ncomp for each box
  Vector<Box>       boxes;
  long long         totalsize; //=numpts*ncomp

  //calculate offset stuff on this proc and allocate and output to
  //local buffer
  blockLocalOffsets(sizes,  totalsize, boxes,  a_data, a_comps, a_outputGhost);

  //allocate buffer --- this totalsize is the local size
  size_t type_size = H5Tget_size(types[0]);

  //the argument is the number of chars in the buffer length
  void* buffer = mallocMT(totalsize*type_size);

  //write to the buffer
  blockWriteToBuffer(buffer, sizes, a_data, a_comps, a_outputGhost);

  //write attributes
  HDF5HeaderData info;
  info.m_int["comps"] = a_comps.size();
  info.m_string["objectType"] = name(dummy); //this is a function call
  std::string group = a_handle.getGroup();

  a_handle.setGroup(group+"/"+a_name+"_attributes");
  info.writeToFile(a_handle);
  a_handle.setGroup(group);

  //exchange offset information
  //output stuff to file
  ret = blockWriteBufferToFile(a_handle,
                               buffer,
                               a_name,
                               boxes,
                               sizes,
                               types,
                               layout,
                               totalsize);
  if (ret != 0)
    {
      return ret;
    }
//  ret = blockWriteBufferToFile(a_handle, buffer, a_name, boxes, offsets, layout, totalsize);

  //free the buffer memory
  freeMT(buffer);
  buffer = NULL;

  return ret;
}

template <class T> int
blockWriteLevel(HDF5Handle&            a_handle,
                const int&             a_level,
                const LevelData<T>&    a_data,
                const Real&            a_dx,
                const Real&            a_dt,
                const Real&            a_time,
                const Box&             a_domain,
                const int&             a_refRatio,
                const IntVect&         a_outputGhost,
                const Interval&        a_comps)
{
  int error;
  char levelName[10];
  std::string currentGroup = a_handle.getGroup();
  sprintf(levelName, "/level_%i",a_level);
  error = a_handle.setGroup(currentGroup + levelName);
  if (error != 0) return 1;

  HDF5HeaderData meta;
  meta.m_real["dx"] = a_dx;
  meta.m_real["dt"] = a_dt;
  meta.m_real["time"] = a_time;
  meta.m_box["prob_domain"] = a_domain;
  meta.m_int["ref_ratio"] = a_refRatio;

  error = meta.writeToFile(a_handle);
  if (error != 0) return 2;

  error = write(a_handle, a_data.boxLayout());
  if (error != 0) return 3;

  error = blockWrite(a_handle, a_data, std::string("data"), a_outputGhost, a_comps);
  if (error != 0) return 4;

  a_handle.setGroup(currentGroup);

  return 0;
}
//==================================================================//
template <class T> int
blockReadLevel(HDF5Handle&     a_handle,
               const int&      a_level,
               LevelData<T>&   a_data,
               Real&           a_dx,
               Real&           a_dt,
               Real&           a_time,
               Box&            a_domain,
               int&            a_refRatio,
               const Interval& a_comps,
               bool            a_setGhost)
{
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  //unused
  // int nComp = header.m_int["num_components"];

  int error;
  char levelName[10];
  std::string currentGroup = a_handle.getGroup();
  sprintf(levelName, "/level_%i",a_level);
  error = a_handle.setGroup(currentGroup + levelName);
  if (error != 0) return 1;

  HDF5HeaderData meta;
  error = meta.readFromFile(a_handle);
  if (error != 0) return 2;
  a_dx       = meta.m_real["dx"];
  a_dt       = meta.m_real["dt"];
  a_time     = meta.m_real["time"];
  a_domain   = meta.m_box["prob_domain"];
  a_refRatio = meta.m_int["ref_ratio"];

  //read in the boxes and load balance
  Vector<Box> boxes;
  error = read(a_handle, boxes);
  Vector<int> procIDs;
  LoadBalance(procIDs, boxes);

  //make the layout
  DisjointBoxLayout layout(boxes, procIDs);

  layout.close();
  if (error != 0) return 3;
  //do the very complicated read
  error = blockRead(a_handle, a_data, "data", layout, a_comps);

  if (error != 0) return 4;

  //return the handle to its input state
  a_handle.setGroup(currentGroup);

  return 0;
}

//////////
template <class T> int
blockRead(HDF5Handle& a_handle, LevelData<T>& a_data, const std::string& a_name,
          const DisjointBoxLayout& a_layout, const Interval& a_comps, bool a_redefineData)
{
  if (a_redefineData)
    {
      HDF5HeaderData info;
      std::string group = a_handle.getGroup();
      if (a_handle.setGroup(group+"/"+a_name+"_attributes"))
        {
          std::string message = "error opening "
                                +a_handle.getGroup()+"/"+a_name+"_attributes" ;
          MayDay::Warning(message.c_str());
          return 1;
        }
      info.readFromFile(a_handle);
      a_handle.setGroup(group);
      int ncomp =  info.m_int["comps"];
      IntVect ghost = info.m_intvect["ghost"];
      if (a_comps.end() > 0 && ncomp < a_comps.end())
        {
          MayDay::Error("attempt to read component interval that is not available");
        }
      if (a_comps.size() == 0)
        a_data.define(a_layout, ncomp, ghost);
      else
        a_data.define(a_layout, a_comps.size(), ghost);
    }
  return blockRead(a_handle, (BoxLayoutData<T>&)a_data, a_name, a_layout, a_comps, false);

}

template <class T>
int
blockRead(HDF5Handle& a_handle, BoxLayoutData<T>& a_data, const std::string& a_name,
         const BoxLayout& a_layout, const Interval& a_comps, bool a_redefineData)
{

  int ret = read(a_handle, a_data, a_name, a_layout, a_comps, a_redefineData);
  return ret;
}

#endif
