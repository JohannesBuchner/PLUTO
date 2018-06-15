#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_USE_HDF5

#include "CH_HDF5.H"
#include "MayDay.H"
#include <cstdio>
#include "parstream.H"
#include "NamespaceHeader.H"
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

//----------------------------------------------------------
// template <class T>
// ostream& operator<<(ostream& os, const Vector<T>& vec)
// {
//   for (int i=0; i<vec.size(); i++) os<<vec[i]<<" ";
//   return os;
// }

void OffsetBuffer::operator=(const OffsetBuffer& rhs)
{
  if (rhs.index.size() == 0)
    {
      index.resize(0);
      offsets.resize(0);
      return;
    }
  index.resize(rhs.index.size());
  int types = rhs.offsets[0].size();
  offsets.resize(rhs.offsets.size(), Vector<int>(types, 0));
  for (int i=0; i<rhs.index.size(); i++)
    {
      index[i] = rhs.index[i];
      for (int j=0; j<types; j++)
        {
          offsets[i][j] = rhs.offsets[i][j];
        }
    }
}

ostream&  operator<<(ostream& os, const OffsetBuffer& ob)
{
  os << ob.index <<"|"<<ob.offsets<<"\n";
  return os;
}

#include "NamespaceFooter.H"
#include "UsingNamespace.H"
#include "BaseNamespaceHeader.H"

//OffsetBuffer specialization of linearSize
template < >
int linearSize(const OffsetBuffer& a_input)
{
  int r = 0;
  r += 2;   // number of entries from this processor, number of types.
  r += a_input.index.size();// indices for these offset entries
  //size==0 can happen when no boxes to a proc.---dtg
  if (a_input.offsets.size() > 0)
    {
      r += a_input.offsets[0].size()*a_input.offsets.size();  // offsets
    }
  return r * sizeof(int);
}

//OffsetBuffer specialization of linearIn
template < >
void linearIn(OffsetBuffer& a_outputT, const void* const a_inBuf)
{
  // pout() << "in:  ";
  //guaranteed to have at least two ints.
  const int* data = (const int*) a_inBuf;
  int num = data[0],  numTypes = data[1];
  a_outputT.index.resize(data[0]);
  a_outputT.offsets.resize(data[0], Vector<int>(data[1]));
  const int* off = data+2+num;
  for (int i=0; i<num; i++)
    {
      a_outputT.index[i] = data[i+2];
      Vector<int>& offsets = a_outputT.offsets[i];
      for (int j=0; j<numTypes; j++)
        {
          offsets[j] = off[j];
          //  pout() << offsets[j]<<" ";
        }
      off+= numTypes;
    }
  //  pout() << endl;
}

//OffsetBuffer specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const OffsetBuffer& a_inputT)
{
  //  pout() << "out: ";
  int* data = (int*) a_outBuf;
  data[0] =  a_inputT.index.size();
  //size==0 can happen when no boxes to a proc.---dtg
  if (a_inputT.offsets.size() > 0)
    {
      data[1] =  a_inputT.offsets[0].size();
    }
  else
    {
      data[1] =  0;
    }
  data+=2;
  for (int i=0; i< a_inputT.index.size(); ++i, ++data)
    *data =  a_inputT.index[i];

  for (int i=0; i<a_inputT.offsets.size(); ++i)
    {
      const Vector<int>& offset = a_inputT.offsets[i];
      for (int t=0; t < offset.size(); ++t, ++data)
        {
          *data = offset[t];
          //  pout() << *data <<" ";
        }
    }

  // pout() << endl;
}
//Vector<OffsetBuffer>  specialization
template < > int linearSize(const Vector<OffsetBuffer>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<OffsetBuffer>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<OffsetBuffer>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

#include "BaseNamespaceFooter.H"
#include "NamespaceHeader.H"

int write(HDF5Handle& a_handle, const BoxLayout& a_layout, const std::string& name)
{
  CH_assert(a_layout.isClosed());
//  herr_t ret;
//  ch_offset_t offset[1];
  hsize_t  flatdims[1], count[1];

  flatdims[0] = a_layout.size();

  hid_t boxdataspace = H5Screate_simple(1, flatdims, NULL);

  H5E_auto_t efunc; void* edata; // turn auto error messaging off

#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  H5Gunlink(a_handle.groupID(), name.c_str()); 
  //removes a pre-existing dataset.
  H5Eset_auto(efunc, edata);
  
  hid_t boxdataset   = H5Dcreate(a_handle.groupID(), name.c_str(),  
				 a_handle.box_id,
                                 boxdataspace, H5P_DEFAULT);
#else
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  H5Ldelete(a_handle.groupID(), name.c_str(),H5P_DEFAULT); 
  //removes a pre-existing dataset.
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);

  hid_t boxdataset   = H5Dcreate2(a_handle.groupID(), name.c_str(),  
				  a_handle.box_id,
				  boxdataspace, H5P_DEFAULT,
				  H5P_DEFAULT,H5P_DEFAULT);
#endif
  if (boxdataset < 0) return boxdataset;

  int  proc = procID();
  hid_t procdataset, procdataspace;
  createDataset(procdataset, procdataspace, a_handle, "Processors",
                &proc, a_layout.size());

  /*
  // write boxes in parallel
  count[0] = 1;
  hid_t memdataspace = H5Screate_simple(1, count, NULL);
  for (DataIterator it = a_layout.dataIterator(); it.ok(); ++it)
    {
      offset[0] = a_layout.index(it());
      ret = H5Sselect_hyperslab (boxdataspace, H5S_SELECT_SET, offset, NULL,
                                 count, NULL);
      if (ret < 0) abort();
      ret = H5Dwrite(boxdataset, a_handle.box_id, memdataspace, boxdataspace,
                     H5P_DEFAULT, &(a_layout[it()]));
      if (ret < 0) abort();

      ret = H5Sselect_hyperslab (procdataspace, H5S_SELECT_SET, offset, NULL,
                                 count, NULL);
      if (ret < 0) abort();
      ret = H5Dwrite(procdataset, H5T_NATIVE_INT, memdataspace, procdataspace,
                     H5P_DEFAULT, &proc);
      if (ret < 0) abort();

    }
  */

  //instead, write boxes serially from proc 0
  if (proc == 0)
  {
    Vector<Box> vbox(a_layout.size());
    Vector<int> pid(a_layout.size());
    count[0]=a_layout.size();
    int b=0;
    hid_t memdataspace = H5Screate_simple(1, count, NULL);
    for (LayoutIterator it = a_layout.layoutIterator(); it.ok(); ++it)
      {
        vbox[b]=a_layout.get(it());
        pid[b]=a_layout.procID(it());
        b++;
      }
    if (vbox.size() > 0)
      {
        H5Dwrite(boxdataset, a_handle.box_id, memdataspace, boxdataspace,
                 H5P_DEFAULT, &(vbox[0]));
        H5Dwrite(procdataset, H5T_NATIVE_INT, memdataspace, procdataspace,
                 H5P_DEFAULT, &(pid[0]));
      }
    else
      {
        H5Dwrite(boxdataset, a_handle.box_id, memdataspace, boxdataspace,
                 H5P_DEFAULT, NULL);
        H5Dwrite(procdataset, H5T_NATIVE_INT, memdataspace, procdataspace,
                 H5P_DEFAULT, NULL);
      }

    H5Sclose(memdataspace);
  }

  H5Sclose(boxdataspace);
  H5Dclose(boxdataset);
  H5Sclose(procdataspace);
  H5Dclose(procdataset);

  return 0;
}

int read(HDF5Handle& a_handle, Vector<Box>& boxes, const std::string& name)
{
#ifdef H516
  hid_t  boxdataset = H5Dopen(a_handle.groupID(), name.c_str());
#else
  hid_t  boxdataset = H5Dopen2(a_handle.groupID(), name.c_str(), 
			       H5P_DEFAULT);
#endif
  if (boxdataset < 0) return boxdataset;
  hid_t boxdataspace =  H5Dget_space(boxdataset);
  if (boxdataspace < 0) return boxdataspace;
  hsize_t dims[1], maxdims[1];
  H5Sget_simple_extent_dims(boxdataspace, dims, maxdims);

  hid_t memdataspace = H5Screate_simple(1, dims, NULL);

  Box* rawboxes = new Box[dims[0]];
  if (rawboxes == NULL)
    MayDay::Error("out of memory in read(HDF5Handle&, Vector<Box>&, string)");

  herr_t error = H5Dread(boxdataset, a_handle.box_id, memdataspace, boxdataspace,
                         H5P_DEFAULT, rawboxes);
  if (error < 0) return error;

  boxes.resize(dims[0]);
  for (int index = 0; index < dims[0]; ++index)
    {
      rawboxes[index].computeBoxLen();
      boxes[index] = rawboxes[index];
    }

  delete[] rawboxes;
  H5Dclose(boxdataset);
  H5Sclose(boxdataspace);
  H5Sclose(memdataspace);
  return 0;
}

int readBoxes(HDF5Handle& a_handle, Vector<Vector<Box> >& boxes)
{
  int error;
  char levelName[100];
  std::string currentGroup = a_handle.getGroup();

  // we want to start in root group.  currentGroup reset on normal exit
  error = a_handle.setGroup("/");
  if (error != 0) return 1;
  std::string workingGroup = a_handle.getGroup();

  HDF5HeaderData header;
  header.readFromFile(a_handle);
  int numLevels = header.m_int["num_levels"];
  boxes.resize(numLevels);

  for (int i = 0; i<numLevels; i++)
  {
        sprintf(levelName, "/level_%i",i);
        error = a_handle.setGroup(workingGroup + levelName);
        if (error != 0) return 1;
        read(a_handle, boxes[i]);
  }

  a_handle.setGroup(currentGroup); // put handle back to where it stated.

  return 0;
}

int readFArrayBox(HDF5Handle& a_handle,
                                  FArrayBox&  a_fab,
                                  int a_level,
                                  int a_boxNumber,
                                  const Interval& a_components,
                                  const std::string& a_dataName)
{
  char levelName[100];
  herr_t err;
  int error=0;
  hsize_t count[1];
  ch_offset_t offset[1];
  hsize_t flatdims[1];
  hid_t dataspace, dataset, memdataspace;

  std::string currentGroup = a_handle.getGroup();

  // we want to start in root group.  currentGroup reset on normal exit
  error = a_handle.setGroup("/");
  if (error != 0) return 1;
  std::string workingGroup = a_handle.getGroup();

  sprintf(levelName, "/level_%i",a_level);
  error = a_handle.setGroup(workingGroup + levelName);
  if (error != 0) return 1;

  // first, fetch the correct box from file.

  Box box;
  flatdims[0]  = 1;
  memdataspace = H5Screate_simple(1, flatdims, NULL);
  count[0]  = 1;
  offset[0] = a_boxNumber;
#ifdef H516
  dataset   = H5Dopen(a_handle.groupID(),  "boxes");
#else
  dataset   = H5Dopen2(a_handle.groupID(),  "boxes",H5P_DEFAULT);
#endif
  dataspace = H5Dget_space(dataset);
  err =  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                             offset, NULL,
                             count, NULL);
  if (err != 0 ) return 2;
  err = H5Dread(dataset, HDF5Handle::box_id, memdataspace, dataspace,
                H5P_DEFAULT, &box);
  if (err != 0 ) return 3;

  H5Dclose(dataset);
  H5Sclose(dataspace);

  box.computeBoxLen();

  //next, fetch the data offset that we need from the file.
  std::string namebuff1 = a_dataName + ":offsets=0";
#ifdef H516
  dataset = H5Dopen(a_handle.groupID(), namebuff1.c_str());
#else
  dataset = H5Dopen2(a_handle.groupID(), namebuff1.c_str(), H5P_DEFAULT);
#endif
  dataspace = H5Dget_space(dataset);
  err =  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                             offset, NULL,
                             count, NULL);
  if (err != 0 ) return 2;

  long long offsetLong;
  err = H5Dread(dataset, H5T_NATIVE_LLONG, memdataspace, dataspace,
                H5P_DEFAULT, &offsetLong);
  if (err != 0 ) return 3;

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(memdataspace);

  //  OK, now lets go get our FArrayBox data:
  int ncomp = a_components.size();
  if (a_fab.box() == box && ncomp == a_fab.nComp())
  {
    //FArrayBox already defined appropriately for this operation
  }
  else
  {
    a_fab.define(box, ncomp);
  }

  offsetLong += box.numPts() * a_components.begin();
  count[0]    = box.numPts() * ncomp;
  flatdims[0] = count[0];
  offset[0]   = offsetLong;
  memdataspace = H5Screate_simple(1, flatdims, NULL);
  std::string namebuff2 = a_dataName + ":datatype=0";
#ifdef H516
  dataset  = H5Dopen(a_handle.groupID(),  namebuff2.c_str() );
#else
  dataset  = H5Dopen2(a_handle.groupID(),  namebuff2.c_str() , H5P_DEFAULT);
#endif
  dataspace = H5Dget_space(dataset);
  err =  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                             offset, NULL, count, NULL);
  if (err != 0 ) return 2;
  err = H5Dread(dataset, H5T_NATIVE_REAL, memdataspace, dataspace,
                H5P_DEFAULT, a_fab.dataPtr());
  if (err != 0 ) return 3;

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(memdataspace);

  a_handle.setGroup(currentGroup);
  return 0;

}

bool HDF5Handle::initialized = false;
hid_t HDF5Handle::box_id = 0;
hid_t HDF5Handle::intvect_id = 0;
hid_t HDF5Handle::realvect_id = 0;
map<std::string, std::string> HDF5Handle::groups = map<std::string, std::string>();

#ifdef H516
extern "C"
{
  herr_t print_and_abort(void* s)
  {
    herr_t e = H5Eprint((FILE*)s);
    abort();
    return e;
  }
}
#else
extern "C"
{
  herr_t print_and_abort(hid_t stackID, void* s)
  {
    herr_t e = H5Eprint2(stackID,(FILE*)s);
    abort();
    return e;
  }
}
#endif

void HDF5Handle::initialize()
{
  H5dont_atexit();

#ifdef H516
  H5Eset_auto(print_and_abort, stderr);
#else
  H5Eset_auto2(H5E_DEFAULT, print_and_abort, stderr);
#endif

  groups["C"] = "Cell";
  groups["X"] = "XFace";
  groups["Y"] = "YFace";
  groups["Z"] = "ZFace";
  groups["N"] = "Node";

  /*  old composite-of-composite version of box
  intvect_id = H5Tcreate (H5T_COMPOUND, sizeof(IntVect));
  H5Tinsert (intvect_id, "intvecti", HOFFSET(IntVect, vect[0]), H5T_NATIVE_INT);
#if CH_SPACEDIM > 1
  H5Tinsert (intvect_id, "intvectj", HOFFSET(IntVect, vect[1]), H5T_NATIVE_INT);
#endif
#if CH_SPACEDIM > 2
  H5Tinsert (intvect_id, "intvectk", HOFFSET(IntVect, vect[2]), H5T_NATIVE_INT);
#endif

  box_id = H5Tcreate (H5T_COMPOUND, sizeof(Box));
  H5Tinsert (box_id, "smallend", HOFFSET(Box, smallend), intvect_id);
  H5Tinsert (box_id, "bigend",   HOFFSET(Box, bigend), intvect_id);
  */

  IntVect d1;
  RealVect d2;
  Box d3;

  intvect_id = H5Tcreate (H5T_COMPOUND, sizeof(IntVect));
  box_id = H5Tcreate (H5T_COMPOUND, sizeof(Box));
  realvect_id = H5Tcreate(H5T_COMPOUND, sizeof(RealVect));
#if CH_SPACEDIM == 1
  H5Tinsert (intvect_id, "intvecti", CHOFFSET(d1, vect), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_i", CHOFFSET(d3, smallend.vect), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_i", CHOFFSET(d3, bigend.vect),   H5T_NATIVE_INT);
  H5Tinsert (realvect_id, "x", CHOFFSET(d2, vect),   H5T_NATIVE_REAL);
#elif CH_SPACEDIM == 2
  H5Tinsert (intvect_id, "intvecti", CHOFFSET(d1, vect), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectj", CHOFFSET(d1, vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_i", CHOFFSET(d3, smallend.vect), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_j", CHOFFSET(d3, smallend.vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_i", CHOFFSET(d3, bigend.vect),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_j", CHOFFSET(d3, bigend.vect)+sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (realvect_id, "x", CHOFFSET(d2, vect),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "y", CHOFFSET(d2, vect)+sizeof(Real),   H5T_NATIVE_REAL);
#elif CH_SPACEDIM == 3
  H5Tinsert (intvect_id, "intvecti", CHOFFSET(d1, vect), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectj", CHOFFSET(d1, vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectk", CHOFFSET(d1, vect)+2*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_i", CHOFFSET(d3, smallend.vect), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_j", CHOFFSET(d3, smallend.vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_k", CHOFFSET(d3, smallend.vect)+2*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_i", CHOFFSET(d3, bigend.vect),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_j", CHOFFSET(d3, bigend.vect)+sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_k", CHOFFSET(d3, bigend.vect)+2*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (realvect_id, "x", CHOFFSET(d2, vect),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "y", CHOFFSET(d2, vect)+sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "z", CHOFFSET(d2, vect)+2*sizeof(Real),   H5T_NATIVE_REAL);
#elif CH_SPACEDIM == 4
  H5Tinsert (intvect_id, "intvecti", CHOFFSET(d1, vect), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectj", CHOFFSET(d1, vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectk", CHOFFSET(d1, vect)+2*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectl", CHOFFSET(d1, vect)+3*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_i", CHOFFSET(d3, smallend.vect), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_j", CHOFFSET(d3, smallend.vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_k", CHOFFSET(d3, smallend.vect)+2*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_l", CHOFFSET(d3, smallend.vect)+3*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_i", CHOFFSET(d3, bigend.vect),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_j", CHOFFSET(d3, bigend.vect)+sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_k", CHOFFSET(d3, bigend.vect)+2*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_l", CHOFFSET(d3, bigend.vect)+3*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (realvect_id, "x", CHOFFSET(d2, vect),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "y", CHOFFSET(d2, vect)+sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "z", CHOFFSET(d2, vect)+2*sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "u", CHOFFSET(d2, vect)+3*sizeof(Real),   H5T_NATIVE_REAL);
#elif CH_SPACEDIM == 5
  H5Tinsert (intvect_id, "intvecti", CHOFFSET(d1, vect), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectj", CHOFFSET(d1, vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectk", CHOFFSET(d1, vect)+2*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectl", CHOFFSET(d1, vect)+3*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectm", CHOFFSET(d1, vect)+4*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_i", CHOFFSET(d3, smallend.vect), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_j", CHOFFSET(d3, smallend.vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_k", CHOFFSET(d3, smallend.vect)+2*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_l", CHOFFSET(d3, smallend.vect)+3*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_m", CHOFFSET(d3, smallend.vect)+4*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_i", CHOFFSET(d3, bigend.vect),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_j", CHOFFSET(d3, bigend.vect)+sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_k", CHOFFSET(d3, bigend.vect)+2*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_l", CHOFFSET(d3, bigend.vect)+3*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_m", CHOFFSET(d3, bigend.vect)+4*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (realvect_id, "x", CHOFFSET(d2, vect),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "y", CHOFFSET(d2, vect)+sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "z", CHOFFSET(d2, vect)+2*sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "u", CHOFFSET(d2, vect)+3*sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "v", CHOFFSET(d2, vect)+4*sizeof(Real),   H5T_NATIVE_REAL);
#elif CH_SPACEDIM == 6
  H5Tinsert (intvect_id, "intvecti", CHOFFSET(d1, vect), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectj", CHOFFSET(d1, vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectk", CHOFFSET(d1, vect)+2*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectl", CHOFFSET(d1, vect)+3*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectm", CHOFFSET(d1, vect)+4*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (intvect_id, "intvectn", CHOFFSET(d1, vect)+5*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_i", CHOFFSET(d3, smallend.vect), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_j", CHOFFSET(d3, smallend.vect)+sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_k", CHOFFSET(d3, smallend.vect)+2*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_l", CHOFFSET(d3, smallend.vect)+3*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_m", CHOFFSET(d3, smallend.vect)+4*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "lo_n", CHOFFSET(d3, smallend.vect)+5*sizeof(int), H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_i", CHOFFSET(d3, bigend.vect),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_j", CHOFFSET(d3, bigend.vect)+sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_k", CHOFFSET(d3, bigend.vect)+2*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_l", CHOFFSET(d3, bigend.vect)+3*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_m", CHOFFSET(d3, bigend.vect)+4*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (box_id, "hi_n", CHOFFSET(d3, bigend.vect)+5*sizeof(int),   H5T_NATIVE_INT);
  H5Tinsert (realvect_id, "x", CHOFFSET(d2, vect),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "y", CHOFFSET(d2, vect)+sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "z", CHOFFSET(d2, vect)+2*sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "u", CHOFFSET(d2, vect)+3*sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "v", CHOFFSET(d2, vect)+4*sizeof(Real),   H5T_NATIVE_REAL);
  H5Tinsert (realvect_id, "w", CHOFFSET(d2, vect)+5*sizeof(Real),   H5T_NATIVE_REAL);
#else
  undefined_dimension!
#endif

  initialized = true;
}

HDF5Handle::HDF5Handle(): m_isOpen(false)
{
  if (!initialized) initialize();
}

HDF5Handle::HDF5Handle(
        const std::string& a_filename,
        mode a_mode,
        const char *a_globalGroupName)
            :
        m_isOpen(false),
        m_level(-1)
{
  int err = open(a_filename, a_mode, a_globalGroupName);
  if (err < 0 )
  {
    char buf[1024];
    sprintf(buf, "Problem opening file %s ",a_filename.c_str());
    MayDay::Error(buf);
  }
}

HDF5Handle::~HDF5Handle()
{
  CH_assert( !m_isOpen );
}

//#ifdef CH_MPI
//#if ( H5_VERS_MAJOR == 1 && ( H5_VERS_MINOR < 4 || ( H5_VERS_MINOR == 4 && H5_VERS_RELEASE <= 1 ) ) )
//// there's a bug in the HDF include file 'H5FDmpio.h': no extern "C"
//extern "C" {
//herr_t H5Pset_fapl_mpio(hid_t fapl_id, MPI_Comm comm, MPI_Info info);
//};
//#endif
//#endif

int HDF5Handle::open(
        const std::string& a_filename,
        mode a_mode,
        const char *a_globalGroupName)
{
  int ret = 0;
  if (m_isOpen)
    {
      MayDay::Error("Calling 'open' on already open file.  use 'close' on finished files");
    }

  m_filename = a_filename;
  if (!initialized) initialize();
  m_group    = "/";

  hid_t file_access = 0;
  if (a_mode != CREATE_SERIAL)
    {
#ifdef CH_MPI
      file_access = H5Pcreate (H5P_FILE_ACCESS);

#if ( H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 2 )
      H5Pset_mpi(file_access,  Chombo_MPI::comm, MPI_INFO_NULL);
#else
      H5Pset_fapl_mpio(file_access,  Chombo_MPI::comm, MPI_INFO_NULL);
#endif
#else
      file_access = H5P_DEFAULT;
#endif
    }

  switch(a_mode)
  {
  case CREATE:
    m_fileID = H5Fcreate(a_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_access);
    if (m_fileID < 0) return m_fileID;
    break;
  case CREATE_SERIAL:
    m_fileID = H5Fcreate(a_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (m_fileID < 0) return m_fileID;
    break;
  case OPEN_RDONLY:
    m_fileID = H5Fopen(a_filename.c_str(), H5F_ACC_RDONLY, file_access);
    if (m_fileID < 0) return m_fileID;
    break;
  case OPEN_RDWR:
    m_fileID = H5Fopen(a_filename.c_str(), H5F_ACC_RDWR, file_access);
    if (m_fileID < 0) return m_fileID;
    break;
  default:
    MayDay::Error("unrecognized file access mode:  HDF5Handle::open");
  }
  H5Pclose(file_access);
#ifdef H516
  m_currentGroupID = H5Gopen(m_fileID, m_group.c_str());
#else
  m_currentGroupID = H5Gopen2(m_fileID, m_group.c_str(),H5P_DEFAULT);
#endif
  if (m_fileID >= 0 && m_currentGroupID >= 0) m_isOpen = true;

  // Write or read dimension checks and accuracy data

  HDF5HeaderData info;
  char buf[10000];
  hid_t attr, datatype, group;
  size_t codeprecision, fileprecision;
  switch(a_mode)
  {
  case CREATE_SERIAL:
  case CREATE:
#ifdef H516
    group = H5Gcreate(m_fileID, a_globalGroupName, 0);
#else
    group = H5Gcreate2(m_fileID, a_globalGroupName, H5P_DEFAULT,
		       H5P_DEFAULT,H5P_DEFAULT);
#endif
    info.m_int["SpaceDim"] = SpaceDim;
    info.m_real["testReal"] = 0.0;
    info.writeToLocation(group);
    break;
  default:
#ifdef H516
    group =  H5Gopen(m_fileID, a_globalGroupName);
#else
    group =  H5Gopen2(m_fileID, a_globalGroupName, H5P_DEFAULT);
#endif
    if (group < 0)
      {
        MayDay::Warning("This files appears to be missing a 'Chombo_global' section");
      }
    info.readFromLocation(group);
    if (info.m_int["SpaceDim"] == 0)
      {
        sprintf(buf,"File '%s' appears to lack a SpaceDim definition",a_filename.c_str());
        MayDay::Error(buf);
      }
    if (info.m_int["SpaceDim"] != SpaceDim)
      {
        sprintf(buf,"SpaceDim of '%s' does not match code",a_filename.c_str());
        MayDay::Error(buf);
      }

#ifdef H516
    attr = H5Aopen_name(group, "testReal");
#else
    attr = H5Aopen_by_name(group, ".", "testReal", H5P_DEFAULT,H5P_DEFAULT);
#endif
    if (attr < 0) return 1;
    datatype = H5Aget_type(attr);
    fileprecision = H5Tget_precision(datatype);
    codeprecision = H5Tget_precision(H5T_NATIVE_REAL);
    if (fileprecision > codeprecision) ret = 2;
    if (codeprecision > fileprecision)
      {
        sprintf(buf, "code is compiled with Real=%i bits, file %s has Real=%i bits",
                (unsigned int)codeprecision, a_filename.c_str(), (unsigned int)fileprecision);
        MayDay::Warning(buf);
        ret = 2;
      }
    H5Aclose(attr);
    H5Tclose(datatype);
  }
  H5Gclose(group);

  return ret;
}

bool HDF5Handle::isOpen() const
{
  return m_isOpen;
}

void HDF5Handle::setGroupToLevel(int a_level)
{
  char ch[100];
  sprintf(ch, "level_%i", a_level);
  setGroup(ch);
}

int HDF5Handle::setGroup(const std::string& group)
{
  if (!m_isOpen)
    {
      MayDay::Error("cannot access group until file is successfully opened");
    }
  if (m_group == group) return 0;

  int ret = 0;

  ret = H5Gclose(m_currentGroupID);
  if (ret < 0)
    {
      MayDay::Warning(" Error closing old group");
      return ret;
    }

  H5E_auto_t efunc; void* edata; // turn auto error messaging off
#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  m_currentGroupID = H5Gopen(m_fileID, group.c_str());
  if (m_currentGroupID < 0)
    {
      H5Eset_auto(efunc, edata); //turn error messaging back on.
      //open failed, go to group creation
      m_currentGroupID = H5Gcreate(m_fileID, group.c_str(), 0);
    }
  if (m_currentGroupID < 0)
    ret = -1;

  H5Eset_auto(efunc, edata); //turn error messaging back on.
#else
  H5Eget_auto(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  m_currentGroupID = H5Gopen2(m_fileID, group.c_str(), H5P_DEFAULT);
  if (m_currentGroupID < 0)
    {
      H5Eset_auto2(H5E_DEFAULT, efunc, edata); //turn error messaging back on.
      //open failed, go to group creation
      m_currentGroupID = H5Gcreate2(m_fileID, group.c_str(),H5P_DEFAULT, 
				    H5P_DEFAULT,H5P_DEFAULT);
    }
  if (m_currentGroupID < 0)
    ret = -1;

  H5Eset_auto2(H5E_DEFAULT, efunc, edata); //turn error messaging back on.
#endif
  m_group = group;
  return ret;
}

int HDF5Handle::pushGroup(const std::string& grp)
{
  if ( getGroup() == "/" )
  {
    return setGroup( getGroup()                    + grp );
  } else
  {
    return setGroup( getGroup() + std::string("/") + grp );
  }
}

int HDF5Handle::popGroup()
{
  std::string grp( getGroup() );

  // Check for weirdness (probably impossible):
  CH_assert( grp.find( "//" ) == grp.npos );

  std::string::size_type i = grp.rfind('/');
  CH_assert( i != grp.npos );
  CH_assert( grp != "/" ); // Error to call popGroup() if group is "/".

  if ( i == (grp.size()-1) )    // Trailing '/' ... possible though not likely.
  {
    i = grp.rfind('/', i-1);
    CH_assert( i != grp.npos );
    CH_assert( i != 0 );
  }

  if ( i == 0 )  // We're popping down to the root group.
  {
    grp.erase( i+1, grp.npos );
  } else
  {
    grp.erase( i, grp.npos );
  }
  return setGroup( grp );
}

void HDF5Handle::close()
{
  CH_assert( m_isOpen );
  if (m_isOpen)
  {
    if (!(m_currentGroupID < 0)) H5Gclose(m_currentGroupID);
    if (!(m_fileID < 0)) H5Fclose(m_fileID);
    m_isOpen = false;
  }
}

const std::string&  HDF5Handle::getGroup() const
{
  return m_group;
}

const hid_t& HDF5Handle::fileID() const
{
  return m_fileID;
}

const hid_t& HDF5Handle::groupID() const
{
  return m_currentGroupID;
}

//=====================================================================================
  /// writes the current attribute list to the current group in 'file'
int HDF5HeaderData::writeToFile(HDF5Handle& file) const
{
  return writeToLocation(file.groupID());
}

int HDF5HeaderData::writeToLocation(hid_t loc_id) const
{
  H5E_auto_t efunc; void* edata;
#ifdef H516
  H5Eget_auto(&efunc, &edata);
#else
  H5Eget_auto2(H5E_DEFAULT, &efunc, &edata);
#endif
  herr_t  ret;
  char messg[1024];
#define INSERT(Ttype, mapName, H5Ttype)                                   \
  for (map<std::string, Ttype>::const_iterator p = mapName.begin();        \
      p!= mapName.end(); ++p)                                             \
    {                                                                     \
      hid_t aid  = H5Screate(H5S_SCALAR);                                 \
      H5Eset_auto(NULL, NULL);                                            \
      hid_t attr = H5Acreate(loc_id, p->first.c_str(), H5Ttype,           \
                             aid, H5P_DEFAULT);                           \
      if (attr < 0)                                                        \
        {                                                                 \
          H5Adelete(loc_id, p->first.c_str());                            \
          attr = H5Acreate(loc_id, p->first.c_str(), H5Ttype,             \
                           aid, H5P_DEFAULT);                             \
          if (attr < 0)                                                    \
            {                                                             \
              sprintf(messg," Problem writing attribute %s",p->first.c_str());  \
              MayDay::Warning(messg);                                     \
            }                                                             \
        }                                                                 \
      H5Eset_auto(efunc, edata);                                          \
      Ttype tmp = p->second;                                              \
      ret = H5Awrite(attr, H5Ttype, &tmp);                                \
      if (ret < 0) return ret;                                             \
      H5Sclose(aid);                                                      \
      H5Aclose(attr);                                                     \
    }                                                                     \

#define INSERT2(Ttype, mapName, H5Ttype)                                   \
  for (map<std::string, Ttype>::const_iterator p = mapName.begin();        \
      p!= mapName.end(); ++p)                                             \
    {                                                                     \
      hid_t aid  = H5Screate(H5S_SCALAR);                                 \
      H5Eset_auto2(H5E_DEFAULT, NULL, NULL);				\
      hid_t attr = H5Acreate2(loc_id, p->first.c_str(), H5Ttype,           \
			      aid, H5P_DEFAULT, H5P_DEFAULT);			\
      if (attr < 0)                                                        \
        {                                                                 \
          H5Adelete(loc_id, p->first.c_str());                            \
          attr = H5Acreate2(loc_id, p->first.c_str(), H5Ttype,             \
			    aid, H5P_DEFAULT, H5P_DEFAULT);			          if (attr < 0)                                                    \
            {                                                             \
              sprintf(messg," Problem writing attribute %s",p->first.c_str());  \
              MayDay::Warning(messg);                                     \
            }                                                             \
        }                                                                 \
      H5Eset_auto2(H5E_DEFAULT, efunc, edata);				\
      Ttype tmp = p->second;                                              \
      ret = H5Awrite(attr, H5Ttype, &tmp);                                \
      if (ret < 0) return ret;                                             \
      H5Sclose(aid);                                                      \
      H5Aclose(attr);                                                     \
    }                                                                     \

#ifdef H516
    INSERT(Real, m_real, H5T_NATIVE_REAL);
    INSERT(int, m_int, H5T_NATIVE_INT);
    INSERT(IntVect, m_intvect, HDF5Handle::intvect_id);
    INSERT(Box, m_box, HDF5Handle::box_id);
    INSERT(RealVect, m_realvect, HDF5Handle::realvect_id);
#else
    INSERT2(Real, m_real, H5T_NATIVE_REAL);
    INSERT2(int, m_int, H5T_NATIVE_INT);
    INSERT2(IntVect, m_intvect, HDF5Handle::intvect_id);
    INSERT2(Box, m_box, HDF5Handle::box_id);
    INSERT2(RealVect, m_realvect, HDF5Handle::realvect_id);
#endif
    // string is different, of course

    for (map<std::string, std::string>::const_iterator p = m_string.begin();
        p!= m_string.end(); ++p)
    {
      hid_t s_type = H5Tcopy(H5T_C_S1);
      H5Tset_size(s_type, p->second.length()); //extra requirement for strings
      hid_t aid  = H5Screate(H5S_SCALAR);
#ifdef H516
      H5Eset_auto(NULL, NULL);
      hid_t attr = H5Acreate(loc_id, p->first.c_str(), s_type,
                             aid, H5P_DEFAULT);
#else
      H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
      hid_t attr = H5Acreate2(loc_id, p->first.c_str(), s_type,
			      aid, H5P_DEFAULT, H5P_DEFAULT);
#endif
      if (attr < 0)
        {
          H5Adelete(loc_id, p->first.c_str());
#ifdef H516
          attr = H5Acreate(loc_id, p->first.c_str(), s_type,
                           aid, H5P_DEFAULT);
#else
          attr = H5Acreate2(loc_id, p->first.c_str(), s_type,
			    aid, H5P_DEFAULT,H5P_DEFAULT);
#endif
          if (attr < 0)
            {
              sprintf(messg," Problem writing attribute %s",p->first.c_str());
              MayDay::Warning(messg);
            }
        }
#ifdef H516
      H5Eset_auto(efunc, edata);
#else
      H5Eset_auto2(H5E_DEFAULT, efunc, edata);
#endif
      char* tmp = (char*)p->second.c_str();
      ret = H5Awrite(attr, s_type, tmp);
      if (ret < 0) return ret;
      H5Sclose(aid);
      H5Aclose(attr);
      H5Tclose(s_type);
    }

    return 0;
}

  /// read process is add/change, does not remove key-value pairs. reads from current group.
int HDF5HeaderData::readFromFile(HDF5Handle& file)
{
#ifdef H516
  return H5Aiterate(file.groupID(), NULL, HDF5HeaderDataattributeScan , this);
#else
  hsize_t* n=0;
  return H5Aiterate2(file.groupID(), H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, n, HDF5HeaderDataattributeScan , this);
#endif

}

int HDF5HeaderData::readFromLocation(hid_t loc_id)
{
#ifdef H516
   return H5Aiterate(loc_id, NULL, HDF5HeaderDataattributeScan , this);
#else
  hsize_t* n=0;
  return H5Aiterate2(loc_id, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, n, HDF5HeaderDataattributeScan , this);
#endif
}

extern "C"
{
#ifdef H516
  herr_t HDF5HeaderDataattributeScan(hid_t loc_id, const char *name, void *opdata)
#else
  herr_t HDF5HeaderDataattributeScan(hid_t loc_id, const char *name, const H5A_info_t* info, void *opdata)
#endif
{
  herr_t ret = 0;
  HDF5HeaderData& data = *(static_cast<HDF5HeaderData*>(opdata));

#ifdef H516
  hid_t attr   = H5Aopen_name(loc_id, name);
#else
  hid_t attr   = H5Aopen_by_name(loc_id,".", name, H5P_DEFAULT, H5P_DEFAULT );
#endif
  hid_t atype  = H5Aget_type(attr);
  hid_t aclass = H5Tget_class(atype);
  char* buf = NULL;  size_t size = 0;

  switch(aclass)
  {
  case H5T_INTEGER :
    int Ivalue;
    ret  = H5Aread(attr, H5T_NATIVE_INT, &Ivalue);
    if (ret < 0) break;
    data.m_int[name] = Ivalue;
    break;
  case H5T_FLOAT:
    Real Rvalue;
    ret = H5Aread(attr, H5T_NATIVE_REAL, &Rvalue);
    if (ret < 0) break;
    data.m_real[name] = Rvalue;
    break;
  case H5T_STRING:
    size = H5Tget_size(atype);
    buf = new char[size+1];
    ret = H5Aread(attr, atype, buf);
    if (ret < 0) break;
    buf[size] = 0; // for some reason HDF5 is not null terminating strings correctly
    data.m_string[name] = std::string(buf);
    break;
  case H5T_COMPOUND:
    if (strcmp(H5Tget_member_name(atype, 0), "lo_i") == 0)
      {
        Box value;
        ret = H5Aread(attr, HDF5Handle::box_id, &value);
        if (ret < 0) break;
        value.computeBoxLen();
        data.m_box[name] = value;
        break;
      }
    else if (strcmp(H5Tget_member_name(atype, 0), "intvecti") == 0)
      {
        IntVect value;
        ret = H5Aread(attr, HDF5Handle::intvect_id, &value);
        if (ret < 0) break;
        data.m_intvect[name] = value;
        break;
      }
    else if (strcmp(H5Tget_member_name(atype, 0), "x") == 0)
      {
        RealVect value;
        ret = H5Aread(attr, HDF5Handle::realvect_id, &value);
        if (ret < 0) break;
        data.m_realvect[name] = value;
        break;
      }
  default:
    // heel if I know what to do about attributes I don't recognize
    MayDay::Warning("HDF5HeaderData::readFromFile encountered unrecognized attribute");
  }
  delete[] buf;
  H5Tclose(atype);
  H5Aclose(attr);
  return ret;
}
          }
void HDF5HeaderData::clear()
{
  m_real.clear();
  m_int.clear();
  m_string.clear();
  m_intvect.clear();
  m_realvect.clear();
  m_box.clear();
}

template< >
hid_t H5Type(const int* dummy)
{
  return H5T_NATIVE_INT;
}

template< >
hid_t H5Type(const float* dummy)
{
  return H5T_NATIVE_FLOAT;

}

template< >
hid_t H5Type(const double* dummy)
{
  return H5T_NATIVE_DOUBLE;
}

template< >
hid_t H5Type(const Box* dummy)
{
  return HDF5Handle::box_id;
}

template< >
hid_t H5Type(const RealVect* dummy)
{
  return  HDF5Handle::realvect_id;
}

template< >
hid_t H5Type(const IntVect* dummy)
{
  return HDF5Handle::intvect_id;
}

template< >
hid_t H5Type(const long long* dummy)
{
  return H5T_NATIVE_LLONG;
}

void writeDataset(hid_t a_dataset,
                  hid_t a_dataspace,
                  const void* start,
                  ch_offset_t off,
                  hsize_t  count)
{
  ch_offset_t offset[1]   ; offset[0] = off;
  hsize_t  flatdims[1] ; flatdims[0] = count;

  hid_t memdataspace = H5Screate_simple(1, flatdims, NULL);
  H5Sselect_hyperslab (a_dataspace, H5S_SELECT_SET, offset, NULL,
                       flatdims, NULL);
  H5Dwrite(a_dataset, H5Dget_type(a_dataset), memdataspace, a_dataspace,
           H5P_DEFAULT, start);
  H5Sclose(memdataspace);
}

void readDataset(hid_t a_dataset,
                 hid_t a_dataspace,
                 void* start,
                 ch_offset_t off,
                 hsize_t  count)
{
  ch_offset_t offset[1]   ; offset[0] = off;
  hsize_t  flatdims[1] ; flatdims[0] = count;

  hid_t memdataspace = H5Screate_simple(1, flatdims, NULL);
  H5Sselect_hyperslab (a_dataspace, H5S_SELECT_SET, offset, NULL,
                       flatdims, NULL);
  H5Dread(a_dataset, H5Dget_type(a_dataset), memdataspace, a_dataspace,
           H5P_DEFAULT, start);
  H5Sclose(memdataspace);
}

void createData(hid_t& a_dataset,
                hid_t& a_dataspace,
                HDF5Handle& handle,
                const std::string& name,
                hid_t type,
                hsize_t size)
{

  hsize_t  flatdims[1];
  flatdims[0] = size;
  H5E_auto_t efunc; void* edata; // turn auto error messaging off
#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  H5Gunlink(handle.groupID(), name.c_str()); //removes a pre-existing dataset.
  H5Eset_auto(efunc, edata);
  a_dataspace = H5Screate_simple(1, flatdims, NULL);
  a_dataset   = H5Dcreate(handle.groupID(), name.c_str(),  type,
                          a_dataspace, H5P_DEFAULT);
#else
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  H5Ldelete(handle.groupID(), name.c_str(), H5P_DEFAULT); //removes a pre-existing dataset.
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);
  a_dataspace = H5Screate_simple(1, flatdims, NULL);
  a_dataset   = H5Dcreate2(handle.groupID(), name.c_str(),  type,
			   a_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

}

ostream& operator<<(ostream& os, const HDF5HeaderData& data)
{

#define PRINT(Ttype, mapName) \
  for (map<std::string, Ttype>::const_iterator p = data.mapName.begin();        \
      p!= data.mapName.end(); ++p)                   \
  os<<p->first<<" :"<<p->second<<"\n";

  PRINT(Real, m_real);
  PRINT(int, m_int);
  PRINT(std::string, m_string);
  PRINT(IntVect, m_intvect);
  PRINT(Box, m_box);
  return os;
}

void HDF5HeaderData::dump() const
{
  std::cout<<*this<<std::endl;
}

#include "NamespaceFooter.H"
#endif // CH_USE_HDF5
