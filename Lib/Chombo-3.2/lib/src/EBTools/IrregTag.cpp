#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IrregTag.H"
#include "CH_HDF5.H"
#include "NamespaceHeader.H"

TagSet::TagSet(const Vector<int>* l, const Vector<double>* d)
  :
  index(0),
  list(l),
  value(d)
{
}

IrregTag::IrregTag()
  :
  m_closed(false)
{
}

IrregTag::~IrregTag()
{
}

void IrregTag::setTags(const VolIndex& a_vol, const Vector<int>& tags,
                       const Vector<double>& vals)
{
  CH_assert(!m_closed);
  Entry next( a_vol, tags, vals);
  m_accumulator.push_back(next);
}

void IrregTag::close()
{
  m_accumulator.sort();
  m_set.resize(m_accumulator.size());
  std::list<Entry>::const_iterator a_it = m_accumulator.begin();
  for (int i=0; i<m_set.size(); i++, a_it++)
  {
    m_set[i] = *a_it;
    m_ivs|=(*a_it).m_index.gridIndex();
  }
  /*
#ifdef CH_MPI
  Vector<Vector<Entry> > agg;
  gather(agg,m_set, 0);
  if (procID() == 0)
    {
      int size=agg[0].size();
      for (int p=1; p<agg.size(); p++) size+=agg[p].size();
      m_set.resize(size);
      int index=0;
      for (int p=0; p<m_set.size(); p++)
        {
          const Vector<Entry>& v = agg[p];
          for (int i=0; i<v.size(); i++, index++)
            m_set[index] = v[i];
        }
    }
  broadcast(m_set, 0);

  for (int i=0; i<m_set.size(); i++) m_ivs|=(m_set[i].m_index.gridIndex());

#endif
  */
  m_closed = true;
  m_accumulator.clear();
}

TagSet IrregTag::tags(const VolIndex& a_index) const
{
  static Vector<int> emptyVector;
  static Vector<double> emptyDouble;
  static TagSet emptySet(&emptyVector, &emptyDouble);
  if (!m_ivs.contains(a_index.gridIndex()))
    return emptySet;

  TagSet rtn;

  //binary search in sorted list
  int first = 0;
  int last  = m_set.size()-1;
  int index =0;

  while (first <= last)
    {
      index = (first+last)/2;
      const Entry& entry = m_set[index];

      if (a_index == entry.m_index) break;
      if (a_index < entry.m_index)
        {
          last = index-1;
        }
      else
        {
          first = index + 1;
        }
    }

  const Entry& entry = m_set[index];
  CH_assert(a_index == entry.m_index);
  rtn.list = &(entry.m_tags);
  rtn.value = &(entry.m_values);

  return rtn;
}

#ifdef CH_USE_HDF5
void IrregTag::write(const std::string& a_filename)
{
  CH_assert(m_closed);
  HDF5Handle handle(a_filename, HDF5Handle::CREATE);
  write(handle);
}

void IrregTag::write(HDF5Handle& a_handle)
{
  CH_assert(m_closed);
  CH_assert(a_handle.isOpen());

  if (procID() == 0)
    {
      HDF5HeaderData att;
      att.m_int["num_irreg"] = m_set.size();
      att.writeToFile(a_handle);
      int size = 0;
      for (int i=0; i<m_set.size(); i++)
      {
        size += m_set[i].linearSize();
      }

      char* buffer = (char*)malloc(size);
      if (buffer == NULL)
      {
        MayDay::Error("Error in memory allocation in IrregTag::write");
      }
      char* buf=buffer;
      for (int i=0; i<m_set.size(); i++)
        {
          m_set[i].linearOut(buf);
          buf+=m_set[i].linearSize();
        }
      hsize_t  flatdims[1];

      flatdims[0] = size/sizeof(int);
      hid_t dataspace = H5Screate_simple(1, flatdims, NULL);
#ifdef H516
      hid_t dataset   = H5Dcreate(a_handle.groupID(), "IrregTags",
                                  H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
#else
      hid_t dataset   = H5Dcreate2(a_handle.groupID(), "IrregTags",  
				   H5T_NATIVE_INT,
				   dataspace, H5P_DEFAULT,
				   H5P_DEFAULT,H5P_DEFAULT);

#endif
      hid_t memdataspace = H5Screate_simple(1, flatdims, NULL);
      H5Dwrite(dataset, H5T_NATIVE_INT , memdataspace, dataspace,
               H5P_DEFAULT, buffer);

      H5Sclose(dataspace);
      H5Sclose(memdataspace);
      H5Dclose(dataset);
      free(buffer);
    }
}

void IrregTag::read(const std::string& a_filename)
{
  HDF5Handle handle(a_filename, HDF5Handle::OPEN_RDONLY);
  read(handle);

}
void IrregTag::read(HDF5Handle& a_handle)
{

#ifdef H516
  hid_t  dataset = H5Dopen(a_handle.groupID(), "IrregTags");
#else
  hid_t  dataset = H5Dopen2(a_handle.groupID(), "IrregTags",H5P_DEFAULT);
#endif
  CH_assert(dataset > 0 );
  hid_t dataspace =  H5Dget_space(dataset);
  CH_assert(dataspace > 0);

  hsize_t dims[1], maxdims[1];
  H5Sget_simple_extent_dims(dataspace, dims, maxdims);

  hid_t memdataspace = H5Screate_simple(1, dims, NULL);

  char* buffer = (char*)malloc(dims[0]*sizeof(int));
  if (buffer == NULL)
    MayDay::Error("Memory allocation error in  IrregTags::read");

  herr_t error = H5Dread(dataset, H5T_NATIVE_INT, memdataspace, dataspace,
                         H5P_DEFAULT, buffer);
  CH_assert(error >= 0);

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(memdataspace);

  m_ivs.makeEmpty();

  HDF5HeaderData att;
  att.readFromFile(a_handle);
  int size = att.m_int["num_irreg"];
  m_set.resize(size);
  char* buf = buffer;
  for (int i=0; i<size; i++)
    {
      Entry& entry = m_set[i];
      entry.linearIn(buf);
      buf += entry.linearSize();
      m_ivs|=entry.m_index.gridIndex();
    }

  free(buffer);
  m_closed = true;
  m_accumulator.clear();
}
#endif // CH_USE_HDF5

int IrregTag::numVol() const
{
  if (m_closed) return m_set.size();
  return m_accumulator.size();
}

TagSet IrregTag::tags(int a_tagsetIndex) const
{
  const Entry* entry=NULL;
  if (m_closed)
    {
      entry = &(m_set[a_tagsetIndex]);
    }
  else
    {
      MayDay::Error("IrregTag::tags(int) called on open IrregTag");
    }
  TagSet rtn;
  rtn.list = &(entry->m_tags);
  rtn.value = &(entry->m_values);

  return rtn;
}
#include "NamespaceFooter.H"
