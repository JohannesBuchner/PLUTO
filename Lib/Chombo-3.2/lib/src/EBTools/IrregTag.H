#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _IRREGTAG_H_
#define _IRREGTAG_H_

#include "Vector.H"
#include "VolIndex.H"
#include "IntVectSet.H"
#include <list>
#include "SPMD.H"
#include "NamespaceHeader.H"

class HDF5Handle;
class IrregTag;

/**
   class used as return type of IrregTag access functions to
   encapsulate the implementation.  Also provide a lightweight
   return object.
*/

class TagSet
{
public:
  TagSet()
    :
    index(0),
    list(NULL),
    value(NULL)
  {
  }

  TagSet(const Vector<int>* l, const Vector<double>* d);

  void begin()
  {
    index=0;
  }

  bool ok() const
  {
    return index < list->size();
  }

  void operator++()
  {
    index++;
  }

  int operator()() const
  {
    return (*list)[index];
  }

  double val() const
  {
    return (*value)[index];
  }

private:
  int index;
  const Vector<int>* list;
  const Vector<double>* value;
  friend class IrregTag;
};

class IrregTag
{
public:

  /**
   */
  IrregTag();

  /**
   */
  ~IrregTag();

  /*==============================*/

 /**
     /name I/O routines
  */
  /*@{*/

  void read(const std::string& a_filename);
  void read(HDF5Handle& handle);

  /** IrregTag class must be closed before writing
      @see close
  */
  void write(const std::string& a_filename);

  /** IrregTag class must be closed before writing
      @see close
  */
  void write(HDF5Handle& handle);
  /*@}*/

  /*==============================*/

  /**
     /name Construction routines
  */
  /*@{*/

  void setTags(const VolIndex& a_vol, const Vector<int>& tags,
               const Vector<double>& values);

  /**  Called when no more calls to setTags will be performed.
   */
  void close();

  /*@}*/

  /*==============================*/

  /**
     /name Access routines

     IrregTag object must be 'closed' before these operations can
     be performed.
     @see close()
  */
  /*@{*/
  TagSet tags(const VolIndex& a_index) const ;

  int numVol() const ;

  TagSet tags(int a_tagsetIndex) const;

  /*@}*/

private:

  class Entry
  {
  public:
    VolIndex       m_index;
    Vector<int>    m_tags;
    Vector<double> m_values;

    Entry()
    {
    }

    Entry(const VolIndex& index,
          const Vector<int>& tags,
          const Vector<double>& values)
      :
      m_index(index),
      m_tags(tags),
      m_values(values)
    {
    }

    void define(const VolIndex& index,
                const Vector<int>& tags,
                const Vector<double>& values)
    {
      m_index = index;
      m_tags = tags;
      m_values = values;
    }

    bool operator < (const Entry& rhs) const
    {
      return m_index < rhs.m_index;
    }

    int linearSize() const
    {
      return m_index.linearSize() + CH_XD::linearListSize(m_tags) + CH_XD::linearListSize(m_values);
    }

    void linearOut(void* const a_outBuf) const
    {
      char* buf = (char*)a_outBuf;
      buf+=m_index.linearSize();
      m_index.linearOut(a_outBuf);
      CH_XD::linearOut(buf, m_tags);
      buf += CH_XD::linearSize(m_tags);
      CH_XD::linearOut(buf, m_values);
    }

    void linearIn(void* const inBuf)
    {
      const char* buf = (const char*)inBuf;
      buf+=m_index.linearSize();
      m_index.linearIn(inBuf);
      CH_XD::linearIn(m_tags, buf);
      buf += CH_XD::linearSize(m_tags);
      CH_XD::linearIn(m_values, buf);

    }
  };

  bool m_closed;
  IntVectSet m_ivs;
  std::list<Entry> m_accumulator;
  Vector<Entry> m_set;
};

#include "NamespaceFooter.H"
#endif
