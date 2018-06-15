#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IFIDTABLE_H_
#define _IFIDTABLE_H_

#include "IFid.H"
#include <algorithm>
#include "NamespaceHeader.H"

//! \class IFidTable
//! This class maintains a mapping between integer indices and
//! IFids for a given EB index space.
class IFidTable
{
  public:

  //! Creates an empty table.
  IFidTable()
  {
  }

  //! Destructor (not virtual).
  ~IFidTable()
  {
  }

  //! Adds the given IFid to the table, returning the index at which
  //! it may be subsequently retrieved. If the IFid already exists in
  //! the table, the table is unaltered and the index of the existing
  //! entry is returned.
  int add(const IFid& a_id)
  {
    std::vector<IFid>::const_iterator iter =
      find(m_ids.stdVector().begin(), m_ids.stdVector().end(), a_id);
    if (iter == m_ids.stdVector().end()) // Not found!
    {
      int index = m_ids.size();
      m_ids.push_back(a_id);
      return index;
    }
    else
      return iter - m_ids.stdVector().begin(); // Return the existing index.
  }

  //! Retrieves the identifier corresponding to the given index.
  const IFid& operator[](int a_index) const
  {
    CH_assert(a_index >= 0);
    CH_assert(a_index < m_ids.size());
    return m_ids[a_index];
  }

  //! Returns the number of IDs in the table.
  int size() const
  {
    return m_ids.size();
  }

  // Serialization helpers.
  int linearSize() const
  {
    return linearListSize(m_ids);
  }
  void linearIn(const void* const inBuf)
  {
    linearListIn(m_ids, inBuf);
  }
  void linearOut(void* const a_outBuf) const
  {
    linearListOut(a_outBuf, m_ids);
  }

  private:

  // List of IDs.
  Vector<IFid> m_ids;

  // Forbidden.
  IFidTable(const IFidTable&);
  IFidTable& operator=(const IFidTable&);
};

//-----------------------------------------------------------------------
// Template specializations of serialization functions.
//-----------------------------------------------------------------------
template <>
inline int
linearSize(const IFidTable& a_table)
{
  return a_table.linearSize();
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <>
inline void
linearIn(IFidTable& a_table,
         const void* const inBuf)
{
  a_table.linearIn(inBuf);
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
template <>
inline void
linearOut(void* const a_outBuf,
          const IFidTable& a_table)
{
  a_table.linearOut(a_outBuf);
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

#endif
