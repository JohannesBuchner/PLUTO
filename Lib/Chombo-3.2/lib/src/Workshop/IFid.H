#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IFID_H_
#define _IFID_H_

#include "Vector.H"
#include <string>
#include "CH_assert.H"
#include "SPMD.H"
#include "NamespaceHeader.H"

//! \class IFid
//! This class uniquely identifies a boundary represented by a primitive
//! implicit function. We don't envision this being subclassed.
class IFid
{
  public:

  //! Constructs an empty ID. This is needed for serialization.
  IFid()
    :
    m_names()
  {
  }

  //! Constructs an ID from the given string.
  explicit IFid(const std::string& a_name)
    :
    m_names(1, a_name)
  {
  }

  //! Copy constructor.
  IFid(const IFid& a_rhs)
    :
    m_names(a_rhs.m_names)
  {
  }

  //! This factory method creates an ID by aggregating the given list of IDs.
  static IFid aggregate(const Vector<IFid>& a_IDs)
  {
    CH_assert(a_IDs.size() > 0);
    IFid agg(a_IDs[0]);
    for (int i = 1; i < a_IDs.size(); ++i)
    {
      agg.append(a_IDs[i]);
    }
    return agg;
  }

  //! Destructor (not virtual).
  ~IFid()
  {
  }

  //! Assignment operator.
  IFid& operator=(const IFid& a_rhs)
  {
    if (this != &a_rhs)
    {
      m_names = a_rhs.m_names;
    }
    return *this;
  }

  //! Sets the names associated with this identifier.
  void setNames(const Vector<std::string>& a_names)
  {
    m_names = a_names;
  }

  //! Returns the list of names identifying all functions associated with this
  //! IFid.
  const Vector<std::string>& names() const
  {
    return m_names;
  }

  //! Returns true if the given IFid is the same as this one, false otherwise.
  bool operator==(const IFid& a_rhs) const
  {
    if (a_rhs.m_names.size() != m_names.size()) return false;
    for (int i = 0; i < m_names.size(); ++i)
    {
      if (m_names[i] != a_rhs.m_names[i])
        return false;
    }
    return true;
  }

  //! Returns true if the given IFid differs from this one, false otherwise.
  bool operator!=(const IFid& a_rhs) const
  {
    return !(*this == a_rhs);
  }

  private:

  void append(const IFid& a_rhs)
  {
    for (int i = 0; i < a_rhs.m_names.size(); ++i)
      m_names.push_back(a_rhs.m_names[i]);
  }

  // List of names.
  Vector<std::string> m_names;

};

//-----------------------------------------------------------------------
// Template specializations of serialization functions.
//-----------------------------------------------------------------------
template <>
inline int
linearSize(const IFid& a_id)
{
  // An IFid is basically a list of strings.
  return linearSize(a_id.names());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <>
inline void
linearIn(IFid& a_id,
         const void* const inBuf)
{
  Vector<std::string> names;
  linearIn(names, inBuf);
  a_id.setNames(names);
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
template <>
inline void
linearOut(void* const a_outBuf,
          const IFid& a_id)
{
  linearOut(a_outBuf, a_id.names());
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

#endif
