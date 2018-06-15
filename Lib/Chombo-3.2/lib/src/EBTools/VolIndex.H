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

#ifndef _VOLINDEX_H_
#define _VOLINDEX_H_

#include "REAL.H"
#include "IntVect.H"
#include "SPMD.H"
#include "BaseIndex.H"
#include "NamespaceHeader.H"

/// Volume of Fluid Index
/**
   VolIndex is a very lightweight object
   to be used to distinguish between volumes of
   fluid in a cell.  We use VolIndex for
   array indexing.   Non-trivial constructors
   are all private because only EBISBox
   is supposed to make valid volume indicies.
*/
class VolIndex: public BaseIndex
{
public:
  ///
  VolIndex();

  virtual ~VolIndex();

  ///
  //VolIndex& operator= (const VolIndex& a_vofin);

  ///
  inline
  bool operator== (const VolIndex& a_vofin) const;

  ///
  inline
  const IntVect& gridIndex() const;

  ///
  inline
  int cellIndex() const;

  ///
  inline
  bool isDefined() const;

  ///
  // VolIndex(const VolIndex& a_vofin);

  ///
  void define(const VolIndex& a_vofin);

  ///
  inline
  VolIndex(const IntVect& a_ix,const  int& a_vofID);

  ///
  inline
  void define(const IntVect& a_ix,const  int& a_vofID);

  inline
  bool operator<(const VolIndex& rhs) const;

  inline
  bool operator!=(const VolIndex& rhs) const;
  /**
     /name linearization routines
  */
  /*@{*/
  int linearSize() const;
  void linearOut(void* const a_outBuf) const;
  void linearIn(const void* const inBuf);
  /*@}*/

  //  not for users
  static size_t iv_offset, index_offset;
  static int initializeOffsets();

private:

  //cell of the VolIndex
  IntVect m_iv;

  //which VoF in the cell (0 to nVoF-1)
  int m_cellIndex;

  bool m_isDefined;
};

std::ostream& operator<< (std::ostream&       os,
                          const VolIndex& iv);

inline bool VolIndex::operator<(const VolIndex& rhs) const
{
  if (m_iv.lexLT(rhs.m_iv))
  {
    return true;
  }
  if (m_iv == rhs.m_iv)
  {
    return m_cellIndex < rhs.m_cellIndex;
  }
  return false;
}
/*****************************************/
inline
const IntVect&
VolIndex::gridIndex() const
{
  return m_iv;
}

/*****************************************/
inline
int
VolIndex::cellIndex() const
{
  return m_cellIndex;
}
/*****************************************/
inline void
VolIndex::define(const IntVect& a_ix,const int& a_vofID)
{
  m_isDefined = true;
  m_iv = a_ix;
  m_cellIndex = a_vofID;
}
/*****************************************/
inline
VolIndex::VolIndex(const IntVect& a_ix,const int& a_vofID)
  :BaseIndex(),m_iv(a_ix), m_cellIndex(a_vofID), m_isDefined(true)
{
  ;
}
/*****************************************/
inline bool
VolIndex::isDefined() const
{
  return m_isDefined;
}

/*****************************************/
inline
VolIndex::VolIndex()
:m_iv(IntVect::Zero), m_cellIndex(-1), m_isDefined(false)
{
  //set index to bogus number to make undefined
  //ones catchable.

}

inline bool
VolIndex::operator== (const VolIndex& a_vofin) const
{
  return ((m_iv == a_vofin.m_iv)&&
          (m_cellIndex == a_vofin.m_cellIndex));
}

inline bool
VolIndex::operator!=(const VolIndex& rhs) const
{
  return !(*this == rhs);
}

#include "NamespaceFooter.H"
#endif
