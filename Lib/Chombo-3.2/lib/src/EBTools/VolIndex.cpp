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

#include "VolIndex.H"
#include "parstream.H"
#include "CH_HDF5.H"
#include "NamespaceHeader.H"

using std::ostream;

VolIndex vold1;

size_t VolIndex::iv_offset =  CHOFFSET(vold1, m_iv.vect);
size_t VolIndex::index_offset = CHOFFSET(vold1, m_cellIndex);

using std::endl;

int VolIndex::linearSize() const
{
  return (CH_SPACEDIM + 1)*sizeof(int);
}

void VolIndex::linearOut(void* const a_outBuf) const
{
  CH_assert(m_isDefined);
  int* buf = (int*)a_outBuf;
  D_TERM(buf[0]=m_iv[0],; buf[1]=m_iv[1],; buf[2]=m_iv[2]);
  buf[CH_SPACEDIM]=m_cellIndex;

}

void VolIndex::linearIn(const void* const inBuf)
{
  int* buf = (int*)inBuf;
  D_TERM(m_iv[0]=buf[0],; m_iv[1]=buf[1],; m_iv[2]=buf[2]);
  m_cellIndex = buf[CH_SPACEDIM];
  m_isDefined = true;
}

/*****************************************/
ostream&
operator<< (ostream&       os,
            const VolIndex& p)
{
  IntVect iv = p.gridIndex();
  int vofind = p.cellIndex();
  os << D_TERM( "((" << iv[0] , <<
                ',' << iv[1] , <<
                ',' << iv[2])  << ")(" << vofind << "))" ;
  if (os.fail())
    MayDay::Error("operator<<(ostream&,VolIndex&) failed");
  return os;
}
/*****************************************/
int VolIndex::initializeOffsets()
{
  //iv_offset = offsetof(VolIndex, m_iv.vect);
  //index_offset = offsetof(VolIndex, m_cellIndex);
  return 1;
}

static int init =  VolIndex::initializeOffsets();

/*****************************************/
void
VolIndex::define(const VolIndex& a_vofin)
{
  //must be allowed to propogate undefinednitude
  //because of vector<vol>
  //  CH_assert(a_vofin.isDefined());

  m_isDefined = a_vofin.m_isDefined;
  m_iv = a_vofin.m_iv;
  m_cellIndex = a_vofin.m_cellIndex;
}

VolIndex::~VolIndex()
{
}
/*****************************************/
#include "NamespaceFooter.H"
