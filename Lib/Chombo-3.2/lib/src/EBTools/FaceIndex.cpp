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

#include "FaceIndex.H"
#include "MayDay.H"
#include "CH_HDF5.H"
#include "NamespaceHeader.H"

using std::cerr;
using std::cout;
using std::endl;

FaceIndex faced1;

size_t FaceIndex::lo_offset = CHOFFSET(faced1, m_loiv.vect);
size_t FaceIndex::hi_offset = CHOFFSET(faced1, m_hiiv.vect);
size_t FaceIndex::rest_offset = CHOFFSET(faced1, m_loIndex);


/*****************************************/
FaceIndex::FaceIndex()
{
  m_loiv = IntVect::Zero;
  m_hiiv = IntVect::Zero;
  m_loIndex = -7;
  m_hiIndex = -7;
  m_direction = -7;
  m_isDefined = false;
  m_isBoundary = false;

}
/*****************************************/
const bool&
FaceIndex::isDefined() const
{
  return m_isDefined;
}
/******************************/
const int&
FaceIndex::direction() const
{
  return m_direction;
}

/*****************************************/
const bool&
FaceIndex::isBoundary() const
{
  return m_isBoundary;
}
/*****************************************/
FaceIndex::FaceIndex(const FaceIndex& a_facein)
{
  define(a_facein);
}
/*****************************************/
void
FaceIndex::define(const FaceIndex& a_facein)
{
  //must be allowed to propogate undefinednitude
  //because of vector<vol>
  m_isDefined = a_facein.m_isDefined;
  m_isBoundary = a_facein.m_isBoundary;
  m_loiv = a_facein.m_loiv;
  m_hiiv = a_facein.m_hiiv;
  m_loIndex = a_facein.m_loIndex;
  m_hiIndex = a_facein.m_hiIndex;
  m_direction = a_facein.m_direction;
}

/*****************************************/
FaceIndex::FaceIndex(const VolIndex& a_vof1,
                     const VolIndex& a_vof2,
                     const int& a_direction)
{
  define(a_vof1, a_vof2, a_direction);
}

/*****************************************/
FaceIndex::FaceIndex(const VolIndex& a_vof1,
                     const VolIndex& a_vof2)
{
  define(a_vof1, a_vof2);
}

/*****************************************/
void
FaceIndex::define(const VolIndex& a_vof1,
                  const VolIndex& a_vof2,
                  const int& a_direction)
{
  CH_assert((a_direction >= 0) &&
            (a_direction < SpaceDim));
  //one vof or the other has to be in the domain.
  CH_assert((a_vof1.cellIndex() >= 0) || (a_vof2.cellIndex() >= 0));

  //if either vof is outside the domain--shown by a negative
  //cell index, the face is a bounary face
  m_isBoundary =
    ((a_vof1.cellIndex() < 0) || (a_vof2.cellIndex() < 0));
  m_isDefined  = true;

  const IntVect& iv1 = a_vof1.gridIndex();
  const IntVect& iv2 = a_vof2.gridIndex();

  int isign = 1;
  if (iv1[a_direction] < iv2[a_direction])
    {
      m_loiv = iv1;
      m_hiiv = iv2;
      m_loIndex = a_vof1.cellIndex();
      m_hiIndex = a_vof2.cellIndex();
      isign = -1;
    }
  else
    {
      m_loiv = iv2;
      m_hiiv = iv1;
      m_loIndex = a_vof2.cellIndex();
      m_hiIndex = a_vof1.cellIndex();
      isign =  1;
    }

  m_direction = a_direction;
  //this check needs to go away in the land
  //of periodic boundary conditions
  if ((isign*(iv1-iv2)) != BASISV(a_direction))
    {
      cerr << "FaceIndex constructor--bad arguments" << endl;
      cerr << "iv1, iv2, direction =" << iv1 << iv2 << a_direction << endl;
      MayDay::Error("FaceIndex constructor--bad arguments");
    }
}
/*****************************************/
void
FaceIndex::define(const VolIndex& a_vof1,
                  const VolIndex& a_vof2)
{
  const IntVect& iv1 = a_vof1.gridIndex();
  const IntVect& iv2 = a_vof2.gridIndex();

  IntVect div = absolute(iv1 - iv2);
  int direction = -1;

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (div == BASISV(idir))
        {
          direction = idir;
          break;
        }
    }

  if (direction < 0)
    MayDay::Error("a_vof1 and a_vof2 are not neighbors!");

  define(a_vof1, a_vof2, direction);
}

/*****************************************/
int FaceIndex::faceSign(const VolIndex& a_vof) const
{
  if (a_vof.gridIndex() == m_loiv)
    return 1;
  else if (a_vof.gridIndex() == m_hiiv)
    return -1;
  else
    return 0;
}

/*****************************************/
const int&
FaceIndex::cellIndex(const Side::LoHiSide& a_sd) const
{
  CH_assert(isDefined());
  const int* retval;
  if (a_sd == Side::Lo)
    retval=  &m_loIndex;
  else
    retval=  &m_hiIndex;

  return *retval;
}
/*****************************************/
const IntVect&
FaceIndex::gridIndex(const Side::LoHiSide& a_sd) const
{
  CH_assert(isDefined());
  const IntVect* retval;
  if (a_sd == Side::Lo)
    retval=  &m_loiv;
  else
    retval=  &m_hiiv;

  return *retval;
}
/*****************************************/
VolIndex
FaceIndex::getVoF(const Side::LoHiSide& a_sd) const
{
  CH_assert(isDefined());
  VolIndex retval;
  if (a_sd == Side::Lo)
    retval=  VolIndex(m_loiv, m_loIndex);
  else
    retval=  VolIndex(m_hiiv, m_hiIndex);

  return retval;
}
/*****************************************/
FaceIndex::~FaceIndex()
{
}

/*****************************************/
FaceIndex&
FaceIndex::operator= (const FaceIndex& a_facein)
{
  if (&a_facein != this)
    {
      //CH_assert(a_facein.isDefined());
      //must be allowed to propogate undefinednitude
      //because of vector<vol>
      m_isDefined = a_facein.m_isDefined;
      m_isBoundary = a_facein.m_isBoundary;
      m_loiv = a_facein.m_loiv;
      m_hiiv = a_facein.m_hiiv;
      m_loIndex = a_facein.m_loIndex;
      m_hiIndex = a_facein.m_hiIndex;
      m_direction = a_facein.m_direction;
    }
  return *this;
}

/******************************/
bool FaceIndex::operator==(const FaceIndex& a_facein) const
{
  return((m_loiv == a_facein.m_loiv) &&
         (m_hiiv == a_facein.m_hiiv) &&
         (m_loIndex == a_facein.m_loIndex) &&
         (m_hiIndex == a_facein.m_hiIndex));
}

/******************************/
bool FaceIndex::operator!=(const FaceIndex& a_facein) const
{

  return (!(*this == a_facein));
}

template < >
int linearSize<FaceIndex>(const FaceIndex& findex)
{
  return sizeof(FaceIndex);
}
//FaceIndex specialization of linearIn
template < >
void linearIn<FaceIndex>(FaceIndex& a_outputT, const void* const inBuf)
{
  unsigned char* bob = (unsigned char*)inBuf;
  unsigned char* to = (unsigned char*)&a_outputT;

  memcpy(to + FaceIndex::lo_offset, bob, SpaceDim*sizeof(int));
  bob += SpaceDim*sizeof(int);
  memcpy(to + FaceIndex::hi_offset, bob, SpaceDim*sizeof(int));
  bob += SpaceDim*sizeof(int);
  memcpy(to + FaceIndex::rest_offset, bob, sizeof(FaceIndex)-FaceIndex::rest_offset);
}

//FaceIndex specialization of linearOut
template < >
void linearOut<FaceIndex>(void* const a_outBuf, const FaceIndex& a_inputT)
{
  unsigned char* bob = (unsigned char*)a_outBuf;
  const unsigned char* from = (const unsigned char*)&a_inputT;

  memcpy(bob, from + FaceIndex::lo_offset, SpaceDim*sizeof(int));
  bob += SpaceDim*sizeof(int);
  memcpy(bob, from + FaceIndex::hi_offset, SpaceDim*sizeof(int));
  bob += SpaceDim*sizeof(int);
  memcpy(bob,  from + FaceIndex::rest_offset, sizeof(FaceIndex)-FaceIndex::rest_offset);
}


ostream& operator<<( ostream& out, const FaceIndex& fi )
{
  out << "[(" << fi.m_loiv << "):" << fi.m_loIndex << ";"
      << " (" << fi.m_hiiv << "):" << fi.m_hiIndex << "]";
  return out;
}
#include "NamespaceFooter.H"
