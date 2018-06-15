#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTG 10-03-02

#include "CoveredFaceFAB.H"
#include "NamespaceHeader.H"

// first do simple access functions
// ---------------------------------------------------------
int
CoveredFaceFAB::nComp() const
{
  return m_nComp;
}
// ---------------------------------------------------------
int
CoveredFaceFAB::getIndex(int a_index, Side::LoHiSide a_sd) const
{
  int retval;
  if (a_sd == Side::Lo)
    {
      retval = a_index;
    }
  else
    {
      retval = a_index + SpaceDim;;
    }
  return retval;
}
// ---------------------------------------------------------
const IntVectSet&
CoveredFaceFAB::getIVS(int a_dir, Side::LoHiSide a_sd) const
{
  int index = getIndex(a_dir, a_sd);
  return m_sets[index];
}
// ---------------------------------------------------------
BaseIVFAB<Real>&
CoveredFaceFAB::operator() (const int a_dir, Side::LoHiSide a_sd)
{
  int index = getIndex(a_dir, a_sd);
  return m_data[index];
}
// ---------------------------------------------------------
const BaseIVFAB<Real>&
CoveredFaceFAB::operator() (const int a_dir, Side::LoHiSide a_sd)  const
{
  int index = getIndex(a_dir, a_sd);
  return m_data[index];
}
// ---------------------------------------------------------
CoveredFaceFAB::CoveredFaceFAB()
{
  m_nComp = -1;
}
// ---------------------------------------------------------
CoveredFaceFAB::CoveredFaceFAB(const IntVectSet&  a_ivs,
                               const EBISBox&     a_ebisBox,
                               int                a_nComp)
{
  define(a_ivs, a_ebisBox, a_nComp);
}
// ---------------------------------------------------------
void
CoveredFaceFAB::define(const IntVectSet&  a_ivs,
                       const EBISBox&     a_ebisBox,
                       int                a_nComp)
{
  clear();
  m_isDefined = true;
  m_ebisBox = a_ebisBox;
  m_nComp = a_nComp;

  //compute the covered sets.
  for (SideIterator sit; sit.ok(); ++sit)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          int index = getIndex(idir, sit());
          m_sets[index] = a_ivs;
          for (IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
            {
              const IntVect& iv = ivsit();
              Vector<VolIndex> vofs = m_ebisBox.getVoFs(iv);
              bool allVoFsHaveFaces = true;
              for (int ivof = 0; ivof < vofs.size(); ivof++)
                {
                  const VolIndex& vof = vofs[ivof];
                  Vector<FaceIndex> faces =
                    m_ebisBox.getFaces(vof, idir, sit());
                  if (faces.size() == 0)
                    {
                      allVoFsHaveFaces = false;
                    }
                }
              if (allVoFsHaveFaces)
                m_sets[index] -= iv;
            }
        }
    }

  //define the data holders
  for (int index = 0; index < 2*SpaceDim; index++)
    {
      m_data[index].define(m_sets[index], a_ebisBox.getEBGraph(), a_nComp);
    }
}
// ---------------------------------------------------------
bool
CoveredFaceFAB::isDefined() const
{
  return m_isDefined;
}
// ---------------------------------------------------------
CoveredFaceFAB::~CoveredFaceFAB()
{
  clear();
}
// ---------------------------------------------------------
void
CoveredFaceFAB::clear()
{
  // first delete storage
  for (int index=0; index < 2*SpaceDim; index++)
    {
      m_data[index].clear();
      m_sets[index].makeEmpty();
    }

  m_nComp = -1;
}
// ---------------------------------------------------------
void
CoveredFaceFAB::setVal(const Real& val)
{
  CH_assert(m_nComp > 0);
  for (int index = 0; index < 2*SpaceDim; index++)
    {
      m_data[index].setVal(val);
    }
}

// ---------------------------------------------------------
void
CoveredFaceFAB::copy(const Box& Rfrom, const Interval& Cdest,
                     const Box& Rto, const CoveredFaceFAB& src,
                     const Interval& Csrc)
{
  CH_assert(Rfrom == Rto);
  Box R = Rto;
  for (int index=0; index<2*SpaceDim; index++)
    {
      const BaseIVFAB<Real>& srcFab = src.m_data[index];
      m_data[index].copy(R,Cdest,R, srcFab,Csrc);
    }

}
#include "NamespaceFooter.H"
