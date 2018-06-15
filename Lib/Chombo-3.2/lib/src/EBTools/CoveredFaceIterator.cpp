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

#include "CoveredFaceIterator.H"
#include "NamespaceHeader.H"

/********************************/
CoveredFaceIterator::CoveredFaceIterator()
{
  m_isDefined = false;
}
/********************************/
CoveredFaceIterator::~CoveredFaceIterator()
{
}
/********************************/
CoveredFaceIterator::CoveredFaceIterator(const IntVectSet& a_ivs,
                                         const EBISBox& a_ebisBox,
                                         int a_dir, Side::LoHiSide a_side)
{
  define(a_ivs, a_ebisBox, a_dir, a_side);
}
/********************************/
void
CoveredFaceIterator::define(const IntVectSet& a_ivs,
                            const EBISBox& a_ebisBox,
                            int a_dir, Side::LoHiSide a_side)
{
  m_isDefined = true;
  m_vols.resize(0);

  for (IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      Vector<VolIndex> vofs = a_ebisBox.getVoFs(iv);
      for (int ivof = 0; ivof < vofs.size(); ivof++)
        {
          const VolIndex& vof = vofs[ivof];
          Vector<FaceIndex> faces = a_ebisBox.getFaces(vof, a_dir, a_side);
          if (faces.size() == 0)
            {
              m_vols.push_back(vof);
            }
        }
    }
  reset();
}

/********************************/
void
CoveredFaceIterator::reset()
{
  CH_assert(isDefined());
  m_ivol = 0;
}
/********************************/
//keys off if icell >= m_cellVols.size()
/********************************/
void
CoveredFaceIterator::operator++()
{
  CH_assert(isDefined());
  m_ivol++;
}
/********************************/
const VolIndex&
CoveredFaceIterator::operator() () const
{
  CH_assert(isDefined());
  CH_assert(m_ivol < m_vols.size());
  return m_vols[m_ivol];
}
/********************************/
bool
CoveredFaceIterator::ok() const
{
  return (m_ivol < m_vols.size());
}
/********************************/
bool
CoveredFaceIterator::isDefined() const
{
  return m_isDefined;
}
/********************************/
#include "NamespaceFooter.H"
