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

#include "VoFIterator.H"
#include "EBGraph.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

/********************************/
const Vector<VolIndex>
VoFIterator::getVector() const
{
  CH_assert(m_isDefined);
  return m_vols;
}
/********************************/
VoFIterator::VoFIterator()
{
  m_isDefined = false;
}
/********************************/
VoFIterator::~VoFIterator()
{
}
/********************************/
VoFIterator::VoFIterator(const IntVectSet& a_ivs,
                         const EBGraph& a_ebgraph)
{
  define(a_ivs, a_ebgraph);
}
/********************************/
void
VoFIterator::define(const IntVectSet& a_ivs,
                    const EBGraph& a_ebgraph)
{
  CH_TIMELEAF("VoFIterator::define");
  //can't do this because minbox is broken
  //  CH_assert((a_ebgraph.getRegion().contains(a_ivs.minBox())||
  //          a_ivs.isEmpty()));
  m_isDefined = true;
  m_vols.resize(0);
  for (IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
    {
      m_vols.append(a_ebgraph.getVoFs(ivsit()));
    }
  reset();
}
/********************************/
void
VoFIterator::reset()
{
  CH_assert(isDefined());
  m_ivol = 0;
}
/********************************/
void
VoFIterator::operator++()
{
  CH_assert(isDefined());
  m_ivol++;
}
/********************************/
const VolIndex&
VoFIterator::operator() () const
{
  CH_assert(isDefined());
  CH_assert(m_ivol < m_vols.size());
  return m_vols[m_ivol];
}
/********************************/
bool
VoFIterator::ok() const
{
  return (m_ivol < m_vols.size());
}
/********************************/
bool
VoFIterator::isDefined() const
{
  return m_isDefined;
}
/********************************/
#include "NamespaceFooter.H"
