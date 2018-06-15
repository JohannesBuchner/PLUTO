#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CFRegion.H"
#include "CFStencil.H"
#include "NamespaceHeader.H"

void
CFRegion::define(const DisjointBoxLayout& a_grids,
                 const ProblemDomain&     a_domain)
{
  Vector<Box> periodicBoxes;
  CFStencil::buildPeriodicVector(periodicBoxes, a_domain, a_grids);
  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      LayoutData<CFIVS>& lo =  m_loCFIVS[i];
      LayoutData<CFIVS>& hi =  m_hiCFIVS[i];
      lo.define(a_grids);
      hi.define(a_grids);
      for (DataIterator dit(a_grids); dit.ok(); ++dit)
        {
//          lo[dit].define(a_domain, dit(), a_grids, i, Side::Lo);
//          hi[dit].define(a_domain, dit(), a_grids, i, Side::Hi);
          lo[dit].define(a_domain, a_grids.get(dit), periodicBoxes, i, Side::Lo);
          hi[dit].define(a_domain, a_grids.get(dit), periodicBoxes, i, Side::Hi);
        }
    }
  m_defined=true;
}

const CFIVS& CFRegion::loCFIVS(const DataIndex& a_dit, int dir)
{
  CH_assert(m_defined);
  return (m_loCFIVS[dir])[a_dit];
}
const CFIVS& CFRegion::hiCFIVS(const DataIndex& a_dit, int dir)
{
  CH_assert(m_defined);
  return (m_hiCFIVS[dir])[a_dit];
}

const CFRegion& CFRegion::operator=(const CFRegion& a_rhs)
{
  if (a_rhs.m_defined)
    {
      const BoxLayout& grids=(a_rhs.m_loCFIVS[0]).boxLayout();
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          m_loCFIVS[i].define(grids);
          m_hiCFIVS[i].define(grids);
          for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
            {
              m_loCFIVS[i][dit] = a_rhs.m_loCFIVS[i][dit];
              m_hiCFIVS[i][dit] = a_rhs.m_hiCFIVS[i][dit];
            }
        }
      m_defined=true;
    } else
    {
      m_defined=false;
    }
  return *this;
}

void CFRegion::coarsen(int refRatio)
{
  for (DataIterator dit = m_loCFIVS[0].dataIterator(); dit.ok(); ++dit)
    {
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          (m_loCFIVS[i])[dit].coarsen(refRatio);
          (m_hiCFIVS[i])[dit].coarsen(refRatio);
        }
    }
}


#include "NamespaceFooter.H"
