#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeMGInterp.cpp
// petermc, 20 Jun 2001

#include "NodeMGInterp.H"
#include "REAL.H"
#include "DataIterator.H"
#include "NodeLevelMGF_F.H"
#include "NamespaceHeader.H"

// ---------------------------------------------------------
NodeMGInterp::NodeMGInterp() : is_defined(false)
{
}


// ---------------------------------------------------------
NodeMGInterp::~NodeMGInterp()
{
}


// ---------------------------------------------------------
NodeMGInterp::NodeMGInterp(const DisjointBoxLayout& a_grids,
                           int a_numcomps,
                           int a_refRatio,
                           const Box& a_domain)
  :
  is_defined(false)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_numcomps, a_refRatio, probdomain);
}


// ---------------------------------------------------------
NodeMGInterp::NodeMGInterp(const DisjointBoxLayout& a_grids,
                           int a_numcomps,
                           int a_refRatio,
                           const ProblemDomain& a_domain)
  :
  is_defined(false)
{
  define(a_grids, a_numcomps, a_refRatio, a_domain);
}


// ---------------------------------------------------------
void
NodeMGInterp::define(const DisjointBoxLayout& a_grids,
                     int a_numcomps,
                     int a_refRatio,
                     const Box& a_domain)
{
  ProblemDomain probdomain(a_domain);
  define(a_grids, a_numcomps, a_refRatio, probdomain);
}


// ---------------------------------------------------------
void
NodeMGInterp::define(const DisjointBoxLayout& a_grids,
                     int a_numcomps,
                     int a_refRatio,
                     const ProblemDomain& a_domain)
{
  m_refRatio = a_refRatio;
  m_domain = a_domain;
  m_grids = a_grids;

  m_boxRef = Box(IntVect::Zero, (m_refRatio-1)*IntVect::Unit);
  Box corners(IntVect::Zero, IntVect::Unit);
  m_weights.define(corners, m_boxRef.numPts());
  FORT_NODEINTERPMG_GETWEIGHTS(CHF_CONST_INT(m_refRatio),
                               CHF_BOX(m_boxRef),
                               CHF_FRA(m_weights));

  // create the work array
  DisjointBoxLayout coarsenedGrids;
  coarsen(coarsenedGrids, a_grids, m_refRatio);

  m_coarsenedFine.define(coarsenedGrids, a_numcomps);

  is_defined = true;
}



// ---------------------------------------------------------
// interpolate from coarse level to fine level
void
NodeMGInterp::interpToFine(LevelData<NodeFArrayBox>& a_fine,
                           const LevelData<NodeFArrayBox>& a_coarse,
                           bool a_sameGrids) // a_sameGrids default false
{
  CH_assert(is_defined);
  const int nComp = a_fine.nComp();
  CH_assert(a_coarse.nComp() == nComp);

  // Copy a_coarse to m_coarsenedFine on grids of m_coarsenedFine.

  // petermc, 15 Nov 2002:
  // You don't need a Copier for this copyTo, because
  // the grid layout of the destination, m_coarsenedFine,
  // will be contained in that of the source, a_coarse.
  if (! a_sameGrids)
    {
      a_coarse.copyTo(a_coarse.interval(),
                      m_coarsenedFine,
                      m_coarsenedFine.interval() );
    }

  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Box fineBox = m_grids.get(dit());
      Box crseBox = coarsen(fineBox, m_refRatio);

      const FArrayBox& crseFab = (a_sameGrids) ?
        a_coarse[dit()].getFab() : m_coarsenedFine[dit()].getFab();

      FArrayBox& fineFab = a_fine[dit()].getFab();

      FORT_NODEINTERPMG(CHF_FRA(fineFab),
                        CHF_CONST_FRA(crseFab),
                        CHF_BOX(crseBox),
                        CHF_CONST_INT(m_refRatio),
                        CHF_BOX(m_boxRef),
                        CHF_FRA(m_weights));

      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
}


// ---------------------------------------------------------
bool
NodeMGInterp::isDefined() const
{
  return ( is_defined );
}

#include "NamespaceFooter.H"
