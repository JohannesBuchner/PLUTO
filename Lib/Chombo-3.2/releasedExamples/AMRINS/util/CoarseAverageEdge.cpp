#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CoarseAverageEdge.H"
#include "AverageEdgeF_F.H"

// ----------------------------------------------------------
CoarseAverageEdge::CoarseAverageEdge() :m_isDefined(false)
{
}

// ----------------------------------------------------------
CoarseAverageEdge::~CoarseAverageEdge()
{
}

// ----------------------------------------------------------
CoarseAverageEdge::CoarseAverageEdge(const DisjointBoxLayout& a_fineGrids,
                                     int a_nComp, int a_nRef)
  : m_isDefined(false)
{
  define(a_fineGrids, a_nComp, a_nRef);
}

// ----------------------------------------------------------
void
CoarseAverageEdge::define(const DisjointBoxLayout& a_fineGrids,
                          int a_nComp, int a_nRef)
{
  m_nRef = a_nRef;

  DisjointBoxLayout coarsened_fine_domain;
  coarsen(coarsened_fine_domain, a_fineGrids, m_nRef);
  m_coarsenedFineData.define(coarsened_fine_domain, a_nComp);

  m_isDefined = true;
}

// ----------------------------------------------------------
bool
CoarseAverageEdge::isDefined() const
{
  return m_isDefined;
}

// ----------------------------------------------------------
// this function is shamelessly based on the ANAG CoarseAverage
// (cell-centered) version
void
CoarseAverageEdge::averageToCoarse(LevelData<FluxBox>& a_coarseData,
                                   const LevelData<FluxBox>& a_fineData)
{
  CH_assert(isDefined());

  DataIterator dit = a_fineData.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& coarsenedFine = m_coarsenedFineData[dit()];
      const FluxBox& fine = a_fineData[dit()];

      // coarsen from the entire fine grid onto the entire coarse grid
      averageGridData(coarsenedFine, fine);
    }

  m_coarsenedFineData.copyTo(m_coarsenedFineData.interval(),
                             a_coarseData, a_coarseData.interval());
}

// ----------------------------------------------------------
void
CoarseAverageEdge::averageGridData(FluxBox& a_coarsenedFine,
                                   const FluxBox& a_fine) const
{
  for (int dir=0; dir<SpaceDim; dir++)
    {
      FArrayBox& coarseFab = a_coarsenedFine[dir];
      const FArrayBox& fineFab = a_fine[dir];

      const Box& coarseBox = coarseFab.box();

      // set up refinement box
      int boxHi = m_nRef-1;
      IntVect hiVect(D_DECL(boxHi,boxHi,boxHi));
      // don't want to index at all in dir direction --
      // instead, want to just march along edge.
      hiVect.setVal(dir,0);
      IntVect loVect(D_DECL(0,0,0));
      Box refBox(loVect, hiVect);

      FORT_AVERAGEEDGE( CHF_FRA(coarseFab),
                        CHF_CONST_FRA(fineFab),
                        CHF_BOX(coarseBox),
                        CHF_CONST_INT(dir),
                        CHF_CONST_INT(m_nRef),
                        CHF_BOX(refBox));
    }
}
