#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CoarseAverageCodim.H"
#include "AverageCodimF_F.H"
#include "NamespaceHeader.H"


/******************************************************************************/
/**
 *  \file
 *
 *  \brief Non-inline definitions for CoarseAverageCodim class
 *
 *\\*+*************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default constructor
/*--------------------------------------------------------------------*/

CoarseAverageCodim::CoarseAverageCodim()
  :m_isDefined(false),
   m_isAveraged(false)
{
}

/*--------------------------------------------------------------------*/
//  Full constructor
/** \param[in]  a_fineGrids
 *                      Fine grid layout
 *  \param[in]  a_codim The codimension
 *  \param[in]  a_nComp Number of components to average
 *  \param[in]  a_nRef  The refinement ratio
 *//*-----------------------------------------------------------------*/

CoarseAverageCodim::CoarseAverageCodim(const DisjointBoxLayout& a_fineGrids,
                                       const int                a_codim,
                                       const int                a_nComp,
                                       const int                a_nRef)
  :m_isDefined(false),
   m_isAveraged(false)
{
  define(a_fineGrids, a_codim, a_nComp, a_nRef);
}

/*--------------------------------------------------------------------*/
///  Destructor
/*--------------------------------------------------------------------*/

CoarseAverageCodim::~CoarseAverageCodim()
{
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Defines the object
/** \param[in]  a_fineGrids
 *                      Fine grid layout
 *  \param[in]  a_codim The codimension
 *  \param[in]  a_nComp Number of components to average
 *  \param[in]  a_nRef  The refinement ratio
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::define(const DisjointBoxLayout& a_fineGrids,
                           const int                a_codim,
                           const int                a_nComp,
                           const int                a_nRef)
{
  m_nRef = a_nRef;

  DisjointBoxLayout coarsened_fine_domain;
  coarsen(coarsened_fine_domain, a_fineGrids, m_nRef);
  m_coarsenedFineData.define(coarsened_fine_domain, a_nComp, IntVect::Zero,
                             CodimBoxFactory<FArrayBox> (a_codim));

  m_isDefined = true;
  m_isAveraged = false;
}

/*--------------------------------------------------------------------*/
//  Averages fine-level data to internal coarse representation of fine
//  grid
/** \param[in]  a_fineData
 *                      Fine grid data
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::average(const LevelData<CodimBox<FArrayBox> >& a_fineData)
{
  average(a_fineData, arithmetic, m_nRef);
}

/*--------------------------------------------------------------------*/
//  Harmonic averaging of fine data to coarse representation of fine
//  grid
/** \param[in]  a_fineData
 *                      Fine grid data
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::averageHarmonic(const LevelData<CodimBox<FArrayBox> >& a_fineData)
{
  average(a_fineData, harmonic, m_nRef);
}

/*--------------------------------------------------------------------*/
//  Summation of fine data to internal coarse representation of fine
//  grid
/** \param[in]  a_fineData
 *                      Fine grid data
 *  \param[in]  a_refFactor
 *                      Sum of fine values is divided by
 *                      a_refFactor^(CH_SPACEDIM-1).  For sums this
 *                      should be set to one (default).
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::sum(const LevelData<CodimBox<FArrayBox> >& a_fineData,
                        const int                  a_refFactor)
{
  average(a_fineData, arithmetic, a_refFactor);
}

/*--------------------------------------------------------------------*/
//  Obtain averaged results by copying to the destination
/** \param[out] a_coarseData
 *                      Coarse grid data
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::copyTo(LevelData<CodimBox<FArrayBox> >& a_coarseData)
{
  CH_assert(m_isAveraged);

  // If coarseData's DisjointBoxLayout is not a simple coarsenening of
  // the fine one, then it needs to have at least one ghost cell in
  // order to ensure that this copyTo is done correctly. In
  // particular, this is required in order to ensure that we handle
  // the case where the coarse-fine interface is coincident with a
  // coarse-coarse boundary. The other solution to this would be to
  // build a specialized Copier for LevelData<CodimBox<FArrayBox> >, but we're
  // hoping to avoid that for now...

  if ((a_coarseData.ghostVect() == IntVect::Zero) &&
      !(a_coarseData.getBoxes().compatible(m_coarsenedFineData.getBoxes())))
    {
      MayDay::Error("CoarseAverageCodim requires that coarse data which is not "
                    "a coarsenening of the fine grids have at least one ghost "
                    "cell");
    }

  // Copy to the coarse data
  m_coarsenedFineData.copyTo(m_coarsenedFineData.interval(),
                             a_coarseData,
                             a_coarseData.interval());
}

/*--------------------------------------------------------------------*/
//  Averages fine-level data to coarse level
/** \param[out] a_coarseData
 *                      Coarse grid data
 *  \param[in]  a_fineData
 *                      Fine grid data
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::averageToCoarse(LevelData<CodimBox<FArrayBox> >&       a_coarseData,
                                    const LevelData<CodimBox<FArrayBox> >& a_fineData)
{
  computeAverages(a_coarseData, a_fineData, arithmetic);
}

/*--------------------------------------------------------------------*/
//  Averages fine-level data to coarse level using harmonic averaging
/** \param[out] a_coarseData
 *                      Coarse grid data
 *  \param[in]  a_fineData
 *                      Fine grid data
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::averageToCoarseHarmonic(
  LevelData<CodimBox<FArrayBox> >&       a_coarseData,
  const LevelData<CodimBox<FArrayBox> >& a_fineData)
{
  computeAverages(a_coarseData, a_fineData, harmonic);
}

/*--------------------------------------------------------------------*/
//  Utility function to completely determine averages.
/** Called by both averageToCoarse and averageCoarseHarmonic (to
 *  avoid code duplication)
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::computeAverages(LevelData<CodimBox<FArrayBox> >&       a_coarseData,
                                    const LevelData<CodimBox<FArrayBox> >& a_fineData,
                                    const AverageType          a_averageType)
{
  average(a_fineData, a_averageType, m_nRef);
  copyTo(a_coarseData);
}

/*--------------------------------------------------------------------*/
//  Utility for averaging fine-level data to internal coarse
//  representation of fine grid
/** \param[in]  a_fineData
 *                      Fine grid data
 *  \param[in]  a_averageType
 *                      Type of averaging (arithmetic or harmonic)
 *  \param[in]  a_refFactor
 *                      Sum of fine values is divided by
 *                      a_refFactor^(CH_SPACEDIM-1).  Normally this is
 *                      the refinement ratio
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::average(const LevelData<CodimBox<FArrayBox> >& a_fineData,
                            const AverageType          a_averageType,
                            const int                  a_refFactor)
{
  CH_assert(isDefined());
  DataIterator dit = a_fineData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      CodimBox<FArrayBox> & coarsenedFine = m_coarsenedFineData[dit()];
      const CodimBox<FArrayBox> & fine = a_fineData[dit()];

      // Coarsen from the entire fine grid onto the entire coarse grid
      averageGridData(coarsenedFine, fine, a_averageType, a_refFactor);
    }
  m_isAveraged = true;
}

/*--------------------------------------------------------------------*/
//  Averages entire single grid data from fine->coarse
/** \param[out] a_coarsenedFine
 *                      Data on coarsened box aligned with fine box
 *  \param[in]  a_fine  Data on fine box
 *  \param[in]  a_averageType
 *                      Type of averaging (arithmetic or harmonic)
 *  \param[in]  a_refFactor
 *                      Sum of fine values is divided by
 *                      a_refFactor^(CH_SPACEDIM-1).  Normally this is
 *                      the refinement ratio
 *//*-----------------------------------------------------------------*/

void
CoarseAverageCodim::averageGridData(CodimBox<FArrayBox> &         a_coarsenedFine,
                                    const CodimBox<FArrayBox> &   a_fine,
                                    const AverageType a_averageType,
                                    const int         a_refFactor) const
{
  const int codim = a_coarsenedFine.getCodim();
  const int nOr = a_coarsenedFine.getNumOrient();
  int oDir[SpaceDim];
  for (int iOr = 0; iOr != nOr; ++iOr)
    {
      FArrayBox& fabC = a_coarsenedFine.getSequential(iOr);
      const FArrayBox& fabF = a_fine.getSequential(iOr);
      const Box& box = fabC.box();

      // Construct the localized box to iterate over to accumulate the fine
      // average
      IntVect hi((m_nRef-1)*IntVect::Unit);
      // Do not accumulate in orthogonal directions
      a_coarsenedFine.getDirections(iOr, oDir);
      for (int i = 0; i != codim; ++i)
        {
          hi[oDir[i]] = 0;
        }
      Box refBox(IntVect::Zero, hi);

      const int refDim = SpaceDim - codim;

      switch (a_averageType)
        {
        case arithmetic:
          FORT_AVERAGECODIM( CHF_FRA(fabC),
                             CHF_CONST_FRA(fabF),
                             CHF_BOX(box),
                             CHF_CONST_INT(m_nRef),
                             CHF_CONST_INT(a_refFactor),
                             CHF_CONST_INT(refDim),
                             CHF_BOX(refBox));
          break;
        case harmonic:
          FORT_AVERAGECODIMHARMONIC( CHF_FRA(fabC),
                                     CHF_CONST_FRA(fabF),
                                     CHF_BOX(box),
                                     CHF_CONST_INT(m_nRef),
                                     CHF_CONST_INT(a_refFactor),
                                     CHF_CONST_INT(refDim),
                                     CHF_BOX(refBox));
          break;
        }
    }
}

#include "NamespaceFooter.H"
