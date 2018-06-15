#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>

#include "FourthOrderFillPatch.H"
#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderFillPatch::FourthOrderFillPatch()
{
  m_defined = false;
  m_timeInterpDefined = false;
}


//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
FourthOrderFillPatch::~FourthOrderFillPatch()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void FourthOrderFillPatch::define(/// layout at this level
                                  const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                                  /// layout at coarser level
                                  const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                                  /// number of variables
                                  const int&                a_numStates,
                                  /// problem domain on the coarser level
                                  const ProblemDomain&      a_coarseDomain,
                                  /// refinement ratio between this level and the coarser level
                                  const int&                a_refineCoarse,
                                  /// number of layers of ghost cells to fill by interpolation
                                  const int&                a_interpRadius,
                                  /// whether this object is for a fixed time
                                  bool                      a_fixedTime,
                                  /// dimensions that are fixed, not interpolated
                                  Interval                  a_fixedDims)
{
  // Cache data
  m_numStates = a_numStates;
  m_coarseDomain = a_coarseDomain;
  m_refineCoarse = a_refineCoarse;
  m_interpRadius = a_interpRadius;
  m_fixedDims = a_fixedDims;

  m_layout = a_thisDisjointBoxLayout;
  m_coarseLayout = a_coarserDisjointBoxLayout;

  m_refineVect = m_refineCoarse * IntVect::Unit;
  for (int dirf = m_fixedDims.begin(); dirf <= m_fixedDims.end(); dirf++)
    {
      m_refineVect[dirf] = 1;
    }

  ProblemDomain fineDomain = refine(m_coarseDomain, m_refineVect);

  // width of ghost layer of coarse cells with fine cells to fill in
  // int ghostCoarsened = ceil(m_interpRadius / m_refineCoarse);
  int ghostsCoarsened = m_interpRadius / m_refineCoarse;
  if (ghostsCoarsened * m_refineCoarse < m_interpRadius) ghostsCoarsened++;

  // FourthOrderFineInterp m_spaceInterpolator;
  m_spaceInterpolator.define(m_layout, m_numStates,
                             m_refineCoarse, fineDomain, ghostsCoarsened,
                             m_fixedDims);

  if (!a_fixedTime)
    {
      const IntVect& ghostVect =
        m_spaceInterpolator.coarsenedFineData().ghostVect();
      int ghosts = ghostVect[0];
      CH_assert(ghostVect == ghosts * IntVect::Unit);

      // TimeInterpolatorRK4 m_timeInterpolator;
      // we'll need intermediate-time data on ghosted coarsened m_layout
      // (not only ghost cells), and that's exactly what this will give us.
      // FIXME, 4 Oct 2012: time interpolator still isotropic refinement.
      m_timeInterpolator.define(m_layout, m_coarseLayout, fineDomain,
                                m_refineCoarse, m_numStates, ghosts);
      m_timeInterpDefined = true;
    }

  DisjointBoxLayout m_layoutCoarsened;
  coarsen(m_layoutCoarsened, m_layout, m_refineVect);

  // Use m_spaceInterpolator.m_coarseData instead of this.
  // LevelData<FArrayBox> m_coarsenedFineData;
  //  m_coarsenedFineData.define(m_layoutCoarsened,
  //                             m_numStates,
  //                             m_interpRadius * IntVect::Unit);

  // LayoutData<IntVectSet> m_coarsenedGhosts;
  m_coarsenedGhosts.define(m_layoutCoarsened);

  // I copied this code segment from PiecewiseLinearFillPatch::define().
  // We find LayoutData<IntVectSet> m_coarsenedGhosts, using:
  // ProblemDomain m_coarseDomain;
  // DisjointBoxLayout m_layoutCoarsened;
  // int ghostsCoarsened;

  // a box which will determine whether a given box
  // adjoins a periodic boundary
  Box periodicTestBox(m_coarseDomain.domainBox());
  if (m_coarseDomain.isPeriodic())
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (m_coarseDomain.isPeriodic(idir))
            periodicTestBox.grow(idir,-1);
        }
    }

  DataIterator dit = m_layoutCoarsened.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& coarseBox = m_layoutCoarsened[dit];
      Box coarseBoxWithGhosts =
        grow(coarseBox, ghostsCoarsened) & m_coarseDomain;

      // Initially set coarsenedGhostsHere to the whole ghosted coarse box.
      IntVectSet& coarsenedGhostsHere = m_coarsenedGhosts[dit];
      coarsenedGhostsHere = IntVectSet(coarseBoxWithGhosts);

      // Iterate over boxes in coarsened fine layout, and subtract off
      // from the set of coarse cells from which the fine ghost cells
      // will be interpolated.
      LayoutIterator litOther = m_layoutCoarsened.layoutIterator();
      for (litOther.begin(); litOther.ok(); ++litOther)
        {
          const Box& coarseOtherBox = m_layoutCoarsened[litOther];
          coarsenedGhostsHere -= coarseOtherBox;
          // also need to remove periodic images from list of cells
          // to be filled, since they will be filled through exchange
          // as well
          if (m_coarseDomain.isPeriodic()
              && !periodicTestBox.contains(coarseOtherBox)
              && !periodicTestBox.contains(coarseBox))
            {
              ShiftIterator shiftIt = m_coarseDomain.shiftIterator();
              IntVect shiftMult(m_coarseDomain.domainBox().size());
              Box shiftedBox(coarseOtherBox);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect = shiftMult*shiftIt();
                  shiftedBox.shift(shiftVect);
                  coarsenedGhostsHere -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }
    }

  // Everything is defined now.
  m_defined = true;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFillPatch::fillInterp(/// interpolated solution on this level
                                      LevelData<FArrayBox>&         a_fineData,
                                      /// time interpolation coefficient, in [0:1]
                                      Real                          a_timeInterpCoeff,
                                      /// starting coarse data component
                                      int                           a_srcComp,
                                      /// starting fine data component
                                      int                           a_dstComp,
                                      /// number of data components to interpolate
                                      int                           a_numComp)
{
  CH_assert(m_defined);
  CH_assert(m_timeInterpDefined);

  const Interval srcInterval(a_srcComp, a_srcComp + a_numComp-1);
  // Requires that m_timeInterpolator already has data from the
  // coarser level giving us the Taylor polynomial.
  m_timeInterpolator.interpolate(m_spaceInterpolator.coarsenedFineData(),
                                 a_timeInterpCoeff, srcInterval);

  // Interpolate to a_fineData from m_spaceInterpolator.coarsenedFineData(),
  // on the given components.
  fillInterpSpaceFromCoarsened(a_fineData, a_srcComp, a_dstComp, a_numComp);
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFillPatch::fillRK4Intermediate(/// intermediate RK4 solution on this level
                                               LevelData<FArrayBox>&         a_fineData,
                                               /// time interpolation coefficient, in [0:1]
                                               Real                          a_timeInterpCoeff,
                                               /// RK4 stage, in 0:3
                                               int                           a_stage,
                                               /// starting coarse data component
                                               int                           a_srcComp,
                                               /// starting fine data component
                                               int                           a_dstComp,
                                               /// number of data components to interpolate
                                               int                           a_numComp)
{
  CH_assert(m_defined);
  CH_assert(m_timeInterpDefined);

  const Interval srcInterval(a_srcComp, a_srcComp + a_numComp-1);
  // Requires that m_timeInterpolator already has data from the
  // coarser level giving us the Taylor polynomial.
  m_timeInterpolator.intermediate(m_spaceInterpolator.coarsenedFineData(),
                                  a_timeInterpCoeff, a_stage, srcInterval);

  // Interpolate to a_fineData from m_spaceInterpolator.coarsenedFineData(),
  // on the given components.
  fillInterpSpaceFromCoarsened(a_fineData, a_srcComp, a_dstComp, a_numComp);
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFillPatch::fillInterp(/// interpolated solution on this level
                                      LevelData<FArrayBox>&         a_fineData,
                                      /// solution on coarser level
                                      const LevelData<FArrayBox>&   a_coarseData,
                                      /// starting coarse data component
                                      int                           a_srcComp,
                                      /// starting fine data component
                                      int                           a_dstComp,
                                      /// number of data components to interpolate
                                      int                           a_numComp)
{
  CH_assert(m_defined);

  const Interval srcInterval(a_srcComp, a_srcComp + a_numComp-1);
  const Interval dstInterval(a_dstComp, a_dstComp + a_numComp-1);

  a_coarseData.copyTo(srcInterval,
                      m_spaceInterpolator.coarsenedFineData(),
                      dstInterval);

  // Interpolate to a_fineData from m_spaceInterpolator.coarsenedFineData(),
  // on the given components.
  fillInterpSpaceFromCoarsened(a_fineData, a_srcComp, a_dstComp, a_numComp);
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFillPatch::fillInterpSpaceFromCoarsened(/// interpolated solution on this level
                                                        LevelData<FArrayBox>&         a_fineData,

                                                        /// starting coarse data component
                                                        int                           a_srcComp,
                                                        /// starting fine data component
                                                        int                           a_dstComp,
                                                        /// number of data components to interpolate
                                                        int                           a_numComp)
{
  CH_assert(m_defined);

  const Interval srcInterval(a_srcComp, a_srcComp + a_numComp-1);
  // This should be const, but aliasLevelData doesn't let you do that.
  LevelData<FArrayBox>& coarsenedFineData =
    m_spaceInterpolator.coarsenedFineData();
  LevelData<FArrayBox> coarseCompData;
  aliasLevelData(coarseCompData, &coarsenedFineData, srcInterval);

  const Interval dstInterval(a_dstComp, a_dstComp + a_numComp-1);
  LevelData<FArrayBox> fineCompData;
  aliasLevelData(fineCompData, &a_fineData, dstInterval);

  DataIterator dit = m_layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      // Do not fill in all of fineFab:  fill ghost cells only.
      FArrayBox& fineFab = fineCompData[dit];
      const FArrayBox& coarseFab = coarseCompData[dit];
      const IntVectSet& ivs = m_coarsenedGhosts[dit];

      m_spaceInterpolator.interpOnPatch(fineFab, coarseFab, dit, ivs);
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
  //  m_spaceInterpolator.interpToFine(a_fineData,
  //                                   coarsenedFineData);

  // overwrite interpolated data with valid data when present
  // a_fineData.exchange();
}


//////////////////////////////////////////////////////////////////////////////
TimeInterpolatorRK4& FourthOrderFillPatch::getTimeInterpolator()
{
  CH_assert(m_defined);
  CH_assert(m_timeInterpDefined);
  return m_timeInterpolator;
}

#include "NamespaceFooter.H"
