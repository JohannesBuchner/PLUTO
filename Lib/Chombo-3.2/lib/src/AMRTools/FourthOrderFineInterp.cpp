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

#include "FourthOrderFineInterp.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderFineInterp::FourthOrderFineInterp()
{
  m_defined = false;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
FourthOrderFineInterp::~FourthOrderFineInterp()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void FourthOrderFineInterp::define(/// layout at this level
                                   const DisjointBoxLayout&  a_layout,
                                   /// number of variables
                                   const int&                a_numStates,
                                   /// refinement ratio between this level and next coarser level
                                   const int&                a_refineCoarse,
                                   /// problem domain on this level
                                   const ProblemDomain&      a_domain,
                                   /// width of coarse ghost cell layer with fine cells to be filled in
                                   const int&                a_coarseGhostsFill,
                                   /// dimensions that are fixed, not interpolated
                                   Interval                  a_fixedDims)
{
  // Cache data
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_numStates = a_numStates;
  m_coarseGhostsFill = a_coarseGhostsFill;
  m_fixedDims = a_fixedDims;
  // m_ghostVect = a_ghostVect;

  IntVect interpUnit = IntVect::Unit;
  m_refineVect = m_refineCoarse * IntVect::Unit;
  for (int dirf = m_fixedDims.begin(); dirf <= m_fixedDims.end(); dirf++)
    {
      interpUnit[dirf] = 0;
      m_refineVect[dirf] = 1;
    }

  m_layout = a_layout;
  // petermc, 19 Dec 2008:  prevents crash in case of calling define again
  m_layoutCoarsened = DisjointBoxLayout();
  coarsen(m_layoutCoarsened, m_layout, m_refineVect);

  m_maxStencilDist = 2;
  m_patchInterp.define(m_domain, m_refineCoarse, m_maxStencilDist, m_fixedDims);

  // Now find m_stencilHere:
  // For each coarse cell whose fine cells will be filled in,
  // figure out which stencil to use.
  IntVect coarseGhostVect = m_coarseGhostsFill * interpUnit;
  // LevelData< BaseFab<IntVect> > m_stencilHere;
  m_stencilHere.define(m_layoutCoarsened, 1, coarseGhostVect);

  DataIterator dit = a_layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& coarseBox = m_layoutCoarsened[dit];
      // petermc, 9 Jan 2009:  grow by coarseGhostVect for filling in ghosts
      Box coarseGhostedBox = grow(coarseBox, coarseGhostVect);
      m_patchInterp.setCoarseBox(coarseGhostedBox);

      BaseFab<IntVect>& stencilHereFab = m_stencilHere[dit];
      m_patchInterp.setStencil(stencilHereFab);
    }

  // layers of ghost cells needed to interpolate to fine data,
  // which have m_coarseGhostsFill layers of ghost cells
  m_ghostVect = (m_coarseGhostsFill + m_maxStencilDist) * interpUnit;
  // used only in interpToFine function
  m_coarseData.define(m_layoutCoarsened, m_numStates, m_ghostVect);

  // Everything is defined now.
  m_defined = true;
}


//////////////////////////////////////////////////////////////////////////////
const BaseFab<IntVect>& FourthOrderFineInterp::getStencil(const DataIterator&    a_dit) const
{
  return getStencil(a_dit());
}


//////////////////////////////////////////////////////////////////////////////
const BaseFab<IntVect>& FourthOrderFineInterp::getStencil(const DataIndex&    a_dind) const
{
  CH_assert(m_defined);
  return m_stencilHere[a_dind];
}


//////////////////////////////////////////////////////////////////////////////
LevelData<FArrayBox>& FourthOrderFineInterp::coarsenedFineData()
{
  CH_assert(m_defined);
  return m_coarseData;
}


//////////////////////////////////////////////////////////////////////////////
FourthOrderPatchInterp& FourthOrderFineInterp::patchInterp()
{
  CH_assert(m_defined);
  return m_patchInterp;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFineInterp::interpToFine(
                                         /// interpolated solution on this level
                                         LevelData<FArrayBox>&         a_fine,
                                         /// coarse solution
                                         const LevelData<FArrayBox>&   a_coarse)
{
  CH_assert(m_defined);

  // Simple copyTo also fills in ghost cells of m_coarseData.
  a_coarse.copyTo(m_coarseData);

  DataIterator dit = m_layoutCoarsened.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& coarseBox = m_layoutCoarsened[dit];
      const FArrayBox& coarseFab = m_coarseData[dit];
      const BaseFab<IntVect>& stencilHereFab = m_stencilHere[dit];
      FArrayBox& fineFab = a_fine[dit];

      m_patchInterp.setCoarseBox(coarseBox);
      m_patchInterp.interpToFine(fineFab, coarseFab, stencilHereFab);
    }
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFineInterp::interpOnPatch(
                                          /// interpolated solution on this level
                                          FArrayBox&            a_fine,
                                          /// coarse solution
                                          const FArrayBox&      a_coarse,
                                          /// index
                                          const DataIterator&   a_dit)
{
  interpOnPatch(a_fine, a_coarse, a_dit());
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFineInterp::interpOnPatch(
                                          /// interpolated solution on this level
                                          FArrayBox&            a_fine,
                                          /// coarse solution
                                          const FArrayBox&      a_coarse,
                                          /// index
                                          const DataIndex&      a_dind)
{
  CH_assert(m_defined);

  const Box& coarseBox = m_layoutCoarsened[a_dind];
  const BaseFab<IntVect>& stencilHereFab = m_stencilHere[a_dind];

  m_patchInterp.setCoarseBox(coarseBox);
  m_patchInterp.interpToFine(a_fine, a_coarse, stencilHereFab);
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFineInterp::interpOnPatch(
                                          /// interpolated solution on this level
                                          FArrayBox&            a_fine,
                                          /// coarse solution
                                          const FArrayBox&      a_coarse,
                                          /// index
                                          const DataIterator&   a_dit,
                                          /// we fill in fine cells within these coarse cells
                                          const IntVectSet&     a_ivs)
{
  interpOnPatch(a_fine, a_coarse, a_dit(), a_ivs);
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderFineInterp::interpOnPatch(
                                          /// interpolated solution on this level
                                          FArrayBox&            a_fine,
                                          /// coarse solution
                                          const FArrayBox&      a_coarse,
                                          /// index
                                          const DataIndex&      a_dind,
                                          /// we fill in fine cells within these coarse cells
                                          const IntVectSet&     a_ivs)
{
  CH_assert(m_defined);

  // const Box& coarseBox = m_layoutCoarsened[a_dind];
  const BaseFab<IntVect>& stencilHereFab = m_stencilHere[a_dind];

  // not needed
  // m_patchInterp.setCoarseBox(coarseBox);
  m_patchInterp.interpToFine(a_fine, a_coarse, stencilHereFab, a_ivs);
}
