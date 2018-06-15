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

#include "FourthOrderPatchInterp.H"
#include "BoxIterator.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderPatchInterp::FourthOrderPatchInterp()
{
  m_defined = false;
  m_isCoarseBoxSet = false;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
FourthOrderPatchInterp::~FourthOrderPatchInterp()
{
  if (m_defined)
    {
      const Box& stencilBox = m_stencils.box();
      for (BoxIterator bit(stencilBox); bit.ok(); ++bit)
        {
          IntVect offset = bit();
          delete m_stencils(offset, 0);
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void FourthOrderPatchInterp::define(
                                    /// problem domain on this level
                                    const ProblemDomain&      a_domain,
                                    /// refinement ratio between this level and next coarser level
                                    const int&                a_refineCoarse,
                                    /// maximum distance of stencil from domain boundary
                                    const int&        a_maxStencilDist,
                                    /// dimensions that are fixed, not interpolated
                                    Interval               a_fixedDims)
{
  if (m_defined)
    {
      const Box& stencilBox = m_stencils.box();
      for (BoxIterator bit(stencilBox); bit.ok(); ++bit)
        {
          IntVect offset = bit();
          delete m_stencils(offset, 0);
        }
      m_defined = false;
    } 
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_maxStencilDist = a_maxStencilDist;
  m_fixedDims = a_fixedDims;

  IntVect interpUnit = IntVect::Unit;
  m_refineVect = m_refineCoarse * IntVect::Unit;
  for (int dirf = m_fixedDims.begin(); dirf <= m_fixedDims.end(); dirf++)
    {
      interpUnit[dirf] = 0;
      m_refineVect[dirf] = 1;
    }

  m_coarseDomain = coarsen(m_domain, m_refineVect);

  m_degree = 3;
  Box stencilBox(-m_maxStencilDist*interpUnit,
                 m_maxStencilDist*interpUnit);
  m_stencils.define(stencilBox, 1);
  for (BoxIterator bit(stencilBox); bit.ok(); ++bit)
    {
      IntVect offset = bit();
      m_stencils(offset, 0) =
        new FourthOrderInterpStencil(offset, m_refineCoarse, m_degree, m_fixedDims);
    }

  // Everything is defined now.
  m_defined = true;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderPatchInterp::setCoarseBox(const Box& a_coarseBox)
{
  m_coarseBox      = a_coarseBox;
  m_isCoarseBoxSet = true;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderPatchInterp::setStencil(BaseFab<IntVect>&  a_stencil)
{
  CH_assert(m_defined);
  CH_assert(m_isCoarseBoxSet);

  const Box& coarseDomainBox = m_coarseDomain.domainBox();
  const IntVect& coarseDomainLo = coarseDomainBox.smallEnd();
  const IntVect& coarseDomainHi = coarseDomainBox.bigEnd();
  for (BoxIterator bit(m_coarseBox); bit.ok(); ++bit)
    {
      IntVect ivc = bit();
      // Set IntVect stencilHereFab(ivc, 0).
      // Find distance to coarseDomain boundary
      IntVect dist = IntVect::Zero;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (! m_fixedDims.contains(idir))
            {
              if (! m_coarseDomain.isPeriodic(idir))
                {
                  int offLo = coarseDomainLo[idir] - ivc[idir] - 1;
                  int offHi = coarseDomainHi[idir] - ivc[idir] + 1;
                  if (offLo < 0 && offHi >0) // condition means ivc in coarseDomain
                    {
                      if ((offLo >= -m_maxStencilDist) &&
                          (offHi <= m_maxStencilDist))
                        { // both of these:  very narrow domain, you are in trouble
                          MayDay::Error("FourthOrderFineInterp::define bad boxes");
                        }
                      if (offLo >= -m_maxStencilDist) // -1 or -2
                        dist[idir] = offLo;
                      if (offHi <= m_maxStencilDist) // 1 or 2
                        dist[idir] = offHi;
                      // Otherwise, dist[idir] = 0.
                    }
                }
            }
        }
      a_stencil(ivc, 0) = dist;
    }
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderPatchInterp::interpToFine(/// interpolated solution on this level
                                          FArrayBox&                a_fine,
                                          /// coarse solution
                                          const FArrayBox&          a_coarse,
                                          /// stencils
                                          const BaseFab<IntVect>&   a_stencils)
{
  CH_assert(m_defined);
  CH_assert(m_isCoarseBoxSet);
  for (BoxIterator bit(m_coarseBox); bit.ok(); ++bit)
    {
      IntVect ivc = bit();
      const IntVect& stencilIndex = a_stencils(ivc, 0);
      const FourthOrderInterpStencil& stencil =
        *m_stencils(stencilIndex, 0);
      // Using coarseFab, fill fine cells of fineFab within ivc.
      stencil.fillFine(a_fine, a_coarse, ivc);
    }
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderPatchInterp::interpToFine(/// interpolated solution on this level
                                          FArrayBox&                a_fine,
                                          /// coarse solution
                                          const FArrayBox&          a_coarse,
                                          /// stencils
                                          const BaseFab<IntVect>&   a_stencils,
                                          /// we fill in fine cells within these coarse cells
                                          const IntVectSet&         a_ivs)
{
  CH_assert(m_defined);
  // CH_assert(m_isCoarseBoxSet);
  for (IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
    {
      const IntVect& ivc = ivsit();
      const IntVect& stencilIndex = a_stencils(ivc, 0);
      const FourthOrderInterpStencil& stencil =
        *m_stencils(stencilIndex, 0);
      // Using coarseFab, fill fine cells of fineFab within ivc.
      stencil.fillFine(a_fine, a_coarse, ivc);
    }
}
