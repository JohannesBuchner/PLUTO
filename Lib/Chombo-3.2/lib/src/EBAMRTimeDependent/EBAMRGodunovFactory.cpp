#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRLevel.H"
#include "EBAMRGodunov.H"
#include "EBAMRGodunovFactory.H"
#include "NamespaceHeader.H"
/****************/
/****************/
AMRLevel*
EBAMRGodunovFactory::
new_amrlevel() const
{
  EBAMRGodunov* amrg_ptr = new EBAMRGodunov();

  amrg_ptr->CFL(m_cfl);
  amrg_ptr->patchGodunov(m_patchGodunov);
  amrg_ptr->refinementThreshold(m_refineThresh);
  amrg_ptr->tagBufferSize(m_tagBufferSize);
  amrg_ptr->initialDtMultiplier(m_initialDtMultiplier);
  amrg_ptr->domainLength(m_domainLength);
  amrg_ptr->redistRadius(m_redistRad);
  amrg_ptr->verbosity(m_verbosity);
  amrg_ptr->useMassRedistribution(m_useMassRedist);
  amrg_ptr->doRZCoords(m_doRZCoords);
  amrg_ptr->hasSourceTerm(m_hasSourceTerm);
  amrg_ptr->doSmushing(m_doSmushing);
  amrg_ptr->tagAll(m_tagAll);

  return (static_cast <AMRLevel*> (amrg_ptr));
}
/****************/
/****************/
EBAMRGodunovFactory::
~EBAMRGodunovFactory()
{
}
/****************/
/****************/
EBAMRGodunovFactory::
EBAMRGodunovFactory(const Real&                   a_initialDtMultiplier,
                    const Real&                   a_cfl,
                    const int &                   a_redistRad,
                    const RealVect&               a_domainLength,
                    const Real&                   a_refineThresh,
                    const int &                   a_tagBufferSize,
                    const int &                   a_verbosity,
                    const bool&                   a_useMassRedist,
                    const bool&                   a_doSmushing,
                    const bool&                   a_doRZCoords,
                    const bool&                   a_hasSourceTerm,
                    const EBPatchGodunovFactory* const a_patchGodunov,
                    bool                          a_tagAll)
{
  m_tagAll = a_tagAll;

  m_initialDtMultiplier = a_initialDtMultiplier;
  m_cfl                 = a_cfl;
  m_useMassRedist       = a_useMassRedist;
  m_redistRad           = a_redistRad;
  m_domainLength        = a_domainLength;
  m_refineThresh        = a_refineThresh;
  m_tagBufferSize       = a_tagBufferSize;
  m_verbosity           = a_verbosity;
  m_doSmushing          = a_doSmushing;
  m_doRZCoords          = a_doRZCoords;
  m_hasSourceTerm       = a_hasSourceTerm;
  m_patchGodunov        = a_patchGodunov;
}
/****************/
/****************/
#include "NamespaceFooter.H"
