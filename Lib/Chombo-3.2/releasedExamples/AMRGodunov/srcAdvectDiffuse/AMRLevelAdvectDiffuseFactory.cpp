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

#include "AMRLevelAdvectDiffuseFactory.H"

#include "NamespaceHeader.H"

AMRLevelAdvectDiffuseFactory::
AMRLevelAdvectDiffuseFactory(const AdvectPhysics&        a_gphys,
                             AdvectionVelocityFunction   a_advFunc,
                             BCHolder                    a_bcFunc,
                             const Real&                 a_cfl,
                             const Real&                 a_domainLength,
                             const Real&                 a_refineThresh,
                             const int&                  a_tagBufferSize,
                             const Real&                 a_initialDtMultiplier,
                             const bool&                 a_useLimiting,
                             const Real&                 a_nu)
{
  m_bcFunc                 = a_bcFunc;
  m_cfl                    = a_cfl;
  m_domainLength           = a_domainLength;
  m_refineThresh           = a_refineThresh;
  m_tagBufferSize          = a_tagBufferSize;
  m_initialDtMultiplier    = a_initialDtMultiplier;
  m_useLimiting            = a_useLimiting;
  m_nu                     = a_nu;
  m_advFunc                = a_advFunc;
  GodunovPhysics* gphysPtr = a_gphys.new_godunovPhysics();
  m_advPhys = RefCountedPtr<AdvectPhysics>((AdvectPhysics*)gphysPtr);
}

// Virtual constructor
AMRLevel* AMRLevelAdvectDiffuseFactory::new_amrlevel() const
{
  // Create a new AMRLevelAdvectDiffuse
  AMRLevelAdvectDiffuse* amrGodPtr = new AMRLevelAdvectDiffuse(*m_advPhys,
                                                               m_advFunc,
                                                               m_bcFunc,
                                                               m_cfl,
                                                               m_domainLength,
                                                               m_refineThresh,
                                                               m_tagBufferSize,
                                                               m_initialDtMultiplier,
                                                               m_useLimiting,
                                                               m_nu);

  // Return it
  return (static_cast <AMRLevel*> (amrGodPtr));
}

#include "NamespaceFooter.H"
