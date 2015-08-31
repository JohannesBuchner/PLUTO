#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "AMRLevelPlutoFactory.H"

#include "NamespaceHeader.H"

AMRLevelPlutoFactory::AMRLevelPlutoFactory()
{
  m_patchPluto = NULL;
  m_isDefined = false;
}

AMRLevelPlutoFactory::~AMRLevelPlutoFactory()
{
  if (m_patchPluto != NULL)
  {
    delete m_patchPluto;
    m_patchPluto = NULL;
  }

  m_isDefined = false;
}

void AMRLevelPlutoFactory::define(const Real&              a_cfl,
                                  const Real&              a_domainLength,
                                  const int&               a_verbosity,
                                  const Real&              a_refineThresh,
                                  const int&               a_tagBufferSize,
                                  const Real&              a_initialDtMultiplier,
                                  const PatchPluto*  const a_patchPluto)
{
  // Store the CFL number
  m_cfl = a_cfl;

  // Store the physical dimension of the longest side of the domain
  m_domainLength = a_domainLength;

  // Store the verbosity of the object
  m_verbosity = a_verbosity;

  // Store the refinement threshold for gradient
  m_refineThresh = a_refineThresh;

  // Store the tag buffer size
  m_tagBufferSize = a_tagBufferSize;

  // Store the initial dt multiplier
  m_initialDtMultiplier = a_initialDtMultiplier;

  // Delete any existing physics object
  if (m_patchPluto != NULL)
  {
    delete m_patchPluto;
    m_patchPluto = NULL;
  }

  // Store the object that supplies the physics needed by the integrator
  // (used as a factory)
  m_patchPluto = a_patchPluto->new_patchPluto();

  // The object is defined
  m_isDefined = true;
}

AMRLevel* AMRLevelPlutoFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelPluto
  AMRLevelPluto* amrGodPtr = new AMRLevelPluto();

  // Set up new object
  amrGodPtr->defineParams(m_cfl,
                          m_domainLength, 
                          m_verbosity,
                          m_refineThresh,
                          m_tagBufferSize,
                          m_initialDtMultiplier,
                          m_patchPluto);

  // Return it
  return (static_cast <AMRLevel*> (amrGodPtr));
}

// Check that everything is defined
bool AMRLevelPlutoFactory::isDefined() const
{
   return m_isDefined;
}

#include "NamespaceFooter.H"
