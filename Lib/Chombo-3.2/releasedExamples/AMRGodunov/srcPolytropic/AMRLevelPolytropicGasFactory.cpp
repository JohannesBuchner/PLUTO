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

#include "AMRLevelPolytropicGasFactory.H"
#include "AMRLevelPolytropicGas.H"

#include "NamespaceHeader.H"

AMRLevelPolytropicGasFactory::AMRLevelPolytropicGasFactory()
{
  m_godunovPhysics = NULL;
  m_isDefined = false;
}

AMRLevelPolytropicGasFactory::~AMRLevelPolytropicGasFactory()
{
  if (m_godunovPhysics != NULL)
    {
      delete m_godunovPhysics;
      m_godunovPhysics = NULL;
    }

  m_isDefined = false;
}

void AMRLevelPolytropicGasFactory::define(const Real&                 a_cfl,
                                          const Real&                 a_domainLength,
                                          const int&                  a_verbosity,
                                          const Real&                 a_refineThresh,
                                          const int&                  a_tagBufferSize,
                                          const Real&                 a_initialDtMultiplier,
                                          const GodunovPhysics* const a_godunovPhysics,
                                          const int&                  a_normalPredOrder,
                                          const bool&                 a_useFourthOrderSlopes,
                                          const bool&                 a_usePrimLimiting,
                                          const bool&                 a_useCharLimiting,
                                          const bool&                 a_useFlattening,
                                          const bool&                 a_useArtificialViscosity,
                                          const Real&                 a_artificialViscosity,
                                          const bool&                 a_useSourceTerm,
                                          const Real&                 a_sourceTermScaling,
                                          const bool&                 a_highOrderLimiter)
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
  if (m_godunovPhysics != NULL)
    {
      delete m_godunovPhysics;
      m_godunovPhysics = NULL;
    }

  // Store the object that supplies the physics needed by the integrator
  // (used as a factory)
  m_godunovPhysics = a_godunovPhysics->new_godunovPhysics();

  // Store the order of the normal predictor (1 -> PLM, 2 -> PPM)
  m_normalPredOrder = a_normalPredOrder;

  // Store the slope computation parameters
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_usePrimLimiting      = a_usePrimLimiting;
  m_useCharLimiting      = a_useCharLimiting;
  m_useFlattening        = a_useFlattening;

  // Artificial viscosity coefficient must be greater than zero
  CH_assert(!a_useArtificialViscosity || (a_artificialViscosity >= 0.0));

  // Store the artificial viscosity flag and coefficient
  m_useArtificialViscosity = a_useArtificialViscosity;
  m_artificialViscosity    = a_artificialViscosity;

  // Supply a source term to the computation
  m_useSourceTerm     = a_useSourceTerm;
  m_sourceTermScaling = a_sourceTermScaling;

  // Use a high-order limiter?
  m_highOrderLimiter = a_highOrderLimiter;

  // The object is defined
  m_isDefined = true;
}

// Virtual constructor
AMRLevel* AMRLevelPolytropicGasFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelPolytropicGas
  AMRLevelPolytropicGas* amrGodPtr = new AMRLevelPolytropicGas();

  // Define the new object
  amrGodPtr->defineParams(m_cfl,
                          m_domainLength,
                          m_verbosity,
                          m_refineThresh,
                          m_tagBufferSize,
                          m_initialDtMultiplier,
                          m_godunovPhysics,
                          m_normalPredOrder,
                          m_useFourthOrderSlopes,
                          m_usePrimLimiting,
                          m_useCharLimiting,
                          m_useFlattening,
                          m_useArtificialViscosity,
                          m_artificialViscosity,
                          m_useSourceTerm,
                          m_sourceTermScaling,
                          m_highOrderLimiter);

  // Return it
  return (static_cast <AMRLevel*> (amrGodPtr));
}

// Check that everything is defined
bool AMRLevelPolytropicGasFactory::isDefined() const
{
  return m_isDefined;
}

#include "NamespaceFooter.H"
