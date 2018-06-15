#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELPOLYTROPICGASFACTORY_H_
#define _AMRLEVELPOLYTROPICGASFACTORY_H_

#include "AMRLevelFactory.H"

#include "AMRLevelPolytropicGas.H"
#include "NamespaceHeader.H"

/// AMR Godunov factory
/**
 */
class AMRLevelPolytropicGasFactory: public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelPolytropicGasFactory();

  /// Destructor
  /**
   */
  virtual ~AMRLevelPolytropicGasFactory();

  /// Define the parameters of the factory
  /**
   */
  virtual void define(/// CFL number
                      const Real&                 a_cfl,
                      /// physical dimension of the longest side of the domain
                      const Real&                 a_domainLength,
                      /// verbosity level
                      const int&                  a_verbosity,
                      /// undivided gradient size over which a cell will be tagged for refinement
                      const Real&                 a_refineThresh,
                      /// number of buffer cells around each tagged cell that will also be tagged
                      const int&                  a_tagBufferSize,
                      /// CFL number at beginning of calculation
                      const Real&                 a_initialDtMultiplier,
                      /// physics class
                      const GodunovPhysics* const a_godunovPhysics,
                      /// order of the normal predictor: 1 for PLM, 2 for PPM
                      const int&                  a_normalPredOrder,
                      /// whether to use 4th-order slopes
                      const bool&                 a_useFourthOrderSlopes,
                      /// whether to apply slope limiting to primitive variables
                      const bool&                 a_usePrimLimiting,
                      /// whether to apply slope limiting to characteristic variables
                      const bool&                 a_useCharLimiting,
                      /// whether to apply slope flattening
                      const bool&                 a_useFlattening,
                      /// whether to apply artificial viscosity
                      const bool&                 a_useArtificialViscosity,
                      /// artificial viscosity coefficient
                      const Real&                 a_artificialViscosity,
                      /// whether to apply source term
                      const bool&                 a_useSourceTerm,
                      /// scaling of source term
                      const Real&                 a_sourceTermScaling,
                      /// whether to apply 4th-order limiter
                      const bool&                 a_highOrderLimiter);

  /// Create a new AMRLevel (for polytropic gases)
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Check that this object has been defined
  /**
   */
  bool isDefined() const;

protected:
  // CFL number
  Real m_cfl;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;

  // Verbosity of the object
  int m_verbosity;

  // Refinement threshold for gradient
  Real m_refineThresh;

  // Tag buffer size
  int m_tagBufferSize;

  // Initial dt multiplier
  Real m_initialDtMultiplier;

  // Pointer to the object that supplies the physics needed by the integrator
  // (used as a factory here)
  const GodunovPhysics* m_godunovPhysics;

  // Order of the normal predictor (1 -> PLM, 2 -> PPM)
  int m_normalPredOrder;

  // Use 4th other slope computations (otherwise, use 2nd order)
  bool m_useFourthOrderSlopes;

  // Do slope limiting in the primitive or characteristic variables, respect.
  bool m_usePrimLimiting;
  bool m_useCharLimiting;

  // Do slope flattening - MUST BE USING 4th order slopes
  bool m_useFlattening;

  // Apply artificial viscosity of a set value
  bool m_useArtificialViscosity;
  Real m_artificialViscosity;

  // Supply a source term to the computation
  bool m_useSourceTerm;
  Real m_sourceTermScaling;

  // Use a high-order limiter?
  bool m_highOrderLimiter;

  // Has this object been defined
  bool m_isDefined;

private:

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelPolytropicGasFactory&);
  AMRLevelPolytropicGasFactory(const AMRLevelPolytropicGasFactory&);
};

#include "NamespaceFooter.H"

#endif
