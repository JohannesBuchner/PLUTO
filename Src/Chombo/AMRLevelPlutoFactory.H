#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _AMR_LEVEL_GODUNOV_FACTORY_H_
#define _AMR_LEVEL_GODUNOV_FACTORY_H_

#include "AMRLevelFactory.H"
#include "AMRLevelPluto.H"

#include "NamespaceHeader.H"

/// AMR Pluto factory
/**
 */
class AMRLevelPlutoFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  AMRLevelPlutoFactory();

  /// Destructor
  virtual ~AMRLevelPlutoFactory();

  virtual void define(const Real&              a_cfl,
                      const Real&              a_domainLength,
                      const int&               a_verbosity,
                      const Real&              a_refineThresh,
                      const int&               a_tagBufferSize,
                      const Real&              a_initialDtMultiplier,
                      const PatchPluto*  const a_patchPluto);

  /// Create a new AMRLevel
  virtual AMRLevel* new_amrlevel() const;

  /// Check that everything is defined
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
  int  m_tagBufferSize;

  // Initial dt multiplier 
  Real m_initialDtMultiplier;

  // Patch integrator (used as a factory)
  const PatchPluto* m_patchPluto;

  // Has this object been defined
  bool m_isDefined;


private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelPlutoFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelPlutoFactory(const AMRLevelPlutoFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"

#endif
