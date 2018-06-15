#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBPatch_Polytropic_Factory_H_
#define _EBPatch_Polytropic_Factory_H_

#include "EBPatchPolytropic.H"
#include "EBPatchGodunovFactory.H"

#include "NamespaceHeader.H"

///
/**
 */
class EBPatchPolytropicFactory: public EBPatchGodunovFactory
{
public:
  //! Destructor.
  virtual ~EBPatchPolytropicFactory();

  ///
  /**
   */
  EBPatchGodunov* create() const;

  //! Create a new factory, specifying boundary conditions, gas properties, 
  //! and algorithmic parameters.
  //! \param a_bc The boundary condition factory for the desired problem.
  //! \param a_gamma The ratio of specific heats (Cp/Cv) for the ideal gas.
  //! \param a_specHeat The specific heat of the gas, measured at constant pressure.
  //! \param a_useFourthOrderSlopes If set, the EBPatchPolytropic objects created 
  //!                               by this factory will use fourth-order estimates
  //!                               for slopes.
  EBPatchPolytropicFactory(const EBPhysIBCFactory* const a_bc,
                           const Real&                   a_gamma,
                           const Real&                   a_specHeat,
                           const bool&                   a_useFourthOrderSlopes,
                           const bool&                   a_useFlattening,
                           const bool&                   a_useArtificialVisc,
                           const bool&                   a_useLimiting,
                           const bool&                   a_doRZCoords);

  //! Create a new factory, specifying boundary conditions, gas properties, 
  //! and algorithmic parameters. This version omits the specific heat, 
  //! setting it to 1.
  //! \param a_bc The boundary condition factory for the desired problem.
  //! \param a_gamma The ratio of specific heats (Cp/Cv) for the ideal gas.
  //! \param a_specHeat The specific heat of the gas, measured at constant pressure.
  //! \param a_useFourthOrderSlopes If set, the EBPatchPolytropic objects created 
  //!                               by this factory will use fourth-order estimates
  //!                               for slopes.
  EBPatchPolytropicFactory(const EBPhysIBCFactory* const a_bc,
                           const Real&                   a_gamma,
                           const bool&                   a_useFourthOrderSlopes,
                           const bool&                   a_useFlattening,
                           const bool&                   a_useArtificialVisc,
                           const bool&                   a_useLimiting,
                           const bool&                   a_doRZCoords);

protected:
  const EBPhysIBCFactory* m_bcFactoryPtr;
  Real              m_gamma;
  Real              m_specHeat;
  bool              m_useFourthOrderSlopes;
  bool              m_useFlattening;
  bool              m_useArtificialVisc;
  bool              m_doRZCoords;
  bool              m_useLimiting;

private:

  //disallowed for all the usual reasons
  EBPatchPolytropicFactory()
  {
    MayDay::Error("invalid operator");
  }
  void operator=(const EBPatchPolytropicFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchPolytropicFactory(const EBPatchPolytropicFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
