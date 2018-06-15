#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBAMRCNSFACTORY_H_
#define _EBAMRCNSFACTORY_H_

#include "EBAMRCNS.H"
#include "AMRLevelFactory.H"
#include "EBPatchGodunovFactory.H"
#include "EBPatchPolytropicFactory.H"

#include "NamespaceHeader.H"

//! \class EBAMRCNSFactory
//! This factory creates EBAMRCNS compressible Navier-Stokes level solvers.
class EBAMRCNSFactory: public AMRLevelFactory
{

  public:
  //! 
  AMRLevel* new_amrlevel() const;

  //! Create a new factory that assembles EBAMRCNS level solvers.
  //! \param a_params The parameters for the level solvers.
  //! \param a_factory A factory that builds EBPatchGodunov objects for the solvers.
  //! \param a_initialConditions A space-time function that specifies initial conditions.
  EBAMRCNSFactory(const EBAMRCNSParams& a_params,
                  const RefCountedPtr<EBPatchPolytropicFactory>& a_factory,
                  const RefCountedPtr<EBSpaceTimeFunction>& a_initialConditions):
    m_params(a_params),
    m_factory(a_factory),
    m_ICs(a_initialConditions)
  {
  }

  //! Create a new factory that assembles EBAMRCNS level solvers. This 
  //! constructor supports the older-style specification of initial conditions
  //! via boundary conditions. This will probably die soon--use *ONLY* if you 
  //! know what you're doing!!!
  //! \param a_params The parameters for the level solvers.
  //! \param a_factory A factory that builds EBPatchGodunov objects for the solvers.
  EBAMRCNSFactory(const EBAMRCNSParams& a_params,
                  const RefCountedPtr<EBPatchPolytropicFactory>& a_factory):
    m_params(a_params),
    m_factory(a_factory),
    m_ICs(new EBSpaceTimeFunctionIBCAdaptor())
  {
  }

  ///
  virtual ~EBAMRCNSFactory()
  {;}

  protected:

  EBAMRCNSParams m_params;
  RefCountedPtr<EBPatchPolytropicFactory> m_factory;
  RefCountedPtr<EBSpaceTimeFunction> m_ICs;

  private:

  // Disallowed operators.
  EBAMRCNSFactory();

  // The compiler generates the copy constructor and assignment operator.
};

#include "NamespaceFooter.H"
#endif
