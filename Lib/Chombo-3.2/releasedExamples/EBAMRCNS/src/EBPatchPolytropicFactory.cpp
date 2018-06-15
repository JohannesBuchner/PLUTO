#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBPatchPolytropicFactory.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
EBPatchPolytropicFactory::
EBPatchPolytropicFactory(const EBPhysIBCFactory* const a_bcFactoryPtr,
                         const Real&                   a_gamma,
                         const Real&                   a_specHeat,
                         const bool&                   a_useFourthOrderSlopes,
                         const bool&                   a_useFlattening,
                         const bool&                   a_useArtificialVisc,
                         const bool&                   a_useLimiting,
                         const bool&                   a_doRZCoords)
  :EBPatchGodunovFactory()
{
  m_useLimiting          = a_useLimiting;
  m_gamma                = a_gamma;
  m_specHeat             = a_specHeat;
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_bcFactoryPtr         = a_bcFactoryPtr;
  m_useFlattening        = a_useFlattening;
  m_useArtificialVisc    = a_useArtificialVisc;
  m_doRZCoords           = a_doRZCoords;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBPatchPolytropicFactory::
EBPatchPolytropicFactory(const EBPhysIBCFactory* const a_bcFactoryPtr,
                         const Real&                   a_gamma,
                         const bool&                   a_useFourthOrderSlopes,
                         const bool&                   a_useFlattening,
                         const bool&                   a_useArtificialVisc,
                         const bool&                   a_useLimiting,
                         const bool&                   a_doRZCoords)
  :EBPatchGodunovFactory()
{
  m_useLimiting          = a_useLimiting;
  m_gamma                = a_gamma;
  m_specHeat             = 1.0;
  m_useFourthOrderSlopes = a_useFourthOrderSlopes;
  m_bcFactoryPtr         = a_bcFactoryPtr;
  m_useFlattening        = a_useFlattening;
  m_useArtificialVisc    = a_useArtificialVisc;
  m_doRZCoords           = a_doRZCoords;
}
//-----------------------------------------------------------------------

/******************/
EBPatchPolytropicFactory::
~EBPatchPolytropicFactory()
{;}
/******************/
/******************/
EBPatchGodunov*
EBPatchPolytropicFactory::
create() const
{
  EBPatchPolytropic* retval= new EBPatchPolytropic();
  retval->setSlopeParameters(m_useFourthOrderSlopes,
                             m_useFlattening,
                             m_useLimiting);
  retval->artificialViscosity(m_useArtificialVisc);
  retval->setEBPhysIBC(*m_bcFactoryPtr);
  retval->setGamma(m_gamma);
  retval->setSpecHeat(m_specHeat);
  retval->doRZCoords(m_doRZCoords);
  return static_cast<EBPatchGodunov*>(retval);
}
/******************/
/******************/

#include "NamespaceFooter.H"
