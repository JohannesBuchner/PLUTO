#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ModianoIBCFactory.H"
#include "ModianoIBC.H"
#include "PolyGeom.H"

#include "NamespaceHeader.H"

/******************/
/******************/
ModianoIBCFactory::ModianoIBCFactory(const Real&     a_gamma,
                                     const Real&     a_waveAmp,
                                     const Real&     a_waveWidth,
                                     const RealVect& a_center,
                                     const RealVect& a_waveDir,
                                     bool     a_freeStreamOnly,
                                     bool     a_useNegWave)
  :EBPhysIBCFactory()
{
  m_gamma       = a_gamma       ;
  m_waveAmp     = a_waveAmp     ;
  m_waveWidth   = a_waveWidth   ;
  m_waveDir     = a_waveDir     ;
  m_center      = a_center      ;
  m_freeStreamOnly = a_freeStreamOnly;
  m_useNegWave = a_useNegWave;
}

/******************/
/******************/

ModianoIBCFactory::
~ModianoIBCFactory()
{
}
/******************/
/******************/
EBPhysIBC*
ModianoIBCFactory::
create() const
{
  ModianoIBC* retval =
    new ModianoIBC(m_gamma,
                   m_waveAmp,
                   m_waveWidth,
                   m_center,
                   m_waveDir,
                   m_freeStreamOnly,
                   m_useNegWave);
  return static_cast<EBPhysIBC*>(retval);
}
/******************/
/******************/

#include "NamespaceFooter.H"
