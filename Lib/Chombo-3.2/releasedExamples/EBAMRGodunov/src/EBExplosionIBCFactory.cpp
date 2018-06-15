#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBExplosionIBCFactory.H"
#include "EBExplosionIBC.H"

#include "NamespaceHeader.H"

/******************/
/******************/
EBExplosionIBCFactory::
EBExplosionIBCFactory(const Real&     a_gamma,
                      const Real&     a_size,
                      const Real&     a_p0,
                      const Real&     a_r0,
                      const Real&     a_p1,
                      const Real&     a_r1,
                      const RealVect& a_center,
                      const bool&     a_doRZCoords)
  :EBPhysIBCFactory()
{
  m_doRZCoords = a_doRZCoords;
  m_gamma   = a_gamma;
  m_size    = a_size;
  m_p0      = a_p0;
  m_r0      = a_r0;
  m_p1      = a_p1;
  m_r1      = a_r1;
  m_center  = a_center;
}
/******************/
/******************/

EBExplosionIBCFactory::
~EBExplosionIBCFactory()
{
}
/******************/
/******************/
EBPhysIBC*
EBExplosionIBCFactory::
create() const
{
  EBExplosionIBC* retval =
    new EBExplosionIBC(m_gamma,
                       m_size,
                       m_p0,
                       m_r0,
                       m_p1,
                       m_r1,
                       m_center,
                       m_doRZCoords);
  return static_cast<EBPhysIBC*>(retval);
}
/******************/
/******************/

#include "NamespaceFooter.H"
