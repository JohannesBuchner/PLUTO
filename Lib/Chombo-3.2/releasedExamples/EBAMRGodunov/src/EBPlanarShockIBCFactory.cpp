#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBPlanarShockIBCFactory.H"
#include "EBPlanarShockIBC.H"

#include "NamespaceHeader.H"

/******************/
/******************/
EBPlanarShockIBCFactory::
EBPlanarShockIBCFactory(const Real&     a_gamma,
                        const Real&     a_ms,
                        const Real&     a_center,
                        const int&      a_shocknorm,
                        const bool&     a_shockbackward,
                        const bool&     a_doRZCoords)
  :EBPhysIBCFactory()
{
  if (a_gamma < 1.0)
    {
      MayDay::Error("bogus gamma in EBPlanarShock");
    }
  if (a_ms < 1.0)
    {
      MayDay::Error("EBPlanarShock can only handle mach numbers greater than 1");
    }
  m_gamma         = a_gamma;
  m_ms            = a_ms;
  m_center        = a_center;
  m_shocknorm     = a_shocknorm;
  m_shockbackward = a_shockbackward;
  m_doRZCoords    = a_doRZCoords;
}
/******************/
/******************/

EBPlanarShockIBCFactory::
~EBPlanarShockIBCFactory()
{
}
/******************/
/******************/
EBPhysIBC*
EBPlanarShockIBCFactory::
create() const
{
  EBPlanarShockIBC* retval =
    new EBPlanarShockIBC(m_gamma,
                         m_ms,
                         m_center,
                         m_shocknorm,
                         m_shockbackward,
                         m_doRZCoords);
  return static_cast<EBPhysIBC*>(retval);
}
/******************/
/******************/

#include "NamespaceFooter.H"
