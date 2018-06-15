#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "OldPhysIBC.H"
#include "UsingNamespace.H"

// Indicate that define() hasn't been called
OldPhysIBC::OldPhysIBC()
{
  m_isDefined = false;
}

// Define the object
void OldPhysIBC::define(const ProblemDomain& a_domain,
                     const Real&          a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx;
  // set default time to be 0
  m_time = 0.0;

  m_isDefined = true;
}

// Default implementation of artificial viscosity at the boundary does nothing
void OldPhysIBC::artViscBC(FArrayBox&       a_F,
                           const FArrayBox& a_U,
                           const FArrayBox& a_divVel,
                           const int&       a_dir,
                           const Real&      a_time)
{
}

//
void OldPhysIBC::time(Real a_time)
{
  m_time = a_time;
}

//
Real OldPhysIBC::time() const
{
  return m_time;
}
