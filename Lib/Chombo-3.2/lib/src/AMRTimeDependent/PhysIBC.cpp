#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PhysIBC.H"
#include "NamespaceHeader.H"

// Indicate that define() hasn't been called
PhysIBC::PhysIBC()
{
  m_isDefined = false;
}

PhysIBC::~PhysIBC()
{
}

// Define the object
void PhysIBC::define(const ProblemDomain& a_domain,
                     const Real&          a_dx)
{
  m_domain = a_domain;
  m_dx     = a_dx;

  m_isDefined = true;
}

void PhysIBC::getBoundaryFaces(Box&                   a_boundaryBox,
                               const Box&             a_dataFaceBox,
                               const int&             a_dir,
                               const Side::LoHiSide&  a_side)
{
  CH_assert(m_isDefined == true);
  // Return the empty box by default.
  a_boundaryBox = Box();
  if (!m_domain.isPeriodic(a_dir))
    {
      // domainFaces is face-centered in direction a_dir
      Box domainFaces = m_domain.domainBox();
      domainFaces.surroundingNodes(a_dir);
      Box bx = a_dataFaceBox;
      // We need to do this in case a_dataFaceBox has ghosts outside m_domain
      bx &= domainFaces;
      // Determine which side and thus shifting directions
      int lohisign = sign(a_side);
      bx.shiftHalf(a_dir, lohisign);
      // Is there a domain boundary next to this grid?
      if (!m_domain.contains(bx))
        {
          bx &= m_domain;
          if (a_side == Side::Lo)
            {
              a_boundaryBox = bdryLo(bx, a_dir);
            }
          else // (a_side == Side::Hi)
            {
              a_boundaryBox = bdryHi(bx, a_dir);
            }
        }
    }
}

#include "NamespaceFooter.H"
