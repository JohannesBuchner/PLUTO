#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "VelBCHolder.H"

// -------------------------------------------------------------
VelBCHolder::VelBCHolder()
{
}

// --------------------------------------------------------------
VelBCHolder::~VelBCHolder()
{
}

// --------------------------------------------------------------
VelBCHolder::VelBCHolder(const Tuple<BCHolder, SpaceDim>& a_componentBC)
{
  for (int idir=0; idir<SpaceDim; idir++)
    {
      m_componentBC[idir] = a_componentBC[idir];
    }
}

// --------------------------------------------------------------
void
VelBCHolder::applyBCs(LevelData<FArrayBox>& a_state,
                      const DisjointBoxLayout& a_valid,
                      const ProblemDomain& a_domain,
                      const Real a_dx,
                      bool a_homogeneous)
{
  DataIterator dit = a_state.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& stateFab = a_state[dit];
      const Box& bx = a_valid[dit];
      for (int idir=0; idir<SpaceDim; idir++)
        {
          m_componentBC[idir](stateFab, bx,
                              a_domain, a_dx,
                              a_homogeneous);
        }
    }
}

// -------------------------------------------------------------
EdgeVelBCHolder::EdgeVelBCHolder()
{
}

// --------------------------------------------------------------
EdgeVelBCHolder::~EdgeVelBCHolder()
{
}

// --------------------------------------------------------------
EdgeVelBCHolder::EdgeVelBCHolder(const Tuple<BCHolder, SpaceDim>& a_componentBC)
{
  for (int idir=0; idir<SpaceDim; idir++)
    {
      m_componentBC[idir] = a_componentBC[idir];
    }
}

// --------------------------------------------------------------
void
EdgeVelBCHolder::applyBCs(LevelData<FluxBox>& a_state,
                          const DisjointBoxLayout& a_valid,
                          const ProblemDomain& a_domain,
                          const Real a_dx,
                          bool a_homogeneous)
{
  DataIterator dit = a_state.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FluxBox& stateFab = a_state[dit];
      // bx is CELL-centered
      const Box& bx = a_valid[dit];
      for (int idir=0; idir<SpaceDim; idir++)
        {
          // stateFab[idir] is FACE-centered on face idir
          m_componentBC[idir](stateFab[idir], bx,
                              a_domain, a_dx,
                              a_homogeneous);
        }
    }
}
