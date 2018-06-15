#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TensorFineStencilSet.H"
#include "NamespaceHeader.H"

TensorFineStencilSet::TensorFineStencilSet()
  : m_isDefined(false), m_isEmpty(true)
{
}

TensorFineStencilSet::TensorFineStencilSet(const IntVectSet& a_fineIVS,
                                           const ProblemDomain& a_domain,
                                           int a_normalDir)
  : m_isDefined(true), m_isEmpty(false)
{
  define(a_fineIVS, a_domain, a_normalDir);
}



TensorFineStencilSet::~TensorFineStencilSet()
{
}

void
TensorFineStencilSet::define(const IntVectSet& a_fineIVS,
                             const ProblemDomain& a_domain,
                             int a_normalDir)
{
  CH_assert (!a_domain.isEmpty());

  m_normalDir = a_normalDir;

  // set up indexing
  int index = 0;
  for (int tanDir=0; tanDir<SpaceDim; tanDir++)
    {
      if (tanDir != a_normalDir)
        {
          m_tangentialDirections[index] = tanDir;
          index++;
        }
    }



  // are there any cells to deal with here?
  if (!a_fineIVS.isEmpty())
    {
      m_isEmpty = false;

      for (int tanDirIndex=0; tanDirIndex<SpaceDim-1; tanDirIndex++)
        {
          int tanDir = m_tangentialDirections[tanDirIndex];
          m_centeredSet[tanDirIndex] = a_fineIVS;
          // only potentially need to modify anything if domain isn't periodic
          // if periodic in tanDir, then use centered diffs everywhere
          if (!a_domain.isPeriodic(tanDir))
            {
              // this set of operations results in the box over which
              // the centered difference stencil can be used (shrunk
              // by one in non-periodic direction)
              Box centeredStencilBox = a_domain.domainBox();
              centeredStencilBox.grow(tanDir, 1);
              centeredStencilBox &= a_domain;
              centeredStencilBox.grow(tanDir, -1);

              m_centeredSet[tanDirIndex] &= centeredStencilBox;

              Box forwardDiffBox = adjCellLo(a_domain.domainBox(), tanDir);
              forwardDiffBox.shift(tanDir,1);
              m_forwardSet[tanDirIndex] = a_fineIVS;
              m_forwardSet[tanDirIndex] &= forwardDiffBox;

              Box backwardDiffBox = adjCellHi(a_domain.domainBox(),tanDir);
              backwardDiffBox.shift(tanDir,-1);
              m_backwardSet[tanDirIndex] = a_fineIVS;
              m_backwardSet[tanDirIndex] &= backwardDiffBox;
            } // end modifications if domain is not periodic
        } // end loop over tangential directions

    } // end if there is anything to do here

  m_isDefined = true;

}



const IntVectSet &
TensorFineStencilSet::getCenteredStencilSet(int a_tanDir) const
{
  CH_assert (isDefined());
  CH_assert (a_tanDir != m_normalDir);

  int index = -1;

  // probably need to re-think this indexing strategy... (kinda clunky)
  for (int comp=0; comp<SpaceDim-1; comp++)
    if (m_tangentialDirections[comp] == a_tanDir) index = comp;

  return m_centeredSet[index];
}

const IntVectSet&
TensorFineStencilSet::getBackwardStencilSet(int a_tanDir) const
{
  CH_assert (isDefined());
  CH_assert (a_tanDir != m_normalDir);

  int index = -1;

  // probably need to re-think this indexing strategy... (kinda clunky)
  for (int comp=0; comp<SpaceDim-1; comp++)
    if (m_tangentialDirections[comp] == a_tanDir) index = comp;

  return m_backwardSet[index];
}

const IntVectSet&
TensorFineStencilSet::getForwardStencilSet(int a_tanDir) const
{
  CH_assert (isDefined());
  CH_assert (a_tanDir != m_normalDir);

  int index = -1;

  // probably need to re-think this indexing strategy... (kinda clunky)
  for (int comp=0; comp<SpaceDim-1; comp++)
    if (m_tangentialDirections[comp] == a_tanDir) index = comp;

  return m_forwardSet[index];
}

bool
TensorFineStencilSet::isDefined() const
{
  return m_isDefined;
}


bool
TensorFineStencilSet::isEmpty() const
{
  return m_isEmpty;
}
#include "NamespaceFooter.H"
