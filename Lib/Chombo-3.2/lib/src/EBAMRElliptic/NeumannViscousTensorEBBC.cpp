#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "BoxIterator.H"
#include "VoFIterator.H"

#include "EBArith.H"

#include "NeumannViscousTensorEBBC.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"


NeumannViscousTensorEBBC::
NeumannViscousTensorEBBC(
                         const ProblemDomain& a_domain,
                         const EBISLayout&    a_layout,
                         const RealVect&      a_dx,
                         const IntVect*       a_ghostCellsPhi,
                         const IntVect*       a_ghostCellsRhs)
{
}

NeumannViscousTensorEBBC::
~NeumannViscousTensorEBBC()
{
}

void
NeumannViscousTensorEBBC::
applyEBFlux(EBCellFAB&                    a_lphi,
            const EBCellFAB&              a_phi,
            VoFIterator&                  a_vofit,
            const LayoutData<IntVectSet>& a_cfivs,
            const DataIndex&              a_dit,
            const RealVect&               a_probLo,
            const RealVect&               a_dx,
            const Real&                   a_factor,
            const bool&                   a_useHomogeneous,
            const Real&                   a_time)
{
  CH_assert(m_coefSet);
  const EBISBox&   ebisBox = a_phi.getEBISBox();
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      const VolIndex& vof = a_vofit();
      Real grad[SpaceDim][SpaceDim];
      getBoundaryGrad(grad, vof, a_dx, a_probLo, ebisBox);

      Real flux[SpaceDim][SpaceDim];

      getFluxFromGrad(flux, grad, vof, a_dit);

      Real deltaLph[SpaceDim];
      getChangeInSolution(deltaLph, flux, a_dx, vof, a_dit, ebisBox);

      for (int comp = 0; comp < a_lphi.nComp(); comp++)
        {
          a_lphi(vof, comp) += deltaLph[comp];
        }
    }
}

void
NeumannViscousTensorEBBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCFuncEval>();

  m_isFunction = false;
}

void
NeumannViscousTensorEBBCFactory::
setFunction(RefCountedPtr<BaseBCFuncEval> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_isFunction = true;
}

NeumannViscousTensorEBBC*
NeumannViscousTensorEBBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx,
       const IntVect*       a_ghostCellsPhi,
       const IntVect*       a_ghostCellsRhs)
{
  CH_TIME("NeumannViscousTensorEBBC::create");
  NeumannViscousTensorEBBC* fresh;
  fresh = new NeumannViscousTensorEBBC(a_domain,a_layout,a_dx,a_ghostCellsPhi,
                                       a_ghostCellsRhs);

  if (m_isFunction)
    {
      fresh->setFunction(m_func);
    }
  else
    {
      fresh->setValue(m_value);
    }

  return fresh;
}

#include "NamespaceFooter.H"
