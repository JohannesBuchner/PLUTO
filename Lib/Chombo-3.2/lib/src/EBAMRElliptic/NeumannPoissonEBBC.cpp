#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"

#include "NeumannPoissonEBBC.H"
#include "EBStencil.H"
#include "NamespaceHeader.H"
void NeumannPoissonEBBC::getEBFlux(Real&                         a_flux,
                                   const VolIndex&               a_vof,
                                   const LevelData<EBCellFAB>&   a_phi,
                                   const LayoutData<IntVectSet>& a_cfivs,
                                   const DataIndex&              a_dit,
                                   const RealVect&               a_probLo,
                                   const RealVect&               a_dx,
                                   const bool&                   a_useHomogeneous,
                                   const Real&                   a_time,
                                   const pair<int,Real>*         a_cacheHint )
{
  CH_assert(a_phi.nComp() == 1);

  const EBCellFAB& curPhi     = a_phi[a_dit];
  const EBISBox&   curEBISBox = curPhi.getEBISBox();

  if (a_useHomogeneous)
    {
      a_flux = 0.0;
    }
  else
    {
      if (m_dataBased)
        {
          a_flux = (*m_data)[a_dit](a_vof, 0);
        }
      else if (m_isFunction)
        {
          const RealVect& centroid = curEBISBox.bndryCentroid(a_vof);
          const RealVect&   normal = curEBISBox.normal(a_vof);

          Real value = m_flux->value(a_vof,centroid,normal,a_dx,a_probLo,a_dit,a_time,0);
          a_flux = -value;
        }
      else
        {
          if (m_onlyHomogeneous)
            {
              MayDay::Error("NeumannPoissonEBBC::getFaceFlux called with undefined inhomogeneous BC");
            }

          a_flux = m_value;
        }
    }

  const Real& areaFrac = curEBISBox.bndryArea(a_vof);
  a_flux *= areaFrac;
}

NeumannPoissonEBBC::NeumannPoissonEBBC(const ProblemDomain& a_domain,
                                       const EBISLayout&    a_layout,
                                       const RealVect&      a_dx)
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();
  m_dataBased = false;
  m_onlyHomogeneous = true;
  m_isFunction = false;
}

NeumannPoissonEBBC::~NeumannPoissonEBBC()
{
}

void NeumannPoissonEBBC::setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}

void NeumannPoissonEBBC::setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}

void NeumannPoissonEBBC::applyEBFluxPoint(const VolIndex&               a_vof,
                                          EBCellFAB&                    a_lphi,
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
  Real flux = 0.0;

  const EBISBox&   ebisBox = a_phi.getEBISBox();
  const VolIndex& vof = a_vof;

  if (m_dataBased)
    {
      flux = (*m_data)[a_dit](vof, 0);
    }
  else if (m_isFunction)
    {
      const RealVect& centroid = ebisBox.bndryCentroid(vof);
      const RealVect&   normal = ebisBox.normal(vof);

      Real value = m_flux->value(vof,centroid,normal,a_dx,a_probLo,a_dit,a_time,0);
      flux = -value;
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("NeumannPoissonEBBC::getFaceFlux called with undefined inhomogeneous BC");
        }

      flux = m_value;
    }

  const Real& areaFrac = ebisBox.bndryArea(vof);
  flux *= areaFrac;

  Real* lphiPtr = NULL;
  int offset;
  //get ghosted box of lphi
  Box boxLph = a_lphi.getSingleValuedFAB().box();
  //check to see if multi-valued cell or not
  if (ebisBox.numVoFs(vof.gridIndex()) > 1)
    {
      const IntVectSet& ivsLph = ebisBox.getMultiCells(boxLph);
      const EBGraph& ebgraph = ebisBox.getEBGraph();
      BaseIVFAB<Real> baseivfabLph(ivsLph, ebgraph, 1);

      offset = baseivfabLph.getIndex(vof, 0) - baseivfabLph.dataPtr(0);
      Real* multiValuedPtrLph = a_lphi.getMultiValuedFAB().dataPtr(0);
      lphiPtr  = multiValuedPtrLph + offset;
    }
  else
    {
      const IntVect& smallendLph = boxLph.smallEnd();
      IntVect ivLph = vof.gridIndex()  - smallendLph;
      IntVect ncellsLph = boxLph.size();

      offset = ivLph[0] + ivLph[1]*ncellsLph[0] ;
#if CH_SPACEDIM==3
      offset +=  ivLph[2]*ncellsLph[0]*ncellsLph[1];
#endif
      Real* singleValuedPtrLph = a_lphi.getSingleValuedFAB().dataPtr();
      lphiPtr  = singleValuedPtrLph + offset;
    }
  Real& lphi = *lphiPtr;
  lphi += flux * a_factor;
}
void NeumannPoissonEBBC::applyEBFlux(EBCellFAB&                    a_lphi,
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
  CH_TIME("NeumannPoissonEBBC::applyEBFlux");
  CH_assert(a_lphi.nComp() == 1 );
  CH_assert(a_phi.nComp() == 1);

  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      applyEBFluxPoint(a_vofit(),
                       a_lphi,
                       a_phi,
                       a_vofit,
                       a_cfivs,
                       a_dit,
                       a_probLo,
                       a_dx,
                       a_factor,
                       a_useHomogeneous,
                       a_time);
    }
}

NeumannPoissonEBBCFactory::NeumannPoissonEBBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

NeumannPoissonEBBCFactory::~NeumannPoissonEBBCFactory()
{
}

void NeumannPoissonEBBCFactory::setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}

void NeumannPoissonEBBCFactory::setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}

NeumannPoissonEBBC* NeumannPoissonEBBCFactory::create(
                                                      const ProblemDomain& a_domain,
                                                      const EBISLayout&    a_layout,
                                                      const RealVect&      a_dx,
                                                      const IntVect*       a_ghostCellsPhi /*=0*/,
                                                      const IntVect*       a_ghostCellsRhs /*=0*/)
{
  NeumannPoissonEBBC* fresh = new NeumannPoissonEBBC(a_domain,a_layout,a_dx);

  if (!m_onlyHomogeneous)
    {
      if (!m_isFunction)
        {
          fresh->setValue(m_value);
        }
      else
        {
          fresh->setFunction(m_flux);
        }
    }

  return fresh;
}
#include "NamespaceFooter.H"
