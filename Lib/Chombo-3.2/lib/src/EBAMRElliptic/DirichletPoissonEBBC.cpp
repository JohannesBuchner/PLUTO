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

#include "DirichletPoissonEBBC.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

Real g_simulationTime = 0.0;

int DirichletPoissonEBBC::s_velComp = 0;
int DirichletPoissonEBBC::s_leastSquaresRad = 2;
bool DirichletPoissonEBBC::s_areaFracWeighted = false;
bool DirichletPoissonEBBC::s_useQuadrantBasedStencil = true;

DirichletPoissonEBBC::DirichletPoissonEBBC()
{
  m_dataBased = false;
}

DirichletPoissonEBBC::DirichletPoissonEBBC(const ProblemDomain& a_domain,
                                           const EBISLayout&    a_layout,
                                           const RealVect&      a_dx,
                                           const IntVect*       a_ghostCellsPhi /*=0*/,
                                           const IntVect*       a_ghostCellsRhs /*=0*/)
{
  construct(a_domain, a_layout, a_dx, a_ghostCellsPhi, a_ghostCellsRhs);
}

void
DirichletPoissonEBBC::construct(const ProblemDomain& a_domain,
                                const EBISLayout&    a_layout,
                                const RealVect&      a_dx,
                                const IntVect*       a_ghostCellsPhi /*=0*/,
                                const IntVect*       a_ghostCellsRhs /*=0*/)
{
  m_ghostCellsPhi = (*a_ghostCellsPhi) ;
  m_ghostCellsRHS = (*a_ghostCellsRhs) ;

  CH_assert( bool(a_ghostCellsPhi) == bool(a_ghostCellsRhs) ); // !xor

  m_value = 12345.6789;
  m_func = RefCountedPtr<BaseBCValue>();

  m_order = 1;

  m_onlyHomogeneous = true;
  m_isFunction = false;

  m_domain = a_domain;
  m_layout = a_layout;

  m_dx = a_dx;

  m_isDefined = false;
}

DirichletPoissonEBBC::~DirichletPoissonEBBC()
{
}

void DirichletPoissonEBBC::define(const LayoutData<IntVectSet>& a_cfivs,
                                  const Real&                   a_factor)
{
  const DisjointBoxLayout& dbl = m_layout.getDisjointLayout();

  LayoutData<VoFIterator > vofItIrreg;
  vofItIrreg.define(dbl); // vofiterator cache

  m_fluxStencil.define(dbl);
  m_fluxWeight.define(dbl);

  //make the Dirichlet stencils
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& curBox = dbl[dit()];
      const EBISBox& curEBISBox = m_layout[dit()];
      const EBGraph& curEBGraph = curEBISBox.getEBGraph();
      const IntVectSet& cfivsThisBox = a_cfivs[dit()];

      IntVectSet notRegular;
      int nComps = 1;

      notRegular |= curEBISBox.getIrregIVS  (curBox);
      notRegular |= curEBISBox.getMultiCells(curBox);

      vofItIrreg[dit()].define(notRegular,curEBISBox.getEBGraph());

      BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = m_fluxStencil[dit()];
      BaseIVFAB<Real>&       curWeightBaseIVFAB  = m_fluxWeight[dit()];

      curStencilBaseIVFAB.define(notRegular,curEBGraph,nComps);
      curWeightBaseIVFAB.define(notRegular,curEBGraph,nComps);

      for (VoFIterator vofit(notRegular,curEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();

          VoFStencil& curStencil = curStencilBaseIVFAB(vof,0);
          Real areaFrac = curEBISBox.bndryArea(vof);
          Real&       curWeight  = curWeightBaseIVFAB(vof,0);

          if (m_order == 1)
            {
              getFirstOrderStencil(curStencil,curWeight,vof,curEBISBox,m_dx);
            }
          else if (m_order == 2)
            {
              getSecondOrderStencil(curStencil,curWeight,vof,curEBISBox,m_dx,cfivsThisBox);
            }
          else
            {
              MayDay::Error("DirichletPoissonEBBC::define stencil order not 1 or 2");
            }
          //Need to magnify weight in fluxStencil with areafrac*factor, factor = 1/dx;
          //Pass the factor in because m_dx[0] in here may not be same as in EBAMRPoissonOp
          curStencil *= areaFrac*a_factor;
          if (s_areaFracWeighted)
            {
              curStencil *= curEBISBox.areaFracScaling(vof);
            }
        }
    }

  m_isDefined = true;
}

void DirichletPoissonEBBC::setOrder(int a_order)
{
  CH_assert(a_order >= 1 && a_order <= 2);

  if (m_order != a_order)
    {
      m_isDefined = false;
    }

  m_order = a_order;
}

void DirichletPoissonEBBC::setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}

void DirichletPoissonEBBC::setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}

void DirichletPoissonEBBC::applyEBFluxPoint(const VolIndex&               a_vof,
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
  const EBISBox&   ebisBox = a_phi.getEBISBox();
  int comp = 0;

  Real value;
  const VolIndex& vof = a_vof;
  if (m_dataBased)
    {
      value = (*m_data)[a_dit](vof, 0);
    }
  else if (m_isFunction)
    {
      // Compute the bndryCentroid location in physical coordinates
      RealVect point = RealVect::Unit;
      const IntVect& iv = vof.gridIndex();
      point *= 0.5;
      point += iv;
      point += ebisBox.bndryCentroid(vof);
      point *= m_dx;
      point += a_probLo;

      RealVect normal = ebisBox.normal(vof);
      value = m_func->value(point, normal, a_time, s_velComp);
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("DirichletPoissonEBBC::applyEBFlux called with undefined inhomogeneous BC");
        }

      value = m_value;
    }

  const BaseIVFAB<Real>& curWeightBaseIVFAB = m_fluxWeight[a_dit];
  const Real& curWeight = curWeightBaseIVFAB(vof,comp);

  Real flux = curWeight * value;

  Real areaFrac = ebisBox.bndryArea(vof);
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

void DirichletPoissonEBBC::applyEBFlux(EBCellFAB&                    a_lphi,
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
  CH_TIME("DirichletPoissonEBBC::applyEBFlux");
  CH_assert(a_lphi.nComp() == 1 );
  CH_assert(a_phi.nComp() == 1);
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      const VolIndex& vof = a_vofit();
      applyEBFluxPoint(vof, 
                       a_lphi,            
                       a_phi,             
                       a_vofit,           
                       a_cfivs,           
                       a_dit,             
                       a_probLo,          
                       a_dx,              
                       a_factor,          
                       a_useHomogeneous,  
                       a_time)            ;


    }
}
bool
DirichletPoissonEBBC::
getSecondOrderStencil(VoFStencil&          a_stencil,
                      Real&                a_weight,
                      Vector<VoFStencil>&  a_pointStencils,
                      Vector<Real>&        a_distanceAlongLine,
                      const VolIndex&      a_vof,
                      const EBISBox&       a_ebisBox,
                      const RealVect&      a_dx,
                      const IntVectSet&    a_cfivs)
{

  a_stencil.clear();
  bool dropOrder = false;
  EBArith::johanStencil(dropOrder, a_pointStencils, a_distanceAlongLine,
                        a_vof, a_ebisBox, a_dx, a_cfivs);
  if (dropOrder)
    {
      return true;
    }

  //if we got this far, sizes should be at least 2
  CH_assert(a_distanceAlongLine.size() >= 2);
  CH_assert(a_pointStencils.size() >= 2);
  Real x1 = a_distanceAlongLine[0];
  Real x2 = a_distanceAlongLine[1];
  //fit quadratic function to line and find the gradient at the origin
  //grad = (x2*x2*(phi1-phi0)-x1*x1(phi2-phi0))/(x2*x2*x1 - x1*x1*x2);
  Real denom = x2*x2*x1 - x1*x1*x2;
  //not done by reference because i want point stencils checkable externally.
  VoFStencil phi1Sten = a_pointStencils[0];
  VoFStencil phi2Sten = a_pointStencils[1];
  phi1Sten *=-x2*x2/denom;
  phi2Sten *= x1*x1/denom;
  //weight is the multiplier of the inhomogeneous value (phi0)
  a_weight =-x1*x1/denom + x2*x2/denom;
  a_stencil += phi1Sten;
  a_stencil += phi2Sten;
  //if we got this far, we have a second order stencil;
  return false;
}

void DirichletPoissonEBBC::getFirstOrderStencil(VoFStencil&     a_stencil,
                                                Real&           a_weight,
                                                const VolIndex& a_vof,
                                                const EBISBox&  a_ebisBox,
                                                const RealVect& a_dx)
{
  CH_TIME("DirichletPoissonEBBC::getFirstOrderStencil1");
  RealVect normal   = a_ebisBox.normal(  a_vof);
  RealVect centroid = a_ebisBox.centroid(a_vof);
  if (s_useQuadrantBasedStencil)
    {
      EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_vof, a_ebisBox, a_dx, m_domain, 0);
    }
  else
    {
      EBArith::getLeastSquaresGradStenAllVoFsRad(a_stencil, a_weight, normal, centroid,  a_vof, a_ebisBox, a_dx, m_domain, 0, s_leastSquaresRad);
    }
}

void DirichletPoissonEBBC::getSecondOrderStencil(VoFStencil&       a_stencil,
                                                 Real&             a_weight,
                                                 const VolIndex&   a_vof,
                                                 const EBISBox&    a_ebisBox,
                                                 const RealVect&   a_dx,
                                                 const IntVectSet& a_cfivs)
{
  Vector<VoFStencil>  pointStencils;
  Vector<Real>        distanceAlongLine;
  bool needToDropOrder = getSecondOrderStencil(a_stencil, a_weight,
                                               pointStencils, distanceAlongLine,
                                               a_vof, a_ebisBox, a_dx, a_cfivs);
  if (needToDropOrder)
    {
      getFirstOrderStencil(a_stencil,a_weight,a_vof,a_ebisBox,a_dx);
    }
}


DirichletPoissonEBBCFactory::DirichletPoissonEBBCFactory()
{
  m_value = 12345.6789;
  m_func = RefCountedPtr<BaseBCValue>();

  m_order = 1;

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

DirichletPoissonEBBCFactory::~DirichletPoissonEBBCFactory()
{
}

void DirichletPoissonEBBCFactory::setOrder(int a_order)
{
  CH_assert(a_order >= 1 && a_order <= 2);
  m_order = a_order;
}

void DirichletPoissonEBBCFactory::setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}

void DirichletPoissonEBBCFactory::setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}

DirichletPoissonEBBC* DirichletPoissonEBBCFactory::create(const ProblemDomain& a_domain,
                                                          const EBISLayout&    a_layout,
                                                          const RealVect&      a_dx,
                                                          const IntVect*       a_ghostCellsPhi /*=0*/,
                                                          const IntVect*       a_ghostCellsRhs /*=0*/)
{
  CH_TIME("DirichletPoissonEBBC::create");
  DirichletPoissonEBBC* fresh = new DirichletPoissonEBBC(a_domain,a_layout,a_dx,a_ghostCellsPhi,
                                                         a_ghostCellsRhs);

  fresh->setOrder(m_order);

  if (!m_onlyHomogeneous)
    {
      if (!m_isFunction)
        {
          fresh->setValue(m_value);
        }
      else
        {
          fresh->setFunction(m_func);
        }
    }

  return fresh;
}

//deprecated stuff to get EBAMRPoissonOp working again
void DirichletPoissonEBBC::define(const LayoutData<IntVectSet>& a_cfivs)
{
  CH_TIME("DirichletPoissonEBBC::define");
  const DisjointBoxLayout& dbl = m_layout.getDisjointLayout();

  m_fluxStencil.define(dbl);
  m_fluxWeight.define(dbl);

  //make the Dirichlet stencils
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& curBox = dbl[dit()];
      const EBISBox& curEBISBox = m_layout[dit()];
      const EBGraph& curEBGraph = curEBISBox.getEBGraph();
      const IntVectSet& cfivsThisBox = a_cfivs[dit()];

      IntVectSet notRegular;
      int nComps = 1;

      notRegular |= curEBISBox.getIrregIVS  (curBox);
      notRegular |= curEBISBox.getMultiCells(curBox);

      BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = m_fluxStencil[dit()];
      BaseIVFAB<Real>&       curWeightBaseIVFAB  = m_fluxWeight[dit()];

      curStencilBaseIVFAB.define(notRegular,curEBGraph,nComps);
      curWeightBaseIVFAB.define(notRegular,curEBGraph,nComps);

      for (VoFIterator vofit(notRegular,curEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();

          VoFStencil& curStencil = curStencilBaseIVFAB(vof,0);
          Real&       curWeight  = curWeightBaseIVFAB(vof,0);

          if (m_order == 1)
            {
              getFirstOrderStencil(curStencil,curWeight,vof,curEBISBox,m_dx);
            }
          else if (m_order == 2)
            {
              getSecondOrderStencil(curStencil,curWeight,vof,curEBISBox,m_dx,cfivsThisBox);
            }
          else
            {
              MayDay::Error("DirichletPoissonEBBC::define stencil order not 1 or 2");
            }
        }
    }
  m_isDefined = true;
}
//deprecated stuff to get EBAMRPoissonOp working again
void DirichletPoissonEBBC::getEBFlux(Real&                         a_flux,
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
  CH_TIME("DirichletPoissonEBBC::getEBFlux");
  //  CH_assert(a_phi.nComp() == 1);

  const EBCellFAB& curPhi = a_phi[a_dit];
  const EBISBox&   curEBISBox = curPhi.getEBISBox();

  int comp = 0;

  if (!m_isDefined)
    {
      define(a_cfivs);
    }

  // Compute the bndryCentroid location in physical coordinates
  RealVect point = RealVect::Unit;
  const IntVect& iv = a_vof.gridIndex();
  point *= 0.5;
  point += iv;
  point += curEBISBox.bndryCentroid(a_vof);
  point *= a_dx;
  point += a_probLo;

  Real value;

  if (a_useHomogeneous)
    {
      value = 0.0;
    }
  else
    {
      if (m_dataBased)
        {
          value = (*m_data)[a_dit](a_vof, 0);
        }
      else if (m_isFunction)
        {
          RealVect normal = curEBISBox.normal(a_vof);
          value = m_func->value(point, normal, a_time,s_velComp);
        }
      else
        {
          if (m_onlyHomogeneous)
            {
              MayDay::Error("DirichletPoissonEBBC::getFaceFlux called with undefined inhomogeneous BC");
            }

          value = m_value;
        }
    }

  a_flux = 0.0;
  const BaseIVFAB<VoFStencil>& curStencilBaseIVFAB = m_fluxStencil[a_dit];
  const VoFStencil& curStencil = curStencilBaseIVFAB(a_vof,comp);
  for (int i = 0; i < curStencil.size(); i++)
    {
      const VolIndex& curVoF = curStencil.vof(i);
      Real weight = curStencil.weight(i);
      Real phiVal = curPhi(curVoF,comp);
      a_flux += weight * phiVal;
    }

  const BaseIVFAB<Real>& curWeightBaseIVFAB = m_fluxWeight[a_dit];
  const Real& curWeight = curWeightBaseIVFAB(a_vof,comp);

  //areafrac already lives in the stencil
  a_flux += curWeight * value;
}

#include "NamespaceFooter.H"
