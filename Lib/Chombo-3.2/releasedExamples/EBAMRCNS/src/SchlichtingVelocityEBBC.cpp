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
#include "PolyGeom.H"

#include "SchlichtingVelocityEBBC.H"
#include "CH_Timer.H"
#include "EBViscousTensorOp.H"
#include "EBAMRCNSParams.H"

#include "NamespaceHeader.H"

SchlichtingVelocityEBBC::
SchlichtingVelocityEBBC()
{
}


SchlichtingVelocityEBBC::
SchlichtingVelocityEBBC(const SchlichtingParams& a_params,
                        const ProblemDomain& a_domain,
                        const EBISLayout&    a_layout,
                        const RealVect&      a_dx):
  DirichletViscousTensorEBBC(a_domain, a_layout, a_dx, &IntVect::Zero, &IntVect::Zero)
{
  m_params = a_params;
  //this is to set the dirichlet bits correctly for when I use them
  setValue(0.);
}

SchlichtingVelocityEBBC::
~SchlichtingVelocityEBBC()
{
}
/***/
void
SchlichtingVelocityEBBC::
getAverageFluxStencil(VoFStencil& a_aveFluxSten, 
                      const VolIndex&   a_vof, 
                      const DataIndex&  a_dit,
                      int a_idir, int   a_ivar)
{
  a_aveFluxSten.clear();
  int numFaces = 0;
  for(SideIterator sit; sit.ok(); ++sit)
    {
      Vector<FaceIndex> faces = m_eblg.getEBISL()[a_dit].getFaces(a_vof, a_idir, sit());
      for(int iface = 0; iface < faces.size(); iface++)
        {
          numFaces++;
          VoFStencil stenSide;

          const RefCountedPtr<LevelData<EBFluxFAB> >&    etaArg = m_etaOpen;
          const RefCountedPtr<LevelData<EBFluxFAB> >& lambdaArg = m_lambdaOpen;
          const Real                                &     dxArg = m_dx[0];
          const EBLevelGrid                         &   eblgArg = m_eblg;
          const FaceIndex                           &   faceArg = faces[iface];
          

          EBViscousTensorOp::getFluxStencil(stenSide, etaArg, lambdaArg, 
                                            dxArg, eblgArg, faceArg,
                                            a_dit, a_ivar);
          /*
          EBViscousTensorOp::getFluxStencil(stenSide, m_etaOpen, m_lambdaOpen, 
                                            m_dx[0], m_eblg, faces[iface], 
                                            a_dit, a_ivar);
          */

          a_aveFluxSten += stenSide;
        }
    }
  if(numFaces > 0)
    {
      a_aveFluxSten *= (1.0/(Real(numFaces)));
    }
}
/***/
bool
SchlichtingVelocityEBBC::
getSchlichtingStencil(VoFStencil&       a_stencil, 
                      const VolIndex&   a_vof, 
                      const DataIndex&  a_dit,
                      const RealVect&   a_normal, 
                      int a_ivar)
{
  a_stencil.clear();
  Vector<Real> distanceAlongLine;
  Vector<VoFStencil> pointStencils;
  bool dropOrder = false;
  EBArith::johanStencil(dropOrder, pointStencils, distanceAlongLine,
                        a_vof,   m_eblg.getEBISL() [a_dit], 
                        m_dx,  (*m_eblg.getCFIVS())[a_dit]);
  if (dropOrder)
    {
      //just return empty stencil
      return true;
    }

  //if we got this far, sizes should be at least 2
  CH_assert(distanceAlongLine.size() >= 2);
  CH_assert(pointStencils.size() >= 2);
  Real x1 = distanceAlongLine[0];
  Real x2 = distanceAlongLine[1];

  //extrapolated value along line = (x2*phi1 - x1*phi2)/(x2-x1)
  VoFStencil x2phi1 = pointStencils[0];
  x2phi1 *= x2;
  VoFStencil x1phi2 = pointStencils[1];
  x1phi2 *= x1;
  VoFStencil extrapStenc;
  extrapStenc  = x2phi1;
  extrapStenc += x1phi2;
  extrapStenc *= (1.0/(x2-x1));

  //now loop through this mess and for every vof in it calculate flux dot normal and 
  //extrapolate that badboy
  for(int ivof = 0; ivof < extrapStenc.size(); ivof++)
    {
      const VolIndex& vofsten = extrapStenc.vof(ivof);
      const Real&     vofweig = extrapStenc.weight(ivof);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          VoFStencil aveFluxSten;
          getAverageFluxStencil(aveFluxSten, vofsten, a_dit, idir, a_ivar);
          aveFluxSten *= a_normal[idir];
          aveFluxSten *= vofweig;
          a_stencil += aveFluxSten;
        }
    }
  return false;
}
/*****/
void
SchlichtingVelocityEBBC::
define(const LayoutData<IntVectSet>& a_cfivs,
       const Real&                   a_factor)
{
  DirichletViscousTensorEBBC::define(a_cfivs, a_factor);
  const DisjointBoxLayout& dbl = m_ebisl.getDisjointLayout();

  LayoutData<IntVectSet>   irregIVS;
  LayoutData<VoFIterator > vofItIrreg;
  vofItIrreg.define(dbl); // vofiterator cache
  irregIVS.define(dbl);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      irregIVS[dit()] = m_ebisl[dit()].getIrregIVS(dbl.get(dit()));
      vofItIrreg[dit()].define(irregIVS[dit()],m_ebisl[dit()].getEBGraph());
    }

  //the Dirichlet stencils have been made.  set them to our special sauce stencil if not sticky
  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      int numSticky = 0;
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {

          for (vofItIrreg[dit()].reset(); vofItIrreg[dit()].ok(); ++vofItIrreg[dit()])
            {
              const VolIndex& vof = vofItIrreg[dit()]();
              RealVect normal = m_ebisl[dit()].normal(vof);
              RealVect bdCent = m_ebisl[dit()].bndryCentroid(vof);
              bdCent *= m_dx;
              RealVect point  = EBArith::getVofLocation(vof, m_dx, RealVect::Zero);
              point += bdCent;

              bool sticky = m_params.isThisPointSticky(point, normal);
              if(!sticky)
                {
                  m_fluxWeight[ivar][dit()](vof, 0) = 0;
                  bool oneDrop = getSchlichtingStencil(m_fluxStencil[ivar][dit()](vof, 0), 
                                                       vof,  dit(), normal, ivar);
                  if(oneDrop)
                    {
                      pout() << "SchlichtingVelocityEBBC: dropped order at vof " << vof.gridIndex() << endl;
                    }
                }
              else
                {
                  numSticky++;
                }
            }
        }
      pout() << "number of sticky points = " << numSticky << endl;
    }

  m_isDefined = true;
}
//given a normal vector, get spacedim-1 tangential vectors that are
//unit vectors and
void
SchlichtingVelocityEBBC::
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
  //always homogeneous
  return;
}

//homogeneous
RealVect
SchlichtingVelocityEBBC::
getInhomogeneousContribution(const VolIndex&  a_vof,
                             const EBCellFAB& a_phi,
                             const EBISBox&   a_ebisBox,
                             const DataIndex& a_dit,
                             const Real&      a_dx)
{
  return RealVect::Zero;
}

////////////////// factory stuff
SchlichtingVelocityEBBCFactory::SchlichtingVelocityEBBCFactory(const SchlichtingParams& a_params)
{
  m_params = a_params;
}

SchlichtingVelocityEBBCFactory::~SchlichtingVelocityEBBCFactory()
{
}


SchlichtingVelocityEBBC*
SchlichtingVelocityEBBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx,
       const IntVect*       a_ghostCellsPhi,
       const IntVect*       a_ghostCellsRhs)
{
  CH_TIME("SchlichtingVelocityEBBC::create");
  SchlichtingVelocityEBBC* fresh = new SchlichtingVelocityEBBC(m_params, a_domain,a_layout,a_dx);
  return fresh;
}

#include "NamespaceFooter.H"
