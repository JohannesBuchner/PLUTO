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

#include "NeumannViscousTensorDomainBC.H"
#include "TensorCFInterp.H"
#include "ViscousTensorOpF_F.H"
#include "EBViscousTensorOpF_F.H"
#include "EBViscousTensorOp.H"
#include "NamespaceHeader.H"

NeumannViscousTensorDomainBC::
NeumannViscousTensorDomainBC()
{
}

NeumannViscousTensorDomainBC::
~NeumannViscousTensorDomainBC()
{
}

void
NeumannViscousTensorDomainBC::
getFaceFlux(BaseFab<Real>&        a_flux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  // flux is eta(grad B + grad B^T) + lambda* I div B   on the face
  CH_assert(a_phi.nComp() == SpaceDim);


  const Box& box = a_flux.box();
  //remember that box is face centered
  int iside;

  if (a_side == Side::Lo)
    {
      iside = -1;
    }
  else
    {
      iside =  1;
    }
  int ncompGph = SpaceDim*SpaceDim;
  FArrayBox faceGph(box, ncompGph);
  faceGph.setVal(0.);

  if (m_isFunction && (!a_useHomogeneous))
    {
      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          RealVect point = EBArith::getIVLocation(iv,a_dx,a_probLo);
          point[a_idir] -=  0.5 * a_dx[a_idir];//point is now at the face center
          IntVect ivNeigh = iv;
          ivNeigh[a_idir]--;
          const VolIndex vof      = VolIndex(iv,     0);
          const VolIndex vofNeigh = VolIndex(ivNeigh,0);
          const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
          for (int comp = 0; comp < SpaceDim; comp++)
            {
              for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
                {
                  Real deriv = m_func->derivative(point, comp, derivDir);
                  int gradComp = TensorCFInterp::gradIndex(comp, derivDir);
                  faceGph(iv, gradComp) = deriv;
                }
            }
        }
    }
  else
    {
      Real value = 0.0;
      if (!a_useHomogeneous)
        {
          value = m_value;
        }
      for (int comp = 0; comp < SpaceDim; comp++)
        {
          for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
            {
              int gradComp = TensorCFInterp::gradIndex(comp, derivDir);
              if ((derivDir == a_idir) && (comp == a_idir))
                {
                  faceGph.setVal(value, gradComp);
                }
              else
                {
                  faceGph.setVal(0.0,   gradComp);
                }
            }
        }
    }
  getFluxFromGrad(a_flux, faceGph, a_dit, a_idir);
}

void
NeumannViscousTensorDomainBC::
getDivergenceAndGradient(Real&                 a_divergence,
                         Real*                 a_gradient,
                         const int&            a_idir,
                         const FaceIndex&      a_bndryFace,
                         const VolIndex&       a_vof,
                         const EBCellFAB&      a_phi,
                         const RealVect&       a_probLo,
                         const RealVect&       a_dx,
                         const Side::LoHiSide& a_side,
                         const DataIndex&      a_dit,
                         const Real&           a_time,
                         const bool&           a_useHomogeneous)
{
  //get divergence and gradient
  const RealVect  point = EBArith::getFaceLocation(a_bndryFace,a_dx,a_probLo);
  //get gradient stencil and replace normal derivatives
  //with something that takes into account the boundary condition
  for (int icomp = 0; icomp < SpaceDim; icomp++)
    {
      for (int idiff = 0; idiff < SpaceDim; idiff++)
        {
          int gradComp = TensorCFInterp::gradIndex(icomp,idiff);
          a_gradient[gradComp] = 0.0;
          if ((!a_useHomogeneous) && m_isFunction)
            {
              Real deriv = m_func->derivative(point, icomp, idiff);
              a_gradient[gradComp] = deriv;
            }
          else if ((!a_useHomogeneous) && (idiff == a_idir) && (icomp == a_idir))
            {
              a_gradient[gradComp] = m_value;
            }
        }
    }

  //now compute the divergence as a sum of gradients
  a_divergence = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int gradComp = TensorCFInterp::gradIndex(idir, idir);
      a_divergence += a_gradient[gradComp];
    }
}

void
NeumannViscousTensorDomainBC::
getFaceFlux(Real&                 a_faceFlux,
            const VolIndex&       a_vof,
            const int&            a_comp,
            const EBCellFAB&      a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  const EBISBox& ebisBox = a_phi.getEBISBox();

  a_faceFlux = 0.0;
  Vector<FaceIndex> bndryFaces = ebisBox.getFaces(a_vof,a_idir,a_side);
  for (int iface = 0; iface < bndryFaces.size(); iface++)
    {
      Real divergence = 0;
      Real gradient[SpaceDim*SpaceDim];
      const FaceIndex& bndryFace = bndryFaces[iface];
      getDivergenceAndGradient(divergence, gradient, a_idir, bndryFace,
                               a_vof, a_phi, a_probLo, a_dx, a_side, a_dit,
                               a_time, a_useHomogeneous);

      if (a_comp == a_idir)
        {
          Real lambda = (*m_lambda)[a_dit][a_idir](bndryFace, 0);
          a_faceFlux += lambda*divergence;
        }
      Real eta = (*m_eta)[a_dit][a_idir](bndryFace, 0);
      int gradBNormComp = TensorCFInterp::gradIndex(a_comp,a_idir);
      int gradBTranComp = TensorCFInterp::gradIndex(a_comp,a_idir);
      Real gradBNorm = gradient[gradBNormComp];
      Real gradBTran = gradient[gradBTranComp];
      a_faceFlux  += eta*(gradBNorm + gradBTran) ;
      a_faceFlux *= ebisBox.areaFrac(bndryFaces[iface]);
    }
}

NeumannViscousTensorDomainBCFactory::
NeumannViscousTensorDomainBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCFuncEval>();

  m_isFunction = false;
}

NeumannViscousTensorDomainBCFactory::
~NeumannViscousTensorDomainBCFactory()
{
}

void
NeumannViscousTensorDomainBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCFuncEval>();

  m_isFunction = false;
}

void
NeumannViscousTensorDomainBCFactory::
setFunction(RefCountedPtr<BaseBCFuncEval> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
NeumannViscousTensorDomainBC*
NeumannViscousTensorDomainBCFactory::create(const ProblemDomain& a_domain,
                                            const EBISLayout&    a_layout,
                                            const RealVect&      a_dx)
{
  NeumannViscousTensorDomainBC* newBC = new NeumannViscousTensorDomainBC();
  if (m_isFunction)
    {
      newBC->setFunction(m_flux);
    }
  else
    {
      newBC->setValue(m_value);
    }

  return newBC;
}
#include "NamespaceFooter.H"
