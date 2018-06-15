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
#include "EBArith.H"

#include "DirichletViscousTensorDomainBC.H"
#include "TensorCFInterp.H"
#include "ViscousTensorOpF_F.H"
#include "EBViscousTensorOp.H"
#include "EBViscousTensorOpF_F.H"
#include "NamespaceHeader.H"


DirichletViscousTensorDomainBC::
DirichletViscousTensorDomainBC()
{
}

DirichletViscousTensorDomainBC::
~DirichletViscousTensorDomainBC()
{
}

void
DirichletViscousTensorDomainBC::
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

  const FArrayBox& phiFAB = (const FArrayBox& )a_phi;
  const Box& faceBox = a_flux.box();

  int iside;

  if (a_side == Side::Lo)
    {
      iside = -1;
    }
  else
    {
      iside =  1;
    }
  int ncompPhi = SpaceDim;
  int ncompGph = SpaceDim*SpaceDim;
  FArrayBox faceDiv(faceBox, 1);
  FArrayBox faceGph(faceBox, ncompGph);
  FArrayBox phiValu(faceBox, ncompPhi);
  faceDiv.setVal(0.);
  faceGph.setVal(0.);

  if (m_isFunction && (!a_useHomogeneous))
    {
      for (BoxIterator bit(faceBox); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          RealVect point = EBArith::getIVLocation(iv,a_dx,a_probLo);
          point[a_idir] -=  0.5 * a_dx[a_idir];//point is now at the face center
          for (int comp = 0; comp < SpaceDim; comp++)
            {
              Real value = m_func->value(point, comp);
              phiValu(iv, comp) = value;
              for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
                {
                  if (derivDir != a_idir)
                    {
                      Real deriv   = m_func->derivative(point, comp, derivDir);
                      int gradComp = TensorCFInterp::gradIndex(comp, derivDir);
                      faceGph(iv, gradComp) = deriv;
                    }
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
      phiValu.setVal(value);
      faceGph.setVal(0.);
    }

  //want to use the value for normal derivs
  int divDir = a_idir;
  for (int velDir = 0; velDir < SpaceDim; velDir++)
    {
      int gradcomp = TensorCFInterp::gradIndex(velDir,divDir);
      Real dx = a_dx[0];
      FORT_NORMALGRADVISCDIRCH(CHF_FRA1(faceGph, gradcomp),
                               CHF_FRA1(phiFAB,   velDir),
                               CHF_FRA1(phiValu, velDir),
                               CHF_BOX(faceBox),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_INT(iside),
                               CHF_CONST_INT(a_idir));
    }

  getFluxFromGrad(a_flux, faceGph, a_dit, a_idir);
}
/****/
void
DirichletViscousTensorDomainBC::
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

void
DirichletViscousTensorDomainBC::
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
  const EBISBox& ebisBox = a_phi.getEBISBox();
  const RealVect  point = EBArith::getFaceLocation(a_bndryFace,a_dx,a_probLo);
  const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
  //get gradient stencil and replace normal derivatives
  //with something that takes into account the boundary condition
  for (int icomp = 0; icomp < SpaceDim; icomp++)
    {
      for (int idiff = 0; idiff < SpaceDim; idiff++)
        {
          int gradComp = TensorCFInterp::gradIndex(icomp,idiff);
          VoFStencil vofSten;
          EBViscousTensorOp::getGradientStencil(vofSten, icomp , idiff, a_bndryFace, a_dit, a_dx[0], m_eblg);
          a_gradient[gradComp] = applyVoFStencil(vofSten, a_phi, icomp);
        }
    }

  //replace normal grad here with normal der taking into account the boundary condition
  for (int comp = 0; comp < SpaceDim; comp++)
    {
      Real value = 0;
      if (!a_useHomogeneous)
        {
          if (m_isFunction)
            {
              value = m_func->value(point, comp);
            }
          else
            {
              if (comp == a_idir)
                {
                  value = m_value;
                }
              else
                {
                  value = 0.0;
                }
            }
        }
      Side::LoHiSide flipsd = flip(a_side);
      VolIndex vof1 = a_vof;
      Real phi1= a_phi(vof1, comp);
      Real phi2 = 0; //just to get the compiler to shut up
      bool hasVoF2;

      Vector<FaceIndex> faces = ebisBox.getFaces(a_vof, a_idir, flipsd);
      if ((faces.size() ==1) && (!faces[0].isBoundary()))
        {
          hasVoF2 = true;
          VolIndex vof2 = faces[0].getVoF(flipsd);
          phi2 = a_phi(vof2, comp);
        }
      else
        {
          hasVoF2 = false;
        }

      Real normalDer = 0;
      if (hasVoF2)
        {
          normalDer = (1.0/(6.*a_dx[0]))*(9.*(phi1-value) - (phi2-value));
        }
      else
        {
          normalDer = (1.0/(0.5*a_dx[0]))*(phi1-value);
        }

      int gradComp = TensorCFInterp::gradIndex(comp,a_idir);
      a_gradient[gradComp] = normalDer;
    }

  //now compute the divergence as a sum of gradients
  a_divergence = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int gradComp = TensorCFInterp::gradIndex(idir, idir);
      a_divergence += a_gradient[gradComp];
    }
}

DirichletViscousTensorDomainBCFactory::
DirichletViscousTensorDomainBCFactory()
  : m_isFunctional(false),
    m_value(12345.6789),
    m_func(RefCountedPtr<BaseBCFuncEval>())
{
}

DirichletViscousTensorDomainBCFactory::
~DirichletViscousTensorDomainBCFactory()
{
}

void
DirichletViscousTensorDomainBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCFuncEval>();

  m_isFunctional = false;
}

void
DirichletViscousTensorDomainBCFactory::
setFunction(RefCountedPtr<BaseBCFuncEval> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_isFunctional = true;
}

DirichletViscousTensorDomainBC*
DirichletViscousTensorDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  DirichletViscousTensorDomainBC* newBC = new DirichletViscousTensorDomainBC();
  if (m_isFunctional)
    {
      newBC->setFunction(m_func);
    }
  else
    {
      newBC->setValue(m_value);
    }
  return newBC;
}
#include "NamespaceFooter.H"
