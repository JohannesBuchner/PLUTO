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

#include "DirichletViscousTensorEBBC.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

int DirichletViscousTensorEBBC::s_leastSquaresRad = 2;

DirichletViscousTensorEBBC::
DirichletViscousTensorEBBC()
{
}


DirichletViscousTensorEBBC::
DirichletViscousTensorEBBC(
                           const ProblemDomain& a_domain,
                           const EBISLayout&    a_layout,
                           const RealVect&      a_dx,
                           const IntVect*       a_ghostCellsPhi,
                           const IntVect*       a_ghostCellsRhs)
  : m_ghostCellsPhi( *a_ghostCellsPhi ),
    m_ghostCellsRHS( *a_ghostCellsRhs )
{
  CH_assert( bool(a_ghostCellsPhi) == bool(a_ghostCellsRhs) ); // !xor

  m_domain = a_domain;
  m_ebisl = a_layout;

  m_dx = a_dx;

  m_isDefined = false;
}

DirichletViscousTensorEBBC::
~DirichletViscousTensorEBBC()
{
}

void
DirichletViscousTensorEBBC::
define(const LayoutData<IntVectSet>& a_cfivs,
       const Real&                   a_factor)
{
  //factor is 1 here.  ignore.
  CH_assert(m_coefSet);
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

  //make the Dirichlet stencils
  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      m_fluxStencil[ivar].define(dbl);
      m_fluxWeight [ivar].define(dbl);

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {

          m_fluxStencil[ivar][dit()].define(irregIVS[dit()],m_ebisl[dit()].getEBGraph(), 1);
          m_fluxWeight [ivar][dit()].define(irregIVS[dit()],m_ebisl[dit()].getEBGraph(), 1);


          for (vofItIrreg[dit()].reset(); vofItIrreg[dit()].ok(); ++vofItIrreg[dit()])
            {
              Real        fluxWeight[SpaceDim];
              VoFStencil fluxStencil[SpaceDim];
              getFluxStencil(fluxStencil, fluxWeight, dit(),
                             vofItIrreg[dit()](),m_ebisl[dit()],
                             m_dx,a_cfivs[dit()]);

              m_fluxStencil[ivar][dit()](vofItIrreg[dit()](), 0)  = fluxStencil[ivar];
              m_fluxWeight [ivar][dit()](vofItIrreg[dit()](), 0)  =  fluxWeight[ivar];
            }
        }
    }

  m_isDefined = true;
}
//given a normal vector, get spacedim-1 tangential vectors that are
//unit vectors and
void
DirichletViscousTensorEBBC::
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
  if (a_useHomogeneous) return;
  if (!m_isFunction && (Abs(m_value < 1.0e-10))) return;

  CH_assert(m_coefSet);
  const EBISBox&   ebisBox = a_phi.getEBISBox();
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      const VolIndex& vof = a_vofit();
      //want only the difference between the inhomogeneous flux and the homogeneous
      RealVect deltaLph = getInhomogeneousContribution(a_vofit(), a_phi, ebisBox, a_dit, a_dx[0]);
      for (int comp = 0; comp < SpaceDim; comp++)
        {
          a_lphi(vof, comp) += deltaLph[comp];
        }
    }
}

//want only the difference between the inhomogeneous flux and the homogeneous
RealVect
DirichletViscousTensorEBBC::
getInhomogeneousContribution(const VolIndex&  a_vof,
                             const EBCellFAB& a_phi,
                             const EBISBox&   a_ebisBox,
                             const DataIndex& a_dit,
                             const Real&      a_dx)
{
  Real fluxHomog[SpaceDim][SpaceDim];
  Real fluxInhom[SpaceDim][SpaceDim];
  getFlux(fluxHomog, a_vof, a_phi, a_ebisBox, a_dit, a_dx, true);
  getFlux(fluxInhom, a_vof, a_phi, a_ebisBox, a_dit, a_dx, false);
  RealVect normal = a_ebisBox.normal(a_vof);

  //compute f dot n dA/dx
  RealVect retval = RealVect::Zero;

  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          //want only the difference between the inhomogeneous flux and the homogeneous
          Real fluxDiff = fluxInhom[idir][ivar] - fluxHomog[idir][ivar];
          retval[ivar] -= fluxDiff*normal[idir];
        }
    }

  retval  *= (a_ebisBox.bndryArea(a_vof));
  retval  /= a_dx;

  return retval;
}
//want only the difference between the inhomogeneous flux and the homogeneous
void
DirichletViscousTensorEBBC::
getFlux(Real             a_flux[SpaceDim][SpaceDim],
        const VolIndex&  a_vof,
        const EBCellFAB& a_phi,
        const EBISBox&   a_ebisBox,
        const DataIndex& a_dit,
        const Real&      a_dx,
        bool a_homogeneous)
{
  Real gradient[SpaceDim][SpaceDim];
  getGradient(gradient, a_vof, a_phi, a_ebisBox, a_dit, a_dx, a_homogeneous);
  RealVect retval;
  Real lambda = (*m_lambda)[a_dit](a_vof, 0);
  Real    eta =    (*m_eta)[a_dit](a_vof, 0);
  Real divu = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      divu +=  gradient[idir][idir];
    }
  for (int irow = 0; irow < SpaceDim; irow++)
    {
      for (int icol = 0; icol < SpaceDim; icol++)
        {
          Real entry = 0;
          if (irow == icol)
            {
              //first add in lambda*divu I
              entry += lambda*divu;
            }
          entry += eta*(gradient[irow][icol] + gradient[icol][irow]);
          a_flux[irow][icol]  = entry;
        }
    }
}
/****/
void
DirichletViscousTensorEBBC::
getGradient(Real             a_grad[SpaceDim][SpaceDim],
            const VolIndex&  a_vof,
            const EBCellFAB& a_phi,
            const EBISBox&   a_ebisBox,
            const DataIndex& a_dit,
            const Real&      a_dx,
            bool a_homogeneous)
{
  if (m_isFunction && (!a_homogeneous))
    {
      getGradientFunction(a_grad, a_vof, a_phi, a_ebisBox, a_dit, a_dx, a_homogeneous);
    }
  else
    {
      Real weight;
      VoFStencil normalStencil;
      IntVectSet cfivs;
      getNormalStencil(normalStencil, weight,
                       a_vof, a_ebisBox, a_dx*RealVect::Unit, cfivs);
      RealVect normal = a_ebisBox.normal(a_vof);
      RealVect tangents[SpaceDim-1];
      PolyGeom::getTangentVectors(tangents, normal);

      Real Jacobian[SpaceDim][SpaceDim];
      Real Jinverse[SpaceDim][SpaceDim];
      PolyGeom::getJacobianAndInverse(Jacobian, Jinverse, normal, tangents);

      for (int ivar = 0; ivar < SpaceDim; ivar++)
        {
          normalStencil.setAllVariables(ivar);
          Real normalGrad = applyVoFStencil(normalStencil, a_phi, ivar);
          if (!a_homogeneous)
            {
              normalGrad += weight*m_value;
            }
          RealVect normTanGrad = RealVect::Zero;
          normTanGrad[0] = normalGrad;

          //cartesian gradient = jinverse*normalTangent gradient
          RealVect cartesianGrad;
          for (int irow = 0; irow < SpaceDim; irow++)
            {
              cartesianGrad[irow] = 0;
              for (int icol = 0; icol < SpaceDim; icol++)
                {
                  cartesianGrad[irow] += Jinverse[irow][icol]*normTanGrad[icol];
                }
            }
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_grad[ivar][idir] = -cartesianGrad[idir];
            }
        }
    }
}

/****/
void
DirichletViscousTensorEBBC::
getGradInhomOnly(Real             a_grad[SpaceDim][SpaceDim],
                 const Real&      a_weight,
                 const VolIndex&  a_vof,
                 const EBISBox&   a_ebisBox,
                 const Real&      a_dx)
{
  RealVect normal = a_ebisBox.normal(a_vof);
  RealVect tangents[SpaceDim-1];
  PolyGeom::getTangentVectors(tangents, normal);

  Real Jacobian[SpaceDim][SpaceDim];
  Real Jinverse[SpaceDim][SpaceDim];
  PolyGeom::getJacobianAndInverse(Jacobian, Jinverse, normal, tangents);

  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      Real normalGrad;
      RealVect normTanGrad;
      if (!m_isFunction)
        {
          normalGrad = a_weight*m_value;  //JUST DOING INHOM BIT
          normTanGrad = RealVect::Zero;
        }
      else
        {
          MayDay::Error("function case for this not implemented");
        }
      normTanGrad[0] = normalGrad;

      //cartesian gradient = jinverse*normalTangent gradient
      RealVect cartesianGrad;
      for (int irow = 0; irow < SpaceDim; irow++)
        {
          cartesianGrad[irow] = 0;
          for (int icol = 0; icol < SpaceDim; icol++)
            {
              cartesianGrad[irow] += Jinverse[irow][icol]*normTanGrad[icol];
            }
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_grad[ivar][idir] = cartesianGrad[idir];
        }
    }
}
/****/
void
DirichletViscousTensorEBBC::
getGradientStenValue(Real             a_grad[SpaceDim][SpaceDim],
                     const VolIndex&  a_vof,
                     const EBCellFAB& a_phi,
                     const EBISBox&   a_ebisBox,
                     const DataIndex& a_dit,
                     const Real&      a_dx,
                     bool a_homogeneous)
{
  if (m_isFunction)
    {
      MayDay::Error("getGradientStenValue only implemented for now in the value case");
    }
  Real weight;
  VoFStencil normalStencil;
  IntVectSet cfivs;
  getNormalStencil(normalStencil, weight,
                   a_vof, a_ebisBox, a_dx*RealVect::Unit, cfivs);

  VoFStencil gradStencil[SpaceDim][SpaceDim];
  getCartesianGradientStencil(gradStencil, normalStencil,
                              a_dit, a_vof, a_ebisBox, a_dx*RealVect::Unit);
  Real gradInhom[SpaceDim][SpaceDim];
  if (!a_homogeneous)
    {
      getGradInhomOnly(gradInhom, weight, a_vof, a_ebisBox, a_dx);
    }
  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real gradVal = applyVoFStencil(gradStencil[ivar][idir], a_phi, ivar);
          if (!a_homogeneous)
            {
              gradVal += gradInhom[ivar][idir];
            }
          a_grad[ivar][idir] = -gradVal;
        }
    }
}
void
DirichletViscousTensorEBBC::
getGradientFunction(Real             a_grad[SpaceDim][SpaceDim],
                    const VolIndex&  a_vof,
                    const EBCellFAB& a_phi,
                    const EBISBox&   a_ebisBox,
                    const DataIndex& a_dit,
                    const Real&      a_dx,
                    bool a_homogeneous)

{
  if (a_homogeneous)
    {
      MayDay::Error("DirVTEBBC::getgradientfunction not written for homogeneous case");
    }
  //fill a gradient with all cartesian gradients
  RealVect point = RealVect::Unit;
  const IntVect& iv = a_vof.gridIndex();
  point *= 0.5;
  point += iv;
  point += a_ebisBox.bndryCentroid(a_vof);
  point *= a_dx;
  Real gradCart[SpaceDim][SpaceDim], gradNT[SpaceDim][SpaceDim];
  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          gradCart[ivar][idir] = m_func->derivative(point, ivar, idir);
        }
    }

  //transform it to normal-tangent coords
  RealVect normal = a_ebisBox.normal(a_vof);
  RealVect tangents[SpaceDim-1];
  PolyGeom::getTangentVectors(tangents, normal);

  Real Jacobian[SpaceDim][SpaceDim];
  Real Jinverse[SpaceDim][SpaceDim];
  PolyGeom::getJacobianAndInverse(Jacobian, Jinverse, normal, tangents);


  for (int irow = 0; irow < SpaceDim; irow++)
    {
      for (int icol = 0; icol < SpaceDim; icol++)
        {
          Real sum = 0;
          for (int isum = 0; isum < SpaceDim; isum++)
            {
              sum += Jacobian[irow][isum]*gradCart[icol][isum];
            }
          gradNT[irow][icol] = sum;
        }
    }

  //only can use the above stuff for tangential stencils
  Real weight;
  VoFStencil normalStencil;
  IntVectSet cfivs;
  getNormalStencil(normalStencil, weight,
                   a_vof, a_ebisBox, a_dx*RealVect::Unit, cfivs);

  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      normalStencil.setAllVariables(ivar);
      Real normalGrad = applyVoFStencil(normalStencil, a_phi, ivar);
      Real value = m_func->value(point, ivar);
      normalGrad += weight*value;

      RealVect normTanGrad = RealVect::Zero;
      normTanGrad[0] = normalGrad;
      for (int idir = 1; idir < SpaceDim; idir++)
        {
          normTanGrad[idir] = gradNT[idir][ivar];
        }

      //cartesian gradient = jinverse*normalTangent gradient
      RealVect cartesianGrad;
      for (int irow = 0; irow < SpaceDim; irow++)
        {
          cartesianGrad[irow] = 0;
          for (int icol = 0; icol < SpaceDim; icol++)
            {
              cartesianGrad[irow] += Jinverse[irow][icol]*normTanGrad[icol];
            }
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_grad[ivar][idir] = -cartesianGrad[idir];
        }
    }
}
/****/
void
DirichletViscousTensorEBBC::
getJacobianAndInverse(Real a_Jacobian[SpaceDim][SpaceDim],
                      Real a_Jinverse[SpaceDim][SpaceDim],
                      RealVect& a_normal,
                      RealVect a_tangents[SpaceDim-1])
{
  PolyGeom::getJacobianAndInverse(a_Jacobian, a_Jinverse, a_normal, a_tangents);
}

void
DirichletViscousTensorEBBC::
getCartesianGradientStencil(VoFStencil           a_gradStencils[SpaceDim][SpaceDim],
                            VoFStencil &         a_normalStencil,
                            const DataIndex&     a_dit,
                            const VolIndex&      a_vof,
                            const EBISBox&       a_ebisBox,
                            const RealVect&      a_dx)
{
  RealVect normal = a_ebisBox.normal(a_vof);
  RealVect tangents[SpaceDim-1];
  PolyGeom::getTangentVectors(tangents, normal);
  VoFStencil emptyStencil;
  Vector<VoFStencil> normTanGradientStencil(SpaceDim, emptyStencil);
  normTanGradientStencil[0] = a_normalStencil;

  Vector<VoFStencil>   cartesianGradientStencil(SpaceDim);
  Real Jacobian[SpaceDim][SpaceDim];
  Real Jinverse[SpaceDim][SpaceDim];
  getJacobianAndInverse(Jacobian, Jinverse, normal, tangents);
  //for each variable, we have the gradient in the normal-tangent
  //coordinate system.  we need to transform it to cartesian by
  //multiplying by jinverse

  for (int ivar = 0; ivar < SpaceDim; ivar++)
    {
      VoFStencil normalStencil = a_normalStencil;
      normalStencil.setAllVariables(ivar);
      VoFStencil normTanGrad[SpaceDim];
      normTanGrad[0] = normalStencil;

      //cartesian gradient = jinverse*normalTangent gradient
      VoFStencil cartesianGrad[SpaceDim];
      for (int irow = 0; irow < SpaceDim; irow++)
        {
          cartesianGrad[irow].clear();
          for (int icol = 0; icol < SpaceDim; icol++)
            {
              VoFStencil moreSten = normTanGrad[icol];
              moreSten *= Jinverse[irow][icol];
              cartesianGrad[irow] += moreSten;
            }
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_gradStencils[ivar][idir] = cartesianGrad[idir];
        }

    }
}

void
DirichletViscousTensorEBBC::
getFluxStencil(VoFStencil           a_stencils[SpaceDim],
               Real                 a_weights [SpaceDim],
               const DataIndex&     a_dit,
               const VolIndex&      a_vof,
               const EBISBox&       a_ebisBox,
               const RealVect&      a_dx,
               const IntVectSet&    a_cfivs)
{
  VoFStencil gradStencil[SpaceDim][SpaceDim];
  Real weight;
  VoFStencil normalStencil;
  getNormalStencil(normalStencil, weight,
                   a_vof, a_ebisBox, a_dx, a_cfivs);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_weights[idir]=weight;
    }
  getCartesianGradientStencil(gradStencil, normalStencil,
                              a_dit, a_vof, a_ebisBox, a_dx);

  //now get the stencil for the divergence
  VoFStencil divergenceStencil;
  for (int ideriv = 0; ideriv < SpaceDim; ideriv++)
    {
      divergenceStencil += gradStencil[ideriv][ideriv];
    }

  //now get the flux stencils
  //now for all the flux stencils Fij = eta*(gradij +gradji) + lambda*deltaij*divergence
  VoFStencil fluxStencil[SpaceDim][SpaceDim];
  Real lambdaPt = (*m_lambda)[a_dit](a_vof, 0);
  Real    etaPt =    (*m_eta)[a_dit](a_vof, 0);
  for (int icomp = 0; icomp < SpaceDim; icomp++)
    {
      for (int ideriv = 0; ideriv < SpaceDim; ideriv++)
        {
          if (icomp == ideriv)
            {
              fluxStencil[icomp][ideriv] += divergenceStencil;
              fluxStencil[icomp][ideriv] *= lambdaPt;
            }
          VoFStencil gradContrib;
          gradContrib += gradStencil[icomp ][ideriv];
          gradContrib += gradStencil[ideriv][icomp ];
          gradContrib *= etaPt;
          fluxStencil[icomp][ideriv] += gradContrib;
        }
    }
  RealVect normal = a_ebisBox.normal(a_vof);
  //flux through face = sum flux[icomp][ideriv]*normal[deriv]
  for (int icomp = 0; icomp < SpaceDim; icomp++)
    {
      a_stencils[icomp].clear();
      for (int ideriv = 0; ideriv < SpaceDim; ideriv++)
        {
          VoFStencil derivContrib = fluxStencil[icomp][ideriv];
          derivContrib *= normal[ideriv];
          a_stencils[icomp] += derivContrib;
        }
      a_weights[icomp] = weight;
    }
}

void
DirichletViscousTensorEBBC::
getNormalStencil(VoFStencil&          a_stencil,
                 Real&                a_weight,
                 const VolIndex&      a_vof,
                 const EBISBox&       a_ebisBox,
                 const RealVect&      a_dx,
                 const IntVectSet&    a_cfivs)
{
  Vector<VoFStencil>  pointStencils;
  Vector<Real>        distanceAlongLine;
  bool needToDropOrder = getSecondOrderStencil(a_stencil, a_weight,
                                               pointStencils, distanceAlongLine,
                                               a_vof, a_ebisBox, a_dx, a_cfivs);

  //debug force drop
  //needToDropOrder = true;
  //end debug
  if (needToDropOrder)
    {
      getFirstOrderStencil(a_stencil, a_weight, a_vof, a_ebisBox, a_dx, a_cfivs);
    }
}
bool
DirichletViscousTensorEBBC::
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
  //grad = -(x2*x2*(phi1-phi0)-x1*x1(phi2-phi0))/(x2*x2*x1 - x1*x1*x2);
  Real denom = x2*x2*x1 - x1*x1*x2;
  //not done by reference because i want point stencils checkable externally.
  VoFStencil phi1Sten = a_pointStencils[0];
  VoFStencil phi2Sten = a_pointStencils[1];
  phi1Sten *=-x2*x2/denom;
  phi2Sten *= x1*x1/denom;
  //weight is the multiplier of the inhomogeneous value (phi0)
  a_weight = -x1*x1/denom + x2*x2/denom;
  a_stencil += phi1Sten;
  a_stencil += phi2Sten;
  //if we got this far, we have a second order stencil;
  return false;
}

void
DirichletViscousTensorEBBC::
getFirstOrderStencil(VoFStencil&       a_stencil,
                     Real&             a_weight,
                     const VolIndex&   a_vof,
                     const EBISBox&    a_ebisBox,
                     const RealVect&   a_dx,
                     const IntVectSet& a_cfivs)
{
  RealVect normal   = a_ebisBox.normal(  a_vof);
  RealVect centroid = a_ebisBox.centroid(a_vof);
  EBArith::getLeastSquaresGradStenAllVoFsRad(a_stencil, a_weight, normal, centroid,
                                             a_vof, a_ebisBox, a_dx, m_domain, 0, s_leastSquaresRad);
}


DirichletViscousTensorEBBCFactory::DirichletViscousTensorEBBCFactory()
{
  m_value = 12345.6789;
  m_func = RefCountedPtr<BaseBCFuncEval>();

  m_isFunction = false;
}

DirichletViscousTensorEBBCFactory::~DirichletViscousTensorEBBCFactory()
{
}

void
DirichletViscousTensorEBBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCFuncEval>();

  m_isFunction = false;
}

void
DirichletViscousTensorEBBCFactory::
setFunction(RefCountedPtr<BaseBCFuncEval> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_isFunction = true;
}

DirichletViscousTensorEBBC*
DirichletViscousTensorEBBCFactory::
create(
       const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx,
       const IntVect*       a_ghostCellsPhi,
       const IntVect*       a_ghostCellsRhs)
{
  CH_TIME("DirichletViscousTensorEBBC::create");
  DirichletViscousTensorEBBC* fresh;
  fresh = new DirichletViscousTensorEBBC(a_domain,a_layout,a_dx,a_ghostCellsPhi,
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
