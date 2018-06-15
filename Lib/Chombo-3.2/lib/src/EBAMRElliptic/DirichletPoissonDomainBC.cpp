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
#include "PolyGeom.H"

#include "DirichletPoissonDomainBC.H"
#include "DirichletPoissonDomainBCF_F.H"
#include "NamespaceHeader.H"


DirichletPoissonDomainBC::DirichletPoissonDomainBC()
 : m_onlyHomogeneous(true),
   m_isFunctional(false),
   m_value(12345.6789),
   m_func(RefCountedPtr<BaseBCValue>()),
   m_ebOrder(2)
{
}

DirichletPoissonDomainBC::~DirichletPoissonDomainBC()
{
}

void DirichletPoissonDomainBC::setValue(Real a_value)
{
   m_onlyHomogeneous = false;
   m_isFunctional = false;
   m_value = a_value;
   m_func = RefCountedPtr<BaseBCValue>();
}

void DirichletPoissonDomainBC::setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunctional = true;
}

void DirichletPoissonDomainBC::setEBOrder(int a_ebOrder)
{
  CH_assert(a_ebOrder >= 1);
  CH_assert(a_ebOrder <= 2);

  m_ebOrder = a_ebOrder;
}
///
/**
   This never gets called.  InflowOutflowPoissonDomainBC::getFaceVel takes care of it.
*/
void
DirichletPoissonDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{
  CH_assert(a_idir == a_face.direction());
  Real value;
  if (m_isFunctional)
    {
      RealVect pt;
      IntVect iv = a_face.gridIndex(Side::Hi);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (idir != a_face.direction())
            {
              Real ptval = a_dx[a_idir]*(Real(iv[idir]) + 0.5);
              pt[idir] = ptval;
            }
          else
            {
              pt[idir] = a_dx[a_idir]*(Real(iv[idir]));
            }
        }
      RealVect normal = EBArith::getDomainNormal(a_idir, a_side);

      value = m_func->value(pt, normal, a_time,a_icomp);

    }
  else
    {
      value = m_value;
    }
  a_faceFlux = value;
}
///
/**
   This is called by InflowOutflowHelmholtzDomainBC::getFaceFlux,
   which is called by EBAMRPoissonOp::applyDomainFlux for applyOp in reg cells.
   For higher-order Dirichlet boundary conditions on velocity and viscous operator.
*/
void DirichletPoissonDomainBC::getHigherOrderFaceFlux(BaseFab<Real>&        a_faceFlux,
                                                      const BaseFab<Real>&  a_phi,
                                                      const RealVect&       a_probLo,
                                                      const RealVect&       a_dx,
                                                      const int&            a_idir,
                                                      const Side::LoHiSide& a_side,
                                                      const DataIndex&      a_dit,
                                                      const Real&           a_time,
                                                      const bool&           a_useHomogeneous)
{
  for (int comp=0; comp<a_phi.nComp(); comp++)
  {
    const Box& box = a_faceFlux.box();

    int iside;

    if (a_side == Side::Lo)
      {
        iside = 1;
      }
    else
      {
        iside = -1;
      }

    if (a_useHomogeneous)
      {
        Real value = 0.0;

        FORT_SETHODIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
                                    CHF_CONST_FRA1(a_phi,comp),
                                    CHF_CONST_REAL(value),
                                    CHF_CONST_REALVECT(a_dx),
                                    CHF_CONST_INT(a_idir),
                                    CHF_CONST_INT(iside),
                                    CHF_BOX(box));
      }
    else
      {
        if (m_isFunctional)
          {
            Real ihdx;

            ihdx = 2.0 / a_dx[a_idir];

            BoxIterator bit(box);

            for (bit.begin(); bit.ok(); ++bit)
              {
                IntVect iv = bit();
                IntVect ivNeigh = iv;
                ivNeigh[a_idir] += sign(a_side);
                const VolIndex vof      = VolIndex(iv,     0);
                const VolIndex vofNeigh = VolIndex(ivNeigh,0);
                const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
                const RealVect  point = EBArith::getFaceLocation(face,a_dx,a_probLo);
                const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
                Real value = m_func->value(face,a_side,a_dit,point,normal,a_time,comp);
                Real phi0 = a_phi(iv,comp);
                iv[a_idir] += iside;
                Real phi1 = a_phi(iv,comp);
                iv[a_idir] -= iside;
                a_faceFlux(iv,comp) = -iside*(8.0*value + phi1 - 9.0*phi0)/(3.0*a_dx[a_idir]);
              }
          }
        else
          {
            if (m_onlyHomogeneous)
              {
                MayDay::Error("DirichletPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
              }

            Real value = m_value;

            FORT_SETHODIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
                                        CHF_CONST_FRA1(a_phi,comp),
                                        CHF_CONST_REAL(value),
                                        CHF_CONST_REALVECT(a_dx),
                                        CHF_CONST_INT(a_idir),
                                        CHF_CONST_INT(iside),
                                        CHF_BOX(box));
          }
      }
  }
}
///
/**
   This is called by InflowOutflowHelmholtzDomainBC::getFaceFlux,
   which is called by EBAMRPoissonOp::applyOp when EB x domain.
   Used for higher-order velocity boundary conditions in viscous operator.
   This routine iterates over possible multiple faces at a domain face
   and calls getHigherOrderFaceFlux(...) below for the individual face extrapolation
*/
void DirichletPoissonDomainBC::getHigherOrderFaceFlux(Real&                 a_faceFlux,
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
  const ProblemDomain& domainBox = ebisBox.getDomain();
  a_faceFlux = 0.0;
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
  if (faces.size() > 0)
    {
      if (faces.size()==1)
        {
          IntVectSet cfivs;
          FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                           cfivs,
                                                           ebisBox,
                                                           domainBox);
          for (int isten=0; isten < faceSten.size(); isten++)
            {
              const Real& weight = faceSten.weight(isten);
              const FaceIndex& face = faceSten.face(isten);
              Real thisFaceFlux;
              const RealVect centroid = ebisBox.centroid(face);
              getHigherOrderFaceFlux(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
                                     a_side,a_dit,a_time,false,centroid,a_useHomogeneous);

              a_faceFlux += thisFaceFlux*weight;
            }
          a_faceFlux *= ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
          //could have done something like before by adding a division by number of faces but this does not use faceStencil above
          // int numFaces = 0;
          // for (int i = 0; i < faces.size(); i++)
          //   {
          //     const RealVect centroid = ebisBox.centroid(faces[i]);
          //     Real thisFaceFlux;
          //     getHigherOrderFaceFlux(thisFaceFlux,faces[i],a_comp,a_phi,a_probLo,
          //                            a_dx,a_idir,a_side,a_dit,a_time,false,centroid,
          //                            a_useHomogeneous);
          //     a_faceFlux += thisFaceFlux;
          //     numFaces ++;
          //   }
          // if (numFaces > 1){a_faceFlux /= Real(numFaces);}
        }
    }
}

///
/**
   This is called by DirichletPoissonDomain::getHigherOrderFaceflux
   which is called by InflowOutflowHelmholtzDomainBC::getFaceFlux.
   a_useHomogeneous now indicates if using whole stencil or not
   if true, then just get the inhomogeneous value;
   if false, then use whole stencil
   should not matter since this is only called for the velocity at operator constructor time (defineStencils)
*/
void DirichletPoissonDomainBC::getHigherOrderFaceFlux(Real&                 a_faceFlux,
                                                      const FaceIndex&      a_face,
                                                      const int&            a_comp,
                                                      const EBCellFAB&      a_phi,
                                                      const RealVect&       a_probLo,
                                                      const RealVect&       a_dx,
                                                      const int&            a_idir,
                                                      const Side::LoHiSide& a_side,
                                                      const DataIndex&      a_dit,
                                                      const Real&           a_time,
                                                      const bool&           a_useAreaFrac,
                                                      const RealVect&       a_centroid,
                                                      const bool&           a_useHomogeneous)
{
  int iside = -sign(a_side);
  const Real ihdx = 2.0 / a_dx[a_idir];
  const EBISBox& ebisBox = a_phi.getEBISBox();

  Real value = -1.e99;
  if (a_useHomogeneous)
    {
      value = 0.0;
    }
  else if (m_isFunctional)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      value = m_func->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux called with undefined inhomogeneous BC");
        }
      value = m_value;
    }

  const VolIndex& vof = a_face.getVoF(flip(a_side));
  //higher-order Dirichlet bc
  Vector<FaceIndex> facesInsideDomain = ebisBox.getFaces(vof,a_idir,flip(a_side));
  if (facesInsideDomain.size() == 1)
    {
      const VolIndex& vofNextInsideDomain = facesInsideDomain[0].getVoF(flip(a_side));
          a_faceFlux = iside*(9.0*a_phi(vof,a_comp) - a_phi(vofNextInsideDomain,a_comp) - 8.0*value)/(3.0*a_dx[a_idir]);
    }
  else
    {
      a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);
    }
}
///
/**
   This is called by InflowOutflowHelmholtzDomainBC::getFaceFlux,
   which is called by EBAMRPoissonOp::applyOp when EB x domain.
   Used for higher-order velocity boundary conditions in viscous operator.
   This routine iterates over possible multiple faces at a domain face
   and calls getHigherOrderFaceFlux(...) below for the individual face extrapolation
*/
void DirichletPoissonDomainBC::getHigherOrderInhomFaceFlux(Real&                 a_faceFlux,
                                                           const VolIndex&       a_vof,
                                                           const int&            a_comp,
                                                           const EBCellFAB&      a_phi,
                                                           const RealVect&       a_probLo,
                                                           const RealVect&       a_dx,
                                                           const int&            a_idir,
                                                           const Side::LoHiSide& a_side,
                                                           const DataIndex&      a_dit,
                                                           const Real&           a_time)
{
  const EBISBox& ebisBox = a_phi.getEBISBox();
  const ProblemDomain& domainBox = ebisBox.getDomain();
  a_faceFlux = 0.0;
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
  if (faces.size() > 0)
    {
      if (faces.size()==1)
        {
          IntVectSet cfivs;
          FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                           cfivs,
                                                           ebisBox,
                                                           domainBox);
          for (int isten=0; isten < faceSten.size(); isten++)
            {
              const Real& weight = faceSten.weight(isten);
              const FaceIndex& face = faceSten.face(isten);
              Real thisFaceFlux;
              const RealVect centroid = ebisBox.centroid(face);
              getHigherOrderInhomFaceFlux(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
                                          a_side,a_dit,a_time,false,centroid);

              a_faceFlux += thisFaceFlux*weight;
            }
          a_faceFlux *= ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
          //could have done something like before by adding a division by number of faces but this does not use faceStencil above
          // int numFaces = 0;
          // for (int i = 0; i < faces.size(); i++)
          //   {
          //     const RealVect centroid = ebisBox.centroid(faces[i]);
          //     Real thisFaceFlux;
          //     getHigherOrderFaceFlux(thisFaceFlux,faces[i],a_comp,a_phi,a_probLo,
          //                            a_dx,a_idir,a_side,a_dit,a_time,false,centroid,
          //                            a_useHomogeneous);
          //     a_faceFlux += thisFaceFlux;
          //     numFaces ++;
          //   }
          // if (numFaces > 1){a_faceFlux /= Real(numFaces);}
        }
    }
}

///
/**
   This is called by DirichletPoissonDomain::getHigherOrderFaceflux
   which is called by InflowOutflowHelmholtzDomainBC::getFaceFlux.
   a_useHomogeneous now indicates if using whole stencil or not
   if true, then just get the inhomogeneous value;
   if false, then use whole stencil
   should not matter since this is only called for the velocity at operator constructor time (defineStencils)
*/
void DirichletPoissonDomainBC::getHigherOrderInhomFaceFlux(Real&                 a_faceFlux,
                                                           const FaceIndex&      a_face,
                                                           const int&            a_comp,
                                                           const EBCellFAB&      a_phi,
                                                           const RealVect&       a_probLo,
                                                           const RealVect&       a_dx,
                                                           const int&            a_idir,
                                                           const Side::LoHiSide& a_side,
                                                           const DataIndex&      a_dit,
                                                           const Real&           a_time,
                                                           const bool&           a_useAreaFrac,
                                                           const RealVect&       a_centroid)
{
  int iside = -sign(a_side);
  const Real ihdx = 2.0 / a_dx[a_idir];
  const EBISBox& ebisBox = a_phi.getEBISBox();

  Real value = -1.e99;
  if (m_onlyHomogeneous)
    {
      value = 0.0;
    }
  else if (m_isFunctional)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      value = m_func->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("DirichletPoissonDomainBC::getHigherOrderInhomFaceFlux called with undefined inhomogeneous BC");
        }
      value = m_value;
    }

  const VolIndex& vof = a_face.getVoF(flip(a_side));
  //higher-order Dirichlet bc
  Vector<FaceIndex> facesInsideDomain = ebisBox.getFaces(vof,a_idir,flip(a_side));
  if (facesInsideDomain.size() == 1)
    {
      a_faceFlux = iside*(- 8.0*value)/(3.0*a_dx[a_idir]);
    }
  else
    {
      a_faceFlux = iside * ihdx * (- value);
    }
}
///
/**
   This is called by InflowOutflowPoissonDomainBC::getFaceFlux,
   which is called by EBAMRPoissonOp::applyDomainFlux for applyOp in reg cells.
   For Dirichlet boundary conditions on pressure (usually at outflow) in projection.
*/
void DirichletPoissonDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                           const BaseFab<Real>&  a_phi,
                                           const RealVect&       a_probLo,
                                           const RealVect&       a_dx,
                                           const int&            a_idir,
                                           const Side::LoHiSide& a_side,
                                           const DataIndex&      a_dit,
                                           const Real&           a_time,
                                           const bool&           a_useHomogeneous)
{
  for (int comp=0; comp<a_phi.nComp(); comp++)
    {
      const Box& box = a_faceFlux.box();

      int iside;

      if (a_side == Side::Lo)
        {
          iside = 1;
        }
      else
        {
          iside = -1;
        }

      if (a_useHomogeneous)
        {
          Real value = 0.0;

          // FORT_SETHODIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
          //                             CHF_CONST_FRA1(a_phi,comp),
          //                             CHF_CONST_REAL(value),
          //                             CHF_CONST_REALVECT(a_dx),
          //                             CHF_CONST_INT(a_idir),
          //                             CHF_CONST_INT(iside),
          //                             CHF_BOX(box));
          FORT_SETDIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
                                    CHF_CONST_FRA1(a_phi,comp),
                                    CHF_CONST_REAL(value),
                                    CHF_CONST_REALVECT(a_dx),
                                    CHF_CONST_INT(a_idir),
                                    CHF_CONST_INT(iside),
                                    CHF_BOX(box));
        }
      else
        {
          if (m_isFunctional)
            {
              Real ihdx;

              ihdx = 2.0 / a_dx[a_idir];

              BoxIterator bit(box);

              // for (bit.begin(); bit.ok(); ++bit)
              //   {
              //     IntVect iv = bit();
              //     IntVect ivNeigh = iv;
              //     IntVect ivNextNeigh = iv;
              //     ivNeigh[a_idir] += sign(a_side);
              //     ivNextNeigh[a_idir] += 2*sign(a_side);
              //     const VolIndex vof      = VolIndex(iv,     0);
              //     const VolIndex vofNeigh = VolIndex(ivNeigh,0);
              //     const VolIndex vofNextNeigh = VolIndex(ivNextNeigh,0);
              //     const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
              //     const FaceIndex faceNeigh = FaceIndex(vofNeigh,vofNextNeigh,a_idir);
              //     const RealVect  point = EBArith::getFaceLocation(face,a_dx,a_probLo);
              //     const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
              //     Real value = m_func->value(face,a_side,a_dit,point,normal,a_time,comp);
              //     Real phi0 = a_phi(iv,comp);
              //     iv[a_idir] += iside;
              //     Real phi1 = a_phi(iv,comp);
              //     iv[a_idir] -= iside;
              //     a_faceFlux(iv,comp) = -iside*(8.0*value + phi1 - 9.0*phi0)/(3.0*a_dx[a_idir]);
              //   }
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  IntVect ivNeigh = iv;
                  ivNeigh[a_idir] += sign(a_side);
                  const VolIndex vof      = VolIndex(iv,     0);
                  const VolIndex vofNeigh = VolIndex(ivNeigh,0);
                  const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
                  const RealVect  point = EBArith::getFaceLocation(face,a_dx,a_probLo);
                  const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
                  Real value = m_func->value(face,a_side,a_dit,point,normal,a_time,comp);
                  Real phiVal = a_phi(iv,comp);
                  a_faceFlux(iv,comp) = iside * ihdx * (phiVal - value);
                }
            }
          else
            {
              if (m_onlyHomogeneous)
                {
                  MayDay::Error("DirichletPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
                }

              Real value = m_value;

              // FORT_SETHODIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
              //                             CHF_CONST_FRA1(a_phi,comp),
              //                             CHF_CONST_REAL(value),
              //                             CHF_CONST_REALVECT(a_dx),
              //                             CHF_CONST_INT(a_idir),
              //                             CHF_CONST_INT(iside),
              //                             CHF_BOX(box));
              FORT_SETDIRICHLETFACEFLUX(CHF_FRA1(a_faceFlux,comp),
                                        CHF_CONST_FRA1(a_phi,comp),
                                        CHF_CONST_REAL(value),
                                        CHF_CONST_REALVECT(a_dx),
                                        CHF_CONST_INT(a_idir),
                                        CHF_CONST_INT(iside),
                                        CHF_BOX(box));
            }
        }
    }
}
///
/**
   This is called by InflowOutflowPoissonDomainBC::getFaceFlux,
   which is called by EBAMRPoissonOp::applyOp when EB x domain.
   For domain boundary conditions on pressure in projections.
*/
void DirichletPoissonDomainBC::getFaceFlux(Real&                 a_faceFlux,
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
  const ProblemDomain& domainBox = ebisBox.getDomain();
  a_faceFlux = 0.0;
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
  if (faces.size() > 0)
    {
      if (faces.size()==1)
        {
          IntVectSet cfivs;
          FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                           cfivs,
                                                           ebisBox,
                                                           domainBox);
          for (int isten=0; isten < faceSten.size(); isten++)
            {
              const Real& weight = faceSten.weight(isten);
              const FaceIndex& face = faceSten.face(isten);
              Real thisFaceFlux;
              const RealVect centroid = ebisBox.centroid(face);
              getFaceGradPhi(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
                             a_side,a_dit,a_time,false,centroid,a_useHomogeneous);

              a_faceFlux += thisFaceFlux*weight;
            }
          a_faceFlux *= ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
          //could have done something like before by adding a division by number of faces but this does not use faceStencil above
          // int numFaces = 0;
          // for (int i = 0; i < faces.size(); i++)
          //   {
          //     const RealVect centroid = ebisBox.centroid(faces[i]);
          //     Real thisFaceFlux;
          //     getFaceGradPhi(thisFaceFlux,faces[i],a_comp,a_phi,a_probLo,
          //                    a_dx,a_idir,a_side,a_dit,a_time,false,centroid,a_useHomogeneous);
          //     a_faceFlux += thisFaceFlux;
          // numFaces ++;
          //   }
          // if (numFaces > 1){a_faceFlux /= Real(numFaces);}
        }
    }
}

///
/**
   This is called by InflowOutflowPoissonDomainBC::getFaceGradPhi
   when enforcing boundary gradients in projection (called by macEnforceGradientBC).
   if true, then just get the inhomogeneous value;
   if false, then use whole stencil
   does matter because this is not only called for the pressure at operator constructor time (defineStencils)
   but also by macEnforceGradientBC which needs the whole stencil (can stencilize this later)
*/
void DirichletPoissonDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
                                              const FaceIndex&      a_face,
                                              const int&            a_comp,
                                              const EBCellFAB&      a_phi,
                                              const RealVect&       a_probLo,
                                              const RealVect&       a_dx,
                                              const int&            a_idir,
                                              const Side::LoHiSide& a_side,
                                              const DataIndex&      a_dit,
                                              const Real&           a_time,
                                              const bool&           a_useAreaFrac,
                                              const RealVect&       a_centroid,
                                              const bool&           a_useHomogeneous)
{
  int iside = -sign(a_side);
  const Real ihdx = 2.0 / a_dx[a_idir];

  Real value = -1.e99;
  if (a_useHomogeneous)
    {
      value = 0.0;
    }
  else if (m_isFunctional)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      value = m_func->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("DirichletPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
        }
      value = m_value;
    }

  const VolIndex& vof = a_face.getVoF(flip(a_side));

  //higher-order Dirichlet bc
  // Vector<FaceIndex> facesInsideDomain = ebisBox.getFaces(vof,a_idir,flip(a_side));
  // if (facesInsideDomain.size() == 1)
  //   {
  //     const VolIndex& vofNextInsideDomain = facesInsideDomain[0].getVoF(flip(a_side));
  //     a_faceFlux = iside*(9.0*a_phi(vof,a_comp) - a_phi(vofNextInsideDomain,a_comp) - 8.0*value)/(3.0*a_dx[a_idir]);
  //   }
  // else
  //   {
  //     a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);
  //   }

  a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);

}

///
/**
   This is called by InflowOutflowPoissonDomainBC::getFaceFlux,
   which is called by EBAMRPoissonOp::applyOp when EB x domain.
   For domain boundary conditions on pressure in projections.
*/
void DirichletPoissonDomainBC::getInhomFaceFlux(Real&                 a_faceFlux,
                                                const VolIndex&       a_vof,
                                                const int&            a_comp,
                                                const EBCellFAB&      a_phi,
                                                const RealVect&       a_probLo,
                                                const RealVect&       a_dx,
                                                const int&            a_idir,
                                                const Side::LoHiSide& a_side,
                                                const DataIndex&      a_dit,
                                                const Real&           a_time)
{
  const EBISBox& ebisBox = a_phi.getEBISBox();
  const ProblemDomain& domainBox = ebisBox.getDomain();
  a_faceFlux = 0.0;
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
  if (faces.size() > 0)
    {
      if (faces.size()==1)
        {
          IntVectSet cfivs;
          FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                           cfivs,
                                                           ebisBox,
                                                           domainBox);
          for (int isten=0; isten < faceSten.size(); isten++)
            {
              const Real& weight = faceSten.weight(isten);
              const FaceIndex& face = faceSten.face(isten);
              Real thisFaceFlux;
              const RealVect centroid = ebisBox.centroid(face);
              getInhomFaceGradPhi(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
                                  a_side,a_dit,a_time,false,centroid);

              a_faceFlux += thisFaceFlux*weight;
            }
          a_faceFlux *= ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("DirichletPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
          //could have done something like before by adding a division by number of faces but this does not use faceStencil above
          // int numFaces = 0;
          // for (int i = 0; i < faces.size(); i++)
          //   {
          //     const RealVect centroid = ebisBox.centroid(faces[i]);
          //     Real thisFaceFlux;
          //     getFaceGradPhi(thisFaceFlux,faces[i],a_comp,a_phi,a_probLo,
          //                    a_dx,a_idir,a_side,a_dit,a_time,false,centroid,a_useHomogeneous);
          //     a_faceFlux += thisFaceFlux;
          // numFaces ++;
          //   }
          // if (numFaces > 1){a_faceFlux /= Real(numFaces);}
        }
    }
}

///
/**
   This is called by InflowOutflowPoissonDomainBC::getFaceGradPhi
   when enforcing boundary gradients in projection (called by macEnforceGradientBC).
   if true, then just get the inhomogeneous value;
   if false, then use whole stencil
   does matter because this is not only called for the pressure at operator constructor time (defineStencils)
   but also by macEnforceGradientBC which needs the whole stencil (can stencilize this later)
*/
void DirichletPoissonDomainBC::getInhomFaceGradPhi(Real&                 a_faceFlux,
                                                   const FaceIndex&      a_face,
                                                   const int&            a_comp,
                                                   const EBCellFAB&      a_phi,
                                                   const RealVect&       a_probLo,
                                                   const RealVect&       a_dx,
                                                   const int&            a_idir,
                                                   const Side::LoHiSide& a_side,
                                                   const DataIndex&      a_dit,
                                                   const Real&           a_time,
                                                   const bool&           a_useAreaFrac,
                                                   const RealVect&       a_centroid)
{
  int iside = -sign(a_side);
  const Real ihdx = 2.0 / a_dx[a_idir];

  Real value = -1.e99;
  if (m_onlyHomogeneous)
    {
      value = 0.0;
    }
  else if (m_isFunctional)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      value = m_func->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      value = m_value;
    }

  //higher-order Dirichlet bc
  // Vector<FaceIndex> facesInsideDomain = ebisBox.getFaces(vof,a_idir,flip(a_side));
  // if (facesInsideDomain.size() == 1)
  //   {
  //     const VolIndex& vofNextInsideDomain = facesInsideDomain[0].getVoF(flip(a_side));
  //     a_faceFlux = iside*(9.0*a_phi(vof,a_comp) - a_phi(vofNextInsideDomain,a_comp) - 8.0*value)/(3.0*a_dx[a_idir]);
  //   }
  // else
  //   {
  //     a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);
  //   }

  a_faceFlux = iside * ihdx * (- value);

}

void DirichletPoissonDomainBC::getFluxStencil(      VoFStencil&      a_stencil,
                                              const VolIndex&        a_vof,
                                              const int&             a_comp,
                                              const RealVect&        a_dx,
                                              const int&             a_idir,
                                              const Side::LoHiSide&  a_side,
                                              const EBISBox&         a_ebisBox)
{
  const ProblemDomain& domainBox = a_ebisBox.getDomain();
  Vector<FaceIndex> faces = a_ebisBox.getFaces(a_vof, a_idir, a_side);
  if (faces.size() > 0)
    {
      if (faces.size()==1)
        {
          CH_assert(faces[0].isBoundary());
          IntVectSet cfivs;
          FaceStencil faceSten = EBArith::getInterpStencil(faces[0],
                                                           cfivs,
                                                           a_ebisBox,
                                                           domainBox);
          for (int isten=0; isten < faceSten.size(); isten++)
            {
              const Real& weight = faceSten.weight(isten);
              const FaceIndex& face = faceSten.face(isten);

              VoFStencil thisStencil;
              getFluxStencil(thisStencil, face, a_comp, a_dx, a_idir, a_side, a_ebisBox);
              thisStencil *= weight;
              a_stencil += thisStencil;
            }
          a_stencil *= a_ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("DirichletPoissonDomainBC::getFluxStencil has multi-valued faces");
        }
    }
  a_stencil *= 1./a_dx[a_idir];
  a_stencil *= Real(sign(a_side));//this sign is for which side is differenced
}

void DirichletPoissonDomainBC::getFluxStencil(      VoFStencil&      a_stencil,
                                              const FaceIndex&       a_face,
                                              const int&             a_comp,
                                              const RealVect&        a_dx,
                                              const int&             a_idir,
                                              const Side::LoHiSide&  a_side,
                                              const EBISBox&         a_ebisBox)
{
  if (m_ebOrder == 1)//pressure
    {
      getFirstOrderFluxStencil(a_stencil, a_face, a_comp,
                               a_dx, a_idir, a_side,
                               a_ebisBox);
    }
  else if (m_ebOrder == 2)//velocity
    {
      getSecondOrderFluxStencil(a_stencil, a_face, a_comp,
                                a_dx, a_idir, a_side,
                                a_ebisBox);
    }
  else
    {
      MayDay::Error("DirichletPoissonDomainBC::getFluxStencil -- bad BC order");
    }
}

void DirichletPoissonDomainBC::getFirstOrderFluxStencil(      VoFStencil&      a_stencil,
                                                        const FaceIndex&       a_face,
                                                        const int&             a_comp,
                                                        const RealVect&        a_dx,
                                                        const int&             a_idir,
                                                        const Side::LoHiSide&  a_side,
                                                        const EBISBox&         a_ebisBox)
{
  const VolIndex& vof = a_face.getVoF(flip(a_side));
  const Real isign = Real(sign(a_side));//this sign is for the extrapolation direction
  const Real weight = -isign*2.0/a_dx[a_idir];
  a_stencil.add(vof, weight, a_comp);
}

void DirichletPoissonDomainBC::getSecondOrderFluxStencil(      VoFStencil&      a_stencil,
                                                         const FaceIndex&       a_face,
                                                         const int&             a_comp,
                                                         const RealVect&        a_dx,
                                                         const int&             a_idir,
                                                         const Side::LoHiSide&  a_side,
                                                         const EBISBox&         a_ebisBox)
{
  const VolIndex& vof = a_face.getVoF(flip(a_side));
  const Real isign = Real(sign(a_side));//this sign is for the extrapolation direction
  const Real weight = isign*1.0/(3.0*a_dx[a_idir]);
  Vector<FaceIndex> facesInsideDomain = a_ebisBox.getFaces(vof,a_idir,flip(a_side));
  if (facesInsideDomain.size() == 1)
    {
      const VolIndex& vofNextInsideDomain = facesInsideDomain[0].getVoF(flip(a_side));
      a_stencil.add(vofNextInsideDomain,      weight, a_comp);
      a_stencil.add(                vof, -9.0*weight, a_comp);
    }
  else
    {
      a_stencil.add(vof, -6.0*weight, a_comp);
    }
}

DirichletPoissonDomainBCFactory::DirichletPoissonDomainBCFactory()
: m_onlyHomogeneous(true),
  m_isFunctional(false),
  m_value(12345.6789),
  m_func(RefCountedPtr<BaseBCValue>()),
  m_ebOrder(2)
{
}

DirichletPoissonDomainBCFactory::~DirichletPoissonDomainBCFactory()
{
}

void DirichletPoissonDomainBCFactory::setValue(Real a_value)
{
  m_value = a_value;
  m_func = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunctional = false;
}

void DirichletPoissonDomainBCFactory::setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunctional = true;
}

//use this for order of domain boundary
void DirichletPoissonDomainBCFactory::setEBOrder(int a_ebOrder)
{
  CH_assert(a_ebOrder >= 1);
  CH_assert(a_ebOrder <= 2);

  m_ebOrder = a_ebOrder;
}

DirichletPoissonDomainBC*
DirichletPoissonDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  DirichletPoissonDomainBC* newBC = new DirichletPoissonDomainBC();
  newBC->setEBOrder(m_ebOrder);
  if (!m_onlyHomogeneous)
    {
      if (m_isFunctional)
        {
          newBC->setFunction(m_func);
        }
      else
        {
          newBC->setValue(m_value);
        }
    }
  return newBC;
}
#include "NamespaceFooter.H"
