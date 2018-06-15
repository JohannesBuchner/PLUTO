#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DirichletPoissonDomainBCF_F.H"
#include "DirichletConductivityDomainBC.H"
#include "EBArithF_F.H"
#include "NamespaceHeader.H"
/*****/
DirichletConductivityDomainBC::
DirichletConductivityDomainBC()
{
}
/*****/
DirichletConductivityDomainBC::
~DirichletConductivityDomainBC()
{
}

/*****/
void
DirichletConductivityDomainBC::
setValue(Real a_value)
{
   m_onlyHomogeneous = false;
   m_isFunctional = false;
   m_value = a_value;
   m_func = RefCountedPtr<BaseBCValue>();
}

/*****/
int
DirichletConductivityDomainBC::
whichBC(int                  a_idir,
        Side::LoHiSide       a_side)
{
  return 0;
}

/*****/
void
DirichletConductivityDomainBC::
setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunctional = true;
}

/*****/
void
DirichletConductivityDomainBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
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

  //again, following the odd convention of EBAMRPoissonOp
  //(because I am reusing its BC classes),
  //the input flux here is CELL centered and the input box
  //is the box adjacent to the domain boundary on the valid side.
  //because I am not insane (yet) I will just shift the flux's box
  //over and multiply by the appropriate coefficient
  a_faceFlux.shiftHalf(a_idir, -sign(a_side));
  const Box& faceBox = a_faceFlux.box();
  const BaseFab<Real>&   regCoef = (*m_bcoef)[a_dit][a_idir].getSingleValuedFAB();
  int  isrc = 0;
  int  idst = 0;
  int  inum = 1;
  FORT_MULTIPLYTWOFAB(CHF_FRA(a_faceFlux),
                      CHF_CONST_FRA(regCoef),
                      CHF_BOX(faceBox),
                      CHF_INT(isrc),CHF_INT(idst),CHF_INT(inum));

  //shift flux back to cell centered land
  a_faceFlux.shiftHalf(a_idir,  sign(a_side));
}

/*****/
void

DirichletConductivityDomainBC::
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
        }
    }

  Real bcoave = 0;
  Real areaTot = 0;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      Real bcoFace  = (*m_bcoef)[a_dit][a_idir](faces[iface], 0);
      bcoave += areaFrac*bcoFace;
      areaTot += areaFrac;
    }
  if (areaTot > 1.0e-8)
    {
      bcoave /= areaTot;
    }
  a_faceFlux *= bcoave;
}

/*****/
void
DirichletConductivityDomainBC::
getFaceGradPhi(Real&                 a_faceFlux,
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

  if (a_useAreaFrac)
    {
      MayDay::Error("DirichletPoissonDomainBC::getFaceFlux -- useAreaFrac=TRUE");
      a_faceFlux *= ebisBox.areaFrac(a_face);
    }
}

/*****/
void
DirichletConductivityDomainBC::
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
/******/
DirichletConductivityDomainBCFactory::
DirichletConductivityDomainBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

/******/
DirichletConductivityDomainBCFactory::
~DirichletConductivityDomainBCFactory()
{
}
/******/
void
DirichletConductivityDomainBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}
/******/
void
DirichletConductivityDomainBCFactory::
setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
/******/
DirichletConductivityDomainBC*
DirichletConductivityDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  DirichletConductivityDomainBC* newBC = new DirichletConductivityDomainBC();
  if (!m_onlyHomogeneous)
    {
      if (m_isFunction)
        {
          newBC->setFunction(m_flux);
        }
      else
        {
          newBC->setValue(m_value);
        }
    }
  return newBC;
}
#include "NamespaceFooter.H"
