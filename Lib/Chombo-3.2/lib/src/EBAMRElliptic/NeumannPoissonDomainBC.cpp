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
#include "PolyGeom.H"

#include "NeumannPoissonDomainBC.H"
#include "NamespaceHeader.H"

NeumannPoissonDomainBC::NeumannPoissonDomainBC()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

NeumannPoissonDomainBC::~NeumannPoissonDomainBC()
{
}

void NeumannPoissonDomainBC::setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}

void NeumannPoissonDomainBC::setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}

void
NeumannPoissonDomainBC::getFaceVel(Real&                 a_faceFlux,
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
  a_faceFlux =  a_vel[a_idir](a_face, a_icomp);
}

void NeumannPoissonDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                         const BaseFab<Real>&  a_phi,
                                         const RealVect&       a_probLo,
                                         const RealVect&       a_dx,
                                         const int&            a_idir,
                                         const Side::LoHiSide& a_side,
                                         const DataIndex&      a_dit,
                                         const Real&           a_time,
                                         const bool&           a_useHomogeneous)
{
  CH_TIME("NeumannPoissonDomainBC::getFaceFlux");
  CH_assert(a_phi.nComp() == 1);
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

          a_faceFlux.setVal(iside * value);
        }
      else
        {
          if (m_isFunction)
            {
              BoxIterator bit(box);

              for (bit.begin(); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();

                  RealVect point = EBArith::getIVLocation(iv,a_dx,a_probLo);
                  point[a_idir] -= iside * 0.5 * a_dx[a_idir];//point is now at the face center

                  RealVect normal = RealVect::Zero;
                  normal[a_idir] = iside;

                  a_faceFlux(iv,comp) = iside * m_flux->value(iv,a_dit,point,normal,a_time,comp);
                }
            }
          else
            {
              if (m_onlyHomogeneous)
                {
                  MayDay::Error("NeumannPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
                }

              Real value = m_value;

              a_faceFlux.setVal(iside * value);
            }
        }
    }
}

//This extrapolates all gradients to face centers at boundaries and then interpolates to face centroid in the case of eb x domain (called by the operator)
void NeumannPoissonDomainBC::getFaceFlux(Real&                 a_faceFlux,
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
              getFaceFluxGradPhi(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
                                 a_side,a_dit,a_time,false,centroid,a_useHomogeneous);

              a_faceFlux += thisFaceFlux*weight;
            }
          a_faceFlux *= ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("NeumannPoissonDomainBC::getHigherOrderFaceFlux has multi-valued faces (or could be 0 faces)");
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
//called by function above (operator) and below (projection gradient)
void NeumannPoissonDomainBC::getFaceFluxGradPhi(Real&                 a_faceFlux,
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
  const int iside = -sign(a_side);
  Real flux = -1.e99;
  if (a_useHomogeneous)
    {
      flux = 0.0;
    }
  else if (m_isFunction)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      flux = m_flux->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("NeumannPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
        }
      flux = m_value;
    }

  a_faceFlux = iside*flux;

}

//This is for the extrapolation of all gradients to face centers at boundaries in projection
void NeumannPoissonDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
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
  getFaceFluxGradPhi(a_faceFlux,a_face,a_comp,a_phi,a_probLo,
                     a_dx,a_idir,a_side,a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);

}

//This extrapolates all gradients to face centers at boundaries and then interpolates to face centroid in the case of eb x domain (called by the operator)
void NeumannPoissonDomainBC::getInhomFaceFlux(Real&                 a_faceFlux,
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
              getInhomFaceFluxGradPhi(thisFaceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
                                      a_side,a_dit,a_time,false,centroid);

              a_faceFlux += thisFaceFlux*weight;
            }
          a_faceFlux *= ebisBox.areaFrac(faces[0]);
        }
      else
        {
          MayDay::Error("NeumannPoissonDomainBC::getInhomFaceFlux has multi-valued faces (or could be 0 faces)");
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
//called by function above (operator) and below (projection gradient)
void NeumannPoissonDomainBC::getInhomFaceFluxGradPhi(Real&                 a_faceFlux,
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
  const int iside = -sign(a_side);
  Real flux = -1.e99;
  if (m_onlyHomogeneous)
    {
      flux = 0.0;
    }
  else if (m_isFunction)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);
      flux = m_flux->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      flux = m_value;
    }

  a_faceFlux = iside*flux;

}

//This is for the extrapolation of all gradients to face centers at boundaries in projection
void NeumannPoissonDomainBC::getInhomFaceGradPhi(Real&                 a_faceFlux,
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
  getInhomFaceFluxGradPhi(a_faceFlux,a_face,a_comp,a_phi,a_probLo,
                          a_dx,a_idir,a_side,a_dit,a_time,a_useAreaFrac,a_centroid);

}

NeumannPoissonDomainBCFactory::NeumannPoissonDomainBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

NeumannPoissonDomainBCFactory::~NeumannPoissonDomainBCFactory()
{
}

void NeumannPoissonDomainBCFactory::setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}

void NeumannPoissonDomainBCFactory::setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
NeumannPoissonDomainBC*
NeumannPoissonDomainBCFactory::create(const ProblemDomain& a_domain,
                                      const EBISLayout&    a_layout,
                                      const RealVect&      a_dx)
{
  NeumannPoissonDomainBC* newBC = new NeumannPoissonDomainBC();
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
