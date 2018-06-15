#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "NeumannConductivityDomainBC.H"
#include "EBArithF_F.H"
#include "NamespaceHeader.H"
/*****/
NeumannConductivityDomainBC::
NeumannConductivityDomainBC()
{
}
/*****/
NeumannConductivityDomainBC::
~NeumannConductivityDomainBC()
{
}

/*****/
void
NeumannConductivityDomainBC::
setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}

/*****/
int
NeumannConductivityDomainBC::
whichBC(int                  a_idir,
        Side::LoHiSide       a_side)
{
  return 0;
}

/*****/
void
NeumannConductivityDomainBC::
setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}

/*****/
void
NeumannConductivityDomainBC::
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
NeumannConductivityDomainBC::
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

  Real totalMassFlux = 0.0;
  Vector<FaceIndex> faces = ebisBox.getFaces(a_vof,a_idir,a_side);
  for (int i = 0; i < faces.size(); i++)
    {
      const RealVect centroid = ebisBox.centroid(faces[i]);
      Real thisFaceFlux;
      getFaceGradPhi(thisFaceFlux,faces[i],a_comp,a_phi,a_probLo,
                         a_dx,a_idir,a_side,a_dit,a_time,true,centroid,a_useHomogeneous);
      totalMassFlux += thisFaceFlux;
    }

  a_faceFlux = totalMassFlux;

  Real bcoave = 0;
  Real areaTot = 0.0;
  for (int iface = 0; iface < faces.size(); iface++)
    {
      Real areaFrac = ebisBox.areaFrac(faces[iface]);
      areaTot += areaFrac;
      Real bcoFace  = (*m_bcoef)[a_dit][a_idir](faces[iface], 0);
      bcoave += areaFrac*bcoFace;
    }
  if (areaTot > 1.0e-8)
    {
      bcoave /= areaTot;
    }
  a_faceFlux *= bcoave;
}

/*****/
void
NeumannConductivityDomainBC::
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
  const int iside = -sign(a_side);

  const EBISBox& ebisBox = a_phi.getEBISBox();

  Real flux = -1.e99;
  if (a_useHomogeneous)
    {
      flux = 0.0;
    }
  else if (m_isFunction)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      RealVect point  = a_centroid;
      point *= a_dx;
      point[a_face.direction()] = 0.0;//make this not depend on whatever ebisbox is returning for centroid in the face direction.
      point += EBArith::getFaceLocation(a_face,a_dx,a_probLo);
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
  if (a_useAreaFrac)
    {
      a_faceFlux *= ebisBox.areaFrac(a_face);
    }
}

/*****/
void
NeumannConductivityDomainBC::
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
  a_faceFlux =  a_vel[a_idir](a_face, a_icomp);
}
/******/
NeumannConductivityDomainBCFactory::
NeumannConductivityDomainBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

/******/
NeumannConductivityDomainBCFactory::
~NeumannConductivityDomainBCFactory()
{
}
/******/
void
NeumannConductivityDomainBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}
/******/
void
NeumannConductivityDomainBCFactory::
setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
/******/
NeumannConductivityDomainBC*
NeumannConductivityDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
  NeumannConductivityDomainBC* newBC = new NeumannConductivityDomainBC();
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
