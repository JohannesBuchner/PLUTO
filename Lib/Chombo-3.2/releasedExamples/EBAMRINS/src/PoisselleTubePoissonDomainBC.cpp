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

#include "PoisselleTubePoissonDomainBC.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "VoFIterator.H"

#include "UsingNamespace.H"

void
PoisselleTubePoissonDomainBC::
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
  CH_assert(a_vel.nComp() == 1);
  bool isInflow = (a_side==Side::Lo);
  bool isOutflow= (a_side==Side::Hi);

  RealVect point =  EBArith::getFaceLocation(a_face, a_dx, a_probLo);
  RealVect normal = EBArith::getDomainNormal(a_idir, a_side);

  const EBISBox& ebisBox= a_vel[0].getEBISBox();

  if (isOutflow)
    {
      //quadratic extrapolation with Neumann bc
      a_faceFlux = EBArith::extrapFaceVelToOutflow(a_face,a_side,a_idir,ebisBox.getEBGraph(),a_vel[a_idir],0);
    }
  else if (isInflow)
    {
      //bogus values for normal and time because they do not matter
      //input component is always zero.
      //value of face direction corresponds to velocity direciton here
      a_faceFlux = m_bcval[a_idir].value(point, RealVect::Zero, a_time, a_idir);
    }
  else
    {
      MayDay::Error("PoisselleTubePoissonDomainBC::getFaceVel::unexpected case");
    }
}

////
void PoisselleTubePoissonDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                               const BaseFab<Real>&  a_phi,
                                               const RealVect&       a_probLo,
                                               const RealVect&       a_dx,
                                               const int&            a_idir,
                                               const Side::LoHiSide& a_side,
                                               const DataIndex&      a_dit,
                                               const Real&           a_time,
                                               const bool&           a_useHomogeneous)
{
  //for phi: outflow  dirichlet. inflow neumann
  bool isPhiNeum = (a_side == Side::Lo);
  if (isPhiNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}


void PoisselleTubePoissonDomainBC::getFaceFlux(Real&                 a_faceFlux,
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
  bool isPhiNeum = (a_side == Side::Lo);
  if (isPhiNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

void PoisselleTubePoissonDomainBC::getInhomFaceFlux(Real&                 a_faceFlux,
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
  bool isPhiNeum = (a_side == Side::Lo);
  if (isPhiNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                                 a_dx, a_idir, a_side, a_dit,a_time);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                              a_dx, a_idir, a_side, a_dit,a_time);
    }
}

void PoisselleTubePoissonDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
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
  bool isPhiNeum = (a_side == Side::Lo);

  if (isPhiNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceGradPhi(a_faceFlux, a_face, a_comp,a_phi, a_probLo,
                               a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceGradPhi(a_faceFlux, a_face, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);
    }
}

void PoisselleTubePoissonDomainBC::getInhomFaceGradPhi(Real&                 a_faceFlux,
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
  bool isPhiNeum = (a_side == Side::Lo);

  if (isPhiNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getInhomFaceGradPhi(a_faceFlux, a_face, a_comp,a_phi, a_probLo,
                                    a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getInhomFaceGradPhi(a_faceFlux, a_face, a_comp, a_phi, a_probLo,
                                 a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid);
    }
}


void
PoisselleTubeHelmholtzDomainBC::
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
  MayDay::Error("not needed");
}

void PoisselleTubeHelmholtzDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                                 const BaseFab<Real>&  a_phi,
                                                 const RealVect&       a_probLo,
                                                 const RealVect&       a_dx,
                                                 const int&            a_idir,
                                                 const Side::LoHiSide& a_side,
                                                 const DataIndex&      a_dit,
                                                 const Real&           a_time,
                                                 const bool&           a_useHomogeneous)
{
  //vel: outflow is neumann. all others dirichlet.  inflow uses inflow vel as the value
  bool isVelNeum = (a_side == Side::Hi);

  if (isVelNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      RefCountedPtr<BaseBCValue> bcv = RefCountedPtr<BaseBCValue>(new PoisselleTubeBCValue(m_bcval));
      diriBC.setFunction(bcv);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

void PoisselleTubeHelmholtzDomainBC::getFaceFlux(Real&                 a_faceFlux,
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
  //vel: outflow is Neumann. all others Dirichlet
  bool isVelNeum = (a_side == Side::Hi);

  if (isVelNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      RefCountedPtr<BaseBCValue> bcv =RefCountedPtr<BaseBCValue>(new PoisselleTubeBCValue(m_bcval));
      diriBC.setFunction(bcv);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
}

void PoisselleTubeHelmholtzDomainBC::getInhomFaceFlux(Real&                 a_faceFlux,
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
  //vel: outflow is Neumann. all others Dirichlet
  bool isVelNeum = (a_side == Side::Hi);

  if (isVelNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                 a_dx, a_idir, a_side, a_dit,a_time);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      RefCountedPtr<BaseBCValue> bcv =RefCountedPtr<BaseBCValue>(new PoisselleTubeBCValue(m_bcval));
      diriBC.setFunction(bcv);
      diriBC.getHigherOrderInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                         a_dx, a_idir, a_side, a_dit,a_time);
    }
}

void PoisselleTubeHelmholtzDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
                                                    const FaceIndex&       a_face,
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
  MayDay::Error("InOuAnyHelmholtz::getFaceGradPhi: not needed");
}
