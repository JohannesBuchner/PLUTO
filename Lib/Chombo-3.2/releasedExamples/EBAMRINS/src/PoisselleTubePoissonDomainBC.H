#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POISSELLETUBEPOISSONDOMAINBC_H_
#define _POISSELLETUBEPOISSONDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCValue.H"
#include "PoisselleTubeBCValue.H"

#include "UsingNamespace.H"

///
/**
 */
class PoisselleTubePoissonDomainBC: public BaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~PoisselleTubePoissonDomainBC()
  {;}

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  virtual void getInhomFaceFlux(Real&                 a_faceFlux,
                                const VolIndex&       a_vof,
                                const int&            a_comp,
                                const EBCellFAB&      a_phi,
                                const RealVect&       a_probLo,
                                const RealVect&       a_dx,
                                const int&            a_idir,
                                const Side::LoHiSide& a_side,
                                const DataIndex&      a_dit,
                                const Real&           a_time);

  virtual void getFaceGradPhi(Real&                 a_faceFlux,
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
                              const bool&           a_useHomogeneous);

  void getInhomFaceGradPhi(Real&                 a_faceFlux,
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
                           const RealVect&       a_centroid);

  ///
  /**
     This is used by the projections to get
     velocity at domain faces.  Sets velocity to a PoisselleTube  value.
  */
  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow);

  ///
  PoisselleTubePoissonDomainBC(const Vector<PoisselleTubeBCValue>& a_bcval)
  {
    m_bcval = a_bcval;
  }

  ///
  /**
     Return the a face centered velocity that makes a_vof div-free
  */
  static Real getDivFreeFace(const FaceIndex&      a_face,
                             const VolIndex&       a_vof,
                             const int&            a_dir,
                             const Side::LoHiSide& a_side,
                             const EBFluxFAB&      a_vel,
                             const RealVect&       a_dx,
                             const EBISBox&        a_ebisBox);

private:

  ///weak construction forbidden to keep things simple
  PoisselleTubePoissonDomainBC()
  {
  }

  Vector<PoisselleTubeBCValue> m_bcval;

};

class PoisselleTubePoissonDomainBCFactory: public BaseDomainBCFactory
{
public:


  ///
  /**
   */
  PoisselleTubePoissonDomainBCFactory(const Vector<PoisselleTubeBCValue>& a_bcval)
  {
    m_bcval = a_bcval;
  }


  ///
  /**
   */
  virtual ~PoisselleTubePoissonDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual PoisselleTubePoissonDomainBC* create(const ProblemDomain& a_domain,
                                               const EBISLayout&    a_layout,
                                               const RealVect&      a_dx)
  {
    PoisselleTubePoissonDomainBC* newBC = new PoisselleTubePoissonDomainBC(m_bcval);

    return newBC;
  }

private:
  PoisselleTubePoissonDomainBCFactory()
  {
  }

  Vector<PoisselleTubeBCValue> m_bcval;
};



///
/**
 */
class PoisselleTubeHelmholtzDomainBC: public BaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~PoisselleTubeHelmholtzDomainBC()
  {;}

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  virtual void getInhomFaceFlux(Real&                 a_faceFlux,
                                const VolIndex&       a_vof,
                                const int&            a_comp,
                                const EBCellFAB&      a_phi,
                                const RealVect&       a_probLo,
                                const RealVect&       a_dx,
                                const int&            a_idir,
                                const Side::LoHiSide& a_side,
                                const DataIndex&      a_dit,
                                const Real&           a_time);

  virtual void getFaceGradPhi(Real&                 a_faceFlux,
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
                              const bool&           a_useHomogeneous);

  ///
  /**
     This is used by the projections to get
     velocity at domain faces.  Sets velocity to a PoisselleTube  value.
  */
  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow);


  ///
  PoisselleTubeHelmholtzDomainBC(const PoisselleTubeBCValue& a_bcval)
  {
    m_bcval = a_bcval;
  }

private:
  //weak construction forbidden
  PoisselleTubeHelmholtzDomainBC()
  {
  }

  PoisselleTubeBCValue m_bcval;
};

///
/**
 **/
class PoisselleTubeHelmholtzDomainBCFactory: public BaseDomainBCFactory
{
public:

  ///
  PoisselleTubeHelmholtzDomainBCFactory(const PoisselleTubeBCValue& a_bcval)
  {
    m_bcval = a_bcval;
  }


  ///
  /**
   */
  virtual ~PoisselleTubeHelmholtzDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual PoisselleTubeHelmholtzDomainBC* create(const ProblemDomain& a_domain,
                                                 const EBISLayout&    a_layout,
                                                 const RealVect&      a_dx)
  {
    PoisselleTubeHelmholtzDomainBC* newBC = new PoisselleTubeHelmholtzDomainBC(m_bcval);

    return newBC;
  }

private:
  //weak construction
  PoisselleTubeHelmholtzDomainBCFactory()
  {
  }

  PoisselleTubeBCValue m_bcval;
};


#endif
