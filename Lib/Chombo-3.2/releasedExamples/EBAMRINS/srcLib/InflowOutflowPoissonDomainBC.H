#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INFLOWOUTFLOWPOISSONDOMAINBC_H_
#define _INFLOWOUTFLOWPOISSONDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "PoiseuilleInflowBCValue.H"

#include "NamespaceHeader.H"

///
/**
 */
class InflowOutflowPoissonDomainBC: public BaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~InflowOutflowPoissonDomainBC()
  {;}

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const VolIndex&        a_vof,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

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

  ///
  /**
     Implementation of the query of whether the face between a_ivof and a_jvof Dirichlet Domain boundary
     For InflowOutflowPoissonDomainBC, a_phi and a_dit are not used.
     For the moment it is only valid for the pressure calculation
     a_ivof must be inside the valid domain
  */
  virtual bool isDirichletDom(const VolIndex&   a_ivof,
                              const VolIndex&   a_jvof,
                              const EBCellFAB&  a_phi) const;
  ///
  /**
     Elliptic solver inhom flux.
  */
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
     velocity at domain faces.  Sets velocity to a InflowOutflow  value.
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
  InflowOutflowPoissonDomainBC(const int& a_flowDir,
                               const Real& a_inflowVel,
                               const bool& a_doPoiseInflow,
                               const IntVect& a_doSlipWallsHi,
                               const IntVect& a_doSlipWallsLo,
                               const ProblemDomain& a_domain,
                               const RealVect& a_dx,
                               RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                               const bool& a_doWomersleyInflow=false)
  {
    m_flowDir = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_doPoiseInflow = a_doPoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_domain = a_domain;
    m_dx = a_dx;
    m_poiseInflowFunc = a_poiseBCValue;
    m_doWomersleyInflow = a_doWomersleyInflow;
  }

private:

  ///weak construction forbidden to keep things simple
  InflowOutflowPoissonDomainBC()
  {
  }

  int m_flowDir;
  Real m_inflowVel;
  bool m_doPoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  ProblemDomain m_domain;
  RealVect m_dx;
  bool m_doWomersleyInflow;
  RefCountedPtr<PoiseuilleInflowBCValue> m_poiseInflowFunc;

};

class InflowOutflowPoissonDomainBCFactory: public BaseDomainBCFactory
{
public:


  ///
  /**
   */
  InflowOutflowPoissonDomainBCFactory(const int& a_flowDir,
                                      const Real& a_inflowVel,
                                      const bool& a_doPoiseInflow,
                                      const IntVect& a_doSlipWallsHi,
                                      const IntVect& a_doSlipWallsLo,
                                      RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                                      const bool& a_doWomersleyInflow=false)
  {
    m_flowDir   = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_doPoiseInflow = a_doPoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_doWomersleyInflow = a_doWomersleyInflow;
    m_poiseInflowFunc = a_poiseBCValue;
  }


  ///
  /**
   */
  virtual ~InflowOutflowPoissonDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual InflowOutflowPoissonDomainBC* create(const ProblemDomain& a_domain,
                                               const EBISLayout&    a_layout,
                                               const RealVect&      a_dx)
  {
    InflowOutflowPoissonDomainBC* newBC
      = new InflowOutflowPoissonDomainBC(m_flowDir,
                                         m_inflowVel,
                                         m_doPoiseInflow,
                                         m_doSlipWallsHi,
                                         m_doSlipWallsLo,
                                         a_domain,
                                         a_dx,
                                         m_poiseInflowFunc,
                                         m_doWomersleyInflow);

    return newBC;
  }


private:
  InflowOutflowPoissonDomainBCFactory()
  {
  }

  int m_flowDir;
  Real m_inflowVel;
  bool m_doPoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  bool m_doWomersleyInflow;
  RefCountedPtr<PoiseuilleInflowBCValue> m_poiseInflowFunc;

};



///
/**
 */
class InflowOutflowHelmholtzDomainBC: public BaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~InflowOutflowHelmholtzDomainBC()
  {;}

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const VolIndex&        a_vof,
                              const int&             a_comp,
                              const RealVect&        a_dx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

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

  ///
  /**
     Elliptic solver inhom flux.
  */
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

  virtual bool isDirichletDom(const VolIndex&   a_ivof,
                              const VolIndex&   a_jvof,
                              const EBCellFAB&  a_phi) const;
  ///
  /**
     This is used by the projections to get
     velocity at domain faces.  Sets velocity to a InflowOutflow  value.
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
  InflowOutflowHelmholtzDomainBC(const int& a_flowDir,
                                 const Real& a_inflowVel,
                                 const bool& a_doPoiseInflow,
                                 const IntVect& a_doSlipWallsHi,
                                 const IntVect& a_doSlipWallsLo,
                                 const ProblemDomain& a_domain,
                                 const RealVect& a_dx,
                                 RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                                 const bool& a_doWomersleyInflow=false)
  {
    m_flowDir   = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_doPoiseInflow = a_doPoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_domain = a_domain;
    m_dx = a_dx;
    m_poiseInflowFunc = a_poiseBCValue;
    m_doWomersleyInflow = a_doWomersleyInflow;
  }

protected:
  //weak construction forbidden
  InflowOutflowHelmholtzDomainBC()
  {
    MayDay::Error("invalid operator");
  }

  int m_flowDir;
  Real m_inflowVel;
  bool m_doPoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  ProblemDomain m_domain;
  RealVect m_dx;
  bool m_doWomersleyInflow;
  RefCountedPtr<PoiseuilleInflowBCValue> m_poiseInflowFunc;
};

///
/**
 **/
class InflowOutflowHelmholtzDomainBCFactory: public BaseDomainBCFactory
{
public:

  ///
  InflowOutflowHelmholtzDomainBCFactory(const int& a_flowDir,
                                        const Real& a_inflowVel,
                                        const bool& a_doPoiseInflow,
                                        const IntVect& a_doSlipWallsHi,
                                        const IntVect& a_doSlipWallsLo,
                                        RefCountedPtr<PoiseuilleInflowBCValue> a_poiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue> (),
                                        const bool& a_doWomersleyInflow=false)
  {
    m_flowDir   = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_doPoiseInflow = a_doPoiseInflow;
    m_doSlipWallsHi = a_doSlipWallsHi;
    m_doSlipWallsLo = a_doSlipWallsLo;
    m_doWomersleyInflow = a_doWomersleyInflow;
    m_poiseInflowFunc = a_poiseBCValue;
  }


  ///
  /**
   */
  virtual ~InflowOutflowHelmholtzDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual InflowOutflowHelmholtzDomainBC* create(const ProblemDomain& a_domain,
                                                 const EBISLayout&    a_layout,
                                                 const RealVect&      a_dx)
  {
    InflowOutflowHelmholtzDomainBC* newBC
      = new InflowOutflowHelmholtzDomainBC(m_flowDir,
                                           m_inflowVel,
                                           m_doPoiseInflow,
                                           m_doSlipWallsHi,
                                           m_doSlipWallsLo,
                                           a_domain,
                                           a_dx,
                                           m_poiseInflowFunc,
                                           m_doWomersleyInflow);

    return newBC;
  }

protected:
  //weak construction
  InflowOutflowHelmholtzDomainBCFactory()
  {
    MayDay::Error("invalid operator");
  }

  int m_flowDir;
  Real m_inflowVel;
  bool m_doPoiseInflow;
  IntVect m_doSlipWallsHi;
  IntVect m_doSlipWallsLo;
  bool m_doWomersleyInflow;
  RefCountedPtr<PoiseuilleInflowBCValue> m_poiseInflowFunc;
};

#include "NamespaceFooter.H"
#endif
