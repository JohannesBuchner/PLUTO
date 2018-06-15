#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SchlichtingTempDomainBC_H_______
#define _SchlichtingTempDomainBC_H_______

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "EBAMRCNSParams.H"

#include "NamespaceHeader.H"

///
/**
 */
class SchlichtingTempDomainBC: public ConductivityBaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~SchlichtingTempDomainBC()
  {;}

  ///
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
                              const bool&           a_useHomogeneous)
  {
    //this is for projections and such.
    //not used here.
    MayDay::Error();
  }


  ///
  virtual void getFaceVel(Real&                 a_faceFlux,
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
    //this is for projections and such.
    //not used here.
    MayDay::Error();
  }


  ///
  /**
   */
  SchlichtingTempDomainBC(const SchlichtingParams& a_params)
  {
    m_params = a_params;
  }

private:

  ///weak construction forbidden to keep things simple
  SchlichtingTempDomainBC(){;};
  SchlichtingParams m_params;
};

class SchlichtingTempDomainBCFactory: public BaseDomainBCFactory
{
public:


  ///
  /**
   */
  SchlichtingTempDomainBCFactory(const SchlichtingParams& a_params)
  {
    m_params = a_params;
  }


  ///
  /**
   */
  virtual ~SchlichtingTempDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual SchlichtingTempDomainBC* create(const ProblemDomain& a_domain,
                                           const EBISLayout&    a_layout,
                                           const RealVect&      a_dx)
  {
    SchlichtingTempDomainBC* newBC = new SchlichtingTempDomainBC(m_params);

    return newBC;
  }


private:

  //insert incoherent rant about weak construction here
  //(it really does make things simpler)
  SchlichtingTempDomainBCFactory(){;}


  SchlichtingParams m_params;
};

#include "NamespaceFooter.H"
#endif
