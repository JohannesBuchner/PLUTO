#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPlanarShockTemperatureBC_H_
#define _EBPlanarShockTemperatureBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"

#include "NamespaceHeader.H"

///
/**
 */
class EBPlanarShockTemperatureBC: public ConductivityBaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~EBPlanarShockTemperatureBC()
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
  EBPlanarShockTemperatureBC(const int&  a_shockNorm, 
                        const bool& a_shockBackward,   
                        const bool& a_slipWall)
  {
    m_shockNorm       = a_shockNorm;
    m_shockBackward = a_shockBackward;
    m_slipWall      = a_slipWall;
  }

private:
  void
  whereAMI(bool& a_atInflow, 
           bool& a_atOutflow,
           const int&            a_idir, 
           const Side::LoHiSide& a_side);

  ///weak construction forbidden to keep things simple
  EBPlanarShockTemperatureBC(){;};
  int  m_shockNorm;
  bool m_shockBackward;
  bool m_slipWall;

};

class EBPlanarShockTemperatureBCFactory: public BaseDomainBCFactory
{
public:


  ///
  /**
   */
  EBPlanarShockTemperatureBCFactory(const int&  a_shockNorm, 
                               const bool& a_shockBackward, 
                               const bool& a_slipWall)
  {
    m_shockNorm       = a_shockNorm;
    m_shockBackward = a_shockBackward;
    m_slipWall      = a_slipWall;
  }


  ///
  /**
   */
  virtual ~EBPlanarShockTemperatureBCFactory()
  {
  }

  ///
  /**
   */
  virtual EBPlanarShockTemperatureBC* create(const ProblemDomain& a_domain,
                                        const EBISLayout&    a_layout,
                                        const RealVect&      a_dx)
  {
    EBPlanarShockTemperatureBC* newBC = new EBPlanarShockTemperatureBC(m_shockNorm, m_shockBackward,  m_slipWall);

    return newBC;
  }


private:

  //insert incoherent rant about weak construction here
  //(it really does make things simpler)
  EBPlanarShockTemperatureBCFactory(){;}


  int  m_shockNorm;
  bool m_shockBackward;
  bool m_slipWall;

};

#include "NamespaceFooter.H"
#endif
