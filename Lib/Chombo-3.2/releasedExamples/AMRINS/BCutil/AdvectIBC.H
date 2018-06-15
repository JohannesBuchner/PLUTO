#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ADVECTIBC_H_
#define _ADVECTIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "PhysIBC.H"

#include "UsingNamespace.H"

/// IBC for simple advection
/**
   Parameters:

 */
class AdvectIBC:public PhysIBC
{
public:
  /// Null constructor
  /**
   */
  AdvectIBC();

  /// Destructor
  /**
   */
  ~AdvectIBC()
  {
  }

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysIBC object with m_isDefined = false
     (i.e., its define() must be called before it is used) and
     m_isFortranCommonSet set to value of m_isFortranCommonset in the
     current (factory) object.
   */
  PhysIBC *new_physIBC();

  /// set default boundary values
  void setDefaultValues();

  /// Set boundary fluxes
  /**
   */
  void primBC(FArrayBox&            a_WGdnv,
              const FArrayBox&      a_Wextrap,
              const FArrayBox&      a_W,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  void setBdrySlopes(FArrayBox&       a_dW,
                     const FArrayBox& a_W,
                     const int&       a_dir,
                     const Real&      a_time);

  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_U);

  /// set velocity
  void advectionVel(const RealVect& a_advVel);

  ///
  const RealVect& advectionVel() const;

  /// set problem type
  void probType(const int a_probtype);

  /// set boundary value (default is 0)
  void setBoundaryValue(Real a_bcVal, int a_dir,
                        Side::LoHiSide a_hiLo);

  /// access boundary value
  Real getBoundaryValue(int a_dir, Side::LoHiSide a_hiLo) const;

  /// set slope value (default is 0)
  void setSlopeValue(Real a_slopeVal, int a_dir,
                     Side::LoHiSide a_hiLo);

  /// access slope value
  Real getSlopeValue(int a_dir, Side::LoHiSide a_hiLo) const;

  void artViscBC(FArrayBox&       a_F,
                 const FArrayBox& a_U,
                 const FArrayBox& a_divVel,
                 const int&       a_dir,
                 const Real&      a_time);

  /// accessor
  int probType() const;

protected:

  /// advection velocity field
  RealVect m_velocity;

  Real m_bcVal[SpaceDim][2];
  Real m_slopeVal[SpaceDim][2];

  bool m_isBCvalSet;
  bool m_isSlopeValSet;

  /// initial condition type
  int m_probtype;

};

#endif
