#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SMOOTHBUMPIBC_H_
#define _SMOOTHBUMPIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"
#include "SmoothBumpF_F.H"
#include "EBSolidF_F.H"

#include "NamespaceHeader.H"

class SmoothBumpIBC : public EBPhysIBC
{
public:
  virtual ~SmoothBumpIBC();

  SmoothBumpIBC(const Real&     a_gamma,
                const RealVect& a_center,
                const Real&     a_mach,
                const Real&     a_rnot);

  ///
  void define(const ProblemDomain&  a_domain,
              const RealVect&       a_dx);

  ///  For every box in this level, this function is called.
  void fluxBC(EBFluxFAB&            a_flux,
              const EBCellFAB&      a_Wcenter,
              const EBCellFAB&      a_Wextrap,
              const Side::LoHiSide& a_sd,
              const Real&           a_time,
              const EBISBox&        a_ebisBox,
              const DataIndex&      a_dit,
              const Box&            a_box,
              const Box&            a_faceBox,
              const int&            a_dir);

  /// Initialize
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
  */
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir);

protected:
  bool m_isFortranCommonSet;
  bool m_isDefined;
  Real     m_gamma;
  RealVect m_center;
  Real     m_mach;
  Real     m_rnot;

  Real     m_dx;
  ProblemDomain  m_domain;

private:
  SmoothBumpIBC()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const SmoothBumpIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
  SmoothBumpIBC(const SmoothBumpIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }


};

class SmoothBumpIBCFactory: public EBPhysIBCFactory
{
public:

  ///
  virtual ~SmoothBumpIBCFactory()
  {
  }

  SmoothBumpIBCFactory(const Real&     a_gamma,
                       const RealVect& a_center,
                       const Real&     a_mach,
                       const Real&     a_rnot):EBPhysIBCFactory()
  {
    m_gamma  = a_gamma;
    m_center = a_center;
    m_mach   = a_mach;
    m_rnot   = a_rnot;
  }

  ///
  EBPhysIBC* create() const
  {
    SmoothBumpIBC* retval =
      new SmoothBumpIBC(m_gamma,
                        m_center,
                        m_mach,
                        m_rnot);
    return static_cast<EBPhysIBC*>(retval);
  }

protected:
  Real     m_gamma;
  Real     m_mach;
  Real     m_rnot;
  RealVect m_center;

private:
  SmoothBumpIBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const SmoothBumpIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  SmoothBumpIBCFactory(const SmoothBumpIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
