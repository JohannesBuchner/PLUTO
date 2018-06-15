#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef  _EB_EXPLOSIONIBC_Hummena_Hummena___
#define  _EB_EXPLOSIONIBC_Hummena_Hummena___

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "EBPhysIBC.H"
#include "EBExplosionF_F.H"
#include "EBSolidF_F.H"

#include "NamespaceHeader.H"

class EBExplosionIBC : public EBPhysIBC
{
public:
  virtual ~EBExplosionIBC();

  EBExplosionIBC(const Real&     a_gamma,
                 const Real&     a_size,
                 const Real&     a_p0,
                 const Real&     a_r0,
                 const Real&     a_p1,
                 const Real&     a_r1,
                 const RealVect& a_center,
                 const int&      a_doOneDOnly,
                 const RealVect& a_oneDNormal,
                 const RealVect& a_oneDOrigin
                 );

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
  RealVect       m_oneDNormal;
  RealVect       m_oneDOrigin;
  bool           m_isFortranCommonSet;
  bool           m_isDefined;
  Real           m_gamma;
  RealVect       m_center;
  Real           m_size;
  Real           m_p0;
  Real           m_r0;
  Real           m_p1;
  Real           m_r1;
  int            m_doOneDOnly;
  Real           m_dx;
  ProblemDomain  m_domain;

private:
  EBExplosionIBC()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const EBExplosionIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBExplosionIBC(const EBExplosionIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
