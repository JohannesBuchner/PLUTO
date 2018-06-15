#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MODIANOIBC_H_
#define _MODIANOIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "EBPhysIBC.H"
#include "EBFluxFAB.H"

#include "NamespaceHeader.H"

///
/**
   IBC for doing the modiano problem.
 */
class ModianoIBC : public EBPhysIBC
{
public:

  virtual ~ModianoIBC();

  ///
  /**
     waveAmp is the amplitude of the disturbance in the density
     (alpha in Modiano and Colella) . \\
     waveWidth is the width of the wave (w in Modiano Colella) \\
     waveDir is the direction of  the wave.  The wave runs perpendicular to
     the direction of the walls. \\
     gamma is of course the ratio of specific heats for the gas.
   */
  ModianoIBC(const Real&     a_gamma,
             const Real&     a_waveAmp,
             const Real&     a_waveWidth,
             const RealVect& a_center,
             const RealVect& a_waveDir,
             bool  a_freeStreamOnly,
             bool a_negativeWave = false);

  ///
  void define(const ProblemDomain&  a_domain,
              const RealVect&       a_dx);

  //deprecated interface
  void define(const ProblemDomain&  a_domain,
              const Real&       a_dx)
  {
    RealVect dxVect = a_dx*RealVect::Unit;
    define(a_domain, dxVect);
  }

  void  setToExact(LevelData<EBCellFAB>& a_conState,
                   const EBISLayout& a_ebisl,
                   const Real&       a_time) const;

  void  setToExact(EBCellFAB&     a_conState,
                   const EBISBox& a_ebisl,
                   const Real&    a_time) const;


  void  setToExactConsAndPrim(LevelData<EBCellFAB>& a_conState,
                              const EBISLayout& a_ebisl,
                              const Real&       a_time) const;

  void  setToExactConsAndPrim(EBCellFAB&     a_conState,
                              const EBISBox& a_ebisl,
                              const Real&    a_time) const;

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

  ///
  void initialize(EBCellFAB& a_conState,
                  const EBISBox& a_ebisl) const;

  ///
  void
  computeExactFluxes(EBFluxFAB& cflux,
                     BaseIVFAB<Real>  coveredCfluxMinu[SpaceDim],
                     BaseIVFAB<Real>  coveredCfluxPlus[SpaceDim],
                     Vector<VolIndex> coveredFacesMinu[SpaceDim],
                     Vector<VolIndex> coveredFacesPlus[SpaceDim],
                     IntVectSet coveredSetsMinu[SpaceDim],
                     IntVectSet coveredSetsPlus[SpaceDim],
                     BaseIVFAB<Real>& boundaryPress,
                     IntVectSet&      irregIVS,
                     const EBISBox& ebisBox,
                     const Box& a_region,
                     const Box& domainBox,
                     const Real& time);

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

  bool     m_isFortranCommonSet;
  bool     m_isDefined;
  Real     m_gamma;
  Real     m_waveAmp;
  Real     m_waveWidth;
  RealVect m_waveDir;
  RealVect m_center;
  bool m_freeStreamOnly;

  Real m_dx;
  ProblemDomain  m_domain;

  int m_iNegWave;
  int m_idoFreeStreamOnly;
private:
  //forcing strong construction
  ModianoIBC()
  {
    MayDay::Error("invalid operator");
  }
  //disallowed for all the usual reasons
  void operator=(const ModianoIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
  ModianoIBC(const ModianoIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }


};

extern RealVect
mibcGetNormalVector(const RealVect& channelNormal);

#include "NamespaceFooter.H"
#endif
