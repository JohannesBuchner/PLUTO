#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPHYSIBC_H_
#define _EBPHYSIBC_H_

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBFluxFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

///
/**
   Pure virtual base class through which a user specifies the
   initial and boundary conditions for LevelGodunov.
 */
class EBPhysIBC
{
public:
  ///
  EBPhysIBC();

  ///
  virtual ~EBPhysIBC();

  ///
  virtual void define(const ProblemDomain&  a_domain,
                      const RealVect&       a_dx) = 0;

  ///  For every box in this level, this virtual function is called.
  virtual void fluxBC(EBFluxFAB&            a_flux,
                      const EBCellFAB&      a_Wcenter,
                      const EBCellFAB&      a_Wextrap,
                      const Side::LoHiSide& a_sd,
                      const Real&           a_time,
                      const EBISBox&        a_ebisBox,
                      const DataIndex&      a_dit,
                      const Box&            a_box,
                      const Box&            a_faceBox,
                      const int&            a_dir) = 0;

  /// Put initial conserved state into a_conState
  virtual void initialize(LevelData<EBCellFAB>& a_conState,
                          const EBISLayout& a_ebisl) const = 0;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
   */
  virtual void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                              const EBCellFAB& a_primState,
                              const EBISBox&   a_ebisBox,
                              const Box&       a_box,
                              const int&       a_dir) = 0;

};

#include "NamespaceFooter.H"
#endif
