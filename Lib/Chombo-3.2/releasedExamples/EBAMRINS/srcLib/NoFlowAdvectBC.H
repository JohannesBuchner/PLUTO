#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NOFLOWADVECTBC_H_
#define _NOFLOWADVECTBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

///sets flux to zero in all cases
class NoFlowAdvectBC : public EBPhysIBC
{
public:
  virtual ~NoFlowAdvectBC();

  ///
  NoFlowAdvectBC(int a_nComp);

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

  /// Initialize---no op here
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const;

  /// Set boundary slopes --no op here
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir)
  {;}


protected:
  bool           m_isDefined;
  RealVect       m_dx;
  ProblemDomain  m_domain;
  int            m_velComp;

private:
  NoFlowAdvectBC()
  {;}
};


///
/**
 */
class NoFlowAdvectBCFactory: public EBPhysIBCFactory
{
public:

  ///
  NoFlowAdvectBCFactory(int a_velComp);

  ///
  ~NoFlowAdvectBCFactory();


  ///
  EBPhysIBC* create() const;

protected:

  int m_velComp;

private:

  NoFlowAdvectBCFactory()
  {;}
};

#include "NamespaceFooter.H"
#endif
