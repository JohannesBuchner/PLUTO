#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef  _EB_NOFLOWIBC_H_
#define  _EB_NOFLOWIBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

//! \class EBNoFlowIBC
//! This class, for use with EBAMRCNS, ensures that there is no flow across an embedded boundary.
//! Use this only for toy problems where all the action is far from the boundary.
class EBNoFlowIBC: public EBPhysIBC
{
  public:

  //! Create a no-flow boundary condition for an ideal gas.
  //! \param a_gamma The polytropic index of the gas.
  explicit EBNoFlowIBC(Real a_gamma);

  //! Destructor.
  virtual ~EBNoFlowIBC();

  //----------------------------------------------------------
  // The following methods are overridden from the base class.
  //----------------------------------------------------------
  
  void define(const ProblemDomain&  a_domain,
              const RealVect&       a_dx);

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

  
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const;

  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir);

  protected:

  // Domain and cell size for this level.
  RealVect m_dx;
  ProblemDomain m_domain;

  // Set to true if the domain has been defined, false otherwise.
  bool m_isDefined;

  private:

  // Disallowed operators.
  EBNoFlowIBC();
  EBNoFlowIBC(const EBNoFlowIBC& a_rhs);
  EBNoFlowIBC& operator=(const EBNoFlowIBC& a_rhs);
};

//! \class EBNoFlowIBCFactory
//! This factory creates EBNoFlowIBC classes.
class EBNoFlowIBCFactory: public EBPhysIBCFactory
{
  public:

  //! Create a factory that generates no-flow IBCs for a gas with the given 
  //! polytropic index.
  //! \param a_gamma The given polytropic index.
  explicit EBNoFlowIBCFactory(Real a_gamma);

  //! Destructor.
  ~EBNoFlowIBCFactory();

  EBPhysIBC* create() const;

  private:

  // Polytropic index.
  Real m_gamma;

  // Forbidden operations.
  EBNoFlowIBCFactory();
  EBNoFlowIBCFactory(const EBNoFlowIBCFactory&);
  EBNoFlowIBCFactory& operator=(const EBNoFlowIBCFactory&);
};

#include "NamespaceFooter.H"
#endif
