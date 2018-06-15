#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIRICHLETPOISSONEBBC_H_
#define _DIRICHLETPOISSONEBBC_H_

#include "RefCountedPtr.H"

#include "BaseEBBC.H"
#include "BaseBCValue.H"
#include "LayoutData.H"
#include "IntVectSet.H"
#include "EBStencil.H"
#include "NamespaceHeader.H"

///
/**
 */
class DirichletPoissonEBBC: public BaseEBBC
{
public:

  ///
  DirichletPoissonEBBC();

  DirichletPoissonEBBC(const ProblemDomain& a_domain,
                       const EBISLayout&    a_layout,
                       const RealVect&      a_dx,
                       const IntVect*       a_ghostCellsPhi=0,
                       const IntVect*       a_ghostCellsRhs=0);

  void construct(const ProblemDomain& a_domain,
                 const EBISLayout&    a_layout,
                 const RealVect&      a_dx,
                 const IntVect*       a_ghostCellsPhi=0,
                 const IntVect*       a_ghostCellsRhs=0);
  ///
  /**
   */
  virtual ~DirichletPoissonEBBC();

  ///
  /**
   */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor);

  ///
  /**
   */
  virtual void setOrder(int a_order);

  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    return &m_fluxStencil;
  }

  LayoutData<BaseIVFAB<Real> >& getFluxWeight()
  {
    return m_fluxWeight;
  }
  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_func);

  ///
  /**
   */
  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time);
  ///
  /**  deprecated interface
   */
  virtual void getEBFlux(Real&                         a_flux,
                         const VolIndex&               a_vof,
                         const LevelData<EBCellFAB>&   a_phi,
                         const LayoutData<IntVectSet>& a_cfivs,
                         const DataIndex&              a_dit,
                         const RealVect&               a_probLo,
                         const RealVect&               a_dx,
                         const bool&                   a_useHomogeneous,
                         const Real&                   a_time,
                         const pair<int,Real>*         a_cacheHint=0 );


  bool getSecondOrderStencil(VoFStencil&          a_stencil,
                             Real&                a_weight,
                             Vector<VoFStencil>&  a_pointStencils,
                             Vector<Real>&        a_distanceAlongLine,
                             const VolIndex&      a_vof,
                             const EBISBox&       a_ebisBox,
                             const RealVect&      a_dx,
                             const IntVectSet&    a_cfivs);

  void applyEBFluxPoint(const VolIndex&               a_vof,
                        EBCellFAB&                    a_lphi,
                        const EBCellFAB&              a_phi,
                        VoFIterator&                  a_vofit,
                        const LayoutData<IntVectSet>& a_cfivs,
                        const DataIndex&              a_dit,
                        const RealVect&               a_probLo,
                        const RealVect&               a_dx,
                        const Real&                   a_factor,
                        const bool&                   a_useHomogeneous,
                        const Real&                   a_time)          ;


  //hack to let component of velocity get set
  static int s_velComp;
  static int s_leastSquaresRad;
  static bool s_areaFracWeighted;
  static bool s_useQuadrantBasedStencil;

protected:
  virtual void getFirstOrderStencil(VoFStencil&     a_stencil,
                                    Real&           a_weight,
                                    const VolIndex& a_vof,
                                    const EBISBox&  a_ebisBox,
                                    const RealVect& a_dx);

  virtual void getSecondOrderStencil(VoFStencil&     a_stencil,
                                     Real&           a_weight,
                                     const VolIndex& a_vof,
                                     const EBISBox&  a_ebisBox,
                                     const RealVect& a_dx,
                                     const IntVectSet& a_cfivs);

private:

public:
  bool m_isDefined;

  bool m_onlyHomogeneous;
  bool m_isFunction;

  int m_order;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_func;

  ProblemDomain m_domain;
  EBISLayout    m_layout;

  RealVect m_dx;


  //stencils for operator evaluation
  IntVect                   m_ghostCellsPhi;
  IntVect                   m_ghostCellsRHS;

  //the stuff below is to keep EBAMRPoissonOp working in the transition
  virtual void define(const LayoutData<IntVectSet>& a_cfivs);
  // m_fluxWeight is the weight of the Dirichlet BC value used when
  // calculating the EB flux as a linear combination of the stencil
  // cell values and the BC value.
  LayoutData<BaseIVFAB<Real> >       m_fluxWeight;
  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil;

};

///
/**
 */
class DirichletPoissonEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  /**
   */
  DirichletPoissonEBBCFactory();

  ///
  /**
   */
  virtual ~DirichletPoissonEBBCFactory();

  ///
  /**
   */
  virtual void setOrder(int a_order);

  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_func);

  ///
  /**
   */
  virtual DirichletPoissonEBBC* create(const ProblemDomain& a_domain,
                                       const EBISLayout&    a_layout,
                                       const RealVect&      a_dx,
                                       const IntVect*       a_ghostCellsPhi=0,
                                       const IntVect*       a_ghostCellsRhs=0);

protected:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  int m_order;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_func;

};

#include "NamespaceFooter.H"
#endif
