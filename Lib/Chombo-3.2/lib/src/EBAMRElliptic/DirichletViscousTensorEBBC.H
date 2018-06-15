#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DIRICHLETVISCOUSTENSOREBBC_H_
#define _DIRICHLETVISCOUSTENSOREBBC_H_

#include "RefCountedPtr.H"

#include "BaseEBBC.H"
#include "BaseBCFuncEval.H"
#include "LayoutData.H"
#include "IntVectSet.H"
#include "EBStencil.H"
#include "NamespaceHeader.H"
class VoFCache;
class VoFStencil;

///
/**
 */
class DirichletViscousTensorEBBC: public ViscousBaseEBBC
{
public:

  ///
  DirichletViscousTensorEBBC();

  ///
  /**
     If the optional a_ghostCells(Phi,Rhs) args are NULL, we won't use the
     VoF cache.
  */
  DirichletViscousTensorEBBC( const ProblemDomain& a_domain,
                              const EBISLayout&    a_layout,
                              const RealVect&      a_dx,
                              const IntVect*       a_ghostCellsPhi,
                              const IntVect*       a_ghostCellsRhs);

  ///
  /**
   */
  virtual ~DirichletViscousTensorEBBC();

  ///
  /**
   */
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor);

  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    return &m_fluxStencil[ivar];
  }

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
  /**
     public for testing.   gets the gradient of the solution at centroid
     of the embedded boundary given the boundary condition and the solution.
     In both cases, it uses the johansen stencil for the normal gradient.   For the tangential
     gradient, in the function case, it calls the derivative function (in the value case, it sets
     tangential gradients to zero).  Since the derivative function operates in cartesian space
     there is some complication (set all cartesian gradients.   rotate to NT space.   overwrite
     normal gradient.   rotate back to cartesian space).
  **/
  void getGradient(Real             a_grad[SpaceDim][SpaceDim],
                   const VolIndex&  a_vof,
                   const EBCellFAB& a_phi,
                   const EBISBox&   a_ebisBox,
                   const DataIndex& a_dit,
                   const Real&      a_dx,
                   bool a_homogeneous);

  //public to test the stencil.
  void getGradientStenValue(Real             a_grad[SpaceDim][SpaceDim],
                            const VolIndex&  a_vof,
                            const EBCellFAB& a_phi,
                            const EBISBox&   a_ebisBox,
                            const DataIndex& a_dit,
                            const Real&      a_dx,
                            bool a_homogeneous);
protected:
  //get the gradient associated with the inhomogeneous bit.
  void getGradInhomOnly(Real             a_grad[SpaceDim][SpaceDim],
                        const Real&      a_weight,
                        const VolIndex&  a_vof,
                        const EBISBox&   a_ebisBox,
                        const Real&      a_dx);
  //get the stencil for the normal gradient
  void  getNormalStencil(VoFStencil&          a_stencil,
                         Real&                a_weight,
                         const VolIndex&      a_vof,
                         const EBISBox&       a_ebisBox,
                         const RealVect&      a_dx,
                         const IntVectSet&    a_cfivs);

  //get the homogeneous flux stencil
  void getFluxStencil(VoFStencil        a_stencil[SpaceDim],
                      Real              a_weight[ SpaceDim],
                      const DataIndex&  a_dit,
                      const VolIndex&   a_vof,
                      const EBISBox&    a_ebisBox,
                      const RealVect&   a_dx,
                      const IntVectSet& a_cfivs);

  //another testing routine
  void
  getGradientFunction(Real             a_grad[SpaceDim][SpaceDim],
                      const VolIndex&  a_vof,
                      const EBCellFAB& a_phi,
                      const EBISBox&   a_ebisBox,
                      const DataIndex& a_dit,
                      const Real&      a_dx,
                      bool a_homogeneous);

  //return true if you need to drop order
  bool getSecondOrderStencil(VoFStencil&          a_stencil,
                             Real&                a_weight,
                             Vector<VoFStencil>&  a_pointStencil,
                             Vector<Real>&        a_distanceAlongLine,
                             const VolIndex&      a_vof,
                             const EBISBox&       a_ebisBox,
                             const RealVect&      a_dx,
                             const IntVectSet&    a_cfivs);

  void getFirstOrderStencil(VoFStencil&       a_stencil,
                            Real&             a_weight,
                            const VolIndex&   a_vof,
                            const EBISBox&    a_ebisBox,
                            const RealVect&   a_dx,
                            const IntVectSet& a_cfivs);


  void getJacobianAndInverse(Real a_Jacobian[SpaceDim][SpaceDim],
                             Real a_Jinverse[SpaceDim][SpaceDim],
                             RealVect& a_normal,
                             RealVect a_tangents[SpaceDim-1]);


  virtual RealVect getInhomogeneousContribution(const VolIndex&  a_vof,
                                                const EBCellFAB& a_phi,
                                                const EBISBox&   a_ebisBox,
                                                const DataIndex& a_dit,
                                                const Real&      a_dx);

  void getFlux(Real             a_flux[SpaceDim][SpaceDim],
               const VolIndex&  a_vof,
               const EBCellFAB& a_phi,
               const EBISBox&   a_ebisBox,
               const DataIndex& a_dit,
               const Real&      a_dx,
               bool a_homogeneous);


  void getCartesianGradientStencil(VoFStencil           a_gradStencils[SpaceDim][SpaceDim],
                                   VoFStencil &         a_normalStencil,
                                   const DataIndex&     a_dit,
                                   const VolIndex&      a_vof,
                                   const EBISBox&       a_ebisBox,
                                   const RealVect&      a_dx);


  static int s_leastSquaresRad;
protected:

  bool m_isDefined;

  ProblemDomain m_domain;
  EBISLayout    m_ebisl;

  RealVect m_dx;


  //stencils for operator evaluation
  LayoutData<BaseIVFAB<VoFStencil> > m_fluxStencil[CH_SPACEDIM];
  LayoutData<BaseIVFAB<Real> >        m_fluxWeight[CH_SPACEDIM];
  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;


};

///
/**
 */
class DirichletViscousTensorEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  /**
   */
  DirichletViscousTensorEBBCFactory();

  ///
  /**
   */
  virtual ~DirichletViscousTensorEBBCFactory();


  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCFuncEval> a_func);

  ///
  /**
   */
  virtual DirichletViscousTensorEBBC* create(  const ProblemDomain& a_domain,
                                               const EBISLayout&    a_layout,
                                               const RealVect&      a_dx,
                                               const IntVect*       a_ghostCellsPhi,
                                               const IntVect*       a_ghostCellsRhs);

protected:

private:
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCFuncEval>   m_func;
  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;

};



#include "NamespaceFooter.H"
#endif
