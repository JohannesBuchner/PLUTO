#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPOISSONOP_H__
#define _EBPOISSONOP_H__

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include <map>
#include "RefCountedPtr.H"

#include "AMRMultiGrid.H"

#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"

#include "EBStencil.H"

#include "EBLevelDataOps.H"
#include "BaseEBBC.H"
#include "BaseDomainBC.H"
#include "CFIVS.H"
#include "EBFluxRegister.H"
#include "EBMGAverage.H"
#include "EBMGInterp.H"
#include "EBCoarsen.H"
#include "PolyGeom.H"
#include "EBQuadCFInterp.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"


// Weight for weighted Jacobi.  This kind of ought to be a member
// variable that could be set to a different value at run time.  For
// now, this gives the best convergence.
#define EBPOISSONOP_JACOBI_OMEGA 0.67


#if CH_SPACEDIM==2
#define EBPO_NUMSTEN 4
#elif CH_SPACEDIM==3
#define EBPO_NUMSTEN 8
#else
void THIS_IS_AN_ERROR_MESSAGE(void)
{
  THIS_WILL_ONLY_COMPILE_WHEN_CH_SPACEDIM_IS_2_OR_3;
}
#endif
///
/**
   Operator to solve (alpha + beta lapl)phi = rhs.   This follows the AMRLevelOp interface.
*/
class EBPoissonOp: public MGLevelOp<LevelData<EBCellFAB> >
{
public:
  ///
  ~EBPoissonOp();

  ///
  EBPoissonOp();

  ///
  /**
     If you are approaching this operator from this interface, consider backing away and using
     EBPoissonOpFactory to generate these objects.   Really.
     a_eblg,        : grid at this  level \\
     a_eblgCoarMG,  : grid at intermediate multigrid level \\
     a_domainBC,    : domain boundary conditions at this level   \\
     a_ebBC:          eb boundary conditions at this level   \\
     a_dx:            grid spacing at this level \\
     a_origin:        offset to lowest corner of the domain \\
     a_hasCoarserMG:    true  if there is a coarser MultiGrid level. \\
     a_preCondIters:  number of iterations to do for pre-conditioning \\
     a_relaxType:     0 means point Jacobi, 1 is Gauss-Seidel. \\
     a_orderEB:       0 to not do flux interpolation at cut faces. \\
     a_alpha:         coefficent of identity \\
     a_beta:          coefficient of laplacian.\\
     a_ghostCellsPhi:  Number of ghost cells in phi, correction (typically one)\\
     a_ghostCellsRhs:  Number of ghost cells in RHS, residual, lphi (typically zero)\\
     Ghost cell arguments are there for caching reasons.  Once you set them, an error is thrown if
     you send in data that does not match.
  */
  EBPoissonOp(   const EBLevelGrid &                  a_eblg,
                 const EBLevelGrid &                  a_eblgCoarMG,
                 const RefCountedPtr<BaseDomainBC>&   a_domainBC,
                 const RefCountedPtr<BaseEBBC>&       a_ebBC,
                 const RealVect&                      a_dx,
                 const RealVect&                      a_origin,
                 const bool&                          a_hasMGObjects,
                 const int&                           a_numPreCondIters,
                 const int&                           a_relaxType,
                 const int&                           a_orderEB,
                 const Real&                          a_alpha,
                 const Real&                          a_beta,
                 const IntVect&                       a_ghostCellsPhi,
                 const IntVect&                       a_ghostCellsRHS);

  //MGOp operations.  no finer or coarser

  ///
  /**
   */
  virtual void residual(LevelData<EBCellFAB>&       a_residual,
                        const LevelData<EBCellFAB>& a_phi,
                        const LevelData<EBCellFAB>& a_rhs,
                        bool                        a_homogeneousPhysBC=false);

  ///
  /**
   */
  virtual void preCond(LevelData<EBCellFAB>&       a_opPhi,
                       const LevelData<EBCellFAB>& a_phi);


  virtual void
  GSColorAllRegular(LevelData<EBCellFAB>&             a_phi,
                    const LevelData<EBCellFAB>&       a_rhs,
                    const IntVect&                    a_color,
                    const Real&                       a_weight,
                    const bool&                       a_homogeneousPhysBC);

  virtual void
  GSColorAllIrregular(LevelData<EBCellFAB>&             a_phi,
                      const LevelData<EBCellFAB>&       a_rhs,
                      const IntVect&                    a_color,
                      const bool&                       a_homogeneousPhysBC,
                      int icolor);


  ///
  /**
     this is the linearop function.
  */
  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       bool                              a_homogeneousPhysBC,
                       DataIterator& a_dit,
                       bool do_exchange = true);

  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       bool                              a_homogeneousPhysBC)
  {
    DataIterator dit = a_opPhi.dataIterator();
    applyOp(a_opPhi, a_phi, a_homogeneousPhysBC, dit, true);
  }

  ///
  /**
   */
  virtual void create(LevelData<EBCellFAB>&       a_lhs,
                      const LevelData<EBCellFAB>& a_rhs);

  ///
  /**
   */
  virtual void assign(LevelData<EBCellFAB>&       a_lhs,
                      const LevelData<EBCellFAB>& a_rhs);

  ///
  /**
   */
  virtual Real dotProduct(const LevelData<EBCellFAB>& a_1,
                          const LevelData<EBCellFAB>& a_2);

  ///
  /**
   */
  virtual void incr(LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    Real                        a_scale);

  ///
  /**
   */
  virtual void axby(LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    const LevelData<EBCellFAB>& a_y,
                    Real                        a_a,
                    Real                        a_b);

  ///
  /**
   */
  virtual void scale(LevelData<EBCellFAB>& a_lhs,
                     const Real&           a_scale);

  ///
  /**
   */
  virtual Real norm(const LevelData<EBCellFAB>& a_rhs,
                    int                         a_ord);

  ///
  /**
   */
  virtual void setToZero(LevelData<EBCellFAB>& a_lhs);

  ///
  /**
   */
  virtual void setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value);

  ///
  /**
   */
  virtual void createCoarser(LevelData<EBCellFAB>&       a_coarse,
                             const LevelData<EBCellFAB>& a_fine,
                             bool                        a_ghosted);

  ///
  /**
   */
  virtual void relax(LevelData<EBCellFAB>&       a_e,
                     const LevelData<EBCellFAB>& a_residual,
                     int                         a_iterations);

  ///
  /**
     Calculate restricted residual:
     a_resCoarse[2h] = I[h->2h] (a_rhsFine[h] - L[h](a_phiFine[h]))
  */
  virtual void restrictResidual(LevelData<EBCellFAB>&       a_resCoarse,
                                LevelData<EBCellFAB>&       a_phiFine,
                                const LevelData<EBCellFAB>& a_rhsFine);

  ///
  /**
     Correct the fine solution based on coarse correction:
     a_phiThisLevel += I[2h->h] (a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<EBCellFAB>&       a_phiThisLevel,
                                const LevelData<EBCellFAB>& a_correctCoarse);


  static bool nextColor(IntVect& color,
                        const IntVect& limit);

  static int getIndex(const IntVect& a_color);

  static IntVect getColor(const int& a_icolor);

  void levelMulticolorGS(LevelData<EBCellFAB>&       a_phi,
                         const LevelData<EBCellFAB>& a_rhs);

  void colorGS(LevelData<EBCellFAB>&       a_phi,
               const LevelData<EBCellFAB>& a_rhs,
               const IntVect& color, int icolor);



protected:

  void defineStencils();

  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;
  Vector<IntVect> m_colors;


  EBLevelGrid                     m_eblg;
  EBLevelGrid                     m_eblgCoarMG;

  RefCountedPtr<BaseDomainBC>     m_domainBC;
  RefCountedPtr<BaseEBBC>         m_ebBC;

  RealVect                        m_dx;

  RealVect                        m_invDx;
  RealVect                        m_invDx2;
  Real                            m_dxScale;
  Real                            m_alpha;
  Real                            m_beta;
  Real                            m_time;
  RealVect                        m_origin;
  int                             m_orderEB;
  int                             m_numPreCondIters;
  int                             m_relaxType;

  Copier                          m_exchangeCopier;

  //stencils for operator evaluation
  LayoutData<RefCountedPtr<EBStencil> >  m_opEBStencil;
  //stencils for operator evaluation on gauss-seidel colors
  LayoutData<RefCountedPtr<EBStencil> >  m_colorEBStencil[EBPO_NUMSTEN];
  //stencils for lambda*rhs addition to opphi
  LayoutData<RefCountedPtr<EBStencil> >  m_rhsColorEBStencil[EBPO_NUMSTEN];
  //weights that get multiplied by alpha
  LayoutData<BaseIVFAB<Real> >       m_alphaDiagWeight;
  //weights that get multiplied by beta
  LayoutData<BaseIVFAB<Real> >       m_betaDiagWeight;

  //cache the irreg vofiterator
  LayoutData<VoFIterator > m_vofItIrreg;
  LayoutData<VoFIterator > m_vofItIrregColor[EBPO_NUMSTEN];

  LayoutData<VoFIterator > m_vofItIrregDomLo[SpaceDim];
  LayoutData<VoFIterator > m_vofItIrregDomHi[SpaceDim];

  LayoutData<VoFIterator > m_vofItIrregColorDomLo[EBPO_NUMSTEN][SpaceDim];
  LayoutData<VoFIterator > m_vofItIrregColorDomHi[EBPO_NUMSTEN][SpaceDim];

  bool                        m_hasMGObjects;
  //stuff below is only defined if m_hasMGObjects==true
  EBMGAverage                 m_ebAverageMG;
  EBMGInterp                  m_ebInterpMG;
  DisjointBoxLayout           m_dblCoarMG;
  EBISLayout                  m_ebislCoarMG;
  ProblemDomain               m_domainCoarMG;

  /// Get the coefficients (for all cells) to multiply the residual by
  /// before adding it to the previous iterate for Jacobi iteration.
  void getJacobiRelaxCoeff(LevelData<EBCellFAB>& a_relaxCoeff);

private:


  //internally useful but not for general consumption
  //lots of hidden assumptions and all that

  void getOpVoFStencil(VoFStencil&     a_stencil,
                       const EBISBox&  a_ebisbox,
                       const VolIndex& a_vof);

  void getOpVoFStencil(VoFStencil&             a_stencil,
                       const int&              a_dir,
                       const Vector<VolIndex>& a_allMonotoneVoFs,
                       const EBISBox&          a_ebisbox,
                       const VolIndex&         a_vof,
                       const bool&             a_lowOrder);


  void getOpFaceStencil(VoFStencil&             a_stencil,
                        const Vector<VolIndex>& a_allMonotoneVofs,
                        const EBISBox&          a_ebisbox,
                        const VolIndex&         a_vof,
                        int                     a_dir,
                        const Side::LoHiSide&   a_side,
                        const FaceIndex&        a_face,
                        const bool&             a_lowOrder);


  void levelJacobi(LevelData<EBCellFAB>&       a_phi,
                   const LevelData<EBCellFAB>& a_rhs);

  ///
  /** applyOp() on the regular cells.  Only used in line relaxation stuff
      opPhi comes in holding alpha*phi.  this adds on beta*lapl phi*/
  void applyOpRegular( int idir,
                       Box * a_loBox,
                       Box * a_hiBox,
                       int * a_hasLo,
                       int * a_hasHi,
                       Box & a_curOpPhiBox,
                       Box & a_curPhiBox,
                       int a_nComps,
                       BaseFab<Real> & a_curOpPhiFAB,
                       const BaseFab<Real> & a_curPhiFAB,
                       bool a_homogeneousPhysBC,
                       const DataIndex& a_dit,
                       const Real& a_beta);

  ///
  /** applyOp() on the regular cells for all directions
      opPhi comes in holding alpha*phi.  this adds on beta*lapl phi*/
  void applyOpRegularAllDirs(Box * a_loBox,
                             Box * a_hiBox,
                             int * a_hasLo,
                             int * a_hasHi,
                             Box & a_curOpPhiBox,
                             Box & a_curPhiBox,
                             int a_nComps,
                             BaseFab<Real> & a_curOpPhiFAB,
                             const BaseFab<Real> & a_curPhiFAB,
                             bool a_homogeneousPhysBC,
                             const DataIndex& a_dit,
                             const Real& a_beta);

  void applyDomainFlux(Box * a_loBox,
                       Box * a_hiBox,
                       int * a_hasLo,
                       int * a_hasHi,
                       Box & a_curPhiBox,
                       int a_nComps,
                       BaseFab<Real> & a_phiFAB,
                       bool a_homogeneousPhysBC,
                       const DataIndex& a_dit,
                       const Real& a_beta);

  ///
  /**
     This function computes: a_lhs = (1/diagonal)*a_rhs
     Use this function to initialize the preconditioner
  */
  void getInvDiagRHS(LevelData<EBCellFAB>&       a_lhs,
                     const LevelData<EBCellFAB>& a_rhs);

private:
  //copy constructor and operator= disallowed for all the usual reasons
  EBPoissonOp(const EBPoissonOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBPoissonOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};


#include "NamespaceFooter.H"
#endif
