#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtg 2007

#ifndef _RESISTIVITYOP_H_
#define _RESISTIVITYOP_H_

#include "AMRMultiGrid.H"
#include "REAL.H"
#include "Box.H"
#include "LevelDataOps.H"
#include "BCFunc.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "CFIVS.H"
#include "TensorCFInterp.H"
#include "CoarseAverage.H"
#include "AMRTGA.H"
#include "NamespaceHeader.H"

///
/**
   operator L is defined as
   L(B) = alpha B + beta*(divF)
   where
   F = eta((grad B - grad B^T) + I  (div B))
   alpha and beta are constants.  eta is a function of space.
   beta and eta are incorporated into the flux F when getFlux is called.
   This is always 3 variables.  SpaceDim refers to only the derivatives
   that are meaningful in this bizzare application.
*/
class ResistivityOp : public LevelTGAHelmOp<LevelData<FArrayBox>, FluxBox>
{
public:

  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta)
  {
    m_alpha = a_alpha;
    m_beta = a_beta;
    setLambda();
  }

  virtual void diagonalScale(LevelData<FArrayBox>& a_rhs, bool a_kappaWeighted)
  {
    diagonalScale(a_rhs);
  }

  virtual void diagonalScale(LevelData<FArrayBox>& a_rhs)
  {
    //does not apply here
  }

  virtual void divideByIdentityCoef(LevelData<FArrayBox>& a_rhs)
  {
    //does not apply here
  }

  ///
  /**
   */
  virtual ~ResistivityOp()
  {
  }

  //use at thy peril
  virtual void getFlux(FluxBox&             a_flux,
               const LevelData<FArrayBox>&  a_data,
               const Box&       a_grid,
               const DataIndex& a_dit,
               Real             a_scale);

  ///
  /**
   */
  ResistivityOp(const DisjointBoxLayout&                    a_grids,
                const DisjointBoxLayout&                    a_gridsFine,
                const DisjointBoxLayout&                    a_gridsCoar,
                const RefCountedPtr<LevelData<FluxBox> >&   a_eta,
                Real                                        a_alpha,
                Real                                        a_beta,
                int                                         a_refToFine,
                int                                         a_refToCoar,
                const ProblemDomain&                        a_domain,
                const Real&                                 a_dxLevel,
                const Real&                                 a_dxCoar,
                BCFunc                                      a_bc);

  virtual void residual(  LevelData<FArrayBox>& a_lhs,
                          const LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>& a_rhs,
                          bool a_homogeneous = false);

  virtual void preCond(   LevelData<FArrayBox>& a_correction,
                          const LevelData<FArrayBox>& a_residual);

  virtual void applyOp(   LevelData<FArrayBox>& a_lhs,
                          const LevelData<FArrayBox>& a_phi,
                          bool a_homogeneous = false);

  virtual void applyOpNoBoundary( LevelData<FArrayBox>& a_lhs,
                                  const LevelData<FArrayBox>& a_phi);

  virtual void create(    LevelData<FArrayBox>& a_lhs,
                          const LevelData<FArrayBox>& a_rhs);
  virtual void createCoarsened(    LevelData<FArrayBox>& a_lhs,
                                   const LevelData<FArrayBox>& a_rhs,
                                   const int& a_refRat);

  void reflux(const LevelData<FArrayBox>& a_phiFine,
              const LevelData<FArrayBox>& a_phi,
              LevelData<FArrayBox>& residual,
              AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  void getFlux(FArrayBox&       a_flux,
               const FArrayBox& a_data,
               const FArrayBox& a_gradData,
               const FArrayBox& a_etaFace,
               const Box&       a_facebox,
               int              a_dir,
               int ref = 1);

  /// utility function which computes operator after all bc's have been set
  void computeOperatorNoBCs(LevelData<FArrayBox>& a_lhs,
                            const LevelData<FArrayBox>& a_phi);

  virtual void assign(    LevelData<FArrayBox>& a_lhs,
                          const LevelData<FArrayBox>& a_rhs) ;
  virtual Real dotProduct(const LevelData<FArrayBox>& a_1,
                          const LevelData<FArrayBox>& a_2) ;
  virtual void incr( LevelData<FArrayBox>& a_lhs,
                     const LevelData<FArrayBox>& a_x,
                     Real a_scale) ;
  virtual void axby( LevelData<FArrayBox>& a_lhs,
                     const LevelData<FArrayBox>& a_x,
                     const LevelData<FArrayBox>& a_y,
                     Real a, Real b) ;

  virtual void scale(LevelData<FArrayBox>& a_lhs, const Real& a_scale) ;

  virtual Real norm(const LevelData<FArrayBox>& a_x, int a_ord);

  virtual void setToZero( LevelData<FArrayBox>& a_x);
  /*@}*/

  /**
     \name MGLevelOp functions */
  /*@{*/

  virtual void relax(LevelData<FArrayBox>& a_e,
                     const LevelData<FArrayBox>& a_residual,
                     int iterations);

  virtual void createCoarser(LevelData<FArrayBox>& a_coarse,
                             const LevelData<FArrayBox>& a_fine,
                             bool ghosted);
  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<FArrayBox>& a_resCoarse,
                                LevelData<FArrayBox>& a_phiFine,
                                const LevelData<FArrayBox>& a_rhsFine);

  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<FArrayBox>& a_phiThisLevel,
                                const LevelData<FArrayBox>& a_correctCoarse);

  /*@}*/

  /**
     \name AMRLevelOp functions */
  /*@{*/

  /** returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToCoarser()
  {
    return m_refToCoar;
  }

  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMRResidual(LevelData<FArrayBox>& a_residual,
                           const LevelData<FArrayBox>& a_phiFine,
                           const LevelData<FArrayBox>& a_phi,
                           const LevelData<FArrayBox>& a_phiCoarse,
                           const LevelData<FArrayBox>& a_rhs,
                           bool a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** residual assuming no more coarser AMR levels */

  virtual void AMRResidualNC(LevelData<FArrayBox>& a_residual,
                             const LevelData<FArrayBox>& a_phiFine,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_rhs,
                             bool a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMRResidualNF(LevelData<FArrayBox>& a_residual,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_phiCoarse,
                             const LevelData<FArrayBox>& a_rhs,
                             bool a_homogeneousPhysBC);

  /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
      it is assumed that a_resCoarse has already been filled in with the coarse
      version of AMRResidualNF and that this operation is free to overwrite
      in the overlap regions.
  */

  virtual void AMRRestrict(LevelData<FArrayBox>& a_resCoarse,
                           const LevelData<FArrayBox>& a_residual,
                           const LevelData<FArrayBox>& a_correction,
                           const LevelData<FArrayBox>& a_coarseCorrection, 
                           bool a_skip_res = false );

  /** a_correction += I[h->h](a_coarseCorrection) */
  virtual void AMRProlong(LevelData<FArrayBox>& a_correction,
                          const LevelData<FArrayBox>& a_coarseCorrection);

  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(LevelData<FArrayBox>& a_residual,
                                 const LevelData<FArrayBox>& a_correction,
                                 const LevelData<FArrayBox>& a_coarseCorrection);

  ///
  /**
     compute norm over all cells on coarse not covered by finer
  */
  virtual Real AMRNorm(const LevelData<FArrayBox>& a_coarseResid,
                       const LevelData<FArrayBox>& a_fineResid,
                       const int&                  a_refRat,
                       const int&                  a_ord);

  void homogeneousCFInterp(LevelData<FArrayBox>& a_phif);
 ///
  /**
     does homogeneous coarse/fine interpolation for phi
  */
  void
  homogeneousCFInterpPhi(LevelData<FArrayBox>& a_phif,
                         const DataIndex& a_datInd,
                         int a_idir,
                         Side::LoHiSide a_hiorlo);

  ///
  /** does homogeneous coarse/fine interpolation for tangential gradient
      (needs phi ghost cells to be filled in first, so should call
      homogeneousCFInterpPhi first)
  */
  void homogeneousCFInterpTanGrad(LevelData<FArrayBox>& a_tanGrad,
                                  const LevelData<FArrayBox>& a_phi,
                                  const DataIndex& a_datInd,
                                  int a_idir,
                                  Side::LoHiSide a_hiorlo);

  void interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                       const DataIndex& a_datInd,
                       const int a_idir,
                       const Side::LoHiSide a_hiorlo,
                       const IntVectSet& a_interpIVS);
  ///
  /**
     Apply the AMR operator, including coarse-fine matching
  */
  void AMROperator(      LevelData<FArrayBox>& a_LofPhi,
                   const LevelData<FArrayBox>& a_phiFine,
                   const LevelData<FArrayBox>& a_phi,
                   const LevelData<FArrayBox>& a_phiCoarse,
                   bool a_homogeneousDomBC,
                   AMRLevelOp<LevelData<FArrayBox> >*  a_finerOp);

  void AMROperatorNF(      LevelData<FArrayBox>& a_LofPhi,
                     const LevelData<FArrayBox>& a_phi,
                     const LevelData<FArrayBox>& a_phiCoarse,
                     bool a_homogeneousBC);

  virtual void AMROperatorNC(      LevelData<FArrayBox>& a_LofPhi,
                             const LevelData<FArrayBox>& a_phiFine,
                             const LevelData<FArrayBox>& a_phi,
                             bool a_homogeneousBC,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  void cellGrad(FArrayBox&             a_gradPhi,
                const  FArrayBox&      a_phi,
                const Box&             a_grid);

  void cfinterp(const LevelData<FArrayBox>& a_phiFine,
                const LevelData<FArrayBox>& a_phiCoarse);

  virtual void fillGrad(const LevelData<FArrayBox>& a_phiFine);

  void loHiCenterFace(Box&                 a_loBox,
                      int&                 a_hasLo,
                      Box&                 a_hiBox,
                      int&                 a_hasHi,
                      Box&                 a_centerBox,
                      const ProblemDomain& a_eblg,
                      const Box&           a_inBox,
                      const int&           a_dir);

  void getFaceDivAndGrad(FArrayBox&             a_faceDiv,
                         FArrayBox&             a_faceGrad,
                         const FArrayBox&       a_data,
                         const FArrayBox&       a_gradData,
                         const ProblemDomain&   a_domain,
                         const Box&             a_faceBox,
                         const int&             a_faceDir,
                         const Real             a_dx);

  ///take cell centered divergence of the inputs.
  /**
      not part of the operator evaluation.  this is to
      test whether the divergence of b converges at h^2
      as b does if b is analytically divergence-free.
   */
  void divergenceCC(LevelData<FArrayBox>&       a_div,
                    const LevelData<FArrayBox>& a_phi,
                    const LevelData<FArrayBox>* a_phiC);
protected:
  void setLambda();
  RefCountedPtr<LevelData<FluxBox> >         m_eta;
  DisjointBoxLayout                          m_grids;
  Real                                       m_alpha;
  Real                                       m_beta;
  int                                        m_refToCoar;
  int                                        m_refToFine;
  BCFunc                                     m_bc;

  //b is a scalar in 2d, vector in 3d
  static const int        s_nComp; //==3
  static const int        s_nGradComp; //==3*SpaceDim
  Real                    m_dx;
  Real                    m_dxCrse;
  ProblemDomain           m_domain;
  //relaxation coef
  LevelData<FArrayBox>    m_lambda;
  LevelData<FArrayBox>    m_grad;
  LevelDataOps<FArrayBox> m_levelOps;
  Copier                  m_exchangeCopier;
  TensorCFInterp          m_interpWithCoarser;

  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];
  // will need these to do tangential gradient computations
  LayoutData<TensorFineStencilSet> m_hiTanStencilSets[SpaceDim];
  LayoutData<TensorFineStencilSet> m_loTanStencilSets[SpaceDim];
  Vector<IntVect> m_colors;
private:
  ///weak construction is bad
  ResistivityOp()
  {
    MayDay::Error("invalid operator");
  }
  //copy constructor and operator= disallowed for all the usual reasons
  ResistivityOp(const ResistivityOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }
  void operator=(const ResistivityOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

///
/**
   Factory to create ResistivityOps
 */
class ResistivityOpFactory: public AMRLevelOpFactory<LevelData< FArrayBox> >
{
public:
  virtual ~ResistivityOpFactory()
  {
  }

  ResistivityOpFactory(const Vector<DisjointBoxLayout>&                     a_grids,
                       const Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_eta,
                       Real                                                 a_alpha,
                       Real                                                 a_beta,
                       const Vector<int>&                                   a_refRatios,
                       const ProblemDomain&                                 a_domainCoar,
                       const Real&                                          a_dxCoar,
                       BCFunc                                               a_bc);

  ///
  virtual MGLevelOp< LevelData<FArrayBox> >*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int depth,
          bool homoOnly = true);

  ///
  virtual AMRLevelOp< LevelData<FArrayBox> >*
  AMRnewOp(const ProblemDomain& a_indexSpace);

  ///
  virtual int refToFiner(const ProblemDomain&) const;

private:
  Vector<RefCountedPtr<LevelData<FluxBox> > >    m_eta;
  LevelData<FArrayBox>                           m_phic;
  Vector<ProblemDomain>                          m_domains;
  Vector<DisjointBoxLayout>                      m_boxes;
  Vector<Real>                                   m_dx;
  // refinement to next coarser level
  Vector<int>                                    m_refRatios;
  BCFunc  m_bc;
  Real m_alpha, m_beta;

  ///weak construction is bad
  ResistivityOpFactory()
  {
    MayDay::Error("invalid operator");
  }
  //copy constructor and operator= disallowed for all the usual reasons
  ResistivityOpFactory(const ResistivityOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const ResistivityOpFactory& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
