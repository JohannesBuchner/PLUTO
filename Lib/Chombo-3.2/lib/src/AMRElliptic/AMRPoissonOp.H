#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRPOISSONOP_H_
#define _AMRPOISSONOP_H_

#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "QuadCFInterp.H"
#include "LevelFluxRegister.H"
#include "LevelDataOps.H"
#include "AMRMultiGrid.H"
#include "BCFunc.H"
#include "BaseLevelTGA.H"
#include "CFRegion.H"
#include "AMRIO.H"
#include "CornerCopier.H"

#include "NamespaceHeader.H"

///
/**
   Operator for solving (alpha I + beta*Laplacian)(phi) = rho
   over an AMR hierarchy.
*/
class AMRPoissonOp : public LevelTGAHelmOp<LevelData<FArrayBox> , FluxBox>
{
public:

  /**
     \name AMRPoissonOp functions */
  /*@{*/

  ///
  /**
   */
  AMRPoissonOp()
  {
  }

  ///
  /**
   */
  virtual ~AMRPoissonOp()
  {
  }

  ///
  /** full define function for AMRLevelOp with both coarser and finer levels */
  void define(const DisjointBoxLayout&   a_grids,
              const DisjointBoxLayout&   a_gridsFiner,
              const DisjointBoxLayout&   a_gridsCoarser,
              const Real&                a_dxLevel,
              int                        a_refRatio,
              int                        a_refRatioFiner,
              const ProblemDomain&       a_domain,
              BCHolder                   a_bc,
              const Copier&              a_exchange,
              const CFRegion&            a_cfregion);

  /** full define function for AMRLevelOp with finer levels, but no coarser */
  void define(const DisjointBoxLayout&   a_grids,
              const DisjointBoxLayout&   a_gridsFiner,
              const Real&                a_dxLevel,
              int                        a_refRatio, // dummy arg, send in 1
              int                        a_refRatioFiner,
              const ProblemDomain&       a_domain,
              BCHolder                   a_bc,
              const Copier&              a_exchange,
              const CFRegion&            a_cfregion);

  ///
  /**
     define function for AMRLevelOp which has no finer AMR level
  */
  void define(const DisjointBoxLayout&   a_grids,
              const DisjointBoxLayout&   a_baseBAPtr,
              const Real&                a_dxLevel,
              int                        a_refRatio,
              const ProblemDomain&       a_domain,
              BCHolder                   a_bc,
              const Copier&              a_exchange,
              const CFRegion&            a_cfregion);

  ///
  /**
     define function for AMRLevelOp which has no finer or coarser AMR level
   */
  void define(const DisjointBoxLayout&   a_grids,
              const Real&                a_dx,
              const ProblemDomain&       a_domain,
              BCHolder                   a_bc,
              const Copier&              a_exchange,
              const CFRegion&            a_cfregion);

  ///
  /**
     define function for AMRLevelOp which has no finer or coarser AMR level
   */
  void define(const DisjointBoxLayout&   a_grids,
              const Real&                a_dx,
              const ProblemDomain&       a_domain,
              BCHolder                   a_bc);

  /// Full define function that mimics the old PoissonOp.
  /**
      Makes all coarse-fine information and sets internal variables
  */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_baseBAPtr,
              Real                     a_dxLevel,
              int                      a_refRatio,
              const ProblemDomain&     a_domain,
              BCHolder                 a_bc);

  virtual void residual(LevelData<FArrayBox>&       a_lhs,
                        const LevelData<FArrayBox>& a_phi,
                        const LevelData<FArrayBox>& a_rhs,
                        bool                        a_homogeneous = false);

 /// despite what you might think, the "I" here means "Ignore the coarse-fine boundary"
  virtual void residualI(LevelData<FArrayBox>&       a_lhs,
                         const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>& a_rhs,
                         bool                        a_homogeneous = false);

  virtual void preCond(LevelData<FArrayBox>&       a_correction,
                       const LevelData<FArrayBox>& a_residual);

  virtual void applyOp(LevelData<FArrayBox>&       a_lhs,
                       const LevelData<FArrayBox>& a_phi,
                       bool                        a_homogeneous = false);

  /// despite what you might think, the "I" here means "Ignore the coarse-fine boundary"
  virtual void applyOpI(LevelData<FArrayBox>&       a_lhs,
                        const LevelData<FArrayBox>& a_phi,
                        bool                        a_homogeneous = false);

  virtual void applyOpNoBoundary(LevelData<FArrayBox>& a_lhs,
                                 const LevelData<FArrayBox>& a_rhs) ;

  virtual void create(LevelData<FArrayBox>&       a_lhs,
                      const LevelData<FArrayBox>& a_rhs);

  virtual void createCoarsened(LevelData<FArrayBox>&       a_lhs,
                               const LevelData<FArrayBox>& a_rhs,
                               const int&                  a_refRat);

  virtual void assign(LevelData<FArrayBox>&       a_lhs,
                      const LevelData<FArrayBox>& a_rhs);

  virtual void assignLocal(LevelData<FArrayBox>&       a_lhs,
                           const LevelData<FArrayBox>& a_rhs);

  virtual void buildCopier(Copier& a_copier,
                           const LevelData<FArrayBox>& a_lhs,
                           const LevelData<FArrayBox>& a_rhs);

  virtual void assignCopier(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_rhs,
                            const Copier&               a_copier);

  virtual void zeroCovered(LevelData<FArrayBox>& a_lhs,
                           LevelData<FArrayBox>& a_rhs,
                           const Copier&         a_copier);

  virtual Real dotProduct(const LevelData<FArrayBox>& a_1,
                          const LevelData<FArrayBox>& a_2);
  /* multiple dot products (for GMRES) */
  virtual void mDotProduct(const LevelData<FArrayBox>& a_1,
                           const int a_sz,
                           const LevelData<FArrayBox> a_2[],
                           Real a_mdots[]);

  virtual void incr(LevelData<FArrayBox>&       a_lhs,
                    const LevelData<FArrayBox>& a_x,
                    Real                        a_scale);

  virtual void axby(LevelData<FArrayBox>&       a_lhs,
                    const LevelData<FArrayBox>& a_x,
                    const LevelData<FArrayBox>& a_y,
                    Real                        a_a,
                    Real                        a_b);

  virtual void scale(LevelData<FArrayBox>& a_lhs,
                     const Real&           a_scale);

  virtual Real norm(const LevelData<FArrayBox>& a_x,
                    int                         a_ord);

  virtual Real localMaxNorm(const LevelData<FArrayBox>& a_x);

  virtual void setToZero( LevelData<FArrayBox>& a_x);
  /*@}*/

  /**
     \name MGLevelOp functions */
  /*@{*/

  virtual void relax(LevelData<FArrayBox>&       a_e,
                     const LevelData<FArrayBox>& a_residual,
                     int                         a_iterations);

  virtual void createCoarser(LevelData<FArrayBox>&       a_coarse,
                             const LevelData<FArrayBox>& a_fine,
                             bool                        a_ghosted);
  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                LevelData<FArrayBox>&       a_phiFine,
                                const LevelData<FArrayBox>& a_rhsFine);

  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                const LevelData<FArrayBox>& a_correctCoarse);

  /*@}*/

  /**
     \name AMRLevelOp functions */
  /*@{*/

  /** returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToCoarser()
  {
    return m_refToCoarser;
  }

  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMRResidual(LevelData<FArrayBox>&              a_residual,
                           const LevelData<FArrayBox>&        a_phiFine,
                           const LevelData<FArrayBox>&        a_phi,
                           const LevelData<FArrayBox>&        a_phiCoarse,
                           const LevelData<FArrayBox>&        a_rhs,
                           bool                               a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** residual assuming no more coarser AMR levels */

  virtual void AMRResidualNC(LevelData<FArrayBox>&              a_residual,
                             const LevelData<FArrayBox>&        a_phiFine,
                             const LevelData<FArrayBox>&        a_phi,
                             const LevelData<FArrayBox>&        a_rhs,
                             bool                               a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMRResidualNF(LevelData<FArrayBox>&       a_residual,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_phiCoarse,
                             const LevelData<FArrayBox>& a_rhs,
                             bool                        a_homogeneousPhysBC);

  ///
  /**
     Apply the AMR operator, including coarse-fine matching
  */
  virtual void AMROperator(LevelData<FArrayBox>&              a_LofPhi,
                           const LevelData<FArrayBox>&        a_phiFine,
                           const LevelData<FArrayBox>&        a_phi,
                           const LevelData<FArrayBox>&        a_phiCoarse,
                           bool                               a_homogeneousDomBC,
                           AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  ///
  /**
      Apply the AMR operator, including coarse-fine matching
      assume no coarser AMR level
  */
  virtual void AMROperatorNC(LevelData<FArrayBox>&              a_LofPhi,
                             const LevelData<FArrayBox>&        a_phiFine,
                             const LevelData<FArrayBox>&        a_phi,
                             bool                               a_homogeneousBC,
                             AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  ///
  /**
      Apply the AMR operator, including coarse-fine matching.
      assume no finer AMR level
  */
  virtual void AMROperatorNF(LevelData<FArrayBox>&       a_LofPhi,
                             const LevelData<FArrayBox>& a_phi,
                             const LevelData<FArrayBox>& a_phiCoarse,
                             bool                        a_homogeneousBC);

  /**
      a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
      it is assumed that a_resCoarse has already been filled in with the coarse
      version of AMRResidualNF and that this operation is free to overwrite
      in the overlap regions.
  */
  virtual void AMRRestrict(LevelData<FArrayBox>&       a_resCoarse,
                           const LevelData<FArrayBox>& a_residual,
                           const LevelData<FArrayBox>& a_correction,
                           const LevelData<FArrayBox>& a_coarseCorrection,
                           bool a_skip_res = false );

  virtual void AMRRestrictS(LevelData<FArrayBox>&       a_resCoarse,
                            const LevelData<FArrayBox>& a_residual,
                            const LevelData<FArrayBox>& a_correction,
                            const LevelData<FArrayBox>& a_coarseCorrection,
                            LevelData<FArrayBox>&       a_scratch,
                            bool a_skip_res = false );

  /**
      a_correction += I[h->h](a_coarseCorrection)
  */
  virtual void AMRProlong(LevelData<FArrayBox>&       a_correction,
                          const LevelData<FArrayBox>& a_coarseCorrection);

  /**
      optimization of AMRProlong that sends in the existing temporary and copier
  */
  virtual void AMRProlongS(LevelData<FArrayBox>&       a_correction,
                           const LevelData<FArrayBox>& a_coarseCorrection,
                           LevelData<FArrayBox>&       a_temp,
                           const Copier&               a_copier);
  /**
      optimization of AMRProlong that sends in the existing temporary and copier -- higher order
  */
  virtual void AMRProlongS_2(LevelData<FArrayBox>&       a_correction,
                             const LevelData<FArrayBox>& a_coarseCorrection,
                             LevelData<FArrayBox>&       a_temp,
                             const Copier&               a_copier,
                             const Copier&         a_cornerCopier,
                             const AMRLevelOp<LevelData<FArrayBox> >*  a_crsOp
                             );

  /**
      a_residual = a_residual - L(a_correction, a_coarseCorrection)
  */
  virtual void AMRUpdateResidual(LevelData<FArrayBox>&       a_residual,
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

  /// For tga to reset stuff
  /**
  */
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta);

  /// Change boundary conditions
  virtual void setBC(const BCHolder& a_bc);

  /// For tga stuff---in this case a noop
  virtual void        diagonalScale(LevelData<FArrayBox>& a_rhs,
                                    bool a_kappaWeighted)
  {
  }

  /// For tga stuff---in this case a noop
  virtual void divideByIdentityCoef(LevelData<FArrayBox>& a_rhs)
  {
  }

  virtual void fillGrad(const LevelData<FArrayBox>& a_phi)
  {
    ;//does not apply here
  }

  virtual void reflux(const LevelData<FArrayBox>&        a_phiFine,
                      const LevelData<FArrayBox>&        a_phi,
                      LevelData<FArrayBox>&              a_residual,
                      AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  virtual void getFlux(FluxBox&                    a_flux,
                       const LevelData<FArrayBox>& a_data,
                       const Box&                  a_grid,
                       const DataIndex&            a_dit,
                       Real                        a_scale)
  {
    const FArrayBox& data = a_data[a_dit];
    a_flux.define(a_grid, a_data.nComp());
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        getFlux(a_flux[idir], data, a_flux[idir].box(), idir, 1);
        a_flux[idir] *= a_scale;
      }
  }

  virtual void write(const LevelData<FArrayBox>* a_data,
                     const char*                 a_filename);

  /*@}*/

  /// public constants
  Real m_alpha, m_beta, m_aCoef, m_bCoef;

  // needed for homogeneous interpolation
  // set by the factory
  Real  m_dxCrse;

  Vector<IntVect> m_colors;
  static int s_exchangeMode;
  static int s_relaxMode;
  static int s_maxCoarse;
  static int s_prolongType;

  virtual Real dx() const
  {
    return m_dx;
  }

protected:
  Real                    m_dx;
  ProblemDomain           m_domain;

  LevelDataOps<FArrayBox> m_levelOps;

  BCHolder                m_bc;

  CFRegion                m_cfregion;
  Copier                  m_exchangeCopier;
  QuadCFInterp            m_interpWithCoarser;

  LevelFluxRegister       m_levfluxreg;

  DisjointBoxLayout       m_coarsenedMGrids;

  int                     m_refToCoarser;
  int                     m_refToFiner;

  virtual void levelGSRB(LevelData<FArrayBox>&       a_phi,
                         const LevelData<FArrayBox>& a_rhs);

  virtual void levelMultiColor(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs);

  virtual void looseGSRB(LevelData<FArrayBox>&       a_phi,
                         const LevelData<FArrayBox>& a_rhs);

  virtual void overlapGSRB(LevelData<FArrayBox>&       a_phi,
                           const LevelData<FArrayBox>& a_rhs);

  virtual void levelGSRBLazy(LevelData<FArrayBox>&       a_phi,
                             const LevelData<FArrayBox>& a_rhs);

  virtual void levelJacobi(LevelData<FArrayBox>&       a_phi,
                           const LevelData<FArrayBox>& a_rhs);

  virtual void homogeneousCFInterp(LevelData<FArrayBox>& a_phif);

  virtual void homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                                   const DataIndex&      a_datInd,
                                   int                   a_idir,
                                   Side::LoHiSide        a_hiorlo);

  virtual void singleBoxCFInterp(FArrayBox& a_phi);

  virtual void interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                               const DataIndex&      a_datInd,
                               const int             a_idir,
                               const Side::LoHiSide  a_hiorlo,
                               const IntVectSet&     a_interpIVS);

  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_data,
                       const Box&       a_edgebox,
                       int              a_dir,
                       int              a_ref = 1) const ;

  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_data,
                       int              a_dir,
                       int              a_ref = 1) const ;
};

///
/**
   Factory to create AMRPoissonOps
 */
class AMRPoissonOpFactory: public AMRLevelOpFactory<LevelData<FArrayBox> >
{
public:
  virtual ~AMRPoissonOpFactory()
  {
  }

  ///
  /**
     a_coarseDomain is the domain at the coarsest level.
     a_grids is the AMR  hierarchy.
     a_refRatios are the refinement ratios between levels.  The ratio lives
         with the coarser level so a_refRatios[ilev] is the ratio between
         ilev and ilev+1
     a_coarseDx is the grid spacing at the coarsest level.
     a_bc holds the boundary conditions.
     a_alpha is the identity coefficient
     a_beta is the laplacian coefficient.
  */
  void define(const ProblemDomain&             a_coarseDomain,
              const Vector<DisjointBoxLayout>& a_grids,
              const Vector<int>&               a_refRatios,
              const Real&                      a_coarsedx,
              BCHolder                         a_bc,
              Real                             a_alpha = 0.0,
              Real                             a_beta  = 1.0);

  // regular multigrid definition function --deprecated
  void define(const ProblemDomain&     a_domain,
              const DisjointBoxLayout& a_grid,
              const Real&              a_dx,
              BCHolder                 a_bc,
              int                      a_maxDepth = -1,
              Real                     a_alpha = 0.0,
              Real                     a_beta  = 1.0);

  ///
  virtual MGLevelOp<LevelData<FArrayBox> >* MGnewOp(const ProblemDomain& a_FineindexSpace,
                                                     int                  a_depth,
                                                     bool                 a_homoOnly = true);

  ///
  virtual  AMRLevelOp<LevelData<FArrayBox> >* AMRnewOp(const ProblemDomain& a_indexSpace);

  ///
  virtual int refToFiner(const ProblemDomain& a_domain) const;

private:
  Vector<ProblemDomain>     m_domains;
  Vector<DisjointBoxLayout> m_boxes;

  Vector<Real> m_dx;
  Vector<int>  m_refRatios; // refinement to next coarser level

  BCHolder m_bc;

  Real m_alpha;
  Real m_beta;

  Vector<Copier>   m_exchangeCopiers;
  Vector<CFRegion> m_cfregion;
};

#include "NamespaceFooter.H"
#endif
