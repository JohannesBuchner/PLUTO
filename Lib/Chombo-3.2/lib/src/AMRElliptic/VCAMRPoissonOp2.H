#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VCAMRPOISSONOP2_H_
#define _VCAMRPOISSONOP2_H_

#include "AMRPoissonOp.H"
#include "CoefficientInterpolator.H"

#include "NamespaceHeader.H"

///
/**
   Operator for solving variable-coefficient
   (alpha * aCoef(x) * I - beta * Div(bCoef(x) . Grad)) phi = rho
   over an AMR hierarchy.
*/
class VCAMRPoissonOp2 : public AMRPoissonOp
{
public:
  /**
     \name VCAMRPoissonOp2 functions */
  /*@{*/

  ///
  /**
   */
  VCAMRPoissonOp2()
  {
    m_lambdaNeedsResetting = true;
  }

  ///
  /**
   */
  virtual ~VCAMRPoissonOp2()
  {
  }

  ///
  virtual void residualI(LevelData<FArrayBox>&       a_lhs,
                         const LevelData<FArrayBox>& a_phi,
                         const LevelData<FArrayBox>& a_rhs,
                         bool                        a_homogeneous = false);

  ///
  virtual void preCond(LevelData<FArrayBox>&       a_correction,
                       const LevelData<FArrayBox>& a_residual);

  ///
  virtual void applyOpI(LevelData<FArrayBox>&       a_lhs,
                        const LevelData<FArrayBox>& a_phi,
                        bool                        a_homogeneous = false);

  virtual void applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_phi);

  /*@}*/

  /**
<<<<<<< VCAMRPoissonOp2.H
=======
     \name MGLevelOp functions */
  /*@{*/

  /**
>>>>>>> 1.6
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                LevelData<FArrayBox>&       a_phiFine,
                                const LevelData<FArrayBox>& a_rhsFine);

  /*@}*/

  /**
     \name VCAMRPoissonOp2 functions */
  /*@{*/

  /// For tga stuff
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta);

  /// Also calls reset lambda
  virtual void setCoefs(const RefCountedPtr<LevelData<FArrayBox> >& a_aCoef,
                        const RefCountedPtr<LevelData<FluxBox>   >& a_bCoef,
                        const Real& a_alpha,
                        const Real& a_beta);

  /// Should be called before the relaxation parameter is needed.
  virtual void resetLambda();

  /// Compute lambda once alpha, aCoef, beta, bCoef are defined
  virtual void computeLambda();

  virtual void reflux(const LevelData<FArrayBox>&        a_phiFine,
                      const LevelData<FArrayBox>&        a_phi,
                      LevelData<FArrayBox>&              a_residual,
                      AMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

  /*@}*/

  //! This is called on multigrid operators when their AMR operators
  //! are altered.
  void finerOperatorChanged(const MGLevelOp<LevelData<FArrayBox> >& a_operator,
                            int a_coarseningFactor);

  //! Returns identity coefficient data.
  LevelData<FArrayBox>& identityCoef()
  {
    return *m_aCoef;
  }

  //! Sets up a model that modifies b coefficient data when the operator's
  //! time is set.
  //! \param a_bCoefInterpolator A CoefficientInterpolator that will be used
  //!                            to compute the b coefficient at specific
  //!                            times.
  void setBCoefInterpolator(RefCountedPtr<CoefficientInterpolator<LevelData<FluxBox>, LevelData<FArrayBox> > >& a_bCoefInterpolator)
  {
    m_bCoefInterpolator = a_bCoefInterpolator;
  }

  //! Returns the B coefficient.
  LevelData<FluxBox>& BCoef()
  {
    return *m_bCoef;
  }

  // Allows access to the B coefficient interpolator.
  RefCountedPtr<CoefficientInterpolator<LevelData<FluxBox>, LevelData<FArrayBox> > > BCoefInterpolator()
  {
    return m_bCoefInterpolator;
  }

  //! Sets the time centering of the operator. This interpolates b coefficient
  //! data at the given time if an interpolator is set.
  void setTime(Real a_time);

  /// Identity operator spatially varying coefficient storage (cell-centered) --- if you change this call resetLambda()
  RefCountedPtr<LevelData<FArrayBox> > m_aCoef;

  /// Laplacian operator spatially varying coefficient storage (face-centered) --- if you change this call resetLambda()
  RefCountedPtr<LevelData<FluxBox> > m_bCoef;

  /// Reciprocal of the diagonal entry of the operator matrix
  LevelData<FArrayBox> m_lambda;

protected:
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  // Interpolator for b coefficient data.
  RefCountedPtr<CoefficientInterpolator<LevelData<FluxBox>, LevelData<FArrayBox> > > m_bCoefInterpolator;

  // Current time.
  Real m_time;

  // Does the relaxation coefficient need to be reset?
  bool m_lambdaNeedsResetting;

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

  /// computes flux over face-centered a_facebox.
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_data,
                       const FluxBox&   a_bCoef,
                       const Box&       a_facebox,
                       int              a_dir,
                       int              a_ref = 1) const ;
#if 0
  /// utility function which computes operator after all bc's have been set
  void computeOperatorNoBCs(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_phi);
#endif
};

///
/**
   Factory to create VCAMRPoissonOp2s
*/
class VCAMRPoissonOp2Factory: public AMRLevelOpFactory<LevelData<FArrayBox> >
{
public:
  VCAMRPoissonOp2Factory();

  virtual ~VCAMRPoissonOp2Factory()
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
     a_alpha is the identity constant coefficient
     a_beta is the laplacian constant coefficient.
     a_aCoef is the identity spatially varying coefficient
     a_bCoef is the laplacian spatially varying coefficient.
  */
  void define(const ProblemDomain&                           a_coarseDomain,
              const Vector<DisjointBoxLayout>&               a_grids,
              const Vector<int>&                             a_refRatios,
              const Real&                                    a_coarsedx,
              BCHolder                                       a_bc,
              const Real&                                    a_alpha,
              Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
              const Real&                                    a_beta,
              Vector<RefCountedPtr<LevelData<FluxBox> > >&   a_bCoef);

  //! Defines a factory for VCAMRPoissonOp2 which allows the operators to
  //! allocate their own coefficient data. \f$\alpha\f$ and \f$\beta\f$
  //! coefficients are initialized to 1.
  //! \param a_coarseDomain The domain at the coarsest level.
  //! \param a_grids The disjoint box layouts for the various refinement levels.
  //! \param a_refRatios The refinement ratios between levels.
  //! \param a_coarsedx The grid spacing at the coarsest level.
  //! \param a_bc The boundary condition imposed on the solution.
  //! \param a_ghostVect The ghost stencil to use in the created coefficient data.
  void define(const ProblemDomain&                           a_coarseDomain,
              const Vector<DisjointBoxLayout>&               a_grids,
              const Vector<int>&                             a_refRatios,
              const Real&                                    a_coarsedx,
              BCHolder                                       a_bc,
              const IntVect&                                 a_ghostVect);

  ///
  virtual MGLevelOp<LevelData<FArrayBox> >* MGnewOp(const ProblemDomain& a_FineindexSpace,
                                                    int                  a_depth,
                                                    bool                 a_homoOnly = true);

  ///
  virtual AMRLevelOp<LevelData<FArrayBox> >* AMRnewOp(const ProblemDomain& a_indexSpace);

  ///
  virtual int refToFiner(const ProblemDomain& a_domain) const;

  int m_coefficient_average_type;

private:
  void setDefaultValues();

  Vector<ProblemDomain>     m_domains;
  Vector<DisjointBoxLayout> m_boxes;

  Vector<Real> m_dx;
  Vector<int>  m_refRatios; // refinement to next coarser level

  BCHolder m_bc;

  Real m_alpha;
  Real m_beta;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_aCoef;
  Vector<RefCountedPtr<LevelData<FluxBox> > >   m_bCoef;

  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_lambda;

  Vector<Copier>   m_exchangeCopiers;
  Vector<CFRegion> m_cfregion;
};

#include "NamespaceFooter.H"
#endif
