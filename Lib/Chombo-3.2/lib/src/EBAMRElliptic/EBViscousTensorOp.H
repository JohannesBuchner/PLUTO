#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBVISCOUSTENSOROP_H_
#define _EBVISCOUSTENSOROP_H_

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

#include "EBLevelDataOps.H"
#include "BaseEBBC.H"
#include "AMRTGA.H"
#include "BaseDomainBC.H"
#include "CFIVS.H"
#include "EBFastFR.H"
#include "EBMGAverage.H"
#include "EBMGInterp.H"
#include "EBTensorCFInterp.H"
#include "EBCoarsen.H"
#include "PolyGeom.H"
#include "EBQuadCFInterp.H"
#include "EBLevelGrid.H"
#include "EBAMRIO.H"
#include "DivergenceStencil.H"

#include "NamespaceHeader.H"

// Forward declaration.
class EBViscousTensorOp;

//! \class EBViscousTensorOp
//! This class implements an operator that solves the equation
//! alpha a I + (divF) = alpha a I + beta*div(eta(grad B + grad B^T) + lambda* I div B ) = rhs
//! where beta, lambda,  and eta are incorporated into the flux F. It uses
//! the AMRLevelOp interface.
//class EBViscousTensorOp: public AMRLevelOp<LevelData<EBCellFAB> >
class EBViscousTensorOp: public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB >
{
  public:

  static int s_whichLev;
  static int s_step;
  static bool s_setIntersectionsToZero;


  //! Constructs a viscous tensor operator using the given data. This
  //! constructor is for time-independent a and b coefficients.
  EBViscousTensorOp(const EBLevelGrid &                                a_eblgFine,
                    const EBLevelGrid &                                a_eblg,
                    const EBLevelGrid &                                a_eblgCoar,
                    const EBLevelGrid &                                a_eblgCoarMG,
                    const Real&                                        a_alpha,
                    const Real&                                        a_beta,
                    const RefCountedPtr<LevelData<EBCellFAB> >&        a_acoef,
                    const RefCountedPtr<LevelData<EBFluxFAB> >&        a_eta,
                    const RefCountedPtr<LevelData<EBFluxFAB> >&        a_lambda,
                    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_etaIrreg,
                    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_lambdaIrreg,
                    const Real&                                        a_dx,
                    const Real&                                        a_dxCoar,
                    const int&                                         a_refToFine,
                    const int&                                         a_refToCoar,
                    const RefCountedPtr<ViscousBaseDomainBC>&          a_domainBC,
                    const RefCountedPtr<ViscousBaseEBBC>&              a_ebBC,
                    const bool&                                        a_hasMGObjects,
                    const IntVect&                                     a_ghostCellsPhi,
                    const IntVect&                                     a_ghostCellsRHS);

  //! Constructs a viscous tensor operator using the given data. This
  //! constructor is for a time-dependent a coefficient.
  EBViscousTensorOp(const EBLevelGrid &                                a_eblgFine,
                    const EBLevelGrid &                                a_eblg,
                    const EBLevelGrid &                                a_eblgCoar,
                    const EBLevelGrid &                                a_eblgCoarMG,
                    const Real&                                        a_alpha,
                    const Real&                                        a_beta,
                    const RefCountedPtr<LevelData<EBCellFAB> >&        a_acoef0,
                    const RefCountedPtr<LevelData<EBCellFAB> >&        a_acoef1,
                    const RefCountedPtr<LevelData<EBCellFAB> >&        a_acoef,
                    const RefCountedPtr<LevelData<EBFluxFAB> >&        a_eta,
                    const RefCountedPtr<LevelData<EBFluxFAB> >&        a_lambda,
                    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_etaIrreg,
                    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_lambdaIrreg,
                    const Real&                                        a_dx,
                    const Real&                                        a_dxCoar,
                    const int&                                         a_refToFine,
                    const int&                                         a_refToCoar,
                    const RefCountedPtr<ViscousBaseDomainBC>&          a_domainBC,
                    const RefCountedPtr<ViscousBaseEBBC>&              a_ebBC,
                    const bool&                                        a_hasMGObjects,
                    const IntVect&                                     a_ghostCellsPhi,
                    const IntVect&                                     a_ghostCellsRHS);

  //! Destructor.
  virtual ~EBViscousTensorOp();

  ///for tga to reset stuff
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta);

#ifdef CH_USE_HDF5
  virtual void outputLevel(LevelData<EBCellFAB>& a_rhs, string& a_name)
  {
    char filename[1024];
    sprintf(filename, "%s.ebvto.step%d.lev.%d.hdf5",a_name.c_str(), s_step, s_whichLev );
    writeEBLevelname(&a_rhs, filename);
  }


  virtual void outputAMR(Vector< LevelData<EBCellFAB>* >& a_rhs, string& a_name)
  {
    char filename[1024];
    sprintf(filename, "%s.ebvto.step%d.lev.%d.hdf5", a_name.c_str(), s_step, s_whichLev );
    writeEBAMRname(&a_rhs, filename);
  }
#endif

  ///it's tga's world---we just live in it.
  virtual void kappaScale(LevelData<EBCellFAB> & a_rhs)
  {
    EBLevelDataOps::kappaWeight(a_rhs);
  }

  /// compute (tau dot grad u) (for energy equation)
  void getShearStressDotGradU(LevelData<EBCellFAB>      & a_source,
                              const LevelData<EBCellFAB>& a_gradU,
                             int a_level);

  /// compute volfrac(sigma dot grad u) (for energy equation)
  void getKappaDivSigmaU(LevelData<EBCellFAB>      & a_divSigmaU,
                         const LevelData<EBCellFAB>& a_velocity,
                         const LevelData<EBCellFAB>* a_veloCoar,
                         int a_level);

  //average surrounding coeffcients onto faces (simple averages)
  void getCellCenteredCoefficients(LevelData<EBCellFAB>&    a_etaCell,
                                   LevelData<EBCellFAB>& a_lambdaCell);

  //get stress at cell centers given gradient of velocity
  void getCCSigma(LevelData<EBCellFAB>      & a_source,
                  const LevelData<EBCellFAB>& a_gradU,
                  const LevelData<EBCellFAB>& a_eta,
                  const LevelData<EBCellFAB>& a_lambda);
  ///another tgaism
  virtual void diagonalScale(LevelData<EBCellFAB> & a_rhs,
                             bool a_kappaWeighted)
  {
    if (a_kappaWeighted)
      EBLevelDataOps::kappaWeight(a_rhs);
    for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            int isrc = 0; int idst = idir; int inco = 1;
            a_rhs[dit()].mult((*m_acoef)[dit()], isrc, idst, inco);
          }
      }
  }

  virtual void divideByIdentityCoef(LevelData<EBCellFAB> & a_rhs)
  {
    for (DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
      {
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            int isrc = 0; int idst = idir; int inco = 1;
            a_rhs[dit()].divide((*m_acoef)[dit()], isrc, idst, inco);
          }
      }
  }

  ///a leveltgaism
  virtual void fillGrad(const LevelData<EBCellFAB>& a_phi);

  ///another leveltgaism
  virtual void getFlux(EBFluxFAB&                    a_flux,
                       const LevelData<EBCellFAB>&   a_data,
                       const Box&                    a_grid,
                       const DataIndex&              a_dit,
                       Real                          a_scale);

  void getVelDotSigma(LevelData<EBFluxFAB>      & a_velDotSigma,
                      const LevelData<EBFluxFAB>& a_vel,
                      const LevelData<EBFluxFAB>& a_sigma);

  void getFlux(EBFaceFAB&                    a_fluxCentroid,
               const EBCellFAB&              a_phi,
               const Box&                    a_ghostedBox,
               const Box&                    a_fabBox,
               const ProblemDomain&          a_domain,
               const EBISBox&                a_ebisBox,
               const Real&                   a_dx,
               const DataIndex&              a_datInd,
               const int&                    a_idir);

  // This operator may be time-dependent.
  void setTime(Real a_oldTime, Real a_mu, Real a_dt);

  /// apply operator without any boundary or coarse-fine boundary conditions and no finer level
  virtual void applyOpNoBoundary(LevelData<EBCellFAB>&       a_opPhi,
                                 const LevelData<EBCellFAB>& a_phi)
  {
    s_turnOffBCs = true;
    applyOp(a_opPhi, a_phi, true);
    s_turnOffBCs = false;
  }

  //! This is called on multigrid operators when their AMR operators
  //! are altered.
  void finerOperatorChanged(const MGLevelOp<LevelData<EBCellFAB> >& a_operator,
                            int a_coarseningFactor);

  //functions used by the wider world

  ///
  /** a_residual = a_rhs - L(a_phiFine, a_phi)   no coaser AMR level*/
  void AMRResidualNC(LevelData<EBCellFAB>&       a_residual,
                     const LevelData<EBCellFAB>& a_phiFine,
                     const LevelData<EBCellFAB>& a_phi,
                     const LevelData<EBCellFAB>& a_rhs,
                     bool a_homogeneousBC,
                     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);


  ///
  /** apply AMR operator   no coaser AMR level*/
  void AMROperatorNC(LevelData<EBCellFAB>&       a_LofPhi,
                     const LevelData<EBCellFAB>& a_phiFine,
                     const LevelData<EBCellFAB>& a_phi,
                     bool a_homogeneousBC,
                     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

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


  ///
  /**
     this is the linearop function.  CFBC is set to homogeneous.  phic is null
  */
  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       bool                              a_homogeneousPhysBC);

  
  void incrOpWithExternalFlux(EBCellFAB               & a_lhs,
                              const DataIndex         & a_dit,
                              const BaseIVFAB<Real>   & a_ebflux);

  
  /// if the flux is null, it is not used
  void applyOp(LevelData<EBCellFAB>             & a_lhs,
               const LevelData<EBCellFAB>       & a_phi,
               bool                               a_homogeneous,
               const LevelData<BaseIVFAB<Real> >* a_ebflux);
  ///
  /**
   */
  virtual void create(LevelData<EBCellFAB>&       a_lhs,
                      const LevelData<EBCellFAB>& a_rhs);

  ///
  virtual void createCoarsened(LevelData<EBCellFAB>&       a_lhs,
                               const LevelData<EBCellFAB>& a_rhs,
                               const int&                  a_refRat);

  ///
  Real
  AMRNorm(const LevelData<EBCellFAB>& a_coarResid,
          const LevelData<EBCellFAB>& a_fineResid,
          const int& a_refRat,
          const int& a_ord);


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

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToCoarser();

  ///
  /** Refinement ratio between this level and coarser level.
      Returns 1 when there are no coarser AMRLevelOp objects */
  virtual int refToFiner();

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMRResidual(LevelData<EBCellFAB>& a_residual,
                           const LevelData<EBCellFAB>& a_phiFine,
                           const LevelData<EBCellFAB>& a_phi,
                           const LevelData<EBCellFAB>& a_phiCoarse,
                           const LevelData<EBCellFAB>& a_rhs,
                           bool a_homogeneousBC,
                           AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMRResidualNF(LevelData<EBCellFAB>& a_residual,
                             const LevelData<EBCellFAB>& a_phi,
                             const LevelData<EBCellFAB>& a_phiCoarse,
                             const LevelData<EBCellFAB>& a_rhs,
                             bool a_homogeneousBC);


  ///
  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMROperator(LevelData<EBCellFAB>& a_LofPhi,
                           const LevelData<EBCellFAB>& a_phiFine,
                           const LevelData<EBCellFAB>& a_phi,
                           const LevelData<EBCellFAB>& a_phiCoarse,
                           bool a_homogeneousBC,
                           AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMROperatorNF(LevelData<EBCellFAB>& a_LofPhi,
                             const LevelData<EBCellFAB>& a_phi,
                             const LevelData<EBCellFAB>& a_phiCoarse,
                             bool a_homogeneousBC);

  ///
  /** a_resCoarse = I[h-2h] (a_residual - L(a_correction, a_coarseCorrection)) */
  virtual void AMRRestrict(LevelData<EBCellFAB>& a_resCoarse,
                           const LevelData<EBCellFAB>& a_residual,
                           const LevelData<EBCellFAB>& a_correction,
                           const LevelData<EBCellFAB>& a_coarseCorrection, 
                           bool a_skip_res = false );

  ///
  /** a_correction += I[2h->h](a_coarseCorrection) */
  virtual void AMRProlong(LevelData<EBCellFAB>&       a_correction,
                          const LevelData<EBCellFAB>& a_coarseCorrection);

  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                                 const LevelData<EBCellFAB>& a_correction,
                                 const LevelData<EBCellFAB>& a_coarseCorrection);

  ///static for code reuse in bcs
  static void 
  getFaceCenteredFluxStencil(VoFStencil                                & a_fluxStencil,
                             const RefCountedPtr<LevelData<EBFluxFAB> >& a_eta,
                             const RefCountedPtr<LevelData<EBFluxFAB> >& a_lambda,
                             const Real                                & a_dx,
                             const EBLevelGrid                         & a_eblg,
                             const FaceIndex                           & a_face,
                             const DataIndex                           & a_dit,
                             int a_ivar);

  ///static for code reuse in bcs
  static void 
  getFluxStencil(VoFStencil                                & a_fluxStencil,
                 const RefCountedPtr<LevelData<EBFluxFAB> >& a_eta,
                 const RefCountedPtr<LevelData<EBFluxFAB> >& a_lambda,
                 const Real                                & a_dx,
                 const EBLevelGrid                         & a_eblg,
                 const FaceIndex                           & a_face,
                 const DataIndex                           & a_dit,
                 int a_ivar);

  ///static for code reuse in bcs
  static void
  getDivergenceStencil(VoFStencil&      a_divStencil,
                       const FaceIndex& a_face,
                       const DataIndex& a_dit,
                       const Real     & a_dx,
                       const EBLevelGrid& a_eblg);

  ///static for code reuse in bcs
  static void
  getGradientStencil(VoFStencil&  a_gradStencil,
                     int a_ivar,
                     int a_diffDir,
                     const FaceIndex& a_face,
                     const DataIndex& a_dit,
                     const Real     & a_dx,
                     const EBLevelGrid& a_eblg);

  ///
  static void doLazyRelax(bool a_doLazyRelax)
  {
    s_doLazyRelax = a_doLazyRelax;
  }
  static void setForceNoEBCF(bool a_forceNoEBCF)
  {
    s_forceNoEBCF = a_forceNoEBCF;
  }

  //! (Re)define the stencils for the given coefficients.
  void defineStencils();

  protected:

  //internal functions not usable without significant context
  //  void fillGrad(const LevelData<EBCellFAB>&       a_phi);

  void incrOpRegularDir(EBCellFAB&             a_lhs,
                        const EBCellFAB&       a_phi,
                        const bool&            a_homogeneous,
                        const int&             a_dir,
                        const DataIndex&       a_datInd);

  void applyOpIrregular(EBCellFAB&             a_lhs,
                        const EBCellFAB&       a_phi,
                        const bool&            a_homogeneous,
                        const DataIndex&       a_datInd);

  void getFlux(FArrayBox&                    a_flux,
               const FArrayBox&              a_phi,
               const FArrayBox&              a_gradPhi,
               const Box&                    a_faceBox,
               const int&                    a_idir,
               const DataIndex&              a_datInd);
  void
  cfinterp(const LevelData<EBCellFAB>&       a_phi,
           const LevelData<EBCellFAB>&       a_phiCoarse);

  void reflux(const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              LevelData<EBCellFAB>& residual,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  void homogeneousCFInterp(LevelData<EBCellFAB>&   a_phi);

  void getVoFStencil(VoFStencil&      a_vofStencil,
                     const VolIndex&  a_vof,
                     const DataIndex& a_dit,
                     int              a_ivar);

  void getDivergenceStencil(VoFStencil&  a_gradStencil,
                            const FaceIndex& a_face,
                            const DataIndex& a_dit);

  void getGradientStencil(VoFStencil&  a_gradStencil,
                          int a_ivar,
                          int a_diffDir,
                          const FaceIndex& a_face,
                          const DataIndex& a_dit);


  void gsrbColor(LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_lph,
                 const LevelData<EBCellFAB>& a_rhs,
                 const IntVect&              a_color);

  void getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                                  const FaceIndex& a_face,
                                  const DataIndex& a_dit,
                                  int a_ivar);

  void cellGrad(EBCellFAB&             a_gradPhi,
                const  EBCellFAB&      a_phi,
                const Box&             a_grid,
                const DataIndex&       a_datInd);

  void
  incrementFRCoar(EBFastFR& a_fr,
                  const LevelData<EBCellFAB>& a_phiFine,
                  const LevelData<EBCellFAB>& a_phi);

  //simple averaging
  void
  averageToCells(LevelData<EBCellFAB>&      a_cellCoef,
                 const LevelData<EBFluxFAB>& a_faceCoef,
                 const LevelData<BaseIVFAB<Real> >& a_irregCoef);

  void
  incrementFRFine(EBFastFR& a_fr,
                  const LevelData<EBCellFAB>& a_phiFine,
                  const LevelData<EBCellFAB>& a_phi,
                  AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  //member data, bloated though it may be
  static bool s_turnOffBCs;
  static bool s_forceNoEBCF;

  static bool s_doLazyRelax;
  //locations of coarse fine interfaces
  LayoutData<CFIVS>                            m_loCFIVS[CH_SPACEDIM];
  LayoutData<CFIVS>                            m_hiCFIVS[CH_SPACEDIM];
  LayoutData<RefCountedPtr<DivergenceStencil> >               m_divergenceStencil;
  // will need these to do tangential gradient computations
  LayoutData<TensorFineStencilSet> m_hiTanStencilSets[SpaceDim];
  LayoutData<TensorFineStencilSet> m_loTanStencilSets[SpaceDim];

  IntVect m_ghostPhi, m_ghostRHS;
  //grid information
  EBLevelGrid                                  m_eblgFine;
  EBLevelGrid                                  m_eblg;
  EBLevelGrid                                  m_eblgCoar;
  EBLevelGrid                                  m_eblgCoarMG;

  //operator coefficients
  Real                                         m_alpha;
  Real                                         m_beta;
  //weights that get multiplied by alpha
  LayoutData<BaseIVFAB<Real> >       m_alphaDiagWeight;
  //weights that get multiplied by beta
  LayoutData<BaseIVFAB<Real> >       m_betaDiagWeight;
  LevelData<EBCellFAB> m_zeroCoarse;

  //! "Current" (time-interpolated) value of the a coefficient. For a
  //! time-independent a coefficient, this is where the coefficient lives.
  RefCountedPtr<LevelData<EBCellFAB> >         m_acoef;

  //! The a coefficient at the beginning of the time step. This is NULL if
  //! the conductivity operator is time independent.
  RefCountedPtr<LevelData<EBCellFAB> >         m_acoef0;

  //! The a coefficient at the end of the time step. This is NULL if
  //! the conductivity operator is time independent.
  RefCountedPtr<LevelData<EBCellFAB> >         m_acoef1;

  RefCountedPtr<LevelData<EBFluxFAB> >         m_eta;
  RefCountedPtr<LevelData<EBFluxFAB> >         m_lambda;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_etaIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >  m_lambdaIrreg;


  //flux register with finer level
  EBFastFR                                     m_fastFR;

  //grid spacing
  //dxcoar not implied by ref ratio--we could be in some MG level where dx > dxcoar
  //needed for homogeneous cf interp
  Real                                         m_dx;
  Real                                         m_dxCoar;
  Real getSafety();
  void calculateAlphaWeight();
  void calculateRelaxationCoefficient();

  //refinement ratios and whether this object has multigrid objects
  bool                                         m_hasFine;
  bool                                         m_hasCoar;
  int                                          m_refToFine;
  int                                          m_refToCoar;
  bool                                         m_hasMGObjects;

  //ghost cell information needed for ebstencil optimizations
  IntVect                                      m_ghostCellsPhi;
  IntVect                                      m_ghostCellsRHS;

  //stencil to apply irregular operator
  LayoutData<RefCountedPtr<EBStencil> >        m_opEBStencil[CH_SPACEDIM];

  //relaxation coefficent
  LevelData<EBCellFAB>                         m_relCoef;

  //gradient of solution at cell centers
  LevelData<EBCellFAB>                         m_grad;

  LayoutData<VoFIterator >                     m_vofIterIrreg;
  LayoutData<VoFIterator >                     m_vofIterMulti;
  //for domain boundary conditions at ir regular cells
  LayoutData<VoFIterator >                     m_vofIterDomLo[CH_SPACEDIM];
  LayoutData<VoFIterator >                     m_vofIterDomHi[CH_SPACEDIM];

  //for inhomogeneous cf interpolation
  RefCountedPtr<EBTensorCFInterp>              m_interpWithCoarser;

  //averaging operators
  EBMGAverage                                  m_ebAverage;
  EBMGAverage                                  m_ebAverageMG;

  //interpolation operators
  EBMGInterp                                   m_ebInterp;
  EBMGInterp                                   m_ebInterpMG;

  //boundary condition operators
  RefCountedPtr<ViscousBaseDomainBC>           m_domainBC;
  RefCountedPtr<ViscousBaseEBBC>               m_ebBC;
  Vector<IntVect>                              m_colors;

  //scrounged from EBLevelDataOps for performance reasons
  void
  averageCellToFace(EBFaceFAB           &      a_fluxData,
                    const EBCellFAB     &      a_cellData,
                    const Box           &      a_grid,
                    const EBISBox       &      a_ebisBox,
                    const ProblemDomain &      a_domain,
                    const DataIndex& a_dit,
                    int isrc, int idst, int inco,
                    bool a_interpolateToCentroid);


  //scrounged from EBLevelDataOps for performance reasons
  void
  averageCellToFace(LevelData<EBFluxFAB>&         a_fluxData,
                    const LevelData<EBCellFAB>&   a_cellData,
                    const DisjointBoxLayout&      a_grids,
                    const EBISLayout&             a_ebisl,
                    const ProblemDomain&          a_domain,
                    int isrc, int idst, int inco,
                    bool a_interpolateToCentroid);

  
  //scrounged from EBLevelDataOps for performance reasons
  void
  faceCenteredAverageCellsToFaces(EBFaceFAB           &      a_faceData,
                                  const EBCellFAB     &      a_cellData,
                                  const Box           &      ccFluxBox,
                                  const EBISBox       &      a_ebisBox,
                                  const ProblemDomain &      a_domain,
                                  const DataIndex& a_dit,
                                  int isrc, int idst, int inco);

  LayoutData<IntVectSet> m_ivsIrregCCFlux;

private:
  ///weak construction is bad
  EBViscousTensorOp()
  {
    MayDay::Error("invalid operator");
  }

  //copy constructor and operator= disallowed for all the usual reasons
  EBViscousTensorOp(const EBViscousTensorOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBViscousTensorOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};


#include "NamespaceFooter.H"
#endif
