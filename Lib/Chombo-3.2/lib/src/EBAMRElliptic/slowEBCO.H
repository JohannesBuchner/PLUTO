#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SLOWEBCO_H_
#define _SLOWEBCO_H_

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
#include "EBFastFR.H"
#include "EBMGAverage.H"
#include "EBMGInterp.H"
#include "PolyGeom.H"
#include "EBQuadCFInterp.H"
#include "EBLevelGrid.H"
#include "AMRTGA.H"
#include "AMRPoissonOp.H"
#include "CFRegion.H"
#include "NamespaceHeader.H"


///
/**
   Operator to solve (alpha a + beta div (b grad) )phi = rhs.   This follows the AMRLevelOp interface.
*/
class slowEBCO: public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB >
{
public:

  Real getSafety();

  //for tga to reset stuff
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta);

  //another tgaism
  virtual void diagonalScale(LevelData<EBCellFAB> & a_rhs,
                             bool a_kappaWeighted);
  virtual void divideByIdentityCoef(LevelData<EBCellFAB> & a_rhs);

  ///a leveltgaism
  virtual void fillGrad(const LevelData<EBCellFAB>& a_phi)
  {;}

  ///a leveltgaism
  virtual void getFlux(EBFluxFAB&                    a_flux,
                       const LevelData<EBCellFAB>&   a_data,
                       const Box&                    a_grid,
                       const DataIndex&              a_dit,
                       Real                          a_scale);

  void getFlux(EBFaceFAB&                    a_fluxCentroid,
               const EBCellFAB&              a_phi,
               const Box&                    a_ghostedBox,
               const Box&                    a_fabBox,
               const ProblemDomain&          a_domain,
               const EBISBox&                a_ebisBox,
               const Real&                   a_dx,
               const DataIndex&              a_datInd,
               const int&                    a_idir);

  ///
  virtual ~slowEBCO()
  {;}

  ///
  slowEBCO()
  {;}

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
     If you are approaching this operator from this interface, consider backing away and using
     slowEBCOFactory to generate these objects.   Really.
     a_eblgFine,    : grid at finer  level \\
     a_eblg,        : grid at this  level \\
     a_eblgCoar,    : grid at coarser level \\
     a_eblgCoarMG,  : grid at intermediate multigrid level \\
     a_domainBC,    : domain boundary conditions at this level   \\
     a_ebBC:          eb boundary conditions at this level   \\
     a_dx:            grid spacing at this level \\
     a_origin:        offset to lowest corner of the domain \\
     a_refToFine:     refinement ratio to finer level \\
     a_refToCoar:     refinement ratio to coarser level \\
     a_hasFiner:      true if there is a finer AMR level, false otherwise. \\
     a_hasCoarser:    true  if there is a coarser AMR level. \\
     a_hasCoarserMG:    true  if there is a coarser MultiGrid level. \\
     a_preCondIters:  number of iterations to do for pre-conditioning \\
     a_relaxType:     0 means point Jacobi, 1 is Gauss-Seidel. \\
     a_acoef:         coefficent of identity \\
     a_bcoef:          coefficient of gradient.\\
     a_ghostCellsPhi:  Number of ghost cells in phi, correction\\
     a_ghostCellsRhs:  Number of ghost cells in RHS, residual, lphi\\
     Ghost cell arguments are there for caching reasons.  Once you set them, an error is thrown if
     you send in data that does not match.
  */
  slowEBCO(const EBLevelGrid &                                  a_eblgFine,
                   const EBLevelGrid &                                  a_eblg,
                   const EBLevelGrid &                                  a_eblgCoar,
                   const EBLevelGrid &                                  a_eblgCoarMG,
                   const RefCountedPtr<EBQuadCFInterp>&                 a_quadCFI,
                   const RefCountedPtr<ConductivityBaseDomainBC>&       a_domainBC,
                   const RefCountedPtr<ConductivityBaseEBBC>&           a_ebBC,
                   const Real    &                                      a_dx,
                   const Real    &                                      a_dxCoar,
                   const int&                                           a_refToFine,
                   const int&                                           a_refToCoar,
                   const bool&                                          a_hasFine,
                   const bool&                                          a_hasCoar,
                   const bool&                                          a_hasMGObjects,
                   const bool&                                          a_layoutChanged,
                   const Real&                                          a_alpha,
                   const Real&                                          a_beta,
                   const RefCountedPtr<LevelData<EBCellFAB> >&          a_acoef,
                   const RefCountedPtr<LevelData<EBFluxFAB> >&          a_bcoef,
                   const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_bcoIrreg,
                   const IntVect&                                       a_ghostCellsPhi,
                   const IntVect&                                       a_ghostCellsRHS,
                   const int&                                           a_relaxType);

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

  ///
  /**
     This function assumes that coarse-fine boundary condtions have
     been dealt with.
  */
  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       const LevelData<EBCellFAB>* const a_phiCoarse,
                       const bool&                       a_homogeneousPhysBC,
                       const bool&                       a_homogeneousCFBC);

  /// virtual function called by LevelTGA
  virtual void applyOpNoBoundary(LevelData<EBCellFAB>&        a_opPhi,
                                 const LevelData<EBCellFAB>&  a_phi);

  ///
  /**
     this is the linearop function.  CFBC is set to homogeneous.  phic is null
  */
  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       bool                              a_homogeneousPhysBC);

  ///
  /**
   */
  virtual void create(LevelData<EBCellFAB>&       a_lhs,
                      const LevelData<EBCellFAB>& a_rhs);

  ///
  virtual void createCoarsened(LevelData<EBCellFAB>&       a_lhs,
                               const LevelData<EBCellFAB>& a_rhs,
                               const int&                  a_refRat);

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


  virtual void relaxGauSai(LevelData<EBCellFAB>&       a_e,
                          const LevelData<EBCellFAB>& a_residual,
                          int                         a_iterations);

  virtual void relaxPoiJac(LevelData<EBCellFAB>&       a_e,
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

  void reflux(LevelData<EBCellFAB>& a_residual,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  void gsrbColor(LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_lph,
                 const LevelData<EBCellFAB>& a_rhs,
                 const IntVect&              a_color);

  void getDivFStencil(VoFStencil&      a_vofStencil,
                      const VolIndex&  a_vof,
                      const DataIndex& a_dit);

  void getFluxStencil(VoFStencil&      a_fluxStencil,
                      const FaceIndex& a_face,
                      const DataIndex& a_dit);

  void getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                                  const FaceIndex& a_face,
                                  const DataIndex& a_dit);

  void incrOpRegularDir(EBCellFAB&             a_lhs,
                        const EBCellFAB&       a_phi,
                        const bool&            a_homogeneous,
                        const int&             a_dir,
                        const DataIndex&       a_datInd);
  void applyOpIrregular(EBCellFAB&             a_lhs,
                        const EBCellFAB&       a_phi,
                        const bool&            a_homogeneous,
                        const DataIndex&       a_datInd);
protected:

  static bool                     s_turnOffBCs;
  void defineStencils();
  int                             m_relaxType;
  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;

  RefCountedPtr<EBQuadCFInterp>   m_quadCFIWithCoar;

  EBLevelGrid                     m_eblg;
  EBLevelGrid                     m_eblgFine;
  EBLevelGrid                     m_eblgCoar;
  EBLevelGrid                     m_eblgCoarMG;
  EBLevelGrid                     m_eblgCoarsenedFine;

  RefCountedPtr<ConductivityBaseDomainBC>     m_domainBC;
  RefCountedPtr<ConductivityBaseEBBC>         m_ebBC;

  Real                            m_dxFine;
  Real                            m_dx;
  Real                            m_dxCoar;
  RefCountedPtr<LevelData<EBCellFAB> >          m_acoef;
  RefCountedPtr<LevelData<EBFluxFAB> >          m_bcoef;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >   m_bcoIrreg;

  Real                            m_alpha;
  Real                            m_beta;
  //weights that get multiplied by alpha
  LayoutData<BaseIVFAB<Real> >       m_alphaDiagWeight;
  //weights that get multiplied by beta
  LayoutData<BaseIVFAB<Real> >       m_betaDiagWeight;
  int                             m_refToFine;
  int                             m_refToCoar;
  bool                            m_hasFine;
  bool                            m_hasInterpAve;
  bool                            m_hasCoar;

  //restriction object
  EBMGAverage                    m_ebAverage;
  //prolongation object
  EBMGInterp                     m_ebInterp;

  //stencils for operator evaluation
  LayoutData<RefCountedPtr<EBStencil> >  m_opEBStencil;
  //stencils for operator evaluation on gauss-seidel colors
  LevelData<EBCellFAB>       m_relCoef;

  //cache the vofiterators
  //for irregular cell iteration (includes buffer around multivalued cells)
  LayoutData<VoFIterator >                     m_vofIterIrreg;
  LayoutData<VoFIterator >                     m_vofIterMulti;
  //for domain boundary conditions at ir regular cells
  LayoutData<VoFIterator >                     m_vofIterDomLo[CH_SPACEDIM];
  LayoutData<VoFIterator >                     m_vofIterDomHi[CH_SPACEDIM];

  void calculateRelaxationCoefficient();

  // Coarse-fine stencils for homogeneous CFInterp
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  //flux register with finer level
  EBFastFR       m_fastFR;

  //special mg objects for when we do not have
  //a coarser level or when the refinement to coarse
  //is greater than two
  //flag for when we need special MG objects
  bool                        m_hasMGObjects;
  bool                        m_layoutChanged;
  //stuff below is only defined if m_hasMGObjects==true
  EBMGAverage                 m_ebAverageMG;
  EBMGInterp                  m_ebInterpMG;
  DisjointBoxLayout           m_dblCoarMG;
  EBISLayout                  m_ebislCoarMG;
  ProblemDomain               m_domainCoarMG;

  Vector<IntVect> m_colors;

private:

  void incrementFRCoar(EBFastFR&                   a_fluxReg,
                       const LevelData<EBCellFAB>& a_phiFine,
                       const LevelData<EBCellFAB>& a_phi);

  void incrementFRFine(EBFastFR&                   a_fluxReg,
                       const LevelData<EBCellFAB>& a_phiFine,
                       const LevelData<EBCellFAB>& a_phi,
                       AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  void getFlux(FArrayBox&                    a_flux,
               const FArrayBox&              a_phi,
               const Box&                    a_faceBox,
               const int&                    a_idir,
               const Real&                   a_dx,
               const DataIndex&              a_datInd);



  void applyCFBCs(LevelData<EBCellFAB>&             a_phi,
                  const LevelData<EBCellFAB>* const a_phiCoarse,
                  bool a_homogeneousCFBC);

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

  void applyHomogeneousCFBCs(LevelData<EBCellFAB>&   a_phi);

  void applyHomogeneousCFBCs(EBCellFAB&            a_phi,
                             const DataIndex&      a_datInd,
                             int                   a_idir,
                             Side::LoHiSide        a_hiorlo);
private:
  //copy constructor and operator= disallowed for all the usual reasons
  slowEBCO(const slowEBCO& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const slowEBCO& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};


#include "NamespaceFooter.H"
#endif
