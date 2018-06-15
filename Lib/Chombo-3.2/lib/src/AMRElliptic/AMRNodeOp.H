#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRNODEOP_H_
#define _AMRNODEOP_H_

#include "AMRMultiGrid.H"
#include "REAL.H"
#include "Box.H"
#include "LevelDataOps.H"
#include "NodeBCFunc.H"
#include "NodeFArrayBox.H"
#include "NodeCFIVS.H"
#include "NodeQCFI.H"
#include "CoarseAverage.H"
#include "LevelFluxRegister.H"
#include "AMRNodeOpF_F.H"
#include "NodeCoarseAverage.H"
#include "NodeLevelDataOps.H"
#include "NamespaceHeader.H"

///
/**
 */
class AMRNodeOp : public AMRLevelOp<LevelData<NodeFArrayBox> >
{
public:
  ///
  /**
   */
  AMRNodeOp();

  ///
  /**
   */
  virtual ~AMRNodeOp()
  {
  }

  /**
     \name LinearOp functions */
  /*@{*/

  ///
  /**
   */
  void define(const DisjointBoxLayout& a_grids,
              const Real&              a_dx,
              const ProblemDomain&     a_domain,
              NodeBCFunc                   a_bc);

  ///
  /**
     define function for AMRLevelOp which has no finer AMR level
   */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_baseBAPtr,
              const Real&              a_dxLevel,
              int                      a_refRatio,
              const ProblemDomain&     a_domain,
              NodeBCFunc                   a_bc);

  ///
  /** full define function for AMRLevelOp with both coarser and finer levels */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_gridsFiner,
              const DisjointBoxLayout& a_gridsCoarser,
              const Real&              a_dxLevel,
              int                      a_refRatio,
              int                      a_refRatioFiner,
              const ProblemDomain&     a_domain,
              NodeBCFunc                   a_bc);

  ///
  /** full define function for AMRLevelOp with finer levels, but no coarser */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout& a_gridsFiner,
              const Real&              a_dxLevel,
              int                      a_refRatio, //  dummy argument, send in 1
              int                      a_refRatioFiner,
              const ProblemDomain&     a_domain,
              NodeBCFunc                   a_bc);

  void projectFineInterior(LevelData<NodeFArrayBox>&       a_phi,
                           const LevelData<NodeFArrayBox>& a_phiFine);

  virtual void residual(  LevelData<NodeFArrayBox>& a_lhs,
                          const LevelData<NodeFArrayBox>& a_phi,
                          const LevelData<NodeFArrayBox>& a_rhs,
                          bool a_homogeneous = false);

  virtual void preCond(   LevelData<NodeFArrayBox>& a_correction,
                          const LevelData<NodeFArrayBox>& a_residual);

  virtual void applyOpOnly(   LevelData<NodeFArrayBox>& a_lhs,
                              const LevelData<NodeFArrayBox>& a_phi);

  virtual void applyOp(   LevelData<NodeFArrayBox>& a_lhs,
                          const LevelData<NodeFArrayBox>& a_phi,
                          bool a_homogeneous = false);
  virtual void create(    LevelData<NodeFArrayBox>& a_lhs,
                          const LevelData<NodeFArrayBox>& a_rhs);
  virtual void createCoarsened(    LevelData<NodeFArrayBox>& a_lhs,
                                   const LevelData<NodeFArrayBox>& a_rhs,
                                   const int& a_refRat);

  virtual void assign(    LevelData<NodeFArrayBox>& a_lhs,
                          const LevelData<NodeFArrayBox>& a_rhs) ;
  virtual Real dotProduct(const LevelData<NodeFArrayBox>& a_1,
                          const LevelData<NodeFArrayBox>& a_2) ;
  virtual void incr( LevelData<NodeFArrayBox>& a_lhs,
                     const LevelData<NodeFArrayBox>& a_x,
                     Real a_scale) ;
  virtual void axby( LevelData<NodeFArrayBox>& a_lhs, const LevelData<NodeFArrayBox>& a_x,
                     const LevelData<NodeFArrayBox>& a_y,
                     Real a, Real b) ;
  virtual void scale(LevelData<NodeFArrayBox>& a_lhs, const Real& a_scale) ;

  virtual Real norm(const LevelData<NodeFArrayBox>& a_x, int a_ord);

  virtual void setToZero( LevelData<NodeFArrayBox>& a_x);

  virtual void relax(LevelData<NodeFArrayBox>& a_e,
                     const LevelData<NodeFArrayBox>& a_residual,
                     int iterations);

  virtual void createCoarser(LevelData<NodeFArrayBox>& a_coarse,
                             const LevelData<NodeFArrayBox>& a_fine,
                             bool ghosted);

  ///
  /**
     calculate restricted residual
     a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
  */
  virtual void restrictResidual(LevelData<NodeFArrayBox>& a_resCoarse,
                                LevelData<NodeFArrayBox>& a_phiFine,
                                const LevelData<NodeFArrayBox>& a_rhsFine);

  ///
  /**
     correct the fine solution based on coarse correction
     a_phiThisLevel += I[2h->h](a_correctCoarse)
  */
  virtual void prolongIncrement(LevelData<NodeFArrayBox>& a_phiThisLevel,
                                const LevelData<NodeFArrayBox>& a_correctCoarse);

  ///
  virtual int refToCoarser()
  {
    return m_refToCoarser;
  }

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
  virtual void AMRResidual(LevelData<NodeFArrayBox>& a_residual,
                           const LevelData<NodeFArrayBox>& a_phiFine,
                           const LevelData<NodeFArrayBox>& a_phi,
                           const LevelData<NodeFArrayBox>& a_phiCoarse,
                           const LevelData<NodeFArrayBox>& a_rhs,
                           bool a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<NodeFArrayBox> >* a_finerOp);

  ///
  /** residual assuming no more coarser AMR levels */
  virtual void AMRResidualNC(LevelData<NodeFArrayBox>& a_residual,
                             const LevelData<NodeFArrayBox>& a_phiFine,
                             const LevelData<NodeFArrayBox>& a_phi,
                             const LevelData<NodeFArrayBox>& a_rhs,
                             bool a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<NodeFArrayBox> >* a_finerOp);

  ///
  /** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
  virtual void AMRResidualNF(LevelData<NodeFArrayBox>& a_residual,
                             const LevelData<NodeFArrayBox>& a_phi,
                             const LevelData<NodeFArrayBox>& a_phiCoarse,
                             const LevelData<NodeFArrayBox>& a_rhs,
                             bool a_homogeneousPhysBC);

  ///
  /** apply AMR operator */
  virtual void AMROperator(LevelData<NodeFArrayBox>& a_LofPhi,
                           const LevelData<NodeFArrayBox>& a_phiFine,
                           const LevelData<NodeFArrayBox>& a_phi,
                           const LevelData<NodeFArrayBox>& a_phiCoarse,
                           bool a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<NodeFArrayBox> >* a_finerOp);

  ///
  /** apply AMR operator, assuming no more coarser AMR levels */
  virtual void AMROperatorNC(LevelData<NodeFArrayBox>& a_LofPhi,
                             const LevelData<NodeFArrayBox>& a_phiFine,
                             const LevelData<NodeFArrayBox>& a_phi,
                             bool a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<NodeFArrayBox> >* a_finerOp);

  ///
  /** AMR operator, assuming no finer AMR levels */
  virtual void AMROperatorNF(LevelData<NodeFArrayBox>& a_LofPhi,
                             const LevelData<NodeFArrayBox>& a_phi,
                             const LevelData<NodeFArrayBox>& a_phiCoarse,
                             bool a_homogeneousPhysBC);

  ///
  /** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
      it is assumed that a_resCoarse has already been filled in with the coarse
      version of AMRResidualNF and that this operation is free to overwrite
      in the overlap regions.
  */
  virtual void AMRRestrict(LevelData<NodeFArrayBox>& a_resCoarse,
                           const LevelData<NodeFArrayBox>& a_residual,
                           const LevelData<NodeFArrayBox>& a_correction,
                           const LevelData<NodeFArrayBox>& a_coarseCorrection,
                            bool a_skip_res = false);

  /** a_correction += I[h->h](a_coarseCorrection) */
  virtual void AMRProlong(LevelData<NodeFArrayBox>& a_correction,
                          const LevelData<NodeFArrayBox>& a_coarseCorrection);

  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void AMRUpdateResidual(LevelData<NodeFArrayBox>& a_residual,
                                 const LevelData<NodeFArrayBox>& a_correction,
                                 const LevelData<NodeFArrayBox>& a_coarseCorrection);

  ///
  /**
     compute norm over all cells on coarse not covered by finer
  */
  virtual Real AMRNorm(const LevelData<NodeFArrayBox>& a_coarseResid,
                       const LevelData<NodeFArrayBox>& a_fineResid,
                       const int&                  a_refRat,
                       const int&                  a_ord);

  /*@}*/

  Real m_alpha, m_beta;

protected:
  //internal--for code reuse
  void prolongIncrement(LevelData<NodeFArrayBox>&       a_phiThisLevel,
                        const LevelData<NodeFArrayBox>& a_correctCoarse,
                        int a_refRat);

  ProblemDomain m_domainInteriorNodes;
  Real                    m_dx;
  Real                    m_dxCrse;
  ProblemDomain           m_domain;
  NodeLevelDataOps        m_levelOps;
  NodeBCFunc              m_bc;
  LayoutData<NodeCFIVS>   m_loCFIVS[SpaceDim];
  LayoutData<NodeCFIVS>   m_hiCFIVS[SpaceDim];
  Copier                  m_exchangeCopier;
  NodeQCFI                m_interpWithCoarser;
  NodeCoarseAverage       m_averageOpMG;
  int                     m_refToCoarser;
  int                     m_refToFiner;
  DisjointBoxLayout       m_coarsenedFineGrids;
  bool                    m_hasFiner;

  void levelGSRB(LevelData<NodeFArrayBox>& a_e,
                 const LevelData<NodeFArrayBox>& a_residual);

  void homogeneousCFInterp(LevelData<NodeFArrayBox>& a_phif);
  void homogeneousCFInterp(LevelData<NodeFArrayBox>& a_phif,
                           const DataIndex& a_datInd,
                           int a_idir,
                           Side::LoHiSide a_hiorlo);

  /** interior boundary nodes
   */
  LayoutData< Vector<IntVectSet> > m_IVSV;

  /** whether each section of m_IVSV is a complete box
   */
  LayoutData< BitSet > m_IVSVfull;

  // exterior boundary nodes
  LayoutData< Vector<IntVectSet> > m_IVSVext;

  /** interior boundary nodes of the coarsened grids at this level
   */
  LayoutData< Vector<IntVectSet> > m_IVSVcoarsened;

  /** interior boundary nodes of the coarsened grids at next finer level
   */
  LayoutData< Vector<IntVectSet> > m_IVSVcoarsenedFine;

  void setCFIVS(const DisjointBoxLayout& a_grids);
};

class AMRNodeOpFactory: public AMRLevelOpFactory<LevelData<NodeFArrayBox> >
{
public:
  virtual ~AMRNodeOpFactory()
  {
  }

  /// full AMR definition function
  void define(const ProblemDomain& a_coarseDomain,
              const Vector<DisjointBoxLayout>& a_grids,
              const Vector<int> a_refRatios,
              const Real&       a_coarsedx,
              NodeBCFunc a_bc,
              Real   a_alpha = 0.,
              Real   a_beta  = 1.);

  /// regular multigrid definition function
  void define(const ProblemDomain& a_domain,
              const DisjointBoxLayout& a_grid,
              const Real&    a_dx,
              NodeBCFunc a_bc,
              int maxDepth = -1,
              Real   a_alpha = 0.,
              Real   a_beta  = 1.);

  /**
     \name MGLevelOpFactory functions */
  /*@{*/
  virtual AMRNodeOp*
  MGnewOp(const ProblemDomain& a_FineindexSpace,
          int depth,
          bool homoOnly = true);

  /*@}*/

  /**
     \name AMRLevelOpFactory functions */
  /*@{*/
  virtual  AMRNodeOp* AMRnewOp(const ProblemDomain& a_indexSpace);

  virtual int refToFiner(const ProblemDomain&) const;

private:
  Vector<ProblemDomain>     m_domains;
  Vector<DisjointBoxLayout> m_boxes;
  Vector<Real>              m_dx;
  Vector<int>               m_refRatios; // refinement to next coarser level
  NodeBCFunc                m_bc;
  Real m_alpha, m_beta;

};

#include "NamespaceFooter.H"
#endif
