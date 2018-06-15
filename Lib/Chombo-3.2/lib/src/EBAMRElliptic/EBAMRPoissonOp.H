#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRPOISSONOP_H_
#define _EBAMRPOISSONOP_H_

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
#include "EBFastFR.H"
#include "EBMGAverage.H"
#include "EBMGInterp.H"
#include "PolyGeom.H"
#include "EBQuadCFInterp.H"
#include "EBLevelGrid.H"
#include "AMRTGA.H"
#include "EBAMRIO.H"
#include "AMRPoissonOp.H"
#include "CFRegion.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#endif

#if CH_SPACEDIM==2
#define EBAMRPO_NUMSTEN 4
#elif CH_SPACEDIM==3
#define EBAMRPO_NUMSTEN 8
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
class EBAMRPoissonOp: public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB >
{
public:

  ///
  virtual void dumpAMR(Vector<LevelData<EBCellFAB>*>& a_data, string name)
  {
    writeEBAMRname(&a_data, name.c_str());
  }

  virtual void dumpLevel(LevelData<EBCellFAB>& a_data, string name)
  {
    writeEBLevelname(&a_data, name.c_str());
  }


  ///
  /**
     version that does not fill ebislcoar
   */
  static bool getCoarserLayouts(DisjointBoxLayout&       a_dblCoar,
                                ProblemDomain&           a_domainCoar,
                                const DisjointBoxLayout& a_dblFine,
                                const EBISLayout&        a_ebislFine,
                                const ProblemDomain&     a_domainFine,
                                int                      a_refToCoar,
                                const EBIndexSpace*      a_ebisPtr,
                                int                      a_maxBoxSize,
                                bool&                    a_layoutChanged,
                                int                      a_testRef = 2);

  static Real staticMaxNorm(const LevelData<EBCellFAB>& a_rhs, const EBLevelGrid& a_eblg);

  //for tga to reset stuff
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta);

  //another tgaism
  virtual void diagonalScale(LevelData<EBCellFAB>& a_rhs,
                             bool                  a_kappaWeighted = true);
  virtual void divideByIdentityCoef(LevelData<EBCellFAB> & a_rhs)
  {
    //no acoef here.
  }

  virtual void kappaScale(LevelData<EBCellFAB> & a_rhs)
  {
    //since this is a constant coefficient operator with a=1,
    //this means the same thing as diagonal scale
    diagonalScale(a_rhs);

  }

  ///returns m_dx, such function is required by some LinearSolvers
  Real dx() const
  {
    return m_dx[0];
  }

  ///a leveltgaism
  virtual void fillGrad(const LevelData<EBCellFAB>& a_phi)
  {;}

  ///
  /**
     dump stencil as matrix to stdout
   */
  class StencilIndex
  {
  public:
    StencilIndex(VolIndex& a_vof,
                 int&      a_phase)
      {
        m_vof   = a_vof;
        m_phase = a_phase;
      }

    StencilIndex()
    {
      MayDay::Error("No weak construction of StencilIndex class.");
    }

    StencilIndex& operator= (const StencilIndex& a_sin)
    {
      m_vof   =a_sin.m_vof;
      m_phase = a_sin.m_phase;
      return *this;	
    }

    VolIndex vof() const
    {
      return m_vof;
    }

    int phase() const
    {
      return m_phase;
    }

  protected:
    VolIndex m_vof;
    int      m_phase;
  };

  class StencilIndexComparator
  {
  public:
    bool operator() (const StencilIndex& a_s1, const StencilIndex& a_s2) const
    {
      /*
      int p1 = a_s1.phase();
      int p2 = a_s2.phase();
      if (p1 == p2)
        {
          return (a_s1.vof()<a_s2.vof());
        }
      else
        {
          const IntVect& iv1 = a_s1.vof().gridIndex();
          const IntVect& iv2 = a_s2.vof().gridIndex();
          if (iv1 == iv2)
            {
              return (p1<p2);
            }
          else
            {
              return (iv1<iv2);
            }
        }
      */
      const IntVect& iv1 = a_s1.vof().gridIndex();
      const IntVect& iv2 = a_s2.vof().gridIndex();
      for (int idir=0; idir<SpaceDim; ++idir)
        {
          if (iv1[idir] != iv2[idir])
            {
              return (iv1[idir]<iv2[idir]);
            }
        }
      int p1 = a_s1.phase();
      int p2 = a_s2.phase();
      if (p1 == p2)
        {
          int ci1 = a_s1.vof().cellIndex();
          int ci2 = a_s2.vof().cellIndex();
          return (ci1<ci2);
        }
      else
        {
          return (p1<p2);
        }
    }
  };

  /*
    Calculate stencil weights for each VoF. A mapping of each VolIndex
    an integer is also performed. The first output is this mapping, as
    a matrix that can be read into Matlab/Octave. The second output is
    the stencil weight for each VoF pair, in matrix form.
    The volume fraction is an additional output of the mapping matrix.
    The matrices are output in the format:

    M_(RES)(imap, :) = [ gridIndex[0] gridIndex[1] (gridIndex[2]) alphaWeight kappa];
    where 'RES' is the resolution at the EBPoissonOp's level, 'imap' is
    the integer the VoF has been mapper to, 'alphaWeight' is the
    weighting used to get good conditioning of the matrix, and kappa is
    the VoF's volume fraction.

    L_(RES)(imap, jmap) = value;
    where 'RES' is as above, 'imap' is the integer that the VoF whose
    stencil is being computed has been mapped to, 'jmap' is the integer
    for a VoF in that stencil. 'value' is the stencil weight.
   */
  void dumpStencilMatrix();

  /*
    Calculates the stencil for VoFs at the domain boundary.
   */
  void getDomainFluxStencil(      VoFStencil& a_stencil,
                            const VolIndex&   a_vof,
                            const int         a_comp,
                            const DataIndex&  a_dit);
  /*
    A method for testing EBAMRPoissonOp::dumpStencilMatrix. Since it
    uses calls to 'applyOp' to get the stencil weights, it is slow
    but correct.
   */
  void dumpReferenceStencilMatrix();

  virtual void getFlux(EBFluxFAB&                    a_flux,
                       const LevelData<EBCellFAB>&   a_data,
                       const Box&                    a_grid,
                       const DataIndex&              a_dit,
                       Real                          a_scale)
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        Box ghostedBox = a_grid;
        ghostedBox.grow(1);
        ghostedBox.grow(idir,-1);
        ghostedBox &= m_eblg.getDomain();

        getFlux(a_flux[idir], a_data[a_dit], ghostedBox, a_grid,
                m_eblg.getDomain(),
                m_eblg.getEBISL()[a_dit], m_dx, idir);
      }
  }
  EBLevelGrid getEBLG()
  {
    return m_eblg;
  }
  EBLevelGrid getEBLGCoarMG()
  {
    return m_eblgCoarMG;
  }
  static void setOperatorTime(Real a_time)
  {
    s_time = a_time;
  }

  //conforms to AMRTGA
  virtual void setTime(Real a_time)
  {
    setOperatorTime(a_time);
  }
  ///
  virtual ~EBAMRPoissonOp();

  ///
  EBAMRPoissonOp();

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
     EBAMRPoissonOpFactory to generate these objects.   Really.
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
     a_alpha:         coefficent of identity \\
     a_beta:          coefficient of laplacian.\\
     a_ghostCellsPhi:  Number of ghost cells in phi, correction\\
     a_ghostCellsRhs:  Number of ghost cells in RHS, residual, lphi\\
     Ghost cell arguments are there for caching reasons.  Once you set them, an error is thrown if
     you send in data that does not match.
  */
EBAMRPoissonOp(const EBLevelGrid &                  a_eblgFine,
               const EBLevelGrid &                  a_eblg,
               const EBLevelGrid &                  a_eblgCoar,
               const EBLevelGrid &                  a_eblgCoarMG,
               const RefCountedPtr<EBQuadCFInterp>& a_quadCFI,
               const RefCountedPtr<BaseDomainBC>&   a_domainBC,
               const RefCountedPtr<BaseEBBC>&       a_ebBC,
               const RealVect&                      a_dx,
               const RealVect&                      a_dxCoar,
               const RealVect&                      a_origin,
               const int&                           a_refToFine,
               const int&                           a_refToCoar,
               const bool&                          a_hasFine,
               const bool&                          a_hasCoar,
               const bool&                          a_hasMGObjects,
               const bool&                          a_layoutChanged,
               const int&                           a_numPreCondIters,
               const int&                           a_relaxType,
               const Real&                          a_alpha,
               const Real&                          a_beta,
               const IntVect&                       a_ghostCellsPhi,
               const IntVect&                       a_ghostCellsRHS,
               int                                  a_testRef = 2);

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
     For debugging purpose, get the matrix of the op (could be change from time to time due to EBBC update) by calculating 0 - L(phi(i)) for i=0,1,...,n, each time phi(i)=1 for i-th cell and 0 for all other cells
   */

  virtual void getOpMatrix(const LevelData<EBCellFAB>& a_phi,
                           const LevelData<EBCellFAB>& a_rhs);

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
     get EB flux from leveldata if we have one.  otherwise use m_ebBC
     ignore input data in the case of homogeneous Phys BC
  */
  virtual void
  applyOp(LevelData<EBCellFAB>&                    a_opPhi,
          const LevelData<EBCellFAB>&              a_phi,
          const LevelData<EBCellFAB>* const        a_phiCoar,
          const bool&                              a_homogeneousPhysBC,
          const bool&                              a_homogeneousCFBC,
          const LevelData<BaseIVFAB<Real> >* const a_ebFluxBCLD);

  virtual void
  applyOp(LevelData<EBCellFAB>&                    a_opPhi,
          const LevelData<EBCellFAB>&              a_phi,
          const LevelData<EBCellFAB>* const        a_phiCoar,
          DataIterator&                            a_dit,
          const bool&                              a_homogeneousPhysBC,
          const bool&                              a_homogeneousCFBC,
          const LevelData<BaseIVFAB<Real> >* const a_ebFluxBCLD);

  virtual void
  GSColorAllRegular(BaseFab<Real>&                    a_phi,
                    const BaseFab<Real>&              a_rhs,
                    const int&                        a_icolor,
                    const Real&                       a_weight,
                    const bool&                       a_homogeneousPhysBC,
                    const DataIndex&                  a_dit);

  virtual void
  GSColorAllIrregular(EBCellFAB&                        a_phi,
                      const EBCellFAB&                  a_rhs,
                      const int&                        a_icolor,
                      const bool&                       a_homogeneousPhysBC,
                      const DataIndex&                  a_dit);

  virtual void
  GSColorAllRegularClone(LevelData<EBCellFAB>&       a_phi,
                         const LevelData<EBCellFAB>& a_phiOld,
                         const LevelData<EBCellFAB>& a_rhs,
                         const int&                  a_icolor,
                         const Real&                 a_weight,
                         const bool&                 a_homogeneousPhysBC);

  virtual void
  GSColorAllIrregularClone(LevelData<EBCellFAB>&       a_phi,
                           const LevelData<EBCellFAB>& a_phiOld,
                           const LevelData<EBCellFAB>& a_rhs,
                           const int&                  a_icolor,
                           const bool&                 a_homogeneousPhysBC);

  ///
  /**
     this is the linearop function.  CFBC is set to homogeneous.  phic is null
  */
  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       bool                              a_homogeneousPhysBC);

  /// no exchange of cf interp
  void
  applyOpNoCFBCs(LevelData<EBCellFAB>&                    a_opPhi,
                 const LevelData<EBCellFAB>&              a_phi,
                 const LevelData<EBCellFAB>* const        a_phiCoar,
                 DataIterator&                            a_dit,
                 const bool&                              a_homogeneousPhysBC,
                 const bool&                              a_homogeneousCFBC,
                 const LevelData<BaseIVFAB<Real> >* const a_ebFluxBCLD //only non null in multifluid
                 );

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
  virtual Real localMaxNorm(const LevelData<EBCellFAB>& a_rhs);

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

  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection)
      used in multifluid */
  void AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                         const LevelData<EBCellFAB>& a_correction,
                         const LevelData<EBCellFAB>& a_coarseCorrection,
                         const LevelData<BaseIVFAB<Real> >* const a_ebFluxBCLD);


  void reflux(LevelData<EBCellFAB>& a_residual,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);


  void fast_reflux(LevelData<EBCellFAB>& a_residual,
                  const LevelData<EBCellFAB>& a_phiFine,
                  const LevelData<EBCellFAB>& a_phi,
                  AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  void setEBBC(const RefCountedPtr<BaseEBBC>&      a_ebBC);

  //for debugging other operators
  void slowGSRBColor(LevelData<EBCellFAB>&       a_phi,
                     const LevelData<EBCellFAB>& a_lph,
                     const LevelData<EBCellFAB>& a_rhs,
                     const IntVect&              a_color);

  //for debugging other operators
  void levelSlowRelax(LevelData<EBCellFAB>&       a_phi,
                      const LevelData<EBCellFAB>& a_rhs);

  void levelMultiColorGS(LevelData<EBCellFAB>&       a_phi,
                         const LevelData<EBCellFAB>& a_rhs);

  void levelMultiColorGS(LevelData<EBCellFAB>&       a_phi,
                         const LevelData<EBCellFAB>& a_resid,
                         const IntVect& color);

  void colorGS(LevelData<EBCellFAB>&       a_phi,
               const LevelData<EBCellFAB>& a_rhs,
               const int&                  a_icolor);

  void levelGSRB(LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_rhs);

  void levelGSRB(LevelData<EBCellFAB>&       a_phi,
                 const LevelData<EBCellFAB>& a_rhs,
                 const int                   a_color);

  void levelMultiColorGSClone(LevelData<EBCellFAB>&       a_phi,
                              const LevelData<EBCellFAB>& a_rhs);

  void colorGSClone(LevelData<EBCellFAB>&       a_phi,
                    const LevelData<EBCellFAB>& a_phiOld,
                    const LevelData<EBCellFAB>& a_rhs,
                    const int&                  a_icolor);

  static void doLazyRelax(bool a_doLazyRelax);
  static void doEBEllipticLoadBalance(bool a_doEBEllipticLoadBalance);
  static void areaFracWeighted(bool a_areaFracWeighted);

  static void  getDivFStencil(VoFStencil&      a_vofStencil,
                              const VolIndex&  a_vof,
                              const EBISBox&   a_ebisBox,
                              const IntVectSet& a_cfivs,
                              const RealVect&  a_dx,
                              bool doFaceInterp = true);

  void  getDivFStencil(VoFStencil&      a_vofStencil,
                       const VolIndex&  a_vof,
                       const DataIndex& a_dit,
                       bool doFaceInterp);

  //CP recently added
  void getVoFStencil(LayoutData<BaseIVFAB<VoFStencil> > const*& a_vofStencil)
  {
    //I'm defining const* to make sure the content in m_opStencil is not modified
    a_vofStencil = &m_opStencil;
  }

  void getAlphaDiagWeight(LayoutData<BaseIVFAB<Real> > const*& a_alphaDiagWeight)
  {
    a_alphaDiagWeight = &m_alphaDiagWeight;
  }

  void getAlphaBeta(Real& a_alpha, Real& a_beta)
  {
    a_alpha = m_alpha;
    a_beta  = m_beta;
  }
  //CP recently added
  const RefCountedPtr<BaseDomainBC> getDomainBC()
  {
    //I'm defining const* to make sure the content in m_opStencil is not modified
    return m_domainBC;
  }

  void setListValue(const LevelData<EBCellFAB>& a_data, Real a_value);

  void setListValue(EBCellFAB& a_data, const Vector<VolIndex>& a_setList, Real a_value);

  void setRhsSetList(const LayoutData<Vector<VolIndex> >& a_list);

  LayoutData<Vector<VolIndex> >& getRhsSetList()
  {    
    return m_rhsSetList;
  }

  const LayoutData<Vector<VolIndex> >& getListValue()
  {    
    return m_rhsSetList;
  }

  static void getFluxStencil(VoFStencil&       a_fluxStencil,
                             const FaceIndex&  a_face,
                             const EBISBox&    a_ebisBox,
                             const IntVectSet& a_cfivs,
                             const RealVect&   a_dx,
                             bool a_doFaceInterp);

  static void getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                                         const FaceIndex& a_face,
                                         const RealVect&  a_dx);

  void applyCFBCs(LevelData<EBCellFAB>&             a_phi,
                  const LevelData<EBCellFAB>* const a_phiCoarse,
                  bool a_homogeneousCFBC,
                  bool a_doOnlyRegularInterp = false);

  ///
  /**
     This function computes: a_lhs = (1/diagonal)*a_rhs
     It is used to initialize the preconditioner, and by
     MFPoissonOp::levelJacobi.
     Consider using one of the level Gauss-Seidel methods
     instead of monkeying with this.
  */
  void getInvDiagRHS(LevelData<EBCellFAB>&       a_lhs,
                     const LevelData<EBCellFAB>& a_rhs);

  static int                      s_numComps;
  static int                      s_whichComp;

protected:

  void defineMGObjects(const EBLevelGrid& a_eblgCoarMG);
  void defineWithCoarser(const EBLevelGrid& a_eblgCoar, const int& a_refToCoar);
  void defineWithFiner(const EBLevelGrid& a_eblgFine,
                       const int& a_refToFine);

  static bool                     s_turnOffBCs;
  static bool                     s_doEBEllipticLoadBalance;
  static bool                     s_areaFracWeighted;
  void defineStencils();
  void defineEBCFStencils();

  int                             m_testRef;

  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;

  RefCountedPtr<EBQuadCFInterp>   m_quadCFIWithCoar;

  EBLevelGrid                     m_eblg;
  EBLevelGrid                     m_eblgFine;
  EBLevelGrid                     m_eblgCoar;
  EBLevelGrid                     m_eblgCoarMG;
  EBLevelGrid                     m_eblgCoarsenedFine;

  RefCountedPtr<BaseDomainBC>     m_domainBC;
  RefCountedPtr<BaseEBBC>         m_ebBC;
  LayoutData<Vector<VolIndex> >   m_rhsSetList;

  RealVect                        m_dxFine;
  RealVect                        m_dx;
  RealVect                        m_dxCoar;

  RealVect                        m_invDx;
  RealVect                        m_invDx2;
  Real                            m_dxScale;
  Real                            m_alpha;
  Real                            m_aCoef;
  Real                            m_beta;
  Real                            m_bCoef;
  static Real                     s_time;
  static bool                     s_doLazyRelax;
  static bool                     s_doInconsistentRelax;
  static bool                     s_doTrimEdges;
  static bool                     s_doSetListValueResid;  // if this variable is set to true, it will apply setListValue in residual()
  RealVect                        m_origin;
  int                             m_refToFine;
  int                             m_refToCoar;
  bool                            m_hasFine;
  bool                            m_hasInterpAve;
  bool                            m_hasCoar;
  int                             m_numPreCondIters;
  int                             m_relaxType;

  bool                            m_hasEBCF;

  Copier                          m_exchangeCopier;
  //restriction object
  EBMGAverage                    m_ebAverage;
  //prolongation object
  EBMGInterp                     m_ebInterp;

  //stencils for operator evaluation
  LayoutData<BaseIVFAB<VoFStencil> > m_opStencil;
  LayoutData<RefCountedPtr<EBStencil> >  m_opEBStencil;
  LayoutData<RefCountedPtr<EBStencil> >  m_opEBStencilInhomDomLo[SpaceDim];
  LayoutData<RefCountedPtr<EBStencil> >  m_opEBStencilInhomDomHi[SpaceDim];
  //stencils for operator evaluation on gauss-seidel colors
  LayoutData<RefCountedPtr<EBStencil> >  m_colorEBStencil[EBAMRPO_NUMSTEN];
  LayoutData<RefCountedPtr<EBStencil> >  m_colorEBStencilDomLo[EBAMRPO_NUMSTEN][SpaceDim];
  LayoutData<RefCountedPtr<EBStencil> >  m_colorEBStencilDomHi[EBAMRPO_NUMSTEN][SpaceDim];
  LayoutData<RefCountedPtr<EBStencil> >  m_invDiagEBStencil;
  //weights that get multiplied by alpha
  LayoutData<BaseIVFAB<Real> >       m_alphaDiagWeight;
  //weights that get multiplied by beta
  LayoutData<BaseIVFAB<Real> >       m_betaDiagWeight;
  //constant for using in place of weight in EBStencil->apply for inhom dom bcs
  LayoutData<BaseIVFAB<Real> >       m_one;

  //cache the irreg vofiterator
  LayoutData<VoFIterator >     m_vofItIrreg;
  LayoutData<VoFIterator >     m_vofItIrregColor[EBAMRPO_NUMSTEN];

  LayoutData<VoFIterator >     m_vofItIrregDomLo[SpaceDim];
  LayoutData<VoFIterator >     m_vofItIrregDomHi[SpaceDim];

  LayoutData<VoFIterator >     m_vofItIrregColorDomLo[EBAMRPO_NUMSTEN][SpaceDim];
  LayoutData<VoFIterator >     m_vofItIrregColorDomHi[EBAMRPO_NUMSTEN][SpaceDim];
  Vector<Vector<RefCountedPtr<LayoutData<EBCellFAB> > > > m_cacheInhomDomBCLo;
  Vector<Vector<RefCountedPtr<LayoutData<EBCellFAB> > > > m_cacheInhomDomBCHi;
  // LayoutData<BaseIVFAB<Real> > m_cacheIrregFluxDomLo[EBAMRPO_NUMSTEN][SpaceDim];
  // LayoutData<BaseIVFAB<Real> > m_cacheIrregFluxDomHi[EBAMRPO_NUMSTEN][SpaceDim];

  // Coarse-fine stencils for homogeneous CFInterp
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  //flux register with finer level
  EBFastFR       m_fastFR;

  //stuff to make EBCF go faster
  LayoutData< Vector<FaceIndex>  >   m_faceitCoar[2*SpaceDim];
  LayoutData< Vector<VoFStencil> >  m_stencilCoar[2*SpaceDim];


  //special mg objects for when we do not have
  //a coarser level or when the refinement to coarse
  //is greater than two
  //flag for when we need special MG objects
  bool                        m_hasMGObjects;
  bool                        m_layoutChanged;

  Vector<IntVect> m_colors;


  //stuff below is only defined if m_hasMGObjects==true
  EBMGAverage                 m_ebAverageMG;
  EBMGInterp                  m_ebInterpMG;
  DisjointBoxLayout           m_dblCoarMG;
  EBISLayout                  m_ebislCoarMG;
  ProblemDomain               m_domainCoarMG;

private:
  //internally useful but not for general consumption
  //lots of hidden assumptions and all that

  void fast_incrementFRCoar(const LevelData<EBCellFAB>& a_phiFine,
                            const LevelData<EBCellFAB>& a_phi);

  void fast_incrementFRFine(const LevelData<EBCellFAB>& a_phiFine,
                            const LevelData<EBCellFAB>& a_phi,
                            AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  ///this one gets called by base level tga (called by the public getFlux)
  void getFlux(EBFaceFAB&                    a_flux,
               const EBCellFAB&              a_phi,
               const Box&                    a_ghostedBox,
               const Box&                    a_fabBox,
               const ProblemDomain&          a_domainBox,
               const EBISBox&                a_ebisBox,
               const RealVect&               a_dx,
               const int&                    a_idir);

  ///this one is internal (called by refluxing)
  void getFluxEBCF(EBFaceFAB&                    a_flux,
                   const EBCellFAB&              a_phi,
                   const Box&                    a_ghostedBox,
                   Vector<FaceIndex>&            a_faceitEBCF,
                   Vector<VoFStencil>&           a_ebcfsten,
                   const RealVect&               a_dx);

  ///this one is internal (called by refluxing)
  void getFluxRegO(EBFaceFAB&                    a_flux,
                  const EBCellFAB&              a_phi,
                  const Box&                    a_ghostedBox,
                  const RealVect&               a_dx);

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

  Real getRadius(const FaceIndex& a_face, const RealVect& a_centroid);

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

private:
  //copy constructor and operator= disallowed for all the usual reasons
  EBAMRPoissonOp(const EBAMRPoissonOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBAMRPoissonOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};


#include "NamespaceFooter.H"
#endif
