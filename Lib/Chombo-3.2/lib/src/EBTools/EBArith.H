#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBARITH_H_
#define _EBARITH_H_

#include <fstream>

#include "EBCellFAB.H"
#include "EBISBox.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "BoxLayout.H"
#include "BoxLayoutData.H"
#include "Stencils.H"
#include "LayoutData.H"
#include "Tuple.H"
#include "EBFluxFAB.H"
#include "EBCellFAB.H"
#include "NamespaceHeader.H"
///
/**
   enumeration class to ennumerate the different
   extents of norm covereage. sometimes one
   wants to only take a norm over the irregular
   cells.  sometimes one only wants to take a norm
   over the regular cells, and sometimes one wants
   both.
 */
class EBNormType
{
public:

  ///
  /**
     OverOnlyRegular means take the norm only over
     regular cells. \\
     OverOnlyirregular means take the norm only over
     irregular cells. \\
     OverBoth means use all uncovered cells. \\
     All norms exclude covered cells.
   */
  enum  NormMode
  {
    OverOnlyRegular=0,
    OverOnlyIrregular,
    OverBoth
  };
};

///
/**
   class to encapsulate the common arithmatic operations
   for ebfabs
 */
class EBArith
{
public:

  ///
  /**
     Gets the stencil to take the first derivative of  cell centered data.
     no high order one sided stuff (to keep stencil close)
  */
  static int
  getFirstDerivStencilWidthOne(VoFStencil&      a_sten,
                               const VolIndex&  a_vof,
                               const EBISBox&   a_ebisBox,
                               const int&       a_idir,
                               const Real&      a_dx,
                               IntVectSet*    a_cfivsPtr,
                               int ivar);

  ///
  /**
     Least square sten using all cells in a monotone path of radius a_rad
   */
  static void
  getLeastSquaresGradStenAllVoFsRad(VoFStencil&          a_stencil,
                                    Real&                a_weight,
                                    const RealVect&      a_normal,
                                    const RealVect&      a_centroid,
                                    const VolIndex&      a_vof,
                                    const EBISBox&       a_ebisBox,
                                    const RealVect&      a_dx,
                                    const ProblemDomain& a_domain,
                                    int                  a_ivar,
                                    int                  a_rad);

  static int
  getFirstOrderExtrapolationStencil(VoFStencil&     a_stencil,
                                    const RealVect& a_dist,
                                    const RealVect& a_dx,
                                    const VolIndex& a_startVoF,
                                    const EBISBox&  a_ebisBox,
                                    int a_noExtrapThisDir,
                                    IntVectSet*    a_cfivsPtr,
                                    int ivar);
  ///
  /**
   */
  static int orderScript(int icomp, int inorm, int ncomp);

  ///
  /**
   */
  static void
  compareError(Vector<Real>&                            a_orders,
               const Vector< LevelData<EBCellFAB>* >&   a_errorFine,
               const Vector< LevelData<EBCellFAB>* >&   a_errorCoar,
               const Vector< DisjointBoxLayout >&       a_gridsFine,
               const Vector< DisjointBoxLayout >&       a_gridsCoar,
               const Vector< EBISLayout >&              a_ebislFine,
               const Vector< EBISLayout >&              a_ebislCoar,
               const Vector<int>&                       a_refRat,
               const Box&                               a_coarseDom,
               int                                      a_testverbosity,
               fstream*                                 a_fout = NULL,
               Vector<string> names = Vector<string>());

  ///
  static void shrinkIVS(IntVectSet& a_ivs, const int& a_numShrink);

  static void
  timeInterpolate(LevelData<EBCellFAB>&       a_U,
                  const LevelData<EBCellFAB>& a_UOld,
                  const LevelData<EBCellFAB>& a_UNew,
                  const DisjointBoxLayout&    a_grids,
                  const Real&                 a_time,
                  const Real&                 a_told,
                  const Real&                 a_tnew);
  ///
  static Real
  getDiagWeight(  VoFStencil&     a_vofStencil,
                  const VolIndex& a_vof,
                  int a_var = 0);

  ///
  static void getMultiColors(Vector<IntVect>& a_colors);

  ///
  static void
  interpolateCFH(EBCellFAB&                  a_phi,
                 const int  &                a_idir,
                 const Side::LoHiSide&       a_hiorlo,
                 const EBISBox&              a_ebisBox,
                 const Real&                 a_dxfine,
                 const Real&                 a_dxcoar,
                 const IntVectSet&           a_interpIVS);

  ///
  /**
     returns the order of the extrapolation.
     the reason for the last argument is that you might not want the stencil
     to leak over in the noExtrap direction even though you have set a_dist to zero.
     This happens in CF interpolation where you have to really worry about stencil width.
   **/
  static int
  getExtrapolationStencil(VoFStencil&     a_stencil,
                          const RealVect& a_dist,
                          const RealVect& a_dx,
                          const VolIndex& a_startVoF,
                          const EBISBox&  a_ebisBox,
                          int a_noExtrapThisDirection = -1,
                          IntVectSet*    a_cfivsPtr = NULL,
                          int ivar = 0);

  ///
  /**
     returns the order of the extrapolation.
     the reason for the last argument is that you might not want the stencil
     to leak over in the noExtrap direction even though you have set a_dist to zero.
     This happens in CF interpolation where you have to really worry about stencil width.
   **/
  static int
  get1stOrderExtrapolationStencil(VoFStencil&     a_stencil,
                                  const RealVect& a_dist,
                                  const RealVect& a_dx,
                                  const VolIndex& a_startVoF,
                                  const EBISBox&  a_ebisBox,
                                  int a_noExtrapThisDirection = -1,
                                  IntVectSet*    a_cfivsPtr = NULL,
                                  int ivar = 0);

  ///
  /**
   **/
  static void
  getDir1Dir2(int& a_dir1, int& a_dir2, const int& a_dir);

  ///
  /**
     send null ivs if you do not want lohicenter to care about cfivs
   **/
  static void
  loHiCenter(Box&                 a_loBox,
             int&                 a_hasLo,
             Box&                 a_hiBox,
             int&                 a_hasHi,
             Box&                 a_centerBox,
             const ProblemDomain  &   a_eblg,
             const Box&           a_inBox,
             const int&           a_dir,
             IntVectSet*    a_cfivsPtr = NULL);

  ///
  /**
   **/
  static void
  loHi(Box&                 a_loBox,
       int&                 a_hasLo,
       Box&                 a_hiBox,
       int&                 a_hasHi,
       const ProblemDomain& a_eblg,
       const Box&           a_inBox,
       const int&           a_dir);

  ///
  /**
     testRef is the size of the coarsest domain allowed in multigrid. If testRef=2,
     then the coarsest domain in multigrid will be 2x2(x2)
   **/
  static bool getCoarserLayouts(DisjointBoxLayout&       a_dblCoar,
                                ProblemDomain&           a_domainCoar,
                                const DisjointBoxLayout& a_dblFine,
                                const ProblemDomain&     a_domainFine,
                                int                      a_refToCoar,
                                int                      a_maxBoxSize,
                                bool&                    a_layoutChanged,
                                int                      a_testRef = 2);

  ///
  /**
     Gets the stencil to take the first derivative in the given direction
     of  cell centered data.
     Returns the expected order of the derivative.
     When we need them, we prefer
     first derivataves to be second order if at all possible, so we take some pains
     to achieve that
  */
  static int
  getFirstDerivStencil(VoFStencil&      a_sten,
                       const VolIndex&  a_vof,
                       const EBISBox&   a_ebisBox,
                       const int&       a_idir,
                       const Real&      a_dx,
                       IntVectSet*    a_cfivsPtr = NULL,
                       int ivar = 0);


  ///
  /**
     Gets the stencil to take the mixed (in the given directions)
     derivative of  cell centered data.
     Returns the expected order of the derivative.
     When we need them, we usually only need second derivatives to O(h), so
     this just shifts stencil when it has to.
  */
  static int
  getMixedDerivStencil(VoFStencil&      a_sten,
                       const VolIndex&  a_vof,
                       const EBISBox&   a_ebisBox,
                       const int&       a_dir1,
                       const int&       a_dir2,
                       const Real&      a_dx1,
                       const Real&      a_dx2,
                       IntVectSet*    a_cfivsPtr = NULL,
                       int ivar = 0);

  ///
  /**
     Gets the stencil to take the second derivative in the given direction)
     of  cell centered data.
     Returns the expected order of the derivative.
     When we need them, we usually only need second derivatives to O(h), so
     this just shifts stencil when it has to.
  */
  static int
  getSecondDerivStencil(VoFStencil&      a_sten,
                        const VolIndex&  a_vof,
                        const EBISBox&   a_ebisBox,
                        const int&       a_idir,
                        const Real&      a_dx,
                        IntVectSet*    a_cfivsPtr = NULL,
                        int ivar = 0);

  ///
  /**
     a function that could possibly take over the world.
  */
  static void
  getVoFsDir(bool& a_hasClose, VolIndex& a_closeVoF,
             bool& a_hasFar,   VolIndex& a_farVoF,
             const EBISBox& a_ebisBox,
             const VolIndex& a_vof,
             int a_idir, Side::LoHiSide a_sd,
             IntVectSet*    a_cfivsPtr);

  ///
  /**
     define an intvectset of ghost cells which live on
     the coarse fine interface.  includes corner cells.
   */
  static void
  defineCFIVS(LayoutData<IntVectSet>&   a_cfivs,
              const DisjointBoxLayout&  a_grids,
              const ProblemDomain&      a_probDom);

  ///
  static Real
  extrapFaceGradToOutflow(const FaceIndex&      a_bndryFace,
                          const Side::LoHiSide& a_side,
                          const int&            a_idir,
                          const EBGraph&        a_ebGraph,
                          const EBFaceFAB&      a_faceData,
                          const int&            a_comp);

  ///
  static Real
  extrapFaceValueToDomain(const FaceIndex&      a_bndryFace,
                          const Side::LoHiSide& a_side,
                          const int&            a_idir,
                          const EBGraph&        a_ebGraph,
                          const EBFaceFAB&      a_faceData,
                          const int&            a_comp);

  ///
  static Real
  extrapFaceVelToOutflow(const FaceIndex&      a_bndryFace,
                         const Side::LoHiSide& a_side,
                         const int&            a_idir,
                         const EBGraph&        a_ebGraph,
                         const EBFaceFAB&      a_faceData,
                         const int&            a_comp);

  ///
  static Real
  interpolateVel(const EBFaceFAB& a_vel,
                 const EBISBox&   a_ebisBox,
                 const FaceIndex& a_face);

  ///
  static Real
  deInterpolateVel(const Real&      a_centroidVel,
                   const EBFaceFAB& a_vel,
                   const EBISBox&   a_ebisBox,
                   const FaceIndex& a_face);

  ///
  static Real
  getFaceVelForDivFreeCell(const FaceIndex&      a_face,
                           const VolIndex&       a_vof,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const EBFluxFAB&      a_vel,
                           const RealVect&       a_dx,
                           const EBISBox&        a_ebisBox);
  ///
  /**
     evaluates integral(a_divu dV)/integral(dV) over a hierarchy.
     pval = -1--->assumes volfrac already multiplied in
     pval = -2--->multiply by volume fraction
   */
  static void
  meanOverHierarchy(Real&                                   a_mean,
                    const Vector<LevelData<EBCellFAB> *>&   a_divu,
                    const Vector<DisjointBoxLayout>&        a_grids,
                    const Vector<EBISLayout>&               a_ebisl,
                    const Vector<int>&                      a_refRat,
                    const ProblemDomain&                    a_domainCoarsest,
                    const RealVect&                         a_dxCoarsest,
                    const int& a_comp,
                    int pval= -1);

  ///
  /**
   //gets the weight of a vof to interpolate (linear in 2d, bilinear in 3d) to an interpolation point in the plane
   //this is used to get one of the points along the johansen ray
   */
  static Real
  getJohanVoFWeight(const VolIndex& a_whichVoF,     //vof for which we are trying to find weight
                    const RealVect& a_intersectLoc, //location of intersection in real space
                    const IntVect&  a_intersectIV,  //intvect through which ray passes
                    const int&      a_nMaxDir,      //direction of largest normal
                    const RealVect& a_dx);           //grid spacing but you knew that.

  ///version which uses incoming normal, boundary centroid
  /**
     splitting up stuff this way to facillitate multifluid
     which can have multiple normals and boundary centroids per cell.
   */
  static void
  johanStencil(bool&               a_dropOrder,
               Vector<VoFStencil>& a_pointStencils,
               Vector<Real>&       a_distanceAlongLine,
               const RealVect&     a_normal,
               const RealVect&     a_bndryCentroid,
               const VolIndex&     a_vof,
               const EBISBox&      a_ebisBox,
               const RealVect&     a_dx,
               const IntVectSet&   a_cfivs,
               int a_ivar = 0);

  ///
  /**
     Gets the stencils to get the data points along the Johansen ray.   If a_dropOrder
     returns true then the stencil does not exist.  Also returns the distances along the ray
     to each point.
   */
  static void
  johanStencil(bool&               a_dropOrder,
               Vector<VoFStencil>& a_pointStencils,
               Vector<Real>&       a_distanceAlongLine,
               const VolIndex&     a_vof,
               const EBISBox&      a_ebisBox,
               const RealVect&     a_dx,
               const IntVectSet&   a_cfivs,
               int a_ivar = 0);

  ///standard interface
  /**
     gets the normal and calls the other version
   */
  static void
  getLeastSquaresGradSten(VoFStencil&     a_stencil,
                          Real&           a_weight,
                          const VolIndex& a_vof,
                          const EBISBox&  a_ebisBox,
                          const RealVect& a_dx,
                          const ProblemDomain& a_domain,
                          int a_ivar = 0);

  ///version which uses incoming normal, boundary centroid
  /**
     splitting up stuff this way to facillitate multifluid
     which can have multiple normals and boundary centroids per cell.
   */
  static void
  getLeastSquaresGradSten(VoFStencil&     a_stencil,
                          Real&           a_weight,
                          const RealVect& a_normal  ,
                          const RealVect& a_centroid,
                          const VolIndex& a_vof,
                          const EBISBox&  a_ebisBox,
                          const RealVect& a_dx,
                          const ProblemDomain& a_domain,
                          int a_ivar = 0);


  ///version which uses incoming normal, boundary centroid
  /**
     splitting up stuff this way to facillitate multifluid
     which can have multiple normals and boundary centroids per cell.
   */
  static void
  getLeastSquaresGradSten(VoFStencil&     a_stencil,
                          Real&           a_weight,
                          const RealVect& a_normal  ,
                          const RealVect& a_centroid,
                          const IntVect&  a_quadrant,
                          const VolIndex& a_vof,
                          const EBISBox&  a_ebisBox,
                          const RealVect& a_dx,
                          const ProblemDomain& a_domain,
                          int a_ivar = 0);

  ///version which uses incoming normal, boundary centroid, and
  /// all available VoFs.
  /**
     splitting up stuff this way to facillitate multifluid
     which can have multiple normals and boundary centroids per cell.
   */
  static void
  getLeastSquaresGradStenAllVoFs(VoFStencil&          a_stencil,
                                 Real&                a_weight,
                                 const RealVect&      a_normal,
                                 const RealVect&      a_centroid,
                                 const VolIndex&      a_vof,
                                 const EBISBox&       a_ebisBox,
                                 const RealVect&      a_dx,
                                 const ProblemDomain& a_domain,
                                 int                  a_ivar);

  ///version which uses incoming normal, boundary centroid, and
  /// all quadrants.
  /**
     splitting up stuff this way to facillitate multifluid
     which can have multiple normals and boundary centroids per cell.
   */
  static void
  getLeastSquaresGradStenAllQuad(VoFStencil&          a_stencil,
                                 Real&                a_weight,
                                 const RealVect&      a_normal,
                                 const RealVect&      a_centroid,
                                 const VolIndex&      a_vof,
                                 const EBISBox&       a_ebisBox,
                                 const RealVect&      a_dx,
                                 const ProblemDomain& a_domain,
                                 int                  a_ivar,
                                 bool                 a_doSymmetric = false);

  /**
     split matrix computations for least squares into separate
     routine that is independent of stencil size, so that the
     code can be reused
   */
  static void
  calculateWeightingMatrix(RealVect           x0,
                           Vector<RealVect>&  xp,
                           Vector<RealVect>&  weightMatrix,
                           bool&              detZero);

  static void
  calculateWeightingMatrixRed(RealVect           x0,
                              Vector<RealVect>&  xp,
                              IntVect            dimm,
                              Vector<RealVect>&  weightMatrix,
                              bool&              deadRed);

  ///gets the normal of a domain face
  static RealVect  getDomainNormal(int a_idir, Side::LoHiSide a_side);

  ///gets the location in real space of a face center
  static RealVect  getFaceLocation(const FaceIndex& a_face,
                                   const RealVect&  a_dx,
                                   const RealVect&  a_probLo);

  ///get stencils for points along a line (generalized johansen stencil)
  static void dataRayCast(bool&               a_dropOrder,
                          Vector<VoFStencil>& a_pointStencils,
                          Vector<Real>&       a_distanceAlongLine,
                          const RealVect&     a_normal,
                          const RealVect&     a_bndryCentroid,
                          const VolIndex&     a_vof,
                          const EBISBox&      a_ebisBox,
                          const RealVect&     a_dx,
                          const IntVectSet&   a_cfivs,
                          int a_ivar,
                          int a_numPoints);

  ///gets the location in real space of a cell center
  static RealVect  getVofLocation(const VolIndex&  a_vof,
                                  const RealVect&  a_dx,
                                  const RealVect&  a_probLo);

  ///gets the location in real space of a cell center
  static RealVect  getIVLocation(const IntVect&   a_iv,
                                 const RealVect&  a_dx,
                                 const RealVect&  a_probLo);

  ///
  /**
     A very easy to screw up piece of code that was in four places.
   */
  static void
  defineFluxInterpolant(LevelData<BaseIFFAB<Real> >& a_fluxInterpolant,
                        LayoutData<IntVectSet>     & a_irregSetsGrown,
                        const DisjointBoxLayout    & a_dbl,
                        const EBISLayout           & a_ebisl,
                        const ProblemDomain        & a_domain,
                        const int                  & a_ncomp,
                        const int                  & a_faceDir);

  ///
  static void
  interpolateFluxToCentroids(LevelData<EBFluxFAB>&       a_centroidFlux,
                             const LevelData<EBFluxFAB>& a_faceCentFlux,
                             const DisjointBoxLayout&    a_grids,
                             const EBISLayout&           a_ebisl,
                             const ProblemDomain&        a_domain);

  ///
  static void
  interpolateFluxToCentroids(EBFaceFAB&                  a_centroidFlux,
                             const EBFaceFAB&            a_faceCentFlux,
                             const Box&                  a_box,
                             const EBISBox&              a_ebisBox,
                             const ProblemDomain&        a_domain,
                             const int&                  a_idir);

  ///
  static void
  irregNorm(Real& a_ebIrregNorm,
            const BaseIVFAB<Real>& a_ebiError,
            const IntVectSet& a_ivsIrreg,
            const EBISBox& a_ebisBox,
            const int& a_comp,
            const int& a_normtype);

  ///
  static FaceStencil getInterpStencil(const FaceIndex&     a_face,
                                      const IntVectSet&    a_coarseFineIVS,
                                      const EBISBox&       a_ebisBox,
                                      const ProblemDomain& a_domainBox);

  ///
  static void getInterpStencil2D(FaceStencil&          a_sten,
                                 const FaceIndex&      a_face,
                                 const IntVectSet&     a_coarseFineIVS,
                                 const EBISBox&        a_ebisBox,
                                 const ProblemDomain&  a_domainBox) ;

  ///
  static void getInterpStencil3D(FaceStencil&         a_sten,
                                 const FaceIndex&     a_face,
                                 const IntVectSet&    a_coarseFineIVS,
                                 const EBISBox&       a_ebisBox,
                                 const ProblemDomain& a_domainBox);

  ///
  /**
     Return true if there is unique vof associated with a_cell
     that lies in the a_vofsStencil. If so, a_vof = vof. Else return false.
  */
  static bool
  isVoFHere(VolIndex& a_vof,
            const Vector<VolIndex>& a_vofsStencil,
            const IntVect& a_cell );

  ///
  static bool
  isVoFHere(VolIndex& a_vof, int& a_whichVoF,
            const Vector<VolIndex>& a_vofsStencil,
            const IntVect& a_cell );

  ///
  /** Returns all vofs that can be reached from a_vof via a
      monotone path of length <= than radius
  */
  static void
  getAllVoFsInMonotonePath(Vector<VolIndex>& a_vofList,
                           const VolIndex&   a_vof,
                           const EBISBox&    a_ebisBox,
                           const int&        a_redistRad);

  ///get volume of a vof in rz coords
  static void
  getKVolRZ(Real&           a_kvol,
            Real&           a_cellVol,
            const EBISBox&  a_ebisBox,
            const Real&     a_dx,
            const VolIndex& a_vof);

  ///get volume of a vof in rz coords
  static void
  getKVolRZNoDx(Real&           a_kvol,
                Real&           a_cellVolDx3,
                const EBISBox&  a_ebisBox,
                const VolIndex& a_vof);

  ///
  /**
     Return true if there is a unique adjacent face in the
     given direction.  False if either covered or multivalued.
   */
  static
  bool
  getAdjacentFace(FaceIndex& a_adjacentFace,
                  const FaceIndex& a_face,
                  const EBISBox& a_ebisBox,
                  const ProblemDomain& a_domain,
                  const int& a_idir,
                  const Side::LoHiSide& a_side);

  ///
  /**
     compute a stencil for one of our fancy interpolated
     gradients.
   */
  static
  void computeGradFluxStencil(VoFStencil& a_thisStencil,
                              const FaceIndex& a_thisFace,
                              const EBISBox& a_ebisBox,
                              const ProblemDomain& a_domain,
                              const int& a_dir);

  ///
  /**
     compute a stencil for one of our fancy interpolated
     fluxes.
   */
  static
  void computeInterpStencil(FaceStencil& a_thisStencil,
                            const FaceIndex& a_thisFace,
                            const EBISBox& a_ebisBox,
                            const ProblemDomain& a_domain,
                            const int& a_dir);

  ///
  /**
     return lp-norm of component comp of a_src,
     weighted by local volume fraction.
     not normalized by number of points or anything like that.
     if p==0, volume returned is one
     and norm returned is Max(abs(a_src)) over
     uncovered regions.  otherwise,
     returns sum(volfrac*a_src(iv,comp)**p)  of component comp of a_src
     weighted by local volume fraction and
     also returns volume of uncovered regions.
     Only uncovered regions count here.
   */
  static
  void
  volWeightedSum(Real&              a_norm,
                 Real&              a_volume,
                 const EBCellFAB&   a_src,
                 const Box&         a_region,
                 const EBISBox&     a_ebisBox,
                 const RealVect&    a_dx,
                 const IntVectSet&  a_ivsExclude,
                 const int&         a_comp,
                 const int&         a_p,
                 EBNormType::NormMode = EBNormType::OverBoth);


  static void
  getCompVolRZ(Real&           a_compVol,
               const EBISBox&  a_ebisBox,
               const Real&     a_dx,
               const VolIndex& a_vof,
               bool a_verbose = false);
  ///
  /**
     return volume-weighted sum of component comp of a_src,
     weighted by local volume fraction.
     and norm returned is Max(abs(a_src)) over
     uncovered regions.  otherwise,
     returns sum(volfrac*a_src(iv,comp)**p)  of component comp of a_src
     weighted by local volume fraction and
     also returns volume of uncovered regions.
     Only uncovered regions count here.
   */
  static
  void
  computeSum(Real&              a_norm,
             Real&              a_volume,
             const EBCellFAB&   a_src,
             const Box&         a_region,
             const EBISBox&     a_ebisBox,
             const RealVect&    a_dx,
             const IntVectSet&  a_ivsExclude,
             const int&         a_comp,
             EBNormType::NormMode = EBNormType::OverBoth);


  ///
  /**
     return volume-weighted sum of component comp of a_src,
     weighted by local volume fraction.
     and norm returned is Max(abs(a_src)) over
     uncovered regions.  otherwise,
     returns sum(volfrac*a_src(iv,comp)**p)  of component comp of a_src
     weighted by local volume fraction and
     also returns volume of uncovered regions.
     Only uncovered regions count here.
   */
  static
  void
  computeUnweightedSum(Real&             a_norm,
                       Real&             a_volume,
                       const EBCellFAB&  a_src,
                       const Box&        a_region,
                       const EBISBox&    a_ebisBox,
                       const RealVect&   a_dx,
                       const IntVectSet& a_ivsExclude,
                       const int&        a_comp,
                       EBNormType::NormMode = EBNormType::OverBoth);

  ///
  /** return l-p norm of a_src.
      if p==0, v
      norm returned is Max(abs(a_src)) over
      uncovered regions.  otherwise,
      returns 1/vol(sum(volfrac*a_src(iv,comp)**p)^(1/p))
      of component comp of a_src
      weighted by local volume fraction and
      also returns volume of uncovered regions.
      Only uncovered regions count here.
      The data must have the same layout as a_layout with
      the possible exception of ghost cells.

      int pmode = -2;  //norm = (1/v)(sum(phi dv)) ---no absolute values and multiply kappa as you go
      int pmode = -1;  //norm = (1/v)(sum(phi dv)) ---no absolute values and assume kappa already multiplied in

   */
  static
  Real
  norm(const BoxLayoutData<EBCellFAB >& a_dataOne,
       const BoxLayout& a_layout,
       const EBISLayout& a_ebisl,
       const int& comp,
       const int&  p,
       EBNormType::NormMode = EBNormType::OverBoth);

  static
  Real
  norm(const Vector<LevelData<EBCellFAB>* >& a_data,
       const Vector< DisjointBoxLayout > &   a_layout,
       const Vector< EBISLayout >&           a_ebisl,
       const Vector< int >&                  a_refRatio,
       const int& comp,
       const int&  p,
       EBNormType::NormMode = EBNormType::OverBoth);

  ///
  static
  void
  computeCoveredFaces(Vector<VolIndex>&     a_coveredFace,
                      IntVectSet&           a_coveredSets,
                      IntVectSet&           a_irregIVS,
                      const int&            a_idir,
                      const Side::LoHiSide& a_sd,
                      const EBISBox&        a_ebisBox,
                      const Box&            a_region);

  //testing interface
  static
  Real
  norm(Real& volume, const BoxLayoutData<EBCellFAB >& a_dataOne,
       const BoxLayout& a_layout,
       const EBISLayout& a_ebisl,
       const int& comp,
       const int&  p,
       EBNormType::NormMode = EBNormType::OverBoth);

  ///
  /** return l-p norm of a_src.
      if p==0, v
      norm returned is Max(abs(a_src)) over
      uncovered regions.  otherwise,
      returns 1/vol(sum(volfrac*a_src(iv,comp)**p)^(1/p))
      of component comp of a_src
      weighted by local volume fraction and
      also returns volume of uncovered regions.
      Only uncovered regions count here.
      The data must have the same layout as a_layout with
      the possible exception of ghost cells.
   */
  static
  Real
  norm(const EBCellFAB& a_dataOne,
       const Box&       a_grid,
       const EBISBox&   a_ebisl,
       const int&       a_comp,
       const int&       a_p,
       EBNormType::NormMode = EBNormType::OverBoth);

  static
  Real
  norm(Real& volume, const EBCellFAB& a_dataOne,
       const Box&       a_grid,
       const EBISBox&   a_ebisl,
       const int&       a_comp,
       const int&       a_p,
       EBNormType::NormMode = EBNormType::OverBoth);

  ///
  /** return l-p norm of a_src.
      if p==0, v
      norm returned is Max(abs(a_src)) over
      uncovered regions.  otherwise,
      returns sum(volfrac*a_src(iv,comp)**p)  of component comp of a_src
      weighted by local volume fraction and
      also returns volume of uncovered regions.
      Only uncovered regions count here.
      The data must have the same layout as a_layout with
      the possible exception of ghost cells.
   */
  static
  void
  volWeightedSum(Real&                            a_sum,
                 Real&                            a_volume,
                 const BoxLayoutData<EBCellFAB >& a_dataOne,
                 const BoxLayout&                 a_layout,
                 const EBISLayout&                a_ebisl,
                 const RealVect&                  a_dx,
                 const IntVectSet&                a_ivsExclude,
                 const int&                       a_comp,
                 const int&                       a_p,
                 EBNormType::NormMode = EBNormType::OverBoth);

  ///
  /** return the sum of all irregular boundary areas
   */
  static
  void
  sumBndryArea(Real&               a_area,
               const BoxLayout&    a_region,
               const EBISLayout&   a_ebisl);

  ///
  /** return the dotproduct of two leveldatas of ebfabs,
      Only uncovered regions count here.
   */
  static
  Real
  dotProduct(const BoxLayoutData<EBCellFAB >& a_dataOne,
             const BoxLayoutData<EBCellFAB >& a_dataTwo,
             const BoxLayout& a_layout,
             const EBISLayout& a_ebisl,
             const int& a_comp);

  ///
  /** return the dotproduct of two  ebfabs,
      Only uncovered regions count here.
   */
  static
  Real
  dotProduct(const EBCellFAB& a_dataOne,
             const EBCellFAB& a_dataTwo,
             const Box& a_layout,
             const EBISBox& a_ebisBox,
             const int& a_comp);

  ///
  /** Given a VoF, a_vof1, and a cell, a_cell2, determine if there is a vector of
      VoFs in a_cell2 that connects to a_vof1 via a monotone path.
      A vector of VolIndices is returned, even if empty.
   */
  static
  void
  monotonePathVoFToCellMultiVoFs(Vector<VolIndex>& a_vofsInCell,
                                 const Vector<VolIndex>& a_vofsInPath,
                                 const IntVect& a_cell2,
                                 const EBISBox& a_ebisBox);
  ///
  /** Given a VoF, a_vof1, and a cell, a_cell2, determine if there is a single
      VoF in a_cell2 that connects to a_vof1 via a monotone path.  If there is
      one such VoF then TRUE is returned (and the VolIndex is returned,
      a_vof2).  If there is no such VoF or if there are more than one such
      VoF then FALSE is returned (and a_vof2 is unchanged).
   */
  static
  bool
  monotonePathVoFToCellVoF(VolIndex& a_vof2,
                           const VolIndex& a_vof1,
                           const IntVect& a_cell2,
                           const EBISBox& a_ebisBox);
  ///
  /** Returns all vofs that can be reached from a_vof via a
      monotone path of length <= than radius
  */
  static
  void
  getAllVoFsInMonotonePath(Vector<VolIndex>& a_vofs,
                           const IntVect&    a_timesMoved,
                           const IntVect&    a_pathSign,
                           const VolIndex&   a_vof,
                           const EBISBox&    a_ebisBox,
                           const int&        a_radius);

  static VolIndex& getVoFMax()
  {
    return s_vofMax;
  }
  static Real&     getValMax()
  {
    return s_valMax;
  }

  static void setMinVolumeFraction(Real a_minVolFrac)
  {
    s_minVolFrac = a_minVolFrac;
  }

private:
  static VolIndex s_vofMax;
  static Real     s_valMax;
  static Real     s_minVolFrac;

};

extern IntVect ebcoarsen (const IntVect& b,
                          int        refinement_ratio);

extern
Box ebrefine (const Box& b,
              int        refinement_ratio);
extern
ProblemDomain ebrefine (const ProblemDomain& b,
                        int        refinement_ratio);

extern Box ebcoarsen (const Box& b,
                      int        refinement_ratio);
extern ProblemDomain ebcoarsen (const ProblemDomain& b,
                                int        refinement_ratio);

extern void ebrefine(DisjointBoxLayout& output,
                     const DisjointBoxLayout& input,
                     int refinement);

extern void ebcoarsen(DisjointBoxLayout& output,
                      const DisjointBoxLayout& input,
                      int refinement);

#include "NamespaceFooter.H"
#endif
