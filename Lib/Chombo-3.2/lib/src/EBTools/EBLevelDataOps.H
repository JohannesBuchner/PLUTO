#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBLEVELDATAOPS_H_
#define _EBLEVELDATAOPS_H_

#include "IntVect.H"
#include "LevelData.H"
#include "RefCountedPtr.H"

#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBFluxFAB.H"
#include "EBFluxFactory.H"
#include "EBISLayout.H"
#include "NamespaceHeader.H"

#define EBLEVELDATAOPS_INTERIORREGVOFS   0x1
#define EBLEVELDATAOPS_BOUNDARYREGVOFS   0x2
#define EBLEVELDATAOPS_INTERIORIRREGVOFS 0x4
#define EBLEVELDATAOPS_BOUNDARYIRREGVOFS 0x8

#define EBLEVELDATAOPS_ALLVOFS         (EBLEVELDATAOPS_INTERIORREGVOFS   | \
                                        EBLEVELDATAOPS_BOUNDARYREGVOFS   | \
                                        EBLEVELDATAOPS_INTERIORIRREGVOFS | \
                                        EBLEVELDATAOPS_BOUNDARYIRREGVOFS )

#define EBLEVELDATAOPS_REGULARVOFS     (EBLEVELDATAOPS_INTERIORREGVOFS | \
                                        EBLEVELDATAOPS_BOUNDARYREGVOFS )

#define EBLEVELDATAOPS_IRREGULARVOFS   (EBLEVELDATAOPS_INTERIORIRREGVOFS | \
                                        EBLEVELDATAOPS_BOUNDARYIRREGVOFS )

#define EBLEVELDATAOPS_BOUNDARYVOFS    (EBLEVELDATAOPS_BOUNDARYREGVOFS | \
                                        EBLEVELDATAOPS_IRREGULARVOFS   )

///
/**
   Repository of useful EB functions over a level.  All functions are static so just
   use EBLevelDataOps::blah(...);
 */
class EBLevelDataOps
{
public:
  EBLevelDataOps()
  {
  }

  ~EBLevelDataOps()
  {
  }


  ///
  /**
     If a_interpolateToCentroid = true.  first do cell centered average on grown box.   then interpolate to centroids.
     Otherwise just do cell centered averaging.   This does not do anything special about components (the way mac averaging
     does).   Every component is treated as face centered value = average of adjacent cell values.
  ***/
  static void averageCellToFace(EBFaceFAB           &      a_fluxData,
                                const EBCellFAB     &      a_cellData,
                                const Box           &      a_grid,
                                const EBISBox       &      a_ebisl,
                                const ProblemDomain &      a_domain,
                                int isrc, int idst, int inco,
                                bool a_interpolateToCentroid);


  ///
  /**
     straightforward average from cells to faces over components faceval = half(hival + loval).
     at domain boundaries, just use cell centered value.  nothing crazy
     about ghost cells either.
   */
  static void faceCenteredAverageCellsToFaces(EBFaceFAB           &      a_faceData,
                                              const EBCellFAB     &      a_cellData,
                                              const Box           &      ccFluxBox,
                                              const EBISBox       &      a_ebisBox,
                                              const ProblemDomain &      a_domain,
                                              int isrc, int idst, int inco);
  ///
  /**
     If a_interpolateToCentroid = true.  first do cell centered average on grown box.   then interpolate to centroids.
     Otherwise just do cell centered averaging.   This does not do anything special about components (the way mac averaging
     does).   Every component is treated as face centered value = average of adjacent cell values.
  **/
  static void averageCellToFace(LevelData<EBFluxFAB>&         a_fluxData,
                                const LevelData<EBCellFAB>&   a_cellData,
                                const DisjointBoxLayout&      a_grids,
                                const EBISLayout&             a_ebisl,
                                const ProblemDomain&          a_domain,
                                int isrc, int idst, int inco,
                                bool a_interpolateToCentroid);
  ///
  /**
   */
  static void pruneCoveredBoxes(Vector<Box>&              a_boxes,
                                const ProblemDomain&      a_domain,
                                const EBIndexSpace*       a_ebisPtr);

  ///
  /**
   */
  static int  parallelSum(const int&  a_value);

  ///
  /**
   */
  static long long parallelSum(const long long&  a_value);

  ///
  /**
   */
  static Real parallelSum(const Real& a_value);

  ///
  /**
   */
  static int  parallelMin(const int&  a_value);

  ///
  /**
   */
  static int  parallelMax(const int&  a_value);

  ///
  /**
   */
  static Real parallelMin(const Real& a_value);

  ///
  /**
   */
  static Real parallelMax(const Real& a_value);


  ///
  /**
   */
  static void
  averageCellToFaces(LevelData<EBFluxFAB>&         a_fluxData,
                     const LevelData<EBCellFAB>&   a_cellData,
                     const DisjointBoxLayout&      a_grids,
                     const EBISLayout&             a_ebisl,
                     const ProblemDomain&          a_domain,
                     const int&                    a_comp);


  ///
  /**
   */
  static void
  averageCellToFacesMAC(LevelData<EBFluxFAB>&         a_fluxData,
                        const LevelData<EBCellFAB>&   a_cellData,
                        const DisjointBoxLayout&      a_grids,
                        const EBISLayout&             a_ebisl,
                        const ProblemDomain&          a_domain);


  ///
  /**
   */
  static void exchangeAll(LevelData<EBCellFAB>& a_phi);

  ///
  /**
   */
  static void exchangeCorners(LevelData<EBCellFAB>& a_data,
                              const ProblemDomain&  a_domain);

  ///
  /**
   */
  static void exchangeComp(LevelData<EBCellFAB>& a_data,
                           const int&            a_comp);


  ///
  /**
   */
  static void exchangeAll(LevelData<EBFluxFAB>& a_phi);


  ///
  /**
   */
  static void exchangeComp(LevelData<EBFluxFAB>& a_data,
                           const int&            a_comp);


  ///
  /**
   */
  static void setIrregVal(LevelData<EBCellFAB>&    a_data,
                          const DisjointBoxLayout& a_dbl,
                          const EBISLayout&        a_ebisl,
                          const Real&              a_value);

  ///
  /**
   */
  static bool checkNANINF(const LevelData<EBCellFAB>&a_data,
                          const IntVect&             a_iv1 = IntVect::Zero,
                          const IntVect&             a_iv2 = IntVect::Zero,
                          const Real&                a_shift = 0.0);

  ///
  /**
   */
  static void getMaxMin(Real&                       a_maxVal,
                        Real&                       a_minVal,
                        const LevelData<EBCellFAB>& a_data,
                        const int&                  a_comp,
                        const bool&                 a_doAbs=false);
  ///
  /**
   */
  static void setCoveredVal(LevelData<EBCellFAB>&   a_lData,
                            const Real&             a_value);


  ///
  /**
   */
  static void setCoveredVal(LevelData<EBCellFAB>&   a_lData,
                            const int&              a_comp,
                            const Real&             a_value);


  ///
  /**
   */
  static void setCoveredVal(LevelData<EBFluxFAB>&   a_lData,
                            const Real&             a_value);


  ///
  /**
   */
  static void setCoveredVal(LevelData<EBFluxFAB>&   a_lData,
                            const int&              a_comp,
                            const Real&             a_value);

  ///
  /**
   */
  static void averageMultiVofsToRegFAB(LevelData<EBCellFAB>&    a_data,
                                       const DisjointBoxLayout& a_dbl,
                                       const EBISLayout&        a_ebisl);


  ///
  /**
   */
  static void copyToMultiVofsFromRegFAB(LevelData<EBCellFAB>&    a_data,
                                        const DisjointBoxLayout& a_dbl,
                                        const EBISLayout&        a_ebisl);


  ///
  /**
   */
  static void defineLevelData(LevelData<EBCellFAB>&    a_levelData,
                              const EBISLayout&        a_ebisl,
                              const DisjointBoxLayout& a_dbl,
                              const IntVect&           a_ghosts,
                              const int&               a_nComp);


  ///
  /**
   */
  static void defineLevelData(LevelData<EBFluxFAB>&    a_levelData,
                              const EBISLayout&        a_ebisl,
                              const DisjointBoxLayout& a_dbl,
                              const IntVect&           a_ghosts,
                              const int&               a_nComp);


  ///
  /**
   */
  static void setToZero(LevelData<EBCellFAB>& a_result);


  ///
  /**
   */
  static void setToZero(LevelData<EBFluxFAB>& a_result);


  ///
  /**
   */
  static void scale(LevelData<EBFluxFAB>& a_lhs,
                    const Real&           a_scale) ;

  static void scale(LevelData<EBCellFAB>&                  a_lhs,
                    const  LevelData<EBCellFAB>&           a_scale) ;


  ///
  /**
   */
  static void scale(LevelData<EBCellFAB>& a_lhs,
                    const Real&           a_scale) ;


  ///
  /**
   */
  static void scale(LevelData<EBCellFAB>& a_lhs,
                    const Real& a_scale,
                    const int&  a_comp);


  ///
  /**
   */
  static void setVal(LevelData<EBCellFAB>& a_result,
                     const Real&           a_value);

  ///
  /**
   */
  static void setVal(LevelData<BaseIVFAB<Real> >& a_result,
                     const Real&           a_value);

  ///
  /**
   */
  static void setVal(LevelData<EBCellFAB>& a_result,
                     const Real&           a_value,
                     const int&            a_comp);

  ///
  /**
   */
  static void setVal(LevelData<EBFluxFAB>& a_result,
                     const Real&           a_value);

  ///
  /**
   */
  static void assign(LevelData<EBCellFAB>&       a_to,
                     const LevelData<EBCellFAB>& a_from,
                     const Interval&             a_toInterva,
                     const Interval&             a_fromInterval);


  ///
  /**
   */
  static void assign(LevelData<EBCellFAB>& a_lhs,
                     const LevelData<EBCellFAB>& a_rhs);


  ///
  /**
   */
  static void clone(LevelData<EBCellFAB>& a_lhs,
                    const LevelData<EBCellFAB>& a_rhs);


  ///
  /**
   */
  static void assign(LevelData<EBFluxFAB>& a_lhs,
                     const LevelData<EBFluxFAB>& a_rhs);


  ///
  /**
   */
  static void incr( LevelData<EBCellFAB>& a_lhs,
                    const LevelData<EBCellFAB>&   a_rhs,
                    const Real& a_scale);


  ///
  /**
   */
  static void incr( LevelData<EBCellFAB>& a_lhs,
                    const Real& a_scale);


  ///
  /**
   */
  static void axby( LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    const LevelData<EBCellFAB>& a_y,
                    const Real& a,
                    const Real& b);


  ///
  /**
   */
  static void axby( LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    const LevelData<EBCellFAB>& a_y,
                    const Real& a,
                    const Real& b,
                    const int&  a_lhsComp,
                    const int&  a_xComp,
                    const int&  a_yComp);



  ///
  /**
   */
  static void sum(LevelData<EBCellFAB>&       a_result,
                  const LevelData<EBCellFAB>& a_in1,
                  const LevelData<EBCellFAB>& a_in2);


  ///
  /**
   */
  static void addConstant(LevelData<EBCellFAB>& a_data,
                          const Real&           a_constant);


  ///
  /**
   */
  static void power(LevelData<EBCellFAB>& a_result,
                    const Real&           a_exponent,
                    const int&            a_comp);

  ///
  /**
   */
  static void product(LevelData<EBCellFAB>&       a_result,
                      const LevelData<EBCellFAB>& a_in1,
                      const LevelData<EBCellFAB>& a_in2);


  ///
  /**
   */
  static void product(LevelData<EBFluxFAB>&       a_result,
                      const LevelData<EBFluxFAB>& a_in1,
                      const LevelData<EBFluxFAB>& a_in2);


  ///
  /**
   */
  static void product(LevelData<EBCellFAB>&       a_result,
                      const LevelData<EBCellFAB>& a_in1,
                      const LevelData<EBCellFAB>& a_in2,
                      const int&                  a_rComp,
                      const int&                  a_1Comp,
                      const int&                  a_2Comp);


  ///
  /**
   */
  static void invert(LevelData<EBCellFAB>&       a_result,
                     const LevelData<EBCellFAB>& a_in1);


  ///
  /**
   */
  static  void divideVectorByScalar(LevelData<EBCellFAB>&       a_vectorOut,
                                    const LevelData<EBCellFAB>& a_vectorIn,
                                    const LevelData<EBCellFAB>& a_scalar);

  ///
  /**
   */
  static  void divide(LevelData<EBCellFAB>&       a_result,
                      const LevelData<EBCellFAB>& a_in1,
                      const LevelData<EBCellFAB>& a_in2);


  ///
  /**
   */
  static  void divide(LevelData<EBCellFAB>&       a_result,
                      const LevelData<EBCellFAB>& a_in1,
                      const LevelData<EBCellFAB>& a_in2,
                      const int&                  a_rComp,
                      const int&                  a_1Comp,
                      const int&                  a_2Comp);


  //! Multiply each datum in \a a_data by the largest area fraction
  //!  of the faces of the corresponding cell.
  //! \param a_data The data to be multiplied by the maximum area
  //!  fraction.
  static  void areaFracScalingWeight(LevelData<EBCellFAB>& a_data);

  //! Multiply each datum in \a a_data by the area fraction scaling
  //!  of the corresponding cell.
  //! \param a_data The data to be multiplied by the area fraction
  //!  scaling.
  static  void areaFracScalingWeight(EBCellFAB& a_data);

  //! Multiply each datum in \a a_data by the volume fraction of the
  //! corresponding cell.
  //! \param a_data The data to be multiplied by the volume fraction.
  static  void kappaWeight(LevelData<EBCellFAB>& a_data);

  //! Multiply each datum in \a a_data by the volume fraction of the
  //! corresponding cell.
  //! \param a_data The data to be multiplied by the volume fraction.
  static  void kappaWeight(EBCellFAB& a_data);

  //! Scale each datum in \a a_data by the product of \a a_scale with
  //! the volume fraction of the corresponding cell.
  //! \param a_data The data to be multiplied by the volume fraction.
  //! \param a_scale The scale factor in the product.
  static  void kappaScale(LevelData<EBCellFAB>& a_data,
                          const Real&           a_scale);


  ///
  /**
   */
  static  Real kappaNorm(Real&                       a_volume,
                         const LevelData<EBCellFAB>& a_data,
                         int                         a_which,
                         const ProblemDomain&        a_domain,
                         int                         a_p=2);


  ///
  /**
   */
  static  Vector<Real> vectorKappaNorm(Real&                       a_volume,
                                       const LevelData<EBCellFAB>& a_data,
                                       int                         a_which,
                                       const ProblemDomain&        a_domain,
                                       int                         a_p=2);


  ///
  /**
   */
  static  Real noKappaNorm(Real&                       a_volume,
                           const LevelData<EBCellFAB>& a_data,
                           int                         a_which,
                           const ProblemDomain&        a_domain,
                           int                         a_p=2);


  ///
  /**
   */
  static  Real kappaDotProduct(Real&                       a_volume,
                               const LevelData<EBCellFAB>& a_data1,
                               const LevelData<EBCellFAB>& a_data2,
                               int                         a_which,
                               const ProblemDomain&        a_domain);


  ///
  /**
   */
  static  Real noKappaDotProduct(Real&                       a_volume,
                                 const LevelData<EBCellFAB>& a_data1,
                                 const LevelData<EBCellFAB>& a_data2,
                                 int                         a_which,
                                 const ProblemDomain&        a_domain);


  ///
  /**
   */
  static  Real kappaSumLevel(Real&                       a_volume,
                             const LevelData<EBCellFAB>& a_data,
                             int                         a_which,
                             const ProblemDomain&        a_domain);


  ///
  /**
   */
  static  Real noKappaSumLevel(Real&                       a_volume,
                               const LevelData<EBCellFAB>& a_data,
                               int                         a_which,
                               const ProblemDomain&        a_domain);


  ///
  /**
   */
  static  Vector<Real>   vectorSumKappaPow(Real&                a_volume,
                                           const EBCellFAB&     a_data, const Box& curBox,
                                           int                  a_which,
                                           const ProblemDomain& a_domain,
                                           int                  a_p);

  ///
  /**
   */
  static  Real sumKappaPow(Real&                a_volume,
                           const EBCellFAB&     a_data, const Box& a_region,
                           int                  a_which,
                           const ProblemDomain& a_domain,
                           int                  a_p);


  ///
  /**
   */
  static  Real sumNoKappaPow(Real&                a_volume,
                             const EBCellFAB&     a_data,const Box& a_region,
                             int                  a_which,
                             const ProblemDomain& a_domain,
                             int                  a_p);


  ///
  /**
   */
  static  Real sumKappaDotProduct(Real&                a_volume,
                                  const EBCellFAB&     a_data1,
                                  const EBCellFAB&     a_data2,
                                  const Box& a_region,
                                  int                  a_which,
                                  const ProblemDomain& a_domain);

  ///
  /**
   */
  static  Real sumKappaDotProductAllCells(Real&       a_volume,
                                  const EBCellFAB&     a_data1,
                                  const EBCellFAB&     a_data2,
                                  const Box& a_region,
                                  int                  a_which,
                                  const ProblemDomain& a_domain);
  ///
  /**
   */
  static  Real sumNoKappaDotProduct(Real&                a_volume,
                                    const EBCellFAB&     a_data1,
                                    const EBCellFAB&     a_data2,const Box& a_region,
                                    int                  a_which,
                                    const ProblemDomain& a_domain);


  ///
  /**
   */
  static  Real sumKappa(Real&                a_volume,
                        const EBCellFAB&     a_data,const Box& a_region,
                        int                  a_which,
                        const ProblemDomain& a_domain);


  ///
  /**
   */
  static  Real sumNoKappa(Real&                a_volume,
                          const EBCellFAB&     a_data,const Box& a_region,
                          int                  a_which,
                          const ProblemDomain& a_domain);

  ///

  ///
  /**
   */
  static void
  averageCellToFace(EBFaceFAB&             a_faceData,
                    const EBCellFAB&       a_cellData,
                    const EBGraph&         a_ebGraph,
                    const Box&             a_dblBox,
                    const int&             a_ghostFlux,
                    const int&             a_idir,
                    const ProblemDomain&   a_domain,
                    const int&             a_cellComp,
                    const int&             a_faceComp);


  ///
  /**
   */
  static void gatherBroadCast(Real& a_accum, Real& a_volume, const int& a_p);

  ///
  /**
   */
  static void gatherBroadCast(Vector<Real>& a_accum, Real& a_volume, const int& a_p);


  ///
  /**
   */
  static  Real sum(const LevelData<EBCellFAB> &   a_data,
                   const DisjointBoxLayout &      a_grids,
                   const EBISLayout &             a_ebisl,
                   const IntVectSet &             a_ivsExclude,
                   int   a_comp,
                   bool  a_mutiplyByKappa);
private:

};
#include "NamespaceFooter.H"
#endif
