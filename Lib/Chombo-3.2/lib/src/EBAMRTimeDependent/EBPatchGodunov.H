#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHGODUNOV_H_
#define _EBPATCHGODUNOV_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "CH_HDF5.H"
#include "EBCellFAB.H"
#include "EBFluxFAB.H"
#include "EBISBox.H"
#include "BaseIVFAB.H"
#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"
#include "Stencils.H"
#include "AggStencil.H"
#include "NamespaceHeader.H"

///
/**
 */
class EBPatchGodunov
{
public:
  ///
  /**
   */
  EBPatchGodunov();

  ///non-virtual stuff
  void
  incrementWithSource(EBCellFAB&       a_primState,
                      const EBCellFAB& a_source,
                      const Real&      a_scale,
                      const Box&       a_box) ;

  ///
  virtual void
  extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                       const EBCellFAB&        a_primMinu,
                       const EBCellFAB&        a_primPlus,
                       const EBCellFAB&        a_primState,
                       const Vector<VolIndex>& a_coveredFaces,
                       const int&              a_faceDir,
                       const Side::LoHiSide&   a_sd,
                       const Box&              a_box);

  ///
  virtual void
  pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                         const EBCellFAB&        a_primMinu,
                         const EBCellFAB&        a_primPlus,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);

  ///
  virtual void
  pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                         const EBCellFAB&        a_primMinu,
                         const EBCellFAB&        a_primPlus,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);

  ///
  void
  computeCoveredFaces(Vector<VolIndex>&     a_coveredFace,
                      IntVectSet&           a_coveredIVS,
                      IntVectSet&           a_irregIVS,
                      const int&            a_idir,
                      const Side::LoHiSide& a_sd,
                      const Box&            a_region);
  ///
  /**
   */
  const EBPhysIBC*  getEBPhysIBC() const;

  ///
  /**
   */
  Real
  pointLimiter(const Real& a_deltaW1,
               const Real& a_deltaW2);

  ///
  /**
   */
  void
  setEBPhysIBC(const EBPhysIBCFactory& a_bc);

  /// Set parameters for slope computations
  /**
   */
  void setSlopeParameters(bool a_fourthOrderSlopes,
                          bool a_flattening,
                          bool a_useLimiting);

  ///
  /**
     Given left and right one-sided undivided differences /a_dql,a_dqr/,
     apply van Leer limiter $vL$ defined in section
     to a_dq.
     Called by the default implementation of  PatchPolytropic::slope.
  */
  void
  applyLimiter(EBCellFAB&       a_slopePrim,
               const EBCellFAB& a_slopePrimLeft,
               const EBCellFAB& a_slopePrimRigh,
               const int&       a_dir,
               const Box&       a_box);

  virtual void
  doSecondOrderSlopes(EBCellFAB&       a_delta2W,
                      EBCellFAB&       a_deltaWL,
                      EBCellFAB&       a_deltaWR,
                      EBCellFAB&       a_deltaWC,
                      const EBCellFAB& a_primState,
                      const int&       a_dir,
                      const Box&       a_box,
                      bool a_doAggregated = false);

  virtual void
  aggIrregSecondOrderSlopes(EBCellFAB&       a_delta2W,
                            EBCellFAB&       a_deltaWL,
                            EBCellFAB&       a_deltaWR,
                            EBCellFAB&       a_deltaWC,
                            const EBCellFAB& a_primState,
                            const int&       a_dir,
                            const Box&       a_box);


  virtual void
  irregSecondOrderSlopes(EBCellFAB&       a_delta2W,
                         EBCellFAB&       a_deltaWL,
                         EBCellFAB&       a_deltaWR,
                         EBCellFAB&       a_deltaWC,
                         const EBCellFAB& a_primState,
                         const int&       a_dir,
                         const Box&       a_box);

  void
  doFourthOrderSlopes(EBCellFAB& a_delta4W,
                      EBCellFAB& a_deltaWC,
                      const EBCellFAB& a_delta2W,
                      const EBCellFAB& a_deltaWL,
                      const EBCellFAB& a_deltaWR,
                      const EBCellFAB& a_primState,
                      const int& a_dir,
                      const Box& a_box);
  void
  pointGetSlopes(Real&            a_dql,
                 Real&            a_dqr,
                 Real&            a_dqc,
                 bool&            a_hasFacesLeft,
                 bool&            a_hasFacesRigh,
                 const VolIndex&  a_vof,
                 const EBCellFAB& a_primState,
                 const int&       a_dir,
                 const int&       a_ivar,
                 const bool&      a_verbose);

  virtual void
  coveredExtrapSlopes(Real&            a_dqc,
                      const VolIndex&  a_vof,
                      const EBCellFAB& a_primState,
                      const int&       a_dir,
                      const int&       a_ivar);

  ///
  /**
   */
  void
  computeFluxes(EBFluxFAB&          a_flux,
                BaseIVFAB<Real>     a_coveredFluxLo[SpaceDim],
                BaseIVFAB<Real>     a_coveredFluxHi[SpaceDim],
                Vector<VolIndex>    a_coveredFaceLo[SpaceDim],
                Vector<VolIndex>    a_coveredFaceHi[SpaceDim],
                EBCellFAB&          a_primState,
                EBCellFAB           a_slopesPrim[SpaceDim],
                EBCellFAB           a_slopesSeco[SpaceDim],
                const EBCellFAB&    a_flattening,
                const EBCellFAB&    a_consState,
                const EBCellFAB&    a_source,
                const Box&          a_box,
                const DataIndex&    a_dit,
                bool                a_verbose);

  ///
  /**
   */
virtual  void
  extrapolatePrim2D(EBCellFAB           a_primMinu[SpaceDim],
                    EBCellFAB           a_primPlus[SpaceDim],
                    EBCellFAB&          a_primState,
                    EBCellFAB           a_slopesPrim[SpaceDim],
                    EBCellFAB           a_slopesSeco[SpaceDim],
                    const EBCellFAB&    a_flattening,
                    const EBCellFAB&    a_consState,
                    const EBCellFAB&    a_source,
                    const Box&          a_box,
                    const DataIndex&    a_dit,
                    bool                a_verbose);

  ///
  /**
   */
  void
  extrapolatePrim3D(EBCellFAB           a_primMinu[SpaceDim],
                    EBCellFAB           a_primPlus[SpaceDim],
                    EBCellFAB&          a_primState,
                    EBCellFAB           a_slopesPrim[SpaceDim],
                    EBCellFAB           a_slopesSeco[SpaceDim],
                    const EBCellFAB&    a_flattening,
                    const EBCellFAB&    a_consState,
                    const EBCellFAB&    a_source,
                    const Box&          a_box,
                    const DataIndex&    a_dit,
                    bool                a_verbose);

  ///
  /**
     Update the state using flux difference that ignores EB.
     Store fluxes used in this update.
     Store non-conservative divergence.
     Flux coming out of htis this should exist at cell face centers.
   */
  virtual void
  regularUpdate(EBCellFAB&          a_consState,
                EBFluxFAB&          a_flux,
                BaseIVFAB<Real>&    a_ebIrregFlux,
                BaseIVFAB<Real>&    a_nonConservativeDivergence,
                const EBCellFAB&    a_flattening,
                const EBCellFAB&    a_source,
                const Box&          a_box,
                const IntVectSet&   a_ivs,
                const DataIndex&    a_dit,
                bool                a_verbose);

  virtual void
  regularDivergences(EBCellFAB&          a_nonconsdiv,
                     EBCellFAB&          a_consState,
                     EBFluxFAB&          a_flux,
                     BaseIVFAB<Real>&    a_ebIrregFlux,
                     BaseIVFAB<Real>&    a_nonConservativeDivergence,
                     const EBCellFAB&    a_flattening,
                     const EBCellFAB&    a_source,
                     const Box&          a_box,
                     const IntVectSet&   a_ivs,
                     const DataIndex&    a_dit,
                     bool                a_verbose);

  virtual void
  hybridDivergence(EBCellFAB&             a_hybridDiv,
                   EBCellFAB&             a_consState,
                   BaseIVFAB<Real>&       a_massDiff,
                   const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                   const BaseIVFAB<Real>& a_ebIrregFlux,
                   const BaseIVFAB<Real>& a_nonConservativeDivergence,
                   const Box&             a_box,
                   const IntVectSet&      a_ivs);
  ///
  /**
     Update the state at irregular VoFs and compute mass difference and
     the maximum wave speed over the entire box. Flux going into this
     should exist at VoF centroids.
   */
  virtual void
  irregularUpdate(EBCellFAB&             a_consState,
                  Real&                  a_maxWaveSpeed,
                  BaseIVFAB<Real>&       a_massDiff,
                  const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                  const BaseIVFAB<Real>& a_ebIrregFlux,
                  const BaseIVFAB<Real>& a_nonConservativeDivergence,
                  const Box&             a_box,
                  const IntVectSet&      a_ivs);

  ///
  /**
   */
  virtual void
  doNormalDerivativeExtr2D(EBCellFAB              a_primMinu[SpaceDim],
                           EBCellFAB              a_primPlus[SpaceDim],
                           EBFaceFAB              a_fluxOne[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                           EBCellFAB              a_slopesPrim[SpaceDim],
                           EBCellFAB              a_slopesSeco[SpaceDim],
                           const EBCellFAB&       a_flattening,
                           const EBCellFAB&       a_primState,
                           const EBCellFAB&       a_source,
                           const DataIndex&       a_dit,
                           const Box&             a_box);

  ///
  /**
   */
  void
  doNormalDerivativeExtr3D(EBCellFAB              a_primMinu[SpaceDim],
                           EBCellFAB              a_primPlus[SpaceDim],
                           EBFaceFAB              a_fluxOne[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormMinu[SpaceDim],
                           BaseIVFAB<Real>        a_coveredFluxNormPlus[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormMinu[SpaceDim],
                           Vector<VolIndex>       a_coveredFaceNormPlus[SpaceDim],
                           EBCellFAB              a_slopesPrim[SpaceDim],
                           EBCellFAB              a_slopesSeco[SpaceDim],
                           const EBCellFAB&       a_flattening,
                           const EBCellFAB&       a_primState,
                           const EBCellFAB&       a_source,
                           const DataIndex&       a_dit,
                           const Box&             a_box);

  ///
  /**
   */
  virtual void
  finalExtrap2D(EBCellFAB              a_primMinu[SpaceDim],
                EBCellFAB              a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
                const EBFaceFAB        a_fluxOne[SpaceDim],
                const EBCellFAB&       a_primState,
                const EBCellFAB        a_slopesPrim[SpaceDim],
                const EBCellFAB        a_slopesSeco[SpaceDim],
                const Box&             a_box);

  ///
  /**
   */
  void
  finalExtrap3D(EBCellFAB              a_primMinu[SpaceDim],
                EBCellFAB              a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFlux3DMinu[SpaceDim][SpaceDim],
                const BaseIVFAB<Real>  a_coveredFlux3DPlus[SpaceDim][SpaceDim],
                const EBFaceFAB        a_fluxTwoVec[SpaceDim][SpaceDim],
                const EBCellFAB&       a_primState,
                const EBCellFAB        a_slopesPrim[SpaceDim],
                const EBCellFAB        a_slopesSeco[SpaceDim],
                const Box&             a_box);

  ///
  /**
   */
  void
  do111coupling(EBFaceFAB              a_fluxTwoVec[SpaceDim][SpaceDim],
                BaseIVFAB<Real>        a_coveredFlux3DMinu[SpaceDim][SpaceDim],
                BaseIVFAB<Real>        a_coveredFlux3DPlus[SpaceDim][SpaceDim],
                const EBCellFAB        a_primMinu[SpaceDim],
                const EBCellFAB        a_primPlus[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormMinu[SpaceDim],
                const BaseIVFAB<Real>  a_coveredFluxNormPlus[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormMinu[SpaceDim],
                const Vector<VolIndex> a_coveredFaceNormPlus[SpaceDim],
                const EBFaceFAB        a_fluxOne[SpaceDim],
                const EBCellFAB&       a_primState,
                const EBCellFAB        a_slopesPrim[SpaceDim],
                const EBCellFAB        a_slopesSeco[SpaceDim],
                const DataIndex&       a_dit,
                const Box&             a_box);

  ///
  /**
   */
  virtual void
  interpolateFluxToCentroids(BaseIFFAB<Real>              a_centroidFlux[SpaceDim],
                             const BaseIFFAB<Real>* const a_fluxInterpolant[SpaceDim],
                             const IntVectSet&            a_irregIVS);

  ///
  /**
   */
  virtual void
  finalUpdate(EBCellFAB&             a_consState,
              BaseIVFAB<Real>&       a_massDiff,
              const BaseIVFAB<Real>& a_nonConsDivF,
              const BaseIVFAB<Real>& a_conservDivF,
              const IntVectSet&      a_ivs);

  ///
  virtual void
  getFaceDivergence(EBFluxFAB&        a_openDivU,
                    const EBCellFAB&  a_primState,
                    const EBCellFAB   a_slopePrim[SpaceDim],
                    const Box&        a_box,
                    const IntVectSet& a_ivsIrreg);

  ///
  void artificialViscosity(bool a_useArtificialVisc);

  ///
  void
  applyArtificialViscosity(EBFluxFAB&        a_openFlux,
                           BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                           BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                           const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                           const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                           const EBCellFAB&  a_consState,
                           const EBFluxFAB&  a_divVel,
                           const Box&        a_box,
                           const IntVectSet& a_ivsIrreg);

  Real bilinearFunc(const Real  a_WVal[2][2],
                    const Real& a_xd1,
                    const Real& a_xd2);

  Real maxFunction(const Real  a_WVal[2][2],
                    const Real& a_xd1,
                    const Real& a_xd2);

  ///
  /**
   */
  void
  define(const ProblemDomain& a_domain,
         const RealVect& a_dx, bool a_useAgg = false);


  ///deprecated interface
  void
  define(const ProblemDomain& a_domain,
         const Real& a_dx)
  {
    m_useAgg = false;
    RealVect dxVect = a_dx*RealVect::Unit;
    define(a_domain, dxVect);
  }

  ///
  /**
     needs coarse-fine IVS to know where to drop order for interpolation
     virtual in case you need to add anything to definition
   */
  virtual void
  setValidBox(const Box& a_validBox,
              const EBISBox& a_ebisBox,
              const IntVectSet& a_coarseFineIVS,
              const Real& a_time,
              const Real& a_dt);

  void
  computeFlattening(EBCellFAB&       a_flattening,
                    const EBCellFAB& a_primState,
                    const Box&       a_box);

  ///
  /**
     Compute the limited slope /a_dq/ of the primitive variables
     /a_q/ for the components in the interval /a_interval/,
     Calls user-supplied  EBPatchGodunov::applyLimiter.
  */
  virtual void
  slope(EBCellFAB&       a_slopePrim,
        EBCellFAB&       a_slopeNLim,
        const EBCellFAB& a_primState,
        const EBCellFAB& a_flattening,
        const int&       a_dir,
        const Box&       a_box,
        bool a_doAggregated = false);

  ///needs to be virtual because of RZ
  virtual void
  nonconservativeDivergence(EBCellFAB&             a_divF,
                            const EBFluxFAB&       a_flux,
                            const BaseIVFAB<Real>  a_coveredFluxMinu[SpaceDim],
                            const BaseIVFAB<Real>  a_coveredFluxPlus[SpaceDim],
                            const Vector<VolIndex> a_coveredFaceMinu[SpaceDim],
                            const Vector<VolIndex> a_coveredFacePlus[SpaceDim],
                            const Box&             a_box);

  ///virtual in case you want to do something faster than go through constoprim
  virtual void
  updatePrim(EBCellFAB&              a_primMinu,
             EBCellFAB&              a_primPlus,
             const EBFaceFAB&        a_flux,
             const BaseIVFAB<Real>&  a_coveredFluxMinu,
             const BaseIVFAB<Real>&  a_coveredFluxPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_dir,
             const Box&              a_box,
             const Real&             a_scale) ;

  ///virtual because RZ changes this function
  virtual void
  updateCons(EBCellFAB&              a_consState,
             const EBFaceFAB&        a_flux,
             const BaseIVFAB<Real>&  a_coveredFluxMinu,
             const BaseIVFAB<Real>&  a_coveredFluxPlus,
             const Vector<VolIndex>& a_coveredFaceMinu,
             const Vector<VolIndex>& a_coveredFacePlus,
             const int&              a_dir,
             const Box&              a_box,
             const Real&             a_scale);
  ///
  /**
     Returns the interval of component indices in the primitive variable
     EBCellFAB for the velocities.   Only used for artificial visc and flattening
  */
  virtual Interval velocityInterval() const = 0;

  ///needs to be virtual because of RZ
  /**
   */
  virtual void
  consUndividedDivergence(BaseIVFAB<Real>&       a_divF,
                          const BaseIFFAB<Real>  a_centroidFlux[SpaceDim],
                          const BaseIVFAB<Real>& a_ebIrregFlux,
                          const IntVectSet&      a_ivs);

  ///
  /**
   */
  virtual void
  computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                     const EBCellFAB&  a_primState,
                     const EBCellFAB   a_slopePrim[SpaceDim],
                     const IntVectSet& a_irregIVS,
                     const EBCellFAB&  a_source) = 0;

  ///
  /**
     Returns the component index for the pressure. Called only if flattening is used.
  */
  virtual int pressureIndex() const = 0;

  ///
  /**
     Returns the component index for the bulk modulus, used as a
     normalization to measure shock strength in flattening.
     Called only if flattening is used.
  */
  virtual int bulkModulusIndex() const = 0;

  ///
  /**
   */
  virtual ~EBPatchGodunov();

  ///
  /**
   */
  virtual void
  setCoveredConsVals(EBCellFAB& a_consState)= 0;

  ///
  /**
   */
  virtual Real
  getMaxWaveSpeed(const EBCellFAB& a_consState,
                  const Box& a_box) = 0;

  ///
  virtual void
  floorConserved(EBCellFAB& a_consState,
                 const Box& a_box) = 0;

  ///
  virtual void
  floorPrimitives(EBCellFAB& a_primState,
                  const Box& a_box) = 0;

  ///
  virtual void
  floorConserved(BaseIVFAB<Real>&  a_consState,
                 const IntVectSet& a_ivs) = 0;

  ///
  virtual void
  floorPrimitives(BaseIVFAB<Real>&  a_primState,
                  const IntVectSet& a_ivs) = 0;

  ///
  virtual void  getCoveredValuesPrim(Vector<Real>& a_covValues)=0;

  ///
  virtual void  getCoveredValuesCons(Vector<Real>& a_covValues)=0;

  virtual int densityIndex() const = 0;

  ///
  /**
     Return number of components for primitive variables.
  */
  virtual int numPrimitives() const = 0;

  ///
  /**
     Returns number of components for flux variables.
  */
  virtual int numFluxes() const = 0;

  ///
  /**
     Returns number of components for flux variables.
  */
  virtual int numSlopes() const = 0;

  ///
  /**
     Returns number of components for conserved variables.
  */
  virtual int numConserved() const = 0;

  ///
  /**
     Return the names of the variables.  A default
     implementation is provided that puts in generic names.
  */
  virtual Vector<string> stateNames()=0;

  ///
  /**
     Return the names of the variables.  A default
     implementation is provided that puts in generic names.
  */
  virtual Vector<string> primNames()=0;

  virtual void
  normalPred(EBCellFAB&       a_primLo,
             EBCellFAB&       a_primHi,
             const EBCellFAB& a_primState,
             const EBCellFAB& a_slopePrim,
             const Real&      a_scale,
             const int&       a_dir,
             const Box&       a_box) = 0;

  ///
  /**
     Given input left and right states, compute a suitably-upwinded
     flux (e.g. by solving a Riemann problem), as in
  */
  virtual void
  riemann(EBFaceFAB&       a_flux,
          const EBCellFAB& a_primLeft,
          const EBCellFAB& a_primRight,
          const int&       a_dir,
          const Box&       a_box) = 0;

  ///
  /**
     Given input left and right states, compute a suitably-upwinded
     flux (e.g. by solving a Riemann problem).
  */
  virtual  void
  riemann(BaseIVFAB<Real>&        a_coveredFlux,
          const BaseIVFAB<Real>&  a_extendedState,
          const EBCellFAB&        a_primState,
          const Vector<VolIndex>& a_region,
          const int&              a_dir,
          const Side::LoHiSide&   a_sd,
          const Box& a_box)= 0;

  ///rz func.
  virtual void
  setSource(EBCellFAB&       a_source,
              const EBCellFAB& a_consState,
              const Box&       a_box)
  {;}

  ///rz func.
  virtual void
  assembleFluxReg(EBFaceFAB&       a_fluxRegFlux,
                  const EBFaceFAB& a_godunovFlux,
                  const int&       a_idir,
                  const Box&       a_cellBox)
  {;}

  ///rz func.
  virtual void
  assembleFluxIrr(BaseIFFAB<Real>&       a_fluxRegFlux,
                  const BaseIFFAB<Real>& a_godunovFlux,
                  const int&             a_idir,
                  const Box&             a_cellBox,
                  const IntVectSet&      a_set)
  {;}

  ///
  virtual void
  consToPrim(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box,
             int              a_logflag,
             bool             a_verbose = false)=0;

#ifdef CH_USE_HDF5
  virtual void
  expressions(HDF5HeaderData& a_holder)
  {;}
#endif

  ///
  virtual void
  consToPrim(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs) =0;

  ///
  /**
   */
  virtual void
  primToCons(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box) = 0;

  ///
  /**
   */
  virtual void
  primToCons(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs) = 0;

  ///
  /**
     Return true if the application is using flattening.
  */
  virtual bool usesFlattening() const = 0;

  ///
  /**
     Return true if the application is using artificial viscosity.
  */
  virtual bool usesArtificialViscosity() const = 0;

  ///
  /**
     Return true if you are using fourth-order slopes.
     Return false if you are using second-order slopes.
  */
  virtual bool usesFourthOrderSlopes() const = 0;

  ///
  /**
     Returns value of artificial viscosity. Called only if artificial
     viscosity is being used.
  */
  virtual Real artificialViscosityCoefficient() const  = 0;

  ///
  /**
   */
  virtual bool isDefined() const;

  //debugging hooks.  not for the faint of heart
  static void setVerbose(bool a_verbose);
  static void setCurLevel(int a_curLevel);
  static void setCurComp(int a_curComp);
  static void setDoingVel(int a_yesorno);
  static void setDoingAdvVel(int a_yesorno);
  static int getDoingAdvVel();
  static int getDoingVel();
  static int getCurLevel();
  static int getCurComp();
  static bool getVerbose();
  static Real getMaxWaveSpeed();
  static void setMaxWaveSpeed(Real a_maxWaveSpeedIV);
  static IntVect getMaxWaveSpeedIV();
  static void    setMaxWaveSpeedIV(const IntVect& a_maxWaveSpeedIV);

  EBCellFAB& getPrimState()
  {
    return m_primState;
  }
  BaseIVFAB<Real>* getCoveredFluxPlus()
  {
    return m_coveredFluxPlusG4;
  }
  BaseIVFAB<Real>* getCoveredFluxMinu()
  {
    return m_coveredFluxMinuG4;
  }

  Vector<VolIndex>* getCoveredFaceMinu()
  {
    return m_coveredFaceMinuG4;
  }

  Vector<VolIndex>* getCoveredFacePlus()
  {
    return m_coveredFacePlusG4;
  }
  static IntVect s_debugIV;
  static int     s_whichLev;

  ///
  /**
     set to true if the source you will provide
     is in conservative variables.   Default is false
  */
  static void useConservativeSource(bool a_conservativeSource)
  {
    s_conservativeSource = a_conservativeSource;
  }

  //these are the objects I need to make updatePrim less of a dog
  struct
  {
    int     m_offset;
    bool    m_multiValued;
  } typedef pointerOffset_t;

  struct
  {
    pointerOffset_t m_vofOffset;
  } typedef updateStencil_t;

  //gets offsets
  void fillUpdateStencil(EBPatchGodunov::updateStencil_t& a_sten,
                         const VolIndex&                  a_vof);

  //defines and fills cache
  void   cacheEBCF(Vector<Vector<Real> >&       a_cache,
                   const EBCellFAB&             a_input);

  //puts values of cache into output
  void uncacheEBCF(EBCellFAB&                   a_output,
                   const Vector<Vector<Real> >& a_cache);

  //for testing
  void getArgBox(Box a_argBox[SpaceDim]);
  //for testing
  void getEntireBox(Box a_entireBox[SpaceDim])
  {
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        a_entireBox[idir] = m_entireBox[idir];
      }
  }
protected:

  //the wild variety of crap I need to stencilize slopes
  Box m_modBoxOpen[SpaceDim];

  static bool s_conservativeSource;
  //stuff for debugging output
  static int  s_curLevel;
  static int  s_curComp;
  static int  s_doingVel;
  static int  s_doingAdvVel;
  static bool s_verbose;
  static Real s_maxWaveSpeed;
  static IntVect s_maxWaveSpeedIV;
  //set in define()
  ProblemDomain m_domain;
  RealVect m_dx;
  Real m_dxScale;
  bool m_useAgg;
  BaseIFFAB<FaceStencil> m_interpStencils[SpaceDim];
  Box m_validBox;

  EBISBox m_ebisBox;
  Real m_time;
  Real m_dt;
  bool m_isDefined;
  bool m_isBCSet;
  bool m_isBoxSet;
  bool m_isSlopeSet;
  bool m_isArtViscSet;
  bool m_useFourthOrderSlopes;
  bool m_useFlattening;
  bool m_useLimiting;
  bool m_useArtificialVisc;
  Vector<VolIndex> m_irregVoFs;

  Vector<updateStencil_t>  m_updateStencil;
  //set by factory
  EBPhysIBC* m_bc;

  //all stuff sent to UpdatePrim, etc needs to be built with this
  //set of info if any of these stenciling optimizations are to have any hope
  Box                     m_validBoxG4;
  IntVectSet       m_coveredSetsPlusG4[SpaceDim];
  IntVectSet       m_coveredSetsMinuG4[SpaceDim];
  Vector<VolIndex> m_coveredFacePlusG4[SpaceDim];
  Vector<VolIndex> m_coveredFaceMinuG4[SpaceDim];

  BaseIVFAB<Real>  m_extendStatePlusG4[SpaceDim];
  BaseIVFAB<Real>  m_extendStateMinuG4[SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxPlusG4[SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxMinuG4[SpaceDim];

  BaseIVFAB<Real>  m_extendStatePlus3D[SpaceDim][SpaceDim];
  BaseIVFAB<Real>  m_extendStateMinu3D[SpaceDim][SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxPlus3D[SpaceDim][SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxMinu3D[SpaceDim][SpaceDim];

  EBCellFAB        m_primPlus[SpaceDim];
  EBCellFAB        m_primMinu[SpaceDim];
  EBCellFAB        m_primState;
  EBCellFAB        m_primMinuTemp;
  EBCellFAB        m_primPlusTemp;

  EBFluxFAB        m_primGdnv;
  EBFaceFAB        m_fluxOne[SpaceDim];
  EBFaceFAB        m_fluxTwo[SpaceDim][SpaceDim];

  BaseIVFAB<Real>  m_coveredFluxNormMinu[SpaceDim];
  BaseIVFAB<Real>  m_coveredFluxNormPlus[SpaceDim];
  BaseIVFAB<Real>  m_extendStateNormMinu[SpaceDim];
  BaseIVFAB<Real>  m_extendStateNormPlus[SpaceDim];

  //stuff for aggregated slopes
  void setSlopeStuff();

  RefCountedPtr< AggStencil<EBCellFAB, EBCellFAB> > m_slopStenLo[SpaceDim];
  RefCountedPtr< AggStencil<EBCellFAB, EBCellFAB> > m_slopStenHi[SpaceDim];
  Box m_entireBox[SpaceDim];

 struct
  {
    size_t offset;
    int    dataID;
  } typedef access_t;

  struct
  {
    //slopes have all the same box so they will have the same offset and type
    access_t slop_access; //offsets for slopes
    bool hasLo, hasHi;
  } typedef slop_logic_t;

  Vector<slop_logic_t> m_slopVec[SpaceDim];

};
#include "NamespaceFooter.H"
#endif
