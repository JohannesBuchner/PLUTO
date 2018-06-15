#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPATCHADVECT_H_
#define _EBPATCHADVECT_H_

#include "EBPatchGodunov.H"
#include "EBPatchGodunovFactory.H"
#include "NamespaceHeader.H"

#define TOL 1.e-12

///
/**
   An EBPatchGodunov for passive advection with an
   input velocity of a particular function value.
*/
class EBPatchAdvect : public EBPatchGodunov
{
public:

  ///
  EBPatchAdvect();

  ///
  ~EBPatchAdvect();

  ///
  /**
     Return number of components for primitive variables.
  */
  virtual int numPrimitives() const
  {
    return 1;
  }

  ///
  /**
     Returns number of components for flux variables.
  */
  virtual int numFluxes() const
  {
    return 1;
  }

  ///
  /**
     Returns number of components for conserved variables.
  */
  virtual int numConserved() const
  {
    return 1;
  }

  ///
  /**
   */
  virtual int numSlopes() const
  {
    return 1;
  }

  ///
  /**
   */
  void
  normalPred(EBCellFAB&       a_primLo,
             EBCellFAB&       a_primHi,
             const EBCellFAB& a_primState,
             const EBCellFAB& a_slopePrim,
             const Real&      a_scale,
             const int&       a_dir,
             const Box&       a_box);

  ///
  /**
   */
  void
  transversePred(EBCellFAB&       a_primLo,
                 EBCellFAB&       a_primHi,
                 const EBCellFAB& a_primState,
                 const EBCellFAB& a_slopePrim,
                 const Real&      a_scale,
                 const int&       a_dir,
                 const Box&       a_box);

  ///
  /**
   */
  void
  interpolateFluxToCentroids(BaseIFFAB<Real>      a_centroidFlux[SpaceDim],
                             const EBFluxFAB&     a_fluxInterpolant,
                             const IntVectSet&    a_irregIVS);

  ///
  /**
   */
  void
  extrapToCoveredFaces(BaseIVFAB<Real>&        a_extendedPrim,
                       const EBFaceFAB&        a_primFace,
                       const EBCellFAB&        a_primState,
                       const Vector<VolIndex>& a_coveredFaces,
                       const int&              a_faceDir,
                       const Side::LoHiSide&   a_sd,
                       const Box&              a_box);

  ///
  /**
   */
  void
  pointExtrapToCovered2D(Vector<Real>&           a_extrapVal,
                         const EBFaceFAB&        a_primFace,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);

  ///
  /**
   */
  void
  pointExtrapToCovered3D(Vector<Real>&           a_extrapVal,
                         const EBFaceFAB&        a_primFace,
                         const EBCellFAB&        a_primState,
                         const int&              a_faceDir,
                         const VolIndex&         a_vof,
                         const RealVect&         a_normal,
                         const Side::LoHiSide&   a_sd,
                         const int&              a_numPrim);
  ///
  /**
     default does -= u_idir d cons/ dx_idir.
     assumes fluxes are really prims.
  */
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
             const Real&             a_scale);

  void
  upwindSlope(EBCellFAB&       a_slopeUpWi,
              const EBCellFAB& a_primState,
              const int&       a_dir,
              const Box&       a_box);

  void
  pointGetSlopesUpwind(Real&            a_dql,
                       Real&            a_dqr,
                       Real&            a_dqc,
                       bool&            a_hasFacesLeft,
                       bool&            a_hasFacesRigh,
                       const VolIndex&  a_vof,
                       const EBCellFAB& a_primState,
                       const int&       a_dir,
                       const int&       a_ivar,
                       const bool&      a_verbose);

  /// floors if m_isMaxMinSet
  virtual void
  floorPrimitives(EBCellFAB&  a_primState,
                  const Box&  a_box)
  {
    if (m_isMaxMinSet)
      {
        CH_assert(numPrimitives()==1);
        // IntVectSet ivsIrreg = m_ebisBox.getIrregIVS(a_box);
        IntVectSet ivsIrreg(a_box);
        for (VoFIterator vofit(ivsIrreg, m_ebisBox.getEBGraph());vofit.ok(); ++vofit)
          {
            const VolIndex& vof = vofit();
            a_primState(vof, 0)  = Min(a_primState(vof, 0), m_maxVal);
            a_primState(vof, 0)  = Max(a_primState(vof, 0), m_minVal);
          }
      }
  }

  /// floors if m_isMaxMinSet
  virtual void
  floorPrimitives(BaseIVFAB<Real>&   a_primState,
                  const IntVectSet&  a_set)
  {
    if (m_isMaxMinSet)
      {
        CH_assert(numPrimitives()==1);
        for (VoFIterator vofit(a_set, m_ebisBox.getEBGraph());vofit.ok(); ++vofit)
          {
            const VolIndex& vof = vofit();
            a_primState(vof, 0)  = Min(a_primState(vof, 0), m_maxVal);
            a_primState(vof, 0)  = Max(a_primState(vof, 0), m_minVal);
          }
      }
  }

  /// floors if m_isMaxMinSet
  virtual void
  floorConserved(EBCellFAB&  a_consState,
                 const Box&  a_box)
  {
    if (m_isMaxMinSet)
      {
        MayDay::Error("EBPatchAdvect::floorConserved::EBCellFAB:: unexpected case");
      }
  }

  /// floors if m_isMaxMinSet
  virtual void
  floorConserved(BaseIVFAB<Real>&   a_consState,
                 const IntVectSet&  a_set)
  {
    if (m_isMaxMinSet)
      {
        MayDay::Error("EBPatchAdvect::floorConserved::BaseIVFAB:: unexpected case");
      }
  }

  /// default noop
  virtual void  getCoveredValuesPrim(Vector<Real>& a_covValues)
  {;}

  ///default noop
  virtual void  getCoveredValuesCons(Vector<Real>& a_covValues)
  {;}

  ///default noop
  virtual void setCoveredConsVals(EBCellFAB& a_consState)
  {;}

  ///
  /**
     just a copy as default
  */
  virtual void
  consToPrim(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box,
             int              a_logflag,
             bool             a_verbose=false);

  ///
  /**
     just a copy as default
   */
  virtual void
  consToPrim(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs);

  ///
  /**
     just a copy
   */
  virtual void
  consToPrim(BaseIVFAB<Real>&  a_primState,
             const EBCellFAB&  a_conState,
             const IntVectSet& a_ivs);

  ///
  /**
     just a copy
   */
  virtual void
  primToCons(EBCellFAB&       a_primState,
             const EBCellFAB& a_conState,
             const Box&       a_box);

  ///
  /**
     just a copy
   */
  virtual void
  primToCons(BaseIVFAB<Real>&       a_primState,
             const BaseIVFAB<Real>& a_conState,
             const IntVectSet&      a_ivs);

  ///
  /**
     sends back primitive state (NOT FLUX) using upwinding
  */
  virtual void
  riemann(EBFaceFAB&       a_primState,
          const EBCellFAB& a_primLeft,
          const EBCellFAB& a_primRight,
          const int&       a_dir,
          const Box&       a_box);

  ///
  /**
     sends back primitive state (NOT FLUX) using upwinding
  */
  virtual void
  riemann(BaseIVFAB<Real>&        a_coveredFlux,
          const BaseIVFAB<Real>&  a_extendedState,
          const EBCellFAB&        a_primState,
          const Vector<VolIndex>& a_region,
          const int&              a_dir,
          const Side::LoHiSide&   a_sd,
          const Box& a_box);

  ///
  /**
     sets eb flux to zero
   */
  virtual void
  computeEBIrregFlux(BaseIVFAB<Real>&  a_ebIrregFlux,
                     const EBCellFAB&  a_primState,
                     const EBCellFAB   a_slopePrim[SpaceDim],
                     const IntVectSet& a_irregIVS,
                     const EBCellFAB&  a_source);

  ///
  /**
   */
  void
  extrapolatePrim(EBFluxFAB&                       a_flux,
                  EBCellFAB&                       a_primState,
                  EBCellFAB                        a_slopePrim[SpaceDim],
                  EBCellFAB                        a_slopeNLim[SpaceDim],
                  Vector<BaseIVFAB<Real> * >&      a_coveredFluxMinu,
                  Vector<BaseIVFAB<Real> * >&      a_coveredFluxPlus,
                  const Vector<IntVectSet >&       a_coveredSetsMinu,
                  const Vector<IntVectSet >&       a_coveredSetsPlus,
                  const Vector<Vector<VolIndex> >& a_coveredFaceMinu,
                  const Vector<Vector<VolIndex> >& a_coveredFacePlus,
                  const EBCellFAB&                 a_flattening,
                  const EBCellFAB&                 a_consState,
                  const EBCellFAB&                 a_source,
                  const Box&                       a_box,
                  const DataIndex&                 a_dit,
                  bool                             a_verbose);

  ///
  /**
   */
  void
  extrapolateBCG(EBFluxFAB&                       a_flux,
                 // EBCellFAB&                       a_primState,
                 EBCellFAB                        a_slopePrim[SpaceDim],
                 EBCellFAB                        a_slopeNLim[SpaceDim],
                 const EBCellFAB&                 a_flattening,
                 const EBCellFAB&                 a_consState,
                 const EBCellFAB&                 a_source,
                 const Box&                       a_box,
                 const DataIndex&                 a_dit,
                 bool                             a_verbose);

  void
  extrapolateBCG(EBCellFAB           a_primMinu[SpaceDim],
                 EBCellFAB           a_primPlus[SpaceDim],
                 // EBCellFAB&          a_primState,
                 EBCellFAB           a_slopesPrim[SpaceDim],
                 EBCellFAB           a_slopesSeco[SpaceDim],
                 const EBCellFAB&    a_flattening,
                 const EBCellFAB&    a_consState,
                 const EBCellFAB&    a_source,
                 const Box&          a_box,
                 const DataIndex&    a_dit,
                 bool                a_verbose);

  ///
  virtual void
  advectiveDerivative(EBCellFAB                      & a_udotDelRho,
                      const EBFluxFAB                & a_faceRho,
                      const EBFluxFAB                & a_faceVel,
                      const Vector<BaseIVFAB<Real>*> & a_coveredRhoMinu,
                      const Vector<BaseIVFAB<Real>*> & a_coveredRhoPlus,
                      const Vector<BaseIVFAB<Real>*> & a_coveredVelMinu,
                      const Vector<BaseIVFAB<Real>*> & a_coveredVelPlus,
                      const Vector<Vector<VolIndex> >& a_coveredFaceMinu,
                      const Vector<Vector<VolIndex> >& a_coveredFacePlus,
                      const Box                      & a_box);

  ///
  void
  advectiveDerivative(EBCellFAB                      & a_udotDelRho,
                      const EBFluxFAB                & a_faceRho,
                      const EBFluxFAB                & a_faceVel,
                      const Box                      & a_box);

  //Things you do not want to mess with unless you really know what you are doing
  ///
  /**
     always false for now
  */
  virtual bool usesFlattening() const;

  ///
  /**
     always false for now
  */
  virtual bool usesArtificialViscosity() const;


  // virtual void
  // coveredExtrapSlopes(Real&            a_dqc,
  //                     const VolIndex&  a_vof,
  //                     const EBCellFAB& a_primState,
  //                     const int&       a_dir,
  //                     const int&       a_ivar);

  // virtual void
  // coveredExtrapSlopesEBPA(Real&            a_dqc,
  //                         const VolIndex&  a_vof,
  //                         const EBCellFAB& a_primState,
  //                         const Real&      a_wc,
  //                         const int&       a_dir,
  //                         const int&       a_ivar);

  virtual void
  slope(EBCellFAB&       a_slopePrim,
        EBCellFAB&       a_slopeNLim,
        const EBCellFAB& a_primState,
        const EBCellFAB& a_flattening,
        const int&       a_dir,
        const Box&       a_box,
        bool a_doAgg = false);

  virtual void
  doSecondOrderSlopes(EBCellFAB&       a_delta2W,
                      EBCellFAB&       a_deltaWL,
                      EBCellFAB&       a_deltaWR,
                      EBCellFAB&       a_deltaWC,
                      const EBCellFAB& a_primState,
                      const int&       a_dir,
                      const Box&       a_box);

  ///
  /**
     Return true if you are using fourth-order slopes.
     Return false if you are using second-order slopes.
  */
  virtual bool usesFourthOrderSlopes() const;

  ///
  void
  averageVelToCC(EBCellFAB&                        a_normalVel,
                 const EBFluxFAB&                  a_advectionVel,
                 const Vector<BaseIVFAB<Real> * >& a_coveredVeloLo,
                 const Vector<BaseIVFAB<Real> * >& a_coveredVeloHi,
                 const Vector<Vector<VolIndex> >&  a_coveredFaceLo,
                 const Vector<Vector<VolIndex> >&  a_coveredFaceHi,
                 const Box&                        a_box) const;

  virtual void setVelocities(const EBCellFAB& a_normalVel,
                             const EBFluxFAB& a_advectionVel);

  virtual void
  setValidBox(const Box& a_validBox,
              const EBISBox& a_ebisBox,
              const IntVectSet& a_coarseFineIVS,
              const Real& a_time,
              const Real& a_dt);

  virtual void setTimeAndDt(const Real& a_time, const Real& a_dt)
  {
    m_time = a_time;
    m_dt   = a_dt;
  }

  ///sets generic defaults
  virtual Vector<string> stateNames();
  ///sets generic defaults
  virtual Vector<string> primNames();

  ///
  void useLimiting(bool a_limiting);

  // Things that should not be called because they do not apply

  /// should not be called
  Real
  getMaxWaveSpeed(const EBCellFAB& a_consState,
                  const Box& a_box)
  {
    MayDay::Error("should not be called ");
    return -1.0;
  }
  /// should not be called
  void
  getFlux(EBFaceFAB&       a_flux,
          const EBFaceFAB& a_prim,
          const int&       a_dir,
          const Box&       a_box)
  {
    MayDay::Error("should not be called ");
  }

  ///
  /**
     should not be called
  */
  Interval velocityInterval() const
  {
    MayDay::Error("should not be called");
    Interval retval;
    return retval;
  }

  ///
  /**
     only used as a comparison with ivar for flattening
     at covered face extrapolation
     don't want flattening-type stuff to pure advected vars
     also used for tagging cells.
  */
  virtual int densityIndex() const
  {
    return 0;
  }


  ///
  /**
     only used as a comparison with ivar for flattening
     at covered face extrapolation
     don't want flattening-type stuff to pure advected vars
  */
  int pressureIndex() const
  {
    return -1;
  }

  ///     should not be called
  int bulkModulusIndex() const
  {
    MayDay::Error("should not be called");
    return -1;
  }

  ///     should not be called
  Real artificialViscosityCoefficient() const
  {
    MayDay::Error("should not be called");
    return -1.0;
  }

  void setMaxMin(const Real& a_maxVal,
                 const Real& a_minVal)
  {
    m_isMaxMinSet = true;
    m_maxVal = a_maxVal;
    m_minVal = a_minVal;
  }
  EBCellFAB& getPrimState()
  {
    return m_primState;
  }



protected:

  const EBFluxFAB* m_advectionVelPtr;
  const EBCellFAB* m_normalVelPtr;

  bool m_isVelSet;
  bool m_isMaxMinSet;
  Real m_maxVal;
  Real m_minVal;
private:
  //disallowed for all the usual reasons
  void operator=(const EBPatchAdvect& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchAdvect(const EBPatchAdvect& a_input)
  {
    MayDay::Error("invalid operator");
  }

};


///
/**
 */
class EBPatchAdvectFactory: public EBPatchGodunovFactory
{
public:
  ///
  /**
   */
  virtual ~EBPatchAdvectFactory();

  ///
  /**
   */
  EBPatchGodunov* create() const;

  ///
  /**
   */
  EBPatchAdvectFactory(RefCountedPtr<EBPhysIBCFactory>&  a_bc,
                       const bool&                       a_useLimiting,
                       const Real&                       a_maxVal =-1.e99,
                       const Real&                       a_minVal = 1.e99,
                       const bool&                       a_setMaxMin = false);

protected:
  RefCountedPtr<EBPhysIBCFactory>  m_bcFactoryPtr;
  bool                             m_useLimiting;
  bool                             m_setMaxMin;
  Real                             m_maxVal;
  Real                             m_minVal;

private:

  //disallowed for all the usual reasons
  EBPatchAdvectFactory()
  {
    MayDay::Error("invalid operator");
  }
  void operator=(const EBPatchAdvectFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBPatchAdvectFactory(const EBPatchAdvectFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
