#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRCNS_H_
#define _EBAMRCNS_H_

#include "AMR.H"
#include "AMRLevel.H"
#include "EBCellFAB.H"
#include "EBLevelTGA.H"
#include "BaseIVFAB.H"
#include "LevelData.H"
#include "EBLevelRedist.H"
#include "EBCoarToCoarRedist.H"
#include "EBPatchGodunov.H"
#include "EBPatchGodunovFactory.H"
#include "EBBackwardEuler.H"
#include "EBCoarToFineRedist.H"
#include "EBFineToCoarRedist.H"
#include "EBCoarseAverage.H"
#include "EBPWLFineInterp.H"
#include "EBFluxRegister.H"
#include "EBLevelCNS.H"
#include "Box.H"
#include "BiCGStabSolver.H"
#include "EBSimpleSolver.H"
#include "EBAMRCNSParams.H"
#include "IntVectSet.H"
#include "AMRMultiGrid.H"
#include "EBPatchPolytropicFactory.H"
#include "Vector.H"
#include "MomentumTGA.H"
#include "DisjointBoxLayout.H"
#include "EBQuadCFInterp.H"
#include "EBSpaceTimeFunction.H"

#include "NamespaceHeader.H"

//! \class EBAMRCNS
//! This class implements an embedded boundary compressible Navier-Stokes solver with 
//! adaptive mesh refinement. Viscosity and thermal conductivity are represented by 
//! operators that can vary in space and time. The specific heat of the medium, on the 
//! other hand, is assumed to be constant. The solver stores its data in conserved form,
//! and the data is laid out within FABs in the following order: mass density, momentum 
//! components, energy density.
class EBAMRCNS: public AMRLevel
{
public:

  static bool s_noEBCF;

  //! Creates a new embedded boundary AMR-enabled compressive Navier-Stokes solver.
  //! \param a_params A set of parameters to be used by the solver.
  //! \param a_godFactory A factory for creating patches.
  //! \param a_initialConditions A function of space and time that defines 
  //!                            initial conditions for the solver.
  EBAMRCNS(const EBAMRCNSParams& a_params,
           const RefCountedPtr<EBPatchPolytropicFactory>& a_godFactory,
           const RefCountedPtr<EBSpaceTimeFunction>& a_initialConditions);

  //! Destructor.
  virtual ~EBAMRCNS();

  //! Recompute the specific momentum and energy (conserved quantities) 
  //! in the given state with updated values for the velocity and temperature, 
  //! assuming no change in density.
  void updateMomentumAndEnergy(const Vector<LevelData<EBCellFAB>*>& a_newVelocity,
                               const Vector<LevelData<EBCellFAB>*>& a_newTemperature,
                               int a_baseLevel,
                               int a_maxLevel);

  //! Returns a reference to the FAB containing the updated (conserved) state.
  LevelData<EBCellFAB>& getNewState();


  bool convergedToSteadyState();

  //! Returns a reference to the FAB containing the conserved state as it is at the 
  //! beginning of a time step.
  LevelData<EBCellFAB>& getOldState();

  void postInitialGrid(const bool a_restart); 

  ///
  int nGhost() const
  {
    return m_nGhost;
  }
  ///
  int nComp() const
  {
    return m_nComp;
  }

  //! Returns a Vector containing new state data for each of the refinement 
  //! levels in the AMR hierarchy.
  Vector<LevelData<EBCellFAB>*> getNewStates();
  
  //! Returns a Vector containing old state data for each of the refinement 
  //! levels in the AMR hierarchy.
  Vector<LevelData<EBCellFAB>*> getOldStates();

  //! Returns the computed time step size.
  Real getDt() const;

  //! get the LAST dt
  Real getDtOld()
  {
    return m_dtOld;
  }
  //! Returns the layout of the index spaces for the current level.
  EBISLayout getEBISLayout() const;

  //! Returns the layout of the index spaces for the current level.
  EBLevelGrid getEBLG() const
  {
    return m_eblg;
  }

  ///
  void    defineSolvers();

  //! Fills the given FAB with conserved and primitive data. The stride of the data is 
  //! equal to the stride of the conserved variables plus that of the primitive data, 
  //! and the primitive data will appear after the conserved data in the FAB.
  //! \param a_data The FAB to which the conserved and primitive data are to be written.
  //!               \a a_data need not be defined before this method is called.
  void fillConsAndPrim(LevelData<EBCellFAB>& a_data, LevelData<EBCellFAB>& a_state, int logflag) ;

  //! Extracts the velocity from the conserved state.
  //! \param a_velo A FAB that will hold the extracted velocity.
  //! \param a_state The conserved state from which the velocity is extracted.
  void getVelocity(LevelData<EBCellFAB>&       a_velo,
                   const LevelData<EBCellFAB>& a_state) const;

  //! Extracts the bulk temperature from the conserved state.
  //! \param a_temp A FAB that will hold the extracted temperature.
  //! \param a_state The conserved state from which the temperature is extracted.
  void getTemperature(LevelData<EBCellFAB>&       a_temp,
                      const LevelData<EBCellFAB>& a_state) const;

  //! Sums the conserved variable identified by \a_ivar over the entire grid level, placing 
  //! the result in \a a_sumcons.
  //! \param a_sumcons The value of the conserved variable summed over the grid level.
  //! \param a_ivar The index of the conserved variable within its FAB.
  void sumConserved(Real& a_sumcons,
                    const int& a_ivar) const;

  //! Call this to tell the solver that there are no external sources of 
  //! mass, momentum, or energy, and that these quantities should therefore 
  //! be conserved.
  void setConservative(bool a_conservative);

  virtual void syncWithFineLevel();

  void setAirDiffusionCoefficients(const LevelData<EBCellFAB>& a_densCell,
                                   const LevelData<EBCellFAB>& a_tempCell,
                                   const LevelData<EBFluxFAB>& a_densFace,
                                   const LevelData<EBFluxFAB>& a_tempFace);

  void
  setAirThermDiff(const LevelData<EBCellFAB>& a_massDens,
                  const LevelData<EBCellFAB>& a_temperat);

  // ------------------------------------------------
  // AMRLevel methods overridden 
  // ------------------------------------------------
  
  // This appears to override AMRLevel::define, but the second parameter in this case is 
  // a ProblemDomain, which is implicitly casted to a Box, so I'm not sure how the 
  // compiler actually treats it. -JNJ
  virtual void define(AMRLevel*            a_coarser_level_ptr,
                      const ProblemDomain& a_problem_domain,
                      int                  a_level,
                      int                  a_ref_ratio);

  virtual Real advance();

  virtual void postTimeStep();

  virtual void tagCells(IntVectSet& a_tags);

  virtual void tagCellsInit(IntVectSet& a_tags);

  virtual void regrid(const Vector<Box>& a_new_grids);

  virtual void preRegrid(int                         a_base_level,
                         const Vector<Vector<Box> >& a_new_grids);

  virtual void postRegrid(int a_base_level);

  virtual void initialGrid(const Vector<Box>& a_new_grids);

  virtual void initialData();

  virtual void postInitialize();

  virtual Real computeDt();

  virtual Real computeInitialDt();

  ///made public to allow convergence tests
  void hyperbolicSource(LevelData<EBCellFAB>& a_source);

  ///made public to allow convergence tests
  void fluxDivergence(  LevelData<EBCellFAB>& a_divergeF);

  /// set  output to volfrac*((del dot (kappa grad T))  + div(sigma u))
  void kappaEnergySource(LevelData<EBCellFAB>& a_kappaEnergySource,
                         const LevelData<EBCellFAB>& a_velocity,
                         const LevelData<EBCellFAB>* a_veloCoar,
                         const LevelData<EBCellFAB>& a_temperat,
                         const LevelData<EBCellFAB>* a_tempCoar,
                         const LevelData<EBCellFAB>& a_state);//for coefficients

  ///set  output to volfrac*div(sigma)
  void kappaMomentumSource(LevelData<EBCellFAB>& a_kappaDivSigma,
                           const LevelData<EBCellFAB>& a_velocity,
                           const LevelData<EBCellFAB>* a_veloCoar,
                           const LevelData<EBCellFAB>& a_state);//for coefficients

  //if doNormalization = true, averages with neighboring cells
  //otherwise returns kappa * source
  void explicitHyperbolicSource(LevelData<EBCellFAB>&       a_momentSource, 
                                LevelData<EBCellFAB>&       a_energySource,
                                const LevelData<EBCellFAB>& a_state,
                                bool a_doNormalization );

  ///inviscid version of the algorithm
  void explicitAdvance(const LevelData<EBCellFAB>& a_divergeF);

  // advances the solution using the hyperbolic terms.
  //(\kappa \rho^\npo I - \frac{\kappa \dt}{2} L^m)\ubar^\npo = \kappa \rho^\npo \ubold^*
  void  getUStar(LevelData<EBCellFAB>      & a_UStar, 
                 const LevelData<EBCellFAB>& a_UN,
                 const LevelData<EBCellFAB>& a_divergef);

  //this one is explicit and also does reweighting of the 
  //redistribution object in the case of mass weighting
  void hyperbolicRedistribution(LevelData<EBCellFAB>& a_state);

  /// gets divergence of shear stress
  /**
     (\rho^* I - \dt L^m) \ubold^\npo = \dt (\rho \ubold)^*
     L^m(\ubold^\nph) = \rho^*\left(\frac{\ubold^\npo - \ubold^*}{\dt}\right)
     
     Adds divSigma into the momentum of U*
  */
  void  getDivSigma(LevelData<EBCellFAB>& a_divSigma, 
                    LevelData<EBCellFAB>& a_UStar);

  ///
  /**
     $L^d(U^*) =u \cdot \grad \cdot \sigma) + \sigma \cdot (\grad \ubold)$
     (\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
  */
  void  
  getSplitLdOfU(LevelData<EBCellFAB>      & dissFunc, 
                LevelData<EBCellFAB>      & uDotDivSig,  
                LevelData<EBCellFAB>      & UStar, 
                const LevelData<EBCellFAB>& a_divSigma);


  void  
  getSingleLdOfU(LevelData<EBCellFAB>      & a_divSigmaU, 
                 LevelData<EBCellFAB>      & UStar);

  ///
  /**
     (\rho E)^{**} = (\rho E)^* + \dt L^d(U^*)
     mass diff = kappa(1-kappa)*(kappaConsDissFcn - nonConsDissFcn)
  */
  void
  updateEnergyBySplitLdAndRedistribute(LevelData<EBCellFAB>      & a_dissFunc, 
                                       LevelData<EBCellFAB>      & a_uDotDivSig,  
                                       LevelData<EBCellFAB>      & a_kappaConsDissFunc,
                                       LevelData<EBCellFAB>      & a_nonConsDissFunc,
                                       LevelData<EBCellFAB>      & a_UStar);

  void
  updateEnergyBySingleLdAndRedistribute(LevelData<EBCellFAB>      & a_divSigmaU,
                                        LevelData<EBCellFAB>      & a_kappaConsDivSigmaU,
                                        LevelData<EBCellFAB>      & a_nonConsDivSigmaU,
                                        LevelData<EBCellFAB>      & a_UStar);

  void
  floorVariable(LevelData<EBCellFAB>& a_data,
                const EBLevelGrid a_eblg,
                Real a_minVal, int a_ivar);
  //
  /** (\kappa \rho^\npo C_v I - \frac{\kappa \dt} L^k)T^\npo = \kappa \rho^\npo C_v T^{**}
      (\rho E)^{n+1} = (\rho E)^{**} + \dt L^k(T**)
  **/
  void getDivKappaGradT(LevelData<EBCellFAB>& a_divKappaGradT,
                        LevelData<EBCellFAB>& a_Ustar);

  //this includes news so need to call delete
  void
  getCoarserTemperatureAndVelocity(LevelData<EBCellFAB>* & a_tempCoar,
                                   LevelData<EBCellFAB>* & a_veloCoar);


  //put ustar into state
  void finalAdvance(LevelData<EBCellFAB>& a_Ustar);
#ifdef CH_USE_HDF5
  // ---------------------------------------------------
  // HDF5 routines for plotting and writing checkpoints.
  // ---------------------------------------------------
 
  static int s_NewPlotFile;
  static bool s_solversDefined;
  virtual void writePlotHeaderOld    (HDF5Handle& a_handle) const;
  virtual void writePlotLevelOld     (HDF5Handle& a_handle) const;
  virtual void writePlotHeader       (HDF5Handle& a_handle) const;
  virtual void writePlotLevel        (HDF5Handle& a_handle) const;
  virtual void writeCheckpointHeader (HDF5Handle& a_handle) const;
  virtual void writeCheckpointLevel  (HDF5Handle& a_handle) const;
  virtual void readCheckpointHeader  (HDF5Handle& a_handle);
  virtual void readCheckpointLevel   (HDF5Handle& a_handle);
#endif

protected:

  void defineFactories(bool a_atHalfTime);

  void implicitReflux();

  void refluxRedistInteraction();

  void resetWeights();

  // This performs the redistribution of a quantity (mass, momentum, or energy) in m_stateNew. 
  // The quantity is determined by a_interv. If a_putIntoState is set, the redistribution is 
  // done EXPLICITLY, and the redistributed quantity is placed directly into m_stateNew.
  // Otherwise, the quantity is accumulated into m_redisRHS for the later implicit 
  // redistribution.
  void refluxRHSEnergyAndMomentum();
  void explicitReflux(const Interval& a_interv);
  void coarseFineRedistribution(const Interval& a_interv);

  int getFinestLevel();

  void getHalfState(LevelData<EBCellFAB>& a_stateInt);
  EBAMRCNSParams m_params;

  EBLevelGrid    m_eblg;

  LevelData<EBCellFAB>         m_stateOld;
  LevelData<EBCellFAB>         m_stateNew;
  LevelData<EBCellFAB>         m_rhoUDotDele;
  LevelData<EBCellFAB>         m_pDivU;
  LevelData<EBCellFAB>         m_gradU;

  LevelData<EBCellFAB>         m_redisRHS;
  LevelData<BaseIVFAB<Real> >  m_massDiff;
  LayoutData<IntVectSet>       m_sets;

  EBCoarToFineRedist m_ebCoarToFineRedist;
  EBCoarToCoarRedist m_ebCoarToCoarRedist;
  EBFineToCoarRedist m_ebFineToCoarRedist;

  Real m_dtOld;
  int m_nComp;
  int m_nPrim;
  int m_nGhost;
  Vector<string> m_stateNames;
  Vector<string> m_primNames;
  RealVect m_dx;
  Real m_dtNew;
  bool m_hasCoarser, m_hasFiner;
  EBCoarseAverage                                               m_ebCoarseAverage;
  EBPWLFineInterp                                               m_ebFineInterp;
  EBFluxRegister                                                m_divFFluxRegister;
  EBFluxRegister                                                m_veloFluxRegister;
  EBFluxRegister                                                m_tempFluxRegister;
  RefCountedPtr<EBPatchGodunov>                                 m_ebPatchGodunov;
  RefCountedPtr<EBPatchPolytropicFactory>                       m_ebPatchGodunovFactory;
  EBLevelCNS                                                    m_ebLevelGodunov;
  EBLevelRedist                                                 m_ebLevelRedist;

  //! Initial conditions function.
  RefCountedPtr<EBSpaceTimeFunction>                            m_ICs;

  //! Initial values of conserved quantities.
  Real m_originalMass;
  RealVect m_originalMomentum;
  Real m_originalEnergy;

  //! This flag is set if the system is conservative.
  bool m_isConservative;

  //! Beginning-of-step A coefficient for the thermal diffusion equation.
  RefCountedPtr<LevelData<EBCellFAB> >                          m_acoTemp;

  //! End-of-step A coefficient for the viscous diffusion equation.
  RefCountedPtr<LevelData<EBCellFAB> >                          m_acoVelo;

  RefCountedPtr<LevelData<EBFluxFAB> >                          m_eta;
  RefCountedPtr<LevelData<EBFluxFAB> >                          m_bcoTemp;
  RefCountedPtr<LevelData<EBFluxFAB> >                          m_lambda;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >                   m_etaIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >                   m_bcoTempIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >                   m_lambdaIrreg;
  RefCountedPtr<EBQuadCFInterp>                                 m_quadCFI;

  static RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >  s_tempFactory;
  static RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >  s_veloFactory;
  static RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >  s_tempSolver;
  static RefCountedPtr<AMRMultiGrid<     LevelData<EBCellFAB> > >  s_veloSolver;

  //  static BiCGStabSolver<LevelData<EBCellFAB> > s_botSolver;
  static EBSimpleSolver s_botSolver;
  static RefCountedPtr<MomentumBackwardEuler > s_veloIntegratorBE;
  static RefCountedPtr< EBLevelBackwardEuler > s_tempIntegratorBE;

  static RefCountedPtr<MomentumTGA > s_veloIntegratorTGA;
  static RefCountedPtr< EBLevelTGA > s_tempIntegratorTGA;
  void fillCoefficients(const LevelData<EBCellFAB>& a_state);
  void     clearSolvers();

  //(rho I - dt Lv) delta = dt*Dr(Fm)
  void getRefluxDeltaV(Vector<LevelData<EBCellFAB>* >&  a_deltaVelocity, 
                       Vector<LevelData<EBCellFAB>* >&  a_dtRrefluxDivergeM, 
                       int baseLev, int finestLev, Real baseDt);


  //(rho Cv I - dt Lt) delta = dt*Dr(Fe)
  void getRefluxDeltaT(Vector<LevelData<EBCellFAB>* >& a_deltaTemperat,
                       Vector<LevelData<EBCellFAB>* >& a_dtRefluxDivergeE, 
                       int baseLev, int finestLev, Real baseDt);

  //mom += dt*rho*(deltav)
  void incrMomentumByDeltaV(Vector<LevelData<EBCellFAB>* >& a_deltaVelocity, 
                            int baseLev, int finestLev);

  //mom += dt*rho*(deltav)
  void incrMomentumByDeltaV(LevelData<EBCellFAB>& a_state, 
                            LevelData<EBCellFAB>& a_deltaVelocity);

  //ene += dt*rho*Cv*(deltaT)
  void incrEnergyByDeltaT(Vector<LevelData<EBCellFAB>* >& a_deltaTemperat, 
                          int baseLev, int finestLev);
  void incrEnergyByDeltaT(LevelData<EBCellFAB>& a_state,
                          LevelData<EBCellFAB>& a_deltaTemperat);



  void coarseFineIncrement();
  void postTimeStepRefluxRedistDance();

private:

  Real calculateMass() const;

  EBAMRCNS* getCoarserLevel() const;
  EBAMRCNS* getFinerLevel() const;


  void levelSetup();
  //weak construction is an
  //evil that must be rooted out.
  //verily, it is a character flaw.
  EBAMRCNS()
  {
    MayDay::Error("invalid operator");
  }
  //disallowed for all the usual reasons
  void operator=(const EBAMRCNS& a_input)
  {
    MayDay::Error("invalid operator");
  }
  EBAMRCNS(const EBAMRCNS& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
