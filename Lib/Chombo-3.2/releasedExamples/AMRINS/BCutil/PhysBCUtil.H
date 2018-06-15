#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PHYSBCUTIL_H_
#define _PHYSBCUTIL_H_

#include "Tuple.H"
#include "ParmParse.H"
#include "LoHiSide.H"
#include "BCFunc.H"
// include support for *IBC stuff from AMRGodunov, for use with
// PatchGodunov-derived stuff
#include "PhysIBC.H"

#include "UsingNamespace.H"

/// Encapsulates all physical BC's for IAMR solves
/** This class is designed to contain all the information
needed to specify all physical boundary conditions needed
for the IAMR code. This includes BC's for velocity, pressures,
scalars.  While basic solid-wall BCs are coded in to this
class, the idea is that for a more complicated set of physical
boundary conditions (inflow, etc), a problem-dependent class
will be derived from this one which will override the generic
BC's in this class */
class PhysBCUtil
{

public:
  /** this will access the ParmParse database using the
      physBC prefix and use it to define the physical BC types (hi
      and lo) in each direction.  The BC type "custom" will imply
      that the BC is a somewhat nonstandard one, but that the
      derived class will know what to do with it. */
  PhysBCUtil();

  ///
  virtual ~PhysBCUtil();

  ///
  enum BCType
  {
    ///
    bogusBC = -1,
    ///
    SolidWall,
    ///
    Inflow,
    ///
    Outflow,
    ///
    Symmetry,
    ///
    noShear,
    ///
    Custom,
    ///
    NUM_PHYS_BC_TYPES
  };

  /// "virtual" constructor
  virtual PhysBCUtil* newPhysBCUtil() const;

  /// Velocity BC's

  /// returns velocity BC for viscous source
  virtual Tuple<BCHolder, SpaceDim> viscousVelFuncBC() const;

  /// returns velocity BC for tracing (inviscid)
  virtual Tuple<BCHolder, SpaceDim> tracingVelFuncBC() const;

  /// added by petermc, 20 nov 2007
  virtual BCHolder tracingSolveFuncBC(int a_idir) const;

  /// returns edge-centered vel BC for edge-ctrd vel (single component)
  virtual Tuple<BCHolder, SpaceDim> advectionVelFuncBC(bool a_isViscous) const;

  /// returns pre-projection velocity BC
  virtual Tuple<BCHolder, SpaceDim> uStarFuncBC(bool a_isViscous = true) const;

  /// returns BC on cell-centered u-dot-del-u term
  virtual Tuple<BCHolder, SpaceDim> uDelUFuncBC(bool a_isViscous) const;

  /// returns single-component BC for viscous refluxing solves
  virtual BCHolder viscousRefluxBC(int a_dir) const;

  /// returns a BC object compatible with AMRGodunov advection infrastructure
  virtual PhysIBC* advectionVelIBC() const;

  /// Scalar BC's

  /// added by petermc, 5 nov 2007
  virtual BCHolder scalarTraceFuncBC(int a_scalarType) const;

  /// extrapolation BC, order 1 or 2; added by petermc, 7 nov 2007
  virtual BCHolder extrapFuncBC(int a_order) const;

  ///
  virtual BCHolder scalarRefluxSolveBC(int a_scalarType) const;

  /// returns a BC object compatible with AMRGodunov advection infrastructure
  virtual PhysIBC* scalarTraceIBC(int a_scalarType) const;

  virtual PhysIBC* scalarTraceExtrapIBC(int a_scalarType) const;

  /// volume discrepancy BC's
  virtual BCHolder lambdaFuncBC() const;

  /// added by petermc, 12 nov 2007
  virtual BCHolder basicLambdaFuncBC() const;

  ///
  virtual PhysIBC* lambdaTraceIBC() const;

  /// Miscellaneous BC's

  // BC on stream function (used when computing stream func. from vorticity)
  virtual BCFunc streamBC() const;

  // added by petermc, 30 oct 2007
  virtual BCHolder FreestreamCorrFuncBC() const;

  /// sets basic physical BC on pressure (called by public BC functions)
  virtual BCHolder BasicPressureFuncBC(bool a_isHomogeneous) const;

  // added by petermc, 30 oct 2007
  virtual BCHolder LevelPressureFuncBC() const;

  // added by petermc, 30 oct 2007
  virtual BCHolder SyncProjFuncBC() const;

  /// sets basic physical BC on pressure before taking gradients
  virtual BCHolder BasicGradPressureFuncBC() const;

  // added by petermc, 30 oct 2007
  virtual BCHolder gradMacPressureFuncBC() const;

  // added by petermc, 30 oct 2007
  virtual BCHolder gradPiFuncBC() const;

  // added by petermc, 30 oct 2007
  virtual BCHolder gradESyncFuncBC() const;

  // added by petermc, 30 oct 2007
  virtual BCHolder gradELambdaFuncBC() const;

  /// returns single-component BC for viscous solves
  virtual BCHolder viscousSolveFuncBC(int a_dir) const;

  /// returns multi-component BC for viscous solves
  virtual BCHolder viscousFuncBC() const;

  /// miscellaneous utility functions

  /// sets cell spacing
  void Dx(const Real a_dx);

  /// returns cell spacing
  Real Dx() const;

  /// sets time
  void Time(const Real a_time);

  /// returns current time
  Real Time() const;

  ///
  int loBC(int a_dir) const;

  ///
  int hiBC(int a_dir) const;

  ///
  int getBC(int a_dir, const Side::LoHiSide& a_side) const;

  ///  this is in case the BC's have an effect on the timestep
  /** pass in currently computed dt, along with cfl and dx.  If
      effect of boundary conditions requires a decreased timestep,
      the newly reduced timestep is returned. In the default case,
      this just returns a_dt back; however, derived classes may
      actually have an effect.*/
  virtual void computeBoundaryDt(Real& a_dt, Real a_cfl, Real a_dx) const;

protected:

  /// returns a BCHolder for cell-centered velocity
  virtual BCHolder basicCCVelFuncBC(bool a_isHomogeneous,
                                    bool a_isViscous,
                                    int a_velComponent,
                                    const Interval& a_interval) const;

  /// returns a BCHolder for edge-centered velocity
  virtual BCHolder basicECVelFuncBC(bool a_isHomogeneous,
                                    bool a_isViscous,
                                    int a_comp,
                                    const Interval& a_interval) const;

  ///
  virtual BCHolder basicScalarFuncBC(bool a_isHomogeneous,
                                     int  a_scalarType) const;

  ///  interact with ParmParse to set physical BC types
  virtual void setBCs();

  // contains the enumerated values for the lo-end BC
  Tuple<int, SpaceDim> m_loBC;

  // contains the enumerated values for the hi-end BC
  Tuple<int, SpaceDim> m_hiBC;

private:
  /// cell spacing
  Real m_dx;

  /// current time
  Real m_time;

  /// copy constructor
  PhysBCUtil(const PhysBCUtil& rhs);

  ///
  virtual PhysBCUtil& operator=(const PhysBCUtil& rhs);
};

#endif
