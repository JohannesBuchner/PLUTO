#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TIMEINTERPOLATORRK4_H_
#define _TIMEINTERPOLATORRK4_H_

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Time interpolator class using 4th-order Runge-Kutta

/**
 */
class TimeInterpolatorRK4
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  TimeInterpolatorRK4();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~TimeInterpolatorRK4();

  /// Actual constructor.
  /**
     Set up object.
   */
  void define(/// layout at this level
              const DisjointBoxLayout&  a_thisDisjointBoxLayout,
              /// layout at next coarser level
              const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
              /// problem domain on this level
              const ProblemDomain&      a_domain,
              /// refinement ratio between this level and next coarser level
              const int&                a_refineCoarse,
              /// number of variables
              const int&                a_numStates,
              /// layers of ghost cells to be filled in on the coarsened layout at this level
              const int&                a_ghosts);

  /// Set coarse timestep.
  void setDt(const Real&  a_dt);

  /// Update Taylor polynomial coefficients with the coarse solution.
  /**
     Update Taylor polynomial coefficients with coarse solution a_soln.
     Ghost cells, too.

     This function must be called after setDt() and before saveRHS().
   */
  void saveInitialSoln(const LevelData<FArrayBox>&   a_soln);

  /// Update Taylor polynomial coefficients with coarse right-hand side.
  /**
     Update Taylor polynomial coefficients with coarse right-hand side.
     Ghost cells, too.

     This function must be called exactly four times, after saveInitialSoln()
     and before any calls to interpolate().

     The counter m_countRHS keeps track of how many times this is called.
   */
  void saveRHS(const LevelData<FArrayBox>&   a_rhs);

  /// Interpolate in time using the Taylor polynomial.
  /**
     Interpolate in time to a_U on interval a_intvl using the Taylor polynomial.
     Ghost cells, too.
   */
  void interpolate(/// interpolated solution on this level coarsened
                   LevelData<FArrayBox>&   a_U,
                   /// time interpolation coefficient in range [0:1]
                   const Real&             a_timeInterpCoeff,
                   /// interval of a_U to fill in
                   const Interval&         a_intvl);

  /// Set RK4 intermediate values in time using the Taylor polynomial.
  /**
     Set RK4 intermediate values to a_U on interval a_intvl
     using the Taylor polynomial.
     Ghost cells, too.
   */
  void intermediate(/// interpolated solution on this level coarsened
                   LevelData<FArrayBox>&   a_U,
                   /// time interpolation coefficient in range [0:1]
                   const Real&             a_timeInterpCoeff,
                   /// which RK4 stage:  0, 1, 2, 3
                   const int&              a_stage,
                   /// interval of a_U to fill in
                   const Interval&         a_intvl) const;

protected:

  /// Reset this object to use with new data
  void resetData();

  /// Set a_vec = m_dt * (a_c0, a_c1, a_c2, a_c3).
  void setVectorDt(Vector<Real>& a_vec, Real a_c0, Real a_c1, Real a_c2, Real a_c3);

  Vector< Vector<Real> > m_coeffs;

  /// layout for this level, coarsened by m_refToCoarse
  DisjointBoxLayout m_thisCoarsenedLayout;

  /// layout for the coarse level
  DisjointBoxLayout m_coarseLayout;

  /// For copying from rhs on m_coarseLayout to m_rhsCopy on m_thisCoarsenedLayout
  Copier m_copier;

  /// layers of ghost cells around m_thisCoarsenedLayout for m_rhsCopy and m_taylorCoeffs
  int m_ghosts;

  /// number of coefficients in Taylor polynomial:  this is 4
  int m_numCoeffs;

  /// ghost vector around m_thisCoarsenedLayout for m_rhsCopy and m_taylorCoeffs
  IntVect m_ghostVect;

  /// coarse timestep
  Real m_dt;

  /// Copy of rhs on m_thisCoarsenedLayout, to be used within saveRHS
  LevelData<FArrayBox> m_rhsCopy;

  /// coefficients of the third-degree Taylor polynomial on m_thisCoarsenedLayout with ghost vector m_ghostVect; m_numCoeffs*m_numStates components
  LevelData<FArrayBox> m_taylorCoeffs;

  /// difference between f(U1) and f(U2), used for intermediate-value calculations, on m_thisCoarsenedLayout with ghost vector m_ghostVect; m_numCoeffs*m_numStates components
  LevelData<FArrayBox> m_diff12;

  /// whether m_dt has been set
  bool m_gotDt;

  /// whether initial solution has been saved
  bool m_gotInitialSoln;

  /// whether we have the full Taylor polynomial
  bool m_gotFullTaylorPoly;

  /// number of times saveRHS function has been called
  int m_countRHS;

  /// Problem domain - index space for next coarser level
  ProblemDomain m_coarseDomain;

  /// Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  /// Number of variables
  int m_numStates;

  /// define() has been called
  bool m_defined;

private:

  // Disallowed for all the usual reasons
  void operator=(const TimeInterpolatorRK4&);
  TimeInterpolatorRK4(const TimeInterpolatorRK4&);
};

#include "NamespaceFooter.H"
#endif
