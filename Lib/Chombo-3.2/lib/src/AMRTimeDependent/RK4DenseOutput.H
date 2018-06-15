#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RK4DENSEOUTPUT_H_
#define _RK4DENSEOUTPUT_H_


#include "NamespaceHeader.H"

/// This is a standard RK4 implementation
/// that also includes dense output coefficients

// Soln is the type that operators are applied to
// Rhs is the type that operator produce
// EXOp wraps the explicit operator
template <class Soln, class Rhs, class EXOP>
class RK4DenseOutput
{
public:

  RK4DenseOutput<Soln, Rhs, EXOP>() { m_isDefined = false; }

  // This must be called first.
  void define(const Soln& a_state, Real a_dt, bool a_denseOutput = false);

  // Advance one step.
  void advance(Real a_time, Soln& a_state);

  // Return current dense output coefs, 0th power first, etc.
  // NOTE: These are intended to be for <JU> conserved quantities
  void denseOutputCoefs(Vector<Rhs*>& a_interpCoefs);

  // Reset the timestep.
  void resetDt(Real a_dt);

  // Set whether the step starts at 0. and/or ends at 1.
  void start0end1(bool a_start0, bool a_end1);

  // Access to the operators if we're defined already
  // Caller is responsible for making sure they're in a good state 
  // (for example, like have resetDt() called, etc.
  EXOP& getEXOP();

  bool isDefined() const { return m_isDefined; }

  bool hasDenseOutput() const { return m_hasDenseOutput; }

  /// Runge-Kutta coefficients
  static const int  s_nStages = 4;
  static const Real s_c[s_nStages];
  static const Real s_a[s_nStages][s_nStages];
  static const Real s_b[s_nStages];
  static const int  s_nDenseCoefs = 3;
  static const Real s_bstar[s_nDenseCoefs][s_nStages];

protected:
  bool m_isDefined;
  bool m_denseOutput;
  bool m_hasDenseOutput; // If there is dense output data to interpolate
  Real m_dt;
  Real m_time;
  Soln m_phi[s_nStages];
  Rhs m_rhs;
  Rhs m_denseCoefs[s_nDenseCoefs];
  Rhs m_kE;
  EXOP m_opEx;

private:

};

//==============================================

template <class Soln, class Rhs, class EXOP>
void RK4DenseOutput<Soln, Rhs, EXOP>::
define(const Soln& a_state, Real a_dt, bool a_denseOutput)
{
  CH_TIMERS("RK4DenseOutput::define");
  m_dt = a_dt;
  m_denseOutput = a_denseOutput;
  // define Soln types
  for (int stage = 0; stage < s_nStages; stage++)
    m_phi[stage].define(a_state);

  // define Rhs types
  m_kE.define(a_state);
  m_rhs.define(a_state);

  // define opEx
  m_opEx.define(a_state, m_dt);

  // if dense output is requested, need more storage
  for (int coef = 0; m_denseOutput && (coef < s_nDenseCoefs); coef++)
    m_denseCoefs[coef].define(a_state);

  m_hasDenseOutput = false;
  m_isDefined = true;
}

/*
  Get a reference to the implicit-explicit operator
 */
template <class Soln, class Rhs, class EXOP>
EXOP& RK4DenseOutput<Soln, Rhs, EXOP>::
getEXOP()
{
  return m_opEx;
}

/*
  Reset the timestep.
 */
template <class Soln, class Rhs, class EXOP>
void RK4DenseOutput<Soln, Rhs, EXOP>::
resetDt(Real a_dt)
{
  CH_assert(isDefined());

  // Only update everything if dt has changed
  Real reltol = 1e-14;
  if (Abs(m_dt - a_dt) > m_dt*reltol)
  {
    m_dt = a_dt;
    m_hasDenseOutput = false;
    m_opEx.resetDt(m_dt);
  }
}

/*
  Set whether the step starts at 0. and/or ends at 1.
 */
template <class Soln, class Rhs, class EXOP>
void RK4DenseOutput<Soln, Rhs, EXOP>::
start0end1(bool a_start0, bool a_end1)
{
  CH_assert(isDefined());

  int stage0 = -1; // if no stage is at time 0.
  if (a_start0) stage0 = 0;

  int stage1 = -1; // if no stage is at time 1.
  if (a_end1) stage1 = s_nStages-1;

  m_opEx.stage0stage1(stage0, stage1);
}

/*
  Advance solution a_state in time, a_time to a_time + a_dt.
 */
template <class Soln, class Rhs, class EXOP>
void RK4DenseOutput<Soln, Rhs, EXOP>::
advance(Real a_time, Soln& a_state)
{
  CH_TIMERS("RK4DenseOutput::advance");
  CH_assert(isDefined());

  CH_TIMER("RK4DenseOutput::advance - explicit op", t1);

  // Reset the dense output coefs
  if (m_denseOutput)
  {
    m_hasDenseOutput = false;
    for (int icoef=0; icoef < s_nDenseCoefs; ++icoef)
      m_denseCoefs[icoef].zero();
  }

  // Copy a_state into all stages
  for (int stage = 0; stage < s_nStages; stage++)
    m_phi[stage].copy(a_state);

  // Copy a_state values back to itself, to reset any internal state
  // (like fine flux registers for this level)
  a_state.copy(m_phi[0]);

  // For each stage
  for (int stage = 0; stage < s_nStages; stage++)
    {
      Real t = a_time + s_c[stage]*m_dt;
      if (stage > 0)
      {
        // Copy rhs from phi[stage]
        m_rhs.copy(m_phi[stage]);
      }

      // Calculate the operators for this stage
      CH_START(t1);
      m_opEx.explicitOp(m_kE, t, m_phi[stage], stage);
      CH_STOP(t1);
      // Accumulate the explicit op into the stage rhs, final solution
      // NOTE: increment any internal registers
      a_state.increment(m_kE, m_dt*s_b[stage], true);
      for (int k=stage+1; k < s_nStages; ++k)
        m_phi[k].increment(m_kE, m_dt*s_a[k][stage]);

      // This might provide an optimization
      // Accumulate the final solution diff from last stage, explicit op only
      // a_state.increment(m_kE, m_dt*(s_b[stage] - s_aE[s_nStages-1][stage]));

      if (m_denseOutput)
      {
        // pout() << "Stage: " << stage-1 << endl;
        for (int icoef=0; icoef < s_nDenseCoefs; ++icoef)
        {
          m_denseCoefs[icoef].increment(m_kE, m_dt*s_bstar[icoef][stage]);
          
          /*
          LevelData<FArrayBox>& coef = m_denseCoefs[icoef].data();
          DataIterator dit(coef.getBoxes());
          for (dit.begin(); dit.ok(); ++dit)
            pout() << "  Coef[" << icoef+1 << "] = " << coef[dit].min() << endl;
          */
        }
      }
    }

  m_hasDenseOutput = m_denseOutput;
  m_time = a_time;
}

/*
  Return the coefs to interpolate solution, in terms of power of the fraction
  of time between t_old and t_new.
 */
template <class Soln, class Rhs, class EXOP>
void RK4DenseOutput<Soln, Rhs, EXOP>::
denseOutputCoefs(Vector<Rhs*>& a_interpCoefs)
{
  CH_TIMERS("RK4DenseOutput::denseOutputCoefs");

  const int nCoef = s_nDenseCoefs+1;
  CH_assert(m_hasDenseOutput);
  CH_assert(a_interpCoefs.size() == nCoef); 

  for (int icoef=0; icoef < nCoef; ++icoef)
    CH_assert(a_interpCoefs[icoef] != NULL);

  // Copy over the dense coef values

  // First coeficient is just the old state
  a_interpCoefs[0]->copy(m_phi[0]);
  
  // Next coefs are our dense output
  for (int icoef = 1; icoef < nCoef ; ++icoef)
  {
    /*
    LevelData<FArrayBox>& coef = m_denseCoefs[icoef].data();
    for (dit.begin(); dit.ok(); ++dit)
      pout() << "  Coef[" << icoef+1 << "] = " << coef[dit].min() << endl;
    */
    a_interpCoefs[icoef]->copy(m_denseCoefs[icoef-1]);
  }
}


/*
  Static constants for RK4
 */

// Time coefficients for each stage
template <class Soln, class Rhs, class EXOP>
const Real RK4DenseOutput<Soln, Rhs, EXOP>::
s_c[] = { 0.0, 0.5, 0.5, 1.0 };
  
// Stage coefficients - each row is for that stage 
template <class Soln, class Rhs, class EXOP>
const Real RK4DenseOutput<Soln, Rhs, EXOP>::
s_a[][RK4DenseOutput<Soln, Rhs, EXOP>::s_nStages] = {
  {0., 0., 0., 0.},
  {0.5, 0., 0., 0.},
  {0., 0.5, 0., 0.},
  {0., 0., 1.0, 0.},
};

// Final stage coefficients
template <class Soln, class Rhs, class EXOP>
const Real RK4DenseOutput<Soln, Rhs, EXOP>::
s_b[] =
  {0.16666666666666667, 0.33333333333333333, 0.33333333333333333, 0.16666666666666667};

// Coefficients for dense ouput, 4th-order interpolation
template <class Soln, class Rhs, class EXOP>
const Real RK4DenseOutput<Soln, Rhs, EXOP>::
s_bstar[][RK4DenseOutput<Soln, Rhs, EXOP>::s_nStages] = {
  {1.0, 0., 0., 0.},
  {-1.5, 1.0, 1.0, -0.5},
  {0.66666666666666667, -0.66666666666666667, -0.66666666666666667, 0.66666666666666667}
};


#include "NamespaceFooter.H"
#endif 
