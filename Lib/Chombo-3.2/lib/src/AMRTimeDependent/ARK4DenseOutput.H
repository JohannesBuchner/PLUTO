#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ARK4DENSEOUTPUT_H_
#define _ARK4DENSEOUTPUT_H_


#include "NamespaceHeader.H"

/// This is a more flexible, reduced memory version
/// 4th-order additive Runge-Kutta algorithm,
/// that also includes dense output coefficients
/**
   This templated class encapsulates
   the fourth-order additive Runge-Kutta method 
   "ARK4(3)6L[2]SA"
   by Kennedy and Carpenter 2003 Appl. Numer. Math. 44: 139-181

   See also section 3 of Zhang, Johansen, and Colella,
   SIAM J. Sci. Comput. 34, pp. B179-B201.
*/

// Soln is the type that operators are applied to
// Rhs is the type that operator produce
// IMEXOp wraps the explicit and implicit operators
template <class Soln, class Rhs, class IMEXOp>
class ARK4DenseOutput
{
public:

  ARK4DenseOutput<Soln, Rhs, IMEXOp>() { m_isDefined = false; }

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
  IMEXOp& getImExOp();

  bool isDefined() const { return m_isDefined; }

  bool hasDenseOutput() const { return m_hasDenseOutput; }

  /// Runge-Kutta coefficients
  static const int  s_nStages = 6;
  static const Real s_aIdiag;
  static const Real s_c[s_nStages];
  static const Real s_aE[s_nStages][s_nStages];
  static const Real s_aI[s_nStages][s_nStages];
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
  Rhs m_kI;
  IMEXOp m_opImEx;

private:

};

//==============================================

template <class Soln, class Rhs, class IMEXOp>
void ARK4DenseOutput<Soln, Rhs, IMEXOp>::
define(const Soln& a_state, Real a_dt, bool a_denseOutput)
{
  CH_TIMERS("ARK4DenseOutput::define");
  m_dt = a_dt;
  m_denseOutput = a_denseOutput;
  // define Soln types
  for (int stage = 0; stage < s_nStages; stage++)
    m_phi[stage].define(a_state);

  // define Rhs types
  m_kE.define(a_state);
  m_kI.define(a_state);
  m_rhs.define(a_state);

  // define opImEx
  m_opImEx.define(a_state, m_dt, s_aIdiag);

  // if dense output is requested, need more storage
  for (int coef = 0; m_denseOutput && (coef < s_nDenseCoefs); coef++)
    m_denseCoefs[coef].define(a_state);

  m_hasDenseOutput = false;
  m_isDefined = true;
}

/*
  Get a reference to the implicit-explicit operator
 */
template <class Soln, class Rhs, class IMEXOp>
IMEXOp& ARK4DenseOutput<Soln, Rhs, IMEXOp>::
getImExOp()
{
  return m_opImEx;
}

/*
  Reset the timestep.
 */
template <class Soln, class Rhs, class IMEXOp>
void ARK4DenseOutput<Soln, Rhs, IMEXOp>::
resetDt(Real a_dt)
{
  CH_assert(isDefined());

  // Only update everything if dt has changed
  Real reltol = 1e-14;
  if (Abs(m_dt - a_dt) > m_dt*reltol)
  {
    m_dt = a_dt;
    m_hasDenseOutput = false;
    m_opImEx.resetDt(m_dt);
  }
}

/*
  Set whether the step starts at 0. and/or ends at 1.
 */
template <class Soln, class Rhs, class IMEXOp>
void ARK4DenseOutput<Soln, Rhs, IMEXOp>::
start0end1(bool a_start0, bool a_end1)
{
  CH_assert(isDefined());

  int stage0 = -1; // if no stage is at time 0.
  if (a_start0) stage0 = 0;

  int stage1 = -1; // if no stage is at time 1.
  if (a_end1) stage1 = s_nStages-1;

  m_opImEx.stage0stage1(stage0, stage1);
}

/*
  Advance solution a_state in time, a_time to a_time + a_dt.
 */
template <class Soln, class Rhs, class IMEXOp>
void ARK4DenseOutput<Soln, Rhs, IMEXOp>::
advance(Real a_time, Soln& a_state)
{
  CH_TIMERS("ARK4DenseOutput::advance");
  CH_assert(isDefined());

  CH_TIMER("ARK4DenseOutput::advance - explicit op", t1);
  CH_TIMER("ARK4DenseOutput::advance - implicit op", t2);
  CH_TIMER("ARK4DenseOutput::advance - implicit solve", t3);

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
        // Do the solve - copy rhs from phi[stage]
        m_rhs.copy(m_phi[stage]);
        // Solve for m_phi[stage] in
        // (I - s_aIdiag * m_dt * FI) (m_phi[stage]) = m_rhs.
        CH_START(t3);
        m_opImEx.solve(m_phi[stage], m_rhs, t, stage);
        CH_STOP(t3);
      }

      // Calculate the operators for this stage
      CH_START(t1);
      m_opImEx.explicitOp(m_kE, t, m_phi[stage], stage);
      CH_STOP(t1);
      // Accumulate the explicit op into the stage rhs, final solution
      // NOTE: increment any internal registers
      a_state.increment(m_kE, m_dt*s_b[stage], true);
      for (int k=stage+1; k < s_nStages; ++k)
        m_phi[k].increment(m_kE, m_dt*s_aE[k][stage]);

      CH_START(t2);
      m_opImEx.implicitOp(m_kI, t, m_phi[stage], stage);
      CH_STOP(t2);

      // Now accumulate the implicit op into the stage rhs, final solution
      a_state.increment(m_kI, m_dt*s_b[stage], true);
      for (int k=stage+1; k < s_nStages; ++k)
        m_phi[k].increment(m_kI, m_dt*s_aI[k][stage]);

      // This might provide an optimization
      // Accumulate the final solution diff from last stage, explicit op only
      // a_state.increment(m_kE, m_dt*(s_b[stage] - s_aE[s_nStages-1][stage]));

      if (m_denseOutput)
      {
        // pout() << "Stage: " << stage-1 << endl;
        for (int icoef=0; icoef < s_nDenseCoefs; ++icoef)
        {
          m_denseCoefs[icoef].increment(m_kE, m_dt*s_bstar[icoef][stage]);
          m_denseCoefs[icoef].increment(m_kI, m_dt*s_bstar[icoef][stage]);
          
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
template <class Soln, class Rhs, class IMEXOp>
void ARK4DenseOutput<Soln, Rhs, IMEXOp>::
denseOutputCoefs(Vector<Rhs*>& a_interpCoefs)
{
  CH_TIMERS("ARK4DenseOutput::denseOutputCoefs");

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
  Static constants for ARK4
 */

template <class Soln, class Rhs, class IMEXOp>
const Real ARK4DenseOutput<Soln, Rhs, IMEXOp>::
s_aIdiag = 0.25;

// Time coefficients for each stage
template <class Soln, class Rhs, class IMEXOp>
const Real ARK4DenseOutput<Soln, Rhs, IMEXOp>::
s_c[] = { 0.0, 0.5, 0.332, 0.62, 0.85, 1.0 };
  
// Stage coefficients - each row is for that stage 
template <class Soln, class Rhs, class IMEXOp>
const Real ARK4DenseOutput<Soln, Rhs, IMEXOp>::
s_aE[][ARK4DenseOutput<Soln, Rhs, IMEXOp>::s_nStages] = {
  {0., 0., 0., 0., 0., 0.},
  {0.5, 0., 0., 0., 0., 0.},
  {0.221776, 0.110224, 0., 0., 0., 0.},
  {-0.04884659515311857, -0.17772065232640102, 0.8465672474795197, 0., 0., 0.},
  {-0.15541685842491548, -0.3567050098221991, 1.0587258798684427, 0.30339598837867193, 0., 0.},
  { 0.2014243506726763, 0.008742057842904185, 0.15993995707168115, 0.4038290605220775, 0.22606457389066084, 0.}
};

// Implicit stage coefficients
template <class Soln, class Rhs, class IMEXOp>
const Real ARK4DenseOutput<Soln, Rhs, IMEXOp>::
s_aI[][ARK4DenseOutput<Soln, Rhs, IMEXOp>::s_nStages] = {
  {0., 0., 0., 0., 0., 0.},
  {0.25, 0.25, 0., 0., 0., 0.},
  {0.137776, -0.055776, 0.25, 0., 0., 0.},
  {0.14463686602698217, -0.22393190761334475, 0.4492950415863626, 0.25, 0., 0.},
  {0.09825878328356477, -0.5915442428196704, 0.8101210538282996, 0.283164405707806, 0.25, 0.},
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25}
};

// Final stage coefficients
template <class Soln, class Rhs, class IMEXOp>
const Real ARK4DenseOutput<Soln, Rhs, IMEXOp>::
s_b[] =
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25};

// Coefficients for dense ouput, 4th-order interpolation
template <class Soln, class Rhs, class IMEXOp>
const Real ARK4DenseOutput<Soln, Rhs, IMEXOp>::
s_bstar[][ARK4DenseOutput<Soln, Rhs, IMEXOp>::s_nStages] = {
  {0.961753400252887, 0., 0.787405595186356, -2.74544192086633, 3.70351728061223, -1.70723435518514},
  {-1.76418754019038, 0., -0.774504669155511, 9.64023584441292, -12.544886411271, 5.44334277620397},
  {0.960350435099165, 0., 0.173858014493155, -6.21422862823726, 8.56612859966376, -3.48610842101883}
};

#include "NamespaceFooter.H"
#endif 
