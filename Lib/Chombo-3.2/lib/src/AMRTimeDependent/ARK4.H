#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ARK4_H_
#define _ARK4_H_

/// 4th-order additive Runge-Kutta algorithm
/**
   This templated class encapsulates
   the fourth-order additive Runge-Kutta method 
   "ARK4(3)6L[2]SA"
   by Kennedy and Carpenter 2003 Appl. Numer. Math. 44: 139-181

   See also section 3 of Zhang, Johansen, and Colella,
   SIAM J. Sci. Comput. 34, pp. B179-B201.
*/

// Forward declaration
class dXConcept;

// \class XConcept
// Concept of class X, the type of your state data
class XConcept
{
  // allocate space for this
  void define(const XConcept& a_state);

  // copy data to this
  void copy(const XConcept& a_state);

  // increment this by a_factor * a_increment
  void increment(const dXConcept& a_increment, Real a_factor = 1.);

  ~XConcept() {};
};

// \class dXConcept
// Concept of class dX, same type as the operator evaluations
class dXConcept
{
  // allocate space for this
  void define(const XConcept& a_state);

  // initialize with given X
  void init(const XConcept& a_state);

  // increment this by a_factor * a_increment
  void increment(const dXConcept& a_increment, Real a_factor = 1.);

  ~dXConcept() {};
};

// \class FEConcept
// Concept of class FE, which defines your explicit operator
class FEConcept
{
  // allocate space and define infrastructure for this
  void define(const XConcept& a_state, Real a_dt);

  // reset the timestep
  void resetDt(Real a_dt);

  // apply the operator, including multiplication by timestep
  void operator()(dXConcept& a_result, Real a_time, const XConcept& a_state);

  ~FEConcept() {};
};

// \class FIConcept
// Concept of class FI, which defines your implicit operator
class FIConcept
{
  // allocate space and define infrastructure for this
  void define(const XConcept& a_state, Real a_dt, Real a_gamma);

  // reset the timestep
  void resetDt(Real a_dt);

  // apply the operator, including multiplication by timestep
  void operator()(dXConcept& a_result, Real a_time, const XConcept& a_state);

  // solve for a_soln in (I - s_gamma * m_dt * FI(a_time))(a_soln) = a_rhs
  void solve(XConcept& a_soln, const dXConcept& a_rhs, Real a_time);

  ~FIConcept() {};
};

template <class X, class FI, class FE, class dX>
class ARK4
{
public:

  ARK4<X, FI, FE, dX>() { m_isDefined = false; }

  // This must be called first.
  void define(const X& a_state, Real a_dt);

  // Advance one step.
  void advance(Real a_time, X& a_state);

  // Reset the timestep.
  void resetDt(Real a_dt);

  bool isDefined() const { return m_isDefined; }

  /// Runge-Kutta coefficients
  static const int  s_nStages = 6;
  static const Real s_gamma;
  static const Real s_aE[s_nStages][s_nStages];
  static const Real s_aI[s_nStages][s_nStages];
  static const Real s_b[s_nStages];
  static const Real s_c[s_nStages];

protected:

  bool m_isDefined;
  Real m_dt;
  X m_phi[s_nStages];
  dX m_kE[s_nStages-1];
  dX m_kI[s_nStages-1];
  dX m_kEfinal;
  dX m_rhs;
  FI m_fI;
  FE m_fE;

private:

  // Make sure the concepts match the API
  static ARK4<XConcept, FIConcept, FEConcept, dXConcept> testConcept;
};

//==============================================

template <class X, class FI, class FE, class dX>
void ARK4<X, FI, FE, dX>::define(const X& a_state, Real a_dt)
{
  m_dt = a_dt;
  // define X
  for (int stage = 0; stage < s_nStages; stage++)
    {
      m_phi[stage].define(a_state);
    }
  // define dX
  for (int stage = 0; stage < s_nStages-1; stage++)
    {
      m_kE[stage].define(a_state);
      m_kI[stage].define(a_state);
    }
  m_kEfinal.define(a_state);
  m_rhs.define(a_state);
  // define fI
  m_fI.define(a_state, m_dt, s_gamma);
  // define fE
  m_fE.define(a_state, m_dt);
  m_isDefined = true;
}

/*
  Reset the timestep.
 */
template <class X, class FI, class FE, class dX>
void ARK4<X, FI, FE, dX>::resetDt(Real a_dt)
{
  CH_assert(isDefined());

  // Only update everything if dt has changed
  Real reltol = 1e-10;
  if (Abs(m_dt - a_dt) > m_dt*reltol)
  {
    m_dt = a_dt;
    m_fI.resetDt(m_dt);
    m_fE.resetDt(m_dt);
  }
}

/*
  Advance solution a_state in time, a_time to a_time + a_dt.
 */
template <class X, class FI, class FE, class dX>
void ARK4<X, FI, FE, dX>::advance(Real a_time, X& a_state)
{
  CH_assert(isDefined());
  // Set intermediate times t from a_time and m_dt.
  Real t[s_nStages];
  for (int stage = 0; stage < s_nStages; stage++)
    {
      t[stage] = a_time + s_c[stage] * m_dt;
    }

  // Set m_phi[0] := a_state.
  m_phi[0].copy(a_state);
  for (int stage = 1; stage < s_nStages; stage++)
    {
      m_fE(m_kE[stage-1], t[stage-1], m_phi[stage-1]);
      m_fI(m_kI[stage-1], t[stage-1], m_phi[stage-1]);

      // Set m_rhs := a_state +
      //    sum_{j=0:stage-1} ( s_aE[stage][j] * FE(m_phi[j], t[j]) +
      //                        s_aI[stage][j] * FI(m_phi[j]) ).
      m_rhs.init(a_state);
      for (int j = 0; j < stage; j++)
        {
          m_rhs.increment(m_kE[j], s_aE[stage][j]);
          m_rhs.increment(m_kI[j], s_aI[stage][j]);
        }
      // Solve for m_phi[stage] in
      // (I - s_gamma * m_dt * FI) (m_phi[stage]) = m_rhs.
      m_phi[stage].define(a_state);
      m_fI.solve(m_phi[stage], m_rhs, t[stage]);
    }

  // Set a_state := m_phi[s_nStages-1] + 
  //    sum_{j=0:s_nStages-1} ( (s_b[j] - s_aE[s_nStages-1][j]) *
  //                             FE(m_phi[j], t[s_nStages-1]) ).
  // Note these m_fE evaluations are for m_phi[j] at time t[s_nStages-1],
  // not time t[j], so you cannot reuse m_kE[j].
  const int stageFinal = s_nStages-1;
  Real tFinal = t[stageFinal];
  // Calculate the final stage explicit operator
  m_fE(m_kEfinal, tFinal, m_phi[stageFinal]);
  a_state.copy(m_phi[stageFinal]);
  for (int j = 0; j < s_nStages-1; j++)
    {
      a_state.increment(m_kE[j], s_b[j] - s_aE[stageFinal][j]);
    }
  // Add in the final stage explicit operator
  a_state.increment(m_kEfinal, s_b[stageFinal] /* - s_aE, which is 0 */);
}

template<class X, class FI, class FE, class dX>
const Real ARK4<X, FI, FE, dX>::s_gamma = 0.25;

// Time coefficients for each stage
template<class X, class FI, class FE, class dX>
const Real ARK4<X, FI, FE, dX>::s_c[] = { 0.0, 0.5, 0.332, 0.62, 0.85, 1.0 };
  
// Stage coefficients - each row is for that stage 
template<class X, class FI, class FE, class dX>
const Real ARK4<X, FI, FE, dX>::s_aE[][ARK4<X, FI, FE, dX>::s_nStages] = {
  {0., 0., 0., 0., 0., 0.},
  {0.5, 0., 0., 0., 0., 0.},
  {0.221776, 0.110224, 0., 0., 0., 0.},
  {-0.04884659515311857, -0.17772065232640102, 0.8465672474795197, 0., 0., 0.},
  {-0.15541685842491548, -0.3567050098221991, 1.0587258798684427, 0.30339598837867193, 0., 0.},
  { 0.2014243506726763, 0.008742057842904185, 0.15993995707168115, 0.4038290605220775, 0.22606457389066084, 0.}
};

template<class X, class FI, class FE, class dX>
const Real ARK4<X, FI, FE, dX>::s_aI[][ARK4<X, FI, FE, dX>::s_nStages] = {
  {0., 0., 0., 0., 0., 0.},
  {0.25, 0.25, 0., 0., 0., 0.},
  {0.137776, -0.055776, 0.25, 0., 0., 0.},
  {0.14463686602698217, -0.22393190761334475, 0.4492950415863626, 0.25, 0., 0.},
  {0.09825878328356477, -0.5915442428196704, 0.8101210538282996, 0.283164405707806, 0.25, 0.},
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25}
};

template<class X, class FI, class FE, class dX>
const Real ARK4<X, FI, FE, dX>::s_b[] =
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25};

#endif 
