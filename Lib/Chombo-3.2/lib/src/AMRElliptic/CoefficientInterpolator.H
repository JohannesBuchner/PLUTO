#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COEFFICIENTINTERPOLATOR_H__
#define _COEFFICIENTINTERPOLATOR_H__

#include "LevelData.H"
#include "DataIterator.H"
#include "DebugOut.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

// Jeffrey Johnson (jnjohnson@lbl.gov)

//! \class CoefficientInterpolator
//! This abstract base class provides an interface to time-dependent spatial coefficient
//! data for various partial differential equations. One obtains constant access to
//! coefficient data by specifying the desired time centering. Subclasses of this
//! base class define how that coefficient data is computed.
//! \tparam LevelDataType The LevelData type used to store coefficient data.
template <typename LevelData_, typename SolutionLevelData_ = LevelData_>
class CoefficientInterpolator
{
  public:

  typedef LevelData_ LevelDataType;
  typedef SolutionLevelData_ SolutionLevelDataType;

  //! Base class constructor. Called by every subclass.
  //! \param a_numComps The number of components in the interpolated coefficient.
  explicit CoefficientInterpolator(int a_numComps);

  //! Destructor.
  virtual ~CoefficientInterpolator();

  //! Returns the number of components in the interpolated coefficient.
  int numComps() const;

  //! Interpolate the coefficient at the given time, placing the result in the given
  //! LevelData object. This method must be overridden by subclasses.
  //! \param a_result The LevelData object that will store the result.
  //! \param a_time The time at which the coefficient is to be evaluated.
  virtual void interpolate(LevelDataType& a_result,
                           Real a_time);

  //! Interpolate the coefficient at the given time, placing the result in the given
  //! LevelData object. This method assumes that the coefficient depends upon the
  //! solution to the partial differential equation in question, so the solution
  //! is passed into it as an argument.
  //! \param a_result The LevelData object that will store the result.
  //! \param a_solution The solution to the partial differential equation.
  //! \param a_time The time at which the coefficient is to be evaluated.
  virtual void interpolate(LevelDataType& a_result,
                           const SolutionLevelDataType& a_solution,
                           Real a_time);

  //! Returns true if the coefficient depends on the solution to the
  //! partial differential equation (rendering it nonlinear), false
  //! otherwise. By default, the coefficient is assumed not to depend
  //! upon the solution.
  virtual bool dependsUponSolution() const;

  //! Computes the derivative of the coefficient with respect to the
  //! solution at the desired time. By default, this sets \a a_deriv to 0.
  //! \param a_prime The coefficient derivative data will be stored here.
  //! \param a_solution The solution to the partial differential equation.
  //! \param a_time The time at which to compute the coefficient data.
  virtual void interpolatePrime(LevelDataType& a_prime,
                                const SolutionLevelDataType& a_solution,
                                Real a_time);

  //! This virtual void method performs the iterative nonlinear solve
  //! \f$A(\phi) \phi - f(\vec{x}) = 0\f$ for \f$\phi\f$.
  //! \param a_phi The solution to the equation, \f$\phi\f$, will be stored here.
  //! \param a_f The term \f$f(\vec{x})\f$ in the equation.
  //! \param a_time The time at which the equation is solved.
  //! \param a_phi0 The initial estimate for \f$\phi\f$.
  //! \param a_tolerance The threshold for the error in the equation that
  //!                    dictates when iteration should cease.
  virtual void solve(SolutionLevelDataType& a_phi,
                     const SolutionLevelDataType& a_f,
                     Real a_time,
                     const SolutionLevelDataType& a_phi0,
                     Real a_tolerance);

  //! This helper method solves the nonlinear equation
  //! \f$A(\phi) \phi - f(\vec{x}) = 0\f$ for \f$\phi\f$ using Newton-Raphson
  //! iteration. It can be used by subclasses to implement the solve()
  //! virtual method.
  //! \param a_phi The solution to the equation, \f$\phi\f$, will be stored here.
  //! \param a_f The term \f$f(\vec{x})\f$ in the equation.
  //! \param a_time The time at which the equation is solved.
  //! \param a_phi0 The initial estimate for \f$\phi\f$.
  //! \param a_tolerance The threshold for the error in the equation that
  //!                    dictates when iteration should cease.
  void NewtonRaphson(SolutionLevelDataType& a_phi,
                     const SolutionLevelDataType& a_f,
                     Real a_time,
                     const SolutionLevelDataType& a_phi0,
                     Real a_tolerance);

  private:

  // Number of components.
  int m_numComps;

  // Prevents tail-chasing.
  bool m_inCall;

  // Disallowed operators.
  CoefficientInterpolator();
  CoefficientInterpolator(const CoefficientInterpolator&);
  CoefficientInterpolator& operator=(const CoefficientInterpolator&);
};

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
CoefficientInterpolator(int a_numComps)
  :m_numComps(a_numComps),
   m_inCall(false)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
~CoefficientInterpolator()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
int
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
numComps() const
{
  return m_numComps;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
bool
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
dependsUponSolution() const
{
  return false;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
void
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
interpolate(LevelData_& a_result,
            Real a_time)
{
  if (dependsUponSolution())
  {
    MayDay::Error("The solution must be passed to this interpolator!");
  }
  if (m_inCall)
  {
    MayDay::Error("Neither the solution-independent nor -dependent interpolate method has\n"
                  "been defined for this subclass.");
  }

  // It doesn't matter what we pass in for the solution.
  m_inCall = true;
  SolutionLevelData_ phoneyBaloney(a_result.disjointBoxLayout(), 1, a_result.ghostVect());
  interpolate(a_result, phoneyBaloney, a_time);
  m_inCall = false;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
void
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
interpolate(LevelData_& a_result,
            const SolutionLevelData_& a_solution,
            Real a_time)
{
  if (m_inCall)
  {
    MayDay::Error("Neither the solution-independent nor -dependent interpolate method has\n"
                  "been defined for this subclass.");
  }
  if (!dependsUponSolution())
  {
    m_inCall = true;
    interpolate(a_result, a_time);
    m_inCall = false;
  }
  else
  {
    MayDay::Error("The solution-dependent interpolate method has not been defined for this subclass.");
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
void
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
interpolatePrime(LevelData_& a_prime,
                 const SolutionLevelData_& a_solution,
                 Real a_time)
{
  if (!dependsUponSolution())
  {
    for (DataIterator dit = a_prime.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      a_prime[dit()].setVal(0.0);
  }
  else
    MayDay::Error("The derivative w.r.t. the solution has not been implemented for this subclass.");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
void
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
solve(SolutionLevelData_& a_phi,
      const SolutionLevelData_& a_f,
      Real a_time,
      const SolutionLevelData_& a_phi0,
      Real a_tolerance)
{
  MayDay::Error("CoefficientInterpolator::solve() must be implemented for this subclass!");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
template <typename LevelData_, typename SolutionLevelData_>
void
CoefficientInterpolator<LevelData_, SolutionLevelData_>::
NewtonRaphson(SolutionLevelData_& a_phi,
              const SolutionLevelData_& a_f,
              Real a_time,
              const SolutionLevelData_& a_phi0,
              Real a_tolerance)
{
  // How bad is the initial estimate?
  Real maxError = -FLT_MAX;
  LevelData_ F(a_phi0.disjointBoxLayout(), a_phi0.nComp()),
             A(a_phi0.disjointBoxLayout(), a_phi0.nComp());
  interpolate(A, a_phi0, a_time);
  for (DataIterator dit = a_phi0.dataIterator(); dit.ok(); ++dit)
  {
    a_phi[dit()].copy(a_phi0[dit()]); // phi <- phi0.
    F[dit()].copy(A[dit()]);
    F[dit()] *= a_phi[dit()];
    F[dit()] -= a_f[dit()];
    maxError = Max(Abs(F[dit()].max()), Abs(F[dit()].min()));
  }
  if (maxError < a_tolerance)
    return; // Actually, our initial estimate was fine!

  LevelData_ Aprime(a_phi.disjointBoxLayout(), a_phi.nComp());
  while (maxError > a_tolerance)
  {
    // Compute the correction term F/F', apply it to phi, and recompute F.
    // FIXME: This is slow and probably should be moved to Fortran.
    interpolatePrime(Aprime, a_phi, a_time);
    for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      Box box = a_phi.disjointBoxLayout().get(dit());
      const FArrayBox& a = A[dit()];
      const FArrayBox& aprime = Aprime[dit()];
      const FArrayBox& f = F[dit()];
      FArrayBox& phi = a_phi[dit()];
      for (BoxIterator bit(box); bit.ok(); ++bit)
      {
        IntVect i = bit();
        for (int n = 0; n < F.nComp(); ++n)
        {
          if (Abs(f(i, n)) > a_tolerance)
          {
            Real fprime = aprime(i, n) * phi(i, n) + a(i, n);
            phi(i, n) -= f(i, n) / fprime;
          }
        }
      }
    }

    // Recompute F with the corrected value of phi.
    interpolate(A, a_phi, a_time);
    for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
    {
      // Recompute F.
      F[dit()].copy(A[dit()]);
      F[dit()] *= a_phi[dit()];
      F[dit()] -= a_f[dit()];

      // Measure the error again.
      maxError = Max(Abs(F[dit()].max()), Abs(F[dit()].min()));
    }
  }
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"

#endif
