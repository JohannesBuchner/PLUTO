#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <limits>
#include <functional>
#include <cstring>
using std::endl;

#include "REAL.H"
#include "parstream.H"
#include "SPMD.H"
#include "RootSolver.H"

#include "UsingBaseNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testRootSolver();

/// Global variables for handling output:
static const char *pgmname = "testRootSolver" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if (verbose)
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testRootSolver() ;

  if (ret == 0)
    {
      if (verbose)
        pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

struct CrazyFunc1 : public std::unary_function<Real, Real>
{
  // Use following command (substituting p and b with desired values) in gnuplot
  // to see the function
  //   plot [-1:1] (x**p)*(3.5+0.5*sin(30*x)/x)*erf(b*x)
  CrazyFunc1(const int a_p = 0, const Real a_b = 10., const Real a_root = 0.1)
    :
    m_p(a_p),
    m_b(a_b),
    m_root(a_root)
    { }
  Real operator()(const Real& a_x) const
    {
      const Real x = a_x - m_root;
      return std::pow(x, m_p)*(3.5 + 0.5*sin(30.*x)/x)*erf(m_b*x);
    }
private:
  const int  m_p;
  const Real m_b;
  const Real m_root;
};

// Return true if not within specified tolerance
bool checkRootDefaultTol(const Real a_x0, const Real a_root)
{
  return (Abs(a_x0 - a_root) >
          (((Real)4.0)*RootSolver::RootTr<Real>::eps()*Abs(a_root) +
           RootSolver::RootTr<Real>::tolerance()));
}

int trialRootSolver(const int a_p, const Real a_b, const Real a_root)
{
  int numIter;
  const Real xmin = -1.;
  const Real xmax = 1.;
  CrazyFunc1 func(a_p, a_b, a_root);
  const Real x0 = RootSolver::Brent(numIter, func, xmin, xmax);
  int ier = 0;
  if (checkRootDefaultTol(x0, a_root))
    {
      pout() << indent << pgmname
             << " failed to converge within tolerance with\n"
             << indent << "  p     = " << a_p << endl
             << indent << "  b     = " << a_b << endl;
      pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
      pout().precision(std::numeric_limits<Real>::digits10);
      pout() << indent << "  exact = " << a_root << endl;
      pout() << indent << "  x     = " << x0 << endl;
      pout() << indent << "  f(x)  = " << func(x0) << endl;
      pout().setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      pout().precision(6);
      ier = 1;
    }
  return ier;
}

int testRootSolver()
{
  int status = 0;
  for (int iloc = -1; iloc <= 1; ++iloc)
    {
      const Real root = 0.5*iloc;
      // p = -2
      if (trialRootSolver(-2, 10., root)) status += 1;
      // p = 0, b = 10.
      if (trialRootSolver( 0, 10., root)) status += 1;
      // p = 0, b = 1.
      if (trialRootSolver( 0,  1., root)) status += 1;
      // p = 2
      // This one can fail convergence for double -- retry with 2x iterations
      // if it fails
      {
        const int p = 2;
        const Real b = 10.;
        int numIter;
        const Real xmin = -1.;
        const Real xmax = 1.;
        CrazyFunc1 func(p, b, root);
        const Real x0 = RootSolver::Brent(numIter, func, xmin, xmax);
        if (checkRootDefaultTol(x0, root))
          {
            const Real x0 = RootSolver::Brent(
              numIter, func, xmin, xmax,
              RootSolver::RootTr<Real>::tolerance(), 200);
            if (checkRootDefaultTol(x0, root))
              {
                pout() << indent << pgmname
                       << " failed to converge within tolerance with\n"
                       << indent << "  p     = " << p << endl
                       << indent << "  b     = " << b << endl;
                pout().setf(std::ios_base::scientific,
                            std::ios_base::floatfield);
                pout().precision(std::numeric_limits<Real>::digits10);
                pout() << indent << "  exact = " << root << endl;
                pout() << indent << "  x     = " << x0 << endl;
                pout() << indent << "  f(x)  = " << func(x0) << endl;
                pout().setf(std::ios_base::fmtflags(0),
                            std::ios_base::floatfield);
                pout().precision(6);
                status += 1;
              }
          }
      }
      // Test again p=2 by specifying the precision
      {
        const int p = 2;
        const Real b = 10.;
        int numIter;
        const Real xmin = -1.;
        const Real xmax = 1.;
        const int prec = 4;
        CrazyFunc1 func(p, b, root);
        const Real x0 =
          RootSolver::Brent(numIter, func, xmin, xmax, (Real)prec);
        const Real tol = std::pow(10., -std::abs(prec));
        if (Abs(x0 - root) >
            (((Real)4.0)*RootSolver::RootTr<Real>::eps()*Abs(root) +
             tol))
          {
            pout() << indent << pgmname
                   << " failed to converge within tolerance of " << prec
                   << " significant digits with\n"
                   << indent << "  p     = " << p << endl
                   << indent << "  b     = " << b << endl;
            pout().setf(std::ios_base::scientific, std::ios_base::floatfield);
            pout().precision(std::numeric_limits<Real>::digits10);
            pout() << indent << "  exact = " << root << endl;
            pout() << indent << "  x     = " << x0 << endl;
            pout() << indent << "  f(x)  = " << func(x0) << endl;
            pout().setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
            pout().precision(6);
            status += 1;
          }
      }
    }
  return status;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i] ,"-v" ,3 ) == 0)
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if (strncmp( argv[i] ,"-q" ,3 ) == 0)
            {
              verbose = false ;
              // argv[i] = "" ;
            }
          else
            {
              break ;
            }
        }
    }
  return ;
}

