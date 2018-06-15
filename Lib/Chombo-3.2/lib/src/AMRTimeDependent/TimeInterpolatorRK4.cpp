#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>

#include "TimeInterpolatorRK4.H"
#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
TimeInterpolatorRK4::TimeInterpolatorRK4()
{
  m_defined = false;
  resetData();
}


//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
TimeInterpolatorRK4::~TimeInterpolatorRK4()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void TimeInterpolatorRK4::define(/// layout at this level
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
                                 const int&                a_ghosts)
{
  // Cache data
  m_refineCoarse = a_refineCoarse;
  m_coarseDomain = coarsen(a_domain, m_refineCoarse);
  m_numStates = a_numStates;
  m_ghosts = a_ghosts;
  m_ghostVect = m_ghosts * IntVect::Unit;
  m_numCoeffs = 4;

  m_coarseLayout = a_coarserDisjointBoxLayout;
  // petermc, 19 Dec 2008:  prevents crash in case of calling define again
  m_thisCoarsenedLayout = DisjointBoxLayout();
  coarsen(m_thisCoarsenedLayout, a_thisDisjointBoxLayout, m_refineCoarse);

  m_rhsCopy.define(m_thisCoarsenedLayout, m_numStates, m_ghostVect);
  m_taylorCoeffs.define(m_thisCoarsenedLayout,
                        m_numCoeffs * m_numStates,
                        m_ghostVect);
  m_diff12.define(m_thisCoarsenedLayout, m_numStates, m_ghostVect);

  m_copier.define(m_coarseLayout, m_thisCoarsenedLayout,
                  m_coarseDomain, m_ghostVect);
  // Everything is defined now.
  m_defined = true;
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK4::setDt(const Real&  a_dt)
{
  resetData();
  m_dt = a_dt;
  m_gotDt = true;

  m_coeffs.resize(m_numCoeffs);
  // setVectorDt uses m_dt;
  // setVectorDt(vec, a, b, c, d) sets vec := [a, b, c, d]*m_dt
  setVectorDt(m_coeffs[0],  0.,     0.,     0.,     0.);
  setVectorDt(m_coeffs[1],  1.,     0.,     0.,     0.);
  setVectorDt(m_coeffs[2], -3./2.,  1.,     1.,    -1./2.);
  setVectorDt(m_coeffs[3],  2./3., -2./3., -2./3.,  2./3.);
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK4::saveInitialSoln(const LevelData<FArrayBox>&   a_soln)
{
  CH_assert(m_defined);
  CH_assert(m_gotDt);
  CH_assert(!m_gotInitialSoln);
  CH_assert(a_soln.nComp() == m_numStates);

  // First zero out taylorFab.
  DataIterator dit = m_taylorCoeffs.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& taylorFab = m_taylorCoeffs[dit];
      taylorFab.setVal(0.);
    }

  // a_soln is on the coarse layout;
  // m_taylorCoeffs is on the coarsened fine layout.
  // Copy from a_soln to first m_numStates components of m_taylorCoeffs.
  const Interval& srcInt = a_soln.interval();
  a_soln.copyTo(srcInt, m_taylorCoeffs, srcInt, m_copier);

  m_gotInitialSoln = true;
}

//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK4::saveRHS(const LevelData<FArrayBox>&   a_rhs)
{
  CH_assert(m_defined);
  CH_assert(m_gotDt);
  CH_assert(m_gotInitialSoln);
  CH_assert(!m_gotFullTaylorPoly);
  CH_assert(m_countRHS >= 0);
  CH_assert(m_countRHS < 4);
  CH_assert(a_rhs.nComp() == m_numStates);

  // a_rhs is on the coarse layout;
  // m_rhsCopy is on the coarsened fine layout.
  a_rhs.copyTo(m_rhsCopy, m_copier);

  DataIterator dit = m_rhsCopy.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& rhsCopyFab = m_rhsCopy[dit];
      FArrayBox& taylorFab = m_taylorCoeffs[dit];

      for (int ind = 0; ind < m_numCoeffs; ind++)
        {
          Real multFactor = m_coeffs[ind][m_countRHS];
          if (multFactor != 0.)
            {
              taylorFab.plus(rhsCopyFab,
                             multFactor, // multiply rhsCopyFab by this
                             0, // start rhsCopyFab component
                             ind * m_numStates, // start taylorFab component
                             m_numStates); // number of components
            }
        }
      // These are for getting to RK4 intermediates:
      // diff12 = m_dt * (rhs[2] - rhs[1])
      if (m_countRHS == 1)
        {
          FArrayBox& diff12Fab = m_diff12[dit];
          diff12Fab.copy(rhsCopyFab);
        }
      if (m_countRHS == 2)
        {
          FArrayBox& diff12Fab = m_diff12[dit];
          diff12Fab -= rhsCopyFab;
          diff12Fab *= -m_dt;
        }
    }
  m_countRHS++;
  if (m_countRHS == 4)
    {
      m_taylorCoeffs.exchange();
      m_gotFullTaylorPoly = true;
    }
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK4::interpolate(/// interpolated solution on this level coarsened
                                      LevelData<FArrayBox>&   a_U,
                                      /// time interpolation coefficient in range [0:1]
                                      const Real&             a_timeInterpCoeff,
                                      /// interval of a_U to fill in
                                      const Interval&         a_intvl)
{
  CH_assert(m_defined);
  CH_assert(m_gotFullTaylorPoly);
  CH_assert(a_U.nComp() == m_numStates);

  LevelData<FArrayBox> UComp;
  aliasLevelData(UComp, &a_U, a_intvl);

  // For i in 0:m_numCoeffs-1,
  // coeffFirst[i] is index of first component of m_taylorCoeffs
  // that corresponds to a coefficient of t^i.
  Vector<int> coeffFirst(m_numCoeffs);
  int intervalLength = a_intvl.size();
  for (int i = 0; i < m_numCoeffs; i++)
    {
      coeffFirst[i] = a_intvl.begin() + i * m_numStates;
    }

  DataIterator dit = UComp.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& UFab = UComp[dit];
      const FArrayBox& taylorFab = m_taylorCoeffs[dit];

      // Evaluate a0 + a1*t + a2*t^2 + a3*t^3
      // as a0 + t * (a1 + t * (a2 + t * a3)).
      UFab.copy(taylorFab, coeffFirst[m_numCoeffs-1], 0, intervalLength);
      for (int ind = m_numCoeffs - 2; ind >=0; ind--)
        {
          UFab *= a_timeInterpCoeff;
          UFab.plus(taylorFab, coeffFirst[ind], 0, intervalLength);
        }
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK4::intermediate(/// intermediate RK4 solution on this level coarsened
                                       LevelData<FArrayBox>&   a_U,
                                       /// time interpolation coefficient in range [0:1]
                                       const Real&             a_timeInterpCoeff,
                                       /// which RK4 stage:  0, 1, 2, 3
                                       const int&              a_stage,
                                       /// interval of a_U to fill in
                                       const Interval&         a_intvl) const
{
  CH_assert(m_defined);
  CH_assert(m_gotFullTaylorPoly);
  CH_assert(a_U.nComp() == m_numStates);
  CH_assert(a_stage >= 0);
  CH_assert(a_stage < 4);

  Real rinv = 1. / Real(m_refineCoarse);
  Vector<Real> intermCoeffs(4);
  // 0 is coefficient of m_taylorCoeffs[0] = Ucoarse(0)
  // 1 is coefficient of m_taylorCoeffs[1] =           K1
  // 2 is coefficient of m_taylorCoeffs[2] = 1/2 * (-3*K1 + 2*K2 + 2*K3 - K4)
  // 3 is coefficient of m_taylorCoeffs[3] = 2/3 * (   K1 -   K2 -   K3 + K4)
  Real diff12Coeffs;
  // coefficient of m_diff12               =              -   K2 +   K3
  switch (a_stage)
    {
    case 0:
      intermCoeffs[0] = 1.;
      intermCoeffs[1] = a_timeInterpCoeff;
      intermCoeffs[2] = a_timeInterpCoeff * a_timeInterpCoeff;
      intermCoeffs[3] = a_timeInterpCoeff * a_timeInterpCoeff  * a_timeInterpCoeff;
      diff12Coeffs = 0.;
      break;
    case 1:
      intermCoeffs[0] = 1.;
      intermCoeffs[1] = 0.5*rinv + a_timeInterpCoeff;
      intermCoeffs[2] = a_timeInterpCoeff * (rinv + a_timeInterpCoeff);
      intermCoeffs[3] = a_timeInterpCoeff * a_timeInterpCoeff * (1.5*rinv + a_timeInterpCoeff);
      diff12Coeffs = 0.;
      break;
    case 2:
      intermCoeffs[0] = 1.;
      intermCoeffs[1] = 0.5*rinv + a_timeInterpCoeff;
      intermCoeffs[2] = 0.5*rinv*rinv + a_timeInterpCoeff * (rinv + a_timeInterpCoeff);
      intermCoeffs[3] = 0.375*rinv*rinv*rinv + a_timeInterpCoeff * (1.5*rinv*rinv + a_timeInterpCoeff * (1.5*rinv + a_timeInterpCoeff));
      diff12Coeffs = -0.25 * rinv * rinv;
      break;
    case 3:
      intermCoeffs[0] = 1.;
      intermCoeffs[1] = rinv + a_timeInterpCoeff;
      intermCoeffs[2] = rinv*rinv + a_timeInterpCoeff * (2.*rinv + a_timeInterpCoeff);
      intermCoeffs[3] = 0.75*rinv*rinv*rinv + a_timeInterpCoeff * (3.*rinv*rinv + a_timeInterpCoeff * (3.*rinv + a_timeInterpCoeff));
      diff12Coeffs = 0.5 * rinv * rinv;
      break;
    default:
      MayDay::Error("TimeInterpolatorRK4::intermediate must have a_stage in range 0:3");
    }
  LevelData<FArrayBox> UComp;
  aliasLevelData(UComp, &a_U, a_intvl);

  // For i in 0:m_numCoeffs-1,
  // coeffFirst[i] is index of first component of m_taylorCoeffs
  // that corresponds to a coefficient of t^i.
  Vector<int> coeffFirst(m_numCoeffs);
  int intervalLength = a_intvl.size();
  for (int i = 0; i < m_numCoeffs; i++)
    {
      coeffFirst[i] = a_intvl.begin() + i * m_numStates;
    }

  DataIterator dit = UComp.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& UFab = UComp[dit];
      const FArrayBox& taylorFab = m_taylorCoeffs[dit];

      // WAS:  Evaluate a0 + a1*t + a2*t^2 + a3*t^3
      // as a0 + t * (a1 + t * (a2 + t * a3)):
      // that is, set UFab to
      // a3, t*a3 + a2, t*(t*a3 + a2) + a1, t*(t*(t*a3 + a2) + a1) * a0.

      // NEW:  Evaluate a0*c0 + a1*c1 + a2*c2 + a3*c3,
      // where c0, c1, c2, c3 are scalars,
      // and c0 = intermCoeffs[0] = 1.
      UFab.copy(taylorFab, coeffFirst[0], 0, intervalLength);
      for (int ind = 1; ind < 4; ind++)
        {
          UFab.plus(taylorFab, intermCoeffs[ind],
                    coeffFirst[ind], 0, intervalLength);
        }
      const FArrayBox& diff12Fab = m_diff12[dit];
      UFab.plus(diff12Fab, diff12Coeffs, 0, 0, intervalLength);
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK4::resetData()
{
  m_gotDt = false;
  m_gotInitialSoln = false;
  m_gotFullTaylorPoly = false;
  m_countRHS = 0;
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK4::setVectorDt(Vector<Real>& a_vec,
                                      Real a_c0, Real a_c1, Real a_c2, Real a_c3)
{
  a_vec.resize(m_numCoeffs);
  a_vec[0] = a_c0 * m_dt;
  a_vec[1] = a_c1 * m_dt;
  a_vec[2] = a_c2 * m_dt;
  a_vec[3] = a_c3 * m_dt;
}

#include "NamespaceFooter.H"
