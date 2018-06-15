#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRTGA_H_
#define _AMRTGA_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "AMRMultiGrid.H"
#include "NamespaceHeader.H"

//! \class TGAHelmOp
//! This operator is meant to represent the general helmholtz operator that
//! AMRTGA will be solving.  This operator is of the form
//! alpha A(x) phi + beta div B(x) grad phi = rho.
//! AMRTGA needs to reset the constants alpha and beta often.
//! \tparam T The LevelData class that holds solution data for the operator.
template <class T>
class TGAHelmOp : public AMRLevelOp<T>
{
  public:

  //! Base class default constructor. This constructs a time-independent
  //! operator.
  TGAHelmOp()
    :m_isTimeDependent(false),
     m_phonyIdentityCoef(NULL)
  {
  }

  //! Base class constructor which specifies explicitly whether this
  //! operator is time-dependent.
  //! \param a_isTimeDependent If set to true, this Helmholtz operator will
  //!                          be treated as a time-dependent operator.
  explicit TGAHelmOp(bool a_isTimeDependent)
    :m_isTimeDependent(a_isTimeDependent),
     m_phonyIdentityCoef(NULL)
  {
  }

  //! Destructor.
  virtual ~TGAHelmOp()
  {
  }

  //! Sets the scaling constants in the Helmholtz equation.
  //! \param a_alpha The scaling constant that multiplies the identity term.
  //! \param a_beta The scaling constant that multiplies the derivative term.
  virtual void setAlphaAndBeta(const Real& a_alpha,
                               const Real& a_beta) = 0;

  //! Scales the right hand side of the Helmholtz equation by the
  //! identity term in the operator. If you are solving
  //! rho(x) dphi/dt = L(phi), this means multiply by rho (or kappa * rho in
  //! embedded boundary calculations.
  //! \param a_rhs The right hand side of the equation to be scaled.
  //! \param a_kappaWeighted If set to true, \a a_rhs will be scaled by the
  //!                        volume fraction in addition to the identity term.
  virtual void diagonalScale(T& a_rhs, bool a_kappaWeighted)
  {
    MayDay::Error("Two argument version of diagonalScale called but not implemented!");
  }

  /// for eb only.   kappa weight the rhs but do not multiply by the identity coefficient
  virtual void kappaScale(T& a_rhs)
  {
    // default implementation does nothing because only relevant in eb.
  }

  //! Scales the right hand side of the Helmholtz equation by the
  //! identity term in the operator. This version assumes volume fraction
  //! weighting.
  //! \param a_rhs The right hand side of the equation to be scaled.
  virtual void diagonalScale(T& a_rhs)
  {
    diagonalScale(a_rhs, true);
  }

  //! Divides the right hand side of the Helmholtz equation by the
  //! identity coefficient rho(x) in the equation rho(x) dphi/dt = L(phi).
  //! \param a_rhs The right hand side of the equation to be scaled.
  virtual void divideByIdentityCoef(T& a_rhs) = 0;

  //! Apply the differential operator without any boundary or coarse-fine
  //! boundary conditions and no finer level
  //! \param a_ans The result of the application of the operator to \a a_phi.
  //! \param a_phi The data to which the operator is applied.
  virtual void applyOpNoBoundary(T& a_ans, const T& a_phi) = 0;

  //! Sets the time-dependent state of the operator. The default implementation
  //! does nothing and is appropriate for time-independent operators.
  //! \param a_time The time to be used to update the time-dependent operator.
  virtual void setTime(Real a_time)
  {
    ;// most operators are time-independent.
  }

  //! Sets the time-dependent state of the operator. This version of setTime
  //! allows one to linearly interpolate coefficients across an integration
  //! step, since it accepts arguments that define where in the step it is
  //! to be updated. The default implementation calls
  //! setTime(\a a_oldTime + a_mu * a_dt).
  //! \param a_oldTime The time at the beginning of the current step.
  //! \param a_mu The fraction of the current step that has elapsed.
  //! \param a_dt The size of the current step.
  virtual void setTime(Real a_oldTime, Real a_mu, Real a_dt)
  {
    setTime(a_oldTime + a_mu * a_dt);
  }

  //! Returns true if the operator is time-dependent, false otherwise.
  bool isTimeDependent() const
  {
    return m_isTimeDependent;
  }

  //! Allows access to the identity coefficient data for the operator.
  virtual T& identityCoef()
  {
     MayDay::Error("Access to identity coefficient not allowed for this operator");
     // Yes, this is a dereference to a null pointer. No, it shouldn't be
     // called. I'm putting this here to make the compiler shaddup. -JNJ
     return *m_phonyIdentityCoef;
  }

  private:

  // Disallowed operations.
  TGAHelmOp(const TGAHelmOp&);
  TGAHelmOp& operator=(const TGAHelmOp&);

  //! This flag is set if the Helmholtz operator's coefficiens are
  //! time-dependent.
  bool m_isTimeDependent;

  // Phony identity coefficient data to avoid compiler warnings. This
  // should go away soon.
  T* m_phonyIdentityCoef;
};

//! \class LevelTGAHelmOp
//! This subclass of TGAHelmOp proves additional functionality for computing
//! fluxes at refinement level boundaries.
//! \tparam T The LevelData class that holds solution data.
//! \tparam TFlux The FArrayBox class that holds flux data within a LevelData.
template <class T, class TFlux>
class LevelTGAHelmOp : public TGAHelmOp< T >
{
  public:

  //! Base class default constructor. This constructs a time-independent
  //! LevelTGAHelmOp.
  LevelTGAHelmOp()
    :TGAHelmOp<T>::TGAHelmOp(false)
  {
  }

  //! Base class constructor which specifies explicitly whether this
  //! operator is time-dependent.
  //! \param a_isTimeDependent If set to true, this Helmholtz operator will
  //!                          be treated as a time-dependent operator.
  LevelTGAHelmOp(bool a_isTimeIndependent)
    :TGAHelmOp<T>::TGAHelmOp(a_isTimeIndependent)
  {
  }

  //! Destructor.
  virtual ~LevelTGAHelmOp()
  {
  }

  ///  These functions are part of the LevelTGA interface......

  ///
  virtual void fillGrad(const T& a_phi) = 0;

  ///
  virtual void getFlux(TFlux&           a_flux,
                       const T&         a_data,
                       const Box&       a_grid,
                       const DataIndex& a_dit,
                       Real             a_scale) = 0;

};

///
/**
   Template implementation of the TGA algorithm
   to solve the heat equation.
 **/
template <class T>
class AMRTGA
{
public:

  ///
  TGAHelmOp<T>* newOp(const ProblemDomain&                                  a_indexSpace,
                      const AMRLevelOpFactory<T>&                           a_opFact)
  {
    AMRLevelOpFactory<T>& opFact = (AMRLevelOpFactory<T>&) a_opFact;
    TGAHelmOp<T>* retval = (TGAHelmOp<T>*) opFact.AMRnewOp(a_indexSpace);
    return retval;
  }

  ///
  ~AMRTGA();

  ///
  /**
   **/
  AMRTGA(const RefCountedPtr<AMRMultiGrid<T> >     &           a_solver,
         const AMRLevelOpFactory<T>&                           a_factory,
         const ProblemDomain&                                  a_level0Domain,
         const Vector<int>&                                    a_refRat,
         int a_numLevels = -1,
         int a_verbosity = 0);

  ///
  /**
     This advances a parabolic pde from a_phiOld to a_phiNew using TGA
     on a non-moving domain with source term a_source
     If your operator is time-dependent, be sure to send the time.
  **/
  void oneStep(Vector<T*>&             a_phiNew,
               Vector<T*>&             a_phiOld,
               Vector<T*>&             a_source,
               const Real&             a_dt,
               int                     a_lbase,
               int                     a_lmax,
               Real                    a_timeOld = 0);

  ///
  void resetAlphaAndBeta(const Real& a_alpha,
                         const Real& a_beta);

  void setTime(Real time);

protected:
  void solveHelm(Vector<T*>&       a_ans,
                 Vector<T*>&       a_rhs,
                 int               a_lbase,
                 int               a_lmax,
                 Real              a_mu,
                 Real              a_dt);

  void applyHelm(Vector<T*>&       a_ans,
                 Vector<T*>&       a_source,
                 int               a_lbase,
                 int               a_lmax,
                 Real              a_mu,
                 Real              a_dt,
                 bool              a_homogeneousBC);

  void setMu();

  void createData(Vector<T* >&       a_source,
                  int                a_lbase,
                  int                a_lmax);

private:
  //You do not own these operators!!  don't delete it.   the memory is
  //owned by the solver
  Vector<TGAHelmOp<T> * >            m_ops;
  ProblemDomain                                  m_level0Domain;
  Vector<int>                                    m_refRat;
  RefCountedPtr<AMRMultiGrid<T> >                m_solver;
  Real m_mu1, m_mu2, m_mu3, m_mu4;
  int m_verbosity, m_numLevels;
  bool                                           m_dataCreated;
  Vector<T*>                                     m_rhst;
  Vector<T*>                                     m_srct;

  //copy constructor and operator= disallowed for all the usual reasons
  AMRTGA(const AMRTGA<T>& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const AMRTGA<T>& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  /// weak construction is bad.   Ref Counted pointers are your friends.
  AMRTGA()
  {
    MayDay::Error("invalid operator");
  }
};

/*****/
template <class T>
void AMRTGA<T>::
resetAlphaAndBeta(const Real& a_alpha,
                  const Real& a_beta)
{
  Vector<MGLevelOp<T>* > ops = m_solver->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
    {

      TGAHelmOp<T>* helmop = (TGAHelmOp<T>*) ops[iop];
      helmop->setAlphaAndBeta(a_alpha, a_beta);
    }
}

/*****/
template <class T>
void AMRTGA<T>::
setTime(Real a_time)
{
  Vector<MGLevelOp<T>* > ops = m_solver->getAllOperators();
  for (int iop = 0; iop < ops.size(); iop++)
    {

      TGAHelmOp<T>* helmop = (TGAHelmOp<T>*) ops[iop];
      helmop->setTime(a_time);
    }
}

/*****/
template <class T>
AMRTGA<T>::
AMRTGA(const RefCountedPtr<AMRMultiGrid<T> >&                a_solver,
       const AMRLevelOpFactory<T> &                          a_opFact,
       const ProblemDomain&                                  a_level0Domain,
       const Vector<int>&                                    a_refRat,
       int a_numLevels,
       int a_verbosity)
{
  //cast to remove const because base class definition is weird
  m_verbosity = a_verbosity;
  m_level0Domain = a_level0Domain;
  m_refRat = a_refRat;
  m_solver  = a_solver;
  m_numLevels = a_numLevels;
  if (m_numLevels < 0)
    {
      m_numLevels = a_refRat.size();
    }

  m_ops.resize(m_numLevels);
  Vector< AMRLevelOp<T> * >& amrops =  m_solver->getAMROperators();
  for (int ilev = 0; ilev < m_numLevels; ilev++)
    {
      m_ops[ilev] = dynamic_cast<TGAHelmOp<T>* >(amrops[ilev]);
      if (m_ops[ilev]==NULL)
        {
          MayDay::Error("dynamic cast failed---is that operator really a TGAHelmOp?");
        }
    }

  setMu();
  m_dataCreated = false;
}
template <class T>
void AMRTGA<T>::
setMu()
{
  Real tgaEpsilon = 1.e-12;
#ifdef CH_USE_FLOAT
  tgaEpsilon = sqrt(tgaEpsilon);
#endif
  Real a = 2.0 - sqrt(2.0) - tgaEpsilon;
  m_mu1 = (a - sqrt( a*a - 4.0*a + 2.0))/2.0 ;
  m_mu2 = (a + sqrt( a*a - 4.0*a + 2.0))/2.0 ;
  m_mu3 = (1.0 - a);
  m_mu4 = 0.5 - a;
  if (m_verbosity > 4)
    {
      pout() << "   AMRTGA:: epsilon = " << tgaEpsilon << std::endl;
    }
}

template <class T>
void AMRTGA<T>::
createData(Vector<T* >&       a_source,
           int                a_lbase,
           int                a_lmax)
{
  m_dataCreated = true;
  m_rhst.resize(a_source.size(), NULL);
  m_srct.resize(a_source.size(), NULL);
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_rhst[ilev] = new T();
      m_srct[ilev] = new T();
      m_ops[ilev]->create(*m_rhst[ilev], *a_source[ilev]);
      m_ops[ilev]->create(*m_srct[ilev], *a_source[ilev]);
    }
}
/*****/
template <class T>
AMRTGA<T>::~AMRTGA()
{
  for (int ilev = 0; ilev < m_rhst.size(); ilev++)
    {
      if (m_rhst[ilev] != NULL)
        {
          delete m_rhst[ilev];
          delete m_srct[ilev];
          m_rhst[ilev] = NULL;
          m_srct[ilev] = NULL;
        }
    }
}
template <class T>
void AMRTGA<T>::
oneStep(Vector<T* >&       a_phiNew,
        Vector<T* >&       a_phiOld,
        Vector<T* >&       a_source,
        const Real&        a_dt,
        int                a_lbase,
        int                a_lmax,
        Real               a_told)
{
  if (!m_dataCreated)
    {
      createData(a_source, a_lbase, a_lmax);
    }

  if (m_verbosity > 3)
    {
      pout() << "  AMRTGA:: starting mu4 operation" << std::endl;
    }

  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->setToZero(*m_srct[ilev]);
      m_ops[ilev]->incr(*m_srct[ilev], *a_source[ilev], 1.0);

      //this sets srct to S/a
      m_ops[ilev]->divideByIdentityCoef(*m_srct[ilev]);
    }

  //this operation is homogeneous and therefore does not need time set
  //this makes rhs hold (k*a I + mu4 L) S/a
  applyHelm(m_rhst, m_srct, a_lbase, a_lmax, m_mu4, a_dt, true);
  //from here on k is kappa and L is kappa L
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      //this makes rhs hold       dt(k*a I + mu4 L) S/a
      m_ops[ilev]->scale(        *m_rhst[ilev], a_dt);
    }

  if (m_verbosity > 3)
    {
      pout() << "  AMRTGA:: starting mu3 operation" << std::endl;
    }

  //set time to t_old
  setTime(a_told);
  //this makes a_phiNew hold (k*a I + mu3 L) phi^n
  applyHelm(a_phiNew, a_phiOld, a_lbase, a_lmax, m_mu3, a_dt, false);

  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      //this makes rhs hold [(k*a I + mu3 L) phi^n + dt(k*a I + mu4 L) S/a]
      m_ops[ilev]->incr(*m_rhst[ilev], *a_phiNew[ilev], 1.0);
    }

  if (m_verbosity > 2)
    {
      pout() << "  AMRTGA:: starting mu2 operation" << std::endl;
    }
  //set time to t_intermediate
  setTime(a_told + a_dt*(m_mu2 + m_mu3));
  //this makes phinew = phiold if using phiold as init guess
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->assign(*a_phiNew[ilev], *a_phiOld[ilev]);
    }
  solveHelm(a_phiNew, m_rhst, a_lbase, a_lmax, m_mu2, a_dt);
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      //this puts the answer into rhst so we can do the final solve
      m_ops[ilev]->assign(*m_rhst[ilev], *a_phiNew[ilev]);
    }

  //this makes rhst hold k*a[(k*a I - mu2 L)^-1 (rhs)]
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->diagonalScale(*m_rhst[ilev]);
    }

  if (m_verbosity > 2)
    {
      pout() << "  AMRTGA:: starting mu1 operation" << std::endl;
    }
  //this happens at tnew
  setTime(a_told + a_dt);
  //this makes phinew = phiold if using phiold as init guess
  for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
    {
      m_ops[ilev]->assign(*a_phiNew[ilev], *a_phiOld[ilev]);
    }
  solveHelm(a_phiNew, m_rhst, a_lbase, a_lmax, m_mu1, a_dt);
}
/*******/
template <class T>
void AMRTGA<T>::
applyHelm(Vector<T*>&      a_ans,
          Vector<T*>&      a_phi,
          int              a_lbase,
          int              a_lmax,
          Real             a_mu,
          Real             a_dt,
          bool             a_homogeneousBC)
{
  Real factor  = a_mu*a_dt;

  resetAlphaAndBeta(1.0, factor);

  m_solver->computeAMROperator(a_ans,
                               a_phi,
                               a_lmax,
                               a_lbase,
                               a_homogeneousBC);

//   for (int ilev = a_lbase; ilev <= a_lmax; ilev++)
//     {
//       m_ops[ilev]->scale(*a_ans[ilev], -1);
//     }
}
/*******/
template <class T>
void AMRTGA<T>::
solveHelm(Vector<T*>&       a_ans,
          Vector<T*>&       a_rhs,
          int               a_lbase,
          int               a_lmax,
          Real              a_mu,
          Real              a_dt)
{
  Real factor  = -a_mu*a_dt;

  resetAlphaAndBeta(1.0, factor);

  m_solver->solveNoInit(a_ans,
                        a_rhs,
                        a_lmax,
                        a_lbase,
                        false);//zeroPhi
  if (m_solver->m_exitStatus==1 && m_verbosity>3)
    {
      pout() << "AMRTGA:: WARNING: solver exitStatus == 1" << std::endl;
    }
}

#include "NamespaceFooter.H"
#endif
