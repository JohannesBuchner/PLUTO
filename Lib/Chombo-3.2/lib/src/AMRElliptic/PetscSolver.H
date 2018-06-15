#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCSOLVER_H_
#define _PETSCSOLVER_H_

#include "LinearSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "LevelData.H"

#ifdef CH_USE_PETSC
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#endif

#include "NamespaceHeader.H"

///
/**
   Interface to PETSC linear solvers
   BC implementation hooks: addBCdiagValue (formMatrix) and addBCrhsValue (solveprivate)
   m_bccode should be encoded true for cells that are at the domain boundaries at formMatrix time.
 */
///
template <class T>
class PetscSolver : public LinearSolver<T>
{
public:
  PetscSolver();
  virtual ~PetscSolver();
  void destroy();
  ///
  /**
     Set whether the solver uses only homogeneous boundary conditions
   */
  virtual void setHomogeneous(bool a_homogeneous)
  {
    m_homogeneous = a_homogeneous;
  }
  /**
     Set Function F(u) = 0 and Jacobian dF(u)/du for nonlinear solves
  */
#ifdef CH_USE_PETSC
  virtual void setFunctionAndJacobian( PetscErrorCode (*f)(SNES,Vec,Vec,void*),
#if PETSC_VERSION_GE(3,5,0)
				       PetscErrorCode (*j)(SNES,Vec,Mat,Mat,void*)
#else
                                       PetscErrorCode (*j)(SNES,Vec,Mat*,Mat*,MatStructure*,void*)
#endif
                                       )
  {
    m_function = f;
    m_jacobian = j;
  }
  static PetscErrorCode ksp_monitor_pout(KSP ksp, PetscInt it, PetscReal rnorm  ,void *ctx)
  {
    pout() << "      KSP:: iteration = " << it << " residual norm = " <<  rnorm << std::endl;
    return 0;
  }
#endif
  //
  virtual void define( Real a_dx, bool a_homogeneous = true );     // convenience define
  virtual void define( LinearOp<T> *, bool a_homogeneous = true ); // pure virtual in base class
  ///
  /**
     Solve
  */
  virtual void solve( T& a_phi, const T& a_rhs );
  int solve_private(T& a_phi, const T& a_rhs);
  ///
  /**
     Solve a_op*a_phi = a_rhs using the PETSC matrix free functions
     The preconditioner (for which a matrix is formed) need not be 
     the same as the actual operator.
  */
  
  virtual void solve_mfree( T& a_phi, 
                            const T& a_rhs, 
                            LinearOp<T> *a_op );

  int solve_mfree_private(T& a_phi, 
                          const T& a_rhs, 
                          LinearOp<T> *a_op );
  
  LinearOp<T> *m_op_mfree;
  T m_phi_mfree;
  T m_Lphi_mfree;
  bool m_mfree_homogeneous;
  
#ifdef CH_USE_PETSC
  static PetscErrorCode apply_mfree(Mat A, Vec x, Vec f);
#endif
  ///

  /**
     viewResidual
  */
  Real computeResidual( );
  /**
     apply
  */
  virtual int applyOp(T & a_phi, const T & a_rhs );
  /**
     set initial guess non-zero
 */
  void setInitialGuessNonzero( bool b = true )
  {
    m_nz_init_guess = b;
  }

 /**
    set initial guess non-zero
 */
  void setOptionsPrefix( const char prefix[] )
  {
#ifdef CH_USE_PETSC
    if (m_ksp)
      {
        KSPSetOptionsPrefix(m_ksp,prefix);
      }
    else
      {
        strcpy(m_prestring, prefix); // should check len
      }
#endif
  }
#ifdef CH_USE_PETSC
  KSP getKSP() { return m_ksp; }
#endif

  /**
     get an estimate of the number of nnz/row for matrix allocation
  */
  virtual int getNNZPerRow() const
  {
    return 9;
  }
  /**
     do I support having formMatrix precompute exact NNZ/row
  */
  virtual bool supportNNZExact() const
  {
    return false;
  }

  /*
      an interface to apply operations to the rhs
  */
  virtual void rhsOperation(const T& a_rhs)
  {
  }
  ///
  /**
     public member data: whether or not to use inhomogeneous boundary conditions.
   */
  bool m_homogeneous;
  ///
  /**
     member data: grid spacing
   */
public:
  Real m_dx;

  ///
  /**
     handling boundary conditions, turn it into a term to be added to the diag term
     this function coordinates with addBCrhsValue for Dirichlet BC
     for Neumann BC no RHS modification is required
   */
  virtual Real addBCdiagValue(const IntVect& a_iv, const IntVect& a_jv, const T& a_rhs, DataIterator dit, const Real coeff = 1)
  {
    return 0;
  }
public:
  int resetOperator()
  {
#ifdef CH_USE_PETSC
    if (m_defined == 2)
    {
      m_defined = 1;
    }
#endif
    return 0;
  }

  ///
  /**
     Infinity norm
  */
  Real normInfinity( const T& a_phi ) const;

protected:
  ///
  /**
     handling boundary conditions, turn it into a term to be added to the RHS
     this function coordinates with addBCdiagValue for Dirichlet BC, should return a term that is to be added to RHS
   */
  virtual Real addBCrhsValue(const IntVect& a_iv, const T& a_phi, DataIterator dit, const Real& coeff = 1)
  {
    return 0.0;
  }

  /* petsc data */
#ifdef CH_USE_PETSC
public:
  Mat m_mat;
  void *m_ctx; // pointer for nonlnear solver call backs

  Vec m_xx, m_rr, m_bb;
  SNES m_snes;
  KSP m_ksp;
protected:
  PetscInt m_defined;
  PetscErrorCode (*m_function)(SNES,Vec,Vec,void*);
#if PETSC_VERSION_GE(3,5,0)
  PetscErrorCode (*m_jacobian)(SNES,Vec,Mat,Mat,void*);
#else
  PetscErrorCode (*m_jacobian)(SNES,Vec,Mat*,Mat*,MatStructure*,void*);
#endif
#endif
protected:
  char m_prestring[32];
  bool m_null;
  bool m_nz_init_guess;
  LevelData<BaseFab<bool> > m_bccode;
public:
  int setup_solver( const T& a_phi );
  int create_mat_vec( const T& a_phi );
#ifdef CH_USE_PETSC
  LevelData<BaseFab<PetscInt> >  m_gids;
  PetscInt m_gid0;
  PetscErrorCode putPetscInChombo( T& a_phi, const Vec xx );
  PetscErrorCode putChomboInPetsc( Vec out, const T& a_phi );
  virtual int formMatrix( Mat, const T* =0, PetscInt my0=0, PetscInt nloc=0, PetscInt *d=0, PetscInt *o=0) = 0;
#endif
  virtual BaseFab<Real>& getRegFab(T& a_fab, DataIterator& dit) = 0;
  const virtual BaseFab<Real>& getRegFab(const T& a_fab, DataIterator& dit) const = 0;
  const virtual BaseFab<Real>& getRegFab(const T& a_fab, DataIterator& dit, Box& a_box) const = 0;
  void setNull( bool n = true );
};

////////////////////////////////////////////////////////////////////////
// an instantiation of PetscSolver with FArrayBox
////////////////////////////////////////////////////////////////////////
template <class T>
class PetscSolverFAB : public PetscSolver<T>
{
public:
  PetscSolverFAB() : PetscSolver<LevelData<FArrayBox> >()
  {
  }

  BaseFab<Real>& getRegFab(LevelData<FArrayBox>& a_fab, DataIterator& dit)
  {
    BaseFab<Real>& fabref = static_cast<BaseFab<Real>& >(a_fab[dit()]);
    return fabref;
  }

  const BaseFab<Real>& getRegFab(const LevelData<FArrayBox>& a_fab, DataIterator& dit) const
  {
    const BaseFab<Real>& fabref = static_cast<const BaseFab<Real>& >(a_fab[dit()]);
    return fabref;
  }

  const BaseFab<Real>& getRegFab(const LevelData<FArrayBox>& a_fab, DataIterator& a_dit, Box& a_box) const
  {
    const BaseFab<Real>& fabref = static_cast<const BaseFab<Real>& >(a_fab[a_dit()]);
    a_box = fabref.box();
    return fabref;
  }
};

////////////////////////////////////////////////////////////////////////
// an instantiation of PetscSolver with Poisson
//     solve: alpha u - div( beta grad u) = f
//         - alpha and beta are constants
////////////////////////////////////////////////////////////////////////
template <class T>
class PetscSolverPoisson : public PetscSolverFAB<T>
{
public:
  PetscSolverPoisson();
  virtual void define( Real a_dx, bool a_homogeneous = true );
  virtual void define( LinearOp<T> *a_op, bool a_homogeneous = true )
  {
    PetscSolver<T>::define(a_op,a_homogeneous);
  }
#ifdef CH_USE_PETSC
  virtual int formMatrix(Mat, const LevelData<FArrayBox>* =0,PetscInt my0=0, PetscInt nloc=0, PetscInt *d=0, PetscInt *o=0 );
#endif
  Real m_alpha;
  Real m_beta;
};

////////////////////////////////////////////////////////////////////////
// an instantiation of PetscSolver with viscous tensor Poisson
//     solve: alpha a(x) + beta div( eta(x)grad(u) + I * lambda(x) div(u) ) = f
//         - alpha and beta are constants, a, eta and lambda are fields
////////////////////////////////////////////////////////////////////////
template <class T>
class PetscSolverViscousTensor : public PetscSolverFAB<LevelData<FArrayBox> >
{
public:
  PetscSolverViscousTensor();
  virtual void define( LinearOp<T>* a_operator, bool a_homogeneous = true );
  /**
     get an estimate of the number of nnz/row for matrix allocation
  */
  virtual int getNNZPerRow() const
  {
    return 13;
  }
#ifdef CH_USE_PETSC
  virtual int formMatrix(Mat, const LevelData<FArrayBox>* =0, PetscInt my0=0, PetscInt nloc=0, PetscInt *d=0, PetscInt *o=0 );
#endif
protected:
  Real m_alpha, m_beta;
  Real m_dxCrse;
  LevelData<FArrayBox> *m_a;
  LevelData<FluxBox> *m_eta;
  LevelData<FluxBox> *m_lamb;
public:
  void define( Real alpha, Real beta, LevelData<FArrayBox> *a, LevelData<FluxBox> *eta,  LevelData<FluxBox> *lam)
  {
    m_alpha = alpha;
    m_beta = beta;
    m_a = a;
    m_eta = eta;
    m_lamb = lam;
  }
};

#include "NamespaceFooter.H"

#ifdef CH_USE_PETSC
#ifndef CH_EXPLICIT_TEMPLATES
#include "PetscSolverI.H"
#endif // CH_EXPLICIT_TEMPLATES
#endif

#endif /*_PETSCSOLVER_H_*/
