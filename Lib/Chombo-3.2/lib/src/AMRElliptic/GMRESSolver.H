#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// MFA, March 2011

#ifndef _GMRESSOLVER_H_
#define _GMRESSOLVER_H_

#include "LinearSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

///
/**
   Krylov solver using the GMRES algorithm.
      Important parmeters are:
        1) relative tolerance: m_reps(1.e-12)
        2) restart length: m_restrtLen(30)
 */
template <class T>
class GMRESSolver : public LinearSolver<T>
{
public:

  GMRESSolver();

  virtual ~GMRESSolver();

  void clearData();

  virtual void setHomogeneous(bool a_homogeneous)
  {
    m_homogeneous = a_homogeneous;
  }

  ///
  virtual void setConvergenceMetrics(Real a_metric,
                                     Real a_tolerance);

  ///
  /**
     define the solver.   a_op is the linear operator.
     a_homogeneous is whether the solver uses homogeneous boundary
     conditions.
   */
  virtual void define(LinearOp<T>* a_op, bool a_homogeneous=true);

  ///solve the equation.
  virtual void solve(T& a_phi, const T& a_rhs);

  ///
  /**
     public member data: whether the solver is restricted to
     homogeneous boundary conditions
   */
  bool m_homogeneous;

  ///
  /**
     public member data: operator to solve.
   */
  LinearOp<T>* m_op;

  ///
  /**
     private member data: restart length, need to allocate stuff
   */
private:
  int m_restrtLen;
public:
  int restartLen()const
  {
    return m_restrtLen;
  }
  void setRestartLen(int mm);
 ///
  /**
     public member data:  max iterations (eg, >= m_restrtLen)
   */
  int m_imax;
  ///
  /**
     public member data:  how much screen out put the user wants.
     set = 0 for no output.
   */
  int m_verbosity;

  ///
  /**
     public member data:  solver tolerance
   */
  Real m_eps;

  ///
  /**
     public member data:  relative solver tolerance
   */
  Real m_reps;

  ///
  /**
     public member data:
     set = -1 if solver exited for an unknown reason
     set =  1 if solver converged to tolerance
     set =  2 if rho = 0
     set =  3 if max number of restarts was reached
   */
  int m_exitStatus;

  ///
  /**
     public member data:  what the algorithm should consider "close to zero"
   */
  Real m_small;

  int m_normType;

private:
  void allocate();
  void CycleGMRES( T &xx,
                   int &reason, int &itcount, Real &rnorm0,
                   const bool avoidnorms = false);

  void ResidualGMRES( T &a_vv, const T &a_xx,
                      const T &a_bb, T &a_temp );

  void BuildGMRESSoln( Real nrs[], T &a_xx, const int it,
                       const T vv_0[] );

  void UpdateGMRESHessenberg( const int it, bool hapend, Real &res );

  //void ModifiedGramSchmidtOrthogonalization( const int it );
  //void UnmodifiedGramSchmidtOrthogonalization( const int it );
  void TwoUnmodifiedGramSchmidtOrthogonalization( const int it );

  void ApplyAB( T &a_dest, const T &a_xx, T &a_temp ) const;

  // data
  Real *m_data;
  Real *m_hes, *m_hh, *m_d, *m_ee, *m_dd;
  T    *m_work_arr;
};

// *******************************************************
// GMRESSolver Implementation
// *******************************************************

template <class T>
void GMRESSolver<T>::allocate()
{
  int max_k         = m_restrtLen;
  int hh            = (max_k + 2) * (max_k + 1);
  int hes           = (max_k + 1) * (max_k + 1);
  int rs            = (max_k + 2);
  int cc            = (max_k + 1);
  int size          = (hh + hes + rs + 2*cc + 1);

  m_data = new Real[size];
  m_hh = m_data;        // hh
  m_hes = m_data + hh;  // hes_
  m_d = m_hes + hes;    // rs_
  m_ee = m_d + rs;      // cc_
  m_dd = m_ee + cc;     // ss_

  m_work_arr = 0;
}

template <class T>
void GMRESSolver<T>::clearData()
{
  delete [] m_data;
}

template <class T>
void GMRESSolver<T>::setRestartLen(int mm)
{
  CH_assert(mm>0);
  clearData();
  m_restrtLen = mm;
  allocate();
}

template <class T>
GMRESSolver<T>::GMRESSolver()
  :m_homogeneous(false),
   m_op(NULL),
   m_restrtLen(30),
   m_imax(1000),
   m_verbosity(3),
   m_eps(1.0E-50),
   m_reps(1.0E-12),
   m_exitStatus(-1),
   m_small(1.0E-30),
   m_normType(2)
{
  allocate();
}

template <class T>
GMRESSolver<T>::~GMRESSolver()
{
  m_op = NULL;
  clearData();
}

template <class T>
void GMRESSolver<T>::define(LinearOp<T>* a_operator, bool a_homogeneous/*=true*/)
{
  m_homogeneous = a_homogeneous;
  m_op = a_operator;
}

////////////////////////////////////////////////////////////////////////
// GMRES
///////////////////////////////////////////////////////////////////////
#define HH(a,b)  (m_hh        + (b)*(m_restrtLen+2) + (a))
#define HES(a,b) (m_hes       + (b)*(m_restrtLen+1) + (a))
#define CC(a)    (m_ee        + (a))
#define SS(a)    (m_dd        + (a))
#define GRS(a)   (m_d         + (a))

/* vector names */
#define VEC_OFFSET 2
#define VEC_TEMP_RHS       m_work_arr[0]
#define VEC_TEMP_LHS       m_work_arr[1]
#define VEC_VV(i)          m_work_arr[VEC_OFFSET + i]

////////////////////////////////////////////////////////////////////////
// GMRESSolver<T>::solve
///////////////////////////////////////////////////////////////////////
template <class T>
void GMRESSolver<T>::solve( T& a_xx, const T& a_bb )
{
  CH_TIMERS("GMRESSolver::solve");

  CH_TIMER("GMRESSolver::solve::Initialize",timeInitialize);
  CH_TIMER("GMRESSolver::solve::MainLoop",timeMainLoop);
  CH_TIMER("GMRESSolver::solve::Cleanup",timeCleanup);

  CH_START(timeInitialize);

  if (m_verbosity >= 3)
    {
      pout() << "GMRESSolver::solve" << endl;
    }

  const int nwork = VEC_OFFSET + m_restrtLen + 1; // flex = VEC_OFFSET + 2*(m_restrtLen + 1);
  m_work_arr = new T[nwork];
  m_op->create(VEC_TEMP_RHS, a_bb);
  m_op->create(VEC_TEMP_LHS, a_xx);
  for (int i=VEC_OFFSET;i<nwork;i++)
    {
      m_op->create(m_work_arr[i], a_bb);
    }

  int it;
  Real rnorm0 = 0.0;
  T &vv_0 = VEC_VV(0);

  /* Compute the initial (preconditioned) residual (into 'vv_0')*/
  m_op->assign( vv_0, a_bb );
  m_op->setToZero( a_xx );

  CH_STOP(timeInitialize);

  ///
  /**         m_exitStatus
     set = -1 if solver exited for an unknown reason
     set =  1 if solver converged to tolerance
     set =  2 if rho = 0
     set =  3 if max number of restarts was reached
   */

  CH_START(timeMainLoop);

  // doit
  it = 0; m_exitStatus = -1;
  CycleGMRES( a_xx, m_exitStatus, it, rnorm0 );
  // loop for restarts
  while ( m_exitStatus==-1 && it < m_imax )
  {
    if (m_verbosity >= 3)
      {
        pout() << "*";
      }
    ResidualGMRES( vv_0, a_xx, a_bb, VEC_TEMP_RHS );
    CycleGMRES( a_xx, m_exitStatus, it, rnorm0 );
  }
  if (m_exitStatus==-1 && it >= m_imax)
    {
      m_exitStatus = 3;
    }

  CH_STOP(timeMainLoop);

  CH_START(timeCleanup);

  // clean up
  for (int i=0;i<nwork;i++)
    {
      m_op->clear(m_work_arr[i]);
    }
  delete [] m_work_arr; m_work_arr = 0;

  if (m_verbosity >= 3)
    {
      pout() << "GMRESSolver::solve done, status = " << m_exitStatus << endl;
    }

  CH_STOP(timeCleanup);
}

#define CONVERGED(r0,r) (r<r0*m_reps || r<m_eps)

template <class T>
void GMRESSolver<T>::CycleGMRES( T &a_xx,
                 int &a_reason, int &a_itcount, Real &a_rnorm0,
                                 const bool a_avoidnorms/*=false*/ )
{
  CH_assert(m_op != 0);
  Real        res_norm,res,hapbnd,tt;
  int         it;
  bool        hapend = false;
  T          &vv_0 = VEC_VV(0);

  /* scale VEC_VV (the initial residual) */
  //ierr = vv_0->Normalize( &res_norm ); CHKERRQ(ierr);
  res_norm = m_op->norm( vv_0, m_normType );
  res = res_norm;
  *GRS(0) = res_norm;

  /* check for the convergence */
  if ( res == 0. )
  {
    if (m_verbosity >= 3)
      {
        pout() << "GMRESSolver::solve zero residual!!!" << endl;
      }
    a_reason = 1; // 1 == converged,
    return;
  }
  // normilize
  tt = 1./res;
  m_op->scale( vv_0, tt);

  if ( a_itcount == 0 ) a_rnorm0 = res;
  bool conv = CONVERGED( a_rnorm0, res );
  a_reason = conv ? 1 : -1;
  it = 0;
  while ( a_reason == -1 && it < m_restrtLen && a_itcount < m_imax )
  {
    if ( m_verbosity>=4 && (it!=0 || a_itcount==0) )
      {
        pout() << a_itcount << ") GMRES residual = " << res << endl;
      }

    const T &vv_it = VEC_VV(it);
    T &vv_it1 = VEC_VV(it+1);

    // apply AB
    {
      T &Mb = VEC_TEMP_LHS;
      ApplyAB( vv_it1, vv_it, Mb );
    }
    /* update hessenberg matrix and do Gram-Schmidt */
    TwoUnmodifiedGramSchmidtOrthogonalization(it);

    /* vv(i+1) . vv(i+1) */
    //vv_it1->Normalize( &tt ); CHKERRQ(ierr);
    tt = m_op->norm( vv_it1, m_normType );
    /* check for the happy breakdown */
    //CH_assert( *GRS(it) != 0.0 );
    //hapbnd  = Abs( tt / *GRS(it) );
    hapbnd = 1.e-99; // hard wired happy tol!!!
    if (tt < hapbnd)
    {
      if ( m_verbosity>=3 )
        {
          pout() << "Detected happy breakdown, "<<it<<") current hapbnd = "<< tt << endl;
        }
      hapend = true;
    }
    else
      {
        m_op->scale( vv_it1, 1./tt);
      }

    /* save the magnitude */
    *HH(it+1,it) = tt; *HES(it+1,it) = tt;

    UpdateGMRESHessenberg( it, hapend, res );

    // increment
    it++; (a_itcount)++;
    conv = CONVERGED( a_rnorm0, res );
    a_reason = conv ? 1 : -1;
    /* Catch error in happy breakdown and signal convergence and break from loop */
    if ( hapend )
    {
      if ( !conv && m_verbosity>=3)
      {
        pout() << "You reached the happy break down, but convergence was not indicated. Residual norm=" << res << endl;
      }
      break;
    }
  }
  /* Monitor if we know that we will not return for a restart */
  if ( (a_reason!=1 || a_itcount>=m_imax) && m_verbosity>=4 )
    {
      pout() << a_itcount << ") GMRES residual = " << res << endl;
    }

  /*
    Down here we have to solve for the "best" coefficients of the Krylov
    columns, add the solution values together, and possibly unwind the
    preconditioning from the solution
  */
  /* Form the solution (or the solution so far) */
  BuildGMRESSoln( GRS(0), a_xx, it-1, &VEC_VV(0) );
}

template <class T>
void GMRESSolver<T>::ResidualGMRES( T &a_vv, const T &a_xx,
                                    const T &a_bb, T &a_temp_rhs )
{
  CH_assert(m_op != 0);
  m_op->applyOp( a_temp_rhs, a_xx, true); /* t <- Ax */
  m_op->assign( a_vv, a_bb );
  m_op->incr( a_vv, a_temp_rhs, -1.0); /* b - A(B)x */
}
template <class T>
void GMRESSolver<T>::BuildGMRESSoln( Real nrs[], T &a_xx, const int it,
                                     const T vv_0[] )
{
  Real tt;
  int ii,k,j;

  /* Solve for solution vector that minimizes the residual */

  /* If it is < 0, no gmres steps have been performed */
  if (it < 0)
  {
    return;
  }
  if (*HH(it,it) == 0.0)
    {
      pout() << "HH(it,it) is identically zero!!! GRS(it) = " << *GRS(it) << endl;
    }
  if (*HH(it,it) != 0.0)
  {
    nrs[it] = *GRS(it) / *HH(it,it);
  }
  else
  {
    nrs[it] = 0.0;
  }

  for (ii=1; ii<=it; ii++)
  {
    k   = it - ii;
    tt  = *GRS(k);
    for (j=k+1; j<=it; j++) tt  = tt - *HH(k,j) * nrs[j];
    nrs[k]   = tt / *HH(k,k);
  }

  /* Accumulate the correction to the solution of the preconditioned problem in TEMP */
  T &temp = VEC_TEMP_RHS;
  m_op->setToZero(temp);
  //temp->MAXPY( it+1, nrs, vv_0 );
  for (ii=0; ii<it+1; ii++)
    {
      m_op->incr(temp, vv_0[ii], nrs[ii]);
    }
  /* unwind pc */
  /*If we preconditioned on the right, we need to solve for the correction to
    the unpreconditioned problem */
  T &temp_matop = VEC_TEMP_LHS;
  //ierr = pc->Apply( temp, temp_matop );CHKERRQ(ierr);
  m_op->preCond( temp_matop, temp );
  m_op->incr( a_xx, temp_matop, 1.0 );
}
/* GMRESSolver::UpdateGMRESHessenberg *****************************************
 *
 *   INPUT:
 *     - it:
 *     - hapend: good breakdown?
 *     - res: residual (out)
 *
 *   SIDE EFFECTS:
 *     - sets 'nwork_' and allocs 'work_'
 *
 *   RETURN:
 *     - PETSc error code
 */
template <class T>
void GMRESSolver<T>::UpdateGMRESHessenberg( const int it, bool hapend, Real &res )
{
  Real *hh,*cc,*ss,tt;
  int   j;

  hh  = HH(0,it);
  cc  = CC(0);
  ss  = SS(0);

  /* Apply all the previously computed plane rotations to the new column
     of the Hessenberg matrix */
  for (j=1; j<=it; j++)
  {
    tt  = *hh;
    *hh = *cc * tt + *ss * *(hh+1);
    hh++;
    *hh = *cc++ * *hh - (*ss++ * tt);
  }

  /*
    compute the new plane rotation, and apply it to:
     1) the right-hand-side of the Hessenberg system
     2) the new column of the Hessenberg matrix
     thus obtaining the updated value of the residual
  */
  if ( !hapend )
  {
    tt = sqrt( *hh * *hh + *(hh+1) * *(hh+1) );
    if (tt == 0.0)
      {
        pout() << "Your matrix or preconditioner is the null operator\n";
      }
    *cc       = *hh / tt;
    *ss       = *(hh+1) / tt;
    *GRS(it+1) = - (*ss * *GRS(it));
    *GRS(it)   = *cc * *GRS(it);
    *hh       = *cc * *hh + *ss * *(hh+1);
    res      = Abs( *GRS(it+1) );
  }
  else
  {
    /* happy breakdown: HH(it+1, it) = 0, therfore we don't need to apply
       another rotation matrix (so RH doesn't change).  The new residual is
       always the new sine term times the residual from last time (GRS(it)),
       but now the new sine rotation would be zero...so the residual should
       be zero...so we will multiply "zero" by the last residual.  This might
       not be exactly what we want to do here -could just return "zero". */

    res = 0.0;
  }
}
/*
    This is the basic orthogonalization routine using modified Gram-Schmidt.
*/
// template <class T>
// void GMRESSolver<T>::ModifiedGramSchmidtOrthogonalization( const int it )
// {

//   int       ierr,j;
//   double    *hh,*hes,tmp;
//   PromVector_base *vv_1 = VEC_VV(it+1), *vv_j;
//   //const PromMatrix_base * const Amat = alloced_ ? A_ : alias_->A_;

//   PetscFunctionBegin;
//   /* update Hessenberg matrix and do Gram-Schmidt */
//   hh  = HH(0,it);
//   hes = HES(0,it);
//   for (j=0; j<=it; j++) {
//     /* (vv(it+1), vv(j)) */
//     vv_j = VEC_VV(j); ierr = vv_1->Dot( vv_j, hh );CHKERRQ(ierr);
//     *hes++ = *hh;
//     /* vv(it+1) <- vv(it+1) - hh[it+1][j] vv(j) */
//     tmp    = - (*hh++);
//     ierr   = vv_1->AXPY( tmp, vv_j );CHKERRQ(ierr);
//   }
//   PetscFunctionReturn(0);
// }

/*
  This version uses UNMODIFIED Gram-Schmidt.  It is NOT always recommended,
  but it can give MUCH better performance than the default modified form
  when running in a parallel environment.
 */

// int PromSolver::UnmodifiedGramSchmidtOrthogonalization( const int it )
// {
//   int       j,ierr;
//   double    *hh,*hes;
//   PromVector_base *vv_1 = VEC_VV(it+1);
//   const PromVector_base *const*vv_0 = &(VEC_VV(0));
//   //const PromMatrix_base * const Amat = alloced_ ? A_ : alias_->A_;

//   PetscFunctionBegin;
//   /* update Hessenberg matrix and do unmodified Gram-Schmidt */
//   hh  = HH(0,it);
//   hes = HES(0,it);
//   /*
//      This is really a matrix-vector product, with the matrix stored
//      as pointer to rows
//      */
//   ierr = vv_1->MDot( it+1, vv_0, hes );CHKERRQ(ierr);

//   /*
//     This is really a matrix-vector product:
//     [h[0],h[1],...]*[ v[0]; v[1]; ...] subtracted from v[it+1].
//   */

//   for (j=0; j<=it; j++) hh[j] = -hes[j];
//   ierr = vv_1->MAXPY( it+1, hh, vv_0 );CHKERRQ(ierr);
//   for (j=0; j<=it; j++) hh[j] = -hh[j];
//   PetscFunctionReturn(0);
// }

/*
  uses 1 iteration of iterative refinement of UNMODIFIED Gram-Schmidt.
  It can give better performance when running in a parallel
  environment and in some cases even in a sequential environment (because
  MAXPY has more data reuse).

  Care is taken to accumulate the updated HH/HES values.
 */
template <class T>
void GMRESSolver<T>::TwoUnmodifiedGramSchmidtOrthogonalization( const int it )
{
  Real     *hh,*hes,*lhh = 0;
  T        &vv_1 = VEC_VV(it+1);
  const T  *vv_0 = &(VEC_VV(0));

  /* Don't allocate small arrays */
  lhh = new Real[it+1];

  /* update Hessenberg matrix and do unmodified Gram-Schmidt */
  hh  = HH(0,it);
  hes = HES(0,it);

  /* Clear hh and hes since we will accumulate values into them */
  for (int j=0; j<=it; j++)
  {
    hh[j]  = 0.0;
    hes[j] = 0.0;
  }

  for ( int ncnt = 0 ; ncnt < 2 ; ncnt++ )
  {
    /*
       This is really a matrix-vector product, with the matrix stored
       as pointer to rows
    */
    //ierr = vv_1->MDot( it+1, vv_0, lhh ); CHKERRQ(ierr); /* <v,vnew> */
//     for (int j=0; j<=it; j++) // need a vector dot product!!!
//       {
//      lhh[j] = m_op->dotProduct(vv_1, vv_0[j]);
//       }
    m_op->mDotProduct(vv_1, it+1, vv_0, lhh);

    /*
      This is really a matrix vector product:
      [h[0],h[1],...]*[ v[0]; v[1]; ...] subtracted from v[it+1].
    */
    for (int j=0; j<=it; j++) lhh[j] = - lhh[j];
    //ierr = vv_1->MAXPY( it+1, lhh, vv_0 );           CHKERRQ(ierr);
    for (int j=0; j<=it; j++)
      {
        m_op->incr(vv_1, vv_0[j], lhh[j]);
      }
    for (int j=0; j<=it; j++)
    {
      hh[j]  -= lhh[j];     /* hh += <v,vnew> */
      hes[j] += lhh[j];     /* hes += - <v,vnew> */
    }
  }

  delete [] lhh;
}

/* PromSolver::PromPCApplyBAorAB ******************************************
 *
 */
template <class T>
void GMRESSolver<T>::ApplyAB( T &a_dest, const T &a_xx, T &a_tmp_lhs ) const
{
  m_op->preCond( a_tmp_lhs, a_xx );
  m_op->applyOp( a_dest, a_tmp_lhs, true);
}

template <class T>
void GMRESSolver<T>::setConvergenceMetrics(Real a_metric,
                                              Real a_tolerance)
{
  m_eps = a_tolerance;
}

#include "NamespaceFooter.H"
#endif /*_GMRESSOLVER_H_*/
