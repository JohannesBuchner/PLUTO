/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Collection of root-finder algorithms.
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   June 20, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER  60
/* ************************************************************ */
int Brent(double (*Func)(double, void *), void *param, double x1, 
          double x2, double abs_acc, double rel_acc, double *xroot)
/*!
 * Use Brent's method to find the root of \f$ f(x) = 0\f$ known
 * to lie between x1 and x2. 
 * Here the function must be in the form \c double \c Func(x,*p), where 
 * where \c x is the independent variable while \c *p is a (void) pointer 
 * to any data type containing the function parameters.
 * Absolute accuracies can be specified to check for convergence
 * 
 * \param [in] *func    a pointer to a function func(x, *par)
 * \param [in] *param   a void pointer containing the parameters
 * \param [in] x1       the leftmost interval point
 * \param [in] x2       the rightmost interval point
 * \param [in] abs_acc  the desired absolute accuracy.(> 0 if you wish to 
 *                      use it).
 * \param [in] rel_acc  the desired relative accuracy.(> 0 if you wish to 
 *                      use it). 
 * \param [out] xroot   the zero of the function.
 *
 * \return The return value is 0 on success, 1 if the root could
 *         not be found to the desired level of accurat and 2 if the
 *         the initial interval does not contain the root. 
 *************************************************************** */
{
  int iter; 
  double a=x1,b=x2,c=x2,d,e,min1,min2; 
  double fa,fb,fc,p,q,r,s,tol1,xm;

  fa = (*Func)(a, param);
  fb = (*Func)(b, param);
    
  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
    WARNING(
      if(g_stepNumber > 0) print ("! Brent: initial interval does not contain root.\n");
    )
    return 2;   /* Error: initial interval does not contain root */
  }
  fc = fb; 
  for (iter=1; iter<=MAX_ITER; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c  = a; 
      fc = fa; 
      e  = d = b-a;
    } 
    if (fabs(fc) < fabs(fb)) {
      a  = b; 
      b  = c; 
      c  = a; 
      fa = fb; 
      fb = fc; 
      fc = fa;
    } 
    if (abs_acc < 0){
      tol1 = 2.0*2.0e-16*fabs(b);
    }else{
      tol1 = 2.0*2.0e-16*fabs(b) + 0.5*abs_acc; /* EPS = 2.0e-16 machine 
                                                   epsilon for double.    */
    } 
    xm   = 0.5*(c - b); 
/*    if (fabs(xm) <= tol1 || fb == 0.0){ */
    if (fabs(xm) <= 0.5*abs_acc || fabs(xm) <= 0.5*rel_acc*fabs(b)  ||
        fb == 0.0) {
      *xroot = b;
      g_maxRootIter = MAX(g_maxRootIter,iter);
      return 0; 
    }
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb/fa;     /*Attempt inverse quadratic interpolation.*/
      if (a == c) {
       	p = 2.0*xm*s;
	       q = 1.0 - s; 
      } else {
	       q = fa/fc; 
	       r = fb/fc; 
       	p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0)); 
       	q = (q-1.0)*(r-1.0)*(s-1.0);
      } 
      if (p > 0.0) q = -q;	/*Check whether in bounds. */
      p    = fabs(p); 
      min1 = 3.0*xm*q-fabs(tol1*q); 
      min2 = fabs(e*q); 
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	       e = d;
       	d = p/q; 
      } else {
	       d = xm; 
       	e = d;
      } 
    } else {
      d = xm; 
      e = d;
    } 
    a  = b; /* Move last best guess to a.*/
    fa = fb; 
    if (fabs(d) > tol1) b += d;  /*Evaluate new trial root.*/
    else                b += DSIGN(xm)*fabs(tol1); /* SIGN(tol1,xm);  */
    fb=(*Func)(b, param);
  }
  WARNING(
      print ("! Brent: exceed maximum iterations.\n");
  )
  return 1;
}
   
/* ********************************************************************* */   
int Ridder(double (*Func)(double, void *), void *param, 
           double x1, double x2, double abs_acc, double rel_acc, 
           double *xroot)
/*!
 * Use Ridder's method to find the root of \f$ f(x) = 0\f$ known
 * to lie between x1 and x2. 
 * Here the function must be in the form \c double \c Func(x,*p), where 
 * where \c x is the independent variable while \c *p is a (void) pointer 
 * to any data type containing the function parameters.
 * Both relative and absolute accuracies can be specified and convergence
 * is achieved by the condition that is realized first.
 * To choose among the two accuracies, set the other one to a negative
 * number.
 * 
 * \param [in] *func    a pointer to a function func(x, *par)
 * \param [in] *param   a void pointer containing the parameters
 * \param [in] x1       the leftmost interval point
 * \param [in] x2       the rightmost interval point
 * \param [in] abs_acc  the desired absolute accuracy (> 0 if you wish to 
 *                      use it).
 * \param [in] rel_acc  the desired relative accuracy (> 0 if you wish to
 *                      use it).
 * \param [out] xroot   the zero of the function.
 *
 * \return The return value is 0 on success, 1 if the root could
 *         not be found to the desired level of accurat and 2 if the
 *         the initial interval does not contain the root. 
 ************************************************************************ */
{
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
  double delta;

  fl = (*Func)(x1, param);
  fh = (*Func)(x2, param);
  
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl  = x1;
    xh  = x2;
    ans = -1.11e30;
    
    for (j = 1; j <= MAX_ITER; j++) {
      xm = 0.5*(xl + xh);
      fm = (*Func)(xm, param);
      s  = sqrt(fm*fm - fl*fh); 
      if (s == 0.0) {
        *xroot = ans;
        g_maxRootIter = MAX(g_maxRootIter,j);
        return 0;
      }
      xnew = xm + (xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s); /* update */
      delta = fabs(xnew-ans);
      if (delta <= abs_acc || delta <= rel_acc*fabs(xnew)) {
        *xroot = ans;
        g_maxRootIter = MAX(g_maxRootIter,j);
        return 0;
      }
      ans  = xnew;
      fnew = (*Func)(ans, param); 
      if (fnew == 0.0) {
        *xroot = ans;
        g_maxRootIter = MAX(g_maxRootIter,j);
        return 0;
      }
      if (fm*fnew < 0.0) {
        xl  = xm; 
        fl  = fm;
        xh  = ans;
        fh  = fnew;
      } else if (fl*fnew < 0.0) {
        xh = ans;
        fh = fnew;
      } else if (fh*fnew < 0.0) {
        xl = ans;
        fl = fnew;
      } else {
        print ("! Ridder: never get here.\n");
        print ("! xl   = %12.6e, fl   = %12.6e\n",xl,fl);
        print ("! xh   = %12.6e, fh   = %12.6e\n",xh,fh);
        print ("! xm   = %12.6e, fm   = %12.6e\n",xm,fm);
        print ("! xnew = %12.6e, fnew = %12.6e\n",xnew,fnew);
        print ("! sqrt^2 = %12.6e\n",fm*fm - fl*fh);
        QUIT_PLUTO(1);
      }
      
      delta = fabs(xh-xl);
      if (delta <= abs_acc || delta <= rel_acc*fabs(xh)) {
        *xroot = ans;
        g_maxRootIter = MAX(g_maxRootIter,j);
        return 0;
      }
    }
    WARNING(
      print ("! Ridder: exceed maximum iterations.\n");
    )
    return 1;  /* Error: max number of iteration exceeded */
    
  } else {
  
    if (fl == 0.0) {
      *xroot = x1;
      g_maxRootIter = MAX(g_maxRootIter,j);
      return 0;
    }
    if (fh == 0.0) {
      *xroot = x2;
      g_maxRootIter = MAX(g_maxRootIter,j);
      return 0;
    }
    /*  
    WARNING(
      print ("! Ridder: initial interval does not contain root.\n");
      )*/
    return 2;   /* Error: initial interval does not contain root */
  }
}
#undef MAX_ITER

/* ********************************************************************* */
void FDJacobian(int n, double x[], double fv[], double **df, 
		void (*vecfunc)(int , double [] , double [])){
/*
 * Estimate forward difference of the Jacobian.
 *
 *********************************************************************** */

  int i, j;
  double h, temp, *f;
  f = ARRAY_1D(n, double);

  for(j=0;j<n;j++){
    temp=x[j];
    h = EPS_FD_JAC*fabs(temp);
    if (h == 0.0) h = EPS_FD_JAC;
    x[j] = temp + h;
    h = x[j] - temp;
    (*vecfunc)(n,x,f);
    x[j] = temp;
    for(i=0;i<n;i++) df[i][j] = (f[i]-fv[i])/h; 
  }
}

/* ********************************************************************* */
void LineSearch (int n, double xold[], double fold, double g[], 
		 double p[], double x[], double mf[], double *f, double stpmax, 
		 int *check, void (*vecfunc)(int , double [] , double [])){
/*
 * Line search algorithm to minimize functions.
 *
 *********************************************************************** */
  int i;
  double a, alam, alam2, alamin, b, disc, f2, rhs1;
  double rhs2, slope, sum, temp, test, tmplam;
 
  *check = 0;
  sum = 0.0;
  for(i=0;i<n;i++) sum += p[i]*p[i];
  sum = sqrt(sum);
  if (sum > stpmax){
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  }
  slope = 0.0;
  for(i=0;i<n;i++) slope += g[i]*p[i];
  
  if(slope > 0.0){
    print("! LineSearch: Roundoff problem in lnscrh\n");
    QUIT_PLUTO(1);
  }

  test = 0.0;
  for(i=0;i<n;i++){
    temp = fabs(p[i])/MAX(fabs(xold[i]), 1.0);
    if (temp > test) test = temp;
  }
  
  alamin = TOLX/test;
  alam = 1.0;
  for(;;){
    for(i=0;i<n;i++) x[i] = xold[i] + alam*p[i];
    (*vecfunc)(n, x, mf);
    sum = 0.0;
    for (i=0;i< n;i++){
      sum += mf[i] * mf[i];
    }
    *f = 0.5*sum;
    if (alam < alamin){
      for(i=0;i<n;i++) x[i] = xold[i];
      *check = 1;
      return;
    }else if(*f < fold + ALF*alam*slope) return;
    else {
      if (alam == 1.0){
	       tmplam = -slope/(2.0*(*f-fold-slope));
      }else{
       	rhs1 = *f - fold - alam*slope;
       	rhs2 = f2 - fold - alam2*slope;
	       a = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam - alam2);
	       b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2))/(alam - alam2);
	       if (a == 0.0) tmplam = -slope/(2.0*b);
	       else {
	         disc = b*b - 3.0*a*slope;
	         if (disc < 0.0) tmplam = 0.5*alam;
       	  else if (b <= 0.0) tmplam = (-b + sqrt(disc))/(3.0*a);
	         else tmplam = -slope/(b + sqrt(disc));
	       }
	       if (tmplam > 0.5*alam) tmplam = 0.5*alam;
      }
    }
    alam2 = alam;
    f2 = *f;
    alam = MAX(tmplam, 0.1*alam);
  }
}

/* ********************************************************************* */
void Broyden(double x[], int n, int *check, 
             void (*vecfunc)(int, double [], double []))
/* 
 * Broyden multidimensional root finder. 
 *
 *********************************************************************** */
{
  int i, its, j, k, restrt, sing, skip;
  double den, f, fold, stpmax, sum, temp, test; 
  static double *c, *d, *fvcold, *g, *p, **qt, **r, *s, *t, *w, *xold, *fguess;
  
  if (g == NULL){
           
    if(n > MAX_ROOT_EQNS)
      print1 ("! Broyden : Number of equations exceed the maximum limit\n");
    
    c       = ARRAY_1D(MAX_ROOT_EQNS, double);
    d       = ARRAY_1D(MAX_ROOT_EQNS, double);
    fvcold  = ARRAY_1D(MAX_ROOT_EQNS, double);
    g       = ARRAY_1D(MAX_ROOT_EQNS, double);
    p       = ARRAY_1D(MAX_ROOT_EQNS, double);
    qt      = ARRAY_2D(MAX_ROOT_EQNS,MAX_ROOT_EQNS, double);
    r       = ARRAY_2D(MAX_ROOT_EQNS,MAX_ROOT_EQNS,double);
    s       = ARRAY_1D(MAX_ROOT_EQNS, double);
    t       = ARRAY_1D(MAX_ROOT_EQNS, double);
    w       = ARRAY_1D(MAX_ROOT_EQNS, double);
    xold    = ARRAY_1D(MAX_ROOT_EQNS, double); 
    fguess  = ARRAY_1D(MAX_ROOT_EQNS, double); 
  }
  
  /* Use fmin EXPLICTLY routine from math_minimization.c to create fvec. */
  
  (*vecfunc)(n, x, fguess);
  sum = 0.0;
  for (i=0;i< n;i++){
    sum += fguess[i] * fguess[i];
  }
  f = 0.5*sum;

  test = 0.0;
  
  /* Is initial guess the root ? */
  for (i=0; i<n; i++)
    if (fabs(fguess[i]) > test) test = fabs(fguess[i]);
  if (test < 0.01*TOLF) { 
    *check = 0;
    return;
  }
  
  /* Calculate stpmax for line searches */
  sum = 0.0;
  for (i=0;i<n;i++) sum += x[i] * x[i];
  stpmax = STPMX*MAX(sqrt(sum), (double)n);
  restrt = 1; /* To ensure initial computation of Jacobian. */

  /* Start the Iteration Loop */
  for(its = 1; its <= MAXITS; its++){

    if (restrt){
      FDJacobian(n,x,fguess,r,vecfunc); /* Compute Initial Jacobian in r */
      QRDecompose (r, n, c, d, &sing); /* QR Decomposition of r. */
      if (sing){
	       print ("! Broyden: singular Jacobian... Hard Luck !\n");
	       QUIT_PLUTO(1);
      }

      /* Initialize transpose(Q) to Identity matrix. */
      for(i=0;i<n;i++){
	       for(j=0;j<n;j++) qt[i][j] = 0.0;
	       qt[i][i] = 1.0;
      }
      
      /* Compute transpose(Q) explicityly */
      for(k=0; k < n-1; k++) {
	       if (c[k]) {
	         for (j=0; j<n; j++) {
	           sum = 0.0;
	           for(i=k; i<n; i++) {
	             sum += r[i][k]*qt[i][j];
	           }
	           sum /= c[k];
	           for(i=k; i<n; i++) {
	             qt[i][j] -= sum*r[i][k];
	           }
	         }
	       }
      }
      
      /* Create R */
      for(i=0;i<n;i++) {
	       r[i][i] = d[i]; /* Diagonal are in 'd' : QRDecompose */
	       for (j=0;j<i;j++) r[i][j] = 0.0;
      }
      
    } else { /* End of restart and now doing the Broyden Update */
      
      for(i=0; i<n; i++) s[i] = x[i] - xold[i]; /* s = delta(x) */
      for(i=0; i<n; i++) { /* t = R.s */
	       sum = 0.0;
	       for(j=i; j<n; j++) sum += r[i][j] * s[j];
	       t[i] = sum;
      }
      skip = 1;
      for (i=0; i<n; i++){ /* w = delta(F) - B.s */
	       sum = 0.0;
	       for(j=0; j<n; j++) sum += qt[j][i] * t[j];
	       w[i] = fguess[i] - fvcold[i] - sum;
	       if (fabs(w[i]) >= EPS_FD_JAC*(fabs(fguess[i]) + fabs(fvcold[i]))) skip = 0;
	/* No update with noisy components of w */
	       else w[i] = 0.0;
      }
      
      if(!skip) {
	       for(i=0;i<n;i++) { /* t = transpose(Q).w */
	         sum = 0.0;
	         for(j=0;j<n;j++) sum += qt[i][j]*w[j];
	         t[i] = sum;
	       }
	
	       den = 0.0;
	       for(i=0; i<n;i++) den += s[i] * s[i];
	       for(i=0; i<n;i++) s[i] /= den; /* Store s/(s.s) in s */

	       QRUpdate(r,qt,n,t,s); /* Update R and transpose(Q) */
	
	       for(i=0; i<n; i++) {
	         if (r[i][i] == 0.0){
	           print ("! Broyden: R is singular .. Hard Luck ! \n");
	           QUIT_PLUTO(1);
	         }
	         d[i] = r[i][i]; /* Diagonal of R in d */
	       }
      } /* End if(!skip) */
    } /* End of Broyden Update */
    
    /* Compute gradient(f) \approx transpose(Q.R) .F */
    for(i=0;i<n;i++){
      sum = 0.0;
      for(j=0; j<n; j++) sum += qt[i][j] * fguess[j];
      g[i] = sum;
    }

    for(i=n-1;i>=0;i--){
      sum = 0.0;
      for(j=0; j<=i; j++) sum += r[j][i] * g[j];
      g[i] = sum;
    }

    /* Store x and F */
    for(i=0;i<n;i++){
      xold[i] = x[i];
      fvcold[i] = fguess[i];
    }
    fold = f; /* Store Min. Least Sqrs from fmin */
    
    for(i=0;i<n;i++){ /* Compute RHS = -tranpose(Q).F */
      sum = 0.0;
      for(j=0; j<n; j++) sum += qt[i][j] * fguess[j];
      p[i] = -sum;
    }
    
    /* Solve Linear Equation */
    RSolve(r, n , d, p);
    
    /* Do Line Search to get xnew, f and fvec[xnew] */
    /*LineSearch (n, xold, fold, g, p, x, &f, stpmax, check, FMin); */
    LineSearch (n, xold, fold, g, p, x, fguess, &f, stpmax, check, vecfunc);
    
    /* Test for convergenvce on function values. */
    test = 0.0;
    for(i=0;i<n;i++){
      if (fabs(fguess[i]) > test) test = fabs(fguess[i]);
    }
    if (test < TOLF) { /* Function values converged */
      *check = 0;
      return;
    }
    
    if(*check) { /* Line Search Failed */
      if (restrt) {
	       print("! Broyden: Already tried reinitializing \n");
	       return;
      } else {
	       test = 0.0;
	       den = MAX(f, 0.5*n);
	       for(i=0; i<n; i++){
	         temp = fabs(g[i]) * MAX(fabs(x[i]), 1.0)/den;
	         if (temp > test) test = temp;
	       }
	       if (test < TOLMIN) return;
	       else               restrt = 1;
      }
    } else { /* Sucess Step */
      restrt = 0;
      test = 0.0;
      for(i=0; i<n; i++){ /* Test convergence on x */
       	temp = (fabs(x[i] - xold[i])) / MAX(fabs(x[i]), 1.0);
	       if (temp > test) test = temp;
      }
      if (test < TOLX){
	       print("! Brodyen: All is Well !! \n");
	       return;
      }
    }
  }
  print("! Brodyen: MAXITS exceeded in Broyden\n");
  return;
}

