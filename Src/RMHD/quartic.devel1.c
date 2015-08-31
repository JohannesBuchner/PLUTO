/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Solve quartic and cubic equations-
  

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define DEBUG  NO

void PrintSolution (double *z);
int QuarticNewton(double b, double c, double d, double e,   double *z);
void QuarticPrintCoeffs(double b, double c, double d, double e);
double CheckSolution(double b, double c, double d, double e, double x);
double ResolventCubic(double b, double c, double d, double e, double x );
int QuarticSolveNew (double b, double c, double d, double e, double *z);

int QuadraticSolve(double a, double b, double c, double *x);

static int debug_print = DEBUG;

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     These two functions may be moved into tools.c
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

#define swap(x,y) f = x; x = y; y = f;
#define SWAP(x,y) {double _t; _t = x; x = y; y = _t;}

/* ********************************************************************* */
int QuarticSolve (double b, double c, double d, double e, double *z)
/*!
 * Solve a quartic equation in the form 
 * \f[
 *      z^4 + bz^3 + cz^2 + dz + e = 0
 * \f]
 *  For its purpose, it is assumed that \b ALL roots are double. 
 *  This makes things faster.
 *
 * \param [in] b   coefficient of the quartic
 * \param [in] c   coefficient of the quartic
 * \param [in] d   coefficient of the quartic
 * \param [in] e   coefficient of the quartic
 * \param [out] z  vector containing the 
 *                 (double) roots of the quartic
 *   
 * \b Reference:
 *
 *   http://www.1728.com/quartic2.htm 
 * 
 * \return Return 0 on success, 1 if cubic solver fails, 2 if NaN has
 *         been found and 3 if quartic is not satisfied.
 *
 *********************************************************************** */
{
  int  n, j, ifail;
  double b2, f, g, h;
  double a2, a1, a0, u[4];
  double p, q, r, s;
  static double three_256 = 3.0/256.0;
  static double one_64 = 1.0/64.0;
  
#if DEBUG == YES
QuarticSolveNew(b,c,d,e,z);
PrintSolution(z);
exit(1);
#endif

  b2 = b*b;

  f = c - b2*0.375;
  g = d + b2*b*0.125 - b*c*0.5;
  h = e - b2*b2*three_256 + 0.0625*b2*c - 0.25*b*d;
  
  a2 = 0.5*f;
  a1 = (f*f - 4.0*h)*0.0625;
  a0 = -g*g*one_64;

  ifail = CubicSolve(a2, a1, a0, u);

  if (ifail) return 1;

  if (u[1] < 1.e-14){

    p = sqrt(u[2]);
    s = 0.25*b;
    z[0] = z[2] = - p - s;
    z[1] = z[3] = + p - s;
    
  }else{
  
    p = sqrt(u[1]);
    q = sqrt(u[2]);
  
    r = -0.125*g/(p*q);
    s =  0.25*b;
     
    z[0] = - p - q + r - s;
    z[1] =   p - q - r - s;
    z[2] = - p + q - r - s;
    z[3] =   p + q + r - s;

  }  

/* ----------------------------------------
    Sort roots in ascending order.
    Since z[0] < z[3], we first order

     z[0] < z[1] < z[3]

    and then put z[2] into the sequence.
   ---------------------------------------- */

  if      (z[1] < z[0]) {swap(z[1],z[0]);}
  else if (z[1] > z[3]) {swap(z[1],z[3]);}

 
  if (z[2] < z[0]) {
    swap(z[2],z[0]);
    swap(z[2],z[1]);
  }else if (z[2] < z[1]){
    swap(z[2],z[1]);
  }else if (z[2] > z[3]){
    swap(z[2],z[3]);
  }

/*
  for (j = 1; j < 4; j++){
    f = z[j];
    n = j - 1;
    while(n >= 0 && z[n] > f) {
       z[n+1] = z[n];
       n--;
    }
    z[n+1] = f;
  }
*/

  /* ----------------------------------------------
       verify that cmax and cmin satisfy original 
       equation
     ---------------------------------------------- */  

  for (n = 0; n < 4; n++){
    s = e + z[n]*(d + z[n]*(c + z[n]*(b + z[n])));
    if (s != s) {
//      print ("! QuarticSolve: NaN found.\n");
      return 2;
    }
    if (fabs(s) > 1.e-6) {
//      print ("! QuarticSolve: solution does not satisfy f(z) = 0; f(z) = %12.6e\n",s);
      return 3;
    }
  }


  double znew[4];
  ifail = QuarticSolveNew(b,c,d,e,znew);
  if (ifail) {
    debug_print = 1;
    ifail = QuarticSolveNew(b,c,d,e,znew);  
    print ("QuarticSolveNew() has failed\n");
    QUIT_PLUTO(1);
  }
  for (n = 0; n < 4; n++){
    if (znew[n] != znew[n]) {
      print ("! QuarticSolve(): nan in root %d\n",n);
      debug_print = 1;
      ifail = QuarticSolveNew(b,c,d,e,znew);  
      QUIT_PLUTO(1); 
    } 
  }


  if (fabs(z[0]-znew[0]) > 1.e-3 || fabs(z[3]-znew[3]) > 1.e-3){
    debug_print = 1;
    QuarticSolveNew(b,c,d,e,znew);

    printf ("! Different solutions\n");
    print ("  Old: "); PrintSolution(z);
    print ("  Verify[0]: %12.6e\n",CheckSolution(b,c,d,e,z[0]));
    print ("  Verify[3]: %12.6e\n",CheckSolution(b,c,d,e,z[3]));
  
    print ("\n");
    print ("  New: "); PrintSolution(znew);
    print ("  Verify[0]: %12.6e\n",CheckSolution(b,c,d,e,znew[0]));
    print ("  Verify[3]: %12.6e\n",CheckSolution(b,c,d,e,znew[3]));
    print ("\n  Now improving with Newton\n");
    
    
    QuarticNewton(b,c,d,e,znew);
    QuarticNewton(b,c,d,e,znew+3);
    print ("  Improved:"); PrintSolution(znew);
    print ("  Verify[0]: %12.6e\n",CheckSolution(b,c,d,e,znew[0]));
    print ("  Verify[3]: %12.6e\n",CheckSolution(b,c,d,e,znew[3]));

    QuarticPrintCoeffs(b,c,d,e);
    QUIT_PLUTO(1); 
  }
  z[0] = znew[0];
  z[3] = znew[3];

  return(0);
/*  
  printf (" z: %f ; %f ; %f ; %f\n",z[0], z[1], z[2], z[3]);
  */
}
#undef swap
/* ********************************************************************* */
int CubicSolve (double b, double c, double d, double z[])
/*!
 *  Solve a cubic equation in the form 
 *  \f[
 *      z^3 + bz^2 + cz + d = 0
 *  \f]
 *  For its purpose, it is assumed that \b ALL roots are double. 
 *  This makes things faster.
 *
 * \param [in] b coefficient of the cubic
 * \param [in] c coefficient of the cubic
 * \param [in] d coefficient of the cubic
 * \param [out] z vector containing the roots of the cubic.
 *                Roots should be sorted in increasing order.
 *   
 * \b Reference:  http://www.1728.com/cubic2.htm 
 *
 * \return Return 0 on success.
 *
 *********************************************************************** */
{
  double b2, g2;
  double f, g, h;
  double i, i2, j, k, m, n, p;
  static double one_3 = 1.0/3.0, one_27=1.0/27.0;

  b2 = b*b;
 
/*  ----------------------------------------------
     the expression for f should be 
     f = c - b*b/3.0; however, to avoid negative
     round-off making h > 0.0 or g^2/4 - h < 0.0
     we let c --> c(1- 1.1e-16)
    ---------------------------------------------- */

  f  = c*(1.0 - 1.e-16) - b2*one_3;
  g  = b*(2.0*b2 - 9.0*c)*one_27 + d; 
  g2 = g*g;
  i2 = -f*f*f*one_27;
  h  = g2*0.25 - i2;

/* --------------------------------------------
     double roots are possible only when 
   
               h <= 0 
   -------------------------------------------- */

  if (h > 1.e-12){
//    printf ("Only one double root (%12.6e)!\n", h);
  }
  if (i2 < 0.0){
/*
    printf ("i2 < 0.0 %12.6e\n",i2);
    return(1);
*/
    i2 = 0.0;
  }

/* --------------------------------------
       i^2 must be >= g2*0.25
   -------------------------------------- */
  
  i = sqrt(i2);       /*  > 0   */
  j = pow(i, one_3);  /*  > 0   */
  k = -0.5*g/i;

/*  this is to prevent unpleseant situation 
    where both g and i are close to zero       */

  k = (k < -1.0 ? -1.0:k);
  k = (k >  1.0 ?  1.0:k);
  
  k = acos(k)*one_3;       /*  pi/3 < k < 0 */
 
  m = cos(k);              /*   > 0   */
  n = sqrt(3.0)*sin(k);    /*   > 0   */
  p = -b*one_3;

  z[0] = -j*(m + n) + p;
  z[1] = -j*(m - n) + p;
  z[2] =  2.0*j*m + p;

/* ------------------------------------------------------
    Since j, m, n > 0, it should follow that from
    
      z0 = -jm - jn + p
      z1 = -jm + jn + p
      z2 = 2jm + p
      
    z2 is the greatest of the roots, while z0 is the 
    smallest one.
   ------------------------------------------------------ */
      
  return(0);
}










/* ********************************************************************* */
int QuarticSolveNew (double b, double c, double d, double e, double *z)
/*!
 * Solve a quartic equation in the form 
 * \f[
 *      z^4 + bz^3 + cz^2 + dz + e = 0
 * \f]
 *
 * 
 *  For its purpose, it is assumed that \b ALL roots are double. 
 *  This makes things faster.
 *
 * \param [in] b   coefficient of the quartic
 * \param [in] c   coefficient of the quartic
 * \param [in] d   coefficient of the quartic
 * \param [in] e   coefficient of the quartic
 * \param [out] z  vector containing the 
 *                 (double) roots of the quartic
 *   
 * \b Reference:
 *
 *    https://en.wikipedia.org/wiki/Quartic_function
 * \return Return 0 on success.
 *
 *********************************************************************** */
{
  int    status, k;
  double scrh, p, q, del, del0, del1, Q, S, Y;
  double sqp, sqm, phi, sq_del0;
  double D, P;
  double den, bnorm, dnorm;
  double zmin, zmax;
  const double one_third = 1.0/3.0;

#if DEBUG == YES  
/* -- four distinct roots in 3/10, -9/10, 9/10, 5/10 -- */
/*
b = -4.0/5.0;
c = -33.0/50.0;
d = 81.0/125.;
e = -243./2.e3;
*/
/* -- Two simple roots (6/10, -75/100) and one double root (2/10) -- */
/*
b = -0.25;
c = -0.47;
d =  93.0/500.0;
e = -9.0/500.;
*/
/* -- two double roots, 9/10 and -7/10 -- */
/*
e =  3969/1.e4;
d =  63./250.0;
c = -61./50.0;
b = -2.0/5.0;
*/

/* -- One quadrupole root in 0.99 -- */
/*
b = - 99./25.0;
c =   29403./5000.;
d =  -970299./250.e3;
e =   96059601./1.e8;
*/

/* -- Four distinct roots, very close -- */

b = -3.97367591434748e+00;
c =  5.92126911218894e+00;
d = -3.92150958389634e+00;
e =  9.73916387216783e-01;

/*
QuarticNewton (b,c,d,e,z);
print ("\n\n\n Using Newton:\n");
printf ("min and max roots = %12.6e,  %12.6e\n",z[0],z[3]);
exit(1);
*/
#endif

  double b2 = b*b, c2=c*c, bd=b*d;

  del0 = c2 - 3.0*bd + 12.0*e;   /* scales as a^2 */
  del1 = c*(2.0*c2 - 9.0*bd) + 27.0*(b2*e + d*d) - 72.0*c*e; /* scales as a^3 */
  del  = (4.0*del0*del0*del0 - del1*del1)/27.0;  /* >= 0.0, scales as a^6 */   
  D    = 64.0*e - 16.0*c2 + 16.0*b2*c - 16.0*bd - 3.0*b2*b2; /* scales as a^4 */
  P    = 8.0*c - 3.0*b2;
  
  Y    = 0.25*b2 + 2.0*d/(b+1.e-40) - c;

  den = 1.0 + fabs(b) + fabs(c) + fabs(d);
  bnorm = fabs(b)/den;
  dnorm = fabs(d)/den;

//#if DEBUG == YES
  if (debug_print){
    print ("QuarticSolve(): del = %12.6e, del0 = %12.6e D = %12.6e, P = %12.6e, Y = %12.6e\n",
         del,del0, D,P, Y);
  }
/*
   4 distinct roots:     del > 0
   2 simple, 1 double:   del = 0, P < 0, D < 0, del0 != 0
   2 double:             del = 0, P < 0, D = 0
   1 triple:             del = 0, D != 0
   1 quadrupole:         del = 0, del0 = 0, D = 0

*/
//#endif

  
double zero = 1.e-12;
double zero2 = zero*zero;
double zero3 = zero2*zero;


  if (fabs(del) < zero3 && fabs(del0) < zero && fabs(D) < zero2){
if (debug_print) print ("Case A: QuarticSolve(): 1 quadrupole root\n");
    z[0] = z[1] = z[2] = z[3] = -0.25*b;    
  } else if (bnorm < 1.e-12 && dnorm < 1.e-12){  /* Bi-quadratic */
if (debug_print) print ("Case B: QuarticSolve(): Biquadratic\n");
    double a2 = 1.0, b2 = c, c2 = e;
    double x[2];
    QuadraticSolve (a2,b2,c2,x);

    x[0] = MAX(0.0, x[0]);
    x[1] = MAX(0.0, x[1]);

    if (x[0] < 0.0 || x[1] < 0.0){
      print ("QuarticSolve(): case B: cannot continue\n");
      print ("x = %12.6e, %12.6e\n",x[0],x[1]);
      QUIT_PLUTO(1);
    }
    z[0] = -sqrt(x[1]);
    z[1] = -sqrt(x[0]);
    z[2] =  sqrt(x[0]);
    z[3] =  sqrt(x[1]);

  } else if (fabs(Y) < 1.e-14){
if (debug_print) print ("Case C: QuarticSolve(): No need for resolvent cubic\n");
    double x[2];
    double a2, b2, c2, del2;

    a2 = 1.0;
    b2 = 0.5*b;
    del2 = d*d/(b*b) - e;
    del2 = MAX(0.0, del2); /* Prevent small roundoff */
    if (del2 < 0.0){
      print ("! QuarticSolve(): case C: cannot continue, del2 = %12.6e\n",del2);
      QuarticPrintCoeffs(b,c,d,e);
      QUIT_PLUTO(1);
    }
    c2 = d/b + sqrt(del2);
    status = QuadraticSolve(a2,b2,c2,x);
    z[0] = x[0];
    z[3] = x[1];

    c2 = d/b - sqrt(del2);
    status = QuadraticSolve(a2,b2,c2,x);    
    z[1] = x[0];
    z[2] = x[1];

  } else {  /* Non degenerate case. Implies del0 > 0 */

if (debug_print) {
  print ("- QuarticSolve(): case 4: 4 distinct roots, no degeneracy\n");
  print ("  del0 = %12.6e, del1 = %12.6e\n",del0, del1);
}

del0 = MAX(del0,0.0);
    sq_del0 = sqrt(del0);
    scrh    = 0.5*del1/(fabs(del0)*sq_del0); /* This is always less than 1
                                                (in fabs)  since del > 0   */
if (debug_print) print ("  scrh = %8.3e\n",scrh);
scrh = MAX(scrh,-1.0);
scrh = MIN(scrh, 1.0);
    phi = acos(scrh)*one_third;
    p   = c - 3.0*b*b/8.0;
    q   = 0.125*b*b*b - 0.5*b*c + d;  /* this is bY/2 */

    scrh = -2.0*(p - sq_del0*cos(phi))*one_third;

scrh = MAX(scrh,0.0);
    if (scrh < 0.0){
      print ("! Quartic(): canont compute S, scrh = %12.6e\n",scrh);
      return 1;
    }


    S    = 0.5*sqrt(scrh);  /* > 0 */
if (debug_print) print ("  S = %12.6e\n",S);
    if (fabs(S) < 1.e-12){
      print ("! Quartic(): S = %12.6e <= 0 \n",S);
      return 2;
    }
  
    scrh = -4.0*S*S - 2.0*p - q/S;
    scrh = MAX(0.0, scrh);
if (debug_print) print ("  sqp^2 = %12.6e\n",scrh);
    sqm  = sqrt(scrh);
    scrh = -4.0*S*S - 2.0*p + q/S;
    scrh = MAX(0.0, scrh);
if (debug_print) print ("  sqp^2 = %12.6e\n",scrh);
    sqp  = sqrt(scrh);

  /* -- Sort roots by increasing value.
        This ordering has been established on testing only. -- */
  
    z[0] = -0.25*b - S - 0.5*sqp;
    z[1] = -0.25*b - S + 0.5*sqp;
    z[2] = -0.25*b + S - 0.5*sqm;
    z[3] = -0.25*b + S + 0.5*sqm;


if (fabs(CheckSolution(b,c,d,e,z[0]))>1.e-8){
  print ("! QuarticSolve: not correct (0)\n");
  QUIT_PLUTO(1);
}

  }

/* --------------------------------------------------------
     Check and sort
   -------------------------------------------------------- */

  
/*
if (! ( (z[0] <= z[1]) && (z[1] <= z[2]) && (z[2] <= z[3]))){
  print ("! QuarticSolve: roots are not correctly ordered\n");
  print ("! S = %12.6e, sqp = %12.6e, sqm = %12.6e\n",S,sqp,sqm);
  PrintSolution(z);
  QUIT_PLUTO(1);
}
*/

/* Sort roots in ascending order */
/*
  if      (z[1] < z[0]) {SWAP(z[1],z[0]);}
  else if (z[1] > z[3]) {SWAP(z[1],z[3]);}

 
  if (z[2] < z[0]) {
    SWAP(z[2],z[0]);
    SWAP(z[2],z[1]);
  }else if (z[2] < z[1]){
    SWAP(z[2],z[1]);
  }else if (z[2] > z[3]){
    SWAP(z[2],z[3]);
  }
*/
  return 0;  
}

/* ********************************************************************* */
int QuarticNewton(double b, double c, double d, double e,
                  double *z)
/*!
 * Solve the quartic equation using Newton method for multiple roots.
 *
 * \b Reference
 *    - "A first course in Numerical Analysis"
 *      Ralston & Rabinowitz, Eq. [8.6-13] -- [8.6-22]
 *
 *********************************************************************** */
{
  int    k, max_iter=40;
  double x, dx, f, df, d2f;
  double tolx = 1.e-12, tolf = 1.e-18;


  x = *z;
  for (k = 0; k < max_iter; k++){
    f   = e + x*(d + x*(c + x*(b + x)));
    df  = d + x*(2.0*c + x*(3.0*b + 4.0*x));
    d2f = 2.0*c + x*(6.0*b + 12.0*x);
    
    dx = f*df/(df*df - f*d2f); /* search for the root of u=f/df rather than f.
                                  The increment than becomes
                                  u/du = f*fd/(df^2 - f*d2f), see the
                                  discussion before Eq. [8.6-22]. */
    x -= dx;

    if (fabs(dx) < tolx || fabs(f) < tolf) {
      *z = x;

if (debug_print){
  print ("QuarticNewton: k = %d, x = %f, dx = %12.6e,  f = %8.3e\n",k,x, dx,f);
  print ("QuarticNewton, root found in # %d iterations\n",k);
}
      return 0;
    }
  }
  print ("! QuarticNewton: too many steps\n");
  return 1;
  
}
#undef DEBUG

/* ********************************************************************* */
void QuarticPrintCoeffs(double b, double c, double d, double e)
/*
 *
 ********************************************************************** */
{
  print ("! f(x) = %18.12e + x*(%18.12e + x*(%18.12e ",
         e, d, c);
      print ("  + x*(%18.12e + x*%18.12e)))\n", b, 1.0);

  print ("b = %18.14e;\n",b);
  print ("c = %18.14e;\n",c);
  print ("d = %18.14e;\n",d);
  print ("e = %18.14e;\n",e);

  double b3 = -2.0*c;
  double c3 = c*c + b*d - 4.0*e;
  double d3 = -(b*c*d - b*b*e - d*d);

  print ("! Resolvent cubic:\n");
  print ("! g(x) = %18.12e + x*(%18.12e + x*(%18.12e + x))\n", d3, c3, b3);
  
}

/* ********************************************************************* */
double CheckSolution(double b, double c, double d, double e, double x)
/*
 *
 ********************************************************************** */
{
  double f, fp, fm, fmax;
  f =  e + x*(d + x*(c + x*(b + x)));

  x  = 1.0;
  fp = e + x*(d + x*(c + x*(b + x)));

  x  = -1.0;
  fm = e + x*(d + x*(c + x*(b + x)));
   
  fmax = MAX(fp,fm);
  return f/fmax;
}
 
/* ********************************************************************* */
void PrintSolution (double *z)
/*
 *
 ********************************************************************** */
{
  print ("z = %f  %f  %f  %f\n",z[0],z[1],z[2],z[3]);
}

double ResolventCubic(double b, double c, double d, double e, double x )
{
  double b3,c3,d3, f;

  b3 = -2.0*c;
  c3 = c*c + b*d - 4.0*e;
  d3 = -(b*c*d - b*b*e - d*d);
  f  = d3 + x*(c3 + x*(b3 + x));
  return f;
}

/* ********************************************************************* */
int QuadraticSolve(double a, double b, double c, double *x)
/*! 
 *  Solve a quadratic equation in the form
 *
 *   ax^2 + bx + c = 0
 *
 * Return roots in increasing order
 *
 *********************************************************************** */
{
  double del, sb, q;

  del = b*b - 4.0*a*c;
  if (del < 0.0) return 1; /* No real root */

  if (b == 0){
    x[1] = sqrt(c/a);
    x[0] = -x[1];
  }

  sb = (b > 0.0 ? 1.0:-1.0);
  q  = -0.5*(b + sb*sqrt(del));
  x[0] = q/a;
  x[1] = c/q;

  if (x[0] > x[1]) {
    double scrh = x[1];
    x[1] = x[0];
    x[0] = scrh;
  }

  return 0;
}


