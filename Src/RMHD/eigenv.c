/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the eigenvalues for the relativisitc MHD equations.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
int MaxSignalSpeed (double **v, double *a2, double *h,
                    double *cmin, double *cmax, int beg, int end)
/*!
 *  Return the rightmost (cmax) and leftmost (cmin) 
 *  wave propagation speed in the Riemann fan
 *
 *********************************************************************** */
{
  int    i, err;
  double lambda[NFLX];

  for (i = beg; i <= end; i++) {
    err = Magnetosonic (v[i], a2[i], h[i], lambda);
    if (err != 0) return i;
    cmax[i] = lambda[KFASTP];
    cmin[i] = lambda[KFASTM];
  }
  return 0;
}

/* ********************************************************************* */
int Magnetosonic (double *vp, double cs2, double h, double *lambda)
/*!
 * Find the four magneto-sonic speeds (fast and slow) for the 
 * relativistic MHD equations (RMHD).
 * We follow the notations introduced in Eq. (16) in
 *
 * Del Zanna, Bucciantini and Londrillo, A&A, 400, 397--413, 2003 
 * 
 * and also Eq. (57-58) of Mignone & Bodo, MNRAS, 2006 for the
 * degenerate cases.
 *
 * The quartic equation is solved analytically and the solver
 * (function 'QuarticSolve') assumes all roots are real 
 * (which should be always the case here, since eigenvalues 
 * are recovered from primitive variables). The coefficients
 * of the quartic were found through the following MAPLE
 * script:
 * \verbatim
 *   ------------------------------------------------
restart;

u[0] := gamma;
u[x] := gamma*v[x];
b[0] := gamma*vB;
b[x] := B[x]/gamma + b[0]*v[x];
wt   := w + b^2;

#############################################################
"fdZ will be equation (16) of Del Zanna 2003 times w_{tot}";

e2  := c[s]^2 + b^2*(1 - c[s]^2)/wt;
fdZ := (1-e2)*(u[0]*lambda - u[x])^4 + (1-lambda^2)*
       (c[s]^2*(b[0]*lambda - b[x])^2/wt - e2*(u[0]*lambda - u[x])^2);

########################################################
"print the coefficients of the quartic polynomial fdZ";

coeff(fMB,lambda,4);
coeff(fMB,lambda,3);
coeff(fMB,lambda,2);
coeff(fMB,lambda,1);
coeff(fMB,lambda,0);

fdZ := fdZ*wt;  

######################################################
"fMB will be equation (56) of Mignone & Bodo (2006)";

a  := gamma*(lambda-v[x]);
BB := b[x] - lambda*b[0];
fMB := w*(1-c[s]^2)*a^4 - (1-lambda^2)*((b^2+w*c[s]^2)*a^2-c[s]^2*BB^2);

######################################
"check that fdZ = fMB";

simplify(fdZ - fMB);

########################################################
"print the coefficients of the quartic polynomial fMB";

coeff(fMB,lambda,4);
coeff(fMB,lambda,3);
coeff(fMB,lambda,2);
coeff(fMB,lambda,1);
coeff(fMB,lambda,0);


###############################################
"Degenerate case 2: Bx = 0, Eq (58) in MB06";

B[x] := 0;

a2 := w*(c[s]^2 + gamma^2*(1-c[s]^2)) + Q; 
a1 := -2*w*gamma^2*v[x]*(1-c[s]^2);
a0 := w*(-c[s]^2 + gamma^2*v[x]^2*(1-c[s]^2))-Q;

"delc is a good way to evaluate the determinant so that it is > 0";
del  := a1^2 - 4*a2*a0;
delc := 4*(w*c[s]^2 + Q)*(w*(1-c[s]^2)*gamma^2*(1-v[x]^2) + (w*c[s]^2+Q));
simplify(del-delc);

############################################################################
"check the correctness of the coefficients of Eq. (58) in MB06";

Q  := b^2 - c[s]^2*vB^2;
divide(fMB, (lambda-v[x])^2,'fMB2');
simplify(a2 - coeff(fMB2,lambda,2)/gamma^2);
simplify(a1 - coeff(fMB2,lambda,1)/gamma^2);
simplify(a0 - coeff(fMB2,lambda,0)/gamma^2);
     
\endverbatim
 *
 *********************************************************************** */
{
  int   iflag;
  double  vB, vB2, Bm2, vm2, w, w_1;
  double  eps2, one_m_eps2, a2_w;
  double  vx, vx2, u0, u02;
  double  a4, a3, a2, a1, a0;
  double  scrh1, scrh2, scrh;
  double  b2, del, z[4];

#if RMHD_FAST_EIGENVALUES
  Eigenvalues (vp, cs2, h, lambda);
  return 0;
#endif

  scrh   = fabs(vp[VXn])/sqrt(cs2);
  g_maxMach = MAX(scrh, g_maxMach);

  vB  = EXPAND(vp[VX1]*vp[BX1], + vp[VX2]*vp[BX2], + vp[VX3]*vp[BX3]);
  u02 = EXPAND(vp[VX1]*vp[VX1], + vp[VX2]*vp[VX2], + vp[VX3]*vp[VX3]);
  Bm2 = EXPAND(vp[BX1]*vp[BX1], + vp[BX2]*vp[BX2], + vp[BX3]*vp[BX3]);
  vm2 = u02;

  if (u02 >= 1.0){
    print ("! MAGNETOSONIC: |v|= %f > 1\n",u02);
    return 1;
  }

  if (u02 < 1.e-14) {
    
  /* -----------------------------------------------------
      if total velocity is = 0, the eigenvalue equation 
      reduces to a biquadratic:
         
          x^4 + a1*x^2 + a0 = 0
            
      with 
           a0 = cs^2*bx*bx/W > 0,
           a1 = - a0 - eps^2 < 0.
     ----------------------------------------------------- */
       
    w_1  = 1.0/(vp[RHO]*h + Bm2);   
    eps2 = cs2 + Bm2*w_1*(1.0 - cs2);
    a0   = cs2*vp[BXn]*vp[BXn]*w_1;
    a1   = - a0 - eps2;
    del  = a1*a1 - 4.0*a0;
    del  = MAX(del, 0.0);
    del  = sqrt(del); 
    lambda[KFASTP] =  sqrt(0.5*(-a1 + del));
    lambda[KFASTM] = -lambda[KFASTP];
    #if COMPONENTS > 1
     lambda[KSLOWP] =  sqrt(0.5*(-a1 - del));
     lambda[KSLOWM] = -lambda[KSLOWP];
    #endif
    return 0;
  }

  vB2 = vB*vB;
  u02 = 1.0/(1.0 - u02);
  b2  = Bm2/u02 + vB2;
  u0  = sqrt(u02);
  w_1 = 1.0/(vp[RHO]*h + b2);   
  vx  = vp[VXn];
  vx2 = vx*vx;

  if (fabs(vp[BXn]) < 1.e-14){

    w     = vp[RHO]*h;
    scrh1 = b2 + cs2*(w - vB2); /* wc_s^2 + Q */
    scrh2 = w*(1.0 - cs2)*u02;
    del   = 4.0*scrh1*(scrh2*(1.0 - vx2) + scrh1);
      
    a2  = scrh1 + scrh2;
    a1  = -2.0*vx*scrh2;

    lambda[KFASTP] = 0.5*(-a1 + sqrt(del))/a2;
    lambda[KFASTM] = 0.5*(-a1 - sqrt(del))/a2;
    #if COMPONENTS > 1
     lambda[KSLOWP] = lambda[KSLOWM] = vp[VXn];
    #endif

    return 0;

  }else{

    scrh1 = vp[BXn]/u02 + vB*vx;  /* -- this is bx/u0 -- */
    scrh2 = scrh1*scrh1;  
                     
    a2_w       = cs2*w_1;
    eps2       = (cs2*vp[RHO]*h + b2)*w_1;
    one_m_eps2 = u02*vp[RHO]*h*(1.0 - cs2)*w_1;

    /* ---------------------------------------
         Define coefficients for the quartic  
       --------------------------------------- */
    
    scrh = 2.0*(a2_w*vB*scrh1 - eps2*vx);
    a4 = one_m_eps2 - a2_w*vB2 + eps2;
    a3 = - 4.0*vx*one_m_eps2 + scrh;
    a2 =   6.0*vx2*one_m_eps2 + a2_w*(vB2 - scrh2) + eps2*(vx2 - 1.0);
    a1 = - 4.0*vx*vx2*one_m_eps2 - scrh;
    a0 = vx2*vx2*one_m_eps2 + a2_w*scrh2 - eps2*vx2;

    if (a4 < 1.e-12){
      print ("! Magnetosonic: cannot divide by a4\n");
      return 1;
    }

    scrh = 1.0/a4;
     
    a3 *= scrh;
    a2 *= scrh;
    a1 *= scrh;
    a0 *= scrh;
    iflag = QuarticSolve(a3, a2, a1, a0, z);
  
    if (iflag){
      print ("! Magnetosonic: Cannot find max speed [dir = %d]\n", g_dir);
      print ("! rho = %12.6e\n",vp[RHO]);
      print ("! prs = %12.6e\n",vp[PRS]);
      EXPAND(print ("! vx  = %12.6e\n",vp[VX1]);  ,
             print ("! vy  = %12.6e\n",vp[VX2]);  ,
             print ("! vz  = %12.6e\n",vp[VX3]);)
      EXPAND(print ("! Bx  = %12.6e\n",vp[BX1]);  ,
             print ("! By  = %12.6e\n",vp[BX2]);  ,
             print ("! Bz  = %12.6e\n",vp[BX3]);)
      print ("! f(x) = %12.6e + x*(%12.6e + x*(%12.6e ",
               a0*a4, a1*a4, a2*a4);
      print ("  + x*(%12.6e + x*%12.6e)))\n", a3*a4, a4);
      return 1;
    }
/*
if (z[3] < z[2] || z[3] < z[1] || z[3] < z[0] ||
    z[2] < z[1] || z[2] < z[0] || z[1] < z[0]){
  printf ("!Sorting not correct\n");
  printf ("z = %f %f %f %f\n",z[0],z[1],z[2],z[3]);
  QUIT_PLUTO(1);
}
*/
  
    lambda[KFASTM] = z[0];
    lambda[KFASTP] = z[3];
    #if COMPONENTS > 1
     lambda[KSLOWP] = z[2];
     lambda[KSLOWM] = z[1];
    #endif
  }
  return 0;
}

/* ********************************************************************* */
int Eigenvalues(double *v, double cs2, double h, double *lambda)
/*!
 * Compute an approximate expression for the fast magnetosonic speed
 * using the upper-bound estimated outlined by
 * Leismann et al. (A&A 2005, 436, 503), Eq. [27].
 *
 *********************************************************************** */
{
  double vel2, lor2, Bmag2, b2, ca2, om2, vB2;
  double vB, scrh, vl, vx, vx2;

  vel2  = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
  vB    = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
  Bmag2 = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);  

  if (vel2 >= 1.0){
    print ("! Eigenavalues(): |v|^2 = %f > 1\n",vel2);
    return 1;
  }

  vx  = v[VXn];
  vx2 = vx*vx;
  vB2 = vB*vB;
  b2  = Bmag2*(1.0 - vel2) + vB2;
  ca2 = b2/(v[RHO]*h + b2);

  om2  = cs2 + ca2 - cs2*ca2;
  vl   = vx*(1.0 - om2)/(1.0 - vel2*om2);
  scrh = om2*(1.0 - vel2)*((1.0 - vel2*om2) - vx2*(1.0 - om2));
  scrh = sqrt(scrh)/(1.0 - vel2*om2);
  lambda[KFASTM] = vl - scrh;
  lambda[KFASTP] = vl + scrh;

  if (fabs(lambda[KFASTM])>1.0 || fabs(lambda[KFASTP]) > 1.0){
    print ("! Eigenvalues(): vm, vp = %8.3e, %8.3e\n",lambda[KFASTM],lambda[KFASTP]);
    QUIT_PLUTO(1);
  }
  if (lambda[KFASTM] != lambda[KFASTM] || lambda[KFASTP] != lambda[KFASTP]){
    print ("! Eigenvalues(): nan, vm, vp = %8.3e, %8.3e\n",lambda[KFASTM],lambda[KFASTP]);
    QUIT_PLUTO(1);
  }
   
  scrh      = fabs(vx)/sqrt(cs2);
  g_maxMach = MAX(scrh, g_maxMach);
}





#if 1==0   /* this part of the code is beta-testing 
              and is excluded from normal usage  */
#define DEBUG_MODE NO
/* **************************************************************** */
void PRIM_EIGENVECTORS(real *q, real cs2, real h, real *lambda,
                       real **LL, real **RR)
/*
 *
 * PURPOSE
 *
 *   Provide right and left eigenvectors for the primitive
 *   variables (\rho, v, p, B).
 *   We use the procedure outlined in 
 *
 *   "Total Variation Diminishing Scheme for Relativistic MHD"
 *   Balsara, Apj (2001), 132:83.
 *         
 *  The notations are identical, except for the covariant
 *  magnetic field for which we use b^\mu instead of h^\mu
 *  We use square brackets to describe four-vector, always 
 *  intended with an upper subscript, 
 *
 *      v[mu] = v^\mu
 *
 *  !Beware of the distinction between lower and upper index!
 *
 *   u^\mu   = \gamma*(1,v)     u_\mu   = \gamma(-1,v)
 *   phi^\mu = (\lambda,1,0,0)  phi_\mu = (-\lambda,1,0,0)
 *
 *  
 *  
##############################################
#
# Maple script to check the behavior of 
# the eigenvalue equations 
# (Useful to check degeneracies)
#
##############################################

restart;

v[x] := 0;  v[y] := 0; v[z] := 0;
B[x] := BB; B[y] := 0; B[z] := 0;

vB := v[x]*B[x] + v[y]*B[y] + v[z]*B[z];
v2 := v[x]*v[x] + v[y]*v[y] + v[z]*v[z];

u[0] := 1/sqrt(1-v2);
u[1] := u[0]*v[x];
u[2] := u[0]*v[y];
u[3] := u[0]*v[z];

b[0] := u[0]*vB;
b[1] := B[x]/u[0] + u[x]*vB;
b[2] := B[y]/u[0] + u[y]*vB;
b[3] := B[z]/u[0] + u[z]*vB;

b2  := (B[x]^2 + B[y]^2 + B[z]^2)/u[0]^2 + vB^2;

#lambda := b[1]/(u[0]*sqrt(w));
e2 := cs^2*(1-b2/w) + b2/w;
f  := (1-e2)*(u[0]*lambda - u[1])^4 +
     (1-lambda^2)*(cs^2/w*(b[0]*lambda - b[1])^2 -
                     e2*(u[0]*lambda - u[1])^2);
simplify(f);
#################################################
 *
 * LAST MODIFIED
 *
 *   11 Aug 2009 by A. Mignone (mignone@ph.unito.it)
 *
 ****************************************************************** */
{
  int    nv, mu, k, indx[NFLX];
  int    degen1 = 0, degen2[NFLX];
  double vB, Bmag2, eta, b2, E;
  double eta_1, e_P, rho_P, rho_S;
  double a, A, B, betay, betaz;
  double u[4], b[4], l[4], f[4], d[4];
  double la[4], lB[4], fa[4], fB[4];
  double phi[4];
  double s, Ru[10], col[NFLX];
  static double **tmp;

#if DEBUG_MODE == YES
q[RHO] = 1.7;
q[VXn] = 0.;
q[VXt] = 0.0;
q[VXb] = 0.0;
q[BXn] = 0.65723;
q[BXt] = 0.0;
q[BXb] = 0.0;
q[PRS] = 2.;
cs2 = g_gamma*q[PRS]/(q[PRS]+g_gamma/(g_gamma-1.0)*q[RHO]);
h  = 1.0+g_gamma/(g_gamma-1.0)*q[PRS]/q[RHO];

printf ("AA = %12.6e\n",q[BXn]/sqrt(q[BXn]*q[BXn] + q[RHO]*h));
#endif

  if (tmp == NULL) tmp = ARRAY_2D(NFLX, NFLX, double);

  u[0]  = EXPAND(q[VX1]*q[VX1], + q[VX2]*q[VX2], + q[VX3]*q[VX3]);
  u[0]  = 1.0/sqrt(1.0 - u[0]);
  vB    = EXPAND(q[VX1]*q[BX1], + q[VX2]*q[BX2], + q[VX3]*q[BX3]);
  Bmag2 = EXPAND(q[BX1]*q[BX1], + q[BX2]*q[BX2], + q[BX3]*q[BX3]);

  EXPAND(u[1] = u[0]*q[VXn];  ,
         u[2] = u[0]*q[VXt];  ,
         u[3] = u[0]*q[VXb];)

  b[0] = u[0]*vB; 
  EXPAND(b[1] = q[BXn]/u[0] + b[0]*q[VXn];  ,
         b[2] = q[BXt]/u[0] + b[0]*q[VXt];  ,
         b[3] = q[BXb]/u[0] + b[0]*q[VXb];)
  b2  = Bmag2/(u[0]*u[0]) + vB*vB;

/* -- thermodynamics relations -- */

  eta   = q[RHO] + g_gamma/(g_gamma - 1.0)*q[PRS];  
  E     = eta + b2;  /*  -- this is the total enthalpy -- */
  eta_1 = 1.0/eta;

  a     = q[RHO];
  s     = q[PRS]/pow(a,g_gamma);  
  a     = 1.0/g_gamma;
  e_P   = a*q[RHO]/q[PRS] + 1.0/(g_gamma - 1.0);
  rho_P =  a*q[RHO]/q[PRS];
  rho_S = -a*q[RHO]/s;

/* ---------------------------------------
    Find eigenvalues and store their
    content into lambda[k]
   --------------------------------------- */

  MAGNETOSONIC (q, cs2, h, lambda);
  lambda[KALFVP] = (b[1] + u[1]*sqrt(E))/(b[0] + u[0]*sqrt(E));
  lambda[KALFVM] = (b[1] - u[1]*sqrt(E))/(b[0] - u[0]*sqrt(E));
  lambda[KENTRP] = q[VXn];
  #ifdef KDIVB
  #if DIVB_CONTROL != EIGHT_WAVES
    lambda[KDIVB]  = 0.0;
  #endif
  #endif

/* ---------------------------------------
            spot degeneracies 
   --------------------------------------- */

  for (k = 0; k < NFLX; k++) degen2[k] = 0;

  B = sqrt(Bmag2);
  if (fabs(q[BXn]) < 1.e-9*B){  
    degen1 = 1;
  }
  if (fabs(lambda[KALFVM] - lambda[KFASTM]) < 1.e-6) 
    degen2[KFASTM] = 1;

  if (fabs(lambda[KALFVP] - lambda[KFASTP]) < 1.e-6) 
    degen2[KFASTP] = 1;

  if (fabs(lambda[KALFVM] - lambda[KSLOWM]) < 1.e-6) 
    degen2[KSLOWM] = 1;

  if (fabs(lambda[KALFVP] - lambda[KSLOWP]) < 1.e-6) 
    degen2[KSLOWP] = 1;

  

#if DEBUG_MODE == YES
PRINT_STATE(q,lambda,cs2,h);
#endif
  
/* ----------------------------------------
         ~ Magnetosonic waves ~
           ------------------

    Here "a", "B" depends on the wave
    speed. For efficiency, we split 
    l^\mu = l[mu]  and  f^\mu = f[mu] into 

    l[mu] = phi[mu] + la[mu]*a + B*lB[mu]
    f[mu] = fa[mu]*a + fB[mu]*B
   ----------------------------------------  */


  for (mu = 0; mu <= COMPONENTS; mu++) {
    la[mu] = (eta - e_P*b2)*u[mu]*eta_1;
    lB[mu] =  b[mu]*eta_1;
    fa[mu] =  b[mu]*e_P*eta_1;
    fB[mu] = -u[mu]*eta_1;
  }

  /* ----------------------------
          ~ Fast waves ~
     ---------------------------- */

  for (k = 0; k < NFLX; k++){
    if (k != KFASTM && k != KFASTP) continue;

    a = u[0]*(q[VXn] - lambda[k]);
    B = -b[0]*lambda[k] + b[1];
    A = E*a*a - B*B;

#if DEBUG_MODE == YES
printf ("FAST, a = %12.6e, A = %12.6e, B = %12.6e\n",a,A,B);
#endif

    phi[0] = lambda[k];
    EXPAND(phi[1] = 1.0;, phi[2] = 0.0;, phi[3] = 0.0;)

    if (degen2[k]) {
      for (mu = 0; mu <= COMPONENTS; mu++) {
        d[mu] = B*b[mu] - b2*(phi[mu] + a*u[mu]);
      }
    }else {
      for (mu = 0; mu <= COMPONENTS; mu++) {
        l[mu] = phi[mu] + la[mu]*a + lB[mu]*B;
        f[mu] =           fa[mu]*a + fB[mu]*B;
        d[mu] =   a*(B*f[mu] - a*l[mu]) 
                + (B*B - e_P*b2*a*a)*(phi[mu] + 2.0*a*u[mu])*eta_1;
      }
    }

    for (mu = 0; mu <= COMPONENTS; mu++) {
      Ru[  mu] = a*d[mu];
      Ru[4+mu] = B*d[mu] + a*A*f[mu];
    }
    Ru[8] = a*a*A;
    Ru[9] = 0.0;

  /* -- make physical vector -- */

    RR[RHO][k] = rho_P*Ru[8] + rho_S*Ru[9];
    EXPAND(RR[VXn][k] = (-q[VXn]*Ru[0] + Ru[1])/u[0];  ,
           RR[VXt][k] = (-q[VXt]*Ru[0] + Ru[2])/u[0];  ,
           RR[VXb][k] = (-q[VXb]*Ru[0] + Ru[3])/u[0];)
    RR[PRS][k] = Ru[8];

    EXPAND(                                                                ,
           RR[BXt][k] = b[2]*Ru[0] - b[0]*Ru[2] - u[2]*Ru[4] + u[0]*Ru[6];  ,
           RR[BXb][k] = b[3]*Ru[0] - b[0]*Ru[3] - u[3]*Ru[4] + u[0]*Ru[7];)
  }

  /* ----------------------------
          ~ Slow waves ~
     ---------------------------- */

/* -----------------------------------------
    When Bx -> 0 redefine b^\mu
   ----------------------------------------- */

  if (degen1){
    B     = sqrt(Bmag2);
    if (Bmag2 < 1.e-9*q[PRS]){
      betay = betaz = 1.0/sqrt(2.0);
    }else{
      betay = q[BXt]/B;
      betaz = q[BXb]/B;
    }

    A    = q[VXt]*betay + q[VXb]*betaz;
    b[0] =             u[0]*A;
    b[1] =             u[1]*A;
    b[2] = betay/u[0] + u[2]*A;
    b[3] = betaz/u[0] + u[3]*A;

    b2   = (betay*betay + betaz*betaz)/(u[0]*u[0]) + A*A;

#if DEBUG_MODE == YES
printf ("!!!! Degen 1: Bx -> 0,  b^2 = %12.6e\n", b2);
#endif
  } 

  for (k = 0; k < NFLX; k++){
    if (k != KSLOWP && k != KSLOWM) continue;

    if (degen1) { /* -- when Bx->0 also a, A, B -> 0 -- */

      if (k == KSLOWM) {
        for (mu = 0; mu <= COMPONENTS; mu++) {
          Ru[mu]   = -b[mu];
          Ru[mu+4] =  b[mu];
        }
      }else{
        for (mu = 0; mu <= COMPONENTS; mu++) {
          Ru[mu]   = b[mu];
          Ru[mu+4] = b[mu];
        }
      }

      Ru[8] = - b2*sqrt(Bmag2);
      Ru[9] = 0.0;  
    }else{

      a = u[0]*(q[VXn] - lambda[k]);
      B = -b[0]*lambda[k] + b[1];
      A = E*a*a - B*B;

#if DEBUG_MODE == YES
printf ("SLOW, a = %12.6e, A = %12.6e, B = %12.6e\n",a,A,B);
#endif
      phi[0] = lambda[k];
      EXPAND(phi[1] = 1.0;, phi[2] = 0.0;, phi[3] = 0.0;)

      if (degen2[k]) {
        A = 0.0;
        for (mu = 0; mu <= COMPONENTS; mu++) {
          d[mu] = B*b[mu] - b2*(phi[mu] + a*u[mu]);
        }
      }else {
        for (mu = 0; mu <= COMPONENTS; mu++) {
          l[mu] = phi[mu] + la[mu]*a + lB[mu]*B;
          f[mu] =           fa[mu]*a + fB[mu]*B;
          d[mu] =   a*(B*f[mu] - a*l[mu]) 
                + (B*B - e_P*b2*a*a)*(phi[mu] + 2.0*a*u[mu])*eta_1;
        }
      }

      for (mu = 0; mu <= COMPONENTS; mu++) {
        Ru[  mu] = a*d[mu];
        Ru[4+mu] = B*d[mu] +  a*A*f[mu];
      }
      Ru[8] = a*a*A;
      Ru[9] = 0.0;
    }

  /* -- make physical vector -- */

    RR[RHO][k] = rho_P*Ru[8] + rho_S*Ru[9];
    EXPAND(RR[VXn][k] = (-q[VXn]*Ru[0] + Ru[1])/u[0];  ,
           RR[VXt][k] = (-q[VXt]*Ru[0] + Ru[2])/u[0];  ,
           RR[VXb][k] = (-q[VXb]*Ru[0] + Ru[3])/u[0];)
    RR[PRS][k] = Ru[8];

    EXPAND(                                                                ,
           RR[BXt][k] = b[2]*Ru[0] - b[0]*Ru[2] - u[2]*Ru[4] + u[0]*Ru[6];  ,
           RR[BXb][k] = b[3]*Ru[0] - b[0]*Ru[3] - u[3]*Ru[4] + u[0]*Ru[7];)
  }
    
/* -----------------------------------------------------
                    ~ Alfven waves ~
                      ------------

    They are defined through the inner tensor product

    d^a = eps^{abcd}\phi_b u_c b_d

    Since \phi_a = (-lambda,1), the only non-vanishing 
    components, for a=0,1,2,3 are:

    d^0 = eps^{0123} \phi_1 (u_2b_3 - u_3b_2) =  
        =                   (u^2b^3 - u^3b^2)

    d^1 = eps^{1023} \phi_0 (u_2b_3 - u_3b_2) = 
        =            lambda*(u^2b^3 - u^3b^2)

    d^2 =   eps^{2013} \phi_0 (u_1b_3 - u_3b_1)
          + eps^{2103} \phi_1 (u_0b_3 - u_3b_0) =
        =   -lambda*(u_1b_3 - u_3b_1) + B3
       
    d^3 =   eps^{3012} \phi_0 (u_1b_2 - u_2b_1)
          + eps^{3102} \phi_1 (u_0b_2 - u_2b_0) =
        =    lambda*(u_1b_2 - u_2b_1) - B2

    Where 

    \eps^{0123} = \eps^{2013} = \eps{3102} =  1
    \eps^{1023} = \eps^{2103} = \eps{3012} = -1
 
    and B^k = b^ku^0 - u^kb^0 is the lab frame field.
   -----------------------------------------------------  */

  for (k = 0; k < NFLX; k++){
    if (k != KALFVP && k != KALFVM) continue;

    a = u[0]*(q[VXn] - lambda[k]);
    B = -b[0]*lambda[k] + b[1];
    phi[0] = lambda[k];
    EXPAND(phi[1] = 1.0;, phi[2] = 0.0;, phi[3] = 0.0;)

    Ru[0] =            (u[2]*b[3] - u[3]*b[2]);
    Ru[1] =  lambda[k]*(u[2]*b[3] - u[3]*b[2]);
    Ru[2] = -lambda[k]*(u[1]*b[3] - u[3]*b[1])
            +          (u[0]*b[3] - u[3]*b[0]);
    Ru[3] =  lambda[k]*(u[1]*b[2] - u[2]*b[1])
            -          (u[0]*b[2] - u[2]*b[0]);

    s = B/a;
    if (degen1 && k == KALFVP) s =  1.0;
    if (degen1 && k == KALFVM) s = -1.0;

    Ru[4] = Ru[0]*s;
    Ru[5] = Ru[1]*s;
    Ru[6] = Ru[2]*s;
    Ru[7] = Ru[3]*s;
    Ru[8] = Ru[9] = 0.0;

  /* -- make physical vector -- */

    RR[RHO][k] = 0.0;
    EXPAND(RR[VXn][k] = (-q[VXn]*Ru[0] + Ru[1])/u[0];  ,
           RR[VXt][k] = (-q[VXt]*Ru[0] + Ru[2])/u[0];  ,
           RR[VXb][k] = (-q[VXb]*Ru[0] + Ru[3])/u[0];)
    RR[PRS][k] = Ru[8];

    EXPAND(                                                                ,
           RR[BXt][k] = b[2]*Ru[0] - b[0]*Ru[2] - u[2]*Ru[4] + u[0]*Ru[6];  ,
           RR[BXb][k] = b[3]*Ru[0] - b[0]*Ru[3] - u[3]*Ru[4] + u[0]*Ru[7];)
  }

  /* -------------------------
         ~ Entropy wave ~
           ------------
     -------------------------  */

  k = KENTRP;
  RR[RHO][k] = 1.0;

  /* -------------------------
       ~ Divergence wave ~
         ---------------
     -------------------------  */

  #ifdef KDIVB
  #if DIVB_CONTROL != EIGHT_WAVES
   k = KDIVB;
   RR[BXn][k] = 1.0;
  #endif
  #endif
  
/* ----------------------------------------------------
    Now find the left eigenvectors by performing 
    LU decomposition 
   ---------------------------------------------------- */

  for (nv = NFLX; nv--;   ){
  for (k  = NFLX;  k--;   ){
    tmp[nv][k] = RR[nv][k];
  }}

#if DEBUG_MODE == YES
ShowMatrix(RR);
#endif

  if (!LUDecompose(tmp, NFLX, indx, &a)) {
/*    printf ("! EIGENV: Singular eigenvectors\n"); */
    for (nv = NFLX; nv--;   ){
    for (k  = NFLX;  k--;   ){
      RR[nv][k] = LL[k][nv] = (k == nv ? 1.0:0.0);
    }}
    #ifdef KDIVB
     #if DIVB_CONTROL != EIGHT_WAVES
     LL[KDIVB][BXn] = RR[BXn][KDIVB] = 0.0;
     #endif
    #endif
    return;
  }
    
  for (mu = 0; mu < NFLX; mu++){
    for (k = 0; k < NFLX; k++) col[k] = 0.0;
    col[mu] = 1.0;
    LUBackSubst(tmp, NFLX, indx, col);
    for (k = 0; k < NFLX; k++) LL[k][mu] = col[k];
  }

  #ifdef KDIVB
   #if DIVB_CONTROL != EIGHT_WAVES
    LL[KDIVB][BXn] = RR[BXn][KDIVB] = 0.0;
   #endif
  #endif

#if DEBUG_MODE == YES
printf ("-------------------------------\n");
printf ("       INVERSE MATRIX\n");
printf ("-------------------------------\n");
ShowMatrix(LL);
exit(1);
#endif
}
#undef DEBUG_MODE

void PrimToChar (double **Lp, double *v, double *w)
{
  int k, nv;

  for (k = NFLX; k--;   ){
    w[k] = 0.0;
    for (nv = NFLX; nv--;  ){
      w[k] += Lp[k][nv]*v[nv];
    }
  }
}


int PRINT_STATE(double *q, double *lambda, double cs2, double hh)
{

printf ("q[RHO] = %12.6e;\n", q[RHO]);
printf ("q[VXn] = %12.6e;\n", q[VXn]);
printf ("q[VXt] = %12.6e;\n", q[VXt]);
printf ("q[VXb] = %12.6e;\n", q[VXb]);

printf ("q[BXn] = %12.6e;\n", q[BXn]);
printf ("q[BXt] = %12.6e;\n", q[BXt]);
printf ("q[BXb] = %12.6e;\n", q[BXb]);

printf ("q[PRS] = %12.6e;\n", q[PRS]);
printf ("cs2 = %12.6e; hh = %12.6e;\n",cs2,hh);


printf ("fm, %d  %12.6e\n", KFASTM, lambda[KFASTM]);
printf ("am, %d  %12.6e\n", KALFVM, lambda[KALFVM]);
printf ("sm, %d  %12.6e\n", KSLOWM, lambda[KSLOWM]);
printf ("en, %d  %12.6e\n", KENTRP, lambda[KENTRP]);
printf ("sp, %d  %12.6e\n", KSLOWP, lambda[KSLOWP]);
printf ("ap, %d  %12.6e\n", KALFVP, lambda[KALFVP]);
printf ("fp, %d  %12.6e\n", KFASTP, lambda[KFASTP]);

}

#endif /* 1==0 */

