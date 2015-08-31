#include"pluto.h"

/* *********************************************************************** */
double MP5_Reconstruct(double *F, double dx, int j)
/*  
 *
 * 
 *  
 *  
 *
 *
 ************************************************************************* */
{
  double f, d2, d2p, d2m; 
  double dMMm, dMMp;
  double scrh1,scrh2, fmin, fmax; 
  double fAV, fMD, fLC, fUL, fMP;
  static double alpha = 4.0, epsm = 1.e-12;

  f  = 2.0*F[j-2] - 13.0*F[j-1] + 47.0*F[j] + 27.0*F[j+1] - 3.0*F[j+2];
  f /= 60.0;   

  fMP = F[j] + MINMOD(F[j+1]-F[j],alpha*(F[j]-F[j-1]));

  if ((f - F[j])*(f - fMP) <= epsm) return f;

  d2m = F[j-2] + F[j  ] - 2.0*F[j-1];    /* -- Eq. (2.19) -- */
  d2  = F[j-1] + F[j+1] - 2.0*F[j];
  d2p = F[j  ] + F[j+2] - 2.0*F[j+1];    /* -- Eq. (2.19) -- */

  scrh1 = MINMOD(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = MINMOD(d2, d2p);
  dMMp  = MINMOD(scrh1,scrh2);   /* -- Eq. (2.27) -- */

  scrh1 = MINMOD(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = MINMOD(d2, d2m);
  dMMm  = MINMOD(scrh1,scrh2);   /* -- Eq. (2.27) -- */

  fUL = F[j] + alpha*(F[j] - F[j-1]);   /* -- Eq. (2.8) -- */
  fAV = 0.5*(F[j] + F[j+1]);        /* -- Eq. (2.16) -- */
  fMD = fAV - 0.5*dMMp; /* -- Eq. (2.28) -- */
  fLC = 0.5*(3.0*F[j] - F[j-1]) + 4.0/3.0*dMMm;  /* -- Eq. (2.29) -- */

  scrh1 = MIN(F[j], F[j+1]); scrh1 = MIN(scrh1, fMD);
  scrh2 = MIN(F[j], fUL);    scrh2 = MIN(scrh2, fLC);
  fmin  = MAX(scrh1, scrh2);  /* -- Eq. (2.24a) -- */

  scrh1 = MAX(F[j], F[j+1]); scrh1 = MAX(scrh1, fMD);
  scrh2 = MAX(F[j], fUL);    scrh2 = MAX(scrh2, fLC);
  fmax  = MIN(scrh1, scrh2);  /* -- Eq. (2.24b) -- */

  f = Median(f, fmin, fmax); /* -- Eq. (2.26) -- */
  return f;
}

/* ************************************************************* */
double PPM_Reconstruct (double *v, double dx, int i)
/*
 *
 *
 *
 *************************************************************** */
{
  double dvm2, dvm1, dvp1, dvp2;
  double dvc, Sm1, Sp1, Sp2, SM;
  double ap, am, vp, vm;

  dvm2 = v[i-1] - v[i-2];
  dvm1 = v[i]   - v[i-1];
  dvp1 = v[i+1] - v[i];
  dvp2 = v[i+2] - v[i+1];

  dvc = 0.5*(dvm2 + dvm1);
  SM  = 2.0*MINMOD(dvm2, dvm1);
  Sm1 = MINMOD(dvc, SM);

  dvc = 0.5*(dvm1 + dvp1);
  SM  = 2.0*MINMOD(dvm1, dvp1);
  Sp1 = MINMOD(dvc, SM);

  dvc = 0.5*(dvp2 + dvp1);
  SM  = 2.0*MINMOD(dvp2, dvp1);
  Sp2 = MINMOD(dvc, SM);

  vp = 0.5*(v[i] + v[i+1]) - (Sp2 - Sp1)/6.0;
  vm = 0.5*(v[i] + v[i-1]) - (Sp1 - Sm1)/6.0;

  ap = vp - v[i];
  am = vm - v[i];

  if (ap*am >= 0.0) {
    ap = am = 0.0;
  }else{
    if (fabs(ap) >= 2.0*fabs(am)) ap = -2.0*am;
    if (fabs(am) >= 2.0*fabs(ap)) am = -2.0*ap;
  }
  vp = v[i] + ap;
  return vp;
/*  vm = v[i] + am; */
}

/* ************************************************************* */
double LIMO3_Reconstruct(double *v, double dx,  int i)
/*
 *
 *
 *
 *************************************************************** */
{
  double dvp, dvm, r = 1.;
  double a,b,c,q, th, lim;
  double eta, psi, eps = 1.e-12;

  dvp = v[i+1] - v[i];
  dvm = v[i] - v[i-1];

  th  = dvm/(dvp + 1.e-16);

  q = (2.0 + th)/3.0;
  a = MIN(1.5,2.0*th);
  a = MIN(q,a);
  b = MAX(-0.5*th,a);
  c = MIN(q,b);
  psi = MAX(0.0,c);

  eta = r*dx;
  eta = (dvm*dvm + dvp*dvp)/(eta*eta);
  if ( eta <= 1.0 - eps) {
    lim = q;
  }else if (eta >= 1.0 + eps){
    lim = psi;
  }else{
    psi =   (1.0 - (eta - 1.0)/eps)*q
          + (1.0 + (eta - 1.0)/eps)*psi;
    lim = 0.5*psi;
  }
  return (v[i] + 0.5*dvp*lim);
}
/* *********************************************************************** */
double WENOZ_Reconstruct(double *F, double dx,  int j)
/*  
 *
 * 
 *  Provide interface values F_{i+1/2} using 
 *  WENO-Z reconstruction.
 *   
 *
 *
 ************************************************************************* */
#define eps 1.e-6
{
  double IS0, IS1, IS2, IS3, IS4, IS5;
  double a0, a1, a2;
  double f0, f1, f2;
  double w0, w1, w2;
  double sum_a, f;
  static double thirteen_12 = 13.0/12.0;
  double b0, b1, b2,t5;

  a0 = F[j-2] - 2.0*F[j-1] +     F[j]; 
  a1 = F[j-2] - 4.0*F[j-1] + 3.0*F[j]; 
  b0 = thirteen_12*a0*a0 + 0.25*a1*a1;

  a0 = F[j-1] - 2.0*F[j] + F[j+1]; 
  a1 = F[j-1] - F[j+1]; 
  b1 = thirteen_12*a0*a0 + 0.25*a1*a1;

  a0 =     F[j] - 2.0*F[j+1] + F[j+2]; 
  a1 = 3.0*F[j] - 4.0*F[j+1] + F[j+2]; 
  b2 = thirteen_12*a0*a0 + 0.25*a1*a1;

  t5 = fabs(b0-b2);

  a0 = 1.0*(1.0 + t5/(b0 + 1.e-40));
  a1 = 6.0*(1.0 + t5/(b1 + 1.e-40));
  a2 = 3.0*(1.0 + t5/(b2 + 1.e-40));

  sum_a = 1.0/(a0 + a1 + a2);
  w0 = a0*sum_a;
  w1 = a1*sum_a;
  w2 = a2*sum_a;

  f0 = (2.0*F[j-2] - 7.0*F[j-1] + 11.0*F[j]);
  f1 = (   -F[j-1] + 5.0*F[j]   +  2.0*F[j+1]);
  f2 = (2.0*F[j]   + 5.0*F[j+1] -      F[j+2]);

  f = (w0*f0 + w1*f1 + w2*f2)/6.0;
  return(f);
}

/* *********************************************************************** */
double WENO3_Reconstruct(double *F, double dx, int j)
/*  
 *
 * 
 *  Provide interface values F_{i+1/2} using 
 *  third-order WENO reconstruction.
 *   
 *
 *
 ************************************************************************* */
{
  double a0, b0, w0, f0, a1, b1, w1, f1;
  double tau, dx2;

  dx2 = dx*dx;

  b0 = F[j+1] - F[j];
  b1 = F[j-1] - F[j];

  b0 = b0*b0;
  b1 = b1*b1;

/* -- traditional version -- */

/*
  a0 = eps + b0;
  a1 = eps + b1;
  
  a0 = 2.0/(a0*a0);
  a1 = 1.0/(a1*a1);
*/
/* -- improved version -- */

  tau = (F[j+1] - 2.0*F[j] + F[j-1]);
  tau = tau*tau;

  a0 = 2.0*(1.0 + tau/(dx2 + b0));
  a1 = 1.0*(1.0 + tau/(dx2 + b1));

  w0 = a0/(a0 + a1);
  w1 = a1/(a0 + a1);

  f0 =  0.5*( F[j+1] +     F[j]);
  f1 =  0.5*(-F[j-1] + 3.0*F[j]);
  return (w0*f0 + w1*f1);
}
#undef eps


/* *********************************************************************** */
double LIN_Reconstruct(double *F, double dx, int j)
/*  
 *
 * 
 *  Provide interface values F_{i+1/2} using 
 *  simple slope-limited linear reconstruction.
 *   
 *
 *
 ************************************************************************* */
{
  double dFp, dFm, dF;

  dFp = F[j+1] - F[j];
  dFm = F[j] - F[j-1];

  dF = MINMOD(dFp, dFm);
  return (F[j] + 0.5*dF);
}


/* ************************************************************* */
double Median (double a, double b, double c)
/*
 *
 *
 *
 *************************************************************** */
{
  return (a + MINMOD(b-a,c-a));
}
