/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Miscellaneous math functions.
  \author A. Mignone (mignone@ph.unito.it)
          B. Vaidya 
  \date   Oct 4, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
double BesselJ0(double x)
/*!
 * Returns the Bessel function J0(x) for any real x.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  double ax,z;
  double xx,y,ans,ans1,ans2;      /* Accumulate pol.in double prec. */
   
  if ((ax=fabs(x)) < 8.0) {       /* Direct rational function fit. */
    y   = x*x;
    ans1 =  57568490574.0+y*(-13362590354.0+y*(651619640.7
          + y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2 =  57568490411.0+y*(1029532985.0+y*(9494680.718
          + y*(59272.64853+y*(267.8532712+y*1.0))));
    ans = ans1/ans2;
  } else {                                  /*Fitting function (6.5.9).*/
    z    = 8.0/ax;
    y    = z*z;
    xx   = ax-0.785398164;
    ans1 = 1.0 + y*(-0.1098628627e-2+y*(0.2734510407e-4
           + y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = - 0.1562499995e-1+y*(0.1430488765e-3
           + y*(-0.6911147651e-5+y*(0.7621095161e-6
           - y*0.934945152e-7)));
    ans = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}

/* ********************************************************************* */
double BesselJ1(double x)
/*!
 * Returns the Bessel function J1(x) for any real x.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  double ax,z;
  double xx,y,ans,ans1,ans2;    /* Accumulate polynomials in double precision.*/
  if ((ax=fabs(x)) < 8.0) {     /* Direct rational approximation.*/
  
    y=x*x;
    ans1 =   x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
           + y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2 =  144725228442.0+y*(2300535178.0+y*(18583304.74
          + y*(99447.43394+y*(376.9991397+y*1.0))));
    ans  = ans1/ans2;
  } else {                                  /*Fitting function (6.5.9).*/
    z    = 8.0/ax;
    y    = z*z;
    xx   = ax-2.356194491;
    ans1 = 1.0 +y*(0.183105e-2+y*(-0.3516396496e-4
           + y*(0.2457520174e-5+y*(-0.240337019e-6))));
     ans2 = 0.04687499995+y*(-0.2002690873e-3
           + y*(0.8449199096e-5+y*(-0.88228987e-6
           + y*0.105787412e-6)));
     ans  = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    if (x < 0.0) ans = -ans;
  }
  return ans;
}



/* ********************************************************************* */
double BesselI0(double x)
/*!
 * Returns the modified Bessel function I0(x) for any real x.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  double ax,ans;
  double y;    /* Accumulate polynomials in double precision.*/
  if ((ax=fabs(x)) < 3.75) {     /* Polynomial fit. */
     y = x/3.75;
     y *= y;
     ans  = 1.0 + y*(3.5156229+y*(3.0899424 
            + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1
	    + y*0.45813e-2)))));
  } else {                                  
    y    = 3.75/ax;
    ans = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1
	  + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2
	  + y*(-0.2057706e-1 + y*(0.2635537e-1
	  + y*(-0.1647633e-1 + y*0.392377e-2))))))));
  }
  return ans;
}


/* ********************************************************************* */
double BesselI1(double x)
/*!
 * Returns the modified Bessel function I1(x) for any real x.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  double ax, ans1, ans2, ans;
  double y;    /* Accumulate polynomials in double precision.*/
  if ((ax=fabs(x)) < 3.75) {     /* Polynomial fit. */
     y = x/3.75;
     y *= y;
     ans  = ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
	    +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } else {                                  
    y    = 3.75/ax;
    ans1 = 0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
	   -y*0.420059e-2));
    ans2 = 0.39894228+y*(-0.3988024e-1 + y*(-0.362018e-2
	   + y*(0.163801e-2+y*(-0.1031555e-1+y*ans1))));
    ans = ans2*(exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}


/* ********************************************************************* */
double BesselK0(double x)
/*!
 * Returns the modified Bessel function K0(x) for any real x.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  double I0, ans;
  double y;    /* Accumulate polynomials in double precision.*/
  if (x <= 2.0) {     /* Polynomial fit. */
    y = x*x/4.0;
    I0 = BesselI0(x);
    ans = (-log(x/2.0)*I0)+(-0.57721566+y*(0.42278420
	   + y*(0.23069756+y*(0.3488590e-1
	   + y*(0.262698e-2 +y*(0.10750e-3+y*0.74e-5)))))); 
  } else {                                 
    y    = 2.0/x;
    ans = (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
	  + y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
	  + y*(-0.251540e-2+y*0.53208e-3))))));
  }
  return ans;
}

/* ********************************************************************* */
double BesselK1(double x)
/*!
 * Returns the modified Bessel function K1(x) for any real x.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  double I1, ans;
  double y;    /* Accumulate polynomials in double precision.*/
  if (x <= 2.0) {     /* Polynomial fit. */
    y = x*x/4.0;
    I1 = BesselI1(x);
    ans = (log(x/2.0)*I1)+(1.0/x)*(1.0+y*(0.15443144
	   + y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
	   + y*(-0.110404e-2+y*(-0.4686e-4)))))));
  }else{                                 
    y = 2.0/x;
    ans = (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619 +y*(-0.3655620e-1
	  + y*(0.1504268e-1+y*(-0.780353e-2
	  + y*(0.325614e-2+y*(-0.68245e-3)))))));
  }
  return ans;
}

/* ********************************************************************* */
double BesselKn(int n, double x)
/*!
 * Returns the modified Bessel function Kn(x) for positive x and n >= 2.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  int j;
  double bk,bkm,bkp,tox;
  if (n < 2){
    print("!BesselKn: Index n less than 2");
    QUIT_PLUTO(1);
  }  
  tox = 2.0/x;
  bkm = BesselK0(x);
  bk  = BesselK1(x);
  for (j=1;j<n;j++) {
    bkp = bkm + j*tox*bk;
    bkm = bk;
    bk  = bkp;
  }
  return bk;
}


/* ********************************************************************* */
void QuickSort(int *x,int first,int last)
/*!
 * Sort an integer array x in ascending order.
 * Taken from http://www.cquestions.com/2008/01/c-program-for-quick-sort.html
 *
 * \param [in,out] *x     the array to be sorted.
 * \param [in]     first  the starting element of the array
 * \param [in]     last   the final element of the array
 *
 *********************************************************************** */
{
  int pivot,j,temp,i;

  if (first < last){
    pivot = first;
    i     = first;
    j     = last;

    while(i<j){
      while(x[i]<=x[pivot]&&i<last)
        i++;
      while(x[j]>x[pivot])
        j--;
      if(i<j){
        temp=x[i];
        x[i]=x[j];
        x[j]=temp;
      }
    }

    temp=x[pivot];
    x[pivot]=x[j];
    x[j]=temp;
    QuickSort(x,first,j-1);
    QuickSort(x,j+1,last);
  }
}

/* ********************************************************************* */
void SortArray(double *z, int n)
/*!
 * Sort array elements using straight insertion.
 *
 * \param [in,out]  z   pointer to 1D array in double precision
 * \param [in]      n   the number of elements (starting at 0)
 *
 * \b Reference:
 *    -"Numerical Recipes in C", sect. 8.1
 *********************************************************************** */
{
  int    i, j;
  double f;

  for (j = 1; j < n; j++){
    f = z[j];
    i = j - 1;
    while(i >= 0 && z[i] > f) {
       z[i+1] = z[i];
       i--;
    }
    z[i+1] = f;
  }
}

/* ********************************************************************* */
void VectorCartesianComponents(double *v, double x1, double x2, double x3)
/*!
 * Transform the vector \f$ v = (v_{x1}, v_{x2}, v_{x3})\f$ from
 * the chosen coordinate system (= GEOMETRY) to Cartesian components.
 *  
 * \param [in,out]  *v  the original vector.
 *                      On output \c v is replaced by the three Cartesian
 *                      components.
 * \param [in] x1,x2,x3  the grid coordinates
 *
 *********************************************************************** */
{
  double vx1, vx2, vx3;

  vx1 = v[0];
  vx2 = v[1];
  vx3 = v[2];

  #if GEOMETRY == CARTESIAN

  /* nothing to do */

  #elif GEOMETRY == POLAR

   EXPAND(v[0] = vx1*cos(x2) - vx2*sin(x2);   ,
          v[1] = vx1*sin(x2) + vx2*cos(x2);   ,
          v[2] = vx3;)
  
  #elif GEOMETRY == SPHERICAL

   #if DIMENSIONS == 2
    EXPAND(v[0] = vx1*sin(x2) + vx2*cos(x2);   ,
           v[1] = vx1*cos(x2) - vx2*sin(x2);   ,
           v[2] = vx3;)
   #elif DIMENSIONS == 3 
    v[0] = (vx1*sin(x2) + vx2*cos(x2))*cos(x3) - vx3*sin(x3);
    v[1] = (vx1*sin(x2) + vx2*cos(x2))*sin(x3) + vx3*cos(x3);
    v[2] = (vx1*cos(x2) - vx2*sin(x2));
   #endif
  #endif
}
    


