/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Functions for QR decomposition and matrix inversion.
  \author B. Vaidya (bvaidya@unito.it)
  \date   September 11, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ******************************************** */
void RSolve(double **a, int n, double *d, double *b){ 
/*
 * Solves set of n Linear Equations. R.x = b. 
 * Where is R is upper Triangular Matrix stored in a
 * and d. The matrix a and vector d are outputs from 
 * QRDecompse and b[] is the input vector.  
 *
 ********************************************* */
  int i,j;
  double sum;
  
  b[n-1] /= d[n-1];
  for(i=n-2;i>=0;i--){
    sum = 0.0;
    for (j = i+1; j < n; j++){
      sum += a[i][j]*b[j];
    }
    b[i] = (b[i] - sum)/d[i];
  }
}

/* ********************************************************************* */
void QRDecompose (double **a, int n, double *c, double *d, int *sing)
/*   
 * Perform QR decomposition of matrix a. The upper Triangular matrix R 
 * is returned in a except the diagonal elements which are stored in array d.
 * The orthogonal matrix Q is represented as a product of n-1 Householder
 * matrices. Qj = 1.0 - uj \cross uj/cj, where the non-zero components of uj
 * are returned in a[i][j]. 
 * sing returns 1 if there is singularity else on normal exist it returns 0.
 * 
 *********************************************************************** */
{
  int i, j, k;
  double scale, sigma, sum, tau, temp;
  
  *sing = 0;
  for (k = 0; k < n-1; k++) {
    scale = 0.0;
    for (i = k; i < n; i++){
      if ((temp = fabs (a[i][k])) > scale)
        scale = temp;
    }
    
    if (scale == 0.0) {
      *sing = 1;
      c[k] = d[k] = 0.0;
    }else{
     for (i = k; i < n; i++){
       a[i][k] /= scale;
     }
     sum = 0.0;
     for (i = k; i < n; i++) {
       sum  += a[i][k] * a[i][k];
     }

     if (a[k][k] >= 0){
       sigma = fabs(sqrt(sum));
     }else{
       sigma = -fabs(sqrt(sum));
     }

   a[k][k] += sigma;
     c[k] = sigma*a[k][k];
     d[k] = -sigma*scale;
     for (j=k+1;j< n;j++){
       sum = 0.0;
       for (i = k; i < n; i++) {
         sum  += a[i][k] * a[i][j];
       }
       tau = sum /c[k];
       for (i = k; i < n; i++) {
         a[i][j]  -= tau * a[i][k];
        }
      }
    }
  }
  d[n-1] = a[n-1][n-1];
  if (d[n-1] == 0.0) *sing = 1;
}




/* ******************************************** */
void QRSolve(double **a, int n, double *c, double *d, double *b){ 
/*
 * Solves set of n Linear Equations. A.x = b. 
 * The matrix a[4][4] and vector d are outputs from 
 * QRDecompse, thus they represent the fractured R matrix 
 * and b[] is the input RHS vector.  
 *
 ********************************************* */
 int i, j;
 double sum, tau;
 for(j=0;j<n-1;j++){ /* Form Q^T . b */ 
  sum = 0.0;
  for(i=j;i<n;i++){
   sum += a[i][j]*b[i];
   }
  tau = sum/c[j];
  for(i=j;i<n;i++){
   b[i] -= tau*a[i][j];
   }
 }
 RSolve(a, n, d, b);
}

/* ******************************************** */
void QRUpdate(double **r, double **qt, int n, double *u, double *v){ 
/*
 * Solves set of n Linear Equations. R.x = b. 
 * Where is R is upper Triangular Matrix stored in a
 * and d. The matrix a and vector d are outputs from 
 * QRDecompse and b[] is the input vector.  
 *
 ********************************************* */
int i, j, k;
for (k=n-1;k>=0;k--){ /* Find largest k such that u[k] neq 0 */
 if (u[k]) break;
} 

for (i = k-1; i>=0; i--){
 rotate(r,qt, n, i, u[i], -u[i+1]);
 if (u[i] == 0.0) u[i] = fabs(u[i+1]);
 else if (fabs(u[i]) > fabs(u[i+1]))
   u[i] = fabs(u[i])*sqrt(1.0 + (u[i+1]/u[i])*(u[i+1]/u[i]));
 else
   u[i] = fabs(u[i+1])*sqrt(1.0 + (u[i]/u[i+1])*(u[i]/u[i+1]));
}

for (j=0;j<n; j++) r[0][j] += u[0]*v[j];
for (i=0;i<k; i++)
  rotate(r,qt,n,i,r[i][i], -r[i+1][i]);
}

/* ******************************************** */
void rotate(double **r, double **qt, int n, int i,  double a, double b){ 
/*
 * Jacobi rotation.   
 *
 ********************************************* */
int j;
double c, s, w, y, fact;

if (a == 0.0){
  c = 0.0;
  s = (b >= 0.0 ? 1.0 : -1.0);
 }else if (fabs(a) > fabs(b)){
  fact = b/a;
  if (a >= 0.0){
    c = fabs(1.0/sqrt(1.0 + (fact*fact)));
 }else{
    c = -fabs(1.0/sqrt(1.0 + (fact*fact)));
 }
  s = fact * c;
 }else{
  fact = a/b;
  if (b >= 0.0) {
    s = fabs(1.0/sqrt(1.0 + (fact*fact)));
 }else {
    s = -fabs(1.0/sqrt(1.0 + (fact*fact)));
  }
  c = fact * s;
 }

for (j=i;j<n;j++){
 y = r[i][j];
 w = r[i+1][j];
 r[i][j] = c*y - s*w;
 r[i+1][j] = s*y + c*w;
}

for (j=0;j<n;j++){
 y = qt[i][j];
 w = qt[i+1][j];
 qt[i][j] = c*y - s*w;
 qt[i+1][j] = s*y + c*w;
}

}


