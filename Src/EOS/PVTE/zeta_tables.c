#include "pluto.h"

static void MakeZetaTables(double *, double *, int);

/* Vibrational and Rotational Temperatures for molecular hydrogen  */
#define THETA_V  6140.0
#define THETA_R    85.5

/* Switch to choose the molecular hydrogen spin states.
 * 0 : Only Para hydrogen , 
 * 1 : Equilibrium (DEFAULT), 
 * 2 : Ortho to para ratio of 3:1.                                 */
#define ORTHO_PARA_MODE 1

/* ********************************************************************* */
void GetFuncDum(double T, double *funcdum_val)
/*!
 *  Interpolate the value of a function of \c zetaR from the table 
 *  that is used to estimate the value of EH2. 
 * 
 * \param [in]   T             Value of temperature in kelvin.
 * \param [out] *funcdum_val   Pointer to the value of function of zetaR
 * 
 * \return This function has no return value.
 *
 * \b  Reference:\n
 *     D'Angelo, G. et al, ApJ 778 2013.
 *********************************************************************** */
{
  int    klo, khi, kmid;
  int nsteps = 5000;
  double mu, Tmid, dT;
  static double *lnT, *funcdum;
  double y, dy;
  int indx;

/* -------------------------------------------
    Make table on first call
   ------------------------------------------- */

  if (lnT == NULL){    
    lnT    = ARRAY_1D(nsteps, double);
    funcdum = ARRAY_1D(nsteps, double);
    MakeZetaTables(lnT, funcdum, nsteps);
  }
  y = log(T);

/* -------------------------------------------------
    Since the table has regular spacing in log T,
    we divide by the increment to find the nearest
    node in the table.
   ------------------------------------------------- */
  
  if (y > lnT[nsteps-2]) {
    *funcdum_val = funcdum[nsteps-2];
  } else if (y < lnT[0]) {
    *funcdum_val = funcdum[0];
  } else{
    dy   = lnT[1] - lnT[0];
    indx = floor((y - lnT[0])/dy);

    if (indx >= nsteps || indx < 0){
      print1 ("! GetFuncDum: indx out of range, indx = %d\n",indx);
      print1 ("! T = %12.6e\n",T);
      QUIT_PLUTO(1);
    }
    *funcdum_val = (funcdum[indx]*(lnT[indx+1] - y)
                  + funcdum[indx+1]*(y - lnT[indx]))/dy;
  }  
}

/* ********************************************************************* */
void MakeZetaTables(double *lnT, double *funcdum, int nsteps)
/*!
 *  Compute tables from iterative summation involed to estimate 
 *  partition function of parahydrogen (\f$ \zeta_P\f$) and
 *  orthohydrogen (\f$ \zeta_O\f$) and their derivatives with respect
 *  to temperature. Then further estimate \f$\zeta_R\f$ and finally 
 *  the function of \f$ \zeta_R\f$ that goes into the estimation of 
 *  EH2.
 *  The partition function and its derivative can be written as 
 *  \f[ 
 *    \left\{\begin{array}{l}
 *      \zeta  = \DS \sum_i a_ie^{-b_i/T} \\ \noalign{\medskip}
 *      \zeta' = \DS \frac{1}{T^2}\sum_i a_ib_ie^{-b_i/T} 
 *    \end{array}\right.
 *  \f] 
 *  where <tt>i=0,2,4,...</tt> is even for parahydrogen while 
 *  <tt>i=1,3,5,...</tt> is odd for orthoyhdrogen.
 *  One can see that in the limit of zero temperatures, \f$ \zeta_P \to 1,
 *  \, \zeta'_P \to 0\f$ while \f$\zeta_O\to 0,\,\zeta'_O\to0\f$.
 *  Then \f$ \zeta'_O/\zeta_O\f$ is ill-defined in the low temperature limit
 *  since it becomes \c 0/0.
 *  To avoid this, we rewrite \c zetaO by extracting the first term from 
 *  the summation:
 *  \f[
 *    \left\{\begin{array}{l}
 *      \zeta_O  = \DS e^{-b_1/T}\sum_i a_ie^{-\Delta b_i/T} \\ \noalign{\medskip}
 *      \zeta'_O = \DS \frac{1}{T^2}e^{-b_1/T}
 *                   \left[ b_1 \sum_i a_ie^{-\Delta b_i/T} 
 *                             +\sum_i a_i\Delta b_ie^{-\Delta b_i/T}\right]
 *    \end{array}\right.
 *  \f]
 *  where \f$\Delta b_i = b_i-b_1\f$. 
 *  Taking the ratio:
 *  \f[
 *    \frac{\zeta'_O}{\zeta_O} = \frac{1}{T^2}
 *      \frac{   b_1 \sum_i a_ie^{-\Delta b_i/T} 
 *            + \sum_i a_i\Delta b_ie^{-\Delta b_i/T}}
 *           {\sum_i a_ie^{-\Delta b_i/T}}
 *     \qquad\Longrightarrow\qquad
 *    \frac{\zeta'_O}{\zeta_O} - \frac{b_1}{T^2}= \frac{1}{T^2}
 *      \frac{\sum_i a_i\Delta b_ie^{-\Delta b_i/T}}
 *           {\sum_i a_ie^{-\Delta b_i/T}}  \to 0 \quad{\rm as}\quad T\to 0
 *  \f]
 *  The expression on the right (that appears in the computation of 
 *  \f$\zeta'_R\f$) is now well-behaved and it reproduces the
 *  low temperature limit correctly.
 *  You can convince yourself by playing with the following MAPLE script
 *  \code
     restart;
     a[0]  := 1; a[1] := 3;       a[2] := 5;       a[3] := 7;
     b[0] := 0; b[1] := 2*theta; b[2] := 6*theta; b[3] := 12*theta;
 
     zeta[P]  := a[0]*exp(-b[0]*y) + a[2]*exp(-b[2]*y);
     dzeta[P] := y^2*(a[0]*b[0]*exp(-b[0]*y) + a[2]*b[2]*exp(-b[2]*y));
 
     zeta[O]  := a[1]*exp(-b[1]*y) + a[3]*exp(-b[3]*y);
     dzeta[O] := y^2*(a[1]*b[1]*exp(-b[1]*y) + a[3]*b[3]*exp(-b[3]*y));
     rP := dzeta[P]/zeta[P];
     rO := dzeta[O]/zeta[O];
 
     rOminus := rO - b[1]*y^2;
     simplify(rOminus);
 *  \endcode
 *
 *
 *  \param [in]  lnT    Array of Logrithmic values of Gas 
 *                      temperatures from 0.01 K to 10^12 K.
 *  \param [in] nsteps  Number of equal spacings in T. 
 *  \param [out] funcdum The function of zetaR that goes into EH2. 
 *  
 *  \return This function has no return value
 *
 *  \b Reference: \n
 *     - D'Angelo, G. et. al ApJ 778, 2013 (Eq. 23) 
 *
 *********************************************************************** */
{
  int    i,j;
  double Temp0 = 0.01*T_CUT_RHOE;
  double Tmax  = 1.0e12;
  double dy = log(Tmax/Temp0)*(1./nsteps);
  double dT = Temp0*exp(dy); 
  double T, a, b, b1, scrh, inv_T2;
  double zetaP, dzetaP, zetaO, dzetaO, zetaR, dzetaR;
  double dum1, dum2, dum3;
  double alpha, beta, gamma;
  double dzO_zO_m, db, sum1, sum2;

  print1 ("> MakeZetaTables(): generating Zeta tables...\n");

  if (ORTHO_PARA_MODE == 0){ 
    alpha = 1.0; beta = 0.0; gamma = 0.0;
  }else if(ORTHO_PARA_MODE == 2){ 
    alpha = 0.25; beta = 0.75; gamma = 0.0;
  }else{
    alpha = 1.0; beta = 0.0; gamma = 1.0;
  }

  b1 = 2.0*THETA_R;
  for(j = 0; j < nsteps; j++){
    T = Temp0*exp(j*dy);  
    inv_T2 = 1.0/(T*T);
    zetaO = zetaP = dzetaP = dzetaO = 0.0;
    dzO_zO_m = sum1 = sum2 = 0.0;     
    for(i = 0; i <= 10000; i++){
      a = 2*i + 1;
      b = i*(i + 1)*THETA_R;
      if (i%2 == 0){
        scrh    = a*exp(-b/T);
        zetaP  += scrh; 
        dzetaP += scrh*b;
      }else{
        db    = b - b1; 
        scrh  = a*exp(-db/T);
        sum1 += scrh;
        sum2 += scrh*db;
      }
    }
    dzetaP *= inv_T2;

    zetaO  = exp(-b1/T)*sum1;
    dzetaO = exp(-b1/T)*(b1*sum1 + sum2)*inv_T2;

    dzO_zO_m = sum2/sum1*inv_T2; /* = zeta'(O)/zeta(O) - 2*theta/T^2 */

    lnT[j]  = log(T);
    
  /* -----------------------------------------
      Compute table 
     ----------------------------------------- */
    
    scrh   = zetaO*exp(2.0*THETA_R/T);
    
    zetaR  = pow(zetaP,alpha)*pow(scrh,beta) + 3.0*gamma*zetaO;
    dzetaR = (zetaR - 3.0*gamma*zetaO)*(alpha*(dzetaP/zetaP) + 
				          beta*dzO_zO_m)  + 3.0*gamma*dzetaO;
    dum1  = THETA_V/T;                            
    dum2  = dum1*exp(-dum1)/(1.0 - exp(-dum1));  
    dum3  = (T/zetaR)*dzetaR;
    funcdum[j] = 1.5 + dum2 + dum3;
  } 
}

