/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Wave-speeds and characteristic decomposition for the RHD equations.

  This file contains various functions containing Jacobian-related information
  such as characteristic signal speeds, eigenvalues and eigenvector 
  decomposition for the RHD module.
  
  The function MaxSignalSpeed() computes the maximum and minimum 
  characteristic signal velocity for the RHD equations.
  
  The function PrimEigenvectors() calculates left and right eigenvectors
  and the corresponding eigenvalues for the \e primitive form the the
  RHD equations with an adiabatic EoS.

  The function PrimToChar()  compute the matrix-vector multiplcation 
  between the L matrix (containing primitive left eigenvectors) 
  and the vector v. The result containing the characteristic
  variables is stored in the vector w = L.v

  \authors A. Mignone (mignone@ph.unito.it)\n
  \date    Feb 13, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void MaxSignalSpeed (const State *state, double *cmin, double *cmax,
                     int beg, int end)
/*!
 *  Compute the maximum and minimum propagation speeds for the RHD
 *  equations.
 *
 * \param [in]  state  pointer to a State structure
 * \param [out] cmin  1-D array containing the leftmost characteristic
 * \param [out] cmin  1-D array containing the rightmost characteristic
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final   index of computation
 *
 *********************************************************************** */
{
  int     i;
  double  vx, vt2, vel2;
  double  sroot, *q;
  double  *cs2 = state->a2;

  for (i = beg; i <= end; i++) {

    q = state->v[i];
    
    vx   = q[VXn];
    vt2  = EXPAND(0.0, + q[VXt]*q[VXt], + q[VXb]*q[VXb]);
    vel2 = vx*vx + vt2;

    sroot = sqrt(cs2[i]*(1.0 - vx*vx - vt2*cs2[i])*(1.0 - vel2)); 
   
    cmax[i] = (vx*(1.0 - cs2[i]) + sroot)/(1.0 - vel2*cs2[i]);
    cmin[i] = (vx*(1.0 - cs2[i]) - sroot)/(1.0 - vel2*cs2[i]);

    g_maxMach = MAX(fabs(vx)/sqrt(cs2[i]), g_maxMach);
  }
}

/* ******************************************************************** */
void PrimEigenvectors (const State *state, int beg, int end)
/*!
 * Provide eigenvectors eigenvalues of the relativistic
 * equations in primitive form.
 *
 * \param [in,out]  state   Pointer to State structure
 * \param [in]      beg     starting grid index
 * \param [in]      end     final    grid index
 *
 *
 * \Note Only non-zero components of LL and RR are computed; 
 *       RR and LL must be initialized to zero outside.
 *
 *********************************************************************** */
{
  int    i, nv;
  int    kk,ii,jj, ll;
  double   *q, cs2, h, **LL, **RR, *lambda;
  double   cs, g, g2, g3, eta;
  double   D, g2_1d;
  double   a,b, r2, r3, r4;
  double   u, v, w, rp,rm, v2tan;
#if CHECK_EIGENVECTORS == YES  
  double Aw1[NFLX], Aw0[NFLX], AA[NFLX][NFLX];
#endif

  for (i = beg; i <= end; i++){

    q      = state->v[i];
    cs2    = state->a2[i];
    h      = state->h[i];
    LL     = state->Lp[i];
    RR     = state->Rp[i];
    lambda = state->lambda[i];

  /* -- Compute Lorentz factor -- */
     
  #if RECONSTRUCT_4VEL
    g  = EXPAND(q[VXn]*q[VXn], + q[VXt]*q[VXt], + q[VXb]*q[VXb]);
    g  = sqrt(1.0 + g);
    EXPAND(u = q[VXn]/g; ,
           v = q[VXt]/g; ,
           w = q[VXb]/g;)
  #else
    EXPAND(u = q[VXn]; ,
           v = q[VXt]; ,
           w = q[VXb];)
  #endif

  /* -- u,v,w  are now the three dimensional velocities -- */

    cs    = sqrt(cs2);
    v2tan = EXPAND(0.0, + v*v, + w*w);
    eta   = sqrt(1.0 - u*u - cs2*v2tan);
    g2_1d = 1.0/(1.0 - u*u);
    g2    = u*u + v2tan;
    a     = 1.0 - g2*cs2; 
    g2    = 1.0/(1.0 - g2);
    g     = sqrt(g2);
    D     = q[RHO]*g;                  /* Lab - Density */
    rm    = u*(1.0 - cs2) - cs*eta/g;
    rp    = u*(1.0 - cs2) + cs*eta/g;
  
  /* -- Define eigenvalues -- */
  
    lambda[0] = rm/a;
    lambda[1] = rp/a;
    for (kk = 2; kk < NFLX; kk++) lambda[kk] = u;
  
  /* --------------------------------------------------
      Equivalent to 
      
       rp = -cs*(g*eta*u + cs)/(g*D*(1.0 - u*u));
       rm =  cs*(g*eta*u - cs)/(g*D*(1.0 - u*u));
     -------------------------------------------------- */
  
    rp = cs*rp/(q[RHO]*(u*cs - g*eta));
    rm = cs*rm/(q[RHO]*(u*cs + g*eta)); 
  
    a = cs*eta/D;
    b = h*cs2;

  #if RECONSTRUCT_4VEL == NO

  /* -----------------------------------------
          RIGHT Eigenvectors, 3-vel
     ----------------------------------------- */

  /* -- lambda = u - cs -- */

    RR[RHO][0]  = 1.0;
    EXPAND(RR[VXn][0] = -a;           ,
           RR[VXt][0] = v*rm;  ,
           RR[VXb][0] = w*rm;)
    RR[PRS][0] = b;

  /* -- lambda = u + cs -- */

    RR[RHO][1] = 1.0;
    EXPAND(RR[VXn][1] = a;    ,
           RR[VXt][1] = v*rp; ,
           RR[VXb][1] = w*rp;)
    RR[PRS][1] = b;

  /* -- lambda = u -- */

    EXPAND(RR[RHO][2] = 1.0; ,
           RR[VXt][3] = 1.0; ,
           RR[VXb][4] = 1.0;)

  /* -------------------------------------------
            LEFT Eigenvectors, 3-vel
     ------------------------------------------- */

    a = 0.5/a;
    b = 0.5/b;

  /* -- lambda = u - c -- */

    LL[0][VXn] = -a;
    LL[0][PRS] =  b;

  /* -- lambda = u + c -- */

    LL[1][VXn] = a;
    LL[1][PRS] = b;

  /* -- lambda = u -- */

    LL[2][RHO] =  1.0;
    LL[2][PRS] = -2.0*b;

  /* -- lambda = u for transverse vel. -- */

    #if COMPONENTS > 1   
    LL[3][VXn] = u*v*g2_1d; 
    LL[3][VXt] = 1.0;                  
    LL[3][PRS] = v*g2_1d/(g2*q[RHO]*h);
    #endif
    #if COMPONENTS > 2 
    LL[4][VXn] = u*w*g2_1d; 
    LL[4][VXb] = 1.0;
    LL[4][PRS] = w*g2_1d/(g2*q[RHO]*h);
    #endif

  #elif RECONSTRUCT_4VEL == YES

  /* ---------------------------------------------------
          RIGHT Eigenvectors, 4-vel
     --------------------------------------------------- */

    g3 = g2*g;

  /* -- lambda = u - c -- */

    EXPAND(r2 = -(g + g3*u*u)*a + g3*u*rm*v2tan;  ,
           r3 = v*(-g3*u*a + rm*(g + g3*v2tan));  ,
           r4 = w*(-g3*u*a + rm*(g + g3*v2tan)); )

    RR[RHO][0] = 1.0;
    EXPAND(RR[VXn][0] = r2;  ,
           RR[VXt][0] = r3;  ,
           RR[VXb][0] = r4;)
    RR[ENG][0] = b;

  /* -- lambda = u + c -- */

    EXPAND(r2 =  (g + g3*u*u)*a + g3*u*rp*v2tan;  ,
           r3 = v*(g3*u*a + rp*(g + g3*v2tan));   ,
           r4 = w*(g3*u*a + rp*(g + g3*v2tan)); )

    RR[RHO][1] = 1.0;
    EXPAND(RR[VXn][1] = r2;  ,
           RR[VXt][1] = r3;  ,
           RR[VXb][1] = r4;)
    RR[ENG][1] = b;

  /* -- lambda = u -- */

    RR[RHO][2] = 1.0;

  /* -- lambda = u for transverse vel. -- */

    #if COMPONENTS > 1
    EXPAND(RR[VXn][3] = g3*u*v;      ,
           RR[VXt][3] = g + g3*v*v;  ,
           RR[VXb][3] = g3*v*w;)

    #endif
    #if COMPONENTS > 2
    RR[VXn][4] = g3*u*w;  
    RR[VXt][4] = g3*v*w;  
    RR[VXb][4] = g + g3*w*w;
    #endif

  /* --------------------------------------------------
          LEFT Eigenvectors, 4-vel
     -------------------------------------------------- */

  /* -- lambda = u - c -- */

    EXPAND(LL[0][VXn] =-0.5*(1.0 - u*u)/a/g;  ,
           LL[0][VXt] = 0.5*u*v/a/g;    ,
           LL[0][VXb] = 0.5*u*w/a/g;)
    LL[0][ENG] = 0.5/b;

  /* -- lambda = u + c -- */

    EXPAND(LL[1][VXn] =  0.5*(1.0 - u*u)/a/g; ,
           LL[1][VXt] = -0.5*u*v/a/g;      ,
           LL[1][VXb] = -0.5*u*w/a/g;)
    LL[1][ENG] = 0.5/b;

  /* -- lambda = u -- */

    LL[2][RHO] = 1.0;
    LL[2][ENG]  =-1.0/b;

    #if COMPONENTS > 1
    EXPAND(LL[3][VXn] = -0.5*v*(rp - rm)/a/g/g2_1d - u*v/g;  ,
           LL[3][VXt] =  0.5*v*v*(rp - rm)*u/a/g + (1.0 - v*v)/g;   ,
           LL[3][VXb] =  0.5*v*(rp - rm)*u*w/a/g - v*w/g;)
    LL[3][ENG] = -0.5*v*(rm + rp)/b;
    #endif
    #if COMPONENTS > 2
    LL[4][VXn] = -0.5*w*(rp - rm)/a/g/g2_1d - u*w/g;  
    LL[4][VXt] =  0.5*v*(rp - rm)*u*w/a/g - v*w/g;
    LL[4][VXb] =  0.5*w*w*(rp - rm)*u/a/g + (1.0 - w*w)/g;
    LL[4][ENG] = -0.5*w*(rp + rm)/b;
    #endif
  #endif  /* RECONSTRUCT_4VEL == YES */

/* -----------------------------------------
         Check eigenvectors consistency  
   ----------------------------------------- */

#if CHECK_EIGENVECTORS == YES

    for (ii = 0; ii <NFLX; ii++){
    for (jj = 0; jj <NFLX; jj++){
      AA[ii][jj] = 0.0;
      for (kk = 0; kk <NFLX; kk++){
      for (ll = 0; ll <NFLX; ll++){
        AA[ii][jj] += RR[ii][kk]*(kk==ll)*lambda[kk]*LL[ll][jj];
      }}
    }}
  
    PRIM_RHS (q, q, cs2, h, Aw0);
    for (ii = 0; ii <NFLX; ii++){
      Aw1[ii] = 0.0;
      for (jj = 0; jj <NFLX; jj++){
        Aw1[ii] += AA[ii][jj]*q[jj];
      }
    }
  
    a = 0.0; 
    for (ii = 0; ii <NFLX; ii++){
      a += fabs(Aw1[ii] - Aw0[ii]);
    }
  
    if (a > 1.e-1) {
      print ("! Eigenvectors not consistent in EIGENV%12.6e\n",a);
      for (ii = 0; ii < NFLX; ii++){
        for (jj = 0; jj < NFLX; jj++){
          print ("%12.6e   ",AA[ii][jj]);
        }
        print("\n");
      }
  
  
      print ("Aw0:  ");
      for (ii = 0; ii < NFLX; ii++){
        print ("%12.6e ,  ",Aw0[ii]);
      }
      print ("\nAw1:  ");
      for (ii = 0; ii < NFLX; ii++){
        print ("%12.6e ,  ",Aw1[ii]);
      }
      QUIT_PLUTO(1);
    }
  
  /*  Check eigenvectors orthonormality   */
  
    for (ii = 0; ii <NFLX;ii++){
      for (jj = 0; jj <NFLX;jj++){
        a = 0.0;
        for (kk = 0; kk <NFLX;kk++){
          a += LL[kk][ii]*RR[jj][kk];
        }
        if (ii==jj && fabs(a-1.0)>1.e-8) {
          print ("! Eigenvectors not orthogonal!  %d  %d  %12.6e \n",ii,jj,a);
          print ("NSweep: %d\n",g_dir);
          QUIT_PLUTO(1);
        }
        if(ii!=jj && fabs(a)>1.e-8) {
          print ("! Eigenvectors not orthogonal (2) %d  %d  %12.6e !\n",ii,jj,a);
          print ("NSweep: %d\n",g_dir);
          QUIT_PLUTO(1);
        }
      } 
    }
#endif
  } /* End loop on i = beg, end */
}

/* ***************************************************************** */
void PrimToChar (double **L, double *v, double *w)
/*
 *
 *  Compute the matrix-vector product between the
 *  the L matrix (containing primitive left eigenvectors) 
 *  and the vector v. 
 *  For efficiency purpose, multiplication is done 
 *  explicitly, so that only nonzero entries
 *  of the left primitive eigenvectors are considered.
 *  
 *
 ****************************************************************** */
{
  int nv;

  #if RECONSTRUCT_4VEL == NO

   w[0] = L[0][VXn]*v[VXn] + L[0][PRS]*v[PRS];
   w[1] = L[1][VXn]*v[VXn] + L[1][PRS]*v[PRS];
   EXPAND( w[2] = v[RHO] + L[2][PRS]*v[PRS];                    ,
           w[3] = L[3][VXn]*v[VXn] + v[VXt] + L[3][PRS]*v[PRS];   ,
           w[4] = L[4][VXn]*v[VXn] + v[VXb] + L[4][PRS]*v[PRS];)

  #elif RECONSTRUCT_4VEL == YES
 
   w[0] = EXPAND(  L[0][VXn]*v[VXn],
                 + L[0][VXt]*v[VXt],
                 + L[0][VXb]*v[VXb]) + L[0][PRS]*v[PRS];
   w[1] = EXPAND(  L[1][VXn]*v[VXn],
                 + L[1][VXt]*v[VXt],
                 + L[1][VXb]*v[VXb]) + L[1][PRS]*v[PRS];;

   w[2] = v[RHO] + L[2][PRS]*v[PRS];
   
   #if COMPONENTS > 1
    w[3] = EXPAND(  L[3][VXn]*v[VXn],
                  + L[3][VXt]*v[VXt],
                  + L[3][VXb]*v[VXb]) +  L[3][PRS]*v[PRS];

    #if COMPONENTS == 3
    w[4] = L[4][VXn]*v[VXn],
           L[4][VXt]*v[VXt],
           L[4][VXb]*v[VXb] +  L[4][PRS]*v[PRS];
    #endif
   #endif

  #endif

/* ----------------------------------------------- 
     For passive scalars, the characteristic 
     variable is equal to the primitive one, 
     since  l = r = (0,..., 1 , 0 ,....)
    ----------------------------------------------- */		   

#if NSCL > 0 
  NSCL_LOOP(nv) w[nv] = v[nv];
#endif   
}
