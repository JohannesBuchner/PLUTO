#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

/* ************************************************************************* */
//void PatchPluto::computeRefGradient(FArrayBox& gFab, FArrayBox& UFab, const Box& b)
void PatchPluto::computeRefGradient(FArrayBox& gFab, FArrayBox& UFab, 
                                    const FArrayBox& a_dV, const Box& b)
/*
 *
 * PURPOSE
 *
 *   Tag cells for refinement by computing grad[k][j][i]. 
 *   By default a convex combination of the first and second
 *   derivative of the total energy density is used.
 *   alpha = 0 --> triggers refinement towards the 2nd derivative
 *   alpha = 1 --> triggers refinement towards the 1st derivative
 *
 * 
 *
 *************************************************************************** */
{
  CH_assert(m_isDefined);

  int nv, i, j, k;
  int Uib, Uie, Ujb=0, Uje=0, Ukb=0, Uke=0;
  int Gib, Gie, Gjb=0, Gje=0, Gkb=0, Gke=0;

  double rp, rm, r;
  double x, dqx_p, dqx_m, d2qx, den_x;
  double y, dqy_p, dqy_m, d2qy, den_y;
  double z, dqz_p, dqz_m, d2qz, den_z;

  double alpha, qref, gr1, gr2;

  double eps = 0.01;
  double ***UU[NVAR], ***q, ***grad;

  double us[NVAR], vs[NVAR], mu;
  static double **T;

  rp = rm = r = 1.0;

/* -----------------------------------------------
    The solution array U is defined on the box 
    [Uib, Uie] x [Ujb, Uje] x [Ukb, Uke], which 
    differs from that of gFab ([Gib,...Gke]), 
    typically one point larger in each direction. 
   ----------------------------------------------- */
    
  D_EXPAND(Uib = UFab.loVect()[IDIR]; Uie = UFab.hiVect()[IDIR]; ,
           Ujb = UFab.loVect()[JDIR]; Uje = UFab.hiVect()[JDIR]; ,
           Ukb = UFab.loVect()[KDIR]; Uke = UFab.hiVect()[KDIR]; );

  D_EXPAND(Gib = gFab.loVect()[IDIR]; Gie = gFab.hiVect()[IDIR]; ,
           Gjb = gFab.loVect()[JDIR]; Gje = gFab.hiVect()[JDIR]; ,
           Gkb = gFab.loVect()[KDIR]; Gke = gFab.hiVect()[KDIR]; );

  for (nv=0 ; nv<NVAR ; nv ++){
    UU[nv] = ArrayBoxMap(Ukb, Uke, Ujb, Uje, Uib, Uie, UFab.dataPtr(nv));
  }
  grad = ArrayBoxMap(Gkb, Gke, Gjb, Gje, Gib, Gie,gFab.dataPtr(0));

/* -----------------------------------------------
    the parameter alpha controls the bias towards
    1st derivative criterion (alpha = 1) or 2nd
    derivative (alpha = 0)
   ----------------------------------------------- */

  alpha = 0.0;
  #if (EOS != ISOTHERMAL) && (AMR_EN_SWITCH == NO)
   q  = UU[ENG];   
  #else
   q  = UU[RHO];
  #endif

  for (k = Gkb; k <= Gke; k++) { z = k*m_dx + g_domBeg[KDIR];
  for (j = Gjb; j <= Gje; j++) { y = j*m_dx + g_domBeg[JDIR];
  for (i = Gib; i <= Gie; i++) { x = i*m_dx + g_domBeg[IDIR];

    #if GEOMETRY == CYLINDRICAL
     rp = (i+0.5)/(i+1.5);
     rm = (i+0.5)/(i-0.5);
    #endif

    D_EXPAND(dqx_p =   q[k][j][i+1]*rp - q[k][j][i];
             dqx_m = -(q[k][j][i-1]*rm - q[k][j][i]);  ,
             dqy_p =   q[k][j+1][i] - q[k][j][i];
             dqy_m = -(q[k][j-1][i] - q[k][j][i]);     ,
             dqz_p =   q[k+1][j][i] - q[k][j][i];
             dqz_m = -(q[k-1][j][i] - q[k][j][i]);)

  /* --------------------------------------------------------------
      Physical boundary values are not up to date and should be 
      excluded from gradient computation. 
      In this case, left and right derivatives are set equal to 
      each other. This will not trigger refinement in the leftmost 
      and rightmost internal zones (using 2nd derivative) but we 
      really don't care since buffer size will do the job.
     -------------------------------------------------------------- */
      
    D_EXPAND(if (i == 0) dqx_m = dqx_p;  ,
             if (j == 0) dqy_m = dqy_p;  ,
             if (k == 0) dqz_m = dqz_p;)

    D_EXPAND(if (i == m_domain.size(IDIR)-1) dqx_p = dqx_m;  ,
             if (j == m_domain.size(JDIR)-1) dqy_p = dqy_m;  ,
             if (k == m_domain.size(KDIR)-1) dqz_p = dqz_m;)

  /* ---------------------------------------------------
                    New version
     --------------------------------------------------- */
/*
    D_EXPAND(d2qx = dqx_p - dqx_m;  ,
             d2qy = dqy_p - dqy_m;  ,
             d2qz = dqz_p - dqz_m;)

    D_EXPAND(
      den_x = 2.0*fabs(q[k][j][i]) + fabs(q[k][j][i+1]) + fabs(q[k][j][i-1]);
      den_x = fabs(dqx_p) + fabs(dqx_m) + eps*den_x;    ,

      den_y = 2.0*fabs(q[k][j][i]) + fabs(q[k][j+1][i]) + fabs(q[k][j-1][i]);
      den_y = fabs(dqy_p) + fabs(dqy_m) + eps*den_y;    ,

      den_z = 2.0*fabs(q[k][j][i]) + fabs(q[k+1][j][i]) + fabs(q[k-1][j][i]);
      den_z = fabs(dqz_p) + fabs(dqz_m) + eps*den_z;
    )

    gr2  = D_EXPAND(d2qx*d2qx, + d2qy*d2qy, + d2qz*d2qz);
    gr2 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);
    grad[k][j][i] = sqrt(gr2);

  /* ---------------------------------------------------
                      Old version
     --------------------------------------------------- */

 /* -- first derivative -- */

    eps   = 0.05;
    qref = fabs(q[k][j][i]); 
    gr1 = D_EXPAND( fabs(dqx_p + dqx_m)/qref ,
                  + fabs(dqy_p + dqy_m)/qref ,
                  + fabs(dqz_p + dqz_m)/qref);

 /* -- second derivative -- */

    gr2 = D_EXPAND( fabs(dqx_p - dqx_m)/(fabs(dqx_p) + fabs(dqx_m) + eps*qref) ,
                  + fabs(dqy_p - dqy_m)/(fabs(dqy_p) + fabs(dqy_m) + eps*qref) ,
                  + fabs(dqz_p - dqz_m)/(fabs(dqz_p) + fabs(dqz_m) + eps*qref));

    grad[k][j][i] = alpha*gr1 + (1.0 - alpha)*gr2;

  }}}

  for (nv=0 ; nv<NVAR ; nv ++){
    FreeArrayBoxMap(UU[nv], Ukb, Uke, Ujb, Uje, Uib, Uie);
  }
  FreeArrayBoxMap(grad, Gkb, Gke, Gjb, Gje, Gib, Gie);
}

