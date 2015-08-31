#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

static void computeRefVar(double ***UU[], double ***q, double, RBox *Ubox);

#if (EOS != ISOTHERMAL) && (ENTROPY_SWITCH == NO)
 #ifndef CHOMBO_REF_VAR  
  #define CHOMBO_REF_VAR ENG 
 #endif
#else
 #ifndef CHOMBO_REF_VAR  
  #define CHOMBO_REF_VAR RHO
 #endif
#endif

#define REF_CRIT 2   /* 1 == first derivative, 2 == second derivative */

/* ************************************************************************* */
void PatchPluto::computeRefGradient(FArrayBox& gFab, FArrayBox& UFab, 
                                    const FArrayBox& a_dV, const Box& b)
/*!
 * Tag zones for refinement using gradient of the conservative 
 * variables.
 * The gradient is computed by standard finite differences using
 *
 * - REF_CRIT equal to 1 --> compute (normalized) gradient using 1st 
 *                           derivative of the solution;
 * - REF_CRIT equal to 2 --> compute (normalized) gradient using 2nd 
 *                           derivative of the solution (default);
 *                           This approach is based on Lohner (1987).
 *
 * Zones will be flagged for refinement whenever grad[k][j][i] exceeds 
 * the threshold value specified by the 'Refine_thresh' parameter read in
 * pluto.ini.
 *
 * Derivatives are computed using the conserved variable
 * U[CHOMBO_REF_VAR] 
 * where CHOMBO_REF_VAR is taken to be energy density (default).
 * However, by setting CHOMBO_REF_VAR = -1, you can provide your own 
 * physical variable through the function computeRefVar().
 * 
 * \authors C. Zanni   (zanni@oato.inaf.it)\n
 *          A. Mignone (mignone@ph.unito.it)
 * \date    Oct 11, 2012
 *************************************************************************** */
{
  CH_assert(m_isDefined);

  int nv, i, j, k;
  double x1, dqx_p, dqx_m, dqx, d2qx, den_x;
  double x2, dqy_p, dqy_m, dqy, d2qy, den_y;
  double x3, dqz_p, dqz_m, dqz, d2qz, den_z;
  double gr1, gr2, eps = 0.01;
#if CHOMBO_REF_VAR == -1
  double ***UU[NVAR];
#endif
  double ***q, ***grad;
  RBox  Ubox, Gbox;

/* -- check ref criterion -- */

#if REF_CRIT != 1 && REF_CRIT != 2
  print ("! TagCells.cpp: Refinement criterion not valid\n");
  QUIT_PLUTO(1);
#endif

/* -----------------------------------------------------
   1. The solution array U is defined on the box 
      [Uib, Uie] x [Ujb, Uje] x [Ukb, Uke], which 
      differs from that of gFab ([Gib,...Gke]), 
      typically one point larger in each direction. 
   ----------------------------------------------------- */
    
  Ubox.jb = Ubox.je = Ubox.kb = Ubox.ke = 0;
  Gbox.jb = Gbox.je = Gbox.kb = Gbox.ke = 0;

  D_EXPAND(Ubox.ib = UFab.loVect()[IDIR]; Ubox.ie = UFab.hiVect()[IDIR]; ,
           Ubox.jb = UFab.loVect()[JDIR]; Ubox.je = UFab.hiVect()[JDIR]; ,
           Ubox.kb = UFab.loVect()[KDIR]; Ubox.ke = UFab.hiVect()[KDIR]; );

  D_EXPAND(Gbox.ib = gFab.loVect()[IDIR]; Gbox.ie = gFab.hiVect()[IDIR]; ,
           Gbox.jb = gFab.loVect()[JDIR]; Gbox.je = gFab.hiVect()[JDIR]; ,
           Gbox.kb = gFab.loVect()[KDIR]; Gbox.ke = gFab.hiVect()[KDIR]; );

/* --------------------------------------------------------
   2. Input solution array (UFab.dataPtr(nv)) is defined 
      as dV*U/dx^3, where U is an array of conservative 
      variables and dV is the zone volume. 
      To obtain U we must divide by volume.
   -------------------------------------------------------- */

#if CHOMBO_REF_VAR == -1
  FArrayBox tmpU(UFab.box(),NVAR);
  tmpU.copy(UFab);
#else
  FArrayBox tmpU(UFab.box(),1);
  tmpU.copy(UFab,CHOMBO_REF_VAR,0);
#endif 

#if GEOMETRY != CARTESIAN

  #if CHOMBO_REF_VAR == -1

    for (nv = 0; nv < NVAR; nv++) tmpU.divide(a_dV,0,nv);
    #if CHOMBO_CONS_AM == YES
      #if ROTATING_FRAME == YES
        Box curBox = UFab.box();
        for(BoxIterator bit(curBox); bit.ok(); ++bit) {
          const IntVect& iv = bit();
          tmpU(iv,iMPHI) /= a_dV(iv,1);
          tmpU(iv,iMPHI) -= tmpU(iv,RHO)*a_dV(iv,1)*g_OmegaZ;
        }
      #else
        tmpU.divide(a_dV,1,iMPHI);
      #endif
    #endif

  #else

    tmpU.divide(a_dV,0,0);

  #endif

#else

   if (g_stretch_fact != 1.) tmpU /= g_stretch_fact;

#endif // GEOMETRY == CARTESIAN

/* ---------------------------------------------
   3. Set refinement variable
   --------------------------------------------- */

#if CHOMBO_REF_VAR == -1
  for (nv = 0; nv < NVAR; nv++)
    UU[nv] = ArrayBoxMap(Ubox.kb, Ubox.ke,
                         Ubox.jb, Ubox.je,
                         Ubox.ib, Ubox.ie, tmpU.dataPtr(nv));
  q = ArrayBox(Ubox.kb, Ubox.ke, Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
  computeRefVar(UU, q, m_dx, &Ubox);
#else
  q = ArrayBoxMap(Ubox.kb, Ubox.ke,
                  Ubox.jb, Ubox.je,
                  Ubox.ib, Ubox.ie, tmpU.dataPtr(0));
#endif
  
  grad = ArrayBoxMap(Gbox.kb, Gbox.ke, 
                     Gbox.jb, Gbox.je, 
                     Gbox.ib, Gbox.ie, gFab.dataPtr(0));

/* ----------------------------------------------------------------
   4. Main spatial loop for zone tagging based on 1st 
     (REF_CRIT = 1) or 2nd (REF_CRIT = 2) derivative error norm. 
   ---------------------------------------------------------------- */

  BOX_LOOP(&Gbox, k, j, i){
    x3 = (k + 0.5)*m_dx*g_x3stretch + g_domBeg[KDIR];
    x2 = (j + 0.5)*m_dx*g_x2stretch + g_domBeg[JDIR];
#if CHOMBO_LOGR == NO
    x1 = (i + 0.5)*m_dx          + g_domBeg[IDIR];
#else
    double xl = g_domBeg[IDIR] + i*m_dx;
    double xr = xl + m_dx;
    x1 = g_domBeg[IDIR]*0.5*(exp(xr)+exp(xl));
#endif 

    D_EXPAND(dqx_p =    q[k][j][i+1] - q[k][j][i];
             dqx_m = - (q[k][j][i-1] - q[k][j][i]);  ,
             dqy_p =    q[k][j+1][i] - q[k][j][i];
             dqy_m = - (q[k][j-1][i] - q[k][j][i]);  ,
             dqz_p =    q[k+1][j][i] - q[k][j][i];
             dqz_m = - (q[k-1][j][i] - q[k][j][i]);)

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

  /* -----------------------------------------------
         Compute gradient using 1st derivative 
      ---------------------------------------------- */

#if REF_CRIT == 1
    D_EXPAND(dqx = dqx_p + dqx_m;  ,
             dqy = dqy_p + dqy_m;  ,
             dqz = dqz_p + dqz_m;)

    D_EXPAND(den_x = fabs(q[k][j][i+1]) + fabs(q[k][j][i-1]);  ,
             den_y = fabs(q[k][j+1][i]) + fabs(q[k][j-1][i]);  ,
             den_z = fabs(q[k+1][j][i]) + fabs(q[k-1][j][i]);)

    gr1  = D_EXPAND(dqx*dqx, + dqy*dqy, + dqz*dqz);
    gr1 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);

    grad[k][j][i] = sqrt(gr1);
#endif

  /* -----------------------------------------------
         Compute gradient using 2nd derivative 
      ---------------------------------------------- */

#if REF_CRIT == 2
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

    gr2  = D_EXPAND(d2qx*d2qx,   + d2qy*d2qy,   + d2qz*d2qz);
    gr2 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);

    grad[k][j][i] = sqrt(gr2);
#endif
  }

/* --------------------------------------------------------------
   6. Free array
   -------------------------------------------------------------- */
   
  FreeArrayBoxMap(grad, Gbox.kb, Gbox.ke, Gbox.jb, Gbox.je, Gbox.ib, Gbox.ie);

#if CHOMBO_REF_VAR == -1
  for (nv = 0; nv < NVAR; nv++){
    FreeArrayBoxMap(UU[nv], Ubox.kb, Ubox.ke,
                            Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
  }
  FreeArrayBox(q, Ubox.kb, Ubox.jb, Ubox.ib);
#else
  FreeArrayBoxMap(q, Ubox.kb, Ubox.ke, Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
#endif

}

/* ********************************************************************* */
void computeRefVar(double ***UU[], double ***q, double dx, RBox *Ubox)
/*!
 * Compute a user-defined array q(U) function of the conserved
 * variables.
 *
 *
 *********************************************************************** */
{
  int nv, i, j, k;
  double us[NVAR], vs[NVAR];

  BOX_LOOP(Ubox, k, j, i) {
    VAR_LOOP(nv) us[nv] = UU[nv][k][j][i]; 
    q[k][j][i] = us[RHO];
  }

}
#undef REF_CRIT
