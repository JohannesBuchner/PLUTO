#include "pluto.h"

static double alpha, pe;

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double vj, cs, Mach;
  double lor, vy_in, rhob, rhom, pm, pb, hb, hm, eta, vjet;
  double sigma;
  static int first_call = 1;

  g_gamma = 4.0/3.0;

  alpha = 1.0 - g_inputParam[RM]*g_inputParam[RM]/g_inputParam[BETA];
  pe    = 0.5*g_inputParam[BETA]*g_inputParam[BM]*g_inputParam[BM];

  us[RHO] = 1.e3;
  us[VX1] = 0.0;
  us[VX2] = 0.0;
  us[VX3] = 0.0;
  us[BX1] = 0.0;
  us[BX2] = 0.0;
  us[BX3] = 0.0;
  us[PRS] = pe;

  us[AX1] = us[AX2] = us[AX3] = 0.0;

  if (first_call){
    first_call = 0;
    vj   = sqrt(1.0 - 1.0/100.);
    cs   = sqrt(g_gamma*pe/(us[RHO] + pe*g_gamma/(g_gamma-1.0)));
    Mach = vj/cs;
    print ("Mach number = %12.6e\n",Mach);

/* ----- estimate jet velocity --------- */
    
    lor   = 10.0;
    vy_in = sqrt(1.0 - 1.0/lor/lor);
    rhob = g_inputParam[RHO_IN];
    rhom = 1.e3;
    pm   = pb = pe;
    hb   = 1.0 + g_gamma/(g_gamma - 1.0)*pb/rhob;
    hm   = 1.0 + g_gamma/(g_gamma - 1.0)*pm/rhom;

    eta  = rhob*hb/(rhom*hm)*lor*lor;
    vjet = sqrt(eta)/(1.0 + sqrt(eta))*vy_in;
    print (" Estimated jet velocity = %f\n",vjet);

/* ---- estimate sigma (= b^2/rho) ---- */

    sigma = g_inputParam[RM]*g_inputParam[RM]*(0.5 - 2.0*log(g_inputParam[RM]));
    print (" sigma = %f\n",1.0/sigma);

  }
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int     i, j, k, nv;
  double  rho_out, vx_out, vy_out, vz_out;
  double  pr_out,  bx_out, by_out, bz_out;
  double  prof, lor, pr_in, vy_in, bz_in, *r;


  #ifdef STAGGERED_MHD
   #error Boundary not implemented for Staggered MHD
  #endif

  if (side == X2_BEG){
    if (box->vpos == CENTER){
      r = grid->x[IDIR];
      BOX_LOOP(box,k,j,i){

        rho_out = d->Vc[RHO][k][2*JBEG - j - 1][i];            
        EXPAND(vx_out =  d->Vc[VX1][k][2*JBEG - j - 1][i];  ,
               vy_out = -d->Vc[VX2][k][2*JBEG - j - 1][i];  ,
               vz_out =  d->Vc[VX3][k][2*JBEG - j - 1][i];)
        EXPAND(bx_out =  d->Vc[BX1][k][2*JBEG - j - 1][i];  ,
               by_out =  d->Vc[BX2][k][2*JBEG - j - 1][i];  ,
               bz_out = -d->Vc[BX3][k][2*JBEG - j - 1][i];)             
        pr_out  = d->Vc[PRS][k][2*JBEG - j - 1][i];

        prof = (r[i] <= 1.0 ? 1.0 : 0.0);

        lor   = 10.0;
        vy_in = sqrt(1.0 - 1.0/lor/lor);

        if (r[i] <= 1.0) {
          prof = 1.0;
          bz_in = lor*g_inputParam[BM]*g_inputParam[RM]/r[i];
          pr_in = alpha*pe;         
          if (r[i] <= g_inputParam[RM]) {
            bz_in = lor*g_inputParam[BM]*r[i]/g_inputParam[RM];
            pr_in = pe*(alpha + 2.0/g_inputParam[BETA]*(1.0 - r[i]*r[i]/g_inputParam[RM]/g_inputParam[RM]));
          }
        }else{
          prof = 0.0;
        }  

        d->Vc[RHO][k][j][i] = rho_out - (rho_out - g_inputParam[RHO_IN])*prof;
        EXPAND(d->Vc[VX1][k][j][i] = vx_out*(1.0 - prof);              ,
               d->Vc[VX2][k][j][i] = vy_out - (vy_out - vy_in)*prof;   ,
               d->Vc[VX3][k][j][i] = vz_out*(1.0 - prof);)  

        EXPAND(d->Vc[BX1][k][j][i] = bx_out*(1.0 - prof);            ,
               d->Vc[BX2][k][j][i] = by_out*(1.0 - prof);            ,
               d->Vc[BX3][k][j][i] = bz_out - (bz_out - bz_in)*prof; )
        d->Vc[PRS][k][j][i] = pr_out - (pr_out - pr_in)*prof;
      }
    }
  }
}

