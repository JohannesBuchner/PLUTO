#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double x, y, R, Rc, phi, c, s;
  double arg, kp, h;
  double dvx, dvy;

  kp  = g_inputParam[KPAR];
  h   = g_inputParam[HSCALE];  /* -- the scale of the vortex -- */
  Rc  = 1.0;          /* -- the radial position of the vortex -- */

/* --------------------------------------------
              set coordinates 
   -------------------------------------------- */

  R   = x1; 
  phi = x2;
  c = cos(phi); 
  s = sin(phi);

  x = R*c - Rc/sqrt(2.0);
  y = R*s - Rc/sqrt(2.0);

  us[RHO] = 1.0;
  #if EOS == IDEAL
   us[PRS] = 1.0/(g_inputParam[MACH]*g_inputParam[MACH]*g_gamma);
  #else
   g_isoSoundSpeed = 1.0/g_inputParam[MACH];
  #endif
  us[VX1] = 0.0;
  us[VX2] = 1.0/sqrt(R);

  arg = - (x*x + y*y)/(h*h);

  dvx = -y*kp*exp(arg);
  dvy =  x*kp*exp(arg);

  us[VX1] +=  dvx*c + dvy*s;
  us[VX2] += -dvx*s + dvy*c;

  #if ROTATING_FRAME
   g_OmegaZ  = 1.0;
   us[VX2]  -= g_OmegaZ*R;
  #endif

  #if NTRACER == 1
   us[TRC] = (exp(arg) > 1.e-2);
  #endif

}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
  int    i,j,k,nv;
  static int first_call = 1;
  double v[NVAR];
  double *R, *phi, *dVR, *dVphi, *dR, *dphi, vol;
  double **vR, **vphi;
  double E=0.0, w, pw, wmin, pwmin;
  double Ltot, Etot, wtot, pwtot;
  static double voltot;
  FILE *fp;

  R   = grid[IDIR].x;  phi   = grid[JDIR].x;
  dVR = grid[IDIR].dV; dVphi = grid[JDIR].dV;
  dR  = grid[IDIR].dx; dphi  = grid[JDIR].dx;

/* ---------------------------------------------------------
         compute total volume once at the beginning
   --------------------------------------------------------- */

  if (first_call){
    voltot = 0.0;
    DOM_LOOP(k,j,i) voltot += dVR[i]*dVphi[j];
    #ifdef PARALLEL
     MPI_Allreduce (&voltot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Barrier (MPI_COMM_WORLD);
     voltot = vol;
    #endif
  }

/* --------------------------------------------------------
              compute volume averages 
   -------------------------------------------------------- */

  vR = d->Vc[VX1][0]; vphi = d->Vc[VX2][0];

  Etot = Ltot = wtot = pwtot = 0.0;
  wmin = pwmin = 1.e20;
  DOM_LOOP(k,j,i){
    vol = dVR[i]*dVphi[j]/voltot;

    for (nv = NVAR; nv--;   ) v[nv] = d->Vc[nv][k][j][i];

    w  =   (  (vphi[j][i+1] - 1.0/sqrt(R[i+1]))*R[i+1] 
            - (vphi[j][i-1] - 1.0/sqrt(R[i-1]))*R[i-1])/(2.0*R[i]*dR[i]);
    w -=  (vR[j+1][i] - vR[j-1][i])/(2.0*R[i]*dphi[j]);

    pw = w/v[RHO];
    #if EOS == IDEAL
    E  = v[PRS]/(g_gamma-1.0) 
         + 0.5*v[RHO]*(v[VX1]*v[VX1] + v[VX2]*v[VX2]) - v[RHO]/R[i];
     #if PHYSICS == MHD
      E += 0.5*(v[BX1]*v[BX1] + v[BX2]*v[BX2]);
     #endif
    #endif

    Ltot  += R[i]*v[RHO]*v[VX2]*vol;
    Etot  += E*vol;
    wtot  += w*vol;
    pwtot += pw*vol;

    wmin  = MIN(wmin, w);
    pwmin = MIN(pwmin, pw);
  }

  #ifdef PARALLEL
   MPI_Allreduce (&Ltot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Ltot = vol;

   MPI_Allreduce (&Etot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Etot = vol;

   MPI_Allreduce (&wtot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   wtot = vol;

   MPI_Allreduce (&pwtot, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   pwtot = vol;

   MPI_Allreduce (&wmin, &vol, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   wmin = vol;

   MPI_Allreduce (&pwmin, &vol, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   pwmin = vol;

   MPI_Barrier (MPI_COMM_WORLD);
  #endif
/*
printf ("<pm> = %12.6e, <pr> = %12.6e, <pr>/<pm> = %12.6e, <pr*f>/<pm*f> = %12.6e\n",
         pmtot, prtot, prtot/pmtot, beta_num/beta_den);
exit(1);
*/
/* -----------------------------------------------------
         processor #0 does the actual writing 
   ----------------------------------------------------- */

  if (prank == 0){
    static double tpos = -1.0;
    if (g_stepNumber == 0){
      fp = fopen("averages.dat","w");
      fprintf (fp,"#%s %12s  %13s  %13s  %13s  %14s  %12s  %12s\n", 
               "      t","dt","<L>","<E>","<w>","<w/rho>","min(w)","min(pw)");
    }else{
      if (tpos < 0.0){  /* obtain time coordinate of to last entry */
        char   sline[512];
        fp = fopen("averages.dat","r");
        while (fgets(sline, 512, fp))  {}
        sscanf(sline, "%lf\n",&tpos);
        fclose(fp);
      }
      fp = fopen("averages.dat","a");
    }
    if (g_time > tpos){
      fprintf (fp, "%12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n", 
               g_time, g_dt, Ltot, Etot, wtot, pwtot, wmin, pwmin);
    }
    fclose(fp);
  }
  first_call = 0; 
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int    i, j, k, nv;
  double *R, v[256];

  R = grid[IDIR].x;

  if (side == X1_BEG){

    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
        #if EOS == IDEAL
         d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
        #endif
        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IBEG];
        d->Vc[VX2][k][j][i] =  1.0/sqrt(R[i]) + d->Vc[VX2][k][j][IBEG] 
                            - 1.0/sqrt(R[IBEG]);
        #if ROTATING_FRAME == YES
         d->Vc[VX2][k][j][i] += g_OmegaZ*(R[i] - R[IBEG]);
        #endif
        #if NTRACER == 1
         d->Vc[TRC][k][j][i] = 0.0;
        #endif
      }
    }

  } else if (side == X1_END) {

    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];
        #if EOS == IDEAL
         d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND];
        #endif
        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND];
        d->Vc[VX2][k][j][i] =   1.0/sqrt(R[i]) + d->Vc[VX2][k][j][IEND] 
                             - 1.0/sqrt(R[IEND]);
        #if ROTATING_FRAME == YES
         d->Vc[VX2][k][j][i] += g_OmegaZ*(R[i] - R[IEND]);
        #endif
      }
    }
  } 
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = -1.0/(x1*x1);
  g[JDIR] = 0.0; 
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  return 0.0;
}
#endif

