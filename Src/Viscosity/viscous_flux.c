/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Computes viscous fluxes and source terms for the HD/MHD equations. 

  \authors Petros Tzeferacos (petros.tzeferacos@ph.unito.it) \n
           Andrea Mignone
  \date    Sep 18, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ViscousFlux (Data_Arr V, double **ViF, double **ViS, 
                  double **dcoeff, int beg, int end, Grid *grid)
/*!
 *
 *  \param [in]      V  data array containing cell-centered quantities
 *  \param [in,out]  ViF pointer to viscous fluxes
 *  \param [in,out]  ViS pointer to viscous source terms
 *  \param [in,out]  dcoeff  pointer to diffusion coefficient for dt calculation
 *  \param [in]      beg     integer, index for loop beg
 *  \param [in]      end     integer, index for loop end
 *  \param [in]      grid  pointer to array of Grid structures 
 *
 *  \return This function has no return value.
 *
 *  \notes Calculates the stress tensor components (at the i+1/2 face of each cell) and
 *   adds explicit viscous terms to the energy and momentum equation. It is 
 *   called in the during the sweep integrators. The stress tensor is given by
 *
 *  \f[                                                      
 *     T = \left(
 *     \begin{array}{ccc}
 *                     T_{xx} & T_{xy} & T_{xz} \\
 *                     T_{yx} & T_{yy} & T_{yz} \\
 *                     T_{zx} & T_{zy} & T_{zz} 
 *     \end{array}\right)
 * \f]
 *    
 *  where \f$ T_{ij} = T_{ji}\f$ and the components are given by 
 *  \f$ T_{ij} = 2\,m\,e(ij) + (l - 2/3 m) \nabla\cdot {\bf V}\, \delta_{ij} \f$
 *  where \f$ e(ij)= e_{ij}/(h_ih_j)\f$ and  
 *  \f$ e_{ij} = 0.5( V_{i;j} + V_{j;i} ) \f$  whereas m,l are the first and 
 *  second parameter of viscosity respectively.
 *
 ************************************************************************* */
#define D_DX_I(q)  (q[k][j][i + 1] - q[k][j][i])
#define D_DY_J(q)  (q[k][j + 1][i] - q[k][j][i])
#define D_DZ_K(q)  (q[k + 1][j][i] - q[k][j][i])

#define D_DY_I(q)  (  0.25*(q[k][j + 1][i] + q[k][j + 1][i + 1]) \
                    - 0.25*(q[k][j - 1][i] + q[k][j - 1][i + 1]))

#define D_DZ_I(q)  (  0.25*(q[k + 1][j][i] + q[k + 1][j][i + 1])  \
                    - 0.25*(q[k - 1][j][i] + q[k - 1][j][i + 1]))

#define D_DX_J(q)  (  0.25*(q[k][j][i + 1] + q[k][j + 1][i + 1]) \
                    - 0.25*(q[k][j][i - 1] + q[k][j + 1][i - 1]))

#define D_DZ_J(q)  (  0.25*(q[k + 1][j][i] + q[k + 1][j + 1][i]) \
                    - 0.25*(q[k - 1][j][i] + q[k - 1][j + 1][i]))

#define D_DX_K(q)  (  0.25*(q[k][j][i + 1] + q[k + 1][j][i + 1]) \
                    - 0.25*(q[k][j][i - 1] + q[k + 1][j][i - 1]))

#define D_DY_K(q)  (  0.25*(q[k][j + 1][i] + q[k + 1][j + 1][i]) \
                    - 0.25*(q[k][j - 1][i] + q[k + 1][j - 1][i]))                
{
  int i,j,k,n,nv;
  double nu1,nu2;
  double div_v;
  double dx_1, dy_1, dz_1;
  double dx1, dx2, dx3;
  double *x1, *x2, *x3;
  double *x1r, *x2r, *x3r;
  double dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz;
  static double *tau_xx, *tau_xy, *tau_xz,
                *tau_yx, *tau_yy, *tau_yz,
                *tau_zx, *tau_zy, *tau_zz; 
  double ***Vx, ***Vy, ***Vz;
  double scrh;
  double *r, *th, r_1, dr, dr_1, s_1, tan_1, dy;
  double dVxi,dVyi,dVzi;
  double dVxj,dVyj,dVzj;
  double dVxk,dVyk,dVzk;
  double vc[NVAR], vi[NVAR]; /* Center and interface values */
  static double *one_dVr, *one_dmu; /*auxillary volume components for r_1 singularity @ cylindrical and spherical*/
  
  EXPAND(Vx = V[VX1];  ,
         Vy = V[VX2];  ,
         Vz = V[VX3];)
  
  if (tau_xx == NULL){
    tau_xx = ARRAY_1D(NMAX_POINT, double);
    tau_xy = ARRAY_1D(NMAX_POINT, double);
    tau_xz = ARRAY_1D(NMAX_POINT, double);
    tau_yx = ARRAY_1D(NMAX_POINT, double);
    tau_yy = ARRAY_1D(NMAX_POINT, double);
    tau_yz = ARRAY_1D(NMAX_POINT, double);
    tau_zx = ARRAY_1D(NMAX_POINT, double);
    tau_zy = ARRAY_1D(NMAX_POINT, double);
    tau_zz = ARRAY_1D(NMAX_POINT, double);
  }

  #if GEOMETRY != CARTESIAN
   r  = grid[IDIR].x; th = grid[JDIR].x;
   if (one_dVr == NULL){
     one_dVr = ARRAY_1D(NX1_TOT, double);  /* -- intercell (i) and (i+1) volume -- */
     one_dmu = ARRAY_1D(NX2_TOT, double);  /* -- intercell (j) and (j+1) volume -- */
     for (i = 0; i < NX1_TOT - 1; i++){
       one_dVr[i] = r[i+1]*fabs(r[i + 1]) - r[i]*fabs(r[i]);
       one_dVr[i] = 2.0/one_dVr[i];
     }
     for (j = 0; j < NX2_TOT - 1; j++){
       one_dmu[j] = 1.0 - cos(th[j + 1]) - (1.0 - cos(th[j]))*(th[j] > 0.0 ? 1.0:-1.0);
       one_dmu[j] = 1.0/one_dmu[j];
     }
   }
  #endif

/* -- set pointers to coordinates and grid indices --*/

  x1 = grid[IDIR].x; x1r = grid[IDIR].xr; i = g_i; 
  x2 = grid[JDIR].x; x2r = grid[JDIR].xr; j = g_j; 
  x3 = grid[KDIR].x; x3r = grid[KDIR].xr; k = g_k; 

  if (g_dir == IDIR){   

  /* ---------------------------------------------------------
                    Loop on X1 direction 
     --------------------------------------------------------- */

    dy_1 = 1.0/grid[JDIR].dx[j]; 
    dz_1 = 1.0/grid[KDIR].dx[k];
    for (i = beg; i <= end; i++){
      dx_1 = 1.0/grid[IDIR].dx[i];

    /* -- compute face- and cell-centered values */

      VAR_LOOP(nv) {
        vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k][j][i+1]);
        vc[nv] = V[nv][k][j][i];
      }

    /* -- compute viscosity (face center) -- */
    
      Visc_nu(vi, x1r[i], x2[j], x3[k], &nu1, &nu2);

      tau_xx[i] = tau_xy[i]= tau_xz[i]= tau_yx[i]= 
      tau_yy[i] = tau_yz[i]= tau_zx[i]= tau_zy[i]= 
      tau_zz[i] = 0.0;

      dVxi = dVxj = dVxk = dVyi = dVyj = dVyk = dVzi = dVzj = dVzk = 0.0;
      dxVx = dxVy = dxVz = dyVx = dyVy = dyVz = dzVx = dzVy = dzVz = 0.0;

      EXPAND (dVxi = D_DX_I(Vx);, dVyi = D_DX_I(Vy);, dVzi = D_DX_I(Vz);)
      #if DIMENSIONS >= 2
       EXPAND (dVxj = D_DY_I(Vx);, dVyj = D_DY_I(Vy);, dVzj = D_DY_I(Vz);)
       #if DIMENSIONS == 3
        dVxk = D_DZ_I(Vx); dVyk = D_DZ_I(Vy); dVzk = D_DZ_I(Vz);   
       #endif
      #endif

      EXPAND(dxVx = dVxi*dx_1; , dxVy = dVyi*dx_1; , dxVz = dVzi*dx_1; )
      #if DIMENSIONS >= 2 
       EXPAND(dyVx = dVxj*dy_1; , dyVy = dVyj*dy_1; , dyVz = dVzj*dy_1; )
       #if DIMENSIONS == 3
        dzVx = dVxk*dz_1; dzVy = dVyk*dz_1; dzVz = dVzk*dz_1;
       #endif
      #endif

      #if GEOMETRY == CARTESIAN      
       div_v = D_EXPAND(dxVx, + dyVy , + dzVz); 

    /* -- stress tensor components (only some needed for flux computation) -- */

       tau_xx[i] = 2.0*nu1*dxVx + (nu2 - (2.0/3.0)*nu1)*div_v; /*2eta1 dxVx + (eta2 - 2/3 eta1) divV*/
       tau_xy[i] = nu1*(dyVx + dxVy);   /*eta1 (dyVx + dxVy)*/
       tau_xz[i] = nu1*(dzVx + dxVz);   /*eta1 (dzVx + dxVz)*/
       tau_yx[i] = tau_xy[i];                
       tau_zx[i] = tau_xz[i];                
       
    /* -- compute source terms -- */

       EXPAND(ViS[i][MX1] = 0.0; ,
              ViS[i][MX2] = 0.0; ,  
              ViS[i][MX3] = 0.0; )                              
                 
      #elif GEOMETRY == CYLINDRICAL

       r = grid[IDIR].x; dr = grid[IDIR].dx[i];
       dr_1 = 1.0/dr; dz_1 = 0.0;
       r_1 = 1.0/grid[IDIR].x[i];
       
    /* -- calculate the div(v) (trick @ axis for dxVx) -- */

       dxVx = (Vx[k][j][i+1]*r[i+1] - Vx[k][j][i]*fabs(r[i]))*one_dVr[i];
       div_v = D_EXPAND(  dxVx, + dVyj*dy_1, + 0.0); 
       dxVx = dVxi*dx_1;

    /* -- stress tensor components (only some needed for flux computation) -- */
                                
       tau_xx[i] = 2.0*nu1*dxVx + (nu2 - (2.0/3.0)*nu1)*div_v; /* tau_rr = 2 eta1 drVr + (eta2 - 2/3 eta1) divV    */
       tau_xy[i] = nu1*(dyVx + dxVy);   /*tau_rz = eta1 (dzVr + drVz)*/
       EXPAND(tau_xz[i] = nu1*(0.0);,
              tau_xz[i] = nu1*(0.0);,   
              tau_xz[i] = nu1*0.5*(r[i]+r[i+1])*dr_1*((1./r[i+1])*Vz[k][j][i+1] 
                                               - (1./r[i])*Vz[k][j][i]);)   
                               /*tau_rphi = eta1 (1/r dphiVr + drVphi - 1/r Vphi)= eta1 (r dr(1/r Vphi) )*/    
       tau_yx[i] = tau_xy[i]; /*tau_zr*/
       tau_zx[i] = tau_xz[i]; /*tau_phir*/
      
      /* -- we calculate at the center cause we don't need it for flux
            but for src, avoiding 1/r->inf at r=r_f=0 -- */

       Visc_nu(vc, x1[i], x2[j], x3[k], &nu1, &nu2);

       div_v = D_EXPAND( 0.5*(Vx[k][j][i+1]-Vx[k][j][i-1])*dx_1 + Vx[k][j][i]*r_1,
                       + 0.5*(Vy[k][j + 1][i]-Vy[k][j - 1][i])*dy_1, 
                       + 0.0 ); 
 
       tau_zz[i] = 2.0*nu1*r_1*Vx[k][j][i] 
                     + (nu2 - (2.0/3.0)*nu1)*div_v;
       
    /* -- compute source terms -- */

      EXPAND(ViS[i][MX1] = -tau_zz[i]*r_1; ,
             ViS[i][MX2] = 0.0;            ,  
             ViS[i][MX3] = 0.0; ) 

     #elif GEOMETRY == POLAR 
      
      r    = grid[IDIR].xr;    th = grid[JDIR].x;
      dr   = grid[IDIR].dx[i]; dr_1 = 1.0/dr;
      r_1  = 1.0/grid[IDIR].xr[i];

    /* -- calculate div(v) -- */

      div_v =  D_EXPAND(r_1*vi[VX1] + dVxi*dr_1, + r_1*dVyj*dy_1, + dVzk*dz_1); 

    /* -- stress tensor components (only some needed for flux computation) -- */

      tau_xx[i] = 2.0*nu1*dxVx + (nu2 - (2.0/3.0)*nu1)*div_v; /*tau_rr as in cylindrical */
      EXPAND(tau_xy[i] = nu1*dxVy;                              ,
             tau_xy[i] = nu1*(r_1*dyVx + dxVy - r_1*vi[VX2]);   ,
             tau_xy[i] = nu1*(r_1*dyVx + dxVy - r_1*vi[VX2]);)
                      /*tau_rphi = eta1 (1/r dphiVr + drVphi - 1/r Vphi) */
      tau_xz[i] = nu1*(dzVx + dxVz);  /*tau_rz as in cylindrical*/
      tau_yx[i] = tau_xy[i];               /*tau_phir = tau_rphi*/
      EXPAND(tau_yy[i] = 2.0*nu1*r_1*dyVy + (nu2 - (2.0/3.0)*nu1)*div_v;,
             tau_yy[i] = 2.0*nu1*(r_1*dyVy + r_1*vi[VX1])
                         + (nu2 - (2.0/3.0)*nu1)*div_v;,
             tau_yy[i] = 2.0*nu1*(r_1*dyVy + r_1*vi[VX1])
                         + (nu2 - (2.0/3.0)*nu1)*div_v;)
                        /*tau_phiphi = 2 eta1 (1/r dphiVphi + 1/r Vr) + (eta2 - 2/3 eta1) divV*/
      tau_zx[i] = tau_xz[i]; /*tau_zr = tau_rz*/

    /* -- compute source terms -- */

      r_1  = 1.0/grid[IDIR].x[i];
      
      EXPAND(ViS[i][MX1] = -0.5*(tau_yy[i - 1] + tau_yy[i])*r_1; ,
             ViS[i][MX2] = 0.0;                                  ,
             ViS[i][MX3] = 0.0;)
                                            
     #elif GEOMETRY == SPHERICAL

      r = grid[IDIR].xr; th = grid[JDIR].x;
      dr_1 = 1.0/grid[IDIR].dx[i]; r_1  = 1.0/grid[IDIR].xr[i];
      tan_1= 1.0/tan(th[j]); s_1 = 1.0/sin(th[j]);

      /* -- calculate div(v) -- */

      div_v = D_EXPAND(  2.0*r_1*vi[VX1] + dVxi*dr_1 ,
                      + r_1*dVyj*dy_1 + r_1*tan_1*vi[VX2], 
                      + r_1*s_1*dVzk*dz_1); 
       
    /* -- stress tensor components (only some needed for flux computation) -- */
                                
      /* tau_rr = 2eta1 drVr + + (eta2 - 2/3 eta1) divV */
      tau_xx[i] = 2.0*nu1*dxVx + (nu2 - (2.0/3.0)*nu1)*div_v;
      EXPAND(tau_xy[i] = nu1*(r_1*dyVx + dxVy );                   ,
             tau_xy[i] = nu1*(r_1*dyVx + dxVy - r_1*vi[VX2]);      ,
             tau_xy[i] = nu1*(r_1*dyVx + dxVy - r_1*vi[VX2]);)
                      /*tau_rthita = eta1 (1/r dthitaVr + drVthita -1/r Vthita)*/
      EXPAND(tau_xz[i] = nu1*(r_1*s_1*dzVx + dxVz) ;   ,  
             tau_xz[i] = nu1*(r_1*s_1*dzVx + dxVz );   ,   
             tau_xz[i] = nu1*(r_1*s_1*dzVx + dxVz - r_1*vi[VX3]);)
                      /*tau_rphi = eta1 (1/r dthitaVr + drVthita -1/r Vthita)*/
             tau_yx[i] = tau_xy[i]; /*tau_thitar*/
             tau_yy[i] = 2.0*nu1*(r_1*dyVy + r_1*vi[VX1]) 
                        + (nu2 - (2.0/3.0)*nu1)*div_v; 
                        /*tau_thitathita= 2 eta1 (1/r dthitaVthita + 1/r Vr) + (eta2 - 2/3 eta1)divV*/
      tau_zx[i] = tau_xz[i]; /*tau_phir*/
      EXPAND(tau_zz[i] = 2.0*nu1*(r_1*s_1*dzVz + r_1*vi[VX1] ) 
                        + (nu2 - (2.0/3.0)*nu1)*div_v;,
             tau_zz[i] = 2.0*nu1*(r_1*s_1*dzVz + r_1*vi[VX1]
                                         + tan_1*r_1*vi[VX2]) 
                        + (nu2 - (2.0/3.0)*nu1)*div_v;,
             tau_zz[i] = 2.0*nu1*(r_1*s_1*dzVz + r_1*vi[VX1]
                                         + tan_1*r_1*vi[VX2]) 
                        + (nu2 - (2.0/3.0)*nu1)*div_v;)
                        /*tau_phiphi= 2 eta1 (1/rs dphiVphi + 1/r Vr + 1/r cot Vthita) + (eta2 - 2/3 eta1)divV*/

    /* -- compute source terms -- */

      r_1  = 1.0/grid[IDIR].x[i];
      EXPAND(ViS[i][MX1] = -0.5*(  tau_yy[i - 1] + tau_yy[i]
                                 + tau_zz[i - 1] + tau_zz[i])*r_1;,
             ViS[i][MX2] = 0.0;,
             ViS[i][MX3] = 0.0;)       
               
      #endif  /* -- end #if GEOMETRY -- */

    /* -- compute fluxes -- */

      EXPAND(ViF[i][MX1] = tau_xx[i]; ,
             ViF[i][MX2] = tau_xy[i]; ,
             ViF[i][MX3] = tau_xz[i]; )

      #if EOS != ISOTHERMAL
       ViF[i][ENG] = EXPAND(  vi[VX1]*tau_xx[i],
                            + vi[VX2]*tau_yx[i],
                            + vi[VX3]*tau_zx[i]);
      #endif                                               

      dcoeff[i][MX1]  = MAX(nu1, nu2);
      dcoeff[i][MX1] /= vi[RHO];
    }

  }else if (g_dir == JDIR){ 

  /* ---------------------------------------------------------
                    Loop on X2 direction 
     --------------------------------------------------------- */

    dx_1 = 1.0/grid[IDIR].dx[i];
    dz_1 = 1.0/grid[KDIR].dx[k];
    for (j = beg ; j <= end; j++){
      dy_1 = 1.0/grid[JDIR].dx[j];

    /* -- compute face- and cell-centered values */

      VAR_LOOP(nv) {
        vi[nv] = 0.5*(V[nv][k][j+1][i] + V[nv][k][j][i]);
        vc[nv] = V[nv][k][j][i];
      }

    /* -- compute viscosity (face center) -- */
    
      Visc_nu(vi, x1[i], x2r[j], x3[k], &nu1, &nu2);

      tau_xx[j]= tau_xy[j]= tau_xz[j]= tau_yx[j] = tau_yy[j]= 
      tau_yz[j]= tau_zx[j]= tau_zy[j]= tau_zz[j] =0.0;
      
      dVxi = dVxj = dVxk = dVyi = dVyj = dVyk = dVzi = dVzj = dVzk = 0.0;
      dxVx = dxVy = dxVz = dyVx = dyVy = dyVz = dzVx = dzVy = dzVz = 0.0;

    /* ------ CALCULATE DERIVATIVES ---------- */

      EXPAND (dVxi = D_DX_J(Vx);, dVyi = D_DX_J(Vy);, dVzi = D_DX_J(Vz);)
      EXPAND (dVxj = D_DY_J(Vx);, dVyj = D_DY_J(Vy);, dVzj = D_DY_J(Vz);)
      #if DIMENSIONS == 3
       dVxk = D_DZ_J(Vx); dVyk = D_DZ_J(Vy); dVzk = D_DZ_J(Vz);   
      #endif
      
      EXPAND(dxVx = dVxi*dx_1; , dxVy = dVyi*dx_1; , dxVz = dVzi*dx_1; )
      EXPAND(dyVx = dVxj*dy_1; , dyVy = dVyj*dy_1; , dyVz = dVzj*dy_1; )
      #if DIMENSIONS == 3
       dzVx = dVxk*dz_1; dzVy = dVyk*dz_1; dzVz = dVzk*dz_1;
      #endif
     
      #if GEOMETRY == CARTESIAN
   
    /* -- calculate div(v) -- */
       
       div_v = D_EXPAND(dxVx, + dyVy, + dzVz); 
       
    /* -- stress tensor components (only some needed for flux computation) -- */

       tau_xy[j] = nu1*(dyVx + dxVy);
       tau_yx[j] = tau_xy[j];
       tau_yy[j] = 2.0*nu1*dyVy + (nu2 - (2.0/3.0)*nu1)*div_v;
       tau_yz[j] = nu1*(dyVz + dzVy);
       tau_zy[j] = tau_yz[j];
       
    /* -- compute source terms -- */

       EXPAND(ViS[j][MX1] = 0.0; ,
              ViS[j][MX2] = 0.0; ,  
              ViS[j][MX3] = 0.0; )                              

      #elif GEOMETRY == CYLINDRICAL
       
       r = grid[IDIR].x; th = grid[KDIR].x;
       dr = grid[IDIR].dx[i]; dr_1 = 1.0/dr;
       r_1 = 1.0/grid[IDIR].x[i]; dz_1 = 0.0;

     /* -- calculate div(v) -- */
       
       div_v = D_EXPAND(r_1*vi[VX1] + dVxi*dr_1, + dVyj*dy_1, + 0.0); 

    /* -- stress tensor components (only some needed for flux computation) -- */
      
       tau_xy[j] = nu1*(dyVx + dxVy); /*tau_rz = eta1 (dzVr + drVz)*/
       tau_yx[j] = tau_xy[j];  /*tau_zr */
       tau_yy[j] = 2.0*nu1*dyVy + (nu2 - (2.0/3.0)*nu1)*div_v; /*tau_zz = 2eta1 dzVz + (eta2 - 2/3 eta1)divV*/
       EXPAND(tau_yz[j] = 0.0;  ,
              tau_yz[j] = 0.0;  ,
              tau_yz[j] = nu1*(dyVz);) /*tau_zphi = eta1(1/r dphiVz + dzVphi)*/
       tau_zy[j] = tau_yz[j]; /*tau_phiz*/
       
    /* -- compute source terms -- */

       EXPAND(ViS[j][MX1] = 0.0; ,
              ViS[j][MX2] = 0.0; ,  
              ViS[j][MX3] = 0.0; )                              

      #elif GEOMETRY == POLAR 

       r = grid[IDIR].x; th = grid[JDIR].xr; dr   = grid[IDIR].dx[i];
       dr_1 = 1.0/dr; r_1  = 1.0/grid[IDIR].x[i];
       
    /* -- calculate div(v) -- */
       
       div_v = D_EXPAND(r_1*vi[VX1] + dVxi*dr_1, + r_1*dVyj*dy_1, + dVzk*dz_1); 

    /* -- stress tensor components (only some needed for flux computation) -- */

       tau_xy[j] = nu1*(r_1*dyVx + dxVy - r_1*vi[VX2]); 
                    /*tau_rphi = eta1(1/r dphiVr + drVphi -1/r Vphi)*/
       tau_yx[j] = tau_xy[j]; /*tau_phir*/
       tau_yy[j] = 2.0*nu1*(r_1*dyVy + r_1*vi[VX1])
                  + (nu2 - (2.0/3.0)*nu1)*div_v;
                  /*tau_phiphi = 2 eta1 (1/r dphiVphi + 1/rVr) + (eta2 - 2/3 eta1)divV*/
       tau_yz[j] = nu1*(r_1*dyVz + dzVy); /* tau_phiz = eta1(dzVphi + 1/r dphiVz)*/
       tau_zy[j] = tau_yz[j]; /*tau_zphi*/
     
       r_1  = 1.0/grid[IDIR].x[i];

    /* -- compute source terms -- */
       
       EXPAND(ViS[j][MX1] = 0.0; ,
              ViS[j][MX2] = 0.0; ,
              ViS[j][MX3] = 0.0; )
      
      #elif GEOMETRY == SPHERICAL
       
       r = grid[IDIR].x; th = grid[JDIR].xr;
       dr_1 = 1.0/grid[IDIR].dx[i]; r_1  = 1.0/grid[IDIR].x[i];
       tan_1= 1.0/tan(th[j]); s_1 = 1.0/sin(th[j]);

       /*------------------------------------
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Hotfix for tan_1 at the axis: 
        since we have terms tan_1*smth that 
        cannot be treated with the volume
        trick, we set an if condition 
        for this to go to zero, as it is
        the correct behaviour of the  
        term in question there.
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       --------------------------------------*/
       if (fabs(tan(th[j]))< 1.e-12) tan_1 = 0.0;      
       
       /* -- calculate div(v) -- */

       th = grid[JDIR].x;
       
       div_v = D_EXPAND( 2.0*r_1*vi[VX1] + dVxi*dr_1    ,
                       + ( sin(th[j + 1])*Vy[k][j + 1][i] - fabs(sin(th[j]))*Vy[k][j][i])*r_1*one_dmu[j], 
                       + r_1*s_1*dVzk*dz_1 ); 
      
       th = grid[JDIR].xr;

    /* -- stress tensor components (only some needed for flux computation) -- */

      tau_xy[j] = nu1*(r_1*dyVx + dxVy - r_1*vi[VX2]); 
                  /*tau_rthita = eta1(1/r dthitaVr + drVthita -1/r Vthita)*/
      tau_yx[j] = tau_xy[j]; /*tau_thitar*/
      tau_yy[j] = 2.0*nu1*(r_1*dyVy + r_1*vi[VX1]) 
                 + (nu2 - (2.0/3.0)*nu1)*div_v;
                 /*tau_thitathita= 2 eta1 (1/r dthitaVthita + 1/r Vr) + (eta2 - 2/3 eta1)divV*/
       #if DIMENSIONS <= 2                  
        EXPAND(tau_yz[j] = nu1*r_1*dyVz;  ,
               tau_yz[j] = nu1*r_1*dyVz;  ,
               tau_yz[j] = nu1*(r_1*dyVz - tan_1*r_1*vi[VX3]));
       #endif       
       #if DIMENSIONS == 3                  
        tau_yz[j] = nu1*(s_1*r_1*dzVy + r_1*dyVz - tan_1*r_1*vi[VX3]);
       #endif       /*tau_thitaphi= eta1 (1/rs dphiVthita + 1/r dthitaVphi -1/r cot Vphi)*/
       tau_zy[j] = tau_yz[j]; /*tau_phithita*/
       #if DIMENSIONS <= 2                  
        EXPAND(tau_zz[j] = 2.0*nu1*r_1*vi[VX1] 
                          + (nu2 - (2.0/3.0)*nu1)*div_v;    ,
               tau_zz[j] = 2.0*nu1*(r_1*vi[VX1] + tan_1*r_1*vi[VX2]) 
                          + (nu2 - (2.0/3.0)*nu1)*div_v;,
               tau_zz[j] = 2.0*nu1*(r_1*vi[VX1] + tan_1*r_1*vi[VX2]) 
                          + (nu2 - (2.0/3.0)*nu1)*div_v;)
       #endif       
       #if DIMENSIONS == 3                  
        tau_zz[j] = 2.0*nu1*(r_1*s_1*dzVz + r_1*vi[VX1] + tan_1*r_1*vi[VX2]) 
                   + (nu2 - (2.0/3.0)*nu1)*div_v;
       #endif       /*tau_phiphi = 2eta1(1/rs dphiVphi + 1/r Vr + 1/r cot Vthita) + (eta2 -2/3 eta1)divV */
                  
    /* -- compute source terms -- */

       th = grid[JDIR].x; tan_1= 1.0/tan(th[j]);
      
       EXPAND(ViS[j][MX1] = 0.0;                                              ,
              ViS[j][MX2] = 0.5*(  tau_yx[j - 1] + tau_yx[j])*r_1
                                 - tan_1*0.5*(tau_zz[j - 1] + tau_zz[j])*r_1; ,
              ViS[j][MX3] = 0.0; )

      #endif  /* -- end #if GEOMETRY -- */
      
    /* -- compute fluxes -- */

      EXPAND(ViF[j][MX1] = tau_yx[j]; ,
             ViF[j][MX2] = tau_yy[j]; ,
             ViF[j][MX3] = tau_yz[j]; )

      #if EOS != ISOTHERMAL
       ViF[j][ENG] = EXPAND(  vi[VX1]*tau_xy[j],
                            + vi[VX2]*tau_yy[j],
                            + vi[VX3]*tau_zy[j]);
      #endif

      dcoeff[j][MX1]  = MAX(nu1,nu2);
      dcoeff[j][MX1] /= vi[RHO];
    }      

  }else if (g_dir == KDIR){ 
   
  /* ---------------------------------------------------------
                    Loop on X3 direction 
     --------------------------------------------------------- */
     
    dx_1 = 1.0/grid[IDIR].dx[i];
    dy_1 = 1.0/grid[JDIR].dx[j]; 
    for (k = beg; k <= end; k++){
      dz_1 = 1.0/grid[KDIR].dx[k];

    /* -- compute face- and cell-centered values */

      VAR_LOOP(nv) {
        vi[nv] = 0.5*(V[nv][k][j][i] + V[nv][k+1][j][i]);
        vc[nv] = V[nv][k][j][i];
      }

    /* -- compute viscosity (face center) -- */
    
      Visc_nu(vi, x1[i], x2[j], x3r[k], &nu1, &nu2);

      tau_xx[k] = tau_xy[k] = tau_xz[k]= tau_yx[k]= tau_yy[k]= tau_yz[k]= 
      tau_zx[k] = tau_zy[k] = tau_zz[k]=0.0;
      dVxi = dVxj = dVxk = dVyi = dVyj = dVyk = dVzi = dVzj = dVzk = 0.0;

      dVxi = D_DX_K(Vx); dVyi = D_DX_K(Vy); dVzi = D_DX_K(Vz);
      dVxj = D_DY_K(Vx); dVyj = D_DY_K(Vy); dVzj = D_DY_K(Vz);
      dVxk = D_DZ_K(Vx); dVyk = D_DZ_K(Vy); dVzk = D_DZ_K(Vz);
  
      /* ------ CALCULATE DERIVATIVES ---------- */

      dxVx = dVxi*dx_1; dxVy = dVyi*dx_1; dxVz = dVzi*dx_1; 
      dyVx = dVxj*dy_1; dyVy = dVyj*dy_1; dyVz = dVzj*dy_1;
      dzVx = dVxk*dz_1; dzVy = dVyk*dz_1; dzVz = dVzk*dz_1;
         
      #if GEOMETRY == CARTESIAN
        
    /* -- calculate div(v) -- */
       
       div_v = dxVx + dyVy + dzVz;
      
    /* -- stress tensor components (only some needed for flux computation) -- */

       tau_xz[k] = nu1*(dzVx + dxVz);  
       tau_yz[k] = nu1*(dyVz + dzVy);
       tau_zx[k] = tau_xz[k];
       tau_zy[k] = tau_yz[k];
       tau_zz[k] = 2.0*nu1*dzVz + (nu2 - (2.0/3.0)*nu1)*div_v;

    /* -- compute source terms -- */

       EXPAND(ViS[k][MX1] = 0.0; ,
              ViS[k][MX2] = 0.0; ,  
              ViS[k][MX3] = 0.0; )                              

      #elif GEOMETRY == POLAR 

       r = grid[IDIR].x; th = grid[JDIR].x;
       dr   = grid[IDIR].dx[i]; dr_1 = 1.0/dr;
       r_1  = 1.0/grid[IDIR].x[i];

       /* -- calculate the div U -- */
       
       div_v = r_1*vi[VX1] + dVxi*dr_1 + r_1*dVyj*dy_1 + dVzk*dz_1; 

    /* -- stress tensor components (only some needed for flux computation) -- */

       tau_xz[k] = nu1*(dzVx + dxVz);  /*tau_rz = eta1(drVz + dzVr)*/
       tau_yy[k] = 2.0*nu1*(r_1*dyVy + r_1*vi[VX1])
                  + (nu2 - (2.0/3.0)*nu1)*div_v;
                  /*tau_phiphi = 2eta1(1/r dphiVphi + 1/r V) + (eta2- 2/3 eta1)divV*/
       tau_yz[k] = nu1*(r_1*dyVz + dzVy); /*tau_phiz = eta1(dzVphi + 1/r dphiVz)*/
       tau_zx[k] = tau_xz[k]; /* tau_zr   */
       tau_zy[k] = tau_yz[k]; /* tau_zphi */
       tau_zz[k] = 2.0*nu1*dzVz + (nu2 - (2.0/3.0)*nu1)*div_v;
                  /*tau_zz = 2eta1(dzVz) + (eta2- 2/3 eta1)divV*/

       EXPAND(ViS[k][MX1] = 0.0; ,
              ViS[k][MX2] = 0.0; ,  
              ViS[k][MX3] = 0.0; )                              
                   
      #elif GEOMETRY == SPHERICAL
       r    = grid[IDIR].x; th = grid[JDIR].x;
       dr_1 = 1.0/grid[IDIR].dx[i];
       r_1  = 1.0/grid[IDIR].x[i]; tan_1= 1.0/tan(th[j]); s_1 = 1.0/sin(th[j]);

    /* -- calculate div(v) -- */

       div_v =   2.0*r_1*vi[VX1] + dVxi*dr_1 + r_1*dVyj*dy_1 + r_1*tan_1*vi[VX2]
               + r_1*s_1*dVzk*dz_1; 

    /* -- stress tensor components (only some needed for flux computation) -- */

       tau_xz[k] = nu1*(r_1*s_1*dzVx + dxVz - r_1*vi[VX3]);  
                    /*tau_rphi = eta1(drVphi + 1/rs dphiVr - 1/r Vphi)*/ 
       tau_yz[k] = nu1*(s_1*r_1*dzVy + r_1*dyVz - tan_1*r_1*vi[VX3]);
                    /*tau_thitaphi = eta1(1/rs dphiVthita + 1/r dthitaVphi - 1/r cot Vphi)*/ 
       tau_zx[k] = tau_xz[k]; /* tau_phir */
       tau_zy[k] = tau_yz[k]; /* tau_phithita */
       tau_zz[k] = 2.0*nu1*(r_1*s_1*dzVz + r_1*vi[VX1]
                                    + tan_1*r_1*vi[VX2]) 
                  + (nu2 - (2.0/3.0)*nu1)*div_v;
                    /*tau_phiphi = 2eta1(1/rs dphiVphi +1/r Vr +1/r cot Vthita) + (eta2 - 2/3 eta1)divV*/

    /* -- compute source terms -- */

       EXPAND(ViS[k][MX1] = 0.0; ,
              ViS[k][MX2] = 0.0; ,  
              ViS[k][MX3] = 0.0; )                              
                                     
      #endif  /* -- end #if GEOMETRY -- */

    /* -- compute fluxes -- */

      ViF[k][MX1] = tau_zx[k];
      ViF[k][MX2] = tau_zy[k]; 
      ViF[k][MX3] = tau_zz[k]; 

      #if EOS != ISOTHERMAL
       ViF[k][ENG] =   vi[VX1]*tau_xz[k]
                     + vi[VX2]*tau_yz[k]
                     + vi[VX3]*tau_zz[k];
      #endif
      dcoeff[k][MX1]  = MAX(nu1, nu2);
      dcoeff[k][MX1] /= vi[RHO];
    }/*loop*/
  }/*sweep*/

}/*function*/
