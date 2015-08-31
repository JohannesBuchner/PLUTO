#include "pluto.h"

/* ******************************************************************** */
void GetGradient (double ***T, double **gradT, 
                  int beg, int end, Grid *grid)
/*
 *
 * PURPOSE
 *
 *   Compute the gradient of a 3D scalar quantity T in the direction
 *   given by g_dir.
 *   Return a 1D array (dT/dx, dT/dy, dT/dz) along that direction 
 *   computed at cell interfaces, e.g.
 *
 *   if g_dir == IDIR  --> compute 
 *  
 *    [ dT/dl1, dT/dl2, dT/dl3 ] at interface (i+1/2,j,k)
 *
 *   if g_dir == JDIR  --> compute 
 *  
 *    [ dT/dl1, dT/dl2, dT/dl3 ] at interface (i,j+1/2,k)
 * 
 *   if g_dir == KDIR  --> compute 
 *  
 *    [ dT/dl1, dT/dl2, dT/dl3 ] at interface (i,j,k+1/2)
 *   
 *
 *   Here dl1, dl2 and dl3 are the line element in the 
 *   thre directions:
 *
 *   Cartesian:   {dl1, dl2, dl3} = {dx,       dy, dz}
 *   Cylindrical: {dl1, dl2, dl3} = {dr,       dz, - }
 *   Polar:       {dl1, dl2, dl3} = {dr,   r.dphi, dz} 
 *   Spherical:   {dl1, dl2, dl3} = {dr, r.dtheta, r.sin(theta).dphi}
 *   
 * LAST MODIFIED
 *
 *   14 Apr 2011 by T. Matsakos, A. Mignone
 *
 *
 ************************************************************* */
{
  int  i,j,k;
  double *r, *rp;
  double *inv_dx,  *inv_dy,  *inv_dz;
  double *inv_dxi, *inv_dyi, *inv_dzi;
  double dl1, dl2, dl3, theta, r_1, s_1;
  double dx1, dx2, dx3;
  
  D_EXPAND(inv_dx  = grid[IDIR].inv_dx; inv_dxi = grid[IDIR].inv_dxi;  ,
           inv_dy  = grid[JDIR].inv_dx; inv_dyi = grid[JDIR].inv_dxi;  ,
           inv_dz  = grid[KDIR].inv_dx; inv_dzi = grid[KDIR].inv_dxi;)

  r  = grid[IDIR].x;
  rp = grid[IDIR].xr;

  i = g_i; j = g_j; k = g_k;

  if (g_dir == IDIR) {

    #if GEOMETRY == SPHERICAL
     theta = grid[JDIR].x[j];
     s_1   = 1.0/sin(theta);
    #endif
    D_EXPAND(                 ,
      dl2 = dx2 = inv_dy[j];  ,
      dl3 = dx3 = inv_dz[k];
    )
    for (i = beg; i <= end; i++){
      dl1 = inv_dxi[i];  
      #if GEOMETRY == POLAR && DIMENSIONS >= 2
       dl2 = dx2/rp[i]; 
      #elif GEOMETRY == SPHERICAL
       D_EXPAND(               ,
         dl2 = dx2/rp[i];      ,
         dl3 = dx3*s_1/rp[i];)
      #endif
      D_EXPAND( 
        gradT[i][0] = (T[k][j][i+1] - T[k][j][i])*dl1;         ,
        gradT[i][1] = 0.25*(  T[k][j+1][i] + T[k][j+1][i+1]
                            - T[k][j-1][i] - T[k][j-1][i+1])*dl2;  ,
        gradT[i][2] = 0.25*(  T[k+1][j][i] + T[k+1][j][i+1] 
                            - T[k-1][j][i] - T[k-1][j][i+1])*dl3;
      )
    }

  }else if (g_dir == JDIR) {

    r_1  = 1.0/r[i];
    D_EXPAND(
      dl1 = dx1 = inv_dx[i];   ,
                               ,
      dl3 = dx3 = inv_dz[k];
    )
    for (j = beg; j <= end; j++){
      dl2 = inv_dyi[j]; 
      #if GEOMETRY == POLAR
       dl2 *= r_1;
      #elif GEOMETRY == SPHERICAL
       D_EXPAND(               ,
                dl2  *= r_1;   ,
                theta = grid[JDIR].xr[j];
                dl3   = dx3*r_1/sin(theta);)
      #endif
       D_EXPAND( 
        gradT[j][0] = 0.25*(  T[k][j][i+1] + T[k][j+1][i+1]
                            - T[k][j][i-1] - T[k][j+1][i-1])*dl1;   ,
        gradT[j][1] = (T[k][j+1][i] - T[k][j][i])*dl2;              ,
        gradT[j][2] = 0.25*(  T[k+1][j][i] + T[k+1][j+1][i]
                            - T[k-1][j][i] - T[k-1][j+1][i])*dl3;
      )
    }
  
  }else if (g_dir == KDIR){

    dl1 = inv_dx[i];            
    dl2 = inv_dy[j]; 
    #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
     r_1  = 1.0/r[i];
     dl2 *= r_1; 
    #endif
    #if GEOMETRY == SPHERICAL
     theta = grid[JDIR].x[j];
     s_1   = 1.0/sin(theta);
    #endif

    for (k = beg; k <= end; k++){
      dl3 = inv_dzi[k]; 
      #if GEOMETRY == SPHERICAL
       dl3 *= r_1*s_1;
      #endif
      gradT[k][0] = 0.25*(  T[k][j][i+1] + T[k+1][j][i+1]
                          - T[k][j][i-1] - T[k+1][j][i-1])*dl1;
      gradT[k][1] = 0.25*(  T[k][j+1][i] + T[k+1][j+1][i]
                          - T[k][j-1][i] - T[k+1][j-1][i])*dl2;
      gradT[k][2] = (T[k+1][j][i] - T[k][j][i])*dl3;
    }
  }
}
