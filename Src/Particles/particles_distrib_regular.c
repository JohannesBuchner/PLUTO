/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Assign particle coordinates from distribution function

  Given the separable distribution functions
  \f$ f(x) = dN/dx,\,g(y) = dN/dy \f$ we use the cumulative distribution function  
  \f[
    F(x)G(y) = \frac{1}{N_x}\int_{x_b}^{x} f(x)\,dx\quad
               \frac{1}{N_x}\int_{y_b}^{y} g(y)\,dy
    \qquad\mathrm{where}\qquad
    N_x = \int_{x_b}^{x_e} f(x) dx,\,
    N_y = \int_{y_b}^{y_e} g(y) dy,\,
  \f]
  to solve
  \f[
    F(x_i) - R_i = 0
  \f]
   using a Newton-Raphson scheme:
  \f[
      x_i^{(k+1)} = x_i^{(k)} - \frac{F(x_i) - R_i}{f(x)/N} 
  \f]

  \authors A. Mignone (mignone@ph.unito.it)\n

 \b References
    - "Title" \n
      Authors, Journal (year) vol, page

  \date  Oct 15, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER    16

static double Particles_dNdx(double, void *);
static double Particles_dNdy(double, void *);

/* ********************************************************************* */
void Particles_Distrib1(int ip, int Np, double coor[])
/*!
 *  
 *
 *********************************************************************** */
{
  int    i,j,k;
  int    Np_dim;
  double x, y, dx, dy;
  double Nx, rx, Fx, fx;
  double Ny, ry, Gy, gy;
  double Nz, rz, Hz, hz;

  double xb = g_domBeg[IDIR], xe = g_domEnd[IDIR], Lx = xe - xb;
  double yb = g_domBeg[JDIR], ye = g_domEnd[JDIR], Ly = ye - yb;
  double zb = g_domBeg[KDIR], ze = g_domEnd[KDIR], Lz = ze - zb;

  double tol = 1.e-5, ftol = 1.e6;


  Np_dim = ceil(sqrt(Np));
  i      = ip%Np_dim;
  j      = ip/Np_dim;
  
/* -- Generate equally spaced numbers between 0 and 1 -- */

  Nx = GaussQuadrature(Particles_dNdx, NULL, xb, xe, 128, 5);  
  Ny = GaussQuadrature(Particles_dNdy, NULL, yb, ye, 128, 5);  
  rx = (double)(ip+1)/(double)(Np+1);   /* Generate number between 0 and 1 */
//  rx = (double)(i+1)/(double)(Np_dim+1);   /* Generate number between 0 and 1 */
//  ry = (double)(j+1)/(double)(Np_dim+1);   /* Generate number between 0 and 1 */


// $$$$$$$ Check cumulative distribution is computed correctly  $$$$$$$$$$$$$$
static int first_call = 1;
FILE *fp;
if (first_call){
  fp = fopen("distrib.dat","w");
  for (x = xb; x <= xe; x += (xe-xb)/128.){
    fprintf (fp,"%f  %f\n",
              x,GaussQuadrature(Particles_dNdx, NULL, xb, x, 128, 5)/Nx);
    
  }
  fclose(fp);
}
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

print (">> Beg Distrib(). Particle #%d; i,j = %d,%d; rx,ry = %f, %f\n",
        ip, i,j,rx,ry);
  
/* -- Start main iteration loop (x) -- */

  x = xb + rx*Lx;  /* Guess */
/*  
  for (k = 0; k < MAX_ITER; k++){
    Fx  = GaussQuadrature(Particles_dNdx, NULL, xb, x, 128, 5);
    fx  = Particles_dNdx(x, NULL);
    dx  = (Fx - rx*Nx)/fx;
    x  -= dx;
  print ("  iter %d, x = %f, Nx = %f\n",k,x, Fx);
    if (fabs(dx) < tol*Lx || fabs(Fx-rx*Nx) < ftol) {
      break;
    }
  }
*/
/* -- Use bisection method for x -- */

  double xL, xR, funcL, funcR, func;
  xL = xb;
  xR = xe;

  Fx  = GaussQuadrature(Particles_dNdx, NULL, xb, xL, 128, 5);
  fx  = Particles_dNdx(xL, NULL);
  funcL = (Fx - rx*Nx);

  Fx  = GaussQuadrature(Particles_dNdx, NULL, xb, xR, 128, 5);
  fx  = Particles_dNdx(xR, NULL);
  funcR = (Fx - rx*Nx);
    
  for (k = 0; k < MAX_ITER; k++){
    x    = 0.5*(xL+xR);
    Fx   = GaussQuadrature(Particles_dNdx, NULL, xb, x, 128, 5);
    fx   = Particles_dNdx(x, NULL);
    func = (Fx - rx*Nx);
    if (funcL*func < 0.0){
      funcR = func;
      xR    = x;
    }else{
      funcL = func;
      xL    = x;
    }

    if (fabs(xL-xR) < tol*Lx) {
      break;
    }
  }


// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if (first_call) fp = fopen("rx.dat","w");
else            fp = fopen("rx.dat","a");

fprintf (fp,"%d  %f  %f\n",ip,rx,x);
fclose(fp);
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



/* -- Start main iteration loop (x) -- */

  y  = yb + ry*Ly;  /* Guess */
  for (k = 0; k < MAX_ITER; k++){
    Gy  = GaussQuadrature(Particles_dNdy, NULL, yb, y, 128, 5);
    gy  = Particles_dNdy(y, NULL);
    dy  = (Gy - ry*Ny)/gy;
    y  -= dy;
    if (fabs(dy) < tol*Ly || fabs(Gy-ry*Ny) < ftol) {
      break;
    }
  }



print ("<< End Distrib()\n");

  coor[IDIR] = x;
  coor[JDIR] = y*0;
  coor[KDIR] = 0.0;


  
first_call = 0;  
}

double Particles_dNdx(double x, void *param)
{
  return exp(-(x-0.5)*(x-0.5)/0.05);
}

double Particles_dNdy(double y, void *param)
{
  return exp(-y*y/0.05);
}
