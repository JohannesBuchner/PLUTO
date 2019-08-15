/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Read SCvH table.

   - read rho(T,P) and U(T,P)
   - Loop (P/rho, rho){  
        set T(P/rho, rho) by inverting rho(T,P)
     }
   - Loop (T, rho){
        Find P(T,rho), set U(T,rho) = interpolate U(T,P)
     }

  \author A. Mignone (mignone@ph.unito.it)
          B. Vaidya
  \date   9 June, 2014
*/
/* /////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#define P10(a)  pow(10.0, a)

void ReadSCvHTable ()
{
  int    i, j, k, nx, ny;
  double x, xmin, xmax, dx;
  double y, ymin, ymax, dy;
  double rho, p, T, c[11];
  double unit_pressure = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
  Table2D e_tab, rho_tab;
  FILE *fp;

/* --------------------------------------------------------------
   1.  Initialize rho(T,P)
   -------------------------------------------------------------- */

  xmin = 2.1; xmax = 7.06; dx = 0.08;   /* Log(T) */
  ymin = 4.0; ymax = 19.0; dy = 0.2;    /* Log(P) */

  nx = 1 + INT_FLOOR((xmax - xmin)/dx + 0.5); 
  ny = 1 + INT_FLOOR((ymax - ymin)/dy + 0.5);  

printf ("%f %f  nx = %d, ny = %d\n",(xmax-xmin)/dx,(ymax-ymin)/dy,nx,ny);

  InitializeTable2D(&rho_tab, P10(xmin), P10(xmax), nx, 
                              P10(ymin), P10(ymax), ny);
  InitializeTable2D(&e_tab, P10(xmin), P10(xmax), nx, 
                            P10(ymin), P10(ymax), ny);

  fp = fopen("H_TAB_I.A","r");
  if (fp == NULL){
    print ("! File not found\n");
    QUIT_PLUTO(1);
  }

  for (i = 0; i < rho_tab.nx; i++) for (j = 0; j < rho_tab.ny; j++) {
    rho_tab.f[j][i] = -90;  /* --  reset table to low values -- */
    rho_tab.defined[j][i] = 0;  
  }
  
  for (i = 0; i < nx; i++){
    fscanf (fp,"%lf  %d",&x, &ny);
    for (j = 0; j < ny; j++){
      for (k = 0; k < 11; k++) fscanf (fp,"%lf",c+k);
      printf ("(i,j) = (%d, %d), lnT = [%f %f], lnP = [%f %f]\n",
                i,j, x,rho_tab.lnx[i], c[0], rho_tab.lny[j]);
      rho = pow(10.0, c[3]); /* c[3] = Log(rho) in gr/cm^3*/
      rho_tab.f[j][i] = rho;
      rho_tab.defined[j][i] = 1;
        e_tab.f[j][i] = c[5];  /* c[5] = Log(e)   */
    }
  }
  fclose(fp);

  WriteBinaryTable2D("scvh_rho_Tp.bin", &rho_tab);
  WriteBinaryTable2D("scvh_e_Tp.bin",&e_tab);

/* --------------------------------------------------------------
   2. Loop over p/rho and rho to construct T(p, rho)
   -------------------------------------------------------------- */

  int status;
  Table2D Ttab;

  InitializeTable2D(&Ttab, P10(3.2), P10(20.1), 256,     /* p */
                           P10(-12), P10(4.065), 256);  /* rho */

/*
p   = P10(15.0);
rho = P10(-10.0);
status = InverseLookupTable2D(&rho_tab, p, rho, &T);
if (status == 0) Ttab.f[j][i] = T;
else {
  printf ("! Inverse not found, Log(p) = %f, Log(rho) = %f\n",
          log10(p), log10(rho));
  Ttab.f[j][i] = -99999.0;
} 
exit(1);
*/
  for (j = 0; j < Ttab.ny; j++){
  for (i = 0; i < Ttab.nx; i++){
    rho = Ttab.y[j]; 
    p   = Ttab.x[i];
    status = InverseLookupTable2D(&rho_tab, p, rho, &T);
    if (status == 0) Ttab.f[j][i] = T;
    else {
      if (log10(rho) > -3 && log10(rho) < -2){
        printf ("! Inverse not found, Log(p) = %f, Log(rho) = %f\n",
                log10(p), log10(rho));
      }
      Ttab.f[j][i] = -1.0;
    } 
  }}

  WriteBinaryTable2D("scvh_T_prho.bin", &Ttab);
  exit(1);
}
