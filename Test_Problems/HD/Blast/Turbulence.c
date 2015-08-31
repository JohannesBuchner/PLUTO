/* ***************************************************************
 *  Routine that generates a turbulent 
 *  density field.
 *
 *
 * ************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Points_X		128
#define Points_Y		128
#define x_end			1.
#define Norm			.1		/* Normalization of the spectrum power-law */
#define Spectral_Index	-5./3.	/* Index of the spectrum power-law: Kolmogorov -> 5/3 */
#define	IR_Cutoff		2.		/* Needed to avoid an IR divergence; 2 is a good value */
#define Minimum_Lambda	4		/* Grid points corresponding to minimum wavelength */

#define frand()((double)(rand()+1)/(RAND_MAX+1.0))

double **array_2D(int nx, int ny);

int main(void)
{
  int i, j, k, i2, j2, k2, te0, NXLim, NYLim;
  double y_end, xdbl, NXLim_2, NYLim_2;
  double twoPI, x1r, x2r, twopi_a, twopi_b, twopi_c, scrh1, scrh2;
  double xL, xR;
  double **Amplitude, **Phase;
	
  FILE *Out;

/* --------------------------------------------------------
             Grid generation
   -------------------------------------------------------- */

  Out = fopen("grid0.out", "w");
  fprintf(Out, "# GEOMETRY:   CARTESIAN\n");

  fprintf(Out, "%d\n",Points_X);
  for(i2 = 0; i2 < Points_X; i2++){
    xL = -0.5 + (double)(x_end*(i2-0.5)/(Points_X - 1.));
    xR = -0.5 + (double)(x_end*(i2+0.5)/(Points_X - 1.));
    fprintf (Out, "%d   %12.6e  %12.6e \n",i2+1, xL, xR);
  }

  fprintf(Out, "%d\n",Points_Y);
  for(j2 = 0; j2 < Points_Y; j2++){
    xL = -0.5 + (double)(x_end*(j2-0.5)/(Points_Y - 1.));
    xR = -0.5 + (double)(x_end*(j2+0.5)/(Points_Y - 1.));
    fprintf (Out, "%d   %12.6e  %12.6e \n",j2+1, xL, xR);
  }

  fprintf(Out, "%d\n",1);  /* write 3rd coordinates as well */
  fprintf (Out, "1  0.0  1.0\n");
  fclose(Out);

/* --------------------------------------------------------
               Data generation
   -------------------------------------------------------- */

  te0 = time(NULL);
  srand(te0);
  NXLim = 2*((int)(Points_X/2./Minimum_Lambda));
  NYLim = 2*((int)(Points_Y/2./Minimum_Lambda));	
  NXLim_2   = NXLim/2.;
  NYLim_2   = NYLim/2.;
  y_end = x_end*(Points_Y - 1.)/(Points_X - 1.);

  Amplitude = array_2D(NXLim + 1, NYLim + 1);
  Phase     = array_2D(NXLim + 1, NYLim + 1);
  twoPI     = 2.*3.14159265359;
	
  for(i = 0; i <= NXLim; i++) for(j = 0; j <= NYLim; j++){
    x1r = frand();
    x2r = frand();
    xdbl = sqrt(-2.*log(x1r))*cos(twoPI*x2r);
    Phase[i][j] = twoPI*frand();
    scrh1       = IR_Cutoff*IR_Cutoff + pow(twoPI*(i - NXLim_2)/x_end,2.)
                                      + pow(twoPI*(j - NYLim_2)/y_end,2.);
    Amplitude[i][j] = sqrt(Norm*pow(scrh1, Spectral_Index/2.))*xdbl;
  }
  Out = fopen("rho0.dbl", "wb");

  for(j2 = 0; j2 < Points_Y; j2++){
    twopi_b = twoPI/(Points_Y - 1.)*j2;
    for(i2 = 0; i2 < Points_X; i2++){
      twopi_a = twoPI/(Points_X - 1.)*i2;
      xdbl = 0.;
      for(i = NXLim + 1; i--; ){
        scrh1 = twopi_a*(i - NXLim_2); 
        for(j = NYLim + 1; j--; ){
          scrh2 = twopi_b*(j - NYLim_2);
          xdbl += Amplitude[i][j]*cos(scrh1 + scrh2 + Phase[i][j]);
        }
      }
      xdbl = exp(xdbl);
      fwrite (&xdbl, sizeof(double), 1, Out);
    }
  }

  fclose(Out);
  return;
}
    
/* ************************************************************ *    
 *    Functions for memmory Allocation
 * ************************************************************ */    
 
    double **array_2D(int nx, int ny)
    {
	int i;
	double **m;
	
	m = (double **)malloc((size_t)nx*sizeof(double *));
	m[0] = (double *)malloc((size_t)nx*ny*sizeof(double));

	for(i = 1; i < nx; i++) m[i] = m[i - 1] + ny;
	
	return m;
    }
    
    
