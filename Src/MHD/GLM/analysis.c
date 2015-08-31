{
  int i,j,k;
  static int n = 1;
  double ***divB, dB;
  double tot[2048], max;
  double atot, gtot, gmax;
  FILE *fp;

  divB = GLM_GetDivB();
  if (divB == NULL) {
    if (prank == 0){
      fp = fopen("divb.dat","w");
      fprintf (fp,"# \n");
      fprintf (fp,"# time     tot(t)   av_tot(t)    max\n");
      fprintf (fp,"# ------------------------------------\n");
      fclose(fp);
    }
    return;
  }
  tot[n] = 0.0;
  max    = 0.0;
  DOM_LOOP(k,j,i){
    dB      = fabs(divB[k][j][i]);
    tot[n] += dB;
    max     = MAX(dB, max);
  }
  D_EXPAND(
    tot[n] /= (double)grid[IDIR].np_int_glob;  ,
    tot[n] /= (double)grid[JDIR].np_int_glob;  ,
    tot[n] /= (double)grid[KDIR].np_int_glob;
  )

  #ifdef PARALLEL
   MPI_Allreduce (tot + n, &gtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce (&max,    &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   tot[n] = gtot;
   max    = gmax;
  #endif
 
/* -- compute average total -- */

  atot = 0.0;
  for (i = 1; i <= n; i++) atot += tot[i];
  atot /= (double)n;

/* -- write to file -- */

  if (prank == 0){
    printf ("> writing divB to disk (%d)\n",n);
    fp = fopen("divb.dat","a");
    fprintf (fp, "%f  %f  %f  %f\n",g_time, tot[n], atot, max);
    fclose(fp);
  }
  n++;

}

-----------------------------------------------------------
#! /bin/sh

for xx in "0.1" "0.2" "0.4" "0.8" "1.6" "3.2"
do
  echo "Doing $xx"
  ./pluto -par $xx > out.log
  /bin/cp divb.dat divb.$xx.dat
done

 
