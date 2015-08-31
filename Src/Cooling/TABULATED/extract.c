#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main()
{
  int    n, i;
  char   label[64];
  double z[14], L[14], T;
  FILE *f, *fout;

//  f = fopen("cooltab","r");
  f = fopen("cooltab_z03z3.dat","r"); /* newer file */
  if (f == NULL){
    printf ("! File cooltab does not exists\n");
    exit(1);
   }

/* -- read column header -- */
/*
  for (i = 0; i < 14; i++){
    fscanf(f, "%lf  \n", z + i); 
  }
*/
  fscanf(f, "%s  %lf  %lf  %lf  \n", label, z, z+1, z+2); 

  fout = fopen("cooltable.dat","w");
  while (!feof(f)) {
/*
    fscanf (f,"%d  %lf", &n, &T);
    for (i = 0; i < 14; i++){
      fscanf(f, "%lf  \n", L + i); 
    }
*/
    fscanf (f,"%lf  %lf  %lf  %lf", &T, L, L+1, L+2);
    fprintf (fout,"%12.6e   %12.6e\n",T, L[1]*z[1]);
  }
  fclose(f);
  fclose(fout);


  return(0);
}
