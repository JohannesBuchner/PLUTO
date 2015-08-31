/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Provide basic functionality for reading input data files.

  Collects a number of functions for opening, reading and assigning 
  initial conditions from user-supplied data file(s).
  The geometry and dimensions of the input grid can be different from 
  the actual grid employed by PLUTO, as long as the coordinate geometry
  transformation has been implemented.
  The input grid and data files should employ the same format and 
  conventions employed by PLUTO. 
  
  - Gridfile: coordinates should be written using the PLUTO 4.0 grid format.
  - Datafile: variables should be written in sequence in a single binary 
              file using single or double precision. 
              The file extension must be ".flt" or ".dbl" for the former and
              the latter, respectively.
     
    Note that not all of the variables should be present and the input
    array ::get_var specifies which ones are to be searched for.

  The InputDataSet() initialize the module and by assigning values to 
  global variables such as size, geometry and dimensions of the input grid.
  Data values are read through the function InputDataRead() while
  InputDataInterpolate() can be finally used to map input data onto the
  grid employed by PLUTO using bi- or tri-linear interpolation to fill the 
  data array at the desired coordinate location.

  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos 
  \date   Aug 27, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"
#define ID_MAX_NVAR 256 

static int id_nvar; /**< Number of variables to be read on input. */
static int id_var_indx[ID_MAX_NVAR]; /**< The variable index. */
static int id_nx1; /**< Size of input grid in the x1 direction. */
static int id_nx2; /**< Size of input grid in the x2 direction. */
static int id_nx3; /**< Size of input grid in the x3 direction. */

static int id_geom;   /**< Geometry of the input grid. */

static double *id_x1; /**< Array of point coordinates of the x1 input grid. */
static double *id_x2; /**< Array of point coordinates of the x2 input grid. */
static double *id_x3; /**< Array of point coordinates of the x3 input grid. */

static double ***Vin[ID_MAX_NVAR]; /**< An array of 3D data values containing the
                                        initial data file variables. */

/* ********************************************************************* */
void InputDataSet (char *grid_fname, int *get_var)
/*!
 * Initialize access to input data file by assigning values to 
 * grid-related information (geometry, number of points, etc...).
 * This function should be called just once for input-data initialization.
 *
 * \param [in] gname the grid file name
 * \param [in] get_var an array of integers specifying which variables
 *                     have to be read from the input data. 
 * \return Thi function has no return value.
 *
 * The following tasks are performed.
 *********************************************************************** */
{
  int    i, ip, nv, success;
  size_t dsize, dcount;
  char *sub_str, sline[256];
  const char delimiters[] = " \t\r\f\n";
  double xl, xr;
  fpos_t file_pos;
  FILE *fp;

  print1 ("> Input data:\n\n");

/* --------------------------------------------------------------------- */
/*! - Scan grid data file and try to determine the grid geometry 
      (::id_geom). Search for tag "GEOMETRY:" and read the word that
      follows.                                                           */
/* --------------------------------------------------------------------- */

  fp = fopen(grid_fname,"r");
  if (fp == NULL){
    print1 ("! InputDataSet: grid file %s not found\n",grid_fname);
    QUIT_PLUTO(1);
  }
  success = 0;
  while(!success){
    fgets(sline,512,fp);
    sub_str = strtok(sline, delimiters);
    while (sub_str != NULL){
      if (!strcmp(sub_str,"GEOMETRY:")) {
        sub_str = strtok(NULL, delimiters);     
        success = 1;
        break;
      }
      sub_str = strtok(NULL, delimiters);     
    }  
  }
  
  if      (!strcmp(sub_str,"CARTESIAN"))   id_geom = CARTESIAN;
  else if (!strcmp(sub_str,"CYLINDRICAL")) id_geom = CYLINDRICAL;
  else if (!strcmp(sub_str,"POLAR"))       id_geom = POLAR;
  else if (!strcmp(sub_str,"SPHERICAL"))   id_geom = SPHERICAL;
  else{
    print1 ("! InputDataSet: unknown geometry\n");
    QUIT_PLUTO(1);
  }
    
  print1 ("  Input grid file:       %s\n", grid_fname);
  print1 ("  Input grid geometry:   %s\n", sub_str);

/* --------------------------------------------------------------------- */
/*! - Move file pointer until the first line that does not
      begin with a "#".                                                  */
/* --------------------------------------------------------------------- */

  success = 0;
  while(!success){
    fgetpos(fp, &file_pos);
    fgets(sline,512,fp);
    if (sline[0] != '#') success = 1;
  }

  fsetpos(fp, &file_pos);
  
/* --------------------------------------------------------------------- */
/*! - Start reading number of points and grid coordinates. For the
      input x1 direction these are stored inside the module variables
      ::id_nx1 and ::id_x1.                                              */
/* --------------------------------------------------------------------- */
   
  fscanf (fp,"%d \n",&id_nx1);
  id_x1 = ARRAY_1D(id_nx1, double);
  for (i = 0; i < id_nx1; i++){
    fscanf(fp,"%d  %lf %lf\n", &ip, &xl, &xr);
    id_x1[i] = 0.5*(xl + xr);
  }

  fscanf (fp,"%d \n",&id_nx2);
  id_x2 = ARRAY_1D(id_nx2, double);
  for (i = 0; i < id_nx2; i++){
    fscanf(fp,"%d  %lf %lf\n", &ip, &xl, &xr);
    id_x2[i] = 0.5*(xl + xr);
  }

  fscanf (fp,"%d \n",&id_nx3);
  id_x3 = ARRAY_1D(id_nx3, double);
  for (i = 0; i < id_nx3; i++){
    fscanf(fp,"%d  %lf %lf\n", &ip, &xl, &xr);
    id_x3[i] = 0.5*(xl + xr);
  }
  fclose(fp);

/* -- reset grid with 1 point -- */

  if (id_nx1 == 1) id_x1[0] = 0.0;
  if (id_nx2 == 1) id_x2[0] = 0.0;
  if (id_nx3 == 1) id_x3[0] = 0.0;
 
  print1 ("  Input grid extension:  x1 = [%12.3e, %12.3e] (%d points)\n",
             id_x1[0], id_x1[id_nx1-1], id_nx1);
  print1 ("\t\t\t x2 = [%12.3e, %12.3e] (%d points)\n",
             id_x2[0], id_x2[id_nx2-1], id_nx2);
  print1 ("\t\t\t x3 = [%12.3e, %12.3e] (%d points)\n",
             id_x3[0], id_x3[id_nx3-1], id_nx3);

  
/* --------------------------------------------------------------------- */
/*! - Find out how many and which variables we have to read (:id_nvar 
      and ::id_var_indx). 
      Stop counting variables as soon as the first occurrence of "-1" 
      in get_var is encountered                                          */
/* --------------------------------------------------------------------- */
   
  id_nvar = 0;
  for (nv = 0; nv < ID_MAX_NVAR; nv++){
    if (get_var[nv] != -1) {
      id_nvar++;
      id_var_indx[nv] = get_var[nv];
    }else{
      break;
    }
  }
  print1 ("  Number of variables:   %d\n",id_nvar);
}

/* ********************************************************************* */
void InputDataRead (char *data_fname, char *endianity)
/*!
 * Read input data file and store the contents into the local storage
 * array ::Vin. Memory allocation is also done here.  
 * The grid size and number of variables must have 
 * previously set by calling InputDataSet().
 * 
 * \param [in] data_fname the data file name
 * \param [in] endianity  an input string ("little" or "big") giving 
 *                        the byte-order of how the input data file 
 *                        was originally written.
 *                        If an empty string is supplied, no change is 
 *                        made.
 * \return This function has no return value.
 *********************************************************************** */
{
  int   i, j, k, nv, swap_endian=NO;
  size_t dsize, dcount;
  double udbl;
  float  uflt;
  char   ext[] = "   ";
  FILE *fp;

/* ----------------------------------------------------
             Check endianity 
   ---------------------------------------------------- */
  
  if ( (!strcmp(endianity,"big")    &&  IsLittleEndian()) ||
       (!strcmp(endianity,"little") && !IsLittleEndian())) {
    swap_endian = YES;
  }
  
  print1 ("  Input data file:       %s (endianity: %s) \n", 
           data_fname, endianity);
  
/* ------------------------------------------------------
    Get data type from file extensions (dbl or flt).
   ------------------------------------------------------ */

  dcount = strlen(data_fname);
  for (i = 0; i < 3; i++) ext[i] = data_fname[dcount-3+i];

  if (!strcmp(ext,"dbl")){
    print1 ("  Precision:             (double)\n");
    dsize = sizeof(double);
  } else if (!strcmp(ext,"flt")) {
    print1 ("  Precision:\t\t  (single)\n");
    dsize = sizeof(float);
  } else {
    print1 ("! InputDataRead: unsupported data type '%s'\n",ext);
    QUIT_PLUTO(1);
  }
  
/* -------------------------------------------------------
     Read and store data values
   ------------------------------------------------------- */

  fp = fopen(data_fname, "rb");
  if (fp == NULL){
    print1 ("! InputDataRead: file %s does not exist\n");
    QUIT_PLUTO(1);
  }
  for (nv = 0; nv < id_nvar; nv++){
    if (Vin[nv] == NULL) Vin[nv] = ARRAY_3D(id_nx3, id_nx2, id_nx1, double);

    dcount  = 1;

    if (dsize == sizeof(double)){
      for (k = 0; k < id_nx3; k++){ 
      for (j = 0; j < id_nx2; j++){
      for (i = 0; i < id_nx1; i++){
        if (fread (&udbl, dsize, dcount, fp) != dcount){
          print1 ("! InputDataRead: error reading data %d.\n",nv);
          break;
        }
        if (swap_endian) SWAP_VAR(udbl);
        Vin[nv][k][j][i] = udbl;
      }}}
    }else{
      for (k = 0; k < id_nx3; k++){ 
      for (j = 0; j < id_nx2; j++){
      for (i = 0; i < id_nx1; i++){
        if (fread (&uflt, dsize, dcount, fp) != dcount){
          print1 ("! InputDataRead: error reading data %d.\n",nv);
          break;
        }
        if (swap_endian) SWAP_VAR(uflt);
        Vin[nv][k][j][i] = uflt;
      }}}
    }
  }
  fclose(fp);
  print1 ("\n");
}

/* ********************************************************************* */
void InputDataInterpolate (double *vs, double x1, double x2, double x3)
/*!
 * Perform bi- or tri-linear interpolation on external
 * dataset to compute vs[] at the given point {x1,x2,x3}.
 *
 * \param [in] vs interpolated value
 * \param [in] x1 coordinate point at which at interpolates are desired
 * \param [in] x2 coordinate point at which at interpolates are desired
 * \param [in] x3 coordinate point at which at interpolates are desired
 * \return This function has no return value.
 *
 * The function performs the following tasks. 
 *********************************************************************** */
{
  int il = 0, jl = 0, kl = 0;
  int ih, jh, kh;
  int im, jm, km;
  int i, j, k, nv, inv;
  double xx, yy, zz;
  double ***V;

/* --------------------------------------------------------------------- */
/*! - Convert PLUTO coordinates to input grid geometry if necessary.     */
/* --------------------------------------------------------------------- */

  #if GEOMETRY == CARTESIAN
   if (id_geom == GEOMETRY) {  

     /* same coordinate system: nothing to do */
     
   }else if (id_geom == CYLINDRICAL) {  
     double R, z, phi;
     R   = sqrt(x1*x1 + x2*x2);
     phi = atan2(x2,x1);
     if (phi < 0.0) phi += 2.0*CONST_PI;
     z   = x3;

     x1 = R; x2 = z; x3 = phi;
   }else if (id_geom == POLAR) {  
     double R, phi, z;
     R   = sqrt(x1*x1 + x2*x2);
     phi = atan2(x2,x1);
     if (phi < 0.0) phi += 2.0*CONST_PI;
     z   = x3;

     x1 = R; x2 = phi; x3 = z;
   }else if (id_geom == SPHERICAL){
     double r, theta, phi;
     r     = D_EXPAND(x1*x1, + x2*x2, + x3*x3);
     r     = sqrt(r);
     theta = acos(x3/r);
     phi   = atan2(x2,x1);
     if (phi   < 0.0) phi   += 2.0*CONST_PI;
     if (theta < 0.0) theta += 2.0*CONST_PI;
     
     x1 = r; x2 = theta; x3 = phi;
   }else{
     print1 ("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
     QUIT_PLUTO(1);
   }
  #elif GEOMETRY == CYLINDRICAL
   if (id_geom == GEOMETRY) {  

     /* same coordinate system: nothing to do */
     
   }else if (id_geom == SPHERICAL) {  
     double r, theta, phi;
     r     = D_EXPAND(x1*x1, + x2*x2, + 0.0);
     r     = sqrt(r);
     theta = acos(x2/r);
     phi   = 0.0;
     if (theta < 0.0) theta += 2.0*CONST_PI;
     
     x1 = r; x2 = theta; x3 = phi;
   }else{
     print1 ("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
     QUIT_PLUTO(1);
   }
  #elif GEOMETRY == POLAR
   if (id_geom == GEOMETRY) {  

     /* same coordinate system: nothing to do */
     
   }else if (id_geom == CARTESIAN) {  
     double x, y, z;
     x = x1*cos(x2);
     y = x1*sin(x2);
     z = x3;
     
     x1 = x; x2 = y; x3 = z;
   }else{
     print1 ("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
     QUIT_PLUTO(1);
   }   
  #elif GEOMETRY == SPHERICAL
   if (id_geom == GEOMETRY) {  

     /* same coordinate system: nothing to do */
     

   }else if (id_geom == CARTESIAN) {  
     double x, y, z;
     x = x1*sin(x2)*cos(x3);
     y = x1*sin(x2)*sin(x3);
     z = x1*cos(x2);
     
     x1 = x; x2 = y; x3 = z;
   }else{
     print1 ("! InputDataInterpolate: invalid or unsupported coordinate transformation.\n");
     QUIT_PLUTO(1);
   }   
  #endif

/* --------------------------------------------------------------------- */
/*! - Make sure point (x1,x2,x3) does not fall outside input grid range. 
      Limit to input grid edge otherwise.                                */
/* --------------------------------------------------------------------- */
   
  D_EXPAND(if      (x1 < id_x1[0])         x1 = id_x1[0];
           else if (x1 > id_x1[id_nx1-1]) x1 = id_x1[id_nx1-1];  ,
           
           if      (x2 < id_x2[0])         x2 = id_x2[0];
           else if (x2 > id_x2[id_nx2-1]) x2 = id_x2[id_nx2-1];  ,
           
           if      (x3 < id_x3[0])         x3 = id_x3[0];
           else if (x3 > id_x3[id_nx3-1]) x3 = id_x3[id_nx3-1]; )

/* --------------------------------------------------------------------- */
/*! - Use table lookup by binary search to  find the indices 
      il, jl and kl such that grid points of PLUTO fall between 
      [il, il+1], [jl, jl+1], [kl, kl+1].                                */
/* --------------------------------------------------------------------- */

  il = 0;
  ih = id_nx1 - 1;
  while (il != (ih-1)){
    im = (il+ih)/2;
    if   (x1 <= id_x1[im]) ih = im;   
    else                    il = im;
  }
  
  if (id_nx2 > 1){
    jl = 0;
    jh = id_nx2 - 1;
    while (jl != (jh-1)){
      jm = (jl+jh)/2;
      if (x2 <= id_x2[jm]) jh = jm;   
      else                  jl = jm;
    }
  }

  if (id_nx3 > 1){
    kl = 0;
    kh = id_nx3 - 1;
    while (kl != (kh - 1)){
      km = (kl+kh)/2;
      if (x3 <= id_x3[km]) kh = km;   
      else                  kl = km;
    }
  }

/* --------------------------------------------------------------------- */
/*! - Define normalized coordinates between [0,1]:
      - x[il+1] < x1[i] < x[il+1] ==> 0 < xx < 1
      - y[jl+1] < x2[j] < y[jl+1] ==> 0 < yy < 1
      - z[kl+1] < x3[k] < z[kl+1] ==> 0 < zz < 1                            */
/* --------------------------------------------------------------------- */

  xx = yy = zz = 0.0; /* initialize normalized coordinates */

  if (id_nx1 > 1) xx = (x1 - id_x1[il])/(id_x1[il+1] - id_x1[il]);  
  if (id_nx2 > 1) yy = (x2 - id_x2[jl])/(id_x2[jl+1] - id_x2[jl]);  
  if (id_nx3 > 1) zz = (x3 - id_x3[kl])/(id_x3[kl+1] - id_x3[kl]);

/* --------------------------------------------------------------------- */
/*! - Perform bi- or tri-linear interpolation.                           */
/* --------------------------------------------------------------------- */

  for (nv = 0; nv < id_nvar; nv++) { 
    inv = id_var_indx[nv];
    V = Vin[nv];
    vs[inv] =   V[kl][jl][il]*(1.0 - xx)*(1.0 - yy)*(1.0 - zz)
              + V[kl][jl][il+1]*xx*(1.0 - yy)*(1.0 - zz);
    if (id_nx2 > 1){
      vs[inv] +=   V[kl][jl+1][il]*(1.0 - xx)*yy*(1.0 - zz)
                 + V[kl][jl+1][il+1]*xx*yy*(1.0 - zz);
    }
    if (id_nx3 > 1){
     vs[inv] +=   V[kl+1][jl][il]*(1.0 - xx)*(1.0 - yy)*zz
                + V[kl+1][jl][il+1]*xx*(1.0 - yy)*zz
                + V[kl+1][jl+1][il]*(1.0 - xx)*yy*zz
                + V[kl+1][jl+1][il+1]*xx*yy*zz;
    }
  }
}

/* ********************************************************************* */
void InputDataFree (void)
/*!
 * Free memory stored by user-supplied data.
 *
 *********************************************************************** */
{
  int nv;
  for (nv = 0; nv < id_nvar; nv++){
    free ((char *) Vin[nv][0][0]);    
    free ((char *) Vin[nv][0]);
    free ((char *) Vin[nv]);
  }
}

