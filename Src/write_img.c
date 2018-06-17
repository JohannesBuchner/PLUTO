#include "pluto.h"
#ifdef USE_PNG
 #include <png.h>
#endif

static void GetSlice (double ***, Image *, Grid *);
static int GET_SLICE_INDEX (int, double, Grid *);

/* ****************************************************************** */
void WritePPM (double ***Vdbl, char *var_name, char *filename, 
                Grid *grid)

/*
 *
 *
 *
 *
 ******************************************************************** */
{
  int  ic, ir;
  FILE *fl;
  char header[512];
  Image *ppm;

  ppm = GetImage (var_name);
  SetColorMap (ppm->r, ppm->g, ppm->b, ppm->colormap);
  GetSlice (Vdbl, ppm, grid);
  if (prank != 0) return;

  sprintf (header,"P6\n%d %d\n255\n", ppm->ncol, ppm->nrow);
  fl  = fopen (filename,"w");
  fprintf(fl,"%s",header);
  for (ir = 0; ir < ppm->nrow; ir++){
    fwrite (ppm->rgb[ir], sizeof(RGB), ppm->ncol, fl);
  } 
  fclose(fl);
}

#ifdef USE_PNG
/* ****************************************************************** */
void WritePNG (double ***Vdbl, char *var_name, char *filename, 
                Grid *grid)
/*
 *
 *
 *
 *
 ******************************************************************** */
{
  int ic, ir, i;
  png_structp     png_ptr;
  png_infop       info_ptr;
  int backgroundcolour_;
  int bit_depth_;
  int colortype_;
  int compressionlevel_;
  int  indx;
  double filegamma_;
  Image *png;
  unsigned char **image;
  FILE   *fp;
    
  png = GetImage (var_name);
  SetColorMap (png->r, png->g, png->b, png->colormap);
  GetSlice (Vdbl, png, grid);
  if (prank != 0) return;

  image = (png_bytepp)malloc(png->nrow*sizeof(png_bytep));
  for (ir = 0; ir < png->nrow; ir++) {
    image[ir] = (png_bytep)malloc(6*png->ncol*sizeof(png_byte));
  }

  for(ic = 0; ic < png->ncol; ic++){
  for(ir = 0; ir < png->nrow; ir++){
    i = 6*ic;
    image[ir][i]   = png->rgb[ir][ic].r;      /* -- red -- */
    image[ir][i+1] = 0;
    image[ir][i+2] = png->rgb[ir][ic].g;     /* -- green -- */
    image[ir][i+3] = 0;
    image[ir][i+4] = png->rgb[ir][ic].b;     /* -- blue -- */
    image[ir][i+5] = 0;
  }}


 /* -- write --- */

  compressionlevel_ = 6;
  backgroundcolour_ = 0;
  bit_depth_ = 16;
  filegamma_ = 0.;
  colortype_ = 2.; 

  fp = fopen(filename, "wb");
  if(fp == NULL){
    printf(" ! error opening file in writing data\n");
    exit(1);
  }

  png_ptr  = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info_ptr = png_create_info_struct(png_ptr);
  png_init_io(png_ptr, fp);

   /*   if(compressionlevel_ != -2){ */
        png_set_compression_level(png_ptr, compressionlevel_);
   /* }
   else
     {
        png_set_compression_level(png_ptr, PNGWRITER_DEFAULT_COMPRESSION);
   }*/

  png_set_IHDR(png_ptr, info_ptr, png->ncol, png->nrow,
               bit_depth_, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

  if(filegamma_ < 1.0e-1){
    filegamma_ = 0.7;
  }

  png_set_gAMA(png_ptr, info_ptr, filegamma_);
  
   /*
   time_t          gmt;
   png_time        mod_time;
   png_text        text_ptr[5];
   time(&gmt);
   png_convert_from_time_t(&mod_time, gmt);
   png_set_tIME(png_ptr, info_ptr, &mod_time);
   */

  png_write_info(png_ptr, info_ptr);
  png_write_image(png_ptr, image); 
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fp);

  free((char *) image[0]);
  free((char *) image);

}

#endif

/* *********************************************************** */
void GetSlice (double ***Vdbl, Image *image, Grid *grid)
/*
 * 
 * PURPOSE
 *
 *  get a 2D slice from the 3D array Vdbl.
 *  The array is actually written to disk and
 *  re-open for reading immediately after by
 *  proc #0. 
 *  Store its content as 2D rgb structure inside 
 *  image->rgb.
 *   
 *
 * 
 ************************************************************* */
{
  int i, j, k;
  int nx, ny, nz;
  int col_offset, row_offset;
  int offset, ir, ic;
  char filename[256];
  float xflt, slice_min, slice_max;
  static float **slice;
  static RGB **rgb;
  FILE *fl;
  size_t dsize = sizeof(float);

  #if DIMENSIONS == 1
   print ("! PPM output disabled in 1-D\n");
   return;    
  #endif    

/* ------------------------------------------------ 
    Write the whole array in single precision
    Slice are post-processed later 
   ------------------------------------------------ */
  
  fl = FileOpen ("tmp_file.out", SZ_float, "w");
  FileWriteData ((Convert_dbl2flt(Vdbl,1.0, 0))[0][0], dsize, SZ_float, fl, -1); 
  FileClose (fl, SZ_float);
  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
  #endif

  if (prank != 0) return; /* -- rank 0 will do the rest -- */

/* ------------------------------------------------
          get global dimensions
   ------------------------------------------------ */

  nx = grid->gend[IDIR] + 1 - grid->nghost[IDIR];
  ny = grid->gend[JDIR] + 1 - grid->nghost[JDIR];
  nz = grid->gend[KDIR] + 1 - grid->nghost[KDIR];

/* -----------------------------------------
     Allocate memory: make slices big
     enough to contain slices of different 
     sizes.
   ----------------------------------------- */

  if (slice == NULL) {
    ic = MAX(nx,ny); 
    ir = MAX(ny,nz);
    slice = ARRAY_2D(ir, ic, float);
    rgb   = ARRAY_2D(ir, ic, RGB);
  }

/* --------------------------------------------
       set offsets for slice reading
   -------------------------------------------- */

  fl = fopen ("tmp_file.out", "rb");
  if (image->slice_plane == X12_PLANE){

    k = GET_SLICE_INDEX (image->slice_plane, image->slice_coord, grid);
    offset = k*nx*ny*dsize; 
    fseek (fl, offset, SEEK_SET);
    col_offset = 0;
    row_offset = 0;
    image->ncol = nx;
    image->nrow = ny;

  } else if (image->slice_plane == X13_PLANE){

    j = GET_SLICE_INDEX (image->slice_plane, image->slice_coord, grid);
    offset = j*nx*dsize;
    fseek (fl, offset, SEEK_CUR);
    col_offset = 0;
    row_offset = nx*(ny - 1)*dsize; 
    image->ncol = nx;
    image->nrow = nz;

  } else if (image->slice_plane == X23_PLANE){
    
    i = GET_SLICE_INDEX (image->slice_plane, image->slice_coord, grid);
    offset = i*dsize;
    fseek (fl, offset, SEEK_CUR);
    col_offset = (nx - 1)*dsize; 
    row_offset = 0;
    image->ncol = ny;
    image->nrow = nz;

  }

/* -----------------------------------------
            Read slice  
   ----------------------------------------- */

  for (ir = 0; ir < image->nrow; ir++) {
    for (ic = 0; ic < image->ncol; ic++) {     
      fread (&xflt, dsize, 1, fl);
      if (feof(fl)){
         printf (" ! end of file reached\n");
         QUIT_PLUTO(1);
      }
      slice[image->nrow - 1 - ir][ic] = xflt;  /* -- swap row order -- */
      if (col_offset > 0) fseek(fl, col_offset, SEEK_CUR);
    }
    if (row_offset > 0) fseek(fl, row_offset, SEEK_CUR);
  }
  fclose(fl);

/* -----------------------------------------------------------
         Get slice max and min 
   ----------------------------------------------------------- */

  if (fabs(image->max - image->min) < 1.e-8) { /* -- set auto-scale -- */
    slice_min =  1.e38;
    slice_max = -1.e38;
    for (ir = 0; ir < image->nrow; ir++){
    for (ic = 0; ic < image->ncol; ic++){
      slice_min = MIN(slice_min, slice[ir][ic]);
      slice_max = MAX(slice_max, slice[ir][ic]);
    }}
  } else {
    slice_min = image->min;
    slice_max = image->max;
  }

/* -----------------------------------------------------
      Scale the image between 0 and 255.
      Set log/linear scaling. 
      Create {r,g,b} triplet 
   ----------------------------------------------------- */
  
  image->rgb = rgb;
  for (ir = 0; ir < image->nrow; ir++){
  for (ic = 0; ic < image->ncol; ic++){

    if (image->logscale) {
      xflt  = log10(slice[ir][ic]/slice_min)*255.;
      xflt /= log10(slice_max/slice_min);
    }else{
      xflt  = (slice[ir][ic] - slice_min)*255.;
      xflt /= (slice_max - slice_min);
    }

    i = (int) xflt;
    i = MIN (i, 255);
    i = MAX (i, 0);

    image->rgb[ir][ic].r = image->r[i];
    image->rgb[ir][ic].g = image->g[i];
    image->rgb[ir][ic].b = image->b[i];
  }}
}

/* ********************************************************* */
int GET_SLICE_INDEX (int plane, double x, Grid *grid)
/*
 *
 *
 * Find in which cell the coordinate x falls into
 *
 *********************************************************** */
{
  int    dir, i;
  double xr,xl;

  #if DIMENSIONS == 2
   return 0;
  #endif

  if (plane == X12_PLANE) dir = KDIR;
  if (plane == X23_PLANE) dir = IDIR;
  if (plane == X13_PLANE) dir = JDIR;

  for (i = 0; i < grid->np_tot_glob[dir]; i++){
    xl = grid->x_glob[dir][i] - 0.5*grid->dx_glob[dir][i];
    xr = grid->x_glob[dir][i] + 0.5*grid->dx_glob[dir][i];
    if (x >= xl && x <= xr) return MAX(i-IBEG,0);
  }
  return 0; 

}

