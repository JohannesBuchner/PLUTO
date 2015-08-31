#include "pluto.h"
#define MAX_IMAGE_NUMBER 32

static Image image[MAX_IMAGE_NUMBER]; /* -- max number of images is 256 -- */


/* ************************************************************ */ 
void CreateImage (char *var_name)
/*
 * 
 * PUROPOSE:
 *
 *  Create an image associated with 'var_name'. 
 *  Use default 2D values.
 *  The images are collected and stored in this file 
 *  as a repository for later re-use.
 *
 ************************************************************** */
{
  static int icount = 0;

  sprintf (image[icount].basename,"%s",var_name);

  image[icount].slice_plane = X12_PLANE;
  image[icount].slice_coord = 0.0;

  image[icount].max = image[icount].min = 0.0;
  image[icount].logscale = 0;
  image[icount].colormap = "br";

  icount++; 
}

/* ********************************************************** */
Image *GetImage (char *var_name)
/* 
 * 
 * PURPOSE
 *
 *  Get an image structure associated with the variable 
 *  named 'var_name'.
 *
 *
 ************************************************************ */
{
  int indx = -1;

  while (strcmp(image[++indx].basename, var_name)) {
    if (image[indx].colormap == NULL) { /* if colormap is NULL, image does not exist! */
      print1 ("! Error: var '%s' is not associated with a valid image\n"); 
      QUIT_PLUTO(1);
    }
  }
   
  return (image + indx);
}

#undef MAX_IMAGE_NUMBER
