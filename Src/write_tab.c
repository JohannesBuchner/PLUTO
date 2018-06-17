/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Write tabulated 1D or 2D ascii data files.

  WriteTabArray() provides ascii output for 2-D and 1-D data 
  in tabulated multi-column format:
  \verbatim
      .     .         .         .             .
      .     .         .         .             .
      .     .         .         .             .
    x1(i) x2(j)   var_1(i,j)  var_2(i,j)   var_3(i,j) ...
      .     .         .         .             .
      .     .         .         .             .
      .     .         .         .             .
  \endverbatim
 
  Blocks with different x2 coordinates are separated by a blank row
  (in more than 1D).
  One file (with all variables) per output is written to disk. 
  This format does not work in parallel mode and can be used for simple
  serial runs.
 
  \author A. Mignone (mignone@ph.unito.it)
  \date   April 6, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void WriteTabArray (Output *output, char *filename, Grid *grid)
/*!
 * Write tabulated array.
 *
 * \param [in] output    a pointer to the output structure corresponding
 *                       to the tab format
 * \param [in] filename  the output file name
 * \param [in] grid      pointer to an array of Grid structures
 *********************************************************************** */
{
  int  nv, i, j, k;
  FILE *fout;

  #ifdef PARALLEL
   print ("! WriteTabArray: tab output not supported in parallel\n");
   return;
  #endif

  #if DIMENSIONS == 3
   print ("! WriteTabArray: tab output not supported in 3D\n");
   return;
  #endif

/* ------------------------------------------
      dump arrays in ascii format to disk 
   ------------------------------------------ */

  fout = fopen (filename, "w");
  k = 0;
  IDOM_LOOP (i){
    JDOM_LOOP(j){
      fprintf (fout, "%f %f ", grid->x[IDIR][i], grid->x[JDIR][j]);
      for (nv = 0; nv < output->nvar; nv++) {
        if (output->dump_var[nv]) 
          fprintf (fout, "%12.6e ", output->V[nv][k][j][i]);
      }
      fprintf (fout, "\n");   /* newline */
    }
    #if DIMENSIONS > 1
     fprintf (fout, "\n");   /* skip one more empty line in 2D */
    #endif
  }
  fclose (fout);
}
