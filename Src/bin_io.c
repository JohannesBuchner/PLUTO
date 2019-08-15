/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Functions for handling binary I/O.

  This file provides a number of handy functions for opening, reading
  and writing binary files using single or double precision in serial
  or parallel mode.
  It is employed by the following output formats: .dbl, .flt and .vtk.
  
  In parallel mode these functions work as wrappers to the actual
  parallel implementations contained in AL_io.c.
 
  Pointer to data array must be cast into (void *) and are assumed to
  start with index 0 for both cell-centered and staggered data arrays.
  This means that if V3D is a 3D array then:
  - for a cell-centered array, <tt>V = (void *) V3D[0][0]</tt>
  - for a x-staggered array,   <tt>V = (void *)(V3D[0][0]-1)</tt>
  - for a y-staggered array,   <tt>V = (void *)(V3D[0][-1])</tt>
  - for a z-staggered array,   <tt>V = (void *)(V3D[-1][0])</tt>

  \authors A. Mignone (mignone@ph.unito.it)\n
           G. Muscianisi (g.muscianisi@cineca.it)

  \date   Oct 22, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
int FileClose (FILE *fbin, int sz) 
/*!
 * Close file.
 *
 * \param [in] fbin     pointer to the FILE that needs to be closed (serial
 *                      mode only)
 * \param [in] sz       the distributed array descriptor for parallel mode
 *
 * \return Returns 0 on success
 *********************************************************************** */
{
  #ifdef PARALLEL 
   AL_File_close(sz); 
  #else
   fclose (fbin);
  #endif
  return(0);
}

/* ********************************************************************* */
float ***Convert_dbl2flt (double ***Vdbl, double unit, int swap_endian)
/*!
 * Convert the a double-precision 3D array into single precision.
 * Swap endianity if swap_endian == 1.
 *
 * \param [in] Vdbl         pointer to a 3D double precision aray
 * \param [in] unit         a multiplicative constant typically used 
 *                          to write in c.g.s units.
 * \param [in] swap_endian  when set to 1, swap endianity during the 
 *                          conversion.
 * \return a pointer to a 3D array in single precision.
 *********************************************************************** */
{
  int i, j, k;
  float  flt;
  static float ***Vflt;
  
  if (Vflt == NULL) Vflt = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, float);

  if (!swap_endian){
    DOM_LOOP(k,j,i){
      Vflt[k][j][i] = (float)(Vdbl[k][j][i]*unit);
    }
  }else{
    DOM_LOOP(k,j,i){
      flt = (float)(Vdbl[k][j][i]*unit);
      SWAP_VAR(flt);
      Vflt[k][j][i] = flt;
    }
  }

  return (Vflt);
}

/* ********************************************************************* */
FILE *FileOpen (char *filename, int sz, char *mode)  
/*!
 * Open a file for write/reading in binary mode.
 *
 * \param [in] filename a valid file name
 * \param [in] sz       the distributed array descriptor. This parameter
 *                      replaces dsize in parallel mode
 * \param [in] mode     a string giving the opening mode (only "w" or 
 *                      "r" are allowed)
 *
 * \return The pointer to the file.
 *********************************************************************** */
{
  FILE *fp;

/* ----------------------------------------------
    when reading, check if file exists
   ---------------------------------------------- */
   
  if (strcmp(mode,"r") == 0 && prank == 0){
    fp = fopen(filename, "rb");
    if (fp == NULL){
      print ("! FileOpen(): file %s does not exists\n", filename);
      QUIT_PLUTO(1);
    }
    fclose(fp);
  }
  
/* ------------------------------------------------ 
          file exists, keep going
   ------------------------------------------------ */
  
#ifdef PARALLEL
/*  CANNOT USE IT HERE SINCE FILE IS OPENED AND CLOSED FOR APPENDING SEVERAL TIMES
  if (strcmp(mode,"w") == 0){
print ("Deleting file...\n");
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
                  MPI_INFO_NULL, &fh);
    MPI_File_close(&fh);
  }  
*/      
  AL_File_open(filename, sz); 
  return NULL;
#else
  if      (strcmp(mode,"w") == 0) fp = fopen(filename, "wb");
  else if (strcmp(mode,"r") == 0) fp = fopen(filename, "rb");
  if (fp == NULL){
    print ("! Cannot find file %s\n",filename);
    QUIT_PLUTO(1);
  }
  return (fp);
#endif
}

/* ********************************************************************* */
void FileReadData (void *V, size_t dsize, int sz, FILE *fl, int istag,
                   int swap_endian) 
/*!
 * Read a 3D grid data array V[k][j][i] from a binary file. 
 *
 * \param [in] V      pointer to a 3D array,
 *                    V[k][j][i] -->  V[i + NX1*j + NX1*NX2*k]. Must
 *                    start with index 0
 * \param [in] dsize  the size of the each buffer element   
 *                    (sizeof(double) or sizeof(float)). This
 *                    parameter is used only in serial mode.
 * \param [in] sz     the distributed array descriptor. This parameter
 *                    replaces dsize in parallel mode
 * \param [in] fl     a valid FILE pointer
 * \param [in] istag  a flag to identify cell-centered (istag = -1) or
 *                    staggered field data (istag = 0,1 or 2 for staggering 
 *                    in the x1, x2 or x3 directions) 
 * \param [in]  swap_endian  a flag for swapping endianity
 *
 * \return This function has no return value.
 *
 * \remark The data array V is assumed to start \a always with index 0
 *         for  both cell-centered and staggered arrays.
 *********************************************************************** */
{
  int i, j, k;
  int ioff, joff, koff;
  char *Vc;

/* ---------------------------------------
    parallel reading handled by ArrayLib
   --------------------------------------- */

#ifdef PARALLEL 
  AL_Read_array (V, sz, istag);
#else

/* ---------------------------------------
      serial reading
   --------------------------------------- */

  Vc = (char *) V;
  ioff = (istag == 0); 
  joff = (istag == 1); 
  koff = (istag == 2);

  for (k = KBEG; k <= KEND + koff; k++) {
  for (j = JBEG; j <= JEND + joff; j++) {
    i = IBEG + (NX1_TOT + ioff)*(j + (NX2_TOT + joff)*k);
    fread (Vc + i*dsize, dsize, NX1 + ioff, fl);
  }}
#endif

/* -----------------------------------
     now swap endian if necessary
   ----------------------------------- */

  if (swap_endian){
    double *Vd;
    Vd = (double *) V;
    ioff = (istag == 0); 
    joff = (istag == 1); 
    koff = (istag == 2);

    for (k = KBEG-koff; k <= KEND + koff; k++) {
    for (j = JBEG-joff; j <= JEND + joff; j++) {
    for (i = IBEG-ioff; i <= IEND + ioff; i++) {
      SWAP_VAR(Vd[i + (NX1_TOT + ioff)*(j + (NX2_TOT + joff)*k)]);
    }}}     
  }
}

/* ********************************************************************* */
void FileWriteData (void *V, size_t dsize, int sz, FILE *fl, int istag)
/*!
 * Write a 3D grid data array V[k][j][i] to a binary file. 
 *
 * \param [in] V      single pointer to a 3D array,
 *                    V[k][j][i] -->  V[i + NX1*j + NX1*NX2*k]. Must
 *                    start with index 0
 * \param [in] dsize  the size of the each buffer element   
 *                    (sizeof(double) or sizeof(float)). This
 *                    parameter is used only in serial mode.
 * \param [in] sz     the distributed array descriptor. This parameter
 *                    replaces dsize in parallel mode
 * \param [in] fl     a valid FILE pointer
 * \param [in] istag  a flag to identify cell-centered (istag = -1) or
 *                    staggered field data (istag = 0,1 or 2 for staggering 
 *                    in the x1, x2 or x3 directions) 
 *
 * \return This function has no return value.
 *
 * \remark The data array is assumed to start \a always with index 0 for
 *         both cell-centered and staggered arrays.
 *********************************************************************** */
{
  int i, j, k;
  int ioff, joff, koff;
  char *Vc;

#ifdef PARALLEL
  MPI_Barrier (MPI_COMM_WORLD);
  AL_Write_array (V, sz, istag);
  return;
#else
  Vc = (char *) V;
  ioff = (istag == 0); 
  joff = (istag == 1); 
  koff = (istag == 2);

  for (k = KBEG; k <= KEND + koff; k++) {
  for (j = JBEG; j <= JEND + joff; j++) {
    i = IBEG + (NX1_TOT + ioff)*(j + (NX2_TOT + joff)*k);
    fwrite (Vc + i*dsize, dsize, NX1 + ioff, fl);
  }}
#endif
}


/* ********************************************************************* */
int FileDelete (char *fname)
/*!
 * Delete a file (parallel mode only).
 * Useful on some architecture where MPI_File_open will just overwrite
 * the file without erasing it.
 * 
 *
 * \param [in] filename a valid file name
 *
 * \return The pointer to the file.
 *********************************************************************** */
{
#ifdef PARALLEL
  MPI_File fp;

  MPI_File_open(MPI_COMM_WORLD, fname,
                MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
                MPI_INFO_NULL, &fp);
  MPI_File_close(&fp);
#endif
}

/* ********************************************************************* */
void FileWriteArray(void *buf, long int offset, int cnt, size_t size, char *fname)
/*!
 *  Write an array to file using raw binary format.
 *  No AL distributed array descriptor is used here.
 *
 * \param [in]   buf      address of buffer to be written
 * \param [in]   offset   file offset (in bytes)
 * \param [in]   cnt      the number of element to be written
 * \param [in]   size     the sizeof() datatype
 * \param [in]   fname    the file name 
 *********************************************************************** */
{
  char *cbuf = (char *) buf;
  cnt *= size;  /* Convert number of elements to byte */

#ifdef PARALLEL
  MPI_File fhw;
  MPI_Status status;
  MPI_Offset start;   /* Store end of file position */
  
  
  MPI_File_open(MPI_COMM_WORLD, fname, 
                MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fhw);
  MPI_File_get_position(fhw, &start);
    
  MPI_File_set_view (fhw, start+offset, MPI_BYTE, MPI_CHAR, "native", MPI_INFO_NULL);
  MPI_File_write_all(fhw, cbuf, cnt, MPI_CHAR, &status);
  MPI_File_close(&fhw);
#else
  FILE *fp = fopen (fname,"ab"); /* -- Open for appending in binary mode -- */
  fwrite (cbuf, sizeof(char), cnt, fp);
  fclose (fp);    
#endif
  
}

/* ********************************************************************* */
void FileWriteHeader(char *buffer, char fname[], int mode)
/*!
 *  Write or Append an ASCII header to a file.
 *
 * \param [in]   buffer  the input buffer (a string). 
 * \param [in]   fname   the file name
 * \param [in]   mode    an integer specifying if the file has to be
 *                       opened for the first time (mode = -1) or
 *                       appended (mode >= 0)
 * 
 *********************************************************************** */
{
  int nelem = strlen(buffer);
#ifdef PARALLEL
  MPI_File fhw;
  MPI_Status status;
  MPI_Offset start; /* Used to store end of file position */

  if (mode < 0){
  
  /* -- Open file and delete it -- */
  
    MPI_File_open(MPI_COMM_WORLD, fname,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE,
                  MPI_INFO_NULL, &fhw);
    MPI_File_close(&fhw);
    
  /* -- Open file for writing -- */    

    MPI_File_open(MPI_COMM_WORLD, fname, 
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);
    start = 0;
  }else{
    MPI_File_open(MPI_COMM_WORLD, fname, 
                  MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &fhw);
    MPI_File_get_position(fhw, &start);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_set_view(fhw, start, MPI_BYTE, MPI_CHAR, "native", MPI_INFO_NULL);
  if( prank == 0 ){
    MPI_File_write(fhw, buffer, nelem, MPI_CHAR, &status); 
  }
  MPI_File_close(&fhw);
#else
  FILE *fp;
  if (mode < 0){
    fp = fopen (fname,"wb"); /* -- Open for the first time -- */
  }else{
    fp = fopen(fname, "ab"); /* -- for appending in binary mode -- */
  }
  fwrite (buffer, sizeof(char), nelem, fp);
  fclose (fp);    
#endif
}
