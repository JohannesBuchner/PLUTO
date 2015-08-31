/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Memory allocation functions.

  Provides essential functions for allocating and deallocating 
  multi-dimensional arrays.   
  The functions Array1D(), Array2D(), Array3D(), Array4D() can be 
  used to allocate storage for 1-D, 2-D, 3-D and 4-D arrays 
  of any data type with indices starting at 0.
  They are typically invoked from within a corresponding macros
  that handles type casting automatically, e.g., 
  \code
     char *s;
     double **q;
     s = ARRAY_1D(20, char);
     q = ARRAY_2D(30,40,double);
  \endcode
  will allocate memory for a 1D char array with \c 20 elements and a 
  2D double arrays of \c 30x40 elements
  
  The function ArrayBox() can be used to allocate memory for 
  a double precision array with specified index range.

  The functions ArrayMap() can be used to convert a one-dimensional
  array into a 3D array.
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 18, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#define NONZERO_INITIALIZE YES /* Fill arrays to nonsense values to catch
                                  uninitialized values later in the code */

/* ********************************************************************* */
void FreeArray1D (void *v)
/*! 
 * Free memory allocated by the pointer *v.
 *
 *********************************************************************** */
{
  free ((char *)v);
}
/* ********************************************************************* */
void FreeArray2D (void **m)
/*! 
 * Free memory allocated by the double pointer **v.
 *
 *********************************************************************** */
{
  free ((char *) m[0]);
  free ((char *) m);
}
/* ********************************************************************* */
void FreeArray3D (void ***m) 
/*! 
 * Free memory allocated by the pointer ***v.
 *
 *********************************************************************** */
{
  free ((char *) m[0][0]);
  free ((char *) m[0]);
  free ((char *) m);
}
/* ********************************************************************* */
void FreeArray4D (void ****m)
/*! 
 * Free memory allocated by the pointer ****v.
 *
 *********************************************************************** */
{
  free ((char *) m[0][0][0]);
  free ((char *) m[0][0]);
  free ((char *) m[0]);
  free ((char *) m);
}

/* ********************************************************************* */
char *Array1D (int nx, size_t dsize)
/*! 
 * Allocate memory for a 1-D array of any basic data type starting 
 * at 0.
 *
 * \param [in] nx    number of elements to be allocated
 * \param [in] dsize data-type of the array to be allocated
 * 
 * \return A pointer of type (char *) to the allocated memory area. 
 *********************************************************************** */
{
  char *v;
  v = (char *) malloc (nx*dsize);
  PlutoError (!v, "Allocation failure in Array1D");
  g_usedMemory += nx*dsize;

  #if NONZERO_INITIALIZE == YES
   if (dsize==sizeof(double)){
     int i;
     double *q;
     q = (double *) v;
     for (i = nx; i--; ) q[i] = i*1.e18 - 1.e16;
   }
  #endif
  return v;
}

/* ********************************************************************* */
char **Array2D (int nx, int ny, size_t dsize)
/*! 
 * Allocate memory for a 2-D array of any basic data type.
 *
 * \param [in] nx    number of elements in the 2nd dimension
 * \param [in] ny    number of elements in the 1st dimension
 * \param [in] dsize data-type of the array to be allocated
 * 
 * \return A pointer of type (char **) to the allocated memory area 
 *         with index range [0...nx-1][0...ny-1]
 *          
 *********************************************************************** */
{
  int i;
  char **m;
 
  m    = (char **)malloc ((size_t) nx*sizeof(char *));
  PlutoError (!m, "Allocation failure in Array2D (1)");
  m[0] = (char *) malloc ((size_t) nx*ny*dsize);
  PlutoError (!m[0],"Allocation failure in Array2D (2)");
 
  for (i = 1; i < nx; i++) m[i] = m[(i - 1)] + ny*dsize;
 
  g_usedMemory += nx*ny*dsize;

  #if NONZERO_INITIALIZE == YES
   if (dsize==sizeof(double)){
     int j;
     double **q;
     q = (double **)m;
     for (i = nx; i--; ){
     for (j = ny; j--; ){
       q[i][j] = i*j*1.e18 + 1.e16*j;
     }} 
   }
  #endif

  return m;
}

/* ********************************************************************* */
char ***Array3D (int nx, int ny, int nz, size_t dsize)
/*! 
 * Allocate memory for a 3-D array of any basic data type.
 *
 * \param [in] nx    number of elements in the 3rd dimension
 * \param [in] ny    number of elements in the 2nd dimension
 * \param [in] nz    number of elements in the 1st dimension
 * \param [in] dsize data-type of the array to be allocated
 * 
 * \return A pointer of type (char ***) to the allocated memory area 
 *         with index range [0...nx-1][0...ny-1][0...nz-1]
 *          
 *********************************************************************** */
{
  int i, j;
  char ***m;

  m = (char ***) malloc ((size_t) nx*sizeof (char **));
  PlutoError (!m, "Allocation failure in Array3D (1)");

  m[0] = (char **) malloc ((size_t) nx*ny*sizeof(char *));
  PlutoError (!m[0],"Allocation failure in Array3D (2)");

  m[0][0] = (char *) malloc ((size_t) nx*ny*nz*dsize);
  PlutoError (!m[0][0],"Allocation failure in Array3D (3)");

/* ---------------------------
       single subscript: i
   --------------------------- */

  for (i = 1; i < nx; i++) m[i] = m[i - 1] + ny;
  
/* ---------------------------
     double subscript:
     
      (i,0)  (0,j) 
      (i,j)  
   --------------------------- */
  
  for (j = 1; j < ny; j++) m[0][j] = m[0][j - 1] + nz*dsize;
  for (i = 1; i < nx; i++) m[i][0] = m[i - 1][0] + ny*nz*dsize;

  for (i = 1; i < nx; i++) {
  for (j = 1; j < ny; j++) {
    m[i][j] = m[i][j - 1] + nz*dsize;
  }}

  for (j = 0; j < ny; j++){
  for (i = 0; i < nx; i++){
    if (m[i][j] == NULL){
      print ("! Allocation failure in Array3D\n");
      QUIT_PLUTO(1);
    }
  }}
  
  g_usedMemory += nx*ny*nz*dsize;

  #if NONZERO_INITIALIZE == YES
   if (dsize==sizeof(double)){
     double ***q;
     int k;
     q = (double ***)m;
     for (i = nx; i--; ){
     for (j = ny; j--; ){
     for (k = nz; k--; ){
       q[i][j][k] = 1.e18*i + 1.e17*j + 1.e16*k + 1.e15;
     }}}
   }
  #endif

  return m;
}
/* ********************************************************************* */
char ****Array4D (int nx, int ny, int nz, int nv, size_t dsize)
/*! 
 * Allocate memory for a 4-D array of any basic data type.
 *
 * \param [in] nx    number of elements in the 4th dimension
 * \param [in] ny    number of elements in the 3rd dimension
 * \param [in] nz    number of elements in the 2nd dimension
 * \param [in] nv    number of elements in the 1st dimension
 * \param [in] dsize data-type of the array to be allocated
 * 
 * \return A pointer of type (char ****) to the allocated memory area 
 *         with index range [0...nx-1][0...ny-1][0...nz-1][0...nv-1]
 *          
 *********************************************************************** */
{
  int i, j, k;
  char ****m;

  m = (char ****) malloc ((size_t) nx*sizeof (char ***));
  PlutoError (!m, "Allocation failure in Array4D (1)");

  m[0] = (char ***) malloc ((size_t) nx*ny*sizeof (char **));
  PlutoError (!m[0], "Allocation failure in Array4D (2)");

  m[0][0] = (char **) malloc ((size_t) nx*ny*nz*sizeof (char *));
  PlutoError (!m[0][0], "Allocation failure in Array4D (3)");

  m[0][0][0] = (char *) malloc ((size_t) nx*ny*nz*nv*dsize);
  PlutoError (!m[0][0][0], "Allocation failure in Array4D (4)");

/* ---------------------------
       single subscript: i
   --------------------------- */
   
  for (i = 1; i < nx; i++) m[i] = m[i - 1] + ny;

/* ---------------------------
     double subscript:
     
      (i,0)  (0,j) 
      (i,j)  
   --------------------------- */

  for (i = 1; i < nx; i++) {
    m[i][0] = m[i - 1][0] + ny*nz;
  }
  for (j = 1; j < ny; j++) {
    m[0][j] = m[0][j - 1] + nz;
  }

  for (i = 1; i < nx; i++) {
  for (j = 1; j < ny; j++) {
    m[i][j] = m[i][j - 1] + nz;
  }}

/* ---------------------------
     triple subscript:
     
     (i,0,0) (0,j,0) (0,0,k)
     (i,j,0) (i,0,k) (0,j,k)
     (i,j,k)
   --------------------------- */

  for (i = 1; i < nx; i++) {
    m[i][0][0] = m[i - 1][0][0] + ny*nz*nv*dsize;
  }
  for (j = 1; j < ny; j++) {
    m[0][j][0] = m[0][j - 1][0] + nz*nv*dsize;
  }
  
  for (k = 1; k < nz; k++) {
    m[0][0][k] = m[0][0][k - 1] + nv*dsize;
  }
  
  
  for (i = 1; i < nx; i++) {
  for (j = 1; j < ny; j++) {
    m[i][j][0] = m[i][j - 1][0] + nz*nv*dsize;
  }}
   
  for (i = 1; i < nx; i++) {
  for (k = 1; k < nz; k++) {
    m[i][0][k] = m[i][0][k - 1] + nv*dsize;
  }}
  
  for (j = 1; j < ny; j++) {
  for (k = 1; k < nz; k++) {
    m[0][j][k] = m[0][j][k - 1] + nv*dsize;
  }}

  for (i = 1; i < nx; i++) {
    for (j = 1; j < ny; j++) {
      for (k = 1; k < nz; k++) {
        m[i][j][k] = m[i][j][k - 1] + nv*dsize;
      }
    }
  }
      
  g_usedMemory += nx*ny*nz*nv*dsize;

  #if NONZERO_INITIALIZE == YES
   if (dsize==sizeof(double)){
     int l;
     double ****q;
     q = (double ****)m;
     for (i = nx; i--; ){
     for (j = ny; j--; ){
     for (k = nz; k--; ){
     for (l = nv; l--; ){
       q[i][j][k][l] = l*1.e18*i + 1.e17*i*j - 1.e16*k*j + 1.e17;
     }}}}
   }
  #endif

  return m;
}
#undef NONZERO_INITIALIZE

/* ********************************************************************* */
double ***ArrayBox(long int nrl, long int nrh, 
                     long int ncl, long int nch,
                     long int ndl, long int ndh)
/*! 
 * Allocate memory for a 3-D array in double precision with given
 * subscript range [low...high] for each direction.
 * Useful for staggered arrays which do not start at [0].
 *
 * \param [in] nrl   lower bound index for the 3rd dimension
 * \param [in] nrh   upper bound index for the 3rd dimension
 * \param [in] ncl   lower bound index for the 2nd dimension
 * \param [in] nch   upper bound index for the 2nd dimension
 * \param [in] ndl   lower bound index for the 1st dimension
 * \param [in] ndh   upper bound index for the 1st dimension
 * 
 * \return A pointer of type (double ***) to the allocated memory area 
 *         with index range [nrl...nhl][ncl...nch][ndl...ndh].
 *          
 *********************************************************************** */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

/* allocate pointers to pointers to rows */

  t=(double ***) malloc((unsigned int)(nrow*sizeof(double**)));
  if (!t) {
    print ("! ArrayBox: allocation failure (1)\n");
    QUIT_PLUTO(1);
  }
  t -= nrl;

/* allocate pointers to rows and set pointers to them */

  t[nrl]=(double **) malloc((unsigned int)(nrow*ncol*sizeof(double*)));
  if (!t[nrl]) {
    print ("! ArrayBox: allocation failure (2)\n");
    QUIT_PLUTO(1);
  }
  t[nrl] -= ncl;

/* allocate rows and set pointers to them */

  t[nrl][ncl]=(double *) malloc((unsigned int)(nrow*ncol*ndep*sizeof(double)));
  if (!t[nrl][ncl]) {
    print ("! ArrayBox: allocation failure (3)\n");
    QUIT_PLUTO(1);
  }
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

/* return pointer to array of pointers to rows */

  return t;
}
/* ********************************************************************* */
void FreeArrayBox(double ***t, long nrl, long ncl, long ndl)
/*!
 * Free memory allocated by the ::ArrayBox function.
 *
 * \param [in] t   pointer to an allocated memory area
 * \param [in] nrl starting index of the array for the 3rd dimension
 * \param [in] ncl starting index of the array for the 2nd dimension
 * \param [in] ndl starting index of the array for the 1st dimension
 *
 *********************************************************************** */
{
  free((char *) (t[nrl][ncl]+ndl));
  free((char *) (t[nrl]+ncl));
  free((char *) (t+nrl));
}


/* ********************************************************************* */
double ***ArrayMap(int nx, int ny, int nz, double *uptr)
/*!
 * Convert a one dimensional array with (nx*ny*nz) elements into
 * a 3D array with index range [0..nx-1][0..ny-1][0..nz-1].
 * This function is similar, conceptually, to Array3D() except that
 * the memory area is already allocated.
 * 
 * \param [in] uptr  pointer to 1D array
 * \param [in]  nx   number of elements in the 3rd dimension
 * \param [in]  ny   number of elements in the 2nd dimensions
 * \param [in]  nz   number of elements in the 1st dimensions
 *
 * \return A pointer to the 3D array.
 *
 *********************************************************************** */
{
  int i,j;
  double ***t;

/* -- allocate pointers to pointers to rows -- */
   
  t = (double ***) malloc((size_t)((nx)*sizeof(double**)));
  if (!t) {
    print ("! ArrayMap: allocation failure (1) \n");
    QUIT_PLUTO (1);
  }

/* -- allocate pointers to rows and set pointers to them -- */

  t[0] = (double **) malloc((size_t)((nx*ny)*sizeof(double*)));
  if (!t[0]) {
    print ("! ArrayMap: allocation failure (2) \n");
    QUIT_PLUTO(1);
  }

/* -- make the 3D array point to uptr -- */

  t[0][0] = uptr;
  if (!t[0][0]) {
    print ("! ArrayMap: allocation failure (3) \n");
    QUIT_PLUTO (1);
  }

  for(j = 1; j < ny; j++) t[0][j]=t[0][j-1] + nz;

  for(i = 1; i < nx; i++) {
    t[i]    = t[i-1] + ny;
    t[i][0] = t[i-1][0] + ny*nz;
    for(j = 1; j < ny; j++) t[i][j] = t[i][j-1] + nz;
  }
     
  return t;
}

/* ********************************************************************* */
unsigned char ***ArrayCharMap(int nx, int ny, int nz, unsigned char *uptr)
/* ********************************************************************* */
{
  int i,j;
  unsigned char ***t;

/* -- allocate pointers to pointers to rows -- */

  t = (unsigned char ***) malloc((size_t)((nx)*sizeof(unsigned char**)));
  if (!t) {
    print ("! ArrayCharMap: allocation failure (1) \n");
    QUIT_PLUTO (1);
  }

/* -- allocate pointers to rows and set pointers to them -- */

  t[0] = (unsigned char **) malloc((size_t)((nx*ny)*sizeof(unsigned char*)));
  if (!t[0]) {
    print ("! ArrayCharMap: allocation failure (2) \n");
    QUIT_PLUTO(1);
  }

/* -- make the 3D array point to uptr -- */

  t[0][0] = uptr;
  if (!t[0][0]) {
    print ("! ArrayCharMap: allocation failure (3) \n");
    QUIT_PLUTO (1);
  }

  for(j = 1; j < ny; j++) t[0][j]=t[0][j-1] + nz;

  for(i = 1; i < nx; i++) {
    t[i]    = t[i-1] + ny;
    t[i][0] = t[i-1][0] + ny*nz;
    for(j = 1; j < ny; j++) t[i][j] = t[i][j-1] + nz;
  }
  
  return t;
}

/* ********************************************************************* */
void FreeArrayMap(double ***t)
/*!
 *  Free memory allocate with ArrayMap()
 *
 *
 *********************************************************************** */
{
  free((char*) t[0]);
  free((char*) t);
}

/* ********************************************************************* */
void FreeArrayCharMap(unsigned char ***t)
{
  free((char*) t[0]);
  free((char*) t);
}

/* ********************************************************************* */
double ***ArrayBoxMap(int nrl, int nrh,
                          int ncl, int nch,
                          int ndl, int ndh, double *uptr)
/*!
 * Convert a one-dimensional array into a 3D array with given
 * subscript range [low...high] for each direction.
 * 
 *
 *********************************************************************** */
{
  int i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

/* -- allocate pointers to pointers to rows -- */

  t = (double ***) malloc((size_t)((nrow)*sizeof(double**)));
  if (!t) {
    print ("! ArrayBoxMap: allocation failure (1) \n");
    QUIT_PLUTO (1);
  }
  t -= nrl;

/* -- allocate pointers to rows and set pointers to them -- */

  t[nrl] = (double **) malloc((size_t)((nrow*ncol)*sizeof(double*)));
  if (!t[nrl]) {
    print ("! ArrayBoxMap: allocation failure (2) \n");
    QUIT_PLUTO (1);
  }
  t[nrl] -= ncl;

  t[nrl][ncl] = uptr;
  if (!t[nrl][ncl]) {
    print ("! ArrayBoxMap: allocation failure (3) \n");
    QUIT_PLUTO (1);
  }
  t[nrl][ncl] -= ndl;

  for(j = ncl+1; j <= nch; j++) t[nrl][j] = t[nrl][j-1] + ndep;
  for(i = nrl+1; i <= nrh; i++) {
    t[i]      = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for(j = ncl+1; j <= nch; j++) t[i][j] = t[i][j-1] + ndep;
  }

  return t;
}

/* ********************************************************************* */
void FreeArrayBoxMap(double ***t, int nrl, int nrh, 
                                      int ncl, int nch,
                                      int ndl,int ndh)
/*! 
 * Free memory allocated with the ::ArrayBoxMap () function
 *
 *********************************************************************** */
{
  free((char*) (t[nrl]+ncl));
  free((char*) (t+nrl));
}

    
