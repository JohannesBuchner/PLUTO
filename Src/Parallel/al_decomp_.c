/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Miscellaneous of functions act to find the processors 
         distribution

  Find a "maximally cubic" processors distribution in 1D, 2D and 3D.

  \author A. Malagoli (University of Chicago)
  \date Jul 17, 1999
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "al_hidden.h"  /*I "al_hidden.h" I*/

int AL_Decomp_2d_(int nproc,int npdim,int *ldims, int *gdims);
int AL_Decomp_3d_(int nproc,int npdim,int *ldims, int *gdims);

/* ********************************************************************* */
int AL_Auto_Decomp_(int nproc,int npdim,int *ldims, int *gdims)
/*!
 * Find a "maximally cubic" processors distribution.
 *
 * \param [in]  nproc  number of processors
 * \param [in]  npdim  number of parallel dimensions
 * \param [out] ldims  processor decomposition along directions
 * \param [in]  gdims  global array sizes
 ********************************************************************* */
{
  /*
    1D case. We do not need to do anything, really
  */
  if( npdim == 1 ) {
    ldims[0] = nproc;
    return (int) AL_SUCCESS;
  }

  /*
    2D case. Slightly more complicated
  */
  if( npdim == 2 ) {
    if( AL_Decomp_2d_(nproc, npdim, ldims, gdims) == AL_FAILURE) {
      return (int) AL_FAILURE;    
    } else {
      return (int) AL_SUCCESS;
    }
  }

  /*
    3D case. Even more complicated
  */
  if( npdim == 3 ) {
    if( AL_Decomp_3d_(nproc, npdim, ldims, gdims) == AL_FAILURE) {
      return (int) AL_FAILURE;    
    } else {
      return (int) AL_SUCCESS;
    }
  }


return 0;  
}  

/* ********************************************************************* */
int AL_Decomp_3d_(int nproc,int npdim,int *ldims, int *gdims)
/*!
 * Find a "maximally cubic" processors distribution in 3D.
 *
 * \param [in]  nproc  number of processors
 * \param [in]  npdim  number of parallel dimensions
 * \param [out] ldims  processor decomposition along directions
 * \param [in]  gdims  global array sizes
 ********************************************************************* */
{
  int nx, ny, nz, nnz;
  int pow3, powz, myrank;
  int ndim, minxy, maxxy, nproc2, minp2, maxp2, nprocz;
  int nprocx_old, nprocy_old, nprocz_old;
  int nprocx, nprocy;
  register int ipz;
  int l2dims[AL_MAX_DIM], g2dims[AL_MAX_DIM];

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  ndim = npdim;

  nx = gdims[0];
  ny = gdims[1];
  nz = gdims[2];

  pow3 = floor( log((double) nproc)/log(2.0)+0.5 );

#ifdef DEBUG
  printf("AL_Decomp_3d_: pow3 = %d\n",pow3);
#endif

/*   if( pow3 != 0 ) { */
/*     npw3 = pow(2,pow3); */
/*   } else { */
/*     npw3 = 1; */
/*   } */

/*   if( npw3 != nproc ){ */
  if( !AL_POWEROF2(nproc) ) { 
#ifdef DEBUG
    if( myrank == 0 ) printf("AL_Decompose3d: nproc is not a power of two\n");
#endif
   
/*     MPI_Abort(MPI_COMM_WORLD, -14); */
    return (int) AL_FAILURE;
  }

  powz = rint( log((double) nz)/log(2.0) );
  powz = floor( log((double) nz)/log(2.0) +0.5);

  if( powz != 0 ) {
    nnz = pow(2,powz);
  } else {
    nnz = 1;
  }

  maxxy = nx > ny ? nx : ny;
  minxy = nx < ny ? nx : ny;

  for( ipz=0; ipz<=pow3; ipz++){

#ifdef DEBUG
  printf("AL_Decomp_3d_: In Loop ipz = %d\n",ipz);
#endif

    if( ipz != 0 ) {
      nprocz = pow(2,ipz);
    } else {
      nprocz = 1;
    }

    if( ipz != pow3 ) {
      nproc2 = pow(2,(pow3-ipz));
    } else {
      nproc2 = 1;
    }

    g2dims[0] = nx;
    g2dims[1] = ny;

    if( AL_Decomp_2d_(nproc2, 2, l2dims, g2dims) == AL_FAILURE ) 
      return (int) AL_FAILURE;

    nprocx = l2dims[0];
    nprocy = l2dims[1];

    maxp2 = nprocx > nprocy ? nprocx : nprocy;
    minp2 = nprocx < nprocy ? nprocx : nprocy;

    if( nz >= maxxy ) {
      if( nprocz >= maxp2) break;
    }

    if( nz < maxxy && nz > minxy ) {
      if( nprocz < maxp2 && nprocx > minp2) break;
    }

    if( nz <= minxy ) {
      if( nprocz <= minp2 ){
	nprocx_old = nprocx;
	nprocy_old = nprocy;
	nprocz_old = nprocz;
      } else {
	nprocx = nprocx_old;
	nprocy = nprocy_old;
	nprocz = nprocz_old;
      }
    }
	
  }

#ifdef DEBUG
  printf("AL_decompose3d_: %d %d %d\n", nprocx, nprocy, nprocz);
#endif
  ldims[0] = nprocx;
  ldims[1] = nprocy;
  ldims[2] = nprocz;

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Decomp_2d_(int nproc,int npdim,int *ldims, int *gdims)
/*!
 * Find a "maximally cubic" processors distribution in 2D.
 *
 * \param [in]  nproc  number of processors
 * \param [in]  npdim  number of parallel dimensions
 * \param [out] ldims  processor decomposition along directions
 * \param [in]  gdims  global array sizes
 ********************************************************************* */
{
  int nx, ny;
  int np1, np2, nproc1, nproc2, nproc1_old, nproc2_old;
  int np1_old, np2_old, n1, n2;
  int pow2;
  int nprocx, nprocy;
  int myrank;
  register int ip;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if( nproc == 1 ){
    ldims[0] = 1;
    ldims[1] = 1;
    return (int) AL_SUCCESS;
  }

  /* 
     We assume the we have a power of 2 number of grid points
  */
  nx = gdims[0];
  ny = gdims[1];

  if( nx < ny ){
    n1 = nx;
    n2 = ny;
  } else {
    n1 = ny;
    n2 = nx;
  }
  pow2 = floor( log((double) nproc)/log(2.0)+0.5 );

#ifdef DEBUG
  printf("AL_Decomp2d_ [0]: %d %d\n",pow2, nproc);
#endif
 
/*   if( pow2 != 0 ) { */
/*    npw2 = pow(2,pow2); */
/*   } else { */
/*     npw2 = 1; */
/*   } */

/*   if( npw2 != nproc ){ */
  if( ! AL_POWEROF2(nproc) ) {
#ifdef DEBUG
    if( myrank == 0 ) printf("AL_Decompose2d: nproc is not a power of two\n");
#endif
/*     MPI_Abort(MPI_COMM_WORLD, -14); */
    return (int) AL_FAILURE;
  }

/*   powx = floor( log((double) nx)/log(2.0) +0.5); */

/*   if( powx != 0 ) { */
/*     nnx = pow(2,powx); */
/*   } else { */
/*     nnx = 1; */
/*   } */

/*   if( nnx != nx ){ */
  if( ! AL_POWEROF2(nx) ) {
#ifdef DEBUG
    if( myrank == 0 ) printf("AL_Decompose2d: nx is not a power of two: %d\n",nx);
#endif
/*     MPI_Abort(MPI_COMM_WORLD, -14); */
    return (int) AL_FAILURE;
  }

/*   powy = rint ( log((double) ny)/log(2.0) ); */
/*   powy = floor( log((double) ny)/log(2.0) +0.5); */

/*   if( powy != 0 ) { */
/*     nny = pow(2,powy); */
/*   } else { */
/*     nny = 1; */
/*   } */

/*   if( nny != ny ){ */
  if( ! AL_POWEROF2(ny) ) {
#ifdef DEBUG
   if( myrank == 0 ) printf("AL_Decompose2d: ny is not a power of two: %d\n",ny);
#endif
   
/*     MPI_Abort(MPI_COMM_WORLD, -14); */
    return (int) AL_FAILURE;
  }


  for( ip=0; ip<=pow2; ip++){
    if( ip != 0 ) {
      nproc1 = pow(2,ip);
    } else {
      nproc1 = 1;
    }

    if( ip != pow2 ) {
      nproc2 = pow(2,(pow2-ip));
    } else {
      nproc2 = 1;
    }


    np1 = n1/nproc1;
    np2 = n2/nproc2;

#ifdef DEBUG
    printf("AL_decomp2d_ np1, np2: %d %d | %d %d\n",n1,n2,nproc1,nproc2);
#endif


    nproc1_old = nproc1;
    nproc2_old = nproc2;

    np1_old = np1;
    np2_old = np2;

    if( np1 < np2 ) break;

  }
    
  if( nx < ny ){
    nprocx = nproc1_old;
    nprocy = nproc2_old;
    np1    = np1_old;
    np2    = np2_old;
  } else {
    nprocx = nproc2_old;
    nprocy = nproc1_old;
    np1    = np2_old;
    np2    = np1_old;
  }    

#ifdef DEBUG
  printf("AL_Decomp2d_: %d %d | %d %d\n",nprocx, nprocy, gdims[0], gdims[1]);
#endif

  ldims[0] = nprocx;
  ldims[1] = nprocy;

  return (int) AL_SUCCESS;
}  
