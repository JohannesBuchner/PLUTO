/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Miscellaneous functions to obtain informations about the 
         properties of the distributed array  
  
  \author A. Malagoli (University of Chicago)
  \date Jul 17, 1999
*/  
/* ///////////////////////////////////////////////////////////////////// */
#include "al_hidden.h"  /*I "al_hidden.h" I*/

/* 
   The SZ structure stack is defined and maintained
   in al_szptr_.c
   Here we include an external reference to it in
   order to be able to make internal references to it.
*/
extern SZ *sz_stack[AL_MAX_ARRAYS];
extern int stack_ptr[AL_MAX_ARRAYS];

/* ********************************************************************* */
int AL_Get_size(int isz, int *size)
/*!
 * Get the size of the communicator associated with a distributed array
 *
 * \param [in] isz   Integer pointer to the input array descriptor (input)
 * \param [in] size  Integer pointer to size (output)
 *********************************************************************** */
{
  int nproc;
  MPI_Comm comm;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_comm: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  comm = s->comm;

  MPI_Comm_size(comm, &nproc);

  *size = nproc;

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_rank(int isz, int *rank)
/*!
 * Get the rank of this node in the communicator associated with a 
 * distributed array
 *
 * \param [in] isz   Integer pointer to the input array descriptor 
 * \param [out] rank Integer pointer to rank 
 *
 * \return  If the current node does not belong to the communicator, 
 *          AL_UNDEFINED (same as MPI_UNDEFINED) is returned.
 *          Otherwise AL_SUCCESS is returned.
 *********************************************************************** */
{

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_comm: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  *rank = sz_stack[isz]->rank;

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Get_comm(int isz, MPI_Comm *comm)
/*!
 * Get the communicator for a distributed array
 *
 * \param [in]  isz   Integer pointer to the input array descriptor 
 * \param [out] comm  Pointer to communicator 
 *********************************************************************** */
{
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_comm: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  *comm = s->comm;

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_cart_comm(int isz, MPI_Comm *cart_comm)
/*!
 * Get the cartesian communicator for a distributed array
 *
 * \param [in]  isz         Integer pointer to the input array descriptor 
 * \param [out] cart_comm   Pointer to cartesian communicator
 *********************************************************************** */
{
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_comm: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  *cart_comm = s->cart_comm;

  return (int) AL_SUCCESS;
}



/* ********************************************************************* */
int AL_Get_dimensions(int isz, int *ndim)
/*!
 * Get the dimensions of the distributed array
 *
 * \param [in]  isz   Integer pointer to the input array descriptor 
 * \param [out] ndim  Pointer to integer number of dimensions 
 *********************************************************************** */
{
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_dimensions: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  *ndim = s->ndim;

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Get_type(int isz, AL_Datatype *type)
/*!
 * Get the basic data type of a distributed array
 *
 * \param [in]  isz   Integer pointer to the input array descriptor
 * \param [out] type  Pointer to datatype 
 *********************************************************************** */
{
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_type: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  *type = s->type;

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Get_type_size(int isz, int *tsize)
/*!
 * Get size of the basic data type of a distributed array
 *
 * \param [in]  isz    Integer pointer to the input array descriptor
 * \param [out] tsize  Integer pointer to size of datatype. 
 *********************************************************************** */
{
  int itsize;
  AL_Datatype type;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_type_size: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  type = s->type;

  MPI_Type_size( type, &itsize);

  *tsize = itsize;
  
  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_buffsize(int sz_ptr, int *buffsize)
/*!
 * Get the local buffer size for a distributed array, including 
 * the number of elements of the ghsot regions.
 *
 * \param [in]  sz_ptr    Integer pointer to the input array descriptor 
 * \param [out] buffsize  Buffer size [in number of elements] 
 *********************************************************************** */
{
  int isz;
  SZ *s;

  isz = sz_ptr;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_buffsize: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  *buffsize = s->buffsize;

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_global_dim(int isz, int *gdims)
/*!
 * Get the global dimensions of the distributed array
 *
 * \param [in]  isz    Integer pointer to the input array descriptor
 * \param [out] gdims  Array of integers with global dimensions 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_global_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ gdims[i] = s->arrdim[i];}

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_local_dim(int isz, int *ldims)
/*!
 * Get the local dimensions of a distributed array WITHOUT GHOST POINTS
 *
 * \param [in]  isz     Integer pointer to the input array descriptor
 * \param [out] ldims   Array of integers with local dimensions 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_local_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ ldims[i] = s->larrdim[i];}

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Get_local_dim_gp(int sz_ptr, int *ldims_gp)
/*!
 * Get the local dimensions of a distributed array WITH GHOST POINTS
 *
 * \param [in]   sz_ptr    Integer pointer to the input array descriptor
 * \param [out]  ldims_gp  Array of integers with local dimensions 
 *                         (including ghost points)
 *********************************************************************** */
{
  register int i;
  int isz;
  int ndim;
  SZ *s;

  isz = sz_ptr;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_local_dim_gp: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ ldims_gp[i] = s->larrdim_gp[i];}

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_offsets(int sz_ptr, int *offset)
/*!
 * Get the offsets for the for loop index calculation
 *
 * NOTE: 
 * The typical multidimensional C for loop will look like, for example:
 *
 *  ... ;
 *
 *  AL_Get_offsets(sz_ptr, offset);
 *
 *  for(k=kbeg; k<=kend; k++){
 *    koff = k*offset[2];
 *      for(j=jbeg; j<=jend; j++){
 *         joff = j*offset[1];
 *         for(i=ibeg; i<=iend; i++){
 *            ioff = i*offset[0]+joff+koff;
 *            a[ioff] = ... ;
 *         }
 *      }
 *  }
 *  ... ;
 * 
 *  Where a[ioff] is a generic 3D array.
 *
 * \param [in]   sz_ptr   Integer pointer to the input array descriptor
 * \param [out]  offset   Array of integers with offsets
 *********************************************************************** */
{
  register int i;
  int isz;
  int ndim;
  SZ *s;

  isz = sz_ptr;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_offsets: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    offset[i] = s->offset[i];
  }

  return (int) AL_SUCCESS;
}



/* ********************************************************************* */
int AL_Get_stride(int sz_ptr, int *stride)
/*!
 * Get the strides for the for loop index calculation
 * NOTE: 
 * The typical multidimensional C for loop will look like, for example:
 *
 * ... ;
 *
 * AL_Get_stride(sz_ptr, stride);
 *
 * for(k=kbeg; k<=kend; k++){
 *   koff = k*stride[2];
 *     for(j=jbeg; j<=jend; j++){
 *        joff = j*stride[1];
 *        for(i=ibeg; i<=iend; i++){
 *           ioff = i*stride[0]+joff+koff;
 *           a[ioff] = ... ;
 *        }
 *     }
 * }
 * ... ;
 * 
 * Where a[ioff] is a generic 3D array.
 *
 * \param [in]  sz_ptr  Integer pointer to the input array descriptor 
 * \param [out] stride  Array of strides 
 *********************************************************************** */
{
  register int i;
  int isz;
  int ndim;
  SZ *s;

  isz = sz_ptr;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_offsets: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    stride[i] = s->stride[i];
  }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_cart_sizes(int isz, int *cart_sizes)
/*!
 * Get the sizes of the cartesian communicator
 *
 * \param [in]  isz          Integer pointer to the input array descriptor
 * \param [out] cart_sizes   Array of integers with cartesian sizes
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_cart_sizes: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    cart_sizes[i] = s->lsize[i];
  }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_parallel_dim(int isz, int *pardims)
/*!
 * Get the parallel dimensions of a distributed array
 *
 * \param [in]   isz       Integer pointer to the input array descriptor
 * \param [out]  pardims   Array of integers with parallel dimensions 
 *                         [AL_TRUE|AL_FALSE] 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_parallel_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    pardims[i] = s->isparallel[i];
  }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_periodic_dim(int isz, int *periods)
/*!
 * Get the periodic dimensions of the distributed arrays
 *
 * \param [in]  isz      Integer pointer to the input array descriptor
 * \param [out] periods  Array of integers with periodic dimensions 
 *                       [AL_TRUE|AL_FALSE] 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_periodic_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    periods[i] = s->isperiodic[i];
  }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_staggered_dim(int isz, int *stagger)
/*!
 * Get the staggered dimensions of the distributed array
 *
 * \param [in]  isz       Integer pointer to the input array descriptor
 * \param [out] stagger   Array of integers with staggered dimensions 
 *                        [AL_TRUE|AL_FALSE] 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_staggered_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    stagger[i] = s->isstaggered[i];
  }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_ghosts(int isz, int *ghosts)
/*!
 * Get the ghost points of the distributed array
 *
 * \param [in]  isz    Integer pointer to the input array descriptor
 * \param [out] ghost  Array of integers with size of ghost points 
 *                     in each dimension 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_ghosts: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    ghosts[i] = s->bg[i];
  }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_lbounds(int isz, int *lbeg, int *lend, int *gp, int style)
/*!
 * Get the local indexes for the local portion of a distributed array
 * 
 * \param [in]   isz   Integer pointer to the input array descriptor
 * \param [out]  lbeg  Array of ndim integers containing the start points
 * \param [out]  lend  Array of ndim integers containing the end points
 * \param [out]  gp    Array of ndim integers containing the ghost points
 * \param [in]   style Index style: AL_C_INDEXES or AL_FORTRAN_INDEXES
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_lbounds: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    lbeg[i] = s->lbeg[i];
    lend[i] = s->lend[i];
    gp[i] = s->bg[i];
  }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_bounds(int isz, int *beg, int *end, int *gp, int style)
/*!
 * Get the global indexes for the local portion of a distributed array 
 *
 * \param [in]   isz    Integer pointer to the input array descriptor
 * \param [out]  beg    Array of ndim integers containing the start points
 * \param [out]  end    Array of ndim integers containing the end points
 * \param [out]  gp     Array of ndim integers containing the ghost points
 * \param [in]   style  Index style: AL_C_INDEXES or AL_FORTRAN_INDEXES
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_gbounds: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    beg[i] = s->beg[i];
    end[i] = s->end[i];
    gp[i] = s->bg[i];
  }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Get_gbounds(int isz, int *gbeg, int *gend, int *gp, int style)
/*!
 * Get the global bounds of the global array
 *
 * \param [in]   isz     Integer pointer to the input array descriptor
 * \param [out]  gbeg    Array of ndim integers containing the start points
 * \param [out]  gend    Array of ndim integers containing the end points
 * \param [out]  gp      Array of ndim integers containing the ghost points
 * \param [in]   style   Index style: AL_C_INDEXES or AL_FORTRAN_INDEXES
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Get_gbounds: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];
  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    gbeg[i] = s->bg[i];
    gend[i] = s->arrdim[i]+s->bg[i]-1;
    gp[i] = s->bg[i];
  }

  return (int) AL_SUCCESS;
}
