/* //////////////////////////////////////////////////////////////////// */
/*!
  \file 
  \brief Miscellaneous functions to define the properties of the 
         distributed array by setting various parameters via the 
         isz descriptor

  \author  A. Malagoli (University of Chicago)
  \date Jul 17, 1999
*/
/* //////////////////////////////////////////////////////////////////// */
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
int AL_Set_comm(MPI_Comm comm, int isz)
/*! 
 * Set the communicator for a distributed array
 *
 * \param [in]   comm MPI communicator the array is associated with
 * \param [out]  isz  Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the communicator is set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */
{
  int myrank, nproc;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_comm: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nproc);

  s->comm = comm;
  s->size = nproc;
  /*
    Patch for the case in which this node does not
    belong to this communicator 
  */
  if( nproc > 0 ) {
    s->rank = myrank;
  } else { s->rank = MPI_UNDEFINED; }

  return (int) AL_SUCCESS;
}


/* ********************************************************************* */
int AL_Set_dimensions(int ndim, int isz)
/*!
 * Set the dimensions of a distributed array
 *
 * \param [in]   ndim Number of dimensions (1 to 5)
 * \param [out]  isz  Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the dimensions are set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */
{
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_dimensions: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  s->ndim = ndim;

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Set_type(AL_Datatype type, int nelem, int isz)
/*!
 * Set the basic data type of a distributed array. The data types are identical to the MPI datatypes, and they can be defined as AL_<type> or MPI_<type> (e.g. AL_FLOAT or MPI_FLOAT).
 *
 * \param [in]  type  Datatype (AL_Datatype or MPI_Datatype)
 * \param [in]  nelem Number of type 'type' elements in array elements 
 * \param [out] isz   Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the type is set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */
{
  SZ *s;
  AL_Datatype ivector;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_type: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  if( nelem <= 1 ){
    s->type = type;
  } else {
    MPI_Type_contiguous( nelem, type, &ivector);
    MPI_Type_commit(&ivector);
    s->type = ivector;
  }

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Set_global_dim(int *gdims, int isz)
/*!
 * Set the global dimensions of a distributed array
 *
 * \param [in]  gdims Array of integers with global dimensions 
 * \param [out] isz   Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the global dimension are set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */
{
  register int i;
  int ndim;

  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_global_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ s->arrdim[i] = gdims[i] ;}

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Set_local_dim(int *ldims, int isz)
/*!
 * Set the local dimensions of a distributed array
 *
 * \param [in]   ldims Array of integers with local dimensions 
 * \param [out]  isz   Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the local dimension are set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_local_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    s->larrdim[i] = ldims[i];
    s->lbeg[i] = 0;
    s->lend[i] = ldims[i]-1;
  }

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Set_parallel_dim(int *pardims, int isz)
/*!
 * Set the parallel dimensions of a distributed array
 *
 * \param [in]  pardims Array of integers with parallel dimensions [AL_TRUE|AL_FALSE] 
 * \param [out] isz     Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the parallel dimension are set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_parallel_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    s->isparallel[i] = pardims[i];
  }

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Set_periodic_dim(int *periods, int isz)
/*!
 * Set the periodic dimensions of a distributed array
 *
 * \param [in]  periods Array of integers with periodic dimensions [AL_TRUE|AL_FALSE] 
 * \param [out] isz     Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the periodic dimension are set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_periodic_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    if( periods[i] != AL_TRUE && periods[i] != AL_FALSE ){
      printf("Warning: periods has illegal values in AL_Set_periodic_dim\n");
    }
    s->isperiodic[i] = periods[i];
  }

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Set_staggered_dim(int *stagger, int isz)
/*!
 * Set the staggered dimensions of a distributed array
 *
 * \param [in]  stagger Array of integers with staggered dimensions [AL_TRUE|AL_FALSE] 
 * \param [out] isz     Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the staggered dimension are set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_staggered_dim: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    if( stagger[i] != AL_TRUE && stagger[i] != AL_FALSE ){
      printf("Warning: stagger has illegal values in AL_Set_staggered_dim\n");
    }
    s->isstaggered[i] = stagger[i];
  }

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Set_ghosts(int *ghosts, int isz)
/*!
 * Set the ghost points of a distributed array
 *
 * \param [in]  ghost Array of integers with size of ghost points 
                      in each dimension 
 * \param [out] isz   Integer pointer to the input array descriptor
 *
 * \return  AL_SUCCESS if the ghost points are set correctly, 
 *          AL_FAILURE otherwise. 
 *********************************************************************** */ 
{
  register int i;
  int ndim;
  SZ *s;

  /* 
     Check that isz points to an allocated SZ
  */
  if( stack_ptr[isz] == AL_STACK_FREE ){
    printf("AL_Set_ghosts: wrong SZ pointer\n");
    return (int) AL_FAILURE;
  }

  /*
    Get the SZ structure isz is pointing at
  */
  s = sz_stack[isz];

  ndim = s->ndim;
  for(i=0;i<ndim;i++){ 
    s->bg[i] = ghosts[i];
    s->eg[i] = ghosts[i];
  }

  return (int) AL_SUCCESS;
}

