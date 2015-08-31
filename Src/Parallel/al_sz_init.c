/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Allocate and initialize a distributed array descriptor

  \author A. Malagoli (University of Chicago)
  \date Jul 28, 1999 
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
int AL_Sz_init(MPI_Comm comm, int *sz_ptr)
/*!
 * Allocate and initialize a distributed array descriptor
 *
 * \param [in]  comm    MPI communicator the array is associated with
 * \param [out] sz_ptr  Integer pointer to the array descriptor
 *
 * \return  AL_SUCCESS if the array descriptor is correct initialized, 
 *          AL_FAILURE otherwise. 
 *************************************************************************/
{
  int myrank, nproc;
  register int i;

  /*
    If this process does not belong to the communicator,
    or if the communicator is MPI_COMM_NULL, return a
    AL_UNDEFINED pointer
  */
  if( comm == MPI_COMM_NULL ) {
    *sz_ptr = AL_UNDEFINED;
    return AL_FAILURE;
  }

  /*
    Allocate the SZ structure 
  */
  *sz_ptr = AL_Allocate_sz_();

  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nproc);

  /*
    Set the initial entries
  */
  sz_stack[*sz_ptr]->compiled = AL_FALSE; /* The array is not compiled yet */
  sz_stack[*sz_ptr]->comm = comm;
  sz_stack[*sz_ptr]->rank = myrank;
  sz_stack[*sz_ptr]->size = nproc;

  /* 
     Set some defaults 
  */
  sz_stack[*sz_ptr]->type = MPI_DOUBLE;
  for(i=0;i<AL_MAX_DIM;i++){
    sz_stack[*sz_ptr]->isparallel[i] = AL_TRUE;
    sz_stack[*sz_ptr]->isperiodic[i] = AL_TRUE;
    sz_stack[*sz_ptr]->isstaggered[i] = AL_FALSE;
    sz_stack[*sz_ptr]->larrdim_gp[i]=1;
    sz_stack[*sz_ptr]->larrdim[i]=1;
    sz_stack[*sz_ptr]->arrdim[i]=1;
    sz_stack[*sz_ptr]->beg[i]=0;
    sz_stack[*sz_ptr]->end[i]=0;
    sz_stack[*sz_ptr]->lbeg[i]=0;
    sz_stack[*sz_ptr]->lend[i]=0;
    sz_stack[*sz_ptr]->bg[i]=0;
    sz_stack[*sz_ptr]->eg[i]=0;
    sz_stack[*sz_ptr]->offset[i]=1;
    sz_stack[*sz_ptr]->stride[i]=1;
  }

  sz_stack[*sz_ptr]->begs = NULL;
  sz_stack[*sz_ptr]->ends = NULL;

  return (int) AL_SUCCESS;
}

