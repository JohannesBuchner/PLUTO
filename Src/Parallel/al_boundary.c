/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Check if the local beginnings or endings of an array
         correspond to the global beginnings or endings.
  
  Check if the local beginnings or endings of an array
  correspond to the global beginnings or endings.

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
int AL_Is_boundary(int sz_ptr, int *is_gbeg, int *is_gend)
/*!
 * This routine is useful when implementing the physical
 * boundary conditions on a problem. It returns two arrays
 * that are set to AL_TRUE or AL_FALSE depending on wether
 * or not the local beginning (ending) address for the array 
 * in each direction is actually a global one.
 *
 * \param [in]  sz_ptr  pointer to a distributed array descriptor (integer)
 * \param [out] is_gbeg int array set to AL_TRUE or AL_FALSE 
                        for global beginning
 * \param [out] is_gend int array set to AL_TRUE or AL_FALSE 
                        for global ending
 *********************************************************************** */
{
  register int i;
  int myrank, nproc,ndims;
  MPI_Comm comm;
  SZ *s;

  /* DIAGNOSTICS
    Check that sz_ptr points to an allocated SZ
  */
  if( stack_ptr[sz_ptr] == AL_STACK_FREE){
    printf("AL_Is_boundary: wrong SZ pointer\n");
  }

  s = sz_stack[sz_ptr];

  myrank = s->rank;
  nproc  = s->size;
  comm   = s->comm;
  ndims  = s->ndim;

  for(i=0;i<ndims;i++){
    if( s->beg[i]  == s->bg[i] ){
      is_gbeg[i] = AL_TRUE;
    } else {
      is_gbeg[i] = AL_FALSE;
    }
    if( s->end[i] == s->arrdim[i]+s->bg[i]-1 ){
      is_gend[i] = AL_TRUE;
    } else {
      is_gend[i] = AL_FALSE;
    }
  }

  /* DIAGNOSTICS */
#ifdef DEBUG
    printf("AL_Is_boundary: completed\n");
#endif

  return (int) AL_SUCCESS;
}
