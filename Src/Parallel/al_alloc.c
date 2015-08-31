/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief ArrayLib allocation routine
 
  ArrayLib allocation routine
  
  \authors A. Malagoli (University of Chicago)
  
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
void *AL_Allocate_array(int sz_ptr)
/*!
 * Allocate the local buffer for a distributed array
 *
 * \param [in] sz_ptr pointer to SZ array (integer)
 *********************************************************************** */
{
  char *a;
  int myrank, nproc;
  int size;
  long int buffsize;
  MPI_Comm comm;
  SZ *s;



  /* DIAGNOSTICS
    Check that sz_ptr points to an allocated SZ
  */
  if( stack_ptr[sz_ptr] == AL_STACK_FREE){
    printf("AL_Allocate_array: wrong SZ pointer\n");
  }

  s = sz_stack[sz_ptr];

  myrank = s->rank;
  nproc  = s->size;
  comm   = s->comm;

  buffsize = s->buffsize;

  MPI_Type_size(s->type, &size);

  if( !(a = (char *)AL_CALLOC_((int) buffsize,size)) ){
    printf("[%d] AL_Allocate_array: allocation error\n",myrank);
    return 0;
  }

  /* DIAGNOSTICS */
#ifdef DEBUG
    printf("AL_Allocate_array: allocated %ld bytes\n",buffsize*size);
#endif

  return (void *)a;
}


