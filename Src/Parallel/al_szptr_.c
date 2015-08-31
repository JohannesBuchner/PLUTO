/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Miscellaneous of internal routines to convert sz pointers
         to integer pointers. 
  
  \author A. Malagoli (University of Chicago)
  \date Jul 17, 1999
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "al_hidden.h"  /*I "al_hidden.h" I*/

/* In order to make life easier with the Fortran interface,
   we do not pass the pointers to the SZ structures directly,
   but we rather allocate statically an array of pointers to
   SZ structures, and pass an integer index that points to
   an entry in this array */
SZ *sz_stack[AL_MAX_ARRAYS];

/* The array of SZ pointers is treated as a stack, and 
   the routines in this file take care of allocating
   and deallocating pointers to this array. So, we maintain
   an array of integers for the bookkeeping, in order to
   define which entries are free, and which entries are
   allocated */
int stack_ptr[AL_MAX_ARRAYS];

/* More bookkeeping variables */
static int stack_top;  /* First unallocated SZ pointer */
static int stack_used; /* Number of used entries in the stack */

/* ********************************************************************* */
int AL_Init_stack_()
/*!
  * Initialize the stack of SZ pointers
  * This routine is called internally by AL_Init.
  ********************************************************************** */
{
  register int i;
  int myrank;

  stack_top  = 0;
  stack_used = 0;

  for( i=0; i<AL_MAX_ARRAYS;i++){ stack_ptr[i] = AL_STACK_FREE ;}

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
#ifdef DEBUG
  printf("AL_Init_stack_: SZ stack initialized\n");
#endif

  return 0;
}


/* ********************************************************************* */
int AL_Allocate_sz_()
/*!
 * Return an integer pointer to an SZ stack entry.
 *********************************************************************** */
{
  register int i;
  int sz_ptr; /* The returned pointer to the SZ array */

  sz_ptr = -1;

  /*
    First, search the stack for the first 
     available pointer 
  */
  for( i=0; i<AL_MAX_ARRAYS; i++){
    if( stack_ptr[i] == AL_STACK_FREE ){
      sz_ptr = i;
      stack_ptr[i] = AL_STACK_USED;
      break;
    }
  }

  /* 
     If the allocation did not fail, then proceed to
     allocate a SZ structure 
  */
  if( sz_ptr != -1 ){
    if( !(sz_stack[i] = (SZ *)malloc(sizeof(SZ))) ){
      printf("AL_Allocate_sz_: Failed to allocate SZ\n");
    }
  }

  /*
    If we were successful, then set sz.compiled to AL_FALSE,
    indicating the need for intialization.
  */
  sz_stack[sz_ptr]->compiled = AL_FALSE;

  return sz_ptr;

} 


/* ********************************************************************* */
int AL_Deallocate_sz_(int sz_ptr)
/*!
 * Deallocate an integer pointer to an SZ stack entry.
 *
 * \param [in] sz_ptr  Integer pointer to an entry in the SZ stack
 *********************************************************************** */
{
  register int i;
  int ndim;
  SZ *s;

  if(sz_ptr<0) return AL_SUCCESS;
  if(stack_ptr[sz_ptr] == AL_STACK_FREE) return AL_SUCCESS;
  
  s = sz_stack[sz_ptr];
  ndim = s->ndim;

  /* 
     Free the MPI data types
  */
  if( s->compiled == AL_TRUE ){
    for( i=0; i<ndim; i++){
      MPI_Type_free(&(s->strided[i]));
      MPI_Type_free(&(s->type_lr[i]));
      MPI_Type_free(&(s->type_rl[i]));
      MPI_Comm_free(&(s->oned_comm[i]));
    }

    MPI_Type_free(&(s->gsubarr));
    MPI_Type_free(&(s->lsubarr));
    MPI_Comm_free(&(s->cart_comm));
  }

  if( (s->begs != NULL) ) free(s->begs);

  /*
    Begin by dellocating the SZ structure 
  */
  free(sz_stack[sz_ptr]); 
  
  /*
    Now we can reset the integer pointer
  */
  stack_ptr[sz_ptr] = AL_STACK_FREE;
  return (int) AL_SUCCESS;
} 

/* ********************************************************************* */
int AL_Valid_ptr(int sz_ptr)
/*!
 * Return AL_TRUE if the input pointer points
 * to an allocated distributed array. 
 * Return AL_FALSE otherwise.
 *
 * \param [in] sz_ptr  Integer pointer to an entry in the SZ stack
 *********************************************************************** */
{
  if(stack_ptr[sz_ptr]!=AL_STACK_USED ) return AL_FALSE;
  return sz_stack[sz_ptr]->compiled;
}
