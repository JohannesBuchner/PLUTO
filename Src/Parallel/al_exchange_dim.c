/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Fill the ghost boundaries along selected dimensions

  Fill the ghost boundaries along selected dimensions

  \author A. Malagoli (University of Chicago)\n
          A. Mignone 
  \date   Oct 28, 2015
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
int AL_Exchange_dim(char *buf, int *dims, int sz_ptr)
/*!
 * Fill the ghost boundaries along selected dimensions
 *
 * \param [in]  buf     pointer to buffer
 * \param [in]  dims    if dims[i]=0, do not perform the exchange in 
 *                     this dimension (array if int)
 * \param [in]  sz_ptr  integer pointer to the distributed array descriptor
 *********************************************************************** */
{
  register int nd;
  int myrank, nproc;
  int ndim, gp, nleft, nright, tag1, tag2;
  int sendb, recvb;
  MPI_Datatype itype;
  MPI_Comm comm;
  MPI_Status status;
  SZ *s;

  /* DIAGNOSTICS
    Check that sz_ptr points to an allocated SZ
  */
  if( stack_ptr[sz_ptr] == AL_STACK_FREE){
    printf("AL_Decompose: wrong SZ pointer\n");
  }

  s = sz_stack[sz_ptr];

  myrank = s->rank;
  nproc = s->size;
  comm = s->comm;
  ndim = s->ndim;

  for(nd=0;nd<ndim;nd++){
    gp = s->bg[nd];
    
    /* If gp=0, do nothing */
    if( gp > 0 && dims[nd] != 0 ){
      nleft = s->left[nd];
      nright = s->right[nd];
      itype = s->type_rl[nd];
      tag1 = s->tag1[nd];

      sendb = s->sendb1[nd];
      recvb = s->recvb1[nd];

      MPI_Sendrecv(&buf[sendb], 1, itype, nleft, tag1,
		   &buf[recvb], 1, itype, nright,tag1,
		   comm, &status);

      nleft = s->left[nd];
      nright = s->right[nd];
      itype = s->type_lr[nd];
      tag2 = s->tag2[nd];


      sendb = s->sendb2[nd];
      recvb = s->recvb2[nd];

      MPI_Sendrecv(&buf[sendb], 1, itype, nright, tag2,
		   &buf[recvb], 1, itype, nleft,tag2,
		   comm, &status);
    }
  }

  /* DIAGNOSTICS */
#ifdef DEBUG
  if(myrank==0) printf("AL_Exchange: filled ghost regions\n"); 
#endif

  return (int) AL_SUCCESS;
}
