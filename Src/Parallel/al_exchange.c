/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file
  \brief Fill the ghost boundaries

  Fill the ghost boundaries

  \author A. Malagoli (University of Chicago)
  \author A. Mignone (mignone@to.astro.it)

  \date Jun 13, 2007
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
int AL_Exchange(void *vbuf, int sz_ptr)
/*!
 * Fill the ghost boundaries
 *
 * \param [in] vbuf   pointer to buffer
 * \param [in] sz_ptr integer pointer to the distributed array descriptor
 *********************************************************************** */
{
  char *buf;
  register int nd;
  int myrank, nproc;
  int ndim, gp, nleft, nright, tag1, tag2;
  int sendb, recvb;
  MPI_Datatype itype;
  MPI_Comm comm;
  MPI_Status status;
  SZ *s;

  buf = (char *) vbuf;

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
    if( gp > 0 ){
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

/* ********************************************************************* */
int AL_Exchange_periods (void *vbuf, int *periods, int sz_ptr)
/*!
 * Same as AL_Exchange, but exchanges periodic
 * boundaries at physical domain in the dim direction
 * only if periods[dim] = 1.
 * If a dimension is not periodic and periods[dim] = 1
 * nothing changes.
 *
 * \param [in] vbuf     pointer to buffer
 * \param [in] periods  
 * \param [in] sz_ptr   integer pointer to the distributed array descriptor
 *********************************************************************** */
{
  char *buf;
  register int nd;
  int myrank, nproc;
  int ndim, gp, nleft, nright, tag1, tag2;
  int sendb, recvb;
  MPI_Datatype itype;
  MPI_Comm comm;
  MPI_Status status;
  SZ *s;
  int is_beg[3], is_end[3];

  buf = (char *) vbuf;

  /* -- DIAGNOSTICS
        Check that sz_ptr points to an allocated SZ
                                                     -- */
  if( stack_ptr[sz_ptr] == AL_STACK_FREE){
    printf("AL_Decompose: wrong SZ pointer\n");
  }

  s = sz_stack[sz_ptr];

  myrank = s->rank;
  nproc = s->size;
  comm = s->comm;
  ndim = s->ndim;

  AL_Is_boundary (sz_ptr, is_beg, is_end);

  for(nd=0;nd<ndim;nd++){
    gp = s->bg[nd];

    /* If gp=0, do nothing */

    if( gp > 0 ){
      nleft = s->left[nd];
      nright = s->right[nd];

      if (is_beg[nd] && periods[nd] == 0) nleft  = MPI_PROC_NULL;
      if (is_end[nd] && periods[nd] == 0) nright = MPI_PROC_NULL;
      
      itype = s->type_rl[nd];
      tag1 = s->tag1[nd];

      sendb = s->sendb1[nd];
      recvb = s->recvb1[nd];

      MPI_Sendrecv(&buf[sendb], 1, itype, nleft, tag1,
		   &buf[recvb], 1, itype, nright,tag1,
		   comm, &status);

      nleft = s->left[nd];
      nright = s->right[nd];

      if (is_beg[nd] && periods[nd] == 0) nleft  = MPI_PROC_NULL;
      if (is_end[nd] && periods[nd] == 0) nright = MPI_PROC_NULL;

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
  printf("AL_Exchange: filled ghost regions\n"); 
#endif

  return (int) AL_SUCCESS;
}



