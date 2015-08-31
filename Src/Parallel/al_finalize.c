/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief ArrayLib Finalization routine 

  ArrayLib Finalization routine 

  \author A. Malagoli (University of Chicago)
  \date Jul 17, 1999
*/
/* ///////////////////////////////////////////////////////////////////// */
#include <al_hidden.h>

/* ********************************************************************* */
int AL_Finalize()
/*!
 * Finalize the AL Tool. 
 * It contains a call to MPI_Finalize()
 *********************************************************************** */
{
  int myrank, nproc, errcode;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  /* Synchronize just in case */
  MPI_Barrier(MPI_COMM_WORLD);

  errcode = MPI_Finalize();

#ifdef DEBUG
  printf("AL_Finalize: Called MPI_Finalize: %d\n",errcode);
#endif

  return errcode;
}

