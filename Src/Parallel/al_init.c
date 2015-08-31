/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ArrayLib Initialization routines 

  ArrayLib Initialization routines 

  \authors A. Malagoli (University of Chicago) 

  \date   Jul 17, 1999
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "al_hidden.h"  /*I "al_hidden.h" I*/

static int AL_initialized = AL_FALSE;

/* ********************************************************************* */
int AL_Init(int *argc, char ***argv)
/*!
 * Initialize the AL Tool. It contains a call to MPI_Init(). 
 *
 * \param [in] argc  integer pointer to number of arguments
 * \param [in] argv  pointer to argv list
 ********************************************************************* */
{
  int myrank, nproc, errcode;
  int flag;

  errcode = MPI_Initialized(&flag);

  if( !flag ) 
    errcode = MPI_Init(argc, argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

#ifdef DEBUG
  printf("AL_Init: Called MPI_init from C: %d\n",errcode);
#endif

  /* Initialize the SZ pointers stack */
  AL_Init_stack_();

  /* Synchronize just in case */
  MPI_Barrier(MPI_COMM_WORLD);

  AL_initialized = AL_TRUE;

  return errcode;
}
  
/* ********************************************************************* */
int AL_Initialized()
/*! Test whether or not AL was initialized 
 * 
 * \return This function returns AL_TRUE if AL is initialized, AL_FALSE otherwise.
 *********************************************************************** */
{
  return AL_initialized;
}


