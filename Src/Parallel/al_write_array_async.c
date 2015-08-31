/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief ArrayLib routines for asynchronous MPI-IO
 
  ArrayLib routines for asynchronous MPI-IO
  
  \authors G. Muscianisi (g.muscianisi@cineca.it)
  
  \date Feb 28, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "al_hidden.h"  /*I "al_hidden.h" I*/
#include <string.h>

/* 
   The SZ structure stack is defined and maintained
   in al_szptr_.c
   Here we include an external reference to it in
   order to be able to make internal references to it.
*/
extern SZ *sz_stack[AL_MAX_ARRAYS];
extern int stack_ptr[AL_MAX_ARRAYS];

/* ********************************************************************* */
int AL_Write_array_begin(void *va, int sz_ptr, int *output_stag, int *output_dump, int output_nvar)
/*!
 * Write a distributed array to a parallel file by using 
 * asynchronous MPI-IO 
 *
 * \param [in] buffer        pointer to the buffer to write
 * \param [in] sz_ptr        integer pointer to the distributed 
 *                           array descriptor
 * \param [in] output_stag   vector sets to -1 for centred variables,
 *                           and sets to 0,1,2 for staggered field in 
 *                           the x,y,z direction
 * \param [in] output_dump   vector sets to 1 if the variable has to 
 *                           be dumped, 0 in the contrary case
 * \param [in] output_nvar   total number of variables in PLUTO simulation
 *********************************************************************** */
{
  char *a;
  register int i;
  SZ *s;

  MPI_File ifp;

  int errcode;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  a = (char *) va;

  s = sz_stack[sz_ptr];
  ifp = s->ifp;

  /* Definition of a "global" filetype given by nvar times gsubarr for the MPI_File_set_view */ 
  MPI_Datatype filetype_gsub_arr; 
  int count; /* number of variable to be dumped */ 
  count=0;
  for (i=0; i<output_nvar; i++){ 
    if (!(output_dump[i]== 0)) count = count+1;
  } 
  errcode = MPI_Type_contiguous(count, s->gsubarr, &filetype_gsub_arr);
  errcode = MPI_Type_commit(&filetype_gsub_arr); 
  
  /* Definition of a "global" datatype for the MPI_File_write_all_begin */ 
  MPI_Datatype dtype_lsub_arr; 

#ifdef DEBUG
  if(myrank==0) printf("count=%d, output_nvar=%d\n",count,output_nvar);
#endif

  int *block_leng, *displ;
  int k=0;

  block_leng= (int*) malloc(count*sizeof(int));
  displ=(int*) malloc(count*sizeof(int));

  for (i=0; i<output_nvar; i++){
    if (!(output_dump[i]==0)){  /* cioè devo scrivere la variabile */ 
     block_leng[k] = 1;
     displ[k] = i;
     k=k+1;
    }
  }
  
  errcode = MPI_Type_indexed(count,block_leng, displ, s->lsubarr, &dtype_lsub_arr);
  errcode = MPI_Type_commit(&dtype_lsub_arr); 
/*
  MPI_Type_get_extent(dtype_lsub_arr,lb,extent);
  printf("myid, %d, lb %d, extent %d of dtype_lsub_arr\n",myrank,lb,extent);
*/
/* DIAGNOSTICS */
#ifdef DEBUG
  if (myrank==0) {
    for (i=0; i<output_nvar; i++) printf("myid, %d, count %d, output_nvar %d, i %d, output_dump[i] %d, block_leng[i] %d, displ[i] %d\n", myrank, count, output_nvar, i, output_dump[i], block_leng[i], displ[i]); 
  }
#endif

  free(block_leng);
  free(displ);
  MPI_Barrier(s->comm);

  /* Setting of the new view. Each procs see all the file at the beginning of the writing */
  errcode = MPI_File_set_view(ifp, 0, MPI_BYTE, filetype_gsub_arr,
                    "native", MPI_INFO_NULL);
  MPI_Type_free(&filetype_gsub_arr);

/* DIAGNOSTICS */
#ifdef DEBUG
  int len;
  char es[MPI_MAX_ERROR_STRING];
  MPI_Error_string(errcode, es, &len);
  printf("myid %d, Errcode from MPI_File_set_view: %d | %s\n", myrank,  errcode, es);
#endif

  errcode = MPI_File_write_all_begin(ifp, va, 1,dtype_lsub_arr);
  MPI_Type_free(&dtype_lsub_arr);
    
/* DIAGNOSTICS */
#ifdef DEBUG
  printf("myid %d\n", myrank);
  MPI_Error_string(errcode, es, &len);
  printf("myid %d, Errcode from MPI_File_write_all_begin: %d | %s\n", myrank, errcode, es);
#endif

  return (int) AL_SUCCESS;
}

/* ********************************************************************* */
int AL_Write_array_end(void *va, int sz_ptr)
/*!
 * Completition of writing of a distributed array to a 
 * parallel file by using asynchronous MPI-IO 
 *
 * \param [in] buffer  pointer to the buffer to write
 * \param [in] sz_ptr  integer pointer to the distributed array descriptor
 *********************************************************************** */
{
  char *a;
  SZ *s;

  MPI_Status status;

  MPI_File ifp;

  int errcode;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  a = (char *) va;

  s = sz_stack[sz_ptr];
  ifp = s->ifp;

  errcode=MPI_File_write_all_end(ifp, a, &status);

  /* DIAGNOSTICS */
#ifdef DEBUG
  int len;
  char es[128];
  if( errcode ){
    MPI_Error_string(errcode, es, &len);
    printf("Errcode from MPI_File_write_all_end: %d | %s\n", errcode, es);
  }
#endif

  MPI_File_sync(ifp);
   
  return (int) AL_SUCCESS;
}
