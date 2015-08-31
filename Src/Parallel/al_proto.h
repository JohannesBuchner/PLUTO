/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ArrayLib function prototypes header file.

  Contains the function prototypes used in the ArrayLib routines.

  \author A. Malagoli (University of Chicago)
  \author A. Mignone (mignone@ph.unito.it)
  \author G. Muscianisi (g.muscianisi@cineca.it)
 
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#ifndef __AL_PROTO
#define __AL_PROTO


#ifdef __cplusplus
extern "C" {
#endif
/* External prototypes */
extern int AL_Init(int *, char ***);
extern int AL_Finalize();
extern int AL_Initialized();
extern int AL_Sz_init(MPI_Comm, int *);
extern int AL_Free(int);
extern int AL_Sz_free(int);
extern int AL_Valid_ptr(int);

extern int AL_Set_comm(MPI_Comm, int);
extern int AL_Set_dimensions(int, int);
extern int AL_Set_type(AL_Datatype, int, int);
extern int AL_Set_global_dim(int *, int);
extern int AL_Set_local_dim(int *, int);
extern int AL_Set_parallel_dim(int *, int);
extern int AL_Set_periodic_dim(int *, int);
extern int AL_Set_staggered_dim(int *, int);
extern int AL_Set_ghosts(int *, int);

extern int AL_Get_size(int, int *);
extern int AL_Get_comm(int, MPI_Comm *);
extern int AL_Get_cart_comm(int, MPI_Comm *);
extern int AL_Get_dimensions(int, int *);
extern int AL_Get_type(int ,AL_Datatype *);
extern int AL_Get_buffsize(int ,int *);
extern int AL_Get_global_dim(int, int *);
extern int AL_Get_local_dim(int, int *);
extern int AL_Get_local_dim_gp(int, int *);
extern int AL_Get_parallel_dim(int , int *);
extern int AL_Get_periodic_dim(int , int * );
extern int AL_Get_staggered_dim(int , int *);
extern int AL_Get_ghosts(int, int *);
extern int AL_Get_offsets(int, int *);
extern int AL_Get_ghosts (int, int *);
extern int AL_Get_lbounds(int, int *, int *, int *, int);
extern int AL_Get_gbounds(int, int *, int *, int *, int);
extern int AL_Get_bounds(int, int *, int *, int *, int);

extern int AL_Is_boundary(int , int *, int *);
extern int AL_Get_stride(int, int *);

extern int AL_Decompose( int, int *, int );
extern int AL_Type_create_subarray(int, int *, int *, int *, int, MPI_Datatype, MPI_Datatype *);

extern void *AL_Allocate_array(int);
extern int AL_Exchange( void *, int);
extern int AL_Exchange_dim(char *, int *, int);
extern int AL_Exchange_periods (void *vbuf, int *periods, int sz_ptr);

extern int AL_File_open(char *, int);
extern long long AL_Get_offset(int);
extern int AL_Set_offset(int, long long);

extern int AL_Write_header(void *, int, AL_Datatype, int);
extern int AL_File_close( int);
extern int AL_Write_common(void *, int, AL_Datatype, int);
extern int AL_Read_common(void *, int, AL_Datatype, int);
extern int AL_Write_array(void *, int, int);
extern int AL_Read_array(void *, int, int);

extern int AL_Write_array_begin(void *, int , int *, int *, int);
extern int AL_Write_array_end(void *, int);
/*
extern int AL_Write_array_begin(void *, int, int, int *, int);
extern int AL_Write_array_end(void *, int); 
 */
/* Internals prototypes */
extern int AL_Init_stack_();
extern int AL_Allocate_sz_();
extern int AL_Deallocate_sz_(int);
extern int AL_Auto_Decomp_(int, int, int *, int *);
extern int AL_Sort_(int, int *, int *);

#ifdef __cplusplus
}
#endif

#endif /* End ifdef __AL_PROTO */
