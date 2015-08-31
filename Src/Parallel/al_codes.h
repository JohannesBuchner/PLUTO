/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ArrayLib codes header file.

  Contains basic macro definitions used in ArrayLib routines.

  \author A. Malagoli (University of Chicago)
  \author G. Muscianisi (g.muscianisi@cineca.it)

  \date   Aug 29, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#ifndef __AL_CODES
#define __AL_CODES

/* Maximum dimension of supported arrays */
#define AL_MAX_DIM     ((int)5)

/* Maximum number of supported arrays */
#define AL_MAX_ARRAYS  ((int)100)

/* Stack indicator values for stack_ptr (in al_szptr_.c) */
#define AL_STACK_FREE  ((int)0)
#define AL_STACK_USED  ((int)1)

/* Success/Fail codes */
#define AL_TRUE  ((int)1)
#define AL_FALSE ((int)0)

/* Success/Fail codes */
#define AL_SUCCESS  ((AL_Const)0)
#define AL_FAILURE  ((AL_Const)1)

/* Codes for array indexing styles: C or Fortran indexes */
#define AL_C_INDEXES ((int) 7)
#define AL_FORTRAN_INDEXES ((int) 8)

/* Codes for the decomposition mode */
#define AL_AUTO_DECOMP ((AL_Const)  9)  /* Internal decomposition */
#define AL_MPI_DECOMP  ((AL_Const) 10)  /* MPI decomposition */
#define AL_USER_DECOMP ((AL_Const) 11)  /* User defined decomposition */

/* Codes for al_subarray routine (derived from the romio function) */
#define AL_ORDER_C            ((int) 56)
#define AL_ORDER_FORTRAN      ((int) 57)
#define AL_DISTRIBUTE_BLOCK   ((int) 121)
#define AL_DISTRIBUTE_CYCLIC  ((int) 122)
#define AL_DISTRIBUTE_NONE    ((int) 123)


/*
  We duplicate here some common MPI parameters, just in
   case we'll decide to write a non-MPI version in the
   future.
*/

#define AL_CHAR             MPI_CHAR
#define AL_UNSIGNED_CHAR    MPI_UNSIGNED_CHAR
#define AL_BYTE             MPI_BYTE
#define AL_SHORT            MPI_SHORT

#define AL_UNSIGNED         MPI_UNSIGNED      
#define AL_LONG             MPI_LONG          
#define AL_UNSIGNED_LONG    MPI_UNSIGNED_LONG  
#define AL_FLOAT            MPI_FLOAT          
#define AL_DOUBLE           MPI_DOUBLE         
#define AL_LONG_DOUBLE      MPI_LONG_DOUBLE    
#define AL_LONG_LONG_INT    MPI_LONG_LONG_INT  

#define AL_PACKED           MPI_PACKED         
#define AL_LB               MPI_LB             
#define AL_UB               MPI_UB             

/* 
   The layouts for the types MPI_DOUBLE_INT etc are simply
   struct { 
       double var;
       int    loc;
   }
   This is documented in the man pages on the various datatypes.   
 */
#define AL_FLOAT_INT        MPI_FLOAT_INT      
#define AL_DOUBLE_INT       MPI_DOUBLE_INT     
#define AL_LONG_INT         MPI_LONG_INT       
#define AL_SHORT_INT        MPI_SHORT_INT      
#define AL__2INT            MPI_2INT           
#define AL_LONG_DOUBLE_INT  MPI_LONG_DOUBLE_INT

/* Communicators */
#define AL_COMM_WORLD       MPI_COMM_WORLD 
#define AL_COMM_SELF        MPI_COMM_SELF 

/* Groups */
#define AL_GROUP_EMPTY      MPI_GROUP_EMPTY

/* Collective operations */

#define AL_MAX              MPI_MAX    
#define AL_MIN              MPI_MIN    
#define AL_SUM              MPI_SUM    
#define AL_PROD             MPI_PROD   
#define AL_LAND             MPI_LAND   
#define AL_BAND             MPI_BAND   
#define AL_LOR              MPI_LOR    
#define AL_BOR              MPI_BOR    
#define AL_LXOR             MPI_LXOR   
#define AL_BXOR             MPI_BXOR   
#define AL_MINLOC           MPI_MINLOC 
#define AL_MAXLOC           MPI_MAXLOC 


/* Define some null objects */
#define AL_COMM_NULL        MPI_COMM_NULL      
#define AL_OP_NULL          MPI_OP_NULL        
#define AL_GROUP_NULL       MPI_GROUP_NULL     
#define AL_DATATYPE_NULL    MPI_DATATYPE_NULL  
#define AL_REQUEST_NULL     MPI_REQUEST_NULL   
#define AL_ERRHANDLER_NULL  MPI_ERRHANDLER_NULL 

/* These are only guesses; make sure you change them in mpif.h as well */
#define AL_MAX_PROCESSOR_NAME    MPI_MAX_PROCESSOR_NAME
#define AL_MAX_ERROR_STRING      MPI_MAX_ERROR_STRING
#define AL_MAX_NAME_STRING       MPI_MAX_NAME_STRING     /*  How long a name do you need ? */

/* Pre-defined constants */
#define AL_UNDEFINED         MPI_UNDEFINED 
#define AL_UNDEFINED_RANK    MPI_UNDEFINED_RANK
#define AL_KEYVAL_INVALID    MPI_KEYVAL_INVALID

#define AL_PROC_NULL         MPI_PROC_NULL
#define AL_ANYSOURCE         MPI_ANY_SOURCE 
#define AL_ANY_TAG           MPI_ANY_TAG

#endif /* End ifdef __AL_CODES */
