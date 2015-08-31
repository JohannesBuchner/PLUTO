/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ArrayLib function definitions header file.

  Contains function definitions used in ArrayLib routines.

  \author A. Malagoli (University of Chicago)
  \author G. Muscianisi (g.muscianisi@cineca.it)

  \date   Aug 29, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#ifndef __AL_DEFS
#define __AL_DEFS


/* Some definitions for more useful generalizations */
#define AL_Const int

/* The datatype type */
#define AL_Datatype MPI_Datatype

/* Amount of grid overlap in staggered meshes */
#define AL_STAGGERED_OVERLAP 1

/*
  AL_ALLOC_ -- AL memory allocation procedure. 
               We define it temporarily as a macro,
               but it should really become its own
               function at some point.
*/
#define AL_ALLOC_(nelem,size) malloc((nelem)*(size))
#define AL_CALLOC_(nelem,size) calloc((nelem),(size))

/*
  AL_FREE_ -- AL memory freeing procedure. 
              We define it temporarily as a macro,
              but it should really become its own
              function at some point.
*/
#define AL_FREE_(ptr) free((ptr))

/*
  AL_POWEROF2 -- Return 1 if the number is a power of 2
                         0 otherwise
*/
#define AL_POWEROF2(x)  ((((x)-1)&(x))==0)

/* 
  AL_ISEVEN -- Return 1 if the number is even, 0 otherwise
  AL_ISODD  -- Return 1 if the number is odd,  0 otherwise
*/
#define AL_ISEVEN ((x)-(x)/2*2)
#define AL_ISODD ((x)/2*2-(x)+1)

/*
  AL_MAX, AL_MIN -- Max and Min definitions
 */
#define AL_ISMAX(a,b)  ((a) > (b) ? (a) : (b))
#define AL_ISMIN(a,b)  ((a) < (b) ? (a) : (b))

/* End ifndef __AL_DEFS */
#endif 
