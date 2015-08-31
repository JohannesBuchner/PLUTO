/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file 
  \brief  Dellocate a distributed array descriptor

  Dellocate a distributed array descriptor
  
  \author A. Malagoli (University of Chicago)
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

/* ********************************************************************* */
int AL_Sz_free(int sz_ptr)
/*!
 * Dellocate a distributed array descriptor
 *
 * \param [in] sz_ptr Integer pointer to the array descriptor
 *********************************************************************** */
{
  /*
    Deallocate the SZ structure 
  */
  return AL_Deallocate_sz_(sz_ptr);

}

