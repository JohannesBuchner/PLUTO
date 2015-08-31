/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file 
  \brief Creates a datatype describing a subarray of a 
         multidimensional array

  Creates a datatype describing a subarray of a 
  multidimensional array
  
  \authors A. Malagoli (University of Chicago)
  
  \date Jul 17, 1999 
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "al_hidden.h"

/* ********************************************************************* */
int AL_Type_create_subarray(int ndims, int *array_of_sizes, 
                             int *array_of_subsizes, int *array_of_starts,
                             int order, MPI_Datatype oldtype, 
                             MPI_Datatype *newtype)
/*!
 * Creates a datatype describing a subarray of a multidimensional array
 * NOTE:
 * This routine has been modified from R. Thakur's ROMIO implementation
 * of MPI_Type_create_subarry. The reason for keeping a local copy is 
 * mainly to reduce this package's dependence on ROMIO. It only affects
 * the non I/O components of the library (e.g. AL_Allgather).
 *
 * \param [in]   ndims               number of array dimensions 
 * \param [in]   array_of_sizes      number of elements of type
 *                                   oldtype in each dimension of the 
 *                                   full array 
 * \param [in]   array_of_subsizes   number of elements of type 
 *                                   oldtype in each dimension of the 
 *                                   subarray 
 * \param [in]   array_of_starts     starting coordinates of the subarray 
 *                                   in each dimension 
 * \param [in]   order               array storage order flag 
 * \param [in]   oldtype             old datatype (handle)
 * \param [out]  newtype             new datatype (handle)
 *********************************************************************** */
{
    MPI_Aint extent, disps[AL_MAX_DIM], size, size_with_aint;
    int i, blklens[AL_MAX_DIM];
    MPI_Datatype tmp1, tmp2, types[AL_MAX_DIM];
    MPI_Offset size_with_offset;

    if (ndims <= 0) {
        printf("MPI_Type_create_subarray: Invalid ndims argument\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (array_of_sizes <= (int *) 0) {
        printf("MPI_Type_create_subarray: array_of_sizes is an invalid address\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (array_of_subsizes <= (int *) 0) {
        printf("MPI_Type_create_subarray: array_of_subsizes is an invalid address\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (array_of_starts <= (int *) 0) {
        printf("MPI_Type_create_subarray: array_of_starts is an invalid address\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (i=0; i<ndims; i++) {
        if (array_of_sizes[i] <= 0) {
            printf("MPI_Type_create_subarray: Invalid value in array_of_sizes\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (array_of_subsizes[i] <= 0) {
            printf("MPI_Type_create_subarray: Invalid value in array_of_subsizes\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (array_of_starts[i] < 0) {
            printf("MPI_Type_create_subarray: Invalid value in array_of_starts\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    /* order argument checked below */

    if (oldtype == MPI_DATATYPE_NULL) {
        printf("MPI_Type_create_subarray: oldtype is an invalid datatype\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Type_extent(oldtype, (MPI_Aint *) &extent);

/* check if MPI_Aint is large enough for size of global array. 
   if not, complain. */

    size_with_aint = extent;
    for (i=0; i<ndims; i++) size_with_aint *= array_of_sizes[i];
    size_with_offset = extent;
    for (i=0; i<ndims; i++) size_with_offset *= array_of_sizes[i];
    if (size_with_aint != size_with_offset) {
        printf("MPI_Type_create_subarray: Can't use an array of this size unless the MPI implementation defines a 64-bit MPI_Aint\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (order == AL_ORDER_FORTRAN) {
      /* dimension 0 changes fastest */
        if (ndims == 1)
            MPI_Type_contiguous(array_of_subsizes[0], oldtype, &tmp1);
        else {
            MPI_Type_vector(array_of_subsizes[1], array_of_subsizes[0],
                            array_of_sizes[0], oldtype, &tmp1);
            
            size = array_of_sizes[0]*extent;
            for (i=2; i<ndims; i++) {
                size *= array_of_sizes[i-1];
                MPI_Type_hvector(array_of_subsizes[i], 1, size, tmp1, &tmp2);
                MPI_Type_free(&tmp1);
                tmp1 = tmp2;
            }
        }
        
        /* add displacement and UB */
        
        disps[1] = array_of_starts[0];
        size = 1;
        for (i=1; i<ndims; i++) {
            size *= array_of_sizes[i-1];
            disps[1] += size*array_of_starts[i];
        }  
        /* rest done below for both Fortran and C order */
    }

    else if (order == AL_ORDER_C) {
        /* dimension ndims-1 changes fastest */
        if (ndims == 1)
            MPI_Type_contiguous(array_of_subsizes[0], oldtype, &tmp1);
        else {
            MPI_Type_vector(array_of_subsizes[ndims-2],
                            array_of_subsizes[ndims-1],
                            array_of_sizes[ndims-1], oldtype, &tmp1);
            
            size = array_of_sizes[ndims-1]*extent;
            for (i=ndims-3; i>=0; i--) {
                size *= array_of_sizes[i+1];
                MPI_Type_hvector(array_of_subsizes[i], 1, size, tmp1, &tmp2);
                MPI_Type_free(&tmp1);
                tmp1 = tmp2;
            }
        }
        
        /* add displacement and UB */
        
        disps[1] = array_of_starts[ndims-1];
        size = 1;
        for (i=ndims-2; i>=0; i--) {
            size *= array_of_sizes[i+1];
            disps[1] += size*array_of_starts[i];
        }
    }
    else {
        printf("MPI_Type_create_subarray: Invalid order argument\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    disps[1] *= extent;
    
    disps[2] = extent;
    for (i=0; i<ndims; i++) disps[2] *= array_of_sizes[i];
    
    disps[0] = 0;
    blklens[0] = blklens[1] = blklens[2] = 1;
    types[0] = MPI_LB;
    types[1] = tmp1;
    types[2] = MPI_UB;
    
    MPI_Type_struct(3, blklens, disps, types, newtype);

    MPI_Type_free(&tmp1);

    return MPI_SUCCESS;
}
