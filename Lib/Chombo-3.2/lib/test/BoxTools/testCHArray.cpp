#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*------------------------------------------------------------------------------
 * Note:  This file is tested in two locations:
 *
 * test/ChomboFortran - with CHARRAYTESTCHF defined
 * test/BoxTools -      with CHARRAYTESTCPP defined
 *
 * If it is edited, it should be copied to the other location with the correct
 * define statements uncommented.  '#include "testCHArray*F_F.H"' has to be
 * explicitly commented out if 'testCHArray*F_F.H' is not available
 *----------------------------------------------------------------------------*/

// #define CHARRAYTESTCHF
#define CHARRAYTESTCPP

#include <iostream>
#include <cstring>
using std::endl;
#include "parstream.H"
#include "BoxIterator.H"
#include "CHArray.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

#ifdef CHARRAYTESTCHF
// #include "testCHArrayF_F.H"
#if CH_SPACEDIM<6
// The macros fail during expansion at SpaceDim=6 since rank 8 is not supported
// #include "testCHArray_rankDimPlus2F_F.H"
#endif
#endif

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testCHArray();

/// Global variables for handling output:
static const char *pgmname = "testCHArray" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if (verbose)
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testCHArray() ;

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


/*******************************************************************************
 *
 * Routine testCHArray
 *
 ******************************************************************************/

int
testCHArray()
{
  int errors = 0;

/*--------------------------------------------------------------------*
 * Parameters
 *--------------------------------------------------------------------*/

  const int dim0 = 8;
  const int dim1 = 7;
  const int dim2 = 6;
  const int dim3 = 5;
  const int dim4 = 4;
  const int dim5 = 3;
  const int dim6 = 2;

  const int szR1 = dim0;
  const int szR2 = szR1*dim1;
  const int szR3 = szR2*dim2;
  const int szR4 = szR3*dim3;
  const int szR5 = szR4*dim4;
  const int szR6 = szR5*dim5;
  const int szR7 = szR6*dim6;

  const int dim0df = 9;
  const int dim1df = 8;
  const int dim2df = 7;
  const int dim3df = 6;
  const int dim4df = 5;
  const int dim5df = 4;
  const int dim6df = 3;

  const int szR1df = dim0df;
  const int szR2df = szR1df*dim1df;
  const int szR3df = szR2df*dim2df;
  const int szR4df = szR3df*dim3df;
  const int szR5df = szR4df*dim4df;
  const int szR6df = szR5df*dim5df;
  const int szR7df = szR6df*dim6df;


/*==============================================================================
 * Generic testing is performed at SpaceDim = 2.  We test ranks 1, 2 and 7
 * to account for specializations.  All ranks are tested for column-ordered
 * range subscripting
 *============================================================================*/

  if (SpaceDim == 2)
    {

#ifdef CHARRAYTESTCPP

      // Memory for allocating on existing address
      const int szChunk = szR7;
      const int szByteChunk = szChunk*sizeof(int);
      int chunk[szChunk];

#endif

/*--------------------------------------------------------------------*
 * Misc testing
 *--------------------------------------------------------------------*/

#ifdef CHARRAYTESTCPP

//--Confirm empty base class optimizations are still used correctly

      {
        const int expectedSize = 2*sizeof(void*);
        const int actualSize =
          sizeof(CHArray<int, 1, ArZeroRow, ArSp::NewArrayAlloc<int> >);
        if (actualSize > expectedSize)
          {
            pout() << indent << pgmname
                   << ": Empty base class optimizations violated"
                   << endl;
            pout() << indent2
                   << "Expected class size: " << expectedSize << endl;
            pout() << indent2
                   << "Actual class size  : " << actualSize << endl;
            ++errors;
          }
      }

/*--Test "array of matrix"
 *
 *  Flags
 *    1 : size
 *    2 : access/storage
 */

      {
        int status = 0;
        CHArray<CHMatrix, 1, ArZeroCol, ArSp::ArrayOfMatrixAlloc> arrOfM;
        // Matrices have dimensions dim0 x dim1
        arrOfM.getAllocator().define(dim0, dim1);
        // Allocate dim2 matrices
        arrOfM.define(dim2);
        status |= (arrOfM.size()*arrOfM(0).size() != szR3)*1;
        // Test access/storage
        Real c = 0.;
        for (int iM = 0; iM != dim2; ++iM)
          {
            CHMatrix &m = arrOfM(iM);
            for (int icol = 0; icol != dim1; ++icol)
              for (int irow = 0; irow != dim0; ++irow)
                {
                  m(irow, icol) = c;
                  c += 0.5;
                }
          }
        c = 0.;
        for (Real *p = arrOfM(0).begin(); p != arrOfM(dim2-1).end(); ++p)
          {
            status |= (*p != c)*2;
            c += 0.5;
          }
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 1 array of matrix failed storage testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

#endif

/*--------------------------------------------------------------------*
 * Row ordered, zero-based subscript testing
 *
 * Flags
 *   1 : size
 *   2 : access/storage
 *   4 : storage/const access
 *   8 : assignment/storage
 *   16: bounds
 *--------------------------------------------------------------------*/

#ifdef CHARRAYTESTCPP

//--Rank 1

      {
        int status = 0;
        // Test constructor
        CHArray<int, 1, ArZeroRow> arr(dim0);
        status |= (arr.size() != szR1)*1;
        // Test access/storage
        int c = 0;
        for (int i0 = 0; i0 != dim0; ++i0)
          {
            arr(i0) = ++c;
          }
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(dim0df);
        status |= (arr.size() != szR1df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i0 = 0; i0 != dim0df; ++i0)
          {
            status |= (arr(i0) != ++c)*4;
          }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk, dim0);
        status |= (arr.size() != szR1)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i0 = 0; i0 != dim0; ++i0)
          {
            arr(i0) += c++;
          }
        c = 0;
        for (int ic = 0; ic != szR1; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      0)*16;
        status |= (arr.upperBound(0) != dim0-1)*16;
        status |= (arr.size(0) != dim0)*16;
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 1 array failed zero-based row storage testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 2

      {
        int status = 0;
        // Test constructor
        CHArray<int, 2, ArZeroRow> arr(dim1, dim0);
        status |= (arr.size() != szR2)*1;
        // Test access/storage
        int c = 0;
        for (int i1 = 0; i1 != dim1; ++i1)
          for (int i0 = 0; i0 != dim0; ++i0)
            {
              arr(i1, i0) = ++c;
            }
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(dim1df, dim0df);
        status |= (arr.size() != szR2df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i1 = 0; i1 != dim1df; ++i1)
          for (int i0 = 0; i0 != dim0df; ++i0)
            {
              status |= (arr(i1, i0) != ++c)*4;
            }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk, dim1, dim0);
        status |= (arr.size() != szR2)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i1 = 0; i1 != dim1; ++i1)
          for (int i0 = 0; i0 != dim0; ++i0)
            {
              arr(i1, i0) += c++;
            }
        c = 0;
        for (int ic = 0; ic != szR2; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      0)*16;
        status |= (arr.lowerBound(1) !=      0)*16;
        status |= (arr.upperBound(0) != dim0-1)*16;
        status |= (arr.upperBound(1) != dim1-1)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 2 array failed zero-based row storage testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 7

      {
        int status = 0;
        // Test constructor
        CHArray<int, 7, ArZeroRow> arr(
          dim6, dim5, dim4, dim3, dim2, dim1, dim0);
        status |= (arr.size() != szR7)*1;
        // Test access/storage
        int c = 0;
        for (int i6 = 0; i6 != dim6; ++i6)
          for (int i5 = 0; i5 != dim5; ++i5)
            for (int i4 = 0; i4 != dim4; ++i4)
              for (int i3 = 0; i3 != dim3; ++i3)
                for (int i2 = 0; i2 != dim2; ++i2)
                  for (int i1 = 0; i1 != dim1; ++i1)
                    for (int i0 = 0; i0 != dim0; ++i0)
                      {
                        arr(i6, i5, i4, i3, i2, i1, i0) = ++c;
                      }
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(dim6df, dim5df, dim4df, dim3df, dim2df, dim1df, dim0df);
        status |= (arr.size() != szR7df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i6 = 0; i6 != dim6df; ++i6)
          for (int i5 = 0; i5 != dim5df; ++i5)
            for (int i4 = 0; i4 != dim4df; ++i4)
              for (int i3 = 0; i3 != dim3df; ++i3)
                for (int i2 = 0; i2 != dim2df; ++i2)
                  for (int i1 = 0; i1 != dim1df; ++i1)
                    for (int i0 = 0; i0 != dim0df; ++i0)
                      {
                        status |=
                          (arr(i6, i5, i4, i3, i2, i1, i0) != ++c)*4;
                      }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk, dim6, dim5, dim4, dim3, dim2, dim1, dim0);
        status |= (arr.size() != szR7)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i6 = 0; i6 != dim6; ++i6)
          for (int i5 = 0; i5 != dim5; ++i5)
            for (int i4 = 0; i4 != dim4; ++i4)
              for (int i3 = 0; i3 != dim3; ++i3)
                for (int i2 = 0; i2 != dim2; ++i2)
                  for (int i1 = 0; i1 != dim1; ++i1)
                    for (int i0 = 0; i0 != dim0; ++i0)
                      {
                        arr(i6, i5, i4, i3, i2, i1, i0) += c++;
                      }
        c = 0;
        for (int ic = 0; ic != szR7; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      0)*16;
        status |= (arr.lowerBound(1) !=      0)*16;
        status |= (arr.lowerBound(2) !=      0)*16;
        status |= (arr.lowerBound(3) !=      0)*16;
        status |= (arr.lowerBound(4) !=      0)*16;
        status |= (arr.lowerBound(5) !=      0)*16;
        status |= (arr.lowerBound(6) !=      0)*16;
        status |= (arr.upperBound(0) != dim0-1)*16;
        status |= (arr.upperBound(1) != dim1-1)*16;
        status |= (arr.upperBound(2) != dim2-1)*16;
        status |= (arr.upperBound(3) != dim3-1)*16;
        status |= (arr.upperBound(4) != dim4-1)*16;
        status |= (arr.upperBound(5) != dim5-1)*16;
        status |= (arr.upperBound(6) != dim6-1)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
        status |= (arr.size(2) != dim2)*16;
        status |= (arr.size(3) != dim3)*16;
        status |= (arr.size(4) != dim4)*16;
        status |= (arr.size(5) != dim5)*16;
        status |= (arr.size(6) != dim6)*16;
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 7 array failed zero-based row storage testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

#endif

/*--------------------------------------------------------------------*
 * Column ordered, zero-based subscript testing
 *
 * Flags
 *   1 : size
 *   2 : access/storage
 *   4 : storage/const access
 *   8 : assignment/storage
 *   16: bounds
 *   32: fortran
 *--------------------------------------------------------------------*/

//--Rank 1

      {
        int status = 0;
        // Test constructor
        CHArray<int, 1, ArZeroCol> arr(dim0);
        status |= (arr.size() != szR1)*1;
        // Test access/storage
        int c = 0;
        for (int i0 = 0; i0 != dim0; ++i0)
          {
            arr(i0) = ++c;
          }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(dim0df);
        status |= (arr.size() != szR1df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i0 = 0; i0 != dim0df; ++i0)
          {
            status |= (arr(i0) != ++c)*4;
          }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk, dim0);
        status |= (arr.size() != szR1)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i0 = 0; i0 != dim0; ++i0)
          {
            arr(i0) += c++;
          }
        c = 0;
        for (int ic = 0; ic != szR1; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      0)*16;
        status |= (arr.upperBound(0) != dim0-1)*16;
        status |= (arr.size(0) != dim0)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 1, ArZeroCol> arrf(dim0);
        FORT_RANK1CHARRAY(
          CHF_ICHARRAY(1, arrf),
          CHF_CONST_ICHARRAY(1, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 1 array failed zero-based column storage testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 2

      {
        int status = 0;
        // Test constructor
        CHArray<int, 2, ArZeroCol> arr(dim0, dim1);
        status |= (arr.size() != szR2)*1;
        // Test access/storage
        int c = 0;
        for (int i1 = 0; i1 != dim1; ++i1)
          for (int i0 = 0; i0 != dim0; ++i0)
            {
              arr(i0, i1) = ++c;
            }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(dim0df, dim1df);
        status |= (arr.size() != szR2df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i1 = 0; i1 != dim1df; ++i1)
          for (int i0 = 0; i0 != dim0df; ++i0)
            {
              status |= (arr(i0, i1) != ++c)*4;
            }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk, dim0, dim1);
        status |= (arr.size() != szR2)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i1 = 0; i1 != dim1; ++i1)
          for (int i0 = 0; i0 != dim0; ++i0)
            {
              arr(i0, i1) += c++;
            }
        c = 0;
        for (int ic = 0; ic != szR2; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      0)*16;
        status |= (arr.lowerBound(1) !=      0)*16;
        status |= (arr.upperBound(0) != dim0-1)*16;
        status |= (arr.upperBound(1) != dim1-1)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 2, ArZeroCol> arrf(dim0, dim1);
        FORT_RANK2CHARRAY(
          CHF_ICHARRAY(2, arrf),
          CHF_CONST_ICHARRAY(2, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 2 array failed zero-based column storage testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 7

      {
        int status = 0;
        // Test constructor
        CHArray<int, 7, ArZeroCol> arr(
          dim0, dim1, dim2, dim3, dim4, dim5, dim6);
        status |= (arr.size() != szR7)*1;
        // Test access/storage
        int c = 0;
        for (int i6 = 0; i6 != dim6; ++i6)
          for (int i5 = 0; i5 != dim5; ++i5)
            for (int i4 = 0; i4 != dim4; ++i4)
              for (int i3 = 0; i3 != dim3; ++i3)
                for (int i2 = 0; i2 != dim2; ++i2)
                  for (int i1 = 0; i1 != dim1; ++i1)
                    for (int i0 = 0; i0 != dim0; ++i0)
                      {
                        arr(i0, i1, i2, i3, i4, i5, i6) = ++c;
                      }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(dim0df, dim1df, dim2df, dim3df, dim4df, dim5df, dim6df);
        status |= (arr.size() != szR7df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i6 = 0; i6 != dim6df; ++i6)
          for (int i5 = 0; i5 != dim5df; ++i5)
            for (int i4 = 0; i4 != dim4df; ++i4)
              for (int i3 = 0; i3 != dim3df; ++i3)
                for (int i2 = 0; i2 != dim2df; ++i2)
                  for (int i1 = 0; i1 != dim1df; ++i1)
                    for (int i0 = 0; i0 != dim0df; ++i0)
                      {
                        status |= (arr(i0, i1, i2, i3, i4, i5, i6) != ++c)*4;
                      }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk, dim0, dim1, dim2, dim3, dim4, dim5, dim6);
        status |= (arr.size() != szR7)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i6 = 0; i6 != dim6; ++i6)
          for (int i5 = 0; i5 != dim5; ++i5)
            for (int i4 = 0; i4 != dim4; ++i4)
              for (int i3 = 0; i3 != dim3; ++i3)
                for (int i2 = 0; i2 != dim2; ++i2)
                  for (int i1 = 0; i1 != dim1; ++i1)
                    for (int i0 = 0; i0 != dim0; ++i0)
                      {
                        arr(i0, i1, i2, i3, i4, i5, i6) += c++;
                      }
        c = 0;
        for (int ic = 0; ic != szR7; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      0)*16;
        status |= (arr.lowerBound(1) !=      0)*16;
        status |= (arr.lowerBound(2) !=      0)*16;
        status |= (arr.lowerBound(3) !=      0)*16;
        status |= (arr.lowerBound(4) !=      0)*16;
        status |= (arr.lowerBound(5) !=      0)*16;
        status |= (arr.lowerBound(6) !=      0)*16;
        status |= (arr.upperBound(0) != dim0-1)*16;
        status |= (arr.upperBound(1) != dim1-1)*16;
        status |= (arr.upperBound(2) != dim2-1)*16;
        status |= (arr.upperBound(3) != dim3-1)*16;
        status |= (arr.upperBound(4) != dim4-1)*16;
        status |= (arr.upperBound(5) != dim5-1)*16;
        status |= (arr.upperBound(6) != dim6-1)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
        status |= (arr.size(2) != dim2)*16;
        status |= (arr.size(3) != dim3)*16;
        status |= (arr.size(4) != dim4)*16;
        status |= (arr.size(5) != dim5)*16;
        status |= (arr.size(6) != dim6)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 7, ArZeroCol> arrf(
          dim0, dim1, dim2, dim3, dim4, dim5, dim6);
        FORT_RANK7CHARRAY(
          CHF_ICHARRAY(7, arrf),
          CHF_CONST_ICHARRAY(7, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 7 array failed zero-based column storage testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

/*--------------------------------------------------------------------*
 * Row ordered, range subscript testing.
 *
 * Testing here is not needed since there is nothing new outside of
 * ArZeroRow testing and ArRangeCol testing (because of the way
 * templates are used).
 *--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*
 * Column ordered, range subscript testing.
 *
 * The ranges are set as
 *   dim 0 = CHRange(-1, dim0-2)
 *   dim 1 = sz (this is equivalent to CHRange(0, dim1-1))
 *   dim 2 = CHRange(-1, dim2-2)
 *   ...
 *
 * Flags
 *   1 : size
 *   2 : access/storage
 *   4 : storage/const access
 *   8 : assignment/storage
 *   16: bounds
 *   32: fortran
 *--------------------------------------------------------------------*/

//--Rank 1

      {
        int status = 0;
        // Test constructor
        CHArray<int, 1, ArRangeCol> arr(
          CHRange(-1,dim0-2));
        status |= (arr.size() != szR1)*1;
        // Test access/storage
        int c = 0;
        for (int i0 = -1; i0 <= dim0-2; ++i0)
          {
            arr(i0) = ++c;
          }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(
          CHRange(-1,dim0df-2));
        status |= (arr.size() != szR1df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i0 = -1; i0 <= dim0df-2; ++i0)
          {
            status |= (arr(i0) != ++c)*4;
          }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk,
                   CHRange(-1,dim0-2));
        status |= (arr.size() != szR1)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i0 = -1; i0 <= dim0-2; ++i0)
          {
            arr(i0) += c++;
          }
        c = 0;
        for (int ic = 0; ic != szR1; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=     -1)*16;
        status |= (arr.upperBound(0) != dim0-2)*16;
        status |= (arr.size(0) != dim0)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 1, ArRangeCol> arrf(
          CHRange(-1,dim0-2));
        FORT_RANK1CHARRAY(
          CHF_ICHARRAY(1, arrf),
          CHF_CONST_ICHARRAY(1, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 1 array failed ranged subscript column storage "
              "testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 2

      {
        int status = 0;
        // Test constructor
        CHArray<int, 2, ArRangeCol> arr(
          CHRange(-1,dim0-2),
          dim1);
        status |= (arr.size() != szR2)*1;
        // Test access/storage
        int c = 0;
        for (int i1 = 0; i1 != dim1; ++i1)
          for (int i0 = -1; i0 <= dim0-2; ++i0)
            {
              arr(i0, i1) = ++c;
            }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(
          CHRange(-1,dim0df-2),
          dim1df);
        status |= (arr.size() != szR2df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i1 = 0; i1 != dim1df; ++i1)
          for (int i0 = -1; i0 <= dim0df-2; ++i0)
            {
              status |= (arr(i0, i1) != ++c)*4;
            }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk,
                   CHRange(-1,dim0-2),
                   dim1);
        status |= (arr.size() != szR2)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i1 = 0; i1 != dim1; ++i1)
          for (int i0 = -1; i0 <= dim0-2; ++i0)
            {
              arr(i0, i1) += c++;
            }
        c = 0;
        for (int ic = 0; ic != szR2; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      -1)*16;
        status |= (arr.lowerBound(1) !=       0)*16;
        status |= (arr.upperBound(0) !=  dim0-2)*16;
        status |= (arr.upperBound(1) !=  dim1-1)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 2, ArRangeCol> arrf(
          CHRange(-1,dim0-2),
          dim1);
        FORT_RANK2CHARRAY(
          CHF_ICHARRAY(2, arrf),
          CHF_CONST_ICHARRAY(2, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 2 array failed ranged subscript column storage "
              "testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 3

      {
        int status = 0;
        // Test constructor
        CHArray<int, 3, ArRangeCol> arr(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2));
        status |= (arr.size() != szR3)*1;
        // Test access/storage
        int c = 0;
        for (int i2 = -1; i2 <= dim2-2; ++i2)
          for (int i1 = 0; i1 != dim1; ++i1)
            for (int i0 = -1; i0 <= dim0-2; ++i0)
              {
                arr(i0, i1, i2) = ++c;
              }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(
          CHRange(-1,dim0df-2),
          dim1df,
          CHRange(-1,dim2df-2));
        status |= (arr.size() != szR3df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i2 = -1; i2 <= dim2df-2; ++i2)
          for (int i1 = 0; i1 != dim1df; ++i1)
            for (int i0 = -1; i0 <= dim0df-2; ++i0)
              {
                status |= (arr(i0, i1, i2) != ++c)*4;
              }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk,
                   CHRange(-1,dim0-2),
                   dim1,
                   CHRange(-1,dim2-2));
        status |= (arr.size() != szR3)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i2 = -1; i2 <= dim2-2; ++i2)
          for (int i1 = 0; i1 != dim1; ++i1)
            for (int i0 = -1; i0 <= dim0-2; ++i0)
              {
                arr(i0, i1, i2) += c++;
              }
        c = 0;
        for (int ic = 0; ic != szR3; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      -1)*16;
        status |= (arr.lowerBound(1) !=       0)*16;
        status |= (arr.lowerBound(2) !=      -1)*16;
        status |= (arr.upperBound(0) !=  dim0-2)*16;
        status |= (arr.upperBound(1) !=  dim1-1)*16;
        status |= (arr.upperBound(2) !=  dim2-2)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
        status |= (arr.size(2) != dim2)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 3, ArRangeCol> arrf(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2));
        FORT_RANK3CHARRAY(
          CHF_ICHARRAY(3, arrf),
          CHF_CONST_ICHARRAY(3, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 3 array failed ranged subscript column storage "
              "testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 4

      {
        int status = 0;
        // Test constructor
        CHArray<int, 4, ArRangeCol> arr(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2),
          CHRange(-1,dim3-2));
        status |= (arr.size() != szR4)*1;
        // Test access/storage
        int c = 0;
        for (int i3 = -1; i3 <= dim3-2; ++i3)
          for (int i2 = -1; i2 <= dim2-2; ++i2)
            for (int i1 = 0; i1 != dim1; ++i1)
              for (int i0 = -1; i0 <= dim0-2; ++i0)
                {
                  arr(i0, i1, i2, i3) = ++c;
                }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(
          CHRange(-1,dim0df-2),
          dim1df,
          CHRange(-1,dim2df-2),
          CHRange(-1,dim3df-2));
        status |= (arr.size() != szR4df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i3 = -1; i3 <= dim3df-2; ++i3)
          for (int i2 = -1; i2 <= dim2df-2; ++i2)
            for (int i1 = 0; i1 != dim1df; ++i1)
              for (int i0 = -1; i0 <= dim0df-2; ++i0)
                {
                  status |= (arr(i0, i1, i2, i3) != ++c)*4;
                }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk,
                   CHRange(-1,dim0-2),
                   dim1,
                   CHRange(-1,dim2-2),
                   CHRange(-1,dim3-2));
        status |= (arr.size() != szR4)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i3 = -1; i3 <= dim3-2; ++i3)
          for (int i2 = -1; i2 <= dim2-2; ++i2)
            for (int i1 = 0; i1 != dim1; ++i1)
              for (int i0 = -1; i0 <= dim0-2; ++i0)
                {
                  arr(i0, i1, i2, i3) += c++;
                }
        c = 0;
        for (int ic = 0; ic != szR4; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      -1)*16;
        status |= (arr.lowerBound(1) !=       0)*16;
        status |= (arr.lowerBound(2) !=      -1)*16;
        status |= (arr.lowerBound(3) !=      -1)*16;
        status |= (arr.upperBound(0) !=  dim0-2)*16;
        status |= (arr.upperBound(1) !=  dim1-1)*16;
        status |= (arr.upperBound(2) !=  dim2-2)*16;
        status |= (arr.upperBound(3) !=  dim3-2)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
        status |= (arr.size(2) != dim2)*16;
        status |= (arr.size(3) != dim3)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 4, ArRangeCol> arrf(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2),
          CHRange(-1,dim3-2));
        FORT_RANK4CHARRAY(
          CHF_ICHARRAY(4, arrf),
          CHF_CONST_ICHARRAY(4, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 4 array failed ranged subscript column storage "
              "testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 5

      {
        int status = 0;
        // Test constructor
        CHArray<int, 5, ArRangeCol> arr(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2),
          CHRange(-1,dim3-2),
          CHRange(-1,dim4-2));
        status |= (arr.size() != szR5)*1;
        // Test access/storage
        int c = 0;
        for (int i4 = -1; i4 <= dim4-2; ++i4)
          for (int i3 = -1; i3 <= dim3-2; ++i3)
            for (int i2 = -1; i2 <= dim2-2; ++i2)
              for (int i1 = 0; i1 != dim1; ++i1)
                for (int i0 = -1; i0 <= dim0-2; ++i0)
                  {
                    arr(i0, i1, i2, i3, i4) = ++c;
                  }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(
          CHRange(-1,dim0df-2),
          dim1df,
          CHRange(-1,dim2df-2),
          CHRange(-1,dim3df-2),
          CHRange(-1,dim4df-2));
        status |= (arr.size() != szR5df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i4 = -1; i4 <= dim4df-2; ++i4)
          for (int i3 = -1; i3 <= dim3df-2; ++i3)
            for (int i2 = -1; i2 <= dim2df-2; ++i2)
              for (int i1 = 0; i1 != dim1df; ++i1)
                for (int i0 = -1; i0 <= dim0df-2; ++i0)
                  {
                    status |= (arr(i0, i1, i2, i3, i4) != ++c)*4;
                  }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk,
                   CHRange(-1,dim0-2),
                   dim1,
                   CHRange(-1,dim2-2),
                   CHRange(-1,dim3-2),
                   CHRange(-1,dim4-2));
        status |= (arr.size() != szR5)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i4 = -1; i4 <= dim4-2; ++i4)
          for (int i3 = -1; i3 <= dim3-2; ++i3)
            for (int i2 = -1; i2 <= dim2-2; ++i2)
              for (int i1 = 0; i1 != dim1; ++i1)
                for (int i0 = -1; i0 <= dim0-2; ++i0)
                  {
                    arr(i0, i1, i2, i3, i4) += c++;
                  }
        c = 0;
        for (int ic = 0; ic != szR5; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      -1)*16;
        status |= (arr.lowerBound(1) !=       0)*16;
        status |= (arr.lowerBound(2) !=      -1)*16;
        status |= (arr.lowerBound(3) !=      -1)*16;
        status |= (arr.lowerBound(4) !=      -1)*16;
        status |= (arr.upperBound(0) !=  dim0-2)*16;
        status |= (arr.upperBound(1) !=  dim1-1)*16;
        status |= (arr.upperBound(2) !=  dim2-2)*16;
        status |= (arr.upperBound(3) !=  dim3-2)*16;
        status |= (arr.upperBound(4) !=  dim4-2)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
        status |= (arr.size(2) != dim2)*16;
        status |= (arr.size(3) != dim3)*16;
        status |= (arr.size(4) != dim4)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 5, ArRangeCol> arrf(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2),
          CHRange(-1,dim3-2),
          CHRange(-1,dim4-2));
        FORT_RANK5CHARRAY(
          CHF_ICHARRAY(5, arrf),
          CHF_CONST_ICHARRAY(5, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 5 array failed ranged subscript column storage "
              "testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 6

      {
        int status = 0;
        // Test constructor
        CHArray<int, 6, ArRangeCol> arr(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2),
          CHRange(-1,dim3-2),
          CHRange(-1,dim4-2),
          CHRange(-1,dim5-2));
        status |= (arr.size() != szR6)*1;
        // Test access/storage
        int c = 0;
        for (int i5 = -1; i5 <= dim5-2; ++i5)
          for (int i4 = -1; i4 <= dim4-2; ++i4)
            for (int i3 = -1; i3 <= dim3-2; ++i3)
              for (int i2 = -1; i2 <= dim2-2; ++i2)
                for (int i1 = 0; i1 != dim1; ++i1)
                  for (int i0 = -1; i0 <= dim0-2; ++i0)
                    {
                      arr(i0, i1, i2, i3, i4, i5) = ++c;
                    }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(
          CHRange(-1,dim0df-2),
          dim1df,
          CHRange(-1,dim2df-2),
          CHRange(-1,dim3df-2),
          CHRange(-1,dim4df-2),
          CHRange(-1,dim5df-2));
        status |= (arr.size() != szR6df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i5 = -1; i5 <= dim5df-2; ++i5)
          for (int i4 = -1; i4 <= dim4df-2; ++i4)
            for (int i3 = -1; i3 <= dim3df-2; ++i3)
              for (int i2 = -1; i2 <= dim2df-2; ++i2)
                for (int i1 = 0; i1 != dim1df; ++i1)
                  for (int i0 = -1; i0 <= dim0df-2; ++i0)
                    {
                      status |= (arr(i0, i1, i2, i3, i4, i5) != ++c)*4;
                    }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk,
                   CHRange(-1,dim0-2),
                   dim1,
                   CHRange(-1,dim2-2),
                   CHRange(-1,dim3-2),
                   CHRange(-1,dim4-2),
                   CHRange(-1,dim5-2));
        status |= (arr.size() != szR6)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i5 = -1; i5 <= dim5-2; ++i5)
          for (int i4 = -1; i4 <= dim4-2; ++i4)
            for (int i3 = -1; i3 <= dim3-2; ++i3)
              for (int i2 = -1; i2 <= dim2-2; ++i2)
                for (int i1 = 0; i1 != dim1; ++i1)
                  for (int i0 = -1; i0 <= dim0-2; ++i0)
                    {
                      arr(i0, i1, i2, i3, i4, i5) += c++;
                    }
        c = 0;
        for (int ic = 0; ic != szR6; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      -1)*16;
        status |= (arr.lowerBound(1) !=       0)*16;
        status |= (arr.lowerBound(2) !=      -1)*16;
        status |= (arr.lowerBound(3) !=      -1)*16;
        status |= (arr.lowerBound(4) !=      -1)*16;
        status |= (arr.lowerBound(5) !=      -1)*16;
        status |= (arr.upperBound(0) !=  dim0-2)*16;
        status |= (arr.upperBound(1) !=  dim1-1)*16;
        status |= (arr.upperBound(2) !=  dim2-2)*16;
        status |= (arr.upperBound(3) !=  dim3-2)*16;
        status |= (arr.upperBound(4) !=  dim4-2)*16;
        status |= (arr.upperBound(5) !=  dim5-2)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
        status |= (arr.size(2) != dim2)*16;
        status |= (arr.size(3) != dim3)*16;
        status |= (arr.size(4) != dim4)*16;
        status |= (arr.size(5) != dim5)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 6, ArRangeCol> arrf(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2),
          CHRange(-1,dim3-2),
          CHRange(-1,dim4-2),
          CHRange(-1,dim5-2));
        FORT_RANK6CHARRAY(
          CHF_ICHARRAY(6, arrf),
          CHF_CONST_ICHARRAY(6, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 6 array failed ranged subscript column storage "
              "testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

//--Rank 7

      {
        int status = 0;
        // Test constructor
        CHArray<int, 7, ArRangeCol> arr(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2),
          CHRange(-1,dim3-2),
          CHRange(-1,dim4-2),
          CHRange(-1,dim5-2),
          CHRange(-1,dim6-2));
        status |= (arr.size() != szR7)*1;
        // Test access/storage
        int c = 0;
        for (int i6 = -1; i6 <= dim6-2; ++i6)
          for (int i5 = -1; i5 <= dim5-2; ++i5)
            for (int i4 = -1; i4 <= dim4-2; ++i4)
              for (int i3 = -1; i3 <= dim3-2; ++i3)
                for (int i2 = -1; i2 <= dim2-2; ++i2)
                  for (int i1 = 0; i1 != dim1; ++i1)
                    for (int i0 = -1; i0 <= dim0-2; ++i0)
                      {
                        arr(i0, i1, i2, i3, i4, i5, i6) = ++c;
                      }
#ifdef CHARRAYTESTCPP
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            status |= (*p != ++c)*2;
          }
        // Test define
        arr.define(
          CHRange(-1,dim0df-2),
          dim1df,
          CHRange(-1,dim2df-2),
          CHRange(-1,dim3df-2),
          CHRange(-1,dim4df-2),
          CHRange(-1,dim5df-2),
          CHRange(-1,dim6df-2));
        status |= (arr.size() != szR7df)*1;
        // Test storage/const access
        c = 0;
        for (int* p = arr.begin(); p != arr.end(); ++p)
          {
            *p = ++c;
          }
        c = 0;
        for (int i6 = -1; i6 <= dim6df-2; ++i6)
          for (int i5 = -1; i5 <= dim5df-2; ++i5)
            for (int i4 = -1; i4 <= dim4df-2; ++i4)
              for (int i3 = -1; i3 <= dim3df-2; ++i3)
                for (int i2 = -1; i2 <= dim2df-2; ++i2)
                  for (int i1 = 0; i1 != dim1df; ++i1)
                    for (int i0 = -1; i0 <= dim0df-2; ++i0)
                      {
                        status |= (arr(i0, i1, i2, i3, i4, i5, i6) != ++c)*4;
                      }
        // Test define on address
        std::memset(chunk, 0, szByteChunk);
        arr.define(chunk,
                   CHRange(-1,dim0-2),
                   dim1,
                   CHRange(-1,dim2-2),
                   CHRange(-1,dim3-2),
                   CHRange(-1,dim4-2),
                   CHRange(-1,dim5-2),
                   CHRange(-1,dim6-2));
        status |= (arr.size() != szR7)*1;
        // Test assignment/storage
        arr = 1;
        c = 0;
        for (int i6 = -1; i6 <= dim6-2; ++i6)
          for (int i5 = -1; i5 <= dim5-2; ++i5)
            for (int i4 = -1; i4 <= dim4-2; ++i4)
              for (int i3 = -1; i3 <= dim3-2; ++i3)
                for (int i2 = -1; i2 <= dim2-2; ++i2)
                  for (int i1 = 0; i1 != dim1; ++i1)
                    for (int i0 = -1; i0 <= dim0-2; ++i0)
                      {
                        arr(i0, i1, i2, i3, i4, i5, i6) += c++;
                      }
        c = 0;
        for (int ic = 0; ic != szR7; ++ic)
          {
            status |= (chunk[ic] != ++c)*8;
          }
        // Test bounds
        status |= (arr.lowerBound(0) !=      -1)*16;
        status |= (arr.lowerBound(1) !=       0)*16;
        status |= (arr.lowerBound(2) !=      -1)*16;
        status |= (arr.lowerBound(3) !=      -1)*16;
        status |= (arr.lowerBound(4) !=      -1)*16;
        status |= (arr.lowerBound(5) !=      -1)*16;
        status |= (arr.lowerBound(6) !=      -1)*16;
        status |= (arr.upperBound(0) !=  dim0-2)*16;
        status |= (arr.upperBound(1) !=  dim1-1)*16;
        status |= (arr.upperBound(2) !=  dim2-2)*16;
        status |= (arr.upperBound(3) !=  dim3-2)*16;
        status |= (arr.upperBound(4) !=  dim4-2)*16;
        status |= (arr.upperBound(5) !=  dim5-2)*16;
        status |= (arr.upperBound(6) !=  dim6-2)*16;
        status |= (arr.size(0) != dim0)*16;
        status |= (arr.size(1) != dim1)*16;
        status |= (arr.size(2) != dim2)*16;
        status |= (arr.size(3) != dim3)*16;
        status |= (arr.size(4) != dim4)*16;
        status |= (arr.size(5) != dim5)*16;
        status |= (arr.size(6) != dim6)*16;
#endif
#ifdef CHARRAYTESTCHF
        // Test assignment in Fortran
        CHArray<int, 7, ArRangeCol> arrf(
          CHRange(-1,dim0-2),
          dim1,
          CHRange(-1,dim2-2),
          CHRange(-1,dim3-2),
          CHRange(-1,dim4-2),
          CHRange(-1,dim5-2),
          CHRange(-1,dim6-2));
        FORT_RANK7CHARRAY(
          CHF_ICHARRAY(7, arrf),
          CHF_CONST_ICHARRAY(7, arr));  // arrf = arr + inc
        c = 0;
        int cinc = 0;
        for (int* p = arrf.begin(); p != arrf.end(); ++p)
          {
            status |= (*p != ((++c) + (++cinc)))*32;
          }
#endif
        if (status)
          {
            pout() << indent << pgmname
                   << ": Rank 7 array failed ranged subscript column storage "
              "testing"
                   << endl;
            pout() << indent2
                   << "flags: " << status << endl;
            ++errors;
          }
      }

    }  // End if SpaceDim == 2


/*==============================================================================
 * SpaceDim specific testing
 *============================================================================*/

//--The IntVect and Box used for testing

  const IntVect dimiv(D_DECL6(dim0, dim1, dim2, dim3, dim4, dim5));
  const IntVect ivlo(D_DECL6(-1, -2, -3, 1, 2, 3));
  const IntVect ivhi = ivlo + dimiv - IntVect::Unit;
  const Box dimbox(ivlo, ivhi);
  const IntVect dimivdf(
    D_DECL6(dim0df, dim1df, dim2df, dim3df, dim4df, dim5df));
  const IntVect ivhidf = ivlo + dimivdf - IntVect::Unit;
  const Box dimboxdf(ivlo, ivhidf);

  int dimc0 = 0;
  int dimc1 = 0;
  int szRc0 = 0;
  int szRc1 = 0;
  int dimc0df = 0;
  int dimc1df = 0;
  int szRc0df = 0;
  int szRc1df = 0;
  switch (SpaceDim)
    {
    case 1:
      dimc0 = dim1;
      dimc1 = dim2;
      szRc0 = szR2;
      szRc1 = szR3;
      dimc0df = dim1df;
      dimc1df = dim2df;
      szRc0df = szR2df;
      szRc1df = szR3df;
      break;
    case 2:
      dimc0 = dim2;
      dimc1 = dim3;
      szRc0 = szR3;
      szRc1 = szR4;
      dimc0df = dim2df;
      dimc1df = dim3df;
      szRc0df = szR3df;
      szRc1df = szR4df;
      break;
    case 3:
      dimc0 = dim3;
      dimc1 = dim4;
      szRc0 = szR4;
      szRc1 = szR5;
      dimc0df = dim3df;
      dimc1df = dim4df;
      szRc0df = szR4df;
      szRc1df = szR5df;
      break;
    case 4:
      dimc0 = dim4;
      dimc1 = dim5;
      szRc0 = szR5;
      szRc1 = szR6;
      dimc0df = dim4df;
      dimc1df = dim5df;
      szRc0df = szR5df;
      szRc1df = szR6df;
      break;
    case 5:
      dimc0 = dim5;
      dimc1 = dim6;
      szRc0 = szR6;
      szRc1 = szR7;
      dimc0df = dim5df;
      dimc1df = dim6df;
      szRc0df = szR6df;
      szRc1df = szR7df;
      break;
    case 6:  // Unable to test two components at SpaceDim = 6
      dimc0 = dim6;
      szRc0 = szR7;
      dimc0df = dim6df;
      szRc0df = szR7df;
      break;
    }

#ifdef CHARRAYTESTCPP

  // Memory for allocating on existing address
  const int szChunk = szR7;
  const int szByteChunk = szChunk*sizeof(Real);
  Real chunk[szChunk];

#endif

/*--------------------------------------------------------------------*
 * Row ordered, zero-based subscript testing with IntVects.
 *
 * Only test constructor/access/storage to test IntVect macros.
 * Remaining features are tested using column-ordered storage
 *
 * Flags
 *   1 : size
 *   2 : access/storage
 *   4 : bounds
 *--------------------------------------------------------------------*/

#ifdef CHARRAYTESTCPP

//--+1 dimension, (continguous components)

  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+1, ArZeroRow> arr(dimiv, dimc0);
    status |= (arr.size() != szRc0)*1;
    // Test access/storage
    Real c = 0.;
    for (BoxIterator bit(Box(IntVect::Zero, dimiv - IntVect::Unit)); bit.ok();
         ++bit)
      {
        const IntVect iv = bit();
        for (int c0 = 0; c0 != dimc0; ++c0)
          {
            arr(iv, c0) = c;
            c += 0.5;
          }
      }
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test bounds
    status |= (arr.lowerBound(0) !=       0)*4;
    status |= (arr.upperBound(0) != dimc0-1)*4;
    status |= (arr.size(0) != dimc0)*4;
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(1+idim) !=             0)*4;
        status |= (arr.upperBound(1+idim) != dimiv[idim]-1)*4;
        status |= (arr.size(1+idim) != dimiv[idim])*4;
      }
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+1 array failed zero-based, row storage, "
          "continuous component IntVect testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }

//--+2 dimensions, (continguous components)

#if CH_SPACEDIM<6
  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+2, ArZeroRow> arr(dimiv, dimc1, dimc0);
    status |= (arr.size() != szRc1)*1;
    // Test access/storage
    Real c = 0.;
    for (BoxIterator bit(Box(IntVect::Zero, dimiv - IntVect::Unit)); bit.ok();
         ++bit)
      {
        const IntVect iv = bit();
        for (int c1 = 0; c1 != dimc1; ++c1)
          for (int c0 = 0; c0 != dimc0; ++c0)
            {
              arr(iv, c1, c0) = c;
              c += 0.5;
            }
      }
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test bounds
    status |= (arr.lowerBound(0) !=       0)*4;
    status |= (arr.lowerBound(1) !=       0)*4;
    status |= (arr.upperBound(0) != dimc0-1)*4;
    status |= (arr.upperBound(1) != dimc1-1)*4;
    status |= (arr.size(0) != dimc0)*4;
    status |= (arr.size(1) != dimc1)*4;
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(2+idim) !=             0)*4;
        status |= (arr.upperBound(2+idim) != dimiv[idim]-1)*4;
        status |= (arr.size(2+idim) != dimiv[idim])*4;
      }
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+2 array failed zero-based, row storage, "
          "continuous component IntVect testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }
#endif  // CH_SPACEDIM<6

//--+1 dimension, (distributed components)

  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+1, ArZeroRow> arr(dimc0, dimiv);
    status |= (arr.size() != szRc0)*1;
    // Test access/storage
    Real c = 0.;
    for (int c0 = 0; c0 != dimc0; ++c0)
      for (BoxIterator bit(Box(IntVect::Zero, dimiv - IntVect::Unit)); bit.ok();
           ++bit)
        {
          arr(c0, bit()) = c;
          c += 0.5;
        }
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test bounds
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(idim) !=             0)*4;
        status |= (arr.upperBound(idim) != dimiv[idim]-1)*4;
        status |= (arr.size(idim) != dimiv[idim])*4;
      }
    status |= (arr.lowerBound(SpaceDim+0) !=       0)*4;
    status |= (arr.upperBound(SpaceDim+0) != dimc0-1)*4;
    status |= (arr.size(SpaceDim+0) != dimc0)*4;
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+1 array failed zero-based, row storage, "
          "distributed component IntVect testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }

//--+2 dimensions, (distributed components)

#if CH_SPACEDIM<6
  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+2, ArZeroRow> arr(dimc1, dimc0, dimiv);
    status |= (arr.size() != szRc1)*1;
    // Test access/storage
    Real c = 0.;
    for (int c1 = 0; c1 != dimc1; ++c1)
      for (int c0 = 0; c0 != dimc0; ++c0)
        for (BoxIterator bit(Box(IntVect::Zero, dimiv - IntVect::Unit));
             bit.ok(); ++bit)
          {
            arr(c1, c0, bit()) = c;
            c += 0.5;
          }
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test bounds
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(idim) !=             0)*4;
        status |= (arr.upperBound(idim) != dimiv[idim]-1)*4;
        status |= (arr.size(idim) != dimiv[idim])*4;
      }
    status |= (arr.lowerBound(SpaceDim+0) !=       0)*4;
    status |= (arr.lowerBound(SpaceDim+1) !=       0)*4;
    status |= (arr.upperBound(SpaceDim+0) != dimc0-1)*4;
    status |= (arr.upperBound(SpaceDim+1) != dimc1-1)*4;
    status |= (arr.size(SpaceDim+0) != dimc0)*4;
    status |= (arr.size(SpaceDim+1) != dimc1)*4;
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+2 array failed zero-based, row storage, "
          "distributed component IntVect testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }
#endif  // CH_SPACEDIM<6

#endif

/*--------------------------------------------------------------------*
 * Column ordered, zero-based subscript testing with IntVects
 *
 * Flags
 *   1 : size
 *   2 : access/storage
 *   4 : storage/const access
 *   8 : assignment/storage
 *   16: bounds
 *   32: fortran
 *--------------------------------------------------------------------*/

//--+1 dimension, (continguous components)

  {
    int status = 0;
    Box box(IntVect::Zero, dimiv - IntVect::Unit);
    // Test constructor
    CHArray<Real, SpaceDim+1, ArZeroCol> arr(dimc0, dimiv);
    status |= (arr.size() != szRc0)*1;
    // Test access/storage
    Real c = 0.;
    for (BoxIterator bit(box); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c0 = 0; c0 != dimc0; ++c0)
          {
            arr(c0, iv) = c;
            c += 0.5;
          }
      }
#ifdef CHARRAYTESTCPP
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test define
    arr.define(dimc0df, dimivdf);
    status |= (arr.size() != szRc0df)*1;
    // Test storage/const access
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        *p = c;
        c += 0.5;
      }
    c = 0.;
    for (BoxIterator bit(Box(IntVect::Zero, dimivdf - IntVect::Unit)); bit.ok();
         ++bit)
      {
        const IntVect iv = bit();
        for (int c0 = 0; c0 != dimc0df; ++c0)
          {
            status |= (arr(c0, iv) != c)*4;
            c += 0.5;
          }
      }
    // Test define on address
    std::memset(chunk, 0, szByteChunk);
    arr.define(chunk, dimc0, dimiv);
    status |= (arr.size() != szRc0)*1;
    // Test assignment/storage
    arr = -0.5;
    c = 0.5;
    for (BoxIterator bit(box); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c0 = 0; c0 != dimc0; ++c0)
          {
            arr(c0, iv) += c;
            c += 0.5;
          }
      }
    c = 0.;
    for (int ic = 0; ic != szRc0; ++ic)
      {
        status |= (chunk[ic] != c)*8;
        c += 0.5;
      }
    // Test bounds
    status |= (arr.lowerBound(0) !=       0)*16;
    status |= (arr.upperBound(0) != dimc0-1)*16;
    status |= (arr.size(0) != dimc0)*16;
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(1+idim) !=             0)*16;
        status |= (arr.upperBound(1+idim) != dimiv[idim]-1)*16;
        status |= (arr.size(1+idim) != dimiv[idim])*16;
      }
#endif
#ifdef CHARRAYTESTCHF
    // Test assignment in Fortran
    CHArray<Real, SpaceDim+1, ArZeroCol> arrf(dimc0, dimiv);
    FORT_RANKCHARRAYSPACEDIMPLUS1CONT(
      CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, arrf),
      CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, arr),
      CHF_BOX(box));  // arrf = arr + inc
    c = 0.;
    Real cinc = 0.5;
    for (Real* p = arrf.begin(); p != arrf.end(); ++p)
      {
        status |= (*p != (c+cinc))*32;
        c += 0.5;
        cinc += 0.5;
      }
#endif
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+1 array failed zero-based, column "
          "storage, continuous component IntVect testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }

//--+2 dimensions, (continguous components)

#if CH_SPACEDIM<6
  {
    int status = 0;
    Box box(IntVect::Zero, dimiv - IntVect::Unit);
    // Test constructor
    CHArray<Real, SpaceDim+2, ArZeroCol> arr(dimc0, dimc1, dimiv);
    status |= (arr.size() != szRc1)*1;
    // Test access/storage
    Real c = 0.;
    for (BoxIterator bit(box); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c1 = 0; c1 != dimc1; ++c1)
          for (int c0 = 0; c0 != dimc0; ++c0)
            {
              arr(c0, c1, iv) = c;
              c += 0.5;
            }
      }
#ifdef CHARRAYTESTCPP
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test define
    arr.define(dimc0df, dimc1df, dimivdf);
    status |= (arr.size() != szRc1df)*1;
    // Test storage/const access
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        *p = c;
        c += 0.5;
      }
    c = 0.;
    for (BoxIterator bit(Box(IntVect::Zero, dimivdf - IntVect::Unit)); bit.ok();
         ++bit)
      {
        const IntVect iv = bit();
        for (int c1 = 0; c1 != dimc1df; ++c1)
          for (int c0 = 0; c0 != dimc0df; ++c0)
            {
              status |= (arr(c0, c1, iv) != c)*4;
              c += 0.5;
            }
      }
    // Test define on address
    std::memset(chunk, 0, szByteChunk);
    arr.define(chunk, dimc0, dimc1, dimiv);
    status |= (arr.size() != szRc1)*1;
    // Test assignment/storage
    arr = -0.5;
    c = 0.5;
    for (BoxIterator bit(box); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c1 = 0; c1 != dimc1; ++c1)
          for (int c0 = 0; c0 != dimc0; ++c0)
            {
              arr(c0, c1, iv) += c;
              c += 0.5;
            }
      }
    c = 0.;
    for (int ic = 0; ic != szRc1; ++ic)
      {
        status |= (chunk[ic] != c)*8;
        c += 0.5;
      }
    // Test bounds
    status |= (arr.lowerBound(0) !=       0)*16;
    status |= (arr.lowerBound(1) !=       0)*16;
    status |= (arr.upperBound(0) != dimc0-1)*16;
    status |= (arr.upperBound(1) != dimc1-1)*16;
    status |= (arr.size(0) != dimc0)*16;
    status |= (arr.size(1) != dimc1)*16;
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(2+idim) !=             0)*16;
        status |= (arr.upperBound(2+idim) != dimiv[idim]-1)*16;
        status |= (arr.size(2+idim) != dimiv[idim])*16;
      }
#endif
#ifdef CHARRAYTESTCHF
    // Test assignment in Fortran
    CHArray<Real, SpaceDim+2, ArZeroCol> arrf(dimc0, dimc1, dimiv);
    FORT_RANKCHARRAYSPACEDIMPLUS2CONT(
      CHF_RCHARRAY(RANK_SPACEDIM_PLUS_2, arrf),
      CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_2, arr),
      CHF_BOX(box));  // arrf = arr + inc
    c = 0.;
    Real cinc = 0.5;
    for (Real* p = arrf.begin(); p != arrf.end(); ++p)
      {
        status |= (*p != (c+cinc))*32;
        c += 0.5;
        cinc += 0.5;
      }
#endif
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+2 array failed zero-based, column "
          "storage, continuous component IntVect testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }
#endif  // CH_SPACEDIM<6

//--+1 dimension, (distributed components)

  {
    int status = 0;
    Box box(IntVect::Zero, dimiv - IntVect::Unit);
    // Test constructor
    CHArray<Real, SpaceDim+1, ArZeroCol> arr(dimiv, dimc0);
    status |= (arr.size() != szRc0)*1;
    // Test access/storage
    Real c = 0.;
    for (int c0 = 0; c0 != dimc0; ++c0)
      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          arr(bit(), c0) = c;
          c += 0.5;
        }
#ifdef CHARRAYTESTCPP
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test define
    arr.define(dimivdf, dimc0df);
    status |= (arr.size() != szRc0df)*1;
    // Test storage/const access
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        *p = c;
        c += 0.5;
      }
    c = 0.;
    for (int c0 = 0; c0 != dimc0df; ++c0)
      for (BoxIterator bit(Box(IntVect::Zero, dimivdf - IntVect::Unit));
           bit.ok(); ++bit)
        {
          status |= (arr(bit(), c0) != c)*4;
          c += 0.5;
        }
    // Test define on address
    std::memset(chunk, 0, szByteChunk);
    arr.define(chunk, dimiv, dimc0);
    status |= (arr.size() != szRc0)*1;
    // Test assignment/storage
    arr = -0.5;
    c = 0.5;
    for (int c0 = 0; c0 != dimc0; ++c0)
      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          arr(bit(), c0) += c;
          c += 0.5;
        }
    c = 0.;
    for (int ic = 0; ic != szRc0; ++ic)
      {
        status |= (chunk[ic] != c)*8;
        c += 0.5;
      }
    // Test bounds
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(idim) !=             0)*16;
        status |= (arr.upperBound(idim) != dimiv[idim]-1)*16;
        status |= (arr.size(idim) != dimiv[idim])*16;
      }
    status |= (arr.lowerBound(SpaceDim+0) !=       0)*16;
    status |= (arr.upperBound(SpaceDim+0) != dimc0-1)*16;
    status |= (arr.size(SpaceDim+0) != dimc0)*16;
#endif
#ifdef CHARRAYTESTCHF
    // Test assignment in Fortran
    CHArray<Real, SpaceDim+1, ArZeroCol> arrf(dimiv, dimc0);
    FORT_RANKCHARRAYSPACEDIMPLUS1DIST(
      CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, arrf),
      CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, arr),
      CHF_BOX(box),
      CHF_CONST_INT(dimc0));  // arrf = arr + inc
    c = 0.;
    Real cinc = 0.5;
    for (Real* p = arrf.begin(); p != arrf.end(); ++p)
      {
        status |= (*p != (c+cinc))*32;
        c += 0.5;
        cinc += 0.5;
      }
#endif
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+1 array failed zero-based, column "
          "storage, distributed component IntVect testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }

//--+2 dimensions, (distributed components)

#if CH_SPACEDIM<6
  {
    int status = 0;
    Box box(IntVect::Zero, dimiv - IntVect::Unit);
    // Test constructor
    CHArray<Real, SpaceDim+2, ArZeroCol> arr(dimiv, dimc0, dimc1);
    status |= (arr.size() != szRc1)*1;
    // Test access/storage
    Real c = 0.;
    for (int c1 = 0; c1 != dimc1; ++c1)
      for (int c0 = 0; c0 != dimc0; ++c0)
        for (BoxIterator bit(box); bit.ok(); ++bit)
            {
              arr(bit(), c0, c1) = c;
              c += 0.5;
            }
#ifdef CHARRAYTESTCPP
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test define
    arr.define(dimivdf, dimc0df, dimc1df);
    status |= (arr.size() != szRc1df)*1;
    // Test storage/const access
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        *p = c;
        c += 0.5;
      }
    c = 0.;
    for (int c1 = 0; c1 != dimc1df; ++c1)
      for (int c0 = 0; c0 != dimc0df; ++c0)
        for (BoxIterator bit(Box(IntVect::Zero, dimivdf - IntVect::Unit));
             bit.ok(); ++bit)
          {
            status |= (arr(bit(), c0, c1) != c)*4;
            c += 0.5;
          }
    // Test define on address
    std::memset(chunk, 0, szByteChunk);
    arr.define(chunk, dimiv, dimc0, dimc1);
    status |= (arr.size() != szRc1)*1;
    // Test assignment/storage
    arr = -0.5;
    c = 0.5;
    for (int c1 = 0; c1 != dimc1; ++c1)
      for (int c0 = 0; c0 != dimc0; ++c0)
        for (BoxIterator bit(box); bit.ok(); ++bit)
          {
            arr(bit(), c0, c1) += c;
            c += 0.5;
          }
    c = 0.;
    for (int ic = 0; ic != szRc1; ++ic)
      {
        status |= (chunk[ic] != c)*8;
        c += 0.5;
      }
    // Test bounds
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(idim) !=             0)*16;
        status |= (arr.upperBound(idim) != dimiv[idim]-1)*16;
        status |= (arr.size(idim) != dimiv[idim])*16;
      }
    status |= (arr.lowerBound(SpaceDim+0) !=       0)*16;
    status |= (arr.lowerBound(SpaceDim+1) !=       0)*16;
    status |= (arr.upperBound(SpaceDim+0) != dimc0-1)*16;
    status |= (arr.upperBound(SpaceDim+1) != dimc1-1)*16;
    status |= (arr.size(SpaceDim+0) != dimc0)*16;
    status |= (arr.size(SpaceDim+1) != dimc1)*16;
#endif
#ifdef CHARRAYTESTCHF
    // Test assignment in Fortran
    CHArray<Real, SpaceDim+2, ArZeroCol> arrf(dimiv, dimc0, dimc1);
    FORT_RANKCHARRAYSPACEDIMPLUS2DIST(
      CHF_RCHARRAY(RANK_SPACEDIM_PLUS_2, arrf),
      CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_2, arr),
      CHF_BOX(box),
      CHF_CONST_INT(dimc0),
      CHF_CONST_INT(dimc1));  // arrf = arr + inc
    c = 0.;
    Real cinc = 0.5;
    for (Real* p = arrf.begin(); p != arrf.end(); ++p)
      {
        status |= (*p != (c+cinc))*32;
        c += 0.5;
        cinc += 0.5;
      }
#endif
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+2 array failed zero-based, column "
          "storage, distributed component IntVect testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }
#endif  // CH_SPACEDIM<6

/*--------------------------------------------------------------------*
 * Row ordered, zero-based subscript testing with Boxes.
 *
 * Only test constructor/access/storage to test Box macros.
 * Remaining features are tested using column-ordered storage
 *
 * Flags
 *   1 : size
 *   2 : access/storage
 *   4 : bounds
 *--------------------------------------------------------------------*/

#ifdef CHARRAYTESTCPP

//--+1 dimension, (continguous components)

  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+1, ArRangeRow> arr(dimbox, dimc0);
    status |= (arr.size() != szRc0)*1;
    // Test access/storage
    Real c = 0.;
    for (BoxIterator bit(dimbox); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c0 = 0; c0 != dimc0; ++c0)
          {
            arr(iv, c0) = c;
            c += 0.5;
          }
      }
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test bounds
    status |= (arr.lowerBound(0) !=       0)*4;
    status |= (arr.upperBound(0) != dimc0-1)*4;
    status |= (arr.size(0) != dimc0)*4;
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(1+idim) != dimbox.smallEnd(idim))*16;
        status |= (arr.upperBound(1+idim) !=   dimbox.bigEnd(idim))*16;
        status |= (arr.size(1+idim) != dimbox.size(idim))*16;
      }
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+1 array failed zero-based, row storage, "
          "continuous component Box testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }

//--+2 dimensions, (continguous components)

#if CH_SPACEDIM<6
  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+2, ArRangeRow> arr(dimbox, dimc1, dimc0);
    status |= (arr.size() != szRc1)*1;
    // Test access/storage
    Real c = 0.;
    for (BoxIterator bit(dimbox); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c1 = 0; c1 != dimc1; ++c1)
          for (int c0 = 0; c0 != dimc0; ++c0)
            {
              arr(iv, c1, c0) = c;
              c += 0.5;
            }
      }
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test bounds
    status |= (arr.lowerBound(0) !=       0)*4;
    status |= (arr.lowerBound(1) !=       0)*4;
    status |= (arr.upperBound(0) != dimc0-1)*4;
    status |= (arr.upperBound(1) != dimc1-1)*4;
    status |= (arr.size(0) != dimc0)*4;
    status |= (arr.size(1) != dimc1)*4;
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(2+idim) != dimbox.smallEnd(idim))*16;
        status |= (arr.upperBound(2+idim) !=   dimbox.bigEnd(idim))*16;
        status |= (arr.size(2+idim) != dimbox.size(idim))*16;
      }
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+2 array failed zero-based, row storage, "
          "continuous component Box testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }
#endif  // CH_SPACEDIM<6

//--+1 dimension, (distributed components)

  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+1, ArRangeRow> arr(dimc0, dimbox);
    status |= (arr.size() != szRc0)*1;
    // Test access/storage
    Real c = 0.;
    for (int c0 = 0; c0 != dimc0; ++c0)
      for (BoxIterator bit(dimbox); bit.ok(); ++bit)
        {
          arr(c0, bit()) = c;
          c += 0.5;
        }
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test bounds
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(idim) != dimbox.smallEnd(idim))*16;
        status |= (arr.upperBound(idim) !=   dimbox.bigEnd(idim))*16;
        status |= (arr.size(idim) != dimbox.size(idim))*16;
      }
    status |= (arr.lowerBound(SpaceDim+0) !=       0)*4;
    status |= (arr.upperBound(SpaceDim+0) != dimc0-1)*4;
    status |= (arr.size(SpaceDim+0) != dimc0)*4;
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+1 array failed zero-based, row storage, "
          "distributed component Box testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }

//--+2 dimensions, (distributed components)

#if CH_SPACEDIM<6
  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+2, ArRangeRow> arr(dimc1, dimc0, dimbox);
    status |= (arr.size() != szRc1)*1;
    // Test access/storage
    Real c = 0.;
    for (int c1 = 0; c1 != dimc1; ++c1)
      for (int c0 = 0; c0 != dimc0; ++c0)
        for (BoxIterator bit(dimbox); bit.ok(); ++bit)
          {
            arr(c1, c0, bit()) = c;
            c += 0.5;
          }
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test bounds
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(idim) != dimbox.smallEnd(idim))*16;
        status |= (arr.upperBound(idim) !=   dimbox.bigEnd(idim))*16;
        status |= (arr.size(idim) != dimbox.size(idim))*16;
      }
    status |= (arr.lowerBound(SpaceDim+0) !=       0)*4;
    status |= (arr.lowerBound(SpaceDim+1) !=       0)*4;
    status |= (arr.upperBound(SpaceDim+0) != dimc0-1)*4;
    status |= (arr.upperBound(SpaceDim+1) != dimc1-1)*4;
    status |= (arr.size(SpaceDim+0) != dimc0)*4;
    status |= (arr.size(SpaceDim+1) != dimc1)*4;
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+2 array failed zero-based, row storage, "
          "distributed component Box testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }
#endif  // CH_SPACEDIM<6

#endif

/*--------------------------------------------------------------------*
 * Column ordered, zero-based subscript testing with Boxes
 *
 * Flags
 *   1 : size
 *   2 : access/storage
 *   4 : storage/const access
 *   8 : assignment/storage
 *   16: bounds
 *   32: fortran
 *--------------------------------------------------------------------*/

//--+1 dimension, (continguous components)

  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+1, ArRangeCol> arr(dimc0, dimbox);
    status |= (arr.size() != szRc0)*1;
    // Test access/storage
    Real c = 0.;
    for (BoxIterator bit(dimbox); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c0 = 0; c0 != dimc0; ++c0)
          {
            arr(c0, iv) = c;
            c += 0.5;
          }
      }
#ifdef CHARRAYTESTCPP
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test define
    arr.define(dimc0df, dimboxdf);
    status |= (arr.size() != szRc0df)*1;
    // Test storage/const access
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        *p = c;
        c += 0.5;
      }
    c = 0.;
    for (BoxIterator bit(dimboxdf); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c0 = 0; c0 != dimc0df; ++c0)
          {
            status |= (arr(c0, iv) != c)*4;
            c += 0.5;
          }
      }
    // Test define on address
    std::memset(chunk, 0, szByteChunk);
    arr.define(chunk, dimc0, dimbox);
    status |= (arr.size() != szRc0)*1;
    // Test assignment/storage
    arr = -0.5;
    c = 0.5;
    for (BoxIterator bit(dimbox); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c0 = 0; c0 != dimc0; ++c0)
          {
            arr(c0, iv) += c;
            c += 0.5;
          }
      }
    c = 0.;
    for (int ic = 0; ic != szRc0; ++ic)
      {
        status |= (chunk[ic] != c)*8;
        c += 0.5;
      }
    // Test bounds
    status |= (arr.lowerBound(0) !=       0)*16;
    status |= (arr.upperBound(0) != dimc0-1)*16;
    status |= (arr.size(0) != dimc0)*16;
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(1+idim) != dimbox.smallEnd(idim))*16;
        status |= (arr.upperBound(1+idim) !=   dimbox.bigEnd(idim))*16;
        status |= (arr.size(1+idim) != dimbox.size(idim))*16;
      }
#endif
#ifdef CHARRAYTESTCHF
    // Test assignment in Fortran
    CHArray<Real, SpaceDim+1, ArRangeCol> arrf(dimc0, dimbox);
    FORT_RANKCHARRAYSPACEDIMPLUS1CONT(
      CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, arrf),
      CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, arr),
      CHF_BOX(dimbox));  // arrf = arr + inc
    c = 0.;
    Real cinc = 0.5;
    for (Real* p = arrf.begin(); p != arrf.end(); ++p)
      {
        status |= (*p != (c+cinc))*32;
        c += 0.5;
        cinc += 0.5;
      }
#endif
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+1 array failed zero-based, column "
          "storage, continuous component, Box testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }

//--+2 dimensions, (continguous components)

#if CH_SPACEDIM<6
  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+2, ArRangeCol> arr(dimc0, dimc1, dimbox);
    status |= (arr.size() != szRc1)*1;
    // Test access/storage
    Real c = 0.;
    for (BoxIterator bit(dimbox); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c1 = 0; c1 != dimc1; ++c1)
          for (int c0 = 0; c0 != dimc0; ++c0)
            {
              arr(c0, c1, iv) = c;
              c += 0.5;
            }
      }
#ifdef CHARRAYTESTCPP
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test define
    arr.define(dimc0df, dimc1df, dimboxdf);
    status |= (arr.size() != szRc1df)*1;
    // Test storage/const access
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        *p = c;
        c += 0.5;
      }
    c = 0.;
    for (BoxIterator bit(dimboxdf); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c1 = 0; c1 != dimc1df; ++c1)
          for (int c0 = 0; c0 != dimc0df; ++c0)
            {
              status |= (arr(c0, c1, iv) != c)*4;
              c += 0.5;
            }
      }
    // Test define on address
    std::memset(chunk, 0, szByteChunk);
    arr.define(chunk, dimc0, dimc1, dimbox);
    status |= (arr.size() != szRc1)*1;
    // Test assignment/storage
    arr = -0.5;
    c = 0.5;
    for (BoxIterator bit(dimbox); bit.ok(); ++bit)
      {
        const IntVect iv = bit();
        for (int c1 = 0; c1 != dimc1; ++c1)
          for (int c0 = 0; c0 != dimc0; ++c0)
            {
              arr(c0, c1, iv) += c;
              c += 0.5;
            }
      }
    c = 0.;
    for (int ic = 0; ic != szRc1; ++ic)
      {
        status |= (chunk[ic] != c)*8;
        c += 0.5;
      }
    // Test bounds
    status |= (arr.lowerBound(0) !=       0)*16;
    status |= (arr.lowerBound(1) !=       0)*16;
    status |= (arr.upperBound(0) != dimc0-1)*16;
    status |= (arr.upperBound(1) != dimc1-1)*16;
    status |= (arr.size(0) != dimc0)*16;
    status |= (arr.size(1) != dimc1)*16;
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(2+idim) != dimbox.smallEnd(idim))*16;
        status |= (arr.upperBound(2+idim) !=   dimbox.bigEnd(idim))*16;
        status |= (arr.size(2+idim) != dimbox.size(idim))*16;
      }
#endif
#ifdef CHARRAYTESTCHF
    // Test assignment in Fortran
    CHArray<Real, SpaceDim+2, ArRangeCol> arrf(dimc0, dimc1, dimbox);
    FORT_RANKCHARRAYSPACEDIMPLUS2CONT(
      CHF_RCHARRAY(RANK_SPACEDIM_PLUS_2, arrf),
      CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_2, arr),
      CHF_BOX(dimbox));  // arrf = arr + inc
    c = 0.;
    Real cinc = 0.5;
    for (Real* p = arrf.begin(); p != arrf.end(); ++p)
      {
        status |= (*p != (c+cinc))*32;
        c += 0.5;
        cinc += 0.5;
      }
#endif
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+2 array failed zero-based, column "
          "storage, continuous component Box testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }
#endif  // CH_SPACEDIM<6

//--+1 dimension, (distributed components)

  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+1, ArRangeCol> arr(dimbox, dimc0);
    status |= (arr.size() != szRc0)*1;
    // Test access/storage
    Real c = 0.;
    for (int c0 = 0; c0 != dimc0; ++c0)
      for (BoxIterator bit(dimbox); bit.ok(); ++bit)
        {
          arr(bit(), c0) = c;
          c += 0.5;
        }
#ifdef CHARRAYTESTCPP
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test define
    arr.define(dimboxdf, dimc0df);
    status |= (arr.size() != szRc0df)*1;
    // Test storage/const access
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        *p = c;
        c += 0.5;
      }
    c = 0.;
    for (int c0 = 0; c0 != dimc0df; ++c0)
      for (BoxIterator bit(dimboxdf); bit.ok(); ++bit)
        {
          status |= (arr(bit(), c0) != c)*4;
          c += 0.5;
        }
    // Test define on address
    std::memset(chunk, 0, szByteChunk);
    arr.define(chunk, dimbox, dimc0);
    status |= (arr.size() != szRc0)*1;
    // Test assignment/storage
    arr = -0.5;
    c = 0.5;
    for (int c0 = 0; c0 != dimc0; ++c0)
      for (BoxIterator bit(dimbox); bit.ok(); ++bit)
        {
          arr(bit(), c0) += c;
          c += 0.5;
        }
    c = 0.;
    for (int ic = 0; ic != szRc0; ++ic)
      {
        status |= (chunk[ic] != c)*8;
        c += 0.5;
      }
    // Test bounds
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(idim) != dimbox.smallEnd(idim))*16;
        status |= (arr.upperBound(idim) !=   dimbox.bigEnd(idim))*16;
        status |= (arr.size(idim) != dimbox.size(idim))*16;
      }
    status |= (arr.lowerBound(SpaceDim+0) !=       0)*16;
    status |= (arr.upperBound(SpaceDim+0) != dimc0-1)*16;
    status |= (arr.size(SpaceDim+0) != dimc0)*16;
#endif
#ifdef CHARRAYTESTCHF
    // Test assignment in Fortran
    CHArray<Real, SpaceDim+1, ArRangeCol> arrf(dimbox, dimc0);
    FORT_RANKCHARRAYSPACEDIMPLUS1DIST(
      CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, arrf),
      CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, arr),
      CHF_BOX(dimbox),
      CHF_CONST_INT(dimc0));  // arrf = arr + inc
    c = 0.;
    Real cinc = 0.5;
    for (Real* p = arrf.begin(); p != arrf.end(); ++p)
      {
        status |= (*p != (c+cinc))*32;
        c += 0.5;
        cinc += 0.5;
      }
#endif
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+1 array failed zero-based, column "
          "storage, distributed component Box testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }

//--+2 dimensions, (distributed components)

#if CH_SPACEDIM<6
  {
    int status = 0;
    // Test constructor
    CHArray<Real, SpaceDim+2, ArRangeCol> arr(dimbox, dimc0, dimc1);
    status |= (arr.size() != szRc1)*1;
    // Test access/storage
    Real c = 0.;
    for (int c1 = 0; c1 != dimc1; ++c1)
      for (int c0 = 0; c0 != dimc0; ++c0)
        for (BoxIterator bit(dimbox); bit.ok(); ++bit)
            {
              arr(bit(), c0, c1) = c;
              c += 0.5;
            }
#ifdef CHARRAYTESTCPP
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        status |= (*p != c)*2;
        c += 0.5;
      }
    // Test define
    arr.define(dimboxdf, dimc0df, dimc1df);
    status |= (arr.size() != szRc1df)*1;
    // Test storage/const access
    c = 0.;
    for (Real* p = arr.begin(); p != arr.end(); ++p)
      {
        *p = c;
        c += 0.5;
      }
    c = 0.;
    for (int c1 = 0; c1 != dimc1df; ++c1)
      for (int c0 = 0; c0 != dimc0df; ++c0)
        for (BoxIterator bit(dimboxdf); bit.ok(); ++bit)
          {
            status |= (arr(bit(), c0, c1) != c)*4;
            c += 0.5;
          }
    // Test define on address
    std::memset(chunk, 0, szByteChunk);
    arr.define(chunk, dimbox, dimc0, dimc1);
    status |= (arr.size() != szRc1)*1;
    // Test assignment/storage
    arr = -0.5;
    c = 0.5;
    for (int c1 = 0; c1 != dimc1; ++c1)
      for (int c0 = 0; c0 != dimc0; ++c0)
        for (BoxIterator bit(dimbox); bit.ok(); ++bit)
          {
            arr(bit(), c0, c1) += c;
            c += 0.5;
          }
    c = 0.;
    for (int ic = 0; ic != szRc1; ++ic)
      {
        status |= (chunk[ic] != c)*8;
        c += 0.5;
      }
    // Test bounds
    for (int idim = 0; idim != SpaceDim; ++idim)
      {
        status |= (arr.lowerBound(idim) != dimbox.smallEnd(idim))*16;
        status |= (arr.upperBound(idim) !=   dimbox.bigEnd(idim))*16;
        status |= (arr.size(idim) != dimbox.size(idim))*16;
      }
    status |= (arr.lowerBound(SpaceDim+0) !=       0)*16;
    status |= (arr.lowerBound(SpaceDim+1) !=       0)*16;
    status |= (arr.upperBound(SpaceDim+0) != dimc0-1)*16;
    status |= (arr.upperBound(SpaceDim+1) != dimc1-1)*16;
    status |= (arr.size(SpaceDim+0) != dimc0)*16;
    status |= (arr.size(SpaceDim+1) != dimc1)*16;
#endif
#ifdef CHARRAYTESTCHF
    // Test assignment in Fortran
    CHArray<Real, SpaceDim+2, ArRangeCol> arrf(dimbox, dimc0, dimc1);
    FORT_RANKCHARRAYSPACEDIMPLUS2DIST(
      CHF_RCHARRAY(RANK_SPACEDIM_PLUS_2, arrf),
      CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_2, arr),
      CHF_BOX(dimbox),
      CHF_CONST_INT(dimc0),
      CHF_CONST_INT(dimc1));  // arrf = arr + inc
    c = 0.;
    Real cinc = 0.5;
    for (Real* p = arrf.begin(); p != arrf.end(); ++p)
      {
        status |= (*p != (c+cinc))*32;
        c += 0.5;
        cinc += 0.5;
      }
#endif
    if (status)
      {
        pout() << indent << pgmname
               << ": Rank SpaceDim+2 array failed zero-based, column "
          "storage, distributed component Box testing"
               << endl;
        pout() << indent2
               << "flags: " << status << endl;
        ++errors;
      }
  }
#endif  // CH_SPACEDIM<6

  return errors;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i] ,"-v" ,3 ) == 0)
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if (strncmp( argv[i] ,"-q" ,3 ) == 0)
            {
              verbose = false ;
              // argv[i] = "" ;
            }
          else
            {
              break ;
            }
        }
    }
  return ;
}
