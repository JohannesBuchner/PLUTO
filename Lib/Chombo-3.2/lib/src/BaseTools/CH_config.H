#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_CONFIG_H_
#define _CH_CONFIG_H_

/*--------------------------------------------------------------------*
 * System characteristics
 *--------------------------------------------------------------------*/

//--Methods for allocating aligned memory

// Defined if memory allocation can be aligned
#undef CHDEF_SYSTEM_HAVE_ALIGNEDALLOC
// Defined if posix_memalign is used to allocate aligned memory
#undef CHDEF_SYSTEM_HAVE_POSIXMEMALIGN

//--Cache line size

// Line size reported in bytes (the default guess is 64 and we only know it for
// sure if a method is also defined)
#define CHDEF_SYSTEM_CACHELINESIZE 64
// Defined if cache information determined through a system file
#undef CHDEF_SYSTEM_HAVE_CACHEINFOFILE
// Defined if cache information determined through a call to sysconf
#undef CHDEF_SYSTEM_HAVE_CACHESYSCONF

//--X86 vector extensions

// Warning: Your program possibly won't compile unless the right instruction
// sets are enabled for the compiler (for GNU and Intel, these usually look
// like -msse, -msse2, -msse3 -mssse3, -msse4.1, -msse4.2, -msse4a, and -mavx)
// We plan to make this automatic later.  Default settings for most compilers
// on modern (i.e., 64 bit) machines are usually around sse3 (as of 2011).

// Available vector extensions are indicated by the bits set in
// 'CHDEF_SYSTEM_X86VECEXT_CPU_BITS' and 'CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS'.
// The former indicates what the CPU supports and the latter which extensions
// have intrinsics defined that work with the compiler/system
#define CHDEF_BIT_SSE    (1 << 0)
#define CHDEF_BIT_SSE2   (1 << 1)
#define CHDEF_BIT_SSE3   (1 << 2)
#define CHDEF_BIT_SSSE3  (1 << 3)
#define CHDEF_BIT_SSE4D1 (1 << 4)
#define CHDEF_BIT_SSE4D2 (1 << 5)
#define CHDEF_BIT_SSE4A  (1 << 6)
#define CHDEF_BIT_AVX    (1 << 7)
#define CHDEF_BIT_FMA    (1 << 8)
#define CHDEF_BIT_AVX2   (1 << 9)
#define CHDEF_SYSTEM_X86VECEXT_CPU_BITS 0
#define CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS 0

// For each extension that works with the compiler, these are the header files
// to use.
#undef CHDEF_SYSTEM_X86VECEXT_SSE_H
#undef CHDEF_SYSTEM_X86VECEXT_SSE2_H
#undef CHDEF_SYSTEM_X86VECEXT_SSE3_H
#undef CHDEF_SYSTEM_X86VECEXT_SSSE3_H
#undef CHDEF_SYSTEM_X86VECEXT_SSE4D1_H
#undef CHDEF_SYSTEM_X86VECEXT_SSE4D2_H
#undef CHDEF_SYSTEM_X86VECEXT_SSE4A_H
#undef CHDEF_SYSTEM_X86VECEXT_AVX_H
#undef CHDEF_SYSTEM_X86VECEXT_FMA_H
#undef CHDEF_SYSTEM_X86VECEXT_AVX2_H
// This header, if defined, can be used generally for all extensions supported
// by Intel (basically everything except sse4a).  It's usually <immintrin.h>.
#undef CHDEF_SYSTEM_X86VECEXT_INTEL_H

// Defined if operators were found to work ('+' is tested for all of sse, sse2,
// avx, and avx2 that are available as indicated by the compiler bits)
#undef CHDEF_SYSTEM_X86VECEXT_HAVE_OPERATOR

//--How to do a popcnt

// Warning: Your program possibly won't compile unless the right instruction
// sets are enabled for the compiler (for GNU and Intel, this is -mpopcnt)
// We plan to make this automatic later.

// Type of method 0 = software, 1 = GNU x86 assembly
#define CHDEF_SYSTEM_POPCNT_METHOD 0

//--CUDA capable GPU

// Information about CUDA and GPU installation
// The driver version is parsed as DV/1000 << '.' << (DV%100)/10
#undef CHDEF_SYSTEM_CUDA_DRIVER_VERSION
#undef CHDEF_SYSTEM_CUDA_DEVICE_NAME
#undef CHDEF_SYSTEM_CUDA_COMPUTE_CAPABILITY
// GPU CUDA attributes (all memory is in bytes)
#undef CHDEF_SYSTEM_CUDAATTR_TOTAL_GLOBAL_MEMORY
#undef CHDEF_SYSTEM_CUDAATTR_MULTIPROCESSOR_COUNT
#undef CHDEF_SYSTEM_CUDAATTR_TOTAL_CONSTANT_MEMORY
#undef CHDEF_SYSTEM_CUDAATTR_MAX_SHARED_MEMORY_PER_BLOCK
#undef CHDEF_SYSTEM_CUDAATTR_MAX_REGISTERS_PER_BLOCK
#undef CHDEF_SYSTEM_CUDAATTR_WARP_SIZE
#undef CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_MULTIPROCESSOR
#undef CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_BLOCK
#undef CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_X
#undef CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_Y
#undef CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_Z
#undef CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_X
#undef CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_Y
#undef CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_Z

#endif
