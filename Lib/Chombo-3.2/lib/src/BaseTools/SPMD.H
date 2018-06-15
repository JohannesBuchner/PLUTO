#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// These are the inevitable functions that people
// can't avoid using when writing a SPMD program.
// It is a minimal set from which more elaborate
// functionality can be generated.  As always, a
// user is free to utilize the entire MPI programming
// on their own platform.  The functions are
// assured to work on all machines supported.
//

#ifndef _SPMD_H_
#define _SPMD_H_

#include "REAL.H"
#include <string>
#include "Vector.H"
#include "MayDay.H"
#include "Misc.H"

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "BaseNamespaceHeader.H"

#ifdef CH_MPI
struct Chombo_MPI
{
  static MPI_Comm comm;
};

#else
// this can be changed for debugging parallel code in serial
extern int num_procs ;
#endif

extern long long CH_MAX_MPI_MESSAGE_SIZE;

extern long long CH_MaxMPISendSize;
extern long long CH_MaxMPIRecvSize;

int reportMPIStats();

/// local process ID
/**
    Returns the ID of the locally running process in the range 0 <=
    procID() < numProc().  This has no relation to the operating system
    pid.  There is always a procID() == 0.  */
int procID();

inline int CHprocID()
{
  return procID();
}

/// number of parallel processes
/**
    Returns the number of parallel processes running.
    Always returns at least 1. */
unsigned int numProc();

///all ranks wait here to sync-up
/**
   All MPI ranks wait here to sync-up.  Calls MPI_Barrier(comm)  */
void barrier(void);

template <class T>
int linearSize(const T& inputT);

template <class T>
void linearIn(T& a_outputT, const void* const inBuf);

template <class T>
void linearOut(void* const a_outBuf, const T& inputT);

#ifdef CH_MPI
///gather from multiple processes
/**
  Gathers a_input of every process onto Vector<T> a_outVec on process
  number a_dest.  It is an error if a_dest < 0 or a_dest >= numProc().
  a_outVec is numProc() long with the value of a_input on every
  process in its elements.\\

T must have have the following functions:
  \\

  int linearSize<T>(const T& inputT); \\
  void linearIn<T>(T& a_outputT, const void* const inBuf); \\
  void linearOut<T>(void* a_outBuf, const T& inputT); \\

*/
template <class T>
void gather(Vector<T>& a_outVec, const T& a_input, int a_dest);

// this has to be here so that linearSize<IntVectSet>, linearIn<IntVectSet>,
//  and linearOut<IntVectSet> are defined.
// template < >
// void gather(Vector<IntVectSet>& a_outVec,
//             const IntVectSet& a_input,
//             int a_dest);
#endif

///broadcast to every process
/**
  Broadcasts a_inAndOut to every process from process number a_src.
  It is an error if a_src < 0 or a_src >= numProc().  \\

  T must have have the following functions: \\

  int linearSize<T>(const T& inputT); \\
  void linearIn<T>(T& a_outputT, const void* const inBuf); \\
  void linearOut<T>(void* a_outBuf, const T& inputT); \\ */
template <class T>
void broadcast(T& a_inAndOut,  int a_src);

//integer specialization of linearSize
template < >
int linearSize(const int& a_input);

//u l l specialization of linearSize
template < >
int linearSize(const unsigned long long& a_input);

//integer specialization of linearIn
template < >
void linearIn(int& a_outputT, const void* const inBuf);

//u l l specialization of linearIn
template < >
void linearIn(unsigned long long& a_outputT, const void* const inBuf);

//integer specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const int& a_inputT);

//unsigned long long specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const unsigned long long& a_inputT);

//long integer specialization of linearSize
template < >
int linearSize(const long& a_input);

//long integer specialization of linearIn
template < >
void linearIn(long& a_outputT, const void* const inBuf);

//long integer specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const long& a_inputT);

//unsigned long integer specialization of linearSize
template < >
int linearSize(const unsigned long& a_input);

//unsigned long integer specialization of linearIn
template < >
void linearIn(unsigned long& a_outputT, const void* const inBuf);

//unsigned long integer specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const unsigned long& a_inputT);

//Real specialization of linearSize
template < >
int linearSize(const float& a_input);

template < >
int linearSize(const double& a_input);

//Real specialization of linearIn
template < >
void linearIn(float& a_outputT, const void* const a_inBuf);

template < >
void linearIn(double& a_outputT, const void* const a_inBuf);

//Real specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const float& a_inputT);

template < >
void linearOut(void* const a_outBuf, const double& a_inputT);

// std::string specialization.
template <>
int linearSize(const std::string& a_input);
template <>
void linearIn(std::string& a_outputT, const void* const a_inBuf);
template <>
void linearOut(void* const a_outBuf, const std::string& a_inputT);

//Vector<int>  specialization
template < >
int linearSize(const Vector<int>& a_input);
template < >
void linearIn(Vector<int>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<int>& a_inputT);

//Vector<unsigned long long>  specialization
template < >
int linearSize(const Vector<unsigned long long>& a_input);
template < >
void linearIn(Vector<unsigned long long>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<unsigned long long>& a_inputT);

//Vector<long>  specialization
template < >
int linearSize(const Vector<long>& a_input);
template < >
void linearIn(Vector<long>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<long>& a_inputT);

//Vector<Real>  specialization
template < >
int linearSize(const Vector<float>& a_input);
template < >
void linearIn(Vector<float>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<float>& a_inputT);

template < >
int linearSize(const Vector<double>& a_input);
template < >
void linearIn(Vector<double>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<double>& a_inputT);

//Vector<std::string>  specialization
template < >
int linearSize(const Vector<std::string>& a_input);
template < >
void linearIn(Vector<std::string>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<std::string>& a_inputT);

//Vector<Vector<int> > specialization
template < >
int linearSize(const Vector<Vector<int> >& a_input);
template < >
void linearIn(Vector<Vector<int> >& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<Vector<int> >& a_inputT);

//Vector<T> specialization of linearSize
template <class T>
int linearListSize(const Vector<T>& a_input);

//Vector<T> specialization of linearIn
template <class T>
void linearListIn(Vector<T>& a_outputT, const void* const a_inBuf);

//Vector<T> specialization of linearOut
template <class T>
void linearListOut(void* const a_outBuf, const Vector<T>& a_inputT);

class SerialTask
{
public:
    enum task
    {
      compute=0
    };
};

int GetPID(int rank);
int GetRank(int pid);

// return id of unique processor for special serial tasks
int
uniqueProc(const SerialTask::task& a_task);

#include "BaseNamespaceFooter.H"

#include "SPMDI.H"

#endif
