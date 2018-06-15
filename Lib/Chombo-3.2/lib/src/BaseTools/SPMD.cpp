#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <cstring>
// #extern "C" {      // The #extern "C" might have been here for a reason...
#include <unistd.h>
// }

#include "SPMD.H"
#include "parstream.H"
#include "BaseNamespaceHeader.H"

using std::cout;
using std::endl;

// try a 30 Mbyte max message size and see if that helps.

long long CH_MAX_MPI_MESSAGE_SIZE = 30*1024*1024;
long long CH_MaxMPISendSize = 0;
long long CH_MaxMPIRecvSize  = 0;

int reportMPIStats()
{
  pout()<<"Chombo message size limit:"<< CH_MAX_MPI_MESSAGE_SIZE<<"\n"
        <<"Max send message size:"<<CH_MaxMPISendSize<<"\n"
        <<"Max recv message size:"<<CH_MaxMPIRecvSize<<std::endl;
  return 0;
}

Vector<int> pids;

#ifndef CH_MPI

int procID()
{
  return 0;
}

// reset this to fool serial code into thinking its parallel
int num_procs = 1 ;

int GetRank(int pid)
{
  return 0;
}

int GetPID(int rank)
{
  CH_assert(rank == 0);
  return getpid();
}

unsigned int numProc()
{
  return num_procs;
}

#else // CH_MPI version

int GetPID(int rank)
{
  if (rank<0) return -1;
  if (rank>=pids.size()) return -2;
  return pids[rank];
}

int GetRank(int pid)
{
  for (int i=0; i<pids.size(); ++i)
    {
      if (pids[i]== pid) return i;
    }
  return -1;
}

// this is a global variable which indicates whether
// or not an MPI command should be used to deduce rank.
// needed for applications which switch communicators.
// set g_resetProcID=true to force next procID() call to 
// querry MPI_Comm_rank
bool g_resetProcID;

int procID()
{
  static bool firstCall = true;
  static int lastProcID = 0;
  if (firstCall || g_resetProcID )
  {
    g_resetProcID = false;
    firstCall = false;

//     int proc = getpid();
//     gather(pids, proc, 0);
//     broadcast(pids, 0);
    MPI_Comm_rank(Chombo_MPI::comm, &lastProcID);
  }
  return lastProcID;
}

unsigned int numProc()
{
  static int ret = -1;
  if (ret == -1)
  {
    MPI_Comm_size(Chombo_MPI::comm, &ret);
  }
  return ret;
}

// hopefully static copy of opaque handles
MPI_Comm Chombo_MPI::comm = MPI_COMM_WORLD;

#endif // CH_MPI

/// these should work independent of MPI's existence
template < >
void linearIn<float>(float& a_outputT, const void* const a_inBuf)
{
  //this fandango is to avoid unaligned accesses
  char realBuf[sizeof(float)];
  memcpy(realBuf, a_inBuf, sizeof(float));
  float* buffer = (float*)realBuf;
  a_outputT = *buffer;
}

template < >
void linearOut<float>(void* const a_outBuf, const float& a_inputT)
{
  //this fandango is to avoid unaligned accesses
  char realBuf[sizeof(float)];
  float* realPtr = (float*)realBuf;
  *realPtr = a_inputT;
  memcpy(a_outBuf, realBuf, sizeof(float));
}

template < >
int linearSize<float>(const float& inputfloat)
{
  return sizeof(float);
}

template < >
void linearIn<double>(double& a_outputT, const void* const a_inBuf)
{
  //this fandango is to avoid unaligned accesses
  char realBuf[sizeof(double)];
  memcpy(realBuf, a_inBuf, sizeof(double));
  double* buffer = (double*)realBuf;
  a_outputT = *buffer;
}

template < >
void linearOut<double>(void* const a_outBuf, const double& a_inputT)
{
  //this fandango is to avoid unaligned accesses
  char realBuf[sizeof(double)];
  double* realPtr = (double*)realBuf;
  *realPtr = a_inputT;
  memcpy(a_outBuf, realBuf, sizeof(double));
}

template < >
int linearSize<double>(const double& inputdouble)
{
  return sizeof(double);
}

template < >
void linearIn<int>(int& a_outputT, const void* const a_inBuf)
{
  int* buffer = (int*)a_inBuf;
  a_outputT = *buffer;
}

template < >
void linearIn<unsigned long long>(unsigned long long& a_outputT, const void* const a_inBuf)
{
  unsigned long long* buffer = (unsigned long long*)a_inBuf;
  a_outputT = *buffer;
}

template < >
void linearOut<int>(void* const a_outBuf, const int& a_inputT)
{
  int* buffer = (int*)a_outBuf;
  *buffer = a_inputT;
}

template < >
void linearOut<unsigned long long>(void* const a_outBuf, const unsigned long long& a_inputT)
{
  unsigned long long* buffer = (unsigned long long*)a_outBuf;
  *buffer = a_inputT;
}

template < >
int linearSize<int>(const int& a_input)
{
  return sizeof(int);
}

template < >
int linearSize<unsigned long long>(const unsigned long long& a_input)
{
  return sizeof(unsigned long long);
}

template < >
void linearIn<long>(long& a_outputT, const void* const a_inBuf)
{
  long* buffer = (long*)a_inBuf;
  a_outputT = *buffer;
}

template < >
void linearOut<long>(void* const a_outBuf, const long& a_inputT)
{
  long* buffer = (long*)a_outBuf;
  *buffer = a_inputT;
}

template < >
int linearSize<long>(const long& a_input)
{
  return sizeof(long);
}

template < >
void linearIn<unsigned long>(unsigned long& a_outputT, const void* const a_inBuf)
{
  unsigned long* buffer = (unsigned long*)a_inBuf;
  a_outputT = *buffer;
}

template < >
void linearOut<unsigned long>(void* const a_outBuf, const unsigned long& a_inputT)
{
  unsigned long* buffer = (unsigned long*)a_outBuf;
  *buffer = a_inputT;
}

template < >
int linearSize<unsigned long>(const unsigned long& a_input)
{
  return sizeof(unsigned long);
}

// std::string specialization.
template <>
void linearIn<std::string>(std::string& a_outputT, const void* const a_inBuf)
{
  int* intBuffer = (int*)a_inBuf;
  int length = intBuffer[0];
  if (length > 0)
  {
    const char* charBuffer = (const char*)(&intBuffer[1]);
    a_outputT.assign(charBuffer, length);
  }
  else a_outputT = "";
}
template <>
void linearOut<std::string>(void* const a_outBuf, const std::string& a_inputT)
{
  int* intBuffer = (int*)a_outBuf;
  intBuffer[0] = a_inputT.length();
  if (a_inputT.length() > 0)
  {
    char* charBuffer = (char*)(&intBuffer[1]);
    std::copy(a_inputT.begin(), a_inputT.end(), charBuffer);
  }
}
template <>
int linearSize<std::string>(const std::string& a_input)
{
  // A string is stored as its length + its data.
  return sizeof(int) + a_input.length() * sizeof(char);
}

//Vector<int>  specialization
template < > int linearSize(const Vector<int>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<int>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<int>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

//Vector<unsigned long long>  specialization
template < > int linearSize(const Vector<unsigned long long>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<unsigned long long>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<unsigned long long>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

//Vector<long>  specialization
template < > int linearSize(const Vector<long>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<long>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<long>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

//Vector<Real>  specialization
template < > int linearSize(const Vector<float>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<float>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<float>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

template < > int linearSize(const Vector<double>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<double>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<double>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

//Vector<std::string>  specialization
template < > int linearSize(const Vector<std::string>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<std::string>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<std::string>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

//Vector<Vector<int> >  specialization
template < > int linearSize(const Vector<Vector<int> >& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<Vector<int> >& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<Vector<int> >& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}
// template < >
// void gather<IntVectSet>(Vector<IntVectSet>& a_outVec,
//             const IntVectSet& a_input,
//             int a_dest)
// {
//   Vector<Box> boxes = a_input.boxes();
//   Vector<Vector<Box> > all_boxes;

//   gather (all_boxes, boxes, a_dest);

//   const int num_vecs = all_boxes.size();
//   a_outVec.resize(num_vecs);

//   for (int i = 0; i < num_vecs; ++i)
//   {
//     IntVectSet& ivs = a_outVec[i];
//     const Vector<Box>& vb = all_boxes[i];
//     for (int ibox = 0; ibox < vb.size(); ++ibox)
//     {
//       ivs |= vb[ibox];
//     }
//   }
// }

// return id of unique processor for special serial tasks
int
uniqueProc(const SerialTask::task& a_task)
{
#ifdef NDEBUG
    switch (a_task)
    {
    case SerialTask::compute:
    default:
        return (0);
        //break;  // unreachable break can generate compiler warning
    }
#else
// in mpi, the debugger attaches to process 0
    return (0);
#endif
}

#include "BaseNamespaceFooter.H"
