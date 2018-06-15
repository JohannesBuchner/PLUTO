/***********************************************************
 * This file is part of DEMSort, a high-performance        *
 * distributed-memory external-memory sorting program      *
 *                                                         *
 * Copyright (C) 2010 Mirko Rahn and Johannes Singler      *
 * Karlsruhe Institute of Technology (KIT), Germany        *
 *                                                         *
 * All rights reserved.  Free for academic use.            *
 ***********************************************************/

#ifndef _SORT_H_
#define _SORT_H_

#ifdef CH_MPI
#include <mpi.h>
#endif
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <tr1/memory>
#ifdef _GLIBCXX_PARALLEL
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

#include <Vector.H>
#include "CH_Timer.H"
#include "CH_assert.H"

#include "sort_utils.H"
#include "parstream.H"
#include <queue>

#include "UsingNamespace.H"

#ifdef CH_MPI

#ifdef VTRACE
#include "vt_user.h"
#endif

/*********************************************************************** */

#define RECV MPI_Irecv(buf, buf_size, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG \
                       , env->m_comm, req \
                       )

#define WAIT MPI_Wait(req, stat); PE_NUM sender = stat->MPI_SOURCE

#define ISEND(buf, count, id, tag) \
    { MPI_Request dummy; \
      MPI_Isend(buf, count, MPI_BYTE, id, tag, env->m_comm, &dummy); \
      MPI_Request_free(&dummy); \
    }

#define SAY(id, tag) ISEND(NULL, 0, id, tag)

#define ON_TAG_NEED                          \
    {                                        \
      index_type idx;                        \
      memcpy(&idx, buf, sizeof(index_type)); \
      outbuf[sender] = a_data[idx];          \
    }                                        \
    ISEND(outbuf + sender, sizeof(value_type), sender, TAG_GIVE)

#define UPDATE(i)                                                   \
    if (split_position[i] < 0)                                      \
    {                                                               \
        split_value[i] = less.min_value(); pq.push(i);              \
    }                                                               \
    else if (split_position[i] >= a_num_elements_on_pes[i])         \
    {                                                               \
        split_value[i] = less.max_value(); pq.push(i);              \
    }                                                               \
    else                                                            \
    {                                                               \
        ISEND(split_position + i, sizeof(index_type), i, TAG_NEED); \
        ++open_query;                                               \
    }

/* ************************************************************************* */

#define FOR_ALL_PE(i) for (PE_NUM i = 0; i < env->m_numpes; ++i)

template<class value_type, class Comparator>
index_type* mws_v0( value_type *a_data,
                    int* a_num_elements_on_pes,
                    const index_type a_my_cut,
                    Comparator less,
                    SortComm* env
                    )
{
#ifdef VTRACE
  VT_USER_START("MW-Sel-setup");
#endif
  CH_assert(env->m_numpes > 1);  CH_assert(a_my_cut >= 0);

  index_type total_length = 0;
  FOR_ALL_PE(i)
  {
    total_length += a_num_elements_on_pes[i];
    CH_assert(a_num_elements_on_pes[i] >= 0);
  }
  CH_assert(a_my_cut <= total_length);

  /* *********************************************************************** */

  PE_NUM parent = (env->m_mype - 1) / 2;
  PE_NUM lChild = 2 * env->m_mype + 1;
  PE_NUM rChild = 2 * env->m_mype + 2;
  PE_NUM nChild = 0;

  if (lChild < env->m_numpes)
  {
    ++nChild;
  }
  if (rChild < env->m_numpes)
  {
    ++nChild;
  }

  index_type* split_position = new index_type[env->m_numpes]; // the splitter
  index_type* need = new index_type[env->m_numpes];
  value_type* split_value = new value_type[env->m_numpes];  // split_value[i]==e@i[split_position[i]]
  value_type* outbuf = new value_type[env->m_numpes];

  std::fill(split_position, split_position + env->m_numpes, 0);

  index_type have = 0;
  index_type stepsize = 2;

  while (stepsize < a_my_cut)
    {
      stepsize <<= 1;
    }

  stepsize >>= 1;

  if (a_my_cut >= total_length)
    {
      stepsize = 0;
    }

  /* ********************************************************************** */

  MPI_Datatype mpi_value_type;            //value_type specification for MPI
  MPI_Type_contiguous(sizeof(value_type), MPI_BYTE, &mpi_value_type);
  MPI_Type_commit(&mpi_value_type);

  // initial distribution of a_data[0]

  value_type v0 = (a_num_elements_on_pes[env->m_mype] > 0) ? a_data[0] : less.max_value();

  MPI_Allgather(&v0, 1, mpi_value_type,
                split_value, 1, mpi_value_type,
                env->m_comm
                );

  MPI_Type_free(const_cast<MPI_Datatype*>(&mpi_value_type));

  /* ********************************************************************** */

  typedef IndirectInvertedLess<value_type, Comparator> iiLess;
  typedef std::priority_queue<PE_NUM, std::vector<PE_NUM>, iiLess> pqueue;

  pqueue pq(iiLess(split_value, less));

  FOR_ALL_PE(i)
  {
    pq.push(i);
  }

  unsigned buf_size = std::max(sizeof(value_type), sizeof(index_type));
  char* buf = new char[buf_size];
  MPI_Request* req = new MPI_Request;
  MPI_Status* stat = new MPI_Status;

  /* ********************************************************************** */

  RECV;

  unsigned open_query = 0;
  PE_NUM wChild = nChild;
#ifdef VTRACE
  VT_USER_END("MW-Sel-setup");
  VT_USER_START("MW-Sel-loop-1");
#endif
  while (stepsize > 0)
    {
      while (have < a_my_cut)
        {
          while (open_query > 0)
            {
              WAIT;

              switch (stat->MPI_TAG)
                {
                case TAG_GIVE:
                  --open_query;
                  memcpy(split_value + sender, buf, sizeof(value_type));
                  pq.push(sender);
                  break;
                case TAG_NEED: ON_TAG_NEED;
                  break;
                case TAG_DONE: --wChild;
                  break;
                default: CH_assert(false);
                  break;
                }

              RECV;
            } // open_query > 0

            // go to right in one sequence
            PE_NUM min_index = pq.top();
            pq.pop();

            split_position[min_index] += stepsize;
            have += stepsize;
            if (have < a_my_cut)
              {
                UPDATE(min_index);
              }
        } // have < want

        stepsize /= 2;

        if (stepsize)
          {
            while (!pq.empty()) // clear the queue
              {
                pq.pop();
              }

            FOR_ALL_PE(i)
              {
                if (split_position[i] > 0)
                  {
                    split_position[i] -= stepsize;
                    have -= stepsize;
                    UPDATE(i);
                  }
                else
                  {
                    pq.push(i);
                  }
              }
          }
    }
#ifdef VTRACE
    VT_USER_END("MW-Sel-loop-1");
#endif
    /* *********************************************************************** */

    if (env->m_mype > 0 && wChild == 0)
    {
        SAY(parent, TAG_DONE);
    }

    bool stop = false;
#ifdef VTRACE
    VT_USER_START("MW-Sel-loop-2");
#endif
    while (!stop)
      {
        WAIT;

        switch (stat->MPI_TAG)
          {
          case TAG_NEED: ON_TAG_NEED;
            break;
          case TAG_DONE:
            if (--wChild == 0)
              {
                if (env->m_mype > 0)
                  {
                    SAY(parent, TAG_DONE);
                  }
                else
                  {
                    if (lChild < env->m_numpes)
                    {
                      SAY(lChild, TAG_STOP);
                    }
                    if (rChild < env->m_numpes)
                    {
                      SAY(rChild, TAG_STOP);
                    }
                    stop = true;
                  }
              }
            break;
          case TAG_STOP:
            if (lChild < env->m_numpes)
            {
              SAY(lChild, TAG_STOP);
            }
            if (rChild < env->m_numpes)
            {
              SAY(rChild, TAG_STOP);
            }
            stop = true;
            break;
          case TAG_GIVE:
            --open_query; // ????
            break;
          default: CH_assert(false);
            break;
          }

        RECV;
      }
#ifdef VTRACE
    VT_USER_END("MW-Sel-loop-2");
#endif
    MPI_Cancel(req);
    WAIT;

    /* *********************************************************************** */

    delete[] buf;
    delete req;
    delete stat;
    delete[] need;
    delete[] split_value;
    delete[] outbuf;

    /* *********************************************************************** */

    return split_position;
}

////////////////////////////////////////////////////////////////////////
//  mws_v1 -- v0 _ opt. for exisiting order
////////////////////////////////////////////////////////////////////////
template<class value_type, class Comparator>
index_type* mws_v1( value_type *a_data,
                    int* a_num_elements_on_pes,
                    const index_type a_my_cut,
                    Comparator less,
                    SortComm* env
                    )
{
#ifdef VTRACE
  VT_USER_START("MW-Sel-setup");
#endif
  CH_assert(env->m_numpes > 1);  CH_assert(a_my_cut >= 0);

  index_type total_length = 0;
  FOR_ALL_PE(i)
    {
      total_length += a_num_elements_on_pes[i];
      CH_assert(a_num_elements_on_pes[i] >= 0);
    }
  CH_assert(a_my_cut <= total_length);

  /* *********************************************************************** */

  PE_NUM parent = (env->m_mype - 1) / 2;
  PE_NUM lChild = 2 * env->m_mype + 1;
  PE_NUM rChild = 2 * env->m_mype + 2;
  PE_NUM nChild = 0;

  if (lChild < env->m_numpes)
  {
    ++nChild;
  }
  if (rChild < env->m_numpes)
  {
    ++nChild;
  }

  index_type* split_position = new index_type[env->m_numpes]; // the splitter
  index_type* need = new index_type[env->m_numpes];
  value_type* outbuf = new value_type[env->m_numpes];

  std::fill(split_position, split_position + env->m_numpes, 0);

  index_type stepsize = 2;

  while (stepsize < a_my_cut)
    {
      stepsize <<= 1;
    }

  stepsize >>= 1;

  if (a_my_cut >= total_length)
    {
      stepsize = 0;
    }

  /* ********************************************************************** */

  // initial distribution of a_data[0] and get a_data[n(p)]
  value_type* split_value = new value_type[env->m_numpes];  // split_value[i]==e@i[split_position[i]]
  value_type* elem_top = new value_type[env->m_numpes];
  {
    CH_TIME("MW-Select:AllGather");
    MPI_Datatype mpi_value_type;            //value_type specification for MPI
    MPI_Type_contiguous(sizeof(value_type), MPI_BYTE, &mpi_value_type);
    MPI_Type_commit(&mpi_value_type);

    value_type* elem_tmp = new value_type[2*env->m_numpes];

    value_type v0[2] =
    {
      (a_num_elements_on_pes[env->m_mype] > 0) ? a_data[0] : less.max_value(),
                         a_data[a_num_elements_on_pes[env->m_mype]-1]
    };

    MPI_Allgather( &v0, 2, mpi_value_type,
                   elem_tmp, 2, mpi_value_type,
                   env->m_comm
                   );
    FOR_ALL_PE(i)
    {
      split_value[i] = elem_tmp[2*i];
      elem_top[i] = elem_tmp[2*i + 1];
    }

    delete elem_tmp;

    MPI_Type_free(const_cast<MPI_Datatype*>(&mpi_value_type));
  }

  /* ********************************************************************** */

  typedef IndirectInvertedLess<value_type, Comparator> iiLess;
  typedef std::priority_queue<PE_NUM, std::vector<PE_NUM>, iiLess> pqueue;

  pqueue pq(iiLess(split_value, less));

  FOR_ALL_PE(i)
  {
    pq.push(i);
  }

  // find low procs that can be deleted from process
  index_type have = 0;
  std::vector<bool> mask(env->m_numpes);
  {
    CH_TIME("MW-Select:make-mask");
    FOR_ALL_PE(i)
    {
      mask[i] = false;
    }
    PE_NUM root = pq.top();
    while (true)
      {
        index_type extra = 0; CH_assert(!mask[root]);
        FOR_ALL_PE(i)
          {
            // !mask & bot[i] <= top[root] then have a interferer
            if ( !mask[i] && split_value[i].get_key() <= elem_top[root].get_key() )
              {
                extra += a_num_elements_on_pes[i];
              }
          }
        // see if we can add this
        if ( have + extra <= a_my_cut )
          {
            have += a_num_elements_on_pes[root];
            mask[root] = true;
            split_position[root] = a_num_elements_on_pes[root];
            pq.pop();
            root = pq.top();
            if ( have == a_my_cut ) break;
          }
        else
          {
            break;
          }
      }
  }
  delete[] elem_top; elem_top = 0;
  CH_assert(have <= a_my_cut);

  unsigned buf_size = std::max(sizeof(value_type), sizeof(index_type));
  char* buf = new char[buf_size];
  MPI_Request* req = new MPI_Request;
  MPI_Status* stat = new MPI_Status;

  RECV;

  unsigned open_query = 0;
  PE_NUM wChild = nChild;
#ifdef VTRACE
  VT_USER_END("MW-Sel-setup");
  VT_USER_START("MW-Sel-loop-1");
#endif
  while (stepsize > 0)
    {
      while (have < a_my_cut)
        {
          while (open_query > 0)
            {
              WAIT;

                switch (stat->MPI_TAG)
                  {
                  case TAG_GIVE:
                    --open_query;
                    memcpy(split_value + sender, buf, sizeof(value_type));
                    pq.push(sender);
                    break;
                  case TAG_NEED: ON_TAG_NEED;
                    break;
                  case TAG_DONE: --wChild;
                    break;
                  default: CH_assert(false);
                    break;
                  }

                RECV;
            } // open_query > 0

          // go to right in one sequence
          PE_NUM min_index = pq.top();
          pq.pop();

          split_position[min_index] += stepsize;
          have += stepsize;

          if (have < a_my_cut)
            {
              UPDATE(min_index);
            }
        } // have < want

      stepsize /= 2;

      if (stepsize)
        {
          while (!pq.empty()) // clear the queue
            {
              pq.pop();
            }

          FOR_ALL_PE(i)
            {
              if ( !mask[i] )
                {
                  if (split_position[i] > 0)
                    {
                      split_position[i] -= stepsize;
                      have -= stepsize;

                      UPDATE(i);
                    }
                  else
                    {
                      pq.push(i);
                    }
                }
            }
        }
    }
#ifdef VTRACE
  VT_USER_END("MW-Sel-loop-1");
#endif
  /* *********************************************************************** */

  if (env->m_mype > 0 && wChild == 0)
    {
      SAY(parent, TAG_DONE);
    }

  bool stop = false;
#ifdef VTRACE
  VT_USER_START("MW-Sel-loop-2");
#endif
  while (!stop)
    {
      WAIT;

      switch (stat->MPI_TAG)
        {
        case TAG_NEED: ON_TAG_NEED;
          break;
        case TAG_DONE:
          if (--wChild == 0)
            {
              if (env->m_mype > 0)
                {
                  SAY(parent, TAG_DONE);
                }
              else
                {
                  if (lChild < env->m_numpes)
                  {
                    SAY(lChild, TAG_STOP);
                  }
                  if (rChild < env->m_numpes)
                  {
                    SAY(rChild, TAG_STOP);
                  }
                  stop = true;
                }
            }
          break;
        case TAG_STOP:
          if (lChild < env->m_numpes)
          {
            SAY(lChild, TAG_STOP);
          }
          if (rChild < env->m_numpes)
          {
            SAY(rChild, TAG_STOP);
          }
          stop = true;
          break;
        case TAG_GIVE:
          --open_query; // ????
          break;
        default: CH_assert(false);
          break;
        }

      RECV;
    }
  if (have != a_my_cut)
    {
      std::cout << "\t\t\tERROR [" << env->m_mype << "] have=" << have << " a_my_cut=" << a_my_cut << std::endl;
    }
  CH_assert(have == a_my_cut);

#ifdef VTRACE
  VT_USER_END("MW-Sel-loop-2");
#endif
  MPI_Cancel(req);
  WAIT;

  /* *********************************************************************** */

  delete[] buf;
  delete req;
  delete stat;
  delete[] need;
  delete[] split_value;
  delete[] outbuf;

  /* *********************************************************************** */

  return split_position;
}

/* ************************************************************************* */
////////////////////////////////////////////////////////////////////////
//  mws_v2
////////////////////////////////////////////////////////////////////////
/** Alternative implementation for the interface, using collective operations only.
    With existing order optimization */
template<class value_type, class Comparator>
index_type* mws_v2(value_type *a_data,
                   int* a_num_elements_on_pes,
                   const index_type a_my_cut,
                   Comparator less,
                   SortComm* env
                   )
{
  if ( env->m_numpes == 1)
    {
      index_type* split_position = new index_type[env->m_numpes]; // the splitter
      split_position[0] = 0;
      return split_position;
    }
  CH_assert(a_my_cut >= 0);

  const index_type P = env->m_numpes;

  index_type total_length = 0;
  for (PE_NUM i = 0; i < P; ++i)
    {
      total_length += a_num_elements_on_pes[i];
      CH_assert(a_num_elements_on_pes[i] >= 0);
    }
  CH_assert(a_my_cut <= total_length);

  /* *********************************************************************** */

  //for termination protocol
  PE_NUM parent = (env->m_mype - 1) / 2;
  PE_NUM left_child = 2 * env->m_mype + 1;
  PE_NUM right_child = 2 * env->m_mype + 2;
  PE_NUM num_children = 0;
  if (left_child < P)
    {
      ++num_children;
    }
  if (right_child < P)
    {
      ++num_children;
    }

  index_type* split_positions = new index_type[P];
  index_type* queries = new index_type[P];
  value_type* send_buffer = new value_type[P];

  std::fill(split_positions, split_positions + P, 0); //start at the very left

  index_type stepsize = 2;

  while (stepsize < a_my_cut)
    {
      stepsize <<= 1;
    }
  stepsize >>= 1;
  //stepsize can differ across the PEs, because they select a different rank

  if (a_my_cut == total_length)
    {
      for (PE_NUM p = 0; p < P; ++p)
        split_positions[p] = a_num_elements_on_pes[p];
      stepsize = 0;   //do not select myself, but participate in collective data exchange
    }

  /* ********************************************************************** */
  // initial distribution of first elements
  MPI_Datatype mpi_value_type;            //value_type specification for MPI
  MPI_Type_contiguous(sizeof(value_type), MPI_BYTE, &mpi_value_type);
  MPI_Type_commit(&mpi_value_type);
  value_type* split_values = new value_type[P];  // split_values[i]==e@i[split[i]]
  value_type* elem_top = new value_type[P];
  {
    CH_TIME("MW-Select:AllGather");

    value_type* elem_tmp = new value_type[2*P];  // split_value[i]==e@i[split[i]]

    value_type v0[2] =
    {
      (a_num_elements_on_pes[env->m_mype] > 0) ? a_data[0] : less.max_value(),
                         a_data[a_num_elements_on_pes[env->m_mype]-1]
    };

    MPI_Allgather(&v0, 2, mpi_value_type,
                  elem_tmp, 2, mpi_value_type,
                  env->m_comm
                  );
    FOR_ALL_PE(i)
    {
      split_values[i] = elem_tmp[2*i];
      elem_top[i] = elem_tmp[2*i + 1];
    }

    delete elem_tmp;
  }

  /* ********************************************************************** */
  //initial state
  typedef IndirectInvertedLess<value_type, Comparator> iiless;
  typedef std::priority_queue<PE_NUM, std::vector<PE_NUM>, iiless> pqueue;

  //contains numbers p of the PEs, key is corresponding split_values[p]
  pqueue pq(iiless(split_values, less));
  for (PE_NUM p = 0; p < P; ++p)
    {
      pq.push(p);
    }
  // find low procs that can be deleted from process
  index_type have = 0;
  std::vector<bool> mask(P);
  {
    CH_TIME("MW-Select:make-mask");
    FOR_ALL_PE(i)
    {
      mask[i] = false;
    }
    PE_NUM root = pq.top();
    while (true)
      {
        index_type extra = 0; CH_assert(!mask[root]);
        FOR_ALL_PE(i)
        {
          // !mask & bot[i] <= top[root] then have a interferer
          if ( !mask[i] && split_values[i].get_key() <= elem_top[root].get_key() )
            {
              extra += a_num_elements_on_pes[i];
            }
        }
        // see if we can add this
        if ( have + extra <= a_my_cut )
          {
            have += a_num_elements_on_pes[root];
            mask[root] = true;
            split_positions[root] = a_num_elements_on_pes[root];
            pq.pop();
            root = pq.top();
            if ( have == a_my_cut ) break;
          }
        else
          {
            break;
          }
      }
  }
  CH_assert(have <= a_my_cut);
  delete[] elem_top; elem_top = 0;

  //receive either a query or an answer
  unsigned receive_buffer_size = std::max(sizeof(value_type), sizeof(index_type));
  char* receive_buffer = new char[receive_buffer_size];

  MPI_Request req;
  MPI_Status status;

  /* ********************************************************************** */

#define RECV2 MPI_Irecv(receive_buffer, receive_buffer_size, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, env->m_comm, &req)

#define TELL(id, tag) ISEND(NULL, 0, id, tag)

#define ANSWER_QUERY \
    { \
        index_type idx; \
        memcpy(&idx, receive_buffer, sizeof(index_type)); \
        send_buffer[sender] = a_data[idx]; \
    } \
    ISEND(send_buffer + sender, sizeof(value_type), sender, TAG_ANSWER)

#define CHECK_TERMINATION(condition, done_tag, stop_tag, done, done_before) \
  { \
    if (condition) \
    { \
         if (env->m_mype > 0) \
         { \
             TELL(parent, done_tag); \
         } \
         else \
         { \
             if (left_child < P) \
             { \
                 TELL(left_child, stop_tag); \
             } \
             if (right_child < P) \
             { \
                 TELL(right_child, stop_tag); \
             } \
             done = true; \
         } \
    } \
    done_before = condition; \
}

#define CHECK_FINAL_TERMINATION \
CHECK_TERMINATION(num_waiting_children == 0 && stepsize == 0 && !done_before, TAG_DONE, TAG_STOP, done, done_before)

#define CHECK_SHIFT_RIGHT_TERMINATION \
bool dummy; CHECK_TERMINATION(num_children_done_shift_right == num_children, TAG_SHIFT_RIGHT_DONE, TAG_START_SHIFT_LEFT, shift_left, dummy);

    /* ********************************************************************** */
    // find own split positions and serve queries from others
    //in each round
    // -all PEs shift right until they have excess, exchanging queries and answers*
    // -barrier (using tree-shaped termination protocol)
    // -all PEs shift left once, exchanging updated elements in a collective communication*
    // *PEs only participate as long as they have not found the final split position, however,
    //  they still serve answers to queries and participate in the collective communications
    //barrier (using tree-shaped termination protocol)

  PE_NUM num_waiting_children = num_children;
  bool done = false, done_before = false;

  CHECK_FINAL_TERMINATION;

  RECV2;   //post receive

  while (!done)
    {
      PE_NUM num_children_done_shift_right = 0;

      //do query-answer game until the split position is shifted enough to the right
      while (have < a_my_cut)
        {
          bool query_open = false;

          if (have < a_my_cut && stepsize > 0)
            {
              //find sequence with least split value
              PE_NUM min_seq = pq.top();
              pq.pop();

              //move to the right in this sequence
              split_positions[min_seq] += stepsize;
              have += stepsize;

              if (split_positions[min_seq] < 0)
                {
                  split_values[min_seq] = less.min_value();
                  pq.push(min_seq);
                }
              else if (split_positions[min_seq] >= a_num_elements_on_pes[min_seq])
                {
                  split_values[min_seq] = less.max_value();
                  pq.push(min_seq);
                }
              else
                {
                  //send query for element
                  ISEND(split_positions + min_seq, sizeof(index_type),
                        min_seq, TAG_QUERY
                        );
                  query_open = true;
                }
            }

          //wait for answer, answer other queries in the meantime
          while (query_open)
            {
              MPI_Wait(&req, &status);
              PE_NUM sender = status.MPI_SOURCE;

              switch (status.MPI_TAG)
                {
                case TAG_ANSWER:
                  query_open = false;
                  memcpy(split_values + sender, receive_buffer, sizeof(value_type));
                  pq.push(sender);
                  break;
                case TAG_QUERY:
                  ANSWER_QUERY;
                  break;
                case TAG_SHIFT_RIGHT_DONE:
                  ++num_children_done_shift_right;
                  break;
                case TAG_DONE:  //child is done
                  --num_waiting_children;
                  CHECK_FINAL_TERMINATION;
                  break;
                case TAG_STOP:
                  CH_assert(false);
                  if (left_child < P)
                  {
                    TELL(left_child, TAG_STOP);
                  }
                  if (right_child < P)
                  {
                    TELL(right_child, TAG_STOP);
                  }
                  done = true;
                  break;
                default:
                  CH_assert(false);
                  break;
                }

              RECV2;
            } //num_open_queries > 0
        } //have < rank

      bool shift_left = false;

      CHECK_SHIFT_RIGHT_TERMINATION;

      //wait for the command to shift left (kind of a barrier), answer queries meanwhile
      while (!shift_left)
        {
          MPI_Wait(&req, &status);
          PE_NUM sender = status.MPI_SOURCE;

          switch (status.MPI_TAG)
            {
            case TAG_QUERY:
              ANSWER_QUERY;
              break;
            case TAG_SHIFT_RIGHT_DONE:
              ++num_children_done_shift_right;
              CHECK_SHIFT_RIGHT_TERMINATION;
              break;
            case TAG_START_SHIFT_LEFT:
              if (left_child < P)
              {
                TELL(left_child, TAG_START_SHIFT_LEFT);
              }
              if (right_child < P)
              {
                TELL(right_child, TAG_START_SHIFT_LEFT);
              }
              shift_left = true;
              break;
            case TAG_DONE:
              --num_waiting_children;
              CHECK_FINAL_TERMINATION;
              break;
            case TAG_STOP:
              if (left_child < P)
              {
                TELL(left_child, TAG_STOP);
              }
              if (right_child < P)
              {
                TELL(right_child, TAG_STOP);
              }
              done = true;
              break;
            default:
              CH_assert(false);
              break;
            }

          RECV2;
        }

      if (done)
        {
          break;
        }
      //now everybody is ready to shift left

      stepsize /= 2;

      if (stepsize > 0)
        {
          while (!pq.empty()) // clear the queue
            pq.pop();

          //shift left in each sequence
          for (PE_NUM p = 0; p < P; ++p)
            if (split_positions[p] > 0 && !mask[p] )
              {
                split_positions[p] -= stepsize;
                have -= stepsize;
              }
        }

      //query the elements on all split positions in a collective communication operation
      MPI_Alltoall(split_positions, 1, MPI_INDEX_TYPE,
                   queries, 1, MPI_INDEX_TYPE,
                   env->m_comm);

      //prepare answers
      for (PE_NUM p = 0; p < P; ++p)
        {
          if (queries[p] < 0)
            send_buffer[p] = less.min_value();
          else if (queries[p] >= a_num_elements_on_pes[env->m_mype])
            send_buffer[p] = less.max_value();
          else
            send_buffer[p] = a_data[queries[p]];
        }

      //exchange answers
      MPI_Alltoall(send_buffer, 1, mpi_value_type,
                   split_values, 1, mpi_value_type,
                   env->m_comm);

      //refill priority queue
      if (stepsize > 0)
        {
          for (PE_NUM p = 0; p < P; ++p)
            {
              if ( !mask[p] )
                {
                  pq.push(p);
                }
            }
        }

      CHECK_FINAL_TERMINATION;
    }   // !done
  if (have != a_my_cut)
    {
      std::cout << "\t\t\tERROR [" << env->m_mype << "] have=" << have << " a_my_cut=" << a_my_cut << std::endl;
    }
  CH_assert(have == a_my_cut);

  MPI_Cancel( &req );
  MPI_Wait( &req, &status );

  /* *********************************************************************** */

  MPI_Type_free(const_cast<MPI_Datatype*>(&mpi_value_type));

  delete[] receive_buffer;
  delete[] queries;
  delete[] split_values;
  delete[] send_buffer;

  /* *********************************************************************** */

  return split_positions;
}

////////////////////////////////////////////////////////////////////////
//  mws_v4 -- v0 with optimization for inexact load balance and initial guess
////////////////////////////////////////////////////////////////////////
template<class value_type, class Comparator>
index_type* mws_v4(value_type *a_data,
                   int* a_num_elements_on_pes,
                   const index_type a_want,
                   float a_margin,
                   Comparator less,
                   SortComm* env
                   )
{
#ifdef VTRACE
  VT_USER_START("MW-Sel-setup");
#endif
  const int P = env->m_numpes, mype = env->m_mype;
  CH_assert(env->m_numpes > 1);  CH_assert(a_want >= 0);

  index_type total_length = 0;
  FOR_ALL_PE(i)
  {
    total_length += a_num_elements_on_pes[i];
    CH_assert(a_num_elements_on_pes[i] >= 0);
  }
  CH_assert(a_want <= total_length);

  int goal = (int)(total_length/P); CH_assert(goal>0);
  CH_assert(a_margin>0.0 && a_margin < 1.0);
  const int margin = (int)( (float)goal*a_margin ) + 1;
if (mype==P-1)std::cout << "\t[" << mype << "]mws_v4 margin = " << margin<< ", a_margin = " << a_margin << ", num elems = " << total_length << std::endl;

  /* *********************************************************************** */

  PE_NUM parent = (env->m_mype - 1) / 2;
  PE_NUM lChild = 2 * env->m_mype + 1;
  PE_NUM rChild = 2 * env->m_mype + 2;
  PE_NUM nChild = 0;

  if (lChild < env->m_numpes)
  {
    ++nChild;
  }
  if (rChild < env->m_numpes)
  {
    ++nChild;
  }

  index_type* split_position = new index_type[env->m_numpes]; // the splitter
  index_type* need = new index_type[env->m_numpes];
  value_type* split_value = new value_type[env->m_numpes];  // split_value[i]==e@i[split_position[i]]
  value_type* outbuf = new value_type[env->m_numpes];

  std::fill(split_position, split_position + env->m_numpes, 0);

  index_type have = 0;
  index_type stepsize = 2;

  while (stepsize < a_want)
    {
      stepsize <<= 1;
    }

  stepsize >>= 1;

  if (a_want >= total_length)
    {
      stepsize = 0;
    }

  /* ********************************************************************** */

  MPI_Datatype mpi_value_type;            //value_type specification for MPI
  MPI_Type_contiguous(sizeof(value_type), MPI_BYTE, &mpi_value_type);
  MPI_Type_commit(&mpi_value_type);

  // initial distribution of a_data[0]

  value_type v0 = (a_num_elements_on_pes[env->m_mype] > 0) ? a_data[0] : less.max_value();

  MPI_Allgather(&v0, 1, mpi_value_type,
                split_value, 1, mpi_value_type,
                env->m_comm
                );

  MPI_Type_free(const_cast<MPI_Datatype*>(&mpi_value_type));

  /* ********************************************************************** */

  typedef IndirectInvertedLess<value_type, Comparator> iiLess;
  typedef std::priority_queue<PE_NUM, std::vector<PE_NUM>, iiLess> pqueue;

  pqueue pq(iiLess(split_value, less));

  FOR_ALL_PE(i)
  {
    pq.push(i);
  }

  unsigned buf_size = std::max(sizeof(value_type), sizeof(index_type));
  char* buf = new char[buf_size];
  MPI_Request* req = new MPI_Request;
  MPI_Status* stat = new MPI_Status;

  /* ********************************************************************** */

  RECV;

  unsigned open_query = 0;
  PE_NUM wChild = nChild;
#ifdef VTRACE
  VT_USER_END("MW-Sel-setup");
  VT_USER_START("MW-Sel-loop-1");
#endif
  int ii = 0;
  while (stepsize > 0)
    {
      while ( have < a_want - margin )
        {
          while (open_query > 0)
            {

              WAIT;

              switch (stat->MPI_TAG)
                {
                case TAG_GIVE:
                  --open_query;
                  memcpy(split_value + sender, buf, sizeof(value_type));
                  pq.push(sender);
                  break;
                case TAG_NEED: ON_TAG_NEED;
                  break;
                case TAG_DONE: --wChild;
                  break;
                default: CH_assert(false);
                  break;
                }

              RECV;
            } // open_query > 0

            // go to right in one sequence
            PE_NUM min_index = pq.top();
            pq.pop();

            split_position[min_index] += stepsize;
            have += stepsize;

            if (have < a_want - margin)
            {
                UPDATE(min_index);
            }
        } // have < want

        stepsize /= 2;

        if (stepsize)
        {
            while (!pq.empty()) // clear the queue
            {
                pq.pop();
            }

            FOR_ALL_PE(i)
            {
                if (split_position[i] > 0)
                {
                    split_position[i] -= stepsize;
                    have -= stepsize;

                    UPDATE(i);
                }
                else
                {
                    pq.push(i);
                }
            }
        }
    }
#ifdef VTRACE
    VT_USER_END("MW-Sel-loop-1");
#endif
    /* *********************************************************************** */

    if (env->m_mype > 0 && wChild == 0)
    {
        SAY(parent, TAG_DONE);
    }

    bool stop = false;
#ifdef VTRACE
    VT_USER_START("MW-Sel-loop-2");
#endif
    while (!stop)
    {
        WAIT;

        switch (stat->MPI_TAG)
        {
        case TAG_NEED: ON_TAG_NEED;
            break;
        case TAG_DONE:
            if (--wChild == 0)
            {
                if (env->m_mype > 0)
                {
                    SAY(parent, TAG_DONE);
                }
                else
                {
                    if (lChild < env->m_numpes)
                    {
                        SAY(lChild, TAG_STOP);
                    }
                    if (rChild < env->m_numpes)
                    {
                        SAY(rChild, TAG_STOP);
                    }
                    stop = true;
                }
            }
            break;
        case TAG_STOP:
            if (lChild < env->m_numpes)
            {
                SAY(lChild, TAG_STOP);
            }
            if (rChild < env->m_numpes)
            {
                SAY(rChild, TAG_STOP);
            }
            stop = true;
            break;
        default: CH_assert(false);
            break;
        }

        RECV;
    }
#ifdef VTRACE
    VT_USER_END("MW-Sel-loop-2");
#endif
    MPI_Cancel(req);
    WAIT;

    /* *********************************************************************** */

    delete[] buf;
    delete req;
    delete stat;
    delete[] need;
    delete[] split_value;
    delete[] outbuf;

    /* *********************************************************************** */

    return split_position;
}

#endif // CH_MPI

////////////////////////////////////////////////////////////////////////
// class ParticleVector - public Vector<T>
////////////////////////////////////////////////////////////////////////
template <class T, class Comparator>
class ParticleVector : public Vector<T>
{
public:
  ParticleVector(unsigned int size) : Vector<T>(size)
#ifdef CH_MPI
, m_histo(0)
#endif
  {}
#ifdef CH_MPI
  ~ParticleVector()
  {
    if (m_histo!=0) delete m_histo;
  }
  //
  void global_sort( MPI_Comm a_comm, Comparator a_less, const int a_version = 2 );
#else
  void global_sort( Comparator a_less, const int a_version = 2 );
#endif
private:
  T* memsort( T* a_data,
              int *a_local_len,
              const int a_capacity,
              Comparator a_less,
#ifdef CH_MPI
              MPI_Comm a_comm,
#endif
              const int a_version = 2);
#ifdef CH_MPI
  index_type* mws_lborder( T *a_data,
                           int a_nloc,
                           float a_margin,
                           int a_version,
                           Comparator a_less,
                           SortComm* a_env);

  int* mws_sdispl( T *a_data,
                   int a_nloc,
                   float a_margin,
                   Comparator a_less,
                   SortComm* a_env);
  // data
  Histogram *m_histo;
#endif
};

////////////////////////////////////////////////////////////////////////
// ParticleVector<T,Comparator>::global_sort
////////////////////////////////////////////////////////////////////////
template<class T, class Comparator>
void ParticleVector<T,Comparator>::global_sort(
#ifdef CH_MPI
                                               MPI_Comm a_comm,
#endif
                                               Comparator a_less,
                                               const int a_version )
{
  int orig_size = this->size(), new_size =  this->size();
  T *data = &(*this)[0];

  // sort
  T* ret = memsort( data, &new_size, this->capacity(), a_less,
#ifdef CH_MPI
                    a_comm,
#endif
                    a_version );

  if (ret != 0 )
    {
      CH_assert(new_size <= this->capacity());
      // hack to deal with array that changed size -- can we manually move the stack ptr?
      for (int i=new_size ; i < orig_size ; i++) this->pop_back(); // this call the destructor!
      for (int i=orig_size ; i < new_size ; i++ ) this->push_back( data[i] ); // adding itself (direct access to stack pointer?)
      CH_assert(new_size == this->size());
    }
}

#ifdef CH_MPI

////////////////////////////////////////////////////////////////////////
// ParticleVector<T,Comparator>::mws_lborder
////////////////////////////////////////////////////////////////////////
/** This function searches for the smallest elements in sorted container
 * \param data is a object with the operator[]
 * \param a_nloc is the size of the container on this PE
 * \param less is the compare function
 * \param a_env
 * \return array of P indices. Each entry states the number of elements on the PEs belonging to the sum(i=0..rank-1)num_elements_on_pe@i smallest elements
 * This function is a collective function (all PEs in comm have to call this function)
 * The PE do not have to use the same container_type, but the same value_type and Comparator
 */
template<class T, class Comparator>
index_type* ParticleVector<T,Comparator>::mws_lborder( T *a_data,
                                                       int a_nloc, // sort can change this
                                                       float a_margin,
                                                       int a_version,
                                                       Comparator a_less,
                                                       SortComm* a_env
                                                       )
{
  // Johannes Singler's methods
  //the lengths of all containers
  int* num_el_pe = new int[a_env->m_numpes];
  MPI_Allgather( &a_nloc, 1, MPI_INT, num_el_pe, 1, MPI_INT, a_env->m_comm);

  //total number of Elements, which should be selected
  index_type num_elems = 0;
  for (int id = 0; id < a_env->m_numpes; ++id)
    {
      num_elems += num_el_pe[id];
    }
  index_type start_rank = (a_env->m_mype*num_elems)/a_env->m_numpes;

  // MW-Sel
  index_type* lborder;
  switch(a_version)
    {
    case 0:
      lborder = mws_v0<T, Comparator>
        (a_data, num_el_pe, start_rank, a_less, a_env );
      break;
    case 1:
      lborder = mws_v1<T, Comparator>
        (a_data, num_el_pe, start_rank, a_less, a_env );
      break;
    case 2:
      lborder = mws_v2<T, Comparator>
        (a_data, num_el_pe, start_rank, a_less, a_env );
      break;
    case 4:
      lborder = mws_v4<T, Comparator>
        (a_data, num_el_pe, start_rank, a_margin, a_less, a_env );
      break;
    default: CH_assert(0);
    }

  delete[] num_el_pe;

  return lborder;
}
//
// ParticleVector<T,Comparator>::mws_sdispl
//
template<class T, class Comparator>
int* ParticleVector<T,Comparator>::mws_sdispl( T *a_data,
                                               int a_nloc,
                                               float a_margin,
                                               Comparator a_less,
                                               SortComm* a_env
                                               )
{
  key_type max_v = a_less.max_value().get_key();
  key_type mink = a_less.max_value().get_key(), maxk = a_less.min_value().get_key();
  key_type *keys = new key_type[a_nloc];
  for (int i=0;i<a_nloc;i++)
    {
      key_type k = a_data[i].get_key();
      if (k == max_v) k = max_v - 1; // need have a sentinal value ???
      keys[i] = k;
      if ( k<mink ) mink = k;
      if ( k>maxk ) maxk = k;
    }
  key_type mink2, maxk2;
  MPI_Allreduce( &mink, &mink2, 1, MPI_INDEX_TYPE, MPI_MIN, a_env->m_comm );
  MPI_Allreduce( &maxk, &maxk2, 1, MPI_INDEX_TYPE, MPI_MAX, a_env->m_comm );
  if (m_histo==0)
    {
      m_histo = new Histogram( a_env, mink2, maxk2 );
    }
  else
    {
      m_histo->set_limits(mink2, maxk2);
    }

  int* sdispl = m_histo->mws_sdispl( keys, a_nloc, a_env->m_comm, a_margin );

  delete [] keys;

  return sdispl;
}

#endif // CH_MPI

// ************************************************************************* //

inline index_type sort_messages(index_type count, int max_elems_per_message)
{
    index_type quot = count / max_elems_per_message;
    index_type rem = count % max_elems_per_message;

    return (quot + ((rem > 0) ? 1 : 0));
}

// ************************************************************************* //

//This is the Algorithm for sorting in Memory.
////////////////////////////////////////////////////////////////////////
// ParticleVector<T,Comparator>::memsort
////////////////////////////////////////////////////////////////////////
template<class T, class Comparator>
T* ParticleVector<T,Comparator>::memsort( T* a_data,
                                          int *a_local_len,
                                          const int a_capacity,
                                          Comparator a_less,
#ifdef CH_MPI
                                          MPI_Comm a_comm,
#endif
                                          const int a_version)
{
  const int loc_s_len = *a_local_len;
#ifdef VTRACE
  VT_USER_START("Sort-local-sort");
#endif
  CH_TIMERS("memsort");

  CH_TIMER("LocalSort", t1);
  CH_START(t1);
  std::sort(a_data, a_data + loc_s_len, a_less);
  CH_STOP(t1);
#ifdef VTRACE
  VT_USER_END("Sort-local-sort");
#endif

#ifdef CH_MPI

  SortComm tcomm(a_comm), *env = &tcomm; // work with old code

  if (env->m_numpes == 1) return a_data; // all done

  MultiwaySelectResult *msr;
  float b_marg = 0.05; // load balance criteria
  if ( a_version != 3 )
    {
      CH_TIMER("Sort:merge-multiway-select", t9);
      CH_START(t9);
      index_type* lborder = mws_lborder(a_data, loc_s_len, b_marg, a_version,
                                        a_less, env );
      MPI_Barrier( a_comm ); //for correct sort wall time
      CH_STOP(t9);
      //std::cout << "\t[" << env->m_mype << "] lborder=" << lborder[0] << " " << lborder[1] << " " << lborder[2] << " " << lborder[3] << std::endl;
      CH_TIMER("Sort:Select-result", t2);
      CH_START(t2);
#ifdef VTRACE
      VT_USER_START("Sort-WM-result");
#endif
      msr = new MultiwaySelectResult( lborder, loc_s_len, env );
      delete [] lborder;
#ifdef VTRACE
      VT_USER_END("Sort-WM-result");
      VT_USER_START("Sort-Send-Recv");
#endif
      CH_STOP(t2);
    }
  else
    {
      CH_TIMER("Sort:bin-multiway-select", t9);
      CH_START(t9);
      int* sdispl = mws_sdispl(a_data, loc_s_len, b_marg, a_less, env );
      MPI_Barrier( a_comm ); //for correct sort wall time
      CH_STOP(t9);

      CH_TIMER("Sort:Select-result", t2);
      CH_START(t2);
#ifdef VTRACE
      VT_USER_START("Sort-WM-result");
#endif
      msr = new MultiwaySelectResult( sdispl, env );

#ifdef VTRACE
      VT_USER_END("Sort-WM-result");
      VT_USER_START("Sort-Send-Recv");
#endif
      CH_STOP(t2);
    }

  // debug
  for (int ii=0;ii<env->m_numpes;ii++)
  {
    if (msr->m_rcount[ii]<0)
    {
      std::cout << "\t\t\tERROR [" << env->m_mype << "] rcount=" << msr->m_rcount[0] << " " << msr->m_rcount[1] << " " << msr->m_rcount[2] << " " << msr->m_rcount[3] << std::endl;
      std::cout << "\t\t\tERROR [" << env->m_mype << "] sdispl=" << msr->m_sdispl[0] << " " << msr->m_sdispl[1] << " " << msr->m_sdispl[2] << " " << msr->m_sdispl[3] << " " << msr->m_sdispl[4] << std::endl;
      std::cout << "\t\t\tERROR [" << env->m_mype << "] scount=" << msr->m_scount[0] << " " << msr->m_scount[1] << " " << msr->m_scount[2] << " " << msr->m_scount[3] << std::endl;
      std::cout << "\t\t\tERROR [" << env->m_mype << "] rdispl=" << msr->m_rdispl[0] << " " << msr->m_rdispl[1] << " " << msr->m_rdispl[2] << " " << msr->m_rdispl[3] << std::endl;
      exit(122);
    }
  }

  // compute new local length
  int local_recv_len = 0;
  for (int ii=0;ii<env->m_numpes;ii++) local_recv_len += msr->m_rcount[ii];
  *a_local_len = local_recv_len; // output
  if ( local_recv_len > a_capacity || local_recv_len < 0 )
    {
      std::cout << "\t\t[" << env->m_mype << "]memsort ERROR: new local size = " << local_recv_len  << ", capacity = " << a_capacity  << std::endl;
      exit(12);
    }

  MPI_Barrier(a_comm); //for correct sort wall time

  T* merge_data = new T[local_recv_len];

  std::vector<T*> slab( env->m_numpes + 1 );
  slab[0] = merge_data;
  for (PE_NUM id = 1; id <= env->m_numpes; ++id)
    {
      slab[id] = slab[id - 1] + msr->m_rcount[id - 1];
    }

  CH_TIMER("Sort-send-recv.", t3);
  CH_START(t3);

  MPI_Datatype mpi_value_type; //value_type specification for MPI
  MPI_Type_contiguous(sizeof(T), MPI_BYTE, &mpi_value_type);
  MPI_Type_commit( &mpi_value_type );

  if (true)  // collective send/recv good for unsorted data that requires all-to-all
    {
      MPI_Alltoallv( a_data, msr->m_scount, msr->m_sdispl, mpi_value_type,
                     merge_data, msr->m_rcount, msr->m_rdispl, mpi_value_type,
                     env->m_comm);
    }
  else
    {
      int max_elems_per_message = (1 << 20);
      index_type send = 0;
      index_type recv = 0;
      for (PE_NUM id = 0; id < env->m_numpes; ++id)
        {
          send += sort_messages(msr->m_scount[id], max_elems_per_message);
          recv += sort_messages(msr->m_rcount[id], max_elems_per_message);
        }

      // BSEND stuff
      unsigned bsend_buf_size = loc_s_len*sizeof(T) + send*MPI_BSEND_OVERHEAD;
      char* bsend_buf = new char[bsend_buf_size];
      MPI_Buffer_attach( bsend_buf, bsend_buf_size );

      MPI_Request* req = new MPI_Request[send + recv];
      index_type req_pos = 0;
      for (PE_NUM k = 0,
             id_recv = pe_right(env->m_mype, env),
             id_send = pe_left(env->m_mype, env)

             ; k < env->m_numpes

             ; ++k,
             id_recv = pe_right(id_recv, env),
             id_send = pe_left(id_send, env)
           )
        {
          T* buf_send = a_data + msr->m_sdispl[id_send];
          T* buf_recv = slab[id_recv];

          while (msr->m_scount[id_send] > 0 && msr->m_rcount[id_recv] > 0)
            {
              index_type elems_send = std::min(msr->m_scount[id_send],
                                               max_elems_per_message
                                               );

              //sum += elems_send;

              MPI_Ibsend(buf_send, // BSEND stuff
                         elems_send,
                         mpi_value_type,
                         id_send,
                         TAG_ALLTOALLV_MANUAL,
                         env->m_comm,
                         req + req_pos
                         );

              ++req_pos;
              buf_send += elems_send;
              msr->m_scount[id_send] -= elems_send;

              index_type elems_recv = std::min(msr->m_rcount[id_recv],
                                               max_elems_per_message
                                               );
              MPI_Irecv(buf_recv,
                        elems_recv,
                        mpi_value_type,
                        id_recv,
                        TAG_ALLTOALLV_MANUAL,
                        env->m_comm,
                        req + req_pos
                        );

              ++req_pos;
              buf_recv += elems_recv;
              msr->m_rcount[id_recv] -= elems_recv;
            }

          while (msr->m_scount[id_send] > 0)
            {
              index_type elems_send = std::min(msr->m_scount[id_send],
                                               max_elems_per_message
                                               );

              MPI_Ibsend(buf_send, // BSEND stuff
                         elems_send,
                         mpi_value_type,
                         id_send,
                         TAG_ALLTOALLV_MANUAL,
                         env->m_comm,
                         req + req_pos
                         );

              ++req_pos;
              buf_send += elems_send;
              msr->m_scount[id_send] -= elems_send;
            }

          while (msr->m_rcount[id_recv] > 0)
            {

              index_type elems_recv = std::min(msr->m_rcount[id_recv],
                                               max_elems_per_message
                                               );

              MPI_Irecv(buf_recv,
                        elems_recv,
                        mpi_value_type,
                        id_recv,
                        TAG_ALLTOALLV_MANUAL,
                        env->m_comm,
                        req + req_pos
                        );

              ++req_pos;
              buf_recv += elems_recv;
              msr->m_rcount[id_recv] -= elems_recv;
            }
        }

      MPI_Waitall(req_pos, req, MPI_STATUSES_IGNORE);

      int tt = loc_s_len;
      MPI_Buffer_detach( &bsend_buf, &tt ); // BSEND stuff
      delete bsend_buf;

      delete [] req;
    }

  MPI_Type_free(const_cast<MPI_Datatype*>(&mpi_value_type));

  delete msr;

  CH_STOP(t3);

#ifdef VTRACE
  VT_USER_END("Sort-Send-Recv");
  VT_USER_START("Sort-final-merge");
#endif

  //merging
  CH_TIMER("Final-merge", t5);
  CH_START(t5);
  if (true)
    {
      std::vector<std::pair<T*, T*> > sub_run_iterators;
      index_type length = 0;
      for (unsigned s = 0; s + 1 < slab.size(); ++s)
      {
        sub_run_iterators.push_back(std::make_pair(slab[s],
                                                   slab[s + 1]
                                                   )
                                    );
        length += slab[s + 1] - slab[s];
      }
      CH_assert(local_recv_len==length);

      CH_assert(false);  // special gnu4.4 code commented out for now
      /*
      __gnu_parallel::multiway_merge(sub_run_iterators.begin(),
                                     sub_run_iterators.end(),
                                     a_data,
                                     length,
                                     a_less
                                     );
      */
    }
  else // no merge
    {
      for ( int ii = 0 ; ii < local_recv_len ; ii++ )
        {
          a_data[ii] = merge_data[ii];
        }
    }
  delete [] merge_data;

  CH_STOP(t5);

#endif // CH_MPI

#ifdef VTRACE
  VT_USER_END("Sort-final-merge");
#endif

  return a_data;
}

#endif //*SORT_H_*//
