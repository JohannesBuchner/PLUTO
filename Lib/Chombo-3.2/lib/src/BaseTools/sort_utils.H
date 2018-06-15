/***********************************************************
 * This file is part of DEMSort, a high-performance        *
 * distributed-memory external-memory sorting program      *
 *                                                         *
 * Copyright (C) 2010 Mirko Rahn and Johannes Singler      *
 * Karlsruhe Institute of Technology (KIT), Germany        *
 *                                                         *
 * All rights reserved.  Free for academic use.            *
 ***********************************************************/

#ifndef _SORT_UTILS_H_
#define _SORT_UTILS_H_

// to mark unused parameters explicitly
#ifdef __GNUC__
#define PARAM_UNUSED(x) x __attribute__ ((unused))
#else
#define PARAM_UNUSED(x) x
#endif

typedef unsigned long long key_type;
typedef          long long index_type;

template<class S, class Comparator>
struct IndirectLess
{
  S* data;
  Comparator& less;

  IndirectLess(S* data, Comparator& less) : data(data), less(less)
  {
  }

  inline bool operator () (const index_type& a, const index_type& b)
  {
      return (less(data[a], data[b])
              ? true
              : (less(data[b], data[a])
                 ? false
                 : (a < b) //break ties on position in sequence
                 )
              );
  }
};

template<class S, class Comparator>
struct IndirectInvertedLess
{
  S* data;
  Comparator& less;

  IndirectInvertedLess(S* data, Comparator& less) : data(data), less(less)
  {
  }

  inline bool operator () (const index_type& a, const index_type& b)
  {
      return (less(data[b], data[a])
              ? true
              : (less(data[a], data[b])
                 ? false
                 : (b < a) //break ties on position in sequence
                 )
              );
  }
};

template<class InIt, class Comparator>
struct IndirectIndirectLess
{
  InIt* data;
  Comparator& less;

  IndirectIndirectLess(InIt* data, Comparator& less) : data(data), less(less)
  {
  }

  inline bool operator () (const unsigned& a, const unsigned& b)
  {
      return (less(*(data[a]), *(data[b]))
              ? true
              : (less(*(data[b]), *(data[a]))
                 ? false
                 : (a < b) //break ties on position in sequence
                 )
              );
  }
};

template<class InIt, class Comparator>
struct IndirectIndirectInvertedLess
{
  InIt* data;
  Comparator& less;

  IndirectIndirectInvertedLess(InIt* data, Comparator& less) : data(data), less(less)
  {
  }

  inline bool operator () (const unsigned& a, const unsigned& b)
  {
      return (less(*(data[b]), *(data[a]))
              ? true
              : (less(*(data[a]), *(data[b]))
                 ? false
                 : (b < a) //break ties on position in sequence
                 )
              );
  }
};

#ifdef CH_MPI

#define MPI_INDEX_TYPE        MPI_LONG_LONG_INT
#define MPI_KEY_TYPE          MPI_UNSIGNED_LONG_LONG

// given by MPI_Standard
typedef int PE_NUM;

//template <class T>
class SortComm
{
 public:
 SortComm( MPI_Comm a_comm ) : m_comm(a_comm)
  {
    MPI_Comm_size(a_comm,&m_numpes);
    MPI_Comm_rank(a_comm,&m_mype);
  }
  MPI_Comm m_comm;
  int m_mype;
  int m_numpes;
};

// ** print utils in sort_env.h ******************************************** //

std::string sort_format_pe_number(int i)
{
    std::ostringstream ss;

    ss << std::setw(4) << i;

    return ss.str();
}

// ************************************************************************* //

static std::string sort_format_iam(const SortComm* env)
{
    return sort_format_pe_number(env->m_mype) + ": ";
}

// ************************************************************************* //

std::ostream& sort_pe_out(const SortComm* env)
{
    return std::cout << sort_format_iam(env);
}

std::ostream& sort_pe_err(const SortComm* env)
{
    return std::cerr << sort_format_iam(env) << "ERROR: ";
}

std::ostream& sort_nullstream()  // eats everything, prints nothing
{
    static std::ostream rc(NULL);

    return rc;
}

std::ostream& sort_root_out(const SortComm* env)
{
    return (env->m_mype == 0) ? sort_pe_out(env) : sort_nullstream();
}

std::ostream& sort_root_err(const SortComm* env)
{
    return (env->m_mype == 0) ? sort_pe_err(env) : sort_nullstream();
}

// ************************************************************************* //

void sort_pe_msg(const SortComm* env, const std::string msg)
{
    sort_pe_out(env) << msg << std::endl;
}

void sort_root_msg(const SortComm* env, const std::string msg)
{
    sort_root_out(env) << msg << std::endl;
}

PE_NUM pe_right(PE_NUM id, const SortComm* env)
{
    return (id + 1) % env->m_numpes;
}

PE_NUM pe_left(PE_NUM id, const SortComm* env)
{
    return (id - 1 + env->m_numpes) % env->m_numpes;
}

enum Tag
{
// in makeGlobalRuns, in case of block with filler
    TAG_ALLTOALLV_MANUAL=10001
// in alltoall
    , TAG_EXTERNAL_ALLTOALL=10002

// for msr
    , TAG_NEED=10003,
    TAG_QUERY=10004,
    TAG_GIVE=20005,
    TAG_ANSWER=10006,
    TAG_DONE=10007,
    TAG_STOP=10008,
    TAG_SHIFT_RIGHT_DONE=10009,
    TAG_START_SHIFT_LEFT=10010

    , TAG_CHECK=10011

    , TAG_LEFT_MAX=10012,
    TAG_RIGHT_MIN=10013

    , TAG_GRD=10014

    , TAG_PART_DIST=10015
};

/* ************************************************************************* */
// functions to post-process the result of multiway_select
#include <limits.h>
struct MultiwaySelectResult
{
  int* m_sdispl;
  int* m_scount;
  int* m_rcount;
  int* m_rdispl;
  index_type *m_lborder;

  MultiwaySelectResult(index_type* a_lborder,
                       int a_len,
                       SortComm* a_env
                       )
    : m_lborder(a_lborder)
  {
    m_sdispl = new int[a_env->m_numpes + 1];
    m_scount = new int[a_env->m_numpes];
    m_rcount = new int[a_env->m_numpes];
    m_rdispl = new int[a_env->m_numpes];
    const int P = a_env->m_numpes;
    {
      int *lborder= new int[P];
      for (int ii=0;ii<P; ii++)
        {
          CH_assert(a_lborder[ii]<INT_MAX);
          lborder[ii] = a_lborder[ii];
        }

      MPI_Alltoall( lborder, 1, MPI_INT,
                    m_sdispl, 1, MPI_INT, a_env->m_comm );

      delete [] lborder;
    }
    m_sdispl[P] = a_len;
    //if (a_env->m_mype==0)std::cout << "\t[" << a_env->m_mype <<"] m_sdispl = " << m_sdispl[0] << " " << m_sdispl[1] << " " << m_sdispl[2] << " " << m_sdispl[3] << std::endl;

    for (int ii = 0; ii < P; ++ii)
      {
        m_scount[ii] = m_sdispl[ii + 1] - m_sdispl[ii];
      }
    //std::cout << "\t[" << a_env->m_mype <<"] m_scount = " << m_scount[0] << " " << m_scount[1] << " " << m_scount[2] << " " << m_scount[3] << std::endl;

    MPI_Alltoall(m_scount, 1, MPI_INT,
                 m_rcount, 1, MPI_INT, a_env->m_comm);
    //std::cout << "\t[" << a_env->m_mype <<"] m_rcount = " << m_rcount[0] << " " << m_rcount[1] << " " << m_rcount[2] << " " << m_rcount[3] << std::endl;

    m_rdispl[0] = 0;
    for (int ii = 1; ii < P; ++ii)
      {
        m_rdispl[ii] = m_rdispl[ii - 1] + m_rcount[ii - 1];
      }
    //std::cout << "\t[" << a_env->m_mype <<"] m_rdispl = " << m_rdispl[0] << " " << m_rdispl[1] << " " << m_rdispl[2] << " " << m_rdispl[3] << std::endl;
  }

  // constructor with m_sdispl
  MultiwaySelectResult(int* a_sdisp,
                       SortComm* env
                       )
  {
    m_scount = new int[env->m_numpes];
    m_rcount = new int[env->m_numpes];
    m_rdispl = new int[env->m_numpes];
    m_sdispl = a_sdisp;

    //     MPI_Alltoall(lborder, 1, MPI_INT,
    //           m_sdispl, 1, MPI_INT, env->m_comm);

    //     m_sdispl[env->m_numpes] = a_len;
    std::cout << "\t[" << env->m_mype <<"]MultiwaySelectResult m_sdispl = " << m_sdispl[0] << " " << m_sdispl[1] << " " << m_sdispl[2] << " " << m_sdispl[3] << " " << m_sdispl[4] << std::endl;

    for (int i = 0; i < env->m_numpes; ++i)
      {
        m_scount[i] = m_sdispl[i + 1] - m_sdispl[i];
      }
    //MPI_Barrier(env->m_comm);
    //std::cout << "\t[" << env->m_mype <<"]MultiwaySelectResult m_scount = " << m_scount[0] << " " << m_scount[1] << " " << m_scount[2] << " " << m_scount[3] << std::endl;
    MPI_Alltoall(m_scount, 1, MPI_INT,
                 m_rcount, 1, MPI_INT, env->m_comm);

    std::cout << "\t[" << env->m_mype <<"]MultiwaySelectResult m_rcount = " << m_rcount[0] << " " << m_rcount[1] << " " << m_rcount[2] << " " << m_rcount[3] << std::endl;
    m_rdispl[0] = 0;

    for (int i = 1; i < env->m_numpes; ++i)
      {
        m_rdispl[i] = m_rdispl[i - 1] + m_rcount[i - 1];
      }
    MPI_Barrier(env->m_comm);
    std::cout << "\t[" << env->m_mype <<"]MultiwaySelectResult m_rdispl = " << m_rdispl[0] << " " << m_rdispl[1] << " " << m_rdispl[2] << " " << m_rdispl[3] << std::endl;
  }

  ~MultiwaySelectResult()
  {
    delete[] m_sdispl;
    delete[] m_scount;
    delete[] m_rcount;
    delete[] m_rdispl;
  }
};

/* ************************************************************************* */
// verify

template<typename block_type, typename value_type, class Comparator>
void test_multiway_select(value_type* data,
                          index_type length,
                          MultiwaySelectResult* msr,
                          Comparator less,
                          SortComm* env
                          )
{
    unsigned err = 0;

    index_type* lborder_global = new index_type[env->m_numpes * env->m_numpes];
    index_type* sdispl_global = new index_type[env->m_numpes * env->m_numpes];
    index_type* rdispl_global = new index_type[env->m_numpes * env->m_numpes];
    index_type* scount_global = new index_type[env->m_numpes * env->m_numpes];
    index_type* rcount_global = new index_type[env->m_numpes * env->m_numpes];

    MPI_Allgather(msr->m_lborder, env->m_numpes, MPI_INDEX_TYPE,
                  lborder_global, env->m_numpes, MPI_INDEX_TYPE,
                  env->m_comm
                  );

    CH_assert(0); // this is wrong now -- int not index_type
    MPI_Gather(msr->m_sdispl, env->m_numpes, MPI_INDEX_TYPE,
               sdispl_global, env->m_numpes, MPI_INDEX_TYPE,
               0, env->m_comm
               );

    MPI_Gather(msr->m_rdispl, env->m_numpes, MPI_INDEX_TYPE,
               rdispl_global, env->m_numpes, MPI_INDEX_TYPE,
               0, env->m_comm
               );

    MPI_Gather(msr->m_scount, env->m_numpes, MPI_INDEX_TYPE,
               scount_global, env->m_numpes, MPI_INDEX_TYPE,
               0, env->m_comm
               );

    MPI_Gather(msr->m_rcount, env->m_numpes, MPI_INDEX_TYPE,
               rcount_global, env->m_numpes, MPI_INDEX_TYPE,
               0, env->m_comm
               );

    index_type* rlen = new index_type[env->m_numpes];

    MPI_Gather(&length, 1, MPI_INDEX_TYPE,
               rlen, 1, MPI_INDEX_TYPE,
               0, env->m_comm
               );

    if (env->m_mype == 0)
    {
        /* ******************************************************************* */

        std::cout << "lborder = " << std::endl;

        index_type rsum = 0;

        for (int r = 0; r < env->m_numpes; ++r)
        {
            index_type tsum = 0;

            for (int i = 0; i < env->m_numpes; ++i)
            {
                std::cout << std::setw(11) << lborder_global[r * env->m_numpes + i];

                tsum += lborder_global[r * env->m_numpes + i];
            }

            std::cout << ", sum = " << tsum
                      << " is " << (tsum == rsum ? "okay" : "WRONG")
                      << std::endl;

            err += (tsum == rsum) ? 0 : 1;

            rsum += rlen[r];
        }

        /* ******************************************************************* */

        std::cout << "scount = " << std::endl;

        for (int r = 0; r < env->m_numpes; ++r)
        {
            index_type tsum = 0;

            for (int i = 0; i < env->m_numpes; ++i)
            {
                std::cout << std::setw(11) << scount_global[r * env->m_numpes + i];

                tsum += scount_global[r * env->m_numpes + i];
            }

            std::cout << ", sum = " << tsum
                      << " is " << (tsum == rlen[r] ? "okay" : "WRONG")
                      << std::endl;

            err += (tsum == rlen[r]) ? 0 : 1;
        }

        /* ******************************************************************* */

        std::cout << "rcount = " << std::endl;

        for (int r = 0; r < env->m_numpes; ++r)
        {
            index_type tsum = 0;

            for (int i = 0; i < env->m_numpes; ++i)
            {
                std::cout << std::setw(11) << rcount_global[r * env->m_numpes + i];

                tsum += rcount_global[r * env->m_numpes + i];
            }

            std::cout << ", sum = " << tsum
                      << " is " << (tsum == rlen[r] ? "okay" : "WRONG")
                      << std::endl;

            err += (tsum == rlen[r]) ? 0 : 1;
        }

        /* ******************************************************************* */

        std::cout << "sdispl = " << std::endl;

        for (int r = 0; r < env->m_numpes; ++r)
        {
            err += (sdispl_global[r * env->m_numpes + 0] == 0) ? 0 : 1;

            std::cout << std::setw(11) << sdispl_global[r * env->m_numpes + 0];

            bool asc = true;

            for (int i = 1; i < env->m_numpes; ++i)
            {
                std::cout << std::setw(11) << sdispl_global[r * env->m_numpes + i];

                asc &= (sdispl_global[r * env->m_numpes + i]
                        >= sdispl_global[r * env->m_numpes + i - 1]
                        );
            }

            std::cout << ", ascending: " << (asc ? "okay" : "WRONG") << std::endl;

            err += asc ? 0 : 1;
        }

        /* ******************************************************************* */

        std::cout << "rdispl = " << std::endl;

        for (int r = 0; r < env->m_numpes; ++r)
        {
            err += (rdispl_global[r * env->m_numpes + 0] == 0) ? 0 : 1;

            std::cout << std::setw(11) << rdispl_global[r * env->m_numpes + 0];

            bool asc = true;

            for (int i = 1; i < env->m_numpes; ++i)
            {
                std::cout << std::setw(11) << rdispl_global[r * env->m_numpes + i];

                asc &= (rdispl_global[r * env->m_numpes + i]
                        >= rdispl_global[r * env->m_numpes + i - 1]
                        );
            }

            std::cout << ", ascending: " << (asc ? "okay" : "WRONG") << std::endl;

            err += asc ? 0 : 1;
        }

        /* ******************************************************************* */
    }

    sort_pe_out(env) << "msr local errors: " << err << std::endl;

    unsigned serr;

    MPI_Reduce(&err, &serr, 1, MPI_UNSIGNED, MPI_SUM, 0, env->m_comm);

    err = serr;

    sort_root_out(env) << "msr local error sum: " << err << std::endl;

    /* *********************************************************************** */

    value_type* blval = new value_type[env->m_numpes * env->m_numpes];
    value_type* brval = new value_type[env->m_numpes * env->m_numpes];

    for (int r = 0; r < env->m_numpes; ++r)
    {
        if (lborder_global[r * env->m_numpes + env->m_mype] < length)
        {
            value_type val = data[lborder_global[r * env->m_numpes + env->m_mype]];
            blval[env->m_mype * env->m_numpes + r] = val;
        }
        else
        {
            blval[env->m_mype * env->m_numpes + r] = less.max_value();
        }

        if (lborder_global[r * env->m_numpes + env->m_mype] - 1 < 0)
        {
            brval[env->m_mype * env->m_numpes + r] = less.min_value();
        }
        else
        {
            value_type val = data[lborder_global[r * env->m_numpes + env->m_mype] - 1];
            brval[env->m_mype * env->m_numpes + r] = val;
        }
    }

    MPI_Datatype mpi_value_type;            //value_type specification for MPI
    MPI_Type_contiguous(sizeof(value_type), MPI_BYTE, &mpi_value_type);
    MPI_Type_commit(&mpi_value_type);

    MPI_Gather(&blval[env->m_mype * env->m_numpes], env->m_numpes, mpi_value_type,
               blval, env->m_numpes, mpi_value_type,
               0, env->m_comm
               );
    MPI_Gather(&brval[env->m_mype * env->m_numpes], env->m_numpes, mpi_value_type,
               brval, env->m_numpes, mpi_value_type,
               0, env->m_comm
               );

    MPI_Type_free(const_cast<MPI_Datatype*>(&mpi_value_type));

    if (env->m_mype == 0)
    {
        value_type* vmin = new value_type[env->m_numpes];
        value_type* vmax = new value_type[env->m_numpes];

        std::cout << "blval = " << std::endl;

        for (int r = 0; r < env->m_numpes; ++r)
        {
            vmin[r] = less.max_value();

            for (int k = 0; k < env->m_numpes; ++k)
            {
                std::cout << blval[k * env->m_numpes + r];

                if (less(blval[k * env->m_numpes + r], vmin[r]))
                {
                    vmin[r] = blval[k * env->m_numpes + r];
                }
            }
            std::cout << " => min[" << sort_format_pe_number(r) << "] = " << vmin[r] << std::endl;
        }

        std::cout << "brval = " << std::endl;

        for (int r = 0; r < env->m_numpes; ++r)
        {
            vmax[r] = less.min_value();

            for (int k = 0; k < env->m_numpes; ++k)
            {
                std::cout << brval[k * env->m_numpes + r];

                if (less(vmax[r], brval[k * env->m_numpes + r]))
                {
                    vmax[r] = brval[k * env->m_numpes + r];
                }
            }
            std::cout << " => max[" << sort_format_pe_number(r - 1) << "] = " << vmax[r] << std::endl;
        }

        for (int r = 0; r < env->m_numpes - 1; ++r)
        {
            std::cout << "max[" << sort_format_pe_number(r) << "] <= "
                      << "min[" << sort_format_pe_number(r + 1) << "] iff "
                      << vmax[r + 1] << " <= "
                      << vmin[r + 1] << ": "
                      << (less(vmin[r + 1], vmax[r + 1]) ? "WRONG" : "okay") << std::endl;

            err += less(vmin[r + 1], vmax[r + 1]) ? 1 : 0;
        }

        sort_pe_out(env) << "msr global error count: " << err << std::endl;

        if (err > 0)
        {
          sort_pe_msg(env, "* msr is ERRONEOUS");
          CH_assert(false);
        }
        else
        {
          sort_pe_msg(env, "* msr looks GOOD");
        }

        delete[] vmin;
        delete[] vmax;
    }

    delete[] blval;
    delete[] brval;
    delete[] lborder_global;
    delete[] sdispl_global;
    delete[] rdispl_global;
    delete[] scount_global;
    delete[] rcount_global;
}

// ugly but not used
#include "sort_hist.H"

#endif // CH_MPI

#endif // SORT_UTILS_H
