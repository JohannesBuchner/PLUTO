#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SORT_HIST_H_
#define _SORT_HIST_H_

// this method explicity uses keys instead of a compare method
typedef unsigned long long key_type;

class resolved_range
{
public:
  key_type lowerb;
  key_type upperb;
  int chare_dest;
  int total;
};

class probe
{
public:
  key_type* splitter_guesses;
  int num_keys;
  resolved_range* send_info;
  int num_ranges;
};

struct splitter
{
  index_type upperb_count;
  key_type upperb_key;
  index_type lowerb_count;
  key_type lowerb_key;
  index_type goal;
  int virt;
  bool achieved;
  bool broadcasted;
};

class Histogram
{
public:
  Histogram( SortComm* env, key_type min, key_type max );
  void set_limits( key_type min, key_type max );
  int* mws_sdispl( key_type *, int nloc, MPI_Comm c, float margin );
private:
  key_type m_min_key, m_max_key;
  const int m_num_pes;
  const int m_my_pe;
  int m_num_steps;
  int m_ncalls;
  int m_num_achv;
  int m_probe_size;
  int m_probe_size_max;
  int m_margin_error;
  splitter* m_splitters;
  index_type m_num_total_keys;
  key_type* m_last_probe;
  //std::vector<index_type> m_key_splits;
  //std::vector<key_type> m_probe;
private:
  splitter* mws_private( key_type *, int , MPI_Comm );
  void Init( MPI_Comm a_comm, const int a_nloc, float );
  void RefineProbes( std::vector<index_type> & );
  key_type* get_splitter_guesses(int *ptr_num_guesses, std::vector< splitter* > active_splitters);
};

//need to include .C file in order to have it instantiated when the .h file is included externally
//such instantiation is necessary here since the compiler generates templated code on demand
#include "sort_histI.H"
#endif
