#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SORT_HISTI_H_
#define _SORT_HISTI_H_

#include <iterator>

#include "CH_Timer.H"

////////////////////////////////////////////////////////////////////////
// Histogram::Histogram
////////////////////////////////////////////////////////////////////////
Histogram::Histogram( SortComm* env, key_type min, key_type max )
  : m_min_key(min), m_max_key(max), m_num_pes(env->m_numpes),
    m_my_pe(env->m_mype), m_ncalls(0)
{
  m_probe_size = m_num_pes - 1; // probe size
  m_probe_size_max = 2*m_probe_size; // ?????
  m_last_probe = new key_type[m_probe_size];
  m_splitters = new splitter[m_num_pes-1];
}

void Histogram::set_limits( key_type min, key_type max )
{
  m_min_key = min; m_max_key = max;
}

////////////////////////////////////////////////////////////////////////
// Histogram::Init
////////////////////////////////////////////////////////////////////////
void Histogram::Init( MPI_Comm a_comm,
                      const int a_nloc,
                      float a_margin )
{
  //set up initial probe (evenly spaced guesses in data range)
  index_type nl = a_nloc;
  MPI_Allreduce( &nl, &m_num_total_keys, 1, MPI_INDEX_TYPE, MPI_SUM, a_comm );
  int goal = (int)(m_num_total_keys/m_num_pes); CH_assert(goal>0);
  CH_assert(a_margin>0.0 && a_margin < 1.0);
  m_margin_error = (int)( (float)goal*a_margin ) + 1;

  //set up splitters and some variables
  m_num_steps = 0;
  m_num_achv = 0;

  //initial probe divides the range into even chunks
  if ( m_ncalls++ == 0 )  // init init
    {
      for ( key_type i = 0; i < m_probe_size; i++ )
        {
          m_last_probe[i] = m_min_key + ((i+1)*(m_max_key-m_min_key))/m_num_pes;
        }
    }
  else
    {
      for (int i = 0; i < m_probe_size ; i++ )
        {
          m_last_probe[i] = m_splitters[i].lowerb_key;
        }
    }

  // init splitters
  for (int i = 0; i < m_num_pes-1; i++)
  {
    m_splitters[i].lowerb_count = 0;
    m_splitters[i].upperb_count = -1;
    m_splitters[i].lowerb_key = m_min_key;
    m_splitters[i].upperb_key = m_max_key;
    m_splitters[i].goal = -1; //need to know total number of keys to determine this
    m_splitters[i].achieved = false;
    m_splitters[i].broadcasted = false;
  }
}

////////////////////////////////////////////////////////////////////////
// Histogram::mws_sdispl
//  return the number of elements needed from every processor
////////////////////////////////////////////////////////////////////////
int* Histogram::mws_sdispl( key_type *a_keys,
                                        int a_nlocal,
                                        MPI_Comm a_comm,
                                        float a_margin )
{
  Init( a_comm, a_nlocal, a_margin );

  splitter *splits = mws_private( a_keys,
                                              a_nlocal,
                                              a_comm );

  // return number or elements needed on each processor, set 'a_nlocal'
  int *sdisp = new int[m_num_pes+1];
  for ( int p=0, ki=0 ; ki < a_nlocal ; /*void*/)
    {
      while ( a_keys[ki] < splits[p].lowerb_key ) ki++;
      sdisp[p++] = ki;
    }
  sdisp[m_num_pes] = a_nlocal;

  return sdisp;
}

////////////////////////////////////////////////////////////////////////
// Histogram::mws_private
////////////////////////////////////////////////////////////////////////
splitter *Histogram::mws_private( key_type *a_keys,
                                              int a_nlocal,
                                              MPI_Comm a_comm )
{
  const int root = 0;

#define FACT 10000000000

  while (m_num_achv == m_num_pes)
    {
      m_num_steps++;
      // create local histogram, could use 'int' but be safe
      std::vector<index_type> histogram(m_probe_size+1+1); // +1: use for prefix sum
      {
        std::vector<index_type> loc_histo(m_probe_size+1, 0);
        for (int i=0; i<a_nlocal ; i++)
          {
            key_type key = a_keys[i], tkey1;
            int k = m_probe_size/2, inc = m_probe_size/4;
CH_assert(k<m_probe_size-1 && k >= 0);
            if (inc==0) inc=1;
CH_assert(m_last_probe[m_probe_size-1] > key);
            while ( (tkey1=m_last_probe[k]) > key || m_last_probe[k+1] <= key )
              {
                if ( inc < 4 ) // linear search
                  {
                    while ( m_last_probe[k] > key ) k--;
                    while ( k < m_probe_size && m_last_probe[k+1] <= key ) k++;
if (!(k<=m_probe_size && k>=0))std::cout << "[" << m_my_pe << "] ERROR: k=" << k << " key=" << key << " max key = " << m_max_key << std::endl;
                    CH_assert(k<=m_probe_size && k >= 0);
CH_assert(m_last_probe[k] <= key && m_last_probe[k+1] > key);
                    break;
                  }
                else
                  {
                    if ( tkey1 > key ) k -= inc;
                    else k += inc;
                    inc /= 2;
CH_assert(k<m_probe_size-1 && k >= 0);
                  }
              }
CH_assert( m_last_probe[k] <= key && m_last_probe[k+1] > key);
            loc_histo[k]++;
          }
        MPI_Reduce(&loc_histo[0], &histogram[1], m_probe_size, MPI_INDEX_TYPE, MPI_SUM, root, a_comm);
        histogram[0] = 0; // room for prefix sum
      }

      if (m_my_pe == root)
        {
// std::cout << "[0] " << m_num_steps << ") margin=" << m_margin_error << ", min key=" << m_min_key << ", max key=" << m_max_key
// << ", max HISTO ="
// << *std::max_element(histogram.begin(), histogram.end())
// << std::endl << "\thistogram: ";
//        copy(histogram.begin(), histogram.end(), ostream_iterator<index_type>(cout, " "));
//        std::cout << std::endl << "m_last_probe: " ;
//        copy(m_last_probe.begin(), m_last_probe.end(), ostream_iterator<index_type>(cout, " "));
//        std::cout << std::endl << "prefix sum: ";
          // form prefix sum
          for (int hi=0; hi<m_probe_size+1 ; hi++) histogram[hi+1] += histogram[hi];
          histogram[0] = 0;
          copy(histogram.begin(), histogram.end(), ostream_iterator<index_type>(cout, " "));
          std::cout << std::endl;

          if (m_num_steps==2) exit(13);

          // create new probes
          RefineProbes( histogram );

          // broadcast new probes
//        m_probe_size = new_probes.size();
//        int tt[2] = {m_probe_size, m_num_achv }; CH_assert(sizeof(index_type) == sizeof(key_type));
//        m_last_probe.clear();
//        m_last_probe.reserve( m_probe_size );
//        while ( !new_probes.empty() )
//          {
//            key_type probe = new_probes.front(); new_probes.pop_front();
//            m_last_probe.push_back( probe );
//          }
//        MPI_Bcast ( tt, 2, MPI_INT, root, a_comm );
//        MPI_Bcast ( &m_last_probe[0], m_last_probe_size, MPI_INDEX_TYPE, root, a_comm );

//        cout << "[0]sent " << m_probe_size << " new probes, NUM_ACHV *** " << m_num_achv << " probes: ";
//        copy(m_last_probe.begin(), m_last_probe.end(), ostream_iterator<index_type>(cout, " "));
//        std::cout << std::endl;
        }
      else // receive: done?, size of probe, probe
        {
          int tt[2];
          MPI_Bcast ( tt, 2, MPI_INT, root, a_comm );
          m_probe_size = tt[0]; m_num_achv = tt[1];
          //m_last_probe.resize( m_probe_size );
          MPI_Bcast ( &m_last_probe[0], m_probe_size, MPI_INDEX_TYPE, root, a_comm );
        }
    }

  MPI_Barrier(a_comm);
  if (m_my_pe==m_num_pes-1) std::cout << "[" << m_num_steps << "] done splits" << std::endl;
  exit(11);

  return m_splitters;
}

////////////////////////////////////////////////////////////////////////
// Histogram::RefineProbes
////////////////////////////////////////////////////////////////////////
void Histogram::RefineProbes( std::vector<index_type> &histogram )
{
  //uint64_t* hist_counts = (uint64_t*)msg->getData();
  int lenhist = histogram.size();

  ///if first step, initialize total length and splitter goals
  if (m_num_steps == 1)
  {
    m_num_total_keys = histogram[lenhist-1];
    index_type hop = m_num_total_keys/m_num_pes;
    index_type offset = 0;
    for (int i = 0; i < m_num_pes - 1; i++)
    {
      offset += hop;
      m_splitters[i].goal = offset;
      m_splitters[i].upperb_count = m_num_total_keys;
    }
  }

  int idxExp = 0;
  int idxReal = 0;
  int splitterCount = 0;
  key_type histval;
  bool borderChanged;
  int justachv = 0;

  std::vector< splitter* > active_splitters;
  std::vector<int> toBroadcast;
  index_type sumhist = histogram[0];
  index_type diffval;
  ///logic for determining which splitters have been achieved and bounding the rest by the closest counts and splitter-guesses
  ///this code needs at least some comments or something
  while (idxReal < m_num_pes-1)
  {
    if (m_splitters[idxReal].achieved)
    {
      //if this splitter has been determined
      if (!m_splitters[idxReal].broadcasted && (idxReal == 0 || m_splitters[idxReal-1].achieved))
      {
        //if previous splitter has now been determined, broadcast this range
        toBroadcast.push_back(idxReal);
        m_splitters[idxReal].broadcasted = true;
      }
      //don't need to analyze determined splitter, so move on
      idxReal++;
      continue;
    }
    //figure out the current count and current guess from last probe
    if (idxExp >= lenhist-1)
    {
      if (idxExp >= lenhist) sumhist = m_num_total_keys;
      else sumhist =  histogram[idxExp];
      histval = m_max_key;
    }
    else
    {
      sumhist = histogram[idxExp];
      histval = m_last_probe[idxExp];
    }

    //find how far the current splitter is from the guess
    if (sumhist > m_splitters[idxReal].goal)
      diffval = sumhist - m_splitters[idxReal].goal;
    else
      diffval = m_splitters[idxReal].goal - sumhist;

    if (diffval < m_margin_error)
    {
      //if we achieved the splitter within the threshold
      //cement the key and the count
      m_splitters[idxReal].upperb_key = histval;
      m_splitters[idxReal].lowerb_key = histval;
      m_splitters[idxReal].achieved = true;
      justachv++;
      m_splitters[idxReal].upperb_count = sumhist;
      m_splitters[idxReal].lowerb_count = sumhist;

      //broadcast the range if the previous one has been achieved also
      if (idxReal == 0 || m_splitters[idxReal-1].achieved)
      {
        toBroadcast.push_back(idxReal);
        m_splitters[idxReal].broadcasted = true;
      }
      if (idxReal == m_num_pes-2)
      {
        //if this is the last splitter broadcast the last range as well
        toBroadcast.push_back(m_num_pes-1);
      }
      idxExp++;
      idxReal++;
    }
    else if (m_splitters[idxReal].goal < sumhist)
    { //if we did not hit the splitter within the threshold
      borderChanged = false;
      //check if this guess produced better bounds and update the bounds

      if (idxExp > 0)
      {
        if (histogram[idxExp-1] > m_splitters[idxReal].lowerb_count || (histogram[idxExp-1] == m_splitters[idxReal].lowerb_count && m_last_probe[idxExp-1] > m_splitters[idxReal].lowerb_key))
        {
          borderChanged = true;
          m_splitters[idxReal].lowerb_count = histogram[idxExp-1];
          m_splitters[idxReal].lowerb_key = m_last_probe[idxExp-1];
        }
      }
      if ((sumhist < m_splitters[idxReal].upperb_count) || (sumhist == m_splitters[idxReal].upperb_count && histval < m_splitters[idxReal].upperb_key))
      {
        borderChanged = true;
        m_splitters[idxReal].upperb_count = sumhist;
        m_splitters[idxReal].upperb_key = histval;
      }
      //if the bounds did not tighten, ideally, this should never happen
      if (!borderChanged)
      {
#ifdef VERBOSE
        if (lenhist >= m_num_pes) CkPrintf("borderChanged should never be false %d %d %llu %llu %llu %llu %llu!!!\n",
                                           idxExp, idxReal, m_splitters[idxReal].goal, sumhist-histogram[idxExp-1], sumhist,
                                           m_splitters[idxReal].lowerb_count,  m_splitters[idxReal].upperb_count);
        CkPrintf("The distributions is probably too dense to handle the given margin of error\n");
#endif
        //if there was some guess and it did not tighten bounds
        if ((m_probe_size - m_num_achv) < m_num_pes)
        {
          //just declare the splitter converged, currently we cannot do any better
          m_splitters[idxReal].achieved = true;
          if (m_splitters[idxReal].goal - m_splitters[idxReal].lowerb_count <= m_splitters[idxReal].upperb_count - m_splitters[idxReal].goal)
          {
            m_splitters[idxReal].upperb_key = m_splitters[idxReal].lowerb_key;
            m_splitters[idxReal].upperb_count = m_splitters[idxReal].lowerb_count;
          }
          else
          {
            m_splitters[idxReal].lowerb_key = m_splitters[idxReal].upperb_key;
            m_splitters[idxReal].lowerb_count = m_splitters[idxReal].upperb_count;
          }
          //if (m_splitters[idxReal].lowerb_count != 0 && m_splitters[idxReal-1].achieved){
          if (idxReal == 0 || m_splitters[idxReal-1].achieved)
          {
            toBroadcast.push_back(idxReal);
            m_splitters[idxReal].broadcasted = true;
          }
          if (idxReal == m_num_pes-2)
          {
            toBroadcast.push_back(m_num_pes-1);
          }
          idxReal++;
          continue;
        }
      }
      //continue to next splitter
      active_splitters.push_back(&m_splitters[idxReal]);
      idxReal++;
    }
    else
    {
      sumhist = histogram[idxExp];
      idxExp++;
    }
  }
  m_num_achv+=justachv;
  splitterCount = active_splitters.size();

  //if we have successfully not determined close enough values for all the splitters, define new probe
  probe *new_probe;
  int num_guesses = 0;
  key_type *guesses;
  if (splitterCount > 0)
  {
    guesses = get_splitter_guesses(&num_guesses, active_splitters);
    //new_probe = new (num_guesses, toBroadcast.size()) probe;
    new_probe = new probe[num_guesses];
    memcpy(new_probe->splitter_guesses, guesses, num_guesses*sizeof(key_type));
  }
  else
    new_probe = new probe[num_guesses]; // (num_guesses, toBroadcast.size())

  new_probe->num_keys = num_guesses;
  new_probe->num_ranges = toBroadcast.size();

  ///send the resolved ranges
  for (int i = 0; i < toBroadcast.size(); i++)
  {
    if (toBroadcast[i] == 0)
    {
      new_probe->send_info[i].lowerb = m_min_key;
      new_probe->send_info[i].upperb = m_splitters[toBroadcast[i]].upperb_key;
      new_probe->send_info[i].total = m_splitters[0].upperb_count;
    }
    else if (toBroadcast[i] == m_num_pes-1)
    {
      new_probe->send_info[i].lowerb = m_splitters[toBroadcast[i]-1].upperb_key;
      new_probe->send_info[i].upperb = m_max_key;
      new_probe->send_info[i].total = m_num_total_keys - m_splitters[toBroadcast[i]-1].upperb_count;
    }
    else
    {
      new_probe->send_info[i].lowerb = m_splitters[toBroadcast[i]-1].upperb_key;
      new_probe->send_info[i].upperb = m_splitters[toBroadcast[i]].upperb_key;
      new_probe->send_info[i].total = m_splitters[toBroadcast[i]].upperb_count - m_splitters[toBroadcast[i]-1].upperb_count;
    }
    new_probe->send_info[i].chare_dest = toBroadcast[i];
  }
#ifdef DEBUG
  for (int i = 0; i < new_probe->num_keys; i++)
    CkPrintf("probe[%d] = %llu\n", i, new_probe->splitter_guesses[i]);
#endif
  //CkSetQueueing(new_probe, CK_QUEUEING_LIFO);
  if (splitterCount == 0)
  {
    //if all splitters have been converged to
    //c2 = CmiWallTimer();
    //buckets.SendAll(new_probe);
    //first_use = false;
  }
  else
  {
    //if we still need to convergeto some splitters
    //buckets.CountProbe(new_probe);
  }
  active_splitters.clear();
  toBroadcast.clear();

  delete [] m_last_probe;
  //delete msg;
  if (splitterCount > 0)
  {
    m_last_probe = guesses;
  }
}

////////////////////////////////////////////////////////////////////////
// Histogram::get_splitter_guesses
////////////////////////////////////////////////////////////////////////
key_type*  Histogram::get_splitter_guesses( int *ptr_num_guesses, std::vector< splitter* > active_splitters )
{
  int num_guesses;// = *ptr_num_guesses;
  int splitterCount = active_splitters.size();

  ///figure out how much oversplitting we can do
  int virtSent = m_probe_size;
  if (m_probe_size > splitterCount)
  {
    int numVirt = m_probe_size / splitterCount;
    int extraVirt = m_probe_size % splitterCount;
    for (int i = 0; i < splitterCount; i++)
    {
      if (i < extraVirt)
        active_splitters[i]->virt = numVirt+1;
      else
        active_splitters[i]->virt = numVirt;
    }
  }
  else
  {
    virtSent = splitterCount+1;
    for (int i = 0; i < splitterCount; i++)
    {
      active_splitters[i]->virt = 1;
    }
  }

  ///define a new probe
  key_type *newsps = new key_type[m_probe_size];
  num_guesses = 0;
  int idx = 0;
  int divisors;
  key_type divjump = 0;
  int accvirt = 0;
  int firstidx = 0;
  ///Create guesses out of unresolved splitters
  while (idx < splitterCount)
  {
    if (accvirt == 0)
      firstidx = idx;
    accvirt += active_splitters[idx]->virt;
    //place multiple splitters guesses into ranges that bound multiple unresolved splitters
    if (idx == splitterCount-1 || active_splitters[idx]->upperb_key != active_splitters[idx+1]->upperb_key)
    {
      divisors = std::min(m_probe_size_max, accvirt);
      if (divisors > 0)
        divjump = (key_type)(active_splitters[idx]->upperb_key - active_splitters[idx]->lowerb_key)/(divisors+1);
      if (divjump == 0)
      {
        divjump = 1;
      }
      ///we first figure out how many guesses to place in the range, now we place them
      for (int i = 0; i < divisors; i++)
      {
        if (divjump == 1 && (i+1)*divjump + active_splitters[idx]->lowerb_key > active_splitters[idx]->upperb_key)
          newsps[num_guesses] = active_splitters[idx]->upperb_key;
        else
          newsps[num_guesses] = (i+1)*divjump + active_splitters[idx]->lowerb_key;
        num_guesses++;
      }
      accvirt = 0;
    }
    idx++;
  }
  *ptr_num_guesses = num_guesses;
  return newsps;
}

#endif
