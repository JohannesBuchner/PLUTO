#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
using std::pow;

#include "BoxLayoutData.H"
#include "FluxBox.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "NamespaceHeader.H"

int LinearizationTest = 0;

Real norm(const BoxLayoutData<FArrayBox>& a_layout,
          const Interval& interval,
          const int& p)
{
  const BoxLayout validRegion = a_layout.boxLayout();
  DataIterator it = validRegion.dataIterator();
  Real norm = 0;
  if (p == 0)  // max norm
    {
      for (;it.ok(); ++it)
        {
          norm = Max(norm, a_layout[it()].norm(validRegion[it()],
                                               0,
                                               interval.begin(),
                                               interval.size()));
        }
#     ifdef CH_MPI
      Real recv;
      int result = MPI_Allreduce(&norm, &recv, 1, MPI_CH_REAL,
                                 MPI_MAX, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
      {
        //bark!!!
        MayDay::Error("sorry, but I had a communcation error on norm");
      }
      norm = recv;
#     endif
    }
  else if (p == 1) // abs sum norm
    {
      for (;it.ok(); ++it)
        {
          norm += a_layout[it()].norm(validRegion[it()],
                                      1,
                                      interval.begin(),
                                      interval.size());
        }
#     ifdef CH_MPI
      Real recv;
      int result = MPI_Allreduce(&norm, &recv, 1, MPI_CH_REAL,
                                 MPI_SUM, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
      {
        //bark!!!
        MayDay::Error("sorry, but I had a communcation error on norm");
      }
      norm = recv;
#     endif
    }
  else
    {
      for (;it.ok(); ++it)
        {
          norm+=a_layout[it()].sumPow(validRegion[it()],
                                      p,
                                      interval.begin(),
                                      interval.size());
        }
#     ifdef CH_MPI
      Real recv;
      int result = MPI_Allreduce(&norm, &recv, 1, MPI_CH_REAL,
                                 MPI_SUM, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
      {
        //bark!!!
        MayDay::Error("sorry, but I had a communcation error on norm");
      }
      norm = recv;
#     endif
      Real invpwr = 1.0/p;
      if (p == 2) norm = sqrt(norm);
      else       norm = pow(norm, invpwr);
    }

  return norm;
}

template <>
BaseFab<int>* DefaultDataFactory<BaseFab<int> >::create(const Box& box,
                                                        int ncomps,
                                                        const DataIndex& a_datInd) const
{
  return new BaseFab<int>(box, ncomps);
}

template <>
FArrayBox* DefaultDataFactory<FArrayBox>::create(const Box& box, int ncomps,
                                          const DataIndex& a_datInd) const
{
  return new FArrayBox(box, ncomps);
}

FABAliasDataFactory::FABAliasDataFactory(const LayoutData<Real*>& aliases)
{
  define(aliases);
}

void FABAliasDataFactory::define(const LayoutData<Real*>& aliases)
{
  aliasPtrs.define(aliases.boxLayout());
  DataIterator dit = aliases.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      aliasPtrs[dit()] = aliases[dit()];
    }
}

FArrayBox* FABAliasDataFactory::create(const Box& box, int ncomps, const DataIndex& a_datInd) const
{
/* Brian thinks this check shouldn't be done.
  const Box& b = aliasPtrs.box(a_datInd);
  if (b != box)
    {
      MayDay::Error("Aliased data holder dimensions do not match across LevelData const.");
    }
*/
  FArrayBox* rtn = new FArrayBox(box, ncomps, aliasPtrs[a_datInd]);
  return rtn;
}

//--FABAliasFlBxDataFactory

FABAliasFlBxDataFactory::FABAliasFlBxDataFactory(
  BoxLayoutData<FluxBox>* a_original,
  const Interval&         a_interval,
  const int               a_dir)
  :
  m_origPointer(a_original),
  m_interval(a_interval),
  m_dir(a_dir)
{
}

void
FABAliasFlBxDataFactory::define(
  BoxLayoutData<FluxBox>* a_original,
  const Interval&         a_interval,
  const int               a_dir)
{
  m_origPointer = a_original;
  m_interval    = a_interval;
  m_dir         = a_dir;
}

FArrayBox*
FABAliasFlBxDataFactory::create(const Box&       a_box,
                                int              a_ncomps,
                                const DataIndex& a_dataInd) const
{
  CH_assert(a_ncomps = m_interval.size());
  FluxBox& origFlBx = m_origPointer->operator[](a_dataInd);
  return new FArrayBox(m_interval, origFlBx[m_dir]);
}

#include "NamespaceFooter.H"
