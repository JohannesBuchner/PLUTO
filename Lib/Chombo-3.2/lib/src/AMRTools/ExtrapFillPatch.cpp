#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Box.H"
#include "Interval.H"
#include "SPACE.H"
#include "Tuple.H"
#include "FArrayBox.H"
#include "LayoutData.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "LayoutIterator.H"

#include "ExtrapFillPatch.H"
#include "NamespaceHeader.H"

ExtrapFillPatch::ExtrapFillPatch()
  :
  m_is_defined(false),
  m_extrap_interval(0,-1)
{
}

ExtrapFillPatch::~ExtrapFillPatch()
{
}

ExtrapFillPatch::ExtrapFillPatch(
    const DisjointBoxLayout& a_level_domain,
    const Box& a_problem_domain,
    const Interval& a_extrap_interval
    )
  :
  m_is_defined(false),
  m_extrap_interval(0,-1)
{
  ProblemDomain physDomain(a_problem_domain);
  define(a_level_domain,
         physDomain,
         a_extrap_interval);
}

ExtrapFillPatch::ExtrapFillPatch(
    const DisjointBoxLayout& a_level_domain,
    const ProblemDomain& a_problem_domain,
    const Interval& a_extrap_interval
    )
  :
  m_is_defined(false),
  m_extrap_interval(0,-1)
{
  define(a_level_domain,
         a_problem_domain,
         a_extrap_interval);
}

bool
ExtrapFillPatch::isDefined() const
{
  return ( m_is_defined );
}

void
ExtrapFillPatch::define(
  const DisjointBoxLayout& a_level_domain,
  const Box& a_problem_domain,
  const Interval& a_extrap_interval
  )
{
  ProblemDomain physDomain(a_problem_domain);
  define(a_level_domain, physDomain, a_extrap_interval);
}

void
ExtrapFillPatch::define(
  const DisjointBoxLayout& a_level_domain,
  const ProblemDomain& a_problem_domain,
  const Interval& a_extrap_interval
  )
{
  m_extrap_interval = a_extrap_interval;
  const int interp_radius = m_extrap_interval.begin() - 1;
  const int extrap_width = m_extrap_interval.size();
  const int extrap_radius = m_extrap_interval.end();

  DataIterator dit = a_level_domain.dataIterator();
// check to see if boxes are large enough
#ifndef NDEBUG
  std::vector<bool> distrib_boxes_too_small;
  bool boxes_too_small = false;
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& grid = a_level_domain[dit()];
    if (grid.shortside() < extrap_radius)
    {
      boxes_too_small = true;
    }
  }
  if (boxes_too_small)
  {
    cerr << "ExtrapFillPatch::define: minimum Box dimension must not be less than "
         << extrap_radius << endl;
    MayDay::Abort();
  }
#endif
  IntVectSet valid_region;

  // former, ridiculous version of this algorithm.

 //  LayoutIterator lit = a_level_domain.layoutIterator();
//   for (lit.begin(); lit.ok(); ++lit)
//   {
//     const Box& grid = a_level_domain[lit()];
//     valid_region |= grid;
//   }
//   for (int dir = 0; dir < SpaceDim; ++dir)
//   {
//     m_lo_extrap[dir].define(a_level_domain);
//     m_hi_extrap[dir].define(a_level_domain);
//     for (dit.begin(); dit.ok(); ++dit)
//     {
//       const Box& grid = a_level_domain[dit()];
//       const Box interp_box = grow(grid, interp_radius);
//       const Box lo_extrap_box = adjCellLo(interp_box, dir, extrap_width);
//       const Box hi_extrap_box = adjCellHi(interp_box, dir, extrap_width);
//       IntVectSet& local_lo_extrap = m_lo_extrap[dir][dit()];
//       IntVectSet& local_hi_extrap = m_hi_extrap[dir][dit()];
//       local_lo_extrap |= (lo_extrap_box & a_problem_domain);
//       local_lo_extrap -= valid_region;
//       local_hi_extrap |= (hi_extrap_box & a_problem_domain);
//       local_hi_extrap -= valid_region;
//     }
//   }

  // new version of algorithm that doesn't build GIANT IntVectSets
  LayoutIterator it = a_level_domain.layoutIterator();

  // for periodic case, create a testBox which will only contain
  // cells which are _not_ adjacent to a periodic boundary.
  Box periodicTestBox(a_problem_domain.domainBox());
  if (a_problem_domain.isPeriodic())
    {
      for (int idir=0; idir<SpaceDim; idir++)
        {
          if (a_problem_domain.isPeriodic(idir))
            periodicTestBox.grow(idir,-1);
        }
    }

  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      int buffer = interp_radius;
      if (dir==0) buffer+=extrap_width;
      m_lo_extrap[dir].define(a_level_domain);
      m_hi_extrap[dir].define(a_level_domain);
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& grid = a_level_domain[dit()];
          const Box interp_box = grow(grid, interp_radius);
          const Box lo_extrap_box = adjCellLo(interp_box, dir, extrap_width);
          const Box hi_extrap_box = adjCellHi(interp_box, dir, extrap_width);
          IntVectSet& local_lo_extrap = m_lo_extrap[dir][dit()];
          IntVectSet& local_hi_extrap = m_hi_extrap[dir][dit()];
          local_lo_extrap.define(a_problem_domain & lo_extrap_box);
          local_hi_extrap.define(a_problem_domain & hi_extrap_box);
          for (it.begin(); it.ok(); ++it)
            {
              const Box& otherBox = a_level_domain[it()];
              if (otherBox.bigEnd(0)+buffer <= grid.smallEnd(0))
                {
                  local_lo_extrap -= otherBox;
                  local_hi_extrap -= otherBox;
                }
              // now do periodic image checking -- also need to remove
              // periodic images from list of cells to be filled,
              // since they will be filled through exchange as well
              // only do this if otherBox and either hi_extrap_box
              // or lo_extrap_box are adjacent to a periodic
              // boundary (using the periodicTestBox)
              if (a_problem_domain.isPeriodic()
                  && (!periodicTestBox.contains(lo_extrap_box)
                      || !periodicTestBox.contains(hi_extrap_box))
                  && !periodicTestBox.contains(otherBox))
                {
                  ShiftIterator shiftIt = a_problem_domain.shiftIterator();
                  IntVect shiftMult(a_problem_domain.domainBox().size());
                  Box shiftedBox(otherBox);
                  for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                    {
                      IntVect shiftVect = shiftMult*shiftIt();
                      shiftedBox.shift(shiftVect);
                      local_lo_extrap -= shiftedBox;
                      local_hi_extrap -= shiftedBox;
                      shiftedBox.shift(-shiftVect);
                    }
                  // can only skip to end in non-periodic case
                } else if (otherBox.smallEnd(0) > grid.bigEnd(0)+buffer)
                  {
                    it.end();
                  }
            }
        }
    }
  m_is_defined = true;
}

// fill outer ghost cells by extrapolation from inner ghost cells
void
ExtrapFillPatch::fillExtrap(
  LevelData<FArrayBox>& a_data,
  int a_dir,
  int a_dest_comp,
  int a_num_comp
  )
{
  CH_assert(m_is_defined);

  const int interp_radius = m_extrap_interval.begin() - 1;
  const BoxLayout& level_domain = a_data.boxLayout();
  DataIterator dit = level_domain.dataIterator();
  for (int comp = a_dest_comp; comp < a_dest_comp + a_num_comp; ++comp)
  {
    for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& grid_domain = level_domain[dit()];
// find coordinate (in the a_dir direction) of inner ghost cells for
// source value
      const int interp_src_lo = grid_domain.smallEnd(a_dir) - interp_radius;
      const int interp_src_hi = grid_domain.bigEnd(a_dir)   + interp_radius;
      FArrayBox& data_fab = a_data[dit()];

      IVSIterator lo_ivsit(m_lo_extrap[a_dir][dit()]);
      for (lo_ivsit.begin(); lo_ivsit.ok(); ++lo_ivsit)
      {
        const IntVect& dest_iv = lo_ivsit();
        IntVect src_iv(dest_iv);
        src_iv[a_dir] = interp_src_lo;
        data_fab(dest_iv,comp) = data_fab(src_iv,comp);
      }
      IVSIterator hi_ivsit(m_hi_extrap[a_dir][dit()]);
      for (hi_ivsit.begin(); hi_ivsit.ok(); ++hi_ivsit)
      {
        const IntVect& dest_iv = hi_ivsit();
        IntVect src_iv(dest_iv);
        src_iv[a_dir] = interp_src_hi;
        data_fab(dest_iv,comp) = data_fab(src_iv,comp);
      }
    }
  }
}

/*
void
ExtrapFillPatch::printIntVectSets() const
{
  DataIterator lit = m_lo_extrap[0].boxLayout().dataIterator();
  for (lit.begin(); lit.ok(); ++lit)
  {
    cout << "grid " << lit().intCode() << ": " << endl;
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
      cout << "lo extrap [" << dir << "]: " << endl;
      cout << m_lo_extrap[dir][lit()] << endl;
      cout << "hi extrap [" << dir << "]: " << endl;
      cout << m_hi_extrap[dir][lit()] << endl;
    }
  }
}
*/
#include "NamespaceFooter.H"
