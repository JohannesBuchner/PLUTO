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

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "MayDay.H"
using std::cout;
using std::endl;

#include "PiecewiseLinearFillPatchFace.H"

#include "NamespaceHeader.H"

#ifndef copysign
template <class T>
inline
T
copysign (const T& a,
          const T& b)
{
    return ( b >= 0 ) ? ( ( a >= 0 ) ? a : -a) : ( (a >= 0 ) ? -a : a);
}

#endif

const int PiecewiseLinearFillPatchFace::s_stencil_radius = 1;

PiecewiseLinearFillPatchFace::PiecewiseLinearFillPatchFace()
  :
  m_is_defined(false)
{
}


PiecewiseLinearFillPatchFace::~PiecewiseLinearFillPatchFace()
{
}


PiecewiseLinearFillPatchFace::PiecewiseLinearFillPatchFace(
  const DisjointBoxLayout& a_fine_domain,
  const DisjointBoxLayout& a_coarse_domain,
  int a_num_comps,
  const Box& a_crse_problem_domain,
  int a_ref_ratio,
  int a_interp_radius)
  :
  m_is_defined(false)
{
  ProblemDomain crsephysdomain(a_crse_problem_domain);
  define(a_fine_domain,
         a_coarse_domain,
         a_num_comps,
         crsephysdomain,
         a_ref_ratio,
         a_interp_radius);
}


PiecewiseLinearFillPatchFace::PiecewiseLinearFillPatchFace(
  const DisjointBoxLayout& a_fine_domain,
  const DisjointBoxLayout& a_coarse_domain,
  int a_num_comps,
  const ProblemDomain& a_crse_problem_domain,
  int a_ref_ratio,
  int a_interp_radius)
  :
  m_is_defined(false)
{
  define(a_fine_domain,
         a_coarse_domain,
         a_num_comps,
         a_crse_problem_domain,
         a_ref_ratio,
         a_interp_radius);
}


void
PiecewiseLinearFillPatchFace::define(
  const DisjointBoxLayout& a_fine_domain,
  const DisjointBoxLayout& a_coarse_domain,
  int a_num_comps,
  const Box& a_crse_problem_domain,
  int a_ref_ratio,
  int a_interp_radius
  )
{
  ProblemDomain crsephysdomain(a_crse_problem_domain);
  define(a_fine_domain, a_coarse_domain, a_num_comps, crsephysdomain,
         a_ref_ratio, a_interp_radius);
}


void
PiecewiseLinearFillPatchFace::define(
  const DisjointBoxLayout& a_fine_domain,
  const DisjointBoxLayout& a_coarse_domain,
  int a_num_comps,
  const ProblemDomain& a_crse_problem_domain,
  int a_ref_ratio,
  int a_interp_radius
  )
{
  m_ref_ratio = a_ref_ratio;
  m_interp_radius = a_interp_radius;
  m_crse_problem_domain = a_crse_problem_domain;

  ShiftIterator shiftIt = m_crse_problem_domain.shiftIterator();

  if (a_interp_radius != (a_interp_radius/a_ref_ratio)*a_ref_ratio)
  {
    MayDay::Error("PiecewiseLinearFillPatchFace::define: interp_radius must be integral multiple of nRef");
  }

  const ProblemDomain fine_problem_domain = refine(m_crse_problem_domain,
                                                   m_ref_ratio);

  // quick sanity checks
  CH_assert (a_fine_domain.checkPeriodic(fine_problem_domain));
  if (a_coarse_domain.isClosed())
  {
    CH_assert (a_coarse_domain.checkPeriodic(a_crse_problem_domain));

    //
    // create the work array
    DisjointBoxLayout coarsened_fine_domain;
    coarsen ( coarsened_fine_domain,
              a_fine_domain,
              m_ref_ratio );

    const int coarse_slope_radius =
      (m_interp_radius + m_ref_ratio - 1) / m_ref_ratio;
    // note that we add one cell to the radius to account for
    // face centering of data.  Only really need this cell on the
    // high-side, but until i come up with a better way to do this,
    // do it this way
    const int coarse_ghost_radius = coarse_slope_radius + s_stencil_radius +1;
    // (FM commented out)
    // const IntVect coarse_slope = coarse_slope_radius * IntVect::Unit;
    // m_slopes.define(coarsened_fine_domain,
    //                 a_num_comps,
    //                 coarse_slope);
    // const IntVect coarse_ghost = coarse_ghost_radius * IntVect::Unit;
    // m_coarsened_fine_data.define(coarsened_fine_domain,
    //                              a_num_comps,
    //                              coarse_ghost);
    // // Initialize uninitialized memory - a quick hack around the fact that
    // // the ghost cells in coarsened_fine_domain may extend outside the
    // // problem domain (PC, 7/21/00).
    // {
    //   DataIterator dit = coarsened_fine_domain.dataIterator();
    //   for (dit.begin(); dit.ok(); ++dit)
    //     {
    //       m_coarsened_fine_data[dit()].setVal(-666.666);
    //     }
    // }

    // (FM added) coarse slope/ghost/coarsened_fine_domain
    m_coarse_slope = coarse_slope_radius * IntVect::Unit;
    m_coarse_ghost = coarse_ghost_radius * IntVect::Unit;
    m_coarsened_fine_domain= coarsened_fine_domain;

    // allocate intvectsets -- loop over dimensions to do this
    for (int faceDir = 0; faceDir<SpaceDim; faceDir++)
    {
      m_fine_interp[faceDir].define(a_fine_domain);

      for (int dir = 0; dir < SpaceDim; ++dir)
        {
          m_coarse_centered_interp[dir][faceDir].define(coarsened_fine_domain);
          m_coarse_lo_interp[dir][faceDir].define(coarsened_fine_domain);
          m_coarse_hi_interp[dir][faceDir].define(coarsened_fine_domain);
        }
      // Compute intvectsets. We do this mostly in the coarsened domain, with
      // only an intersection with the fine domain at the end.

      // first, create a box which will determine whether a given box
      // adjoins a periodic boundary
      Box periodicTestBox(m_crse_problem_domain.domainBox());
      if (m_crse_problem_domain.isPeriodic())
        {
          for (int idir=0; idir<SpaceDim; idir++)
            {
              if (m_crse_problem_domain.isPeriodic(idir))
                periodicTestBox.grow(idir,-1);
            }
        }

      // also create a fine testBox as well
      Box periodicFineTestBox(fine_problem_domain.domainBox());
      if (m_crse_problem_domain.isPeriodic())
        {
          for (int idir=0; idir<SpaceDim; idir++)
            {
              if (m_crse_problem_domain.isPeriodic(idir))
                periodicFineTestBox.grow(idir,-1);
            }
        }

      // create regions in which to interpolate, then intersect with borders
      // to form regions for centered, one-sided low, and one-sided high
      // differences, per coordinate direction.

      // (dfm, 3/5/01) I think that this logic holds up for face-centered
      // data as well, but will need to check this more carefully...
      DataIterator dit = coarsened_fine_domain.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& fine_box = a_fine_domain[dit()];
          Box coarsened_fine_facebox =
            coarsen(grow(fine_box, m_interp_radius),m_ref_ratio) & m_crse_problem_domain;
          coarsened_fine_facebox.surroundingNodes(faceDir);
          coarsened_fine_facebox.shiftHalf(faceDir,1);
          IntVectSet coarsened_fine_interp(coarsened_fine_facebox);

          // Iterate over boxes in coarsened fine domain, and subtract off
          // from the set of coarse cells from which the fine ghost cells
          // will be interpolated.

          LayoutIterator other_lit = coarsened_fine_domain.layoutIterator();
          for (other_lit.begin(); other_lit.ok(); ++other_lit)
            {
              Box other_coarsened_box = coarsened_fine_domain.get(other_lit());
              other_coarsened_box.surroundingNodes(faceDir);
              other_coarsened_box.shiftHalf(faceDir,1);
              coarsened_fine_interp -= other_coarsened_box;

              // first check to see whether we need to check periodic images
              if (m_crse_problem_domain.isPeriodic()
                  && !periodicTestBox.contains(other_coarsened_box)
                  && !periodicTestBox.contains(coarsened_fine_facebox))
                {
                  IntVect shiftMult(m_crse_problem_domain.domainBox().size());
                  Box shiftedBox(other_coarsened_box);
                  for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                    {
                      IntVect shiftVect = shiftMult*shiftIt();
                      shiftedBox.shift(shiftVect);
                      coarsened_fine_interp -= shiftedBox;
                      shiftedBox.shift(-shiftVect);
                    } // end loop over shift directions
                } // end if we need to check periodic images
            } // end loop over other coarsened fine boxes


          // Now that we have the coarsened cells required for interpolation,
          // construct IntvectSets specifying the one-sided and centered
          // stencil locations.

          for (int dir = 0; dir < SpaceDim; ++dir)
          {
            IntVectSet& coarse_centered_interp
              = m_coarse_centered_interp[dir][faceDir][dit()];
            coarse_centered_interp = coarsened_fine_interp;

            IntVectSet& coarse_lo_interp = m_coarse_lo_interp[dir][faceDir][dit()];
            coarse_lo_interp = coarse_centered_interp;
            coarse_lo_interp.shift(BASISV(dir));

            IntVectSet& coarse_hi_interp = m_coarse_hi_interp[dir][faceDir][dit()];
            coarse_hi_interp = coarse_centered_interp;
            coarse_hi_interp.shift(-BASISV(dir));

            // We iterate over the coarse grids and subtract them off of the
            // one-sided stencils.
            LayoutIterator coarse_lit = a_coarse_domain.layoutIterator();
            for (coarse_lit.begin();coarse_lit.ok();++coarse_lit)
              {
                Box bx = a_coarse_domain.get(coarse_lit());
                bx.surroundingNodes(faceDir);
                bx.shiftHalf(faceDir,1);
                coarse_lo_interp -= bx;
                coarse_hi_interp -= bx;
                // once again, need to do periodic images too
                if (m_crse_problem_domain.isPeriodic()
                    && !periodicTestBox.contains(bx)
                    && !periodicTestBox.contains(coarsened_fine_facebox))
                {
                  IntVect shiftMult(m_crse_problem_domain.domainBox().size());
                  Box shiftedBox(bx);
                  for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                  {
                    IntVect shiftVect = shiftMult*shiftIt();
                    shiftedBox.shift(shiftVect);
                    coarse_lo_interp -= shiftedBox;
                    coarse_hi_interp -= shiftedBox;
                    shiftedBox.shift(-shiftVect);
                  }
                } // end if necessary to check for periodic images
              } // end loop over coarse grids

            coarse_lo_interp.shift(-BASISV(dir));
            coarse_hi_interp.shift(BASISV(dir));
            coarse_centered_interp -= coarse_lo_interp;
            coarse_centered_interp -= coarse_hi_interp;
          } // end loop over directions

          // Finally, we construct the fine cells that are going to be
          // interpolated by doing the same sort of logic

          IntVectSet& fine_interp = m_fine_interp[faceDir][dit()];

          // create fine box
          Box fine_faceBox(fine_box);
          fine_faceBox.grow(m_interp_radius);
          fine_faceBox &= fine_problem_domain;
          fine_faceBox.surroundingNodes(faceDir);
          fine_faceBox.shiftHalf(faceDir, 1);
          fine_interp.define(fine_faceBox);
          // now loop over fine-grid interiors and remove from IVS
          LayoutIterator fine_lit = a_fine_domain.layoutIterator();
          for (fine_lit.begin();fine_lit.ok();++fine_lit)
            {
              Box bx = a_fine_domain.get(fine_lit());
              bx.surroundingNodes(faceDir);
              bx.shiftHalf(faceDir,1);
              fine_interp -= bx;

              // also may need to do periodic checking here as well
              // this looks a lot like what we do on the coarse level
              // first check to see if we need to check other periodic images
              if (fine_problem_domain.isPeriodic()
                  && !periodicFineTestBox.contains(fine_box)
                  && !periodicTestBox.contains(bx))
                {
                  IntVect shiftMult(fine_problem_domain.domainBox().size());
                  Box shiftedBox(bx);
                  for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                    {
                      IntVect shiftVect = shiftMult*shiftIt();
                      shiftedBox.shift(shiftVect);
                      fine_interp -= shiftedBox;
                      shiftedBox.shift(-shiftVect);
                    } // end loop over shift directions
                } // end if we need to check periodic images

            } // end loop over fine boxes at this level

        } // end outer loop over fine boxes
    } // end loop over face directions
    m_is_defined = true;
  } // end if coarse domain is defined
}

bool
PiecewiseLinearFillPatchFace::isDefined() const
{
  return ( m_is_defined );
}


// fill the interpolation region of the fine level ghost cells
void
PiecewiseLinearFillPatchFace::fillInterp(
  LevelData<FluxBox>& a_fine_data,
  const LevelData<FluxBox>& a_old_coarse_data,
  const LevelData<FluxBox>& a_new_coarse_data,
  Real a_time_interp_coef,
  int a_src_comp,
  int a_dest_comp,
  int a_num_comp
  )
{
  CH_assert (m_is_defined);
  CH_assert (a_time_interp_coef >= 0.);
  CH_assert (a_time_interp_coef <= 1.);

  // need to do exchange to ensure that adjoining grids get handled correctly
  Interval fineComps(a_src_comp, a_src_comp+a_num_comp-1);
  a_fine_data.exchange(fineComps);

  // FM modification: define cfd
  LevelData<FluxBox> m_coarsened_fine_data(m_coarsened_fine_domain,
                                           a_num_comp,
                                           m_coarse_ghost);
//
// time interpolation of coarse level data, to coarsened fine level work array
  timeInterp(m_coarsened_fine_data,
             a_old_coarse_data,
             a_new_coarse_data,
             a_time_interp_coef,
             a_src_comp,
             a_dest_comp,
             a_num_comp);
//
// piecewise contant interpolation, from coarsened fine level work
// array, to fine level
  fillConstantInterp(a_fine_data,
                     m_coarsened_fine_data,
                     a_src_comp,
                     a_dest_comp,
                     a_num_comp);
//
// increment fine level data with per-direction linear terms
  for (int dir = 0; dir < SpaceDim; ++dir)
  {
    // computeSlopes(dir,
    //               a_src_comp,
    //               a_num_comp);
    incrementLinearInterpTangential(a_fine_data,
                                    m_coarsened_fine_data,
                                    dir,
                                    a_src_comp,
                                    a_dest_comp,
                                    a_num_comp);
  }

  incrementLinearInterpNormal(a_fine_data,
                              a_src_comp,
                              a_dest_comp,
                              a_num_comp);
}


//
// time interpolation of coarse level data, to coarsened fine level work array
void
PiecewiseLinearFillPatchFace::timeInterp(
  LevelData<FluxBox>& m_coarsened_fine_data,
  const LevelData<FluxBox>& a_old_coarse_data,
  const LevelData<FluxBox>& a_new_coarse_data,
  Real a_time_interp_coef,
  int a_src_comp,
  int a_dest_comp,
  int a_num_comp
  )
{
  Interval src_interval (a_src_comp,  a_src_comp  + a_num_comp - 1);
  Interval dest_interval(a_dest_comp, a_dest_comp + a_num_comp - 1);

  if ( (a_old_coarse_data.boxLayout().size() == 0)  &&
       (a_new_coarse_data.boxLayout().size() == 0) )
  {
    MayDay::Error ( "PiecewiseLinearFillPatchFace::fillInterp: no old coarse data and no new coarse data" );
  }
  else if ( (a_time_interp_coef == 1.) ||
            (a_old_coarse_data.boxLayout().size() == 0) )
  {
// old coarse data is absent, or fine time level is the new coarse time level
    a_new_coarse_data.copyTo( src_interval,
                              m_coarsened_fine_data,
                              dest_interval
                              );
  }
  else if ( (a_time_interp_coef == 0.) ||
            (a_new_coarse_data.boxLayout().size() == 0) )
  {
// new coarse data is absent, or fine time level is the old coarse time level
    a_old_coarse_data.copyTo( src_interval,
                              m_coarsened_fine_data,
                              dest_interval
                              );
  }
  else
  {
// linearly interpolate between old and new time levels
    a_new_coarse_data.copyTo( src_interval,
                              m_coarsened_fine_data,
                              dest_interval
                              );
    const DisjointBoxLayout&
      coarsened_fine_layout = m_coarsened_fine_data.disjointBoxLayout();
    LevelData<FluxBox>
      tmp_coarsened_fine_data(coarsened_fine_layout,
                              m_coarsened_fine_data.nComp(),
                              m_coarsened_fine_data.ghostVect());

// Initialize uninitialized memory - a quick hack around the fact that
// the ghost cells in coarsened_fine_domain may extend outside the
// problem domain (PC, 7/21/00).
    {
      DataIterator dit = coarsened_fine_layout.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
        tmp_coarsened_fine_data[dit()].setVal(-666.666);
      }
    }
    a_old_coarse_data.copyTo( src_interval,
                              tmp_coarsened_fine_data,
                              dest_interval
                              );
    DataIterator dit = coarsened_fine_layout.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& coarsened_fine_fb = m_coarsened_fine_data[dit()];
      FluxBox& tmp_coarsened_fine_fb = tmp_coarsened_fine_data[dit()];
      for (int dir=0; dir<SpaceDim; dir++)
      {
        coarsened_fine_fb[dir] *= a_time_interp_coef;
        tmp_coarsened_fine_fb[dir] *= (1.- a_time_interp_coef);
        coarsened_fine_fb[dir] += tmp_coarsened_fine_fb[dir];
      }
    }
  }
}


// fill the fine interpolation region piecewise-constantly
void
PiecewiseLinearFillPatchFace::fillConstantInterp(
  LevelData<FluxBox>& a_fine_data,
  const LevelData<FluxBox>& m_coarsened_fine_data,
  int a_src_comp,
  int a_dest_comp,
  int a_num_comp
  )
  const
{
  DataIterator dit = a_fine_data.boxLayout().dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
  {
    FluxBox& fine_flux = a_fine_data[dit()];
    const FluxBox& coarse_flux = m_coarsened_fine_data[dit()];
    for (int faceDir=0; faceDir<SpaceDim; faceDir++)
    {
      FArrayBox& fine_fab = fine_flux[faceDir];
      const FArrayBox& coarse_fab = coarse_flux[faceDir];
      const IntVectSet& local_fine_interp = m_fine_interp[faceDir][dit()];
      IVSIterator ivsit(local_fine_interp);
      for (ivsit.begin(); ivsit.ok(); ++ivsit)
      {
        IntVect fine_iv = ivsit();
        // fine-coarse referencing is cell-centered, so convert
        // edge-coordinate to cell, coarsen, then convert back
        // again
        fine_iv.shift(faceDir,-1);
        IntVect coarse_iv = coarsen(fine_iv, m_ref_ratio);
        coarse_iv.shift(faceDir,+1);
        fine_iv.shift(faceDir,+1);
        int coarse_comp = a_src_comp;
        int fine_comp   = a_dest_comp;
        for (; coarse_comp < a_src_comp + a_num_comp; ++fine_comp, ++coarse_comp)
          fine_fab(fine_iv, fine_comp) = coarse_fab(coarse_iv, coarse_comp);
      }
    }
  }
}

// compute slopes at the coarse interpolation sites in the specified direction
void
PiecewiseLinearFillPatchFace::computeSlopes(FArrayBox& slope_fab,
                                            const FArrayBox& data_fab,
                                            const IntVectSet& local_centered_interp,
                                            const IntVectSet& local_lo_interp,
                                            const IntVectSet& local_hi_interp,
                                            int a_dir,
                                            int a_src_comp,
                                            int a_num_comp) const
{

  // loop over face directions -- don't bother to do normal direction
  // for (int faceDir=0; faceDir<SpaceDim; faceDir++)
  {
    // // Initialize uninitialized memory - a quick hack around the fact that
    // // the ghost cells in coarsened_fine_domain may extend outside the
    // // problem domain (PC, 7/21/00).
    // {
    //   DataIterator dit = m_slopes.boxLayout().dataIterator();
    //   for (dit.begin(); dit.ok(); ++dit)
    //     {
    //       m_slopes[dit()][faceDir].setVal(-666.666);
    //     }
    // }
    // if (faceDir != a_dir)
    {
      // DataIterator dit
      //   = m_coarsened_fine_data.boxLayout().dataIterator();
      // //  for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
      // //  {
      // for (dit.begin(); dit.ok(); ++dit)
      {
        // const FArrayBox& data_fab = m_coarsened_fine_data[dit()][faceDir];
        // FArrayBox& slope_fab = m_slopes[dit()][faceDir];
        // const IntVectSet& local_centered_interp
        //   = m_coarse_centered_interp[a_dir][faceDir][dit()];
        // const IntVectSet& local_lo_interp
        //   = m_coarse_lo_interp[a_dir][faceDir][dit()];
        // const IntVectSet& local_hi_interp
        //   = m_coarse_hi_interp[a_dir][faceDir][dit()];
        //
        // van leer limited central difference
        IVSIterator centered_ivsit(local_centered_interp);
        for (centered_ivsit.begin(); centered_ivsit.ok(); ++centered_ivsit)
        {
          const IntVect& iv = centered_ivsit();
          const IntVect ivlo = iv - BASISV(a_dir);
          const IntVect ivhi = iv + BASISV(a_dir);
          for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
          {
            Real dlo = data_fab(iv,comp)   - data_fab(ivlo,comp);
            Real dhi = data_fab(ivhi,comp) - data_fab(iv,comp);
            Real dcenter = .5 * (dlo + dhi);
            Real dlim = 2.* Min(Abs(dlo),Abs(dhi));
            if (dlo*dhi < 0.) dlim = 0.;
            dlim = copysign(Min(Abs(dcenter),dlim), dcenter);
            slope_fab(iv,comp) = dlim;
          }
        }
        //
        // one-sided difference (low)
        IVSIterator lo_ivsit(local_lo_interp);
        for (lo_ivsit.begin(); lo_ivsit.ok(); ++lo_ivsit)
        {
          const IntVect& iv = lo_ivsit();
          const IntVect ivlo = iv - BASISV(a_dir);
          for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
          {
            Real dlo = data_fab(iv,comp) - data_fab(ivlo,comp);
            slope_fab(iv,comp) = dlo;
          }
        }
        //
        // one-sided difference (high)
        IVSIterator hi_ivsit(local_hi_interp);
        for (hi_ivsit.begin(); hi_ivsit.ok(); ++hi_ivsit)
        {
          const IntVect& iv = hi_ivsit();
          const IntVect ivhi = iv + BASISV(a_dir);
          for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
          {
            Real dhi = data_fab(ivhi,comp) - data_fab(iv,comp);
            slope_fab(iv,comp) = dhi;
          }
        }
      }
    } // end if not normal direction
  } // end loop over face directions
}

// increment the fine interpolation sites with linear term for the
// specified coordinate direction
void
PiecewiseLinearFillPatchFace::incrementLinearInterpTangential(
  LevelData<FluxBox>& a_fine_data,
  const LevelData<FluxBox>& m_coarsened_fine_data,
  int a_dir,
  int a_src_comp,
  int a_dest_comp,
  int a_num_comp
  )
  const
{

  for (int faceDir=0; faceDir<SpaceDim; faceDir++)
  {
    // once again, only do tangential directions
    if (faceDir != a_dir)
    {
      DataIterator dit = a_fine_data.boxLayout().dataIterator();
      //const FArrayBox& slope_fab = m_slopes[dit()][faceDir];

      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox slope_fab(surroundingNodes(grow(m_coarsened_fine_domain[dit],
                                                    m_coarse_slope),
                                               faceDir),
                              a_num_comp);

          computeSlopes(slope_fab,
                        //const FArrayBox& data_fab,
                        m_coarsened_fine_data[dit()][faceDir],
                        //const IntVectSet& local_centered_interp,
                        m_coarse_centered_interp[a_dir][faceDir][dit()],
                        //const IntVectSet& local_lo_interp,
                        m_coarse_lo_interp[a_dir][faceDir][dit()],
                        //const IntVectSet& local_hi_interp,
                        m_coarse_hi_interp[a_dir][faceDir][dit()],
                        a_dir,
                        a_src_comp,
                        a_num_comp);

          FArrayBox& fine_data_fab = a_fine_data[dit()][faceDir];
          const IntVectSet& fine_interp = m_fine_interp[faceDir][dit()];
          IVSIterator ivsit(fine_interp);
          for (ivsit.begin(); ivsit.ok(); ++ivsit)
            {
              const IntVect& fine_iv = ivsit();
              const IntVect coarse_iv = coarsen(fine_iv,m_ref_ratio);
              const int offset = fine_iv[a_dir] - m_ref_ratio * coarse_iv[a_dir];
              Real interp_coef = -.5 + (offset +.5) / m_ref_ratio;
              int coarse_comp = a_src_comp;
              int fine_comp   = a_dest_comp;
              for (; coarse_comp < a_src_comp+a_num_comp; ++coarse_comp,++fine_comp)
                {
                  fine_data_fab(fine_iv,fine_comp)
                    += interp_coef * slope_fab(coarse_iv,coarse_comp);
                }
        }
      }
    } // end if not a normal direction
  } // end loop over face directions
}

// interpolate fine-grid values in the normal directions
// by interpolating between fine-grid edges which overlie
// coarse grid edges
void
PiecewiseLinearFillPatchFace::incrementLinearInterpNormal(
  LevelData<FluxBox>& a_fine_data,
  int a_src_comp,
  int a_dest_comp,
  int a_num_comp
  )
  const
{
  DataIterator dit= a_fine_data.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          // loop over fine locations, picking out which faces
          // do not overlie coarse-grid faces

          //    Box interiorFaceBox(IntVect::TheZeroVector(),
          //                    (m_ref_ratio-1)*IntVect::TheUnitVector());
          //interiorFaceBox.surroundingNodes(faceDir);
          //interiorFaceBox.grow(faceDir, -1);

          IntVect hiShift(IntVect::Zero);
          hiShift.setVal(faceDir,m_ref_ratio);


      FArrayBox& fine_data_fab = a_fine_data[dit()][faceDir];
      const IntVectSet& fine_interpIVS = m_fine_interp[faceDir][dit()];
      IVSIterator ivsit(fine_interpIVS);
      for (ivsit.begin(); ivsit.ok(); ++ivsit)
      {
        // check to see whether this face is a refined coarse-grid
        // face
        const IntVect& fine_iv = ivsit();
        int loEdgeComp = fine_iv[faceDir];
        // this looks a bit weird in order to be correct for
        // negative indices
        loEdgeComp = (loEdgeComp < 0) ? (-abs(loEdgeComp+1)/m_ref_ratio-1) :
                                         loEdgeComp/m_ref_ratio;
        loEdgeComp *= m_ref_ratio;
        if (loEdgeComp != fine_iv[faceDir])
        {
          // define fine-grid faces between which we're interpolating
          IntVect loVect(fine_iv);
          loVect.setVal(faceDir,loEdgeComp);
          IntVect hiVect = loVect + hiShift;
          Real fraction = fine_iv[faceDir] - loEdgeComp;
          fraction = fraction/m_ref_ratio;
          for (int comp=a_dest_comp;comp<a_dest_comp+a_num_comp; ++comp)
          {
            fine_data_fab(fine_iv, comp)
              = (1.0-fraction)*fine_data_fab(loVect,comp);
            fine_data_fab(fine_iv,comp)
              += fraction*fine_data_fab(hiVect,comp);
          }
        } // end if fine face doesn't overlie a coarse face
      } // end loop over fine faces to be filled
    } // end loop over fine grids
  } // end loop over face directions
}

void
PiecewiseLinearFillPatchFace::printIntVectSets() const
{

  for (int faceDir=0; faceDir<SpaceDim; faceDir++)
  {
    cout << "face direction = " << faceDir << endl;

    DataIterator lit = m_fine_interp[faceDir].boxLayout().dataIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      cout << "grid " << lit().intCode() << ": " << endl;
      cout << "fine ivs" << endl;
      cout << m_fine_interp[faceDir][lit()] << endl;

      for (int dir = 0; dir < SpaceDim; ++dir)
        {
          cout << "coarse centered ivs [" << dir << "]: " << endl;
          cout << m_coarse_centered_interp[dir][faceDir][lit()] << endl;
          cout << "coarse lo ivs [" << dir << "]: " << endl;
          cout << m_coarse_lo_interp[dir][faceDir][lit()] << endl;
          cout << "coarse hi ivs [" << dir << "]: " << endl;
          cout << m_coarse_hi_interp[dir][faceDir][lit()] << endl;
        }

    }
  } // end loop over face directions
}


#include "NamespaceFooter.H"

