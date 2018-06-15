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
#include <vector>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "BoxIterator.H"
#include "PiecewiseLinearFillPatch.H"
#include "UsingNamespace.H"
#include "FABView.H"
#include "DebugDump.H"

/// Prototypes:
int
testPWL();

void
parseTestOptions(int argc ,char* argv[]) ;

void
initData(LevelData<FArrayBox>& data,
         const IntVect& ghost,
         Real dx,
         Real time_interp_coef);

/// Global variables for handling output:
static const char* pgmname = "testPWL" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  // Make this test pass automatically if DIM > 3, to avoid FAIL message.
  int status = (SpaceDim <= 3) ? testPWL() : 0;
  pout() << indent << pgmname ;
  if ( status == 0 )
    {
      pout() << " passed." << endl;
    }
  else
    {
      pout() << " failed with result code " << status << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status ;
}

int
testPWL()
{
  const int ref_ratio = 2;
  if (verbose)
  {
    pout() << indent2 << "refinement ratio = " << ref_ratio << endl;
  }
  const Box coarse_problem_domain(IntVect(D_DECL6(  0,  0, 0,  0,  0, 0)),
                                  IntVect(D_DECL6( 11, 11, 3, 11, 11, 3)));
  if (verbose)
  {
    pout() << indent2 << "coarse problem domain = " << coarse_problem_domain << endl;
  }
  const Box fine_problem_domain = refine(coarse_problem_domain, ref_ratio);
  if (verbose)
  {
    pout() << indent2 << "fine problem domain = " << fine_problem_domain << endl;
  }
  const int num_coarse_grids = 2;
  if (verbose)
  {
    pout() << indent2 << "num coarse grids = " << num_coarse_grids << endl;
  }
  Vector<int> assign;
  Vector<Box> coarse_boxes(num_coarse_grids);
  coarse_boxes[0] = Box( IntVect( D_DECL6 (0,  0, 0,  0,  0, 0) ),
                         IntVect( D_DECL6 (5, 11, 3,  5, 11, 3) ));
  coarse_boxes[1] = Box(IntVect( D_DECL6 ( 6,  0, 0,  6,  0, 0) ),
                        IntVect( D_DECL6 (11, 11, 3, 11, 11, 3) ));
  LoadBalance(assign, coarse_boxes);
  const DisjointBoxLayout coarse_grids(coarse_boxes, assign);

  if (verbose)
  {
    pout() << indent2 << "coarse grids: " << endl;
    LayoutIterator coarse_lit = coarse_grids.layoutIterator();
    for (coarse_lit.begin(); coarse_lit.ok(); ++coarse_lit)
    {
      pout() << indent2 << coarse_lit().intCode() << ": " << coarse_grids[coarse_lit()] << endl;
    }
  }
  int num_fine_grids = 3;
  if (verbose)
  {
    pout() << indent2 << "num fine grids = " << num_fine_grids << endl;
  }

  Vector<Box> fine_boxes(num_fine_grids);
#if (CH_SPACEDIM == 1)
  // 1d requires a bit of special handling
  fine_boxes[0] = Box(IntVect( 2 ),
                      IntVect( 5 ));
  fine_boxes[1] = Box(IntVect( 6 ),
                      IntVect( 11 ));
  fine_boxes[2] = Box(IntVect( 20 ),
                      IntVect( 23 ));
#else
  fine_boxes[0] = Box(IntVect( D_DECL6 (  2,  4, 0,  2,  4, 0) ),
                      IntVect( D_DECL6 (  7, 11, 7,  7, 11, 7) ));
  fine_boxes[1] = Box(IntVect( D_DECL6 (  6, 12, 0,  6, 12, 0) ),
                      IntVect( D_DECL6 ( 23, 17, 7, 23, 17, 7) ));
  fine_boxes[2] = Box(IntVect( D_DECL6 ( 20,  6, 0, 20,  6, 0) ),
                      IntVect( D_DECL6 ( 23, 11, 7, 23, 11, 7) ));
#endif

  LoadBalance(assign, fine_boxes);
  const DisjointBoxLayout fine_grids(fine_boxes, assign);

  if (verbose)
  {
    pout() << indent2 << "fine grids: " << endl;
    LayoutIterator fine_lit = fine_grids.layoutIterator();
    for (fine_lit.begin(); fine_lit.ok(); ++fine_lit)
    {
      pout() << indent2 << fine_lit().intCode() << ": " << fine_grids[fine_lit()] << endl;
    }
  }
  const int num_comps = 1;
  if (verbose)
  {
    pout() << indent2 << "num comps = " << num_comps << endl;
  }
  const int interp_radius = 1;
  if (verbose)
  {
    pout() << indent2 << "interp radius = " << interp_radius << endl;
  }
  PiecewiseLinearFillPatch pwl(fine_grids,
                               coarse_grids,
                               num_comps,
                               coarse_problem_domain,
                               ref_ratio,
                               interp_radius);
// uncomment the following functions for more verbose deubigging output
//  pwl.printIntVectSets();
//  pwl.writeIntVectSetsPlotfiles();

  const IntVect ghost_vect = interp_radius * IntVect::Unit;

  LevelData<FArrayBox> fine_data(fine_grids, num_comps, ghost_vect);
  LevelData<FArrayBox> old_coarse_data(coarse_grids, num_comps);
  LevelData<FArrayBox> new_coarse_data(coarse_grids, num_comps);

  const Real time_interp_coef = .25;
  if (verbose)
  {
    pout() << indent2 << "time interp coef = " << time_interp_coef << endl;
  }
  const Real fine_dx = 1.;
  initData(fine_data,
           IntVect::Zero,
           fine_dx,
           time_interp_coef);
  const Real coarse_dx = ref_ratio;
  initData(old_coarse_data,
           IntVect::Zero,
           coarse_dx,
           0.);
  initData(new_coarse_data,
           IntVect::Zero,
           coarse_dx,
           1.);
//
// test piecewise linear interpolation
  const int src_comp = 0;
  const int dest_comp = 0;

  pwl.fillInterp(fine_data,
                 old_coarse_data,
                 new_coarse_data,
                 time_interp_coef,
                 src_comp,
                 dest_comp,
                 num_comps);

  fine_data.exchange(fine_data.interval());
//
// compare with exact
  LevelData<FArrayBox> fine_errs(fine_grids, num_comps, ghost_vect);
  initData(fine_errs,
           ghost_vect,
           fine_dx,
           time_interp_coef);
  DataIterator fine_dit = fine_grids.dataIterator();
  Real err_sum = 0.;
  Real data_sum = 0.;
  for (fine_dit.begin(); fine_dit.ok(); ++fine_dit)
  {
    FArrayBox& fine_err_fab = fine_errs[fine_dit()];
    const FArrayBox& fine_data_fab = fine_data[fine_dit()];
    fine_err_fab -= fine_data_fab;
    fine_err_fab /= fine_data_fab;
    const Box domain_interior_box = fine_err_fab.box() & fine_problem_domain;
    err_sum += fine_err_fab.norm(domain_interior_box, 1, 0, num_comps);
    data_sum += fine_data_fab.norm(domain_interior_box, 1, 0, num_comps);
  }
  if (verbose)
  {
    pout() << indent2 << "data  sum = " << data_sum << endl;
    pout() << indent2 << "error sum = " << err_sum << endl;
  }
#if defined(CH_USE_DOUBLE)
  Real err_sum_tol = 1e-10;
#elif defined(CH_USE_FLOAT)
  Real err_sum_tol = 1e-4;
#endif
  if (err_sum < err_sum_tol)
    {
      if (verbose)
        pout() << indent << pgmname << ": piecewise linear coarse-fine interpolation test passed." << endl;
      return 0 ;
    }
  else
    {
      pout() << indent << pgmname << ": piecewise linear coarse-fine interpolation test FAILED!!!" << endl;
    if (verbose)
    {
      pout() << indent2 << "error sum tolerance = " << err_sum_tol << endl;
    }
    return 1 ;
  }
}

void
initData(LevelData<FArrayBox>& data,
         const IntVect& ghost,
         Real dx,
         Real time_interp_coef)
{
  std::vector<Real> linear_coefs(SpaceDim);
  linear_coefs[0] = 2.;
#if (CH_SPACEDIM > 1)
  linear_coefs[1] = .02;
#if (CH_SPACEDIM > 2)
  linear_coefs[2] = .0002;
#endif
#endif
  const int num_comps = data.nComp();
// set fine level data
  const BoxLayout& grids = data.boxLayout();
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    FArrayBox& fab = data[dit()];
    fab.setVal(-666.666);
  }
  for (int comp = 0; comp < num_comps; ++comp)
  {
    for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& fab = data[dit()];
      const Box box = grow(grids[dit()], ghost);
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit)
      {
        const IntVect& iv = bit();
        Real& val = fab(iv,comp);
        val = 0.;
        for (int dir = 0; dir < SpaceDim; ++dir)
        {
          val += dx * linear_coefs[dir] * (iv[dir]+.5);
        }
        val += comp;
      }
    }
  }
}

///
// Parse the standard test options (-p -q) out of the command line.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
  {
    if ( argv[i][0] == '-' ) //if it is an option
    {
      // compare 3 chars to differentiate -x from -xx
      if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
      {
        verbose = true ;
        // argv[i] = "" ;
      }
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // argv[i] = "" ;
      }
    }
  }
  return ;
}
