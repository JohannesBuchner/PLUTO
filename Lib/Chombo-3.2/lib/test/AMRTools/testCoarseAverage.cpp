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
#include "CoarseAverage.H"
#include "UsingNamespace.H"
#include "FABView.H"

/// Prototypes:
int
testAverage();

void
parseTestOptions(int argc ,char* argv[]) ;

void
initData(LevelData<FArrayBox>& data,
         Real dx);

void
initDataHarmonic(LevelData<FArrayBox>& data,
                 Real dx);


void
resetCoveredRegions(LevelData<FArrayBox>& crseData,
                    const DisjointBoxLayout& fine_grids,
                    Real value,
                    int nRef);

/// Global variables for handling output:
static const char* pgmname = "testAverage" ;
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

  int status = testAverage() ;
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
testAverage()
{
  int returnVal = 0;
  int ref_ratio = 2;

  const Box coarse_domain(IntVect(D_DECL6(  0,  0, 0,  0,  0, 0)),
                          IntVect(D_DECL6( 11, 11, 3, 11, 11, 3)));
  const ProblemDomain coarse_problem_domain(coarse_domain);
  if (verbose)
  {
    pout() << indent2 << "coarse problem domain = " << coarse_domain << endl;
  }

  const int num_coarse_grids = 2;
  if (verbose)
  {
    pout() << indent2 << "num coarse grids = " << num_coarse_grids << endl;
  }
  Vector<int> assign;
  Vector<Box> coarse_boxes(num_coarse_grids);
  coarse_boxes[0] = Box( IntVect( D_DECL6 (0,  0, 0, 0, 0, 0) ),
                         IntVect( D_DECL6 (5,  9, 3, 5, 9, 3) ));
  coarse_boxes[1] = Box(IntVect( D_DECL6 ( 6,  2, 0,  6,  2, 0) ),
                        IntVect( D_DECL6 (11, 11, 3, 11, 11, 3) ));
  LoadBalance(assign, coarse_boxes);
  const DisjointBoxLayout coarse_grids(coarse_boxes, assign,
                                       coarse_problem_domain);

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


  const int num_comps = 1;
  if (verbose)
  {
    pout() << indent2 << "num comps = " << num_comps << endl;
  }


  // loop over different refinement ratios; want to catch nRef = 2, 4, and
  // at least one other one, since they're implemented separately.
  int nRefinements = 3;
  for (int pass=0; pass<nRefinements; pass++)
    {

      if (verbose)
        {
          pout() << indent2 << "Refinement ratio = " << ref_ratio << endl;
        }


      Box fine_domain = refine(coarse_domain, ref_ratio);
      ProblemDomain fine_problem_domain(fine_domain);

      LoadBalance(assign, fine_boxes);
      const DisjointBoxLayout fine_grids(fine_boxes, assign,
                                         fine_problem_domain);

      if (verbose)
        {
          pout() << indent2 << "fine problem domain = " << fine_domain << endl;

          pout() << indent2 << "fine grids: " << endl;
          LayoutIterator fine_lit = fine_grids.layoutIterator();
          for (fine_lit.begin(); fine_lit.ok(); ++fine_lit)
            {
              pout() << indent2 << fine_lit().intCode() << ": " << fine_grids[fine_lit()] << endl;
            }
        }

      const IntVect ghost_vect =  IntVect::Unit;
      LevelData<FArrayBox> fine_data(fine_grids, num_comps, ghost_vect);
      LevelData<FArrayBox> coarse_data(coarse_grids, num_comps);

      const Real fine_dx = 1.;
      initData(fine_data,
               fine_dx);

      const Real coarse_dx = ref_ratio;

      initData(coarse_data,
               coarse_dx);

      Real testVal = 666.6;
      resetCoveredRegions(coarse_data, fine_grids,
                          testVal, ref_ratio);


      CoarseAverage averager(fine_grids,
                             coarse_grids,
                             num_comps,
                             ref_ratio);

      //
      // test arithmetic averaging
      averager.averageToCoarse(coarse_data, fine_data);

      //
      // compare with exact
      LevelData<FArrayBox> crse_errs(coarse_grids, num_comps, IntVect::Zero);
      initData(crse_errs,
               coarse_dx);

      DataIterator crse_dit = coarse_grids.dataIterator();
      Real err_sum = 0.;
      Real data_sum = 0.;
      for (crse_dit.begin(); crse_dit.ok(); ++crse_dit)
        {
          FArrayBox& crse_err_fab = crse_errs[crse_dit()];
          const FArrayBox& crse_data_fab = coarse_data[crse_dit()];
          crse_err_fab -= crse_data_fab;
          crse_err_fab /= crse_data_fab;
          const Box domain_interior_box = crse_err_fab.box() & coarse_problem_domain;
          err_sum += crse_err_fab.norm(domain_interior_box, 1, 0, num_comps);
          data_sum += crse_data_fab.norm(domain_interior_box, 1, 0, num_comps);
        }
      if (verbose)
        {
          pout() << indent2 << "arithmetic data  sum = " << data_sum << endl;
          pout() << indent2 << "arithmetic error sum = " << err_sum << endl;
        }
#if defined(CH_USE_DOUBLE)
      Real err_sum_tol = 1e-10;
#elif defined(CH_USE_FLOAT)
      Real err_sum_tol = 1e-4;
#endif
      if (err_sum < err_sum_tol)
        {
          if (verbose)
            pout() << indent << pgmname << ": CoarseAverage arithmetic average test passed." << endl;
        }
      else
        {
          returnVal += 1;
          pout() << indent << pgmname << ": CoarseAverage arithmetic average test FAILED!" << endl;
          if (verbose)
            {
              pout() << indent2 << "error sum tolerance = " << err_sum_tol << endl;
            }
          return returnVal ;
        }

      ///-----------------------------------------------------------
      // now test harmonic averaging

      initDataHarmonic(fine_data,
                       fine_dx);

      initDataHarmonic(coarse_data,
                       coarse_dx);

      testVal = 666.6;
      resetCoveredRegions(coarse_data, fine_grids,
                          testVal, ref_ratio);


      //
      // test arithmetic averaging
      averager.averageToCoarseHarmonic(coarse_data, fine_data);

      //
      // compare with exact
      initDataHarmonic(crse_errs,
                       coarse_dx);

      err_sum = 0.;
      data_sum = 0.;
      for (crse_dit.begin(); crse_dit.ok(); ++crse_dit)
        {
          FArrayBox& crse_err_fab = crse_errs[crse_dit()];
          const FArrayBox& crse_data_fab = coarse_data[crse_dit()];
          crse_err_fab -= crse_data_fab;
          crse_err_fab /= crse_data_fab;
          const Box domain_interior_box = crse_err_fab.box() & coarse_problem_domain;
          err_sum += crse_err_fab.norm(domain_interior_box, 1, 0, num_comps);
          data_sum += crse_data_fab.norm(domain_interior_box, 1, 0, num_comps);
        }
      if (verbose)
        {
          pout() << indent2 << "harmonic data  sum = " << data_sum << endl;
          pout() << indent2 << "harmonic error sum = " << err_sum << endl;
        }

      if (err_sum < err_sum_tol)
        {
          if (verbose)
            pout() << indent << pgmname << ": CoarseAverage harmonic average test passed." << endl;
        }
      else
        {
          returnVal += 2;
          pout() << indent << pgmname << ": CoarseAverage harmonic average test FAILED!" << endl;
          if (verbose)
            {
              pout() << indent2 << "error sum tolerance = " << err_sum_tol << endl;
            }
          return returnVal ;
        }
      //-----------------------------------------------------------
      // refine fine_boxes for next pass
      for (int n=0; n<fine_boxes.size(); n++)
        {
          fine_boxes[n].refine(2);
        }
      ref_ratio *= 2;
    } // end loop over passes through different nRef

  return returnVal ;

}

void
initData(LevelData<FArrayBox>& data,
         Real dx)
{
  std::vector<Real> linear_coefs(SpaceDim);
  linear_coefs[0] = 2.;
#if (CH_SPACEDIM > 1)
  linear_coefs[1] = 2;
#if (CH_SPACEDIM > 2)
  linear_coefs[2] = 2;
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
      const Box box = fab.box();
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


void
initDataHarmonic(LevelData<FArrayBox>& data,
                 Real dx)

{
  std::vector<Real> linear_coefs(SpaceDim);
  linear_coefs[0] = 2.;
#if (CH_SPACEDIM > 1)
  linear_coefs[1] = 2.;
#if (CH_SPACEDIM > 2)
  linear_coefs[2] = 2.;
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
      const Box box = fab.box();
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
        val = 1.0/val;
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


void
resetCoveredRegions(LevelData<FArrayBox>& crseData,
                    const DisjointBoxLayout& fine_grids,
                    Real value,
                    int nRef)
{

  int ncomp = crseData.nComp();
  DataIterator crseDit = crseData.dataIterator();
  for (crseDit.begin(); crseDit.ok(); ++crseDit)
    {
      FArrayBox& crseFab = crseData[crseDit];
      LayoutIterator fineLit = fine_grids.layoutIterator();
      for (fineLit.begin(); fineLit.ok(); ++fineLit)
        {
          Box testBox(fine_grids[fineLit]);
          testBox.coarsen(nRef);
          testBox &= crseFab.box();
          if (!testBox.isEmpty())
            {
              crseFab.setVal(value, testBox, 0, ncomp);
            }
        } // end loop over fine boxes
    } // end loop over coarse boxes
}
