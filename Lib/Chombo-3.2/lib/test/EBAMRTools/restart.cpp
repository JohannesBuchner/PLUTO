#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// Purpose: see comments just above main() function.
//

#include "SphereIF.H"
#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"

#include "EBISLayout.H"
#include "EBAMRIO.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "DebugDump.H"
#include "SPMD.H"
#include "CH_Attach.H"
#include "parstream.H"
#include "restart.H" // In this directory
#include <sstream>
#include <map>
#include <vector>
#include <cstdio>
#include <cerrno>
#include <fstream>

#include "UsingNamespace.H"


/** Produces a bunch of hdf5 files.  Then exec's itself in a mode that causes it
 *  to restart itself from one of those hdf5 files.  Continues for a few time
 *  steps, then stops and compares the last hdf5 file it produced to the hdf5
 *  file from the corresponding time step made in the first run.  If the files
 *  are identical then the test is considered a success.
 *
 *  The test works well enough under mpirun; it runs far enough to report on
 *  its success or failure, but prints some MPI-shutdown-related error messages
 *  after that.
 *  In parallel environments other than mpirun, the program will detect that
 *  there's no mpirun on the PATH, and abort with varying degrees of elegance.
*/
int
main(int argc, char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  {

#ifdef CH_USE_HDF5
    EBRestart::InputParams inputs( argc, argv, "restart.inputs" );
    // All functions and classes specific to this file and restart.H live in
    // their own namespace -- EBRestart.  This is to make it easy to
    // distinguish, in the code below, between Chombo library things, and
    // things specific to this program.

    if ( inputs.test_mode == 0 )
    {
      inputs.first_step = 0;
      EBRestart::mainFunc( inputs );
      //exit(0);
    } else
    {
      //--- test mode ---//
      std::string txtdump( "/tmp/dump.out" );

      if ( inputs.test_pass == 1 )   // Produces inputs.nsteps hdf5 files,
      {                             // then launches second pass.
        int saved_first_step = inputs.first_step;
        inputs.first_step = 0;

        EBRestart::mainFunc( inputs );

        // Save h5dump of last step, and delete all hdf5 files except the
        // inputs.first_step-1th one.
        if ( procID() == 0 )
        {
          std::ostringstream cmd;
          cmd << "h5dump " << EBRestart::filename( inputs.nsteps-1 )
              << " > " << txtdump << ".pass1";
          int ret;
          ret = system( cmd.str().c_str() );
          for ( int s=0; s<inputs.nsteps; ++s )
          {
            if ( s != saved_first_step-1 )
            {
              unlink( EBRestart::filename(s).c_str() );
            }
          }
        }
        // Restart from hdf5 file made at the inputs.first_step-1th step.
        // Proceed with calculation until last step, then compare its hdf5
        // file to the (h5dump of the one) saved from the first test pass.
        typedef char* pchar;
        char** newargv;
        int newargvOffset;
#ifdef CH_MPI
        newargvOffset = 5; // Room for ("mpirun", "-np" "n")
#else
        newargvOffset = 0;
#endif
        newargv = new pchar[argc + 2 + newargvOffset];
        for ( int a=0; a<argc; ++a )
        {
          newargv[a+newargvOffset] = argv[a];
        }
        newargv[argc+newargvOffset] = (char *)"test_pass=2";
        newargv[argc+1+newargvOffset] = 0;

#ifdef CH_MPI
        char hostname[90];
        gethostname(hostname, 90);
        std::ofstream mach;
        mach.open("restart.mach");
        for (int i=0; i<numProc(); i++)
        {
          mach <<hostname<<std::endl;
        }
        mach.close();
        char npbuf[10];
        sprintf( npbuf, "%d", numProc() );
        newargv[1] = "-machinefile";
        newargv[2] = "restart.mach";
        newargv[3] = "-np";
        newargv[4] = npbuf;
        newargv[0] = "mpirun";
#else
        newargv[0] = argv[0];
#endif
        if ( procID() == 0 )
        {
#ifdef CH_MPI
          // This isn't quite a totally perfect cleanup -- there'll be some
          // net_recv errors -- but it's good enough to let things run far
          // to get a chance to compare the hdf5 output of the first and second
          // passes, and thus to determine if this test succeeded.
          // I tried, instead of execvp, fork followed by execvp in the child,
          // but that worked even worse; the second pass didn't get very far
          // at all.  Then I tried starting the second pass with system(): that
          // worked about as well as the execvp strategy.
          // Maybe, to get a totally perfect cleanup, we need to simulate the
          // scoping trick.
#ifdef CH_Linux
          MPI_Finalize(); // Not sure how this is working with this Finalize here
#endif
          if ( 0 != system("which mpirun > /dev/null") )
          {
            pout() << "You don't have mpirun.  This program only works on "
                      << "systems that do.  Sorry.\n";
            exit(0);
          }
#endif
          if ( inputs.verbose )
          {
            pout() << "Restarting at step " << saved_first_step << '\n';
          }
#ifdef CH_PROFILE
          pout() << "CH_PROFILE=TRUE...exiting now (OK).\n";
          exit(0); // execvp() will cause "Profiling timer expired" error.
#endif
          // Only do this if running on Linux platforms.
          // May not be working properly on others
#ifdef CH_Linux
          execvp( newargv[0], newargv );
          perror( "execv" );
#endif
        }
      }
      else
      {
        setPoutBaseName("restart");

        // Second pass, launched by execvp call above.
        // Will compare output against txtdump created above.
        EBRestart::mainFunc( inputs );

        if ( procID() == 0 )
        {
          std::ostringstream cmd;
          cmd << "h5dump " << EBRestart::filename( inputs.nsteps-1 )
              << " > " << txtdump << ".pass2";
          int ret;
          ret = system( cmd.str().c_str() );

          cmd.str("");
          cmd << "cmp " << txtdump << ".pass1 " << txtdump << ".pass2";
          if ( 0 == system( cmd.str().c_str() ) )
          {
            pout() << "restart test passed\n";
          } else
          {
            pout() << "restart test failed.  diff " << txtdump << ".pass1 "
                 << txtdump << ".pass2" << '\n';
            MayDay::Abort("restart test failed");
          }
        }
      }
    }
#endif // CH_USE_HDF5
  } // end scoping trick
#ifdef CH_MPI
    MPI_Finalize();
#endif
}


/**
 * Loops, producing a sequence of hdf5 files step0.hdf5, step1.hdf5, etc.
 * If handed first_step=n, n>0, on the command line, will load step<n>.hdf5
 * and start from that step of the loop.
 *
*/
void
EBRestart::mainFunc( const InputParams& a_inputs )
{
    EBRestart::GeomParams geometry( a_inputs );

    int eekflag = 0;
    RefCountedPtr<EBIndexSpace> ebIndexSpace( new EBIndexSpace );
    RefCountedPtr<EBIndexSpace> ebIndexSpace_old( new EBIndexSpace );

    Vector<LevelData<EBCellFAB>* > ebvector(a_inputs.nlevs);
    Vector<LevelData<EBCellFAB>* > ebvector_old(a_inputs.nlevs);
    for ( int l=0;l<a_inputs.nlevs;++l )
    {
      ebvector[l] = new LevelData<EBCellFAB>;
      ebvector_old[l] = new LevelData<EBCellFAB>;
    }
    CH_assert(eekflag == 0);

    Vector<DisjointBoxLayout> grids(a_inputs.nlevs);

    EBRestart::CheckSumMap checksums; // by step, by level
    checksums[a_inputs.first_step-1] = EBRestart::CheckSumVect(a_inputs.nlevs);

    if ( a_inputs.first_step > 0 ) // Restart.
    {
#ifndef CH_USE_HDF5
      MayDay::Error("inputs.first_step>0 but CH_USE_HDF undefined;"
                    "can't restart without hdf5.");
#endif
      for ( int s=0; s<a_inputs.first_step-1; ++s )
      {
        EBRestart::evolveGeomParams( geometry );
      }
      int eekflag = EBRestart::makeGeometry(*ebIndexSpace, geometry);
      CH_assert( eekflag==0 );

      EBRestart::loadHDF( a_inputs.first_step-1,
                          a_inputs.ghosts,
                          ebIndexSpace,
                          ebvector,
                          grids, geometry, a_inputs );
      EBRestart::updateChecksums( ebvector, a_inputs.nlevs,
                                  checksums[a_inputs.first_step-1] );
      EBRestart::evolveGeomParams( geometry );
    }

    for ( int step=a_inputs.first_step; step<a_inputs.nsteps; ++step )
    {
      checksums[step] = std::vector<EBRestart::CheckSum>(a_inputs.nlevs);
      EBRestart::stepOnce( step,
                           ebIndexSpace,
                           ebIndexSpace_old,
                           ebvector,
                           ebvector_old,
                           grids,
                           a_inputs.ghosts,
                           geometry,
                           a_inputs,
                           checksums[step-1],
                           checksums[step] );
    }

    // Memory cleanup
    for (int l=0; l<a_inputs.nlevs; l++)
    {
      delete ebvector[l];
      delete ebvector_old[l];
    }
}


/*  Dumps the field data and the VOF data, but not the geometry (which must be
 *  recomputed when restarting).
*/
void
EBRestart::dumpHDF( std::string                          a_filename,
                    Vector<std::string>                  a_compNames,
                    int                                  a_ghosts,
                    const Box&                           a_domainCoar,
                    const EBIndexSpace*                  a_ebIndexSpace,
                    const Vector<LevelData<EBCellFAB>*>& a_ebvector,
                    const Vector<DisjointBoxLayout>&     a_grids,
                    int                                  a_nlevs,
                    int                                  a_ncomps,
                    Real                                 a_dx,
                    int                                  a_refRatio )
{
#ifdef CH_USE_HDF5
  RealVect dxvec( RealVect::Unit*a_dx*pow(double(a_refRatio),
                                          double(a_nlevs-1)) );
  HDF5Handle handle;

  createEBFile(handle, a_filename, a_nlevs,
               Vector<int>(a_nlevs,a_refRatio),
               a_domainCoar, dxvec, a_ghosts*IntVect::Unit);
  writeCellCenteredNames(handle, a_compNames);

  for ( int lev=0;lev<a_nlevs;++lev )
  {
    writeCellCentered(handle, lev, a_ebvector[lev]);
  }

  handle.close();
#endif // CH_USE_HDF5
}




/** Essential for restart.
 *  Loads the field data and the VOF data, but not the geometry (which must be
 *  recomputed).
*/
void
EBRestart::loadHDF( int                                   a_step,
                    int                                   a_ghosts,
                    const EBIndexSpace*                   a_ebIndexSpace,
                    const Vector<LevelData<EBCellFAB>* >& a_ebvector,
                    Vector<DisjointBoxLayout>&            a_grids,
                    const EBRestart::GeomParams&          a_geometry,
                    const EBRestart::InputParams&         a_inputs )
{
#ifdef CH_USE_HDF5
  EBRestart::initEBGrids( a_grids, a_ebvector, a_ebIndexSpace, a_ghosts,
                          a_inputs, a_geometry);
  HDF5Handle handle;

  handle.open(filename(a_step), HDF5Handle::OPEN_RDONLY);
  for ( int lev=0; lev<a_inputs.nlevs; ++lev )
  {
    LevelData<EBCellFAB>& field( *a_ebvector[lev] );
    readCellCentered(handle, lev, a_ebIndexSpace, a_ghosts, &field);
  }
  handle.close();
#endif // CH_USE_HDF5
}


/** Called once per time step, from mainFunc(). */
void
EBRestart::stepOnce(
  int                             a_step,
  RefCountedPtr<EBIndexSpace>&    a_ebIndexSpace,
  RefCountedPtr<EBIndexSpace>&    a_ebIndexSpace_old,
  Vector<LevelData<EBCellFAB>* >& a_ebvector,
  Vector<LevelData<EBCellFAB>* >& a_ebvector_old,
  Vector<DisjointBoxLayout>&      a_grids,
  int                             a_ghosts,
  GeomParams&                     a_geometry,
  const InputParams&              a_inputs,
  const EBRestart::CheckSumVect&  a_prevChecksums,
  EBRestart::CheckSumVect&        a_currChecksums
)
{
  int nlevs = a_inputs.nlevs;
  int refRatio = a_inputs.refratio;

  int eekflag = EBRestart::makeGeometry(*a_ebIndexSpace, a_geometry);
  if (eekflag != 0)
  {
    pout() << "Problem in EBRestart::makeGeometry(), eekflag = "
           << eekflag << endl;
  }

  if (a_step > 0)
  {
    // Same grids, new topology
    for ( int l=nlevs-1; l>=0; --l )
    {
      EBISLayout ebisl;
      EBRestart::makeEBISL(ebisl, a_grids[l],
                           coarsenPow(a_geometry.domain,refRatio,nlevs-1-l),
                           a_ebIndexSpace, a_ghosts );
      EBCellFactory factory(ebisl);
      a_ebvector[l]->define(a_grids[l], a_inputs.ncomps, IntVect::Unit*a_ghosts,
                          factory);
    }
    a_ebvector.swap( a_ebvector_old );
  }

  //
  // Regrid
  //
  EBRestart::initEBGrids( a_grids, a_ebvector, a_ebIndexSpace, a_ghosts,
                          a_inputs, a_geometry );

  if ( a_step == 0 )
  {
    EBRestart::fillData( a_ebvector, nlevs, a_inputs.ncomps, a_prevChecksums );
    EBRestart::updateChecksums( a_ebvector, nlevs, a_currChecksums );
  } else
  {
    // Perform second type of remap.  constant EBIndexSpace, but
    // different grids.
    //
    // Coarse grids in this case remain identical.
    a_ebvector_old[0]->copyTo(a_ebvector[0]->interval(),*a_ebvector[0],
                            a_ebvector_old[0]->interval());

    EBRestart::fillData( a_ebvector, nlevs, a_inputs.ncomps, a_prevChecksums );
    EBRestart::updateChecksums( a_ebvector, nlevs, a_currChecksums );
  }


#ifdef CH_USE_HDF5
  // Generate some component names.
  Vector<std::string> compNames(a_inputs.ncomps);
  ostringstream compName;
  for ( int c=0; c<a_inputs.ncomps; ++c )
  {
    compName << "velocity" << c; // Silly, but easy to spot.
    compNames[c] = compName.str();
    compName.str(""); // clears it
  }

  EBRestart::dumpHDF( filename(a_step),
                      compNames,
                      a_ghosts,
                      EBRestart::coarsenPow(a_geometry.domain, refRatio,
                                            nlevs-1),
                      a_ebIndexSpace,
                      a_ebvector,
                      a_grids,
                      nlevs,
                      a_inputs.ncomps,
                      a_geometry.dx,
                      refRatio );

  EBRestart::dumpOldStyle( std::string("oldstyle_") + filename(a_step),
                           compNames,
                           a_ghosts,
                           EBRestart::coarsenPow(a_geometry.domain, refRatio,
                                                 nlevs-1),
                           a_ebIndexSpace,
                           a_ebvector,
                           a_grids,
                           nlevs,
                           a_inputs.ncomps,
                           a_geometry.dx,
                           refRatio );
#endif

  //
  // Change the geometry -- for next step.
  //
  EBRestart::evolveGeomParams( a_geometry );
  a_ebIndexSpace_old.swap( a_ebIndexSpace );
  a_ebvector.swap( a_ebvector_old );
  if ( a_inputs.verbose )
  {
    pout() << "step " << a_step << std::endl;
  }
}


/** InputParams are the things in restart.inputs.
 *  This is the only place we use ParmParse.
*/
EBRestart::InputParams::InputParams(
    int argc, char** argv, const char* a_infilename )
  : n_cell(SpaceDim),
    prob_lo(SpaceDim, 1.0),
    ncomps(0),
    ghosts(0),
    first_step(0),
    test_mode(0),
    test_pass(1),
    verbose(0)
{
  ParmParse pp(argc-1, argv+1, 0, a_infilename);
  pp.getarr("n_cell",n_cell,0,SpaceDim);
  CH_assert(n_cell.size() == SpaceDim);
  pp.getarr("prob_lo",prob_lo,0,SpaceDim);
  pp.get("prob_hi",prob_hi);
  Vector<Real> vcenter;
  pp.getarr("center",vcenter,0,SpaceDim);
  center = RealVect( D_DECL(vcenter[0],vcenter[1],vcenter[2]) );

  pp.get("radius",radius);
  pp.get("nlevs",nlevs);
  pp.get("ncomps",ncomps);
  pp.query("ghosts",ghosts);
  pp.get("refratio",refratio);

  pp.get("nsteps",nsteps);
  pp.get("first_step",first_step);
  CH_assert( first_step < nsteps );

  pp.get("blocking_factor",blocking_factor);
  pp.get("buffer_size",buffer_size);
  pp.get("max_size",max_size);
  pp.get("test_mode",test_mode);

  pp.query("test_pass",test_pass);
  CH_assert( (test_pass == 1) || (test_pass == 2) );

  pp.query("verbose", verbose);
}


EBRestart::GeomParams::GeomParams( const EBRestart::InputParams& a_inputs )
    : domain(Box(IntVect::Zero, IntVect(a_inputs.n_cell) - 1)),
      dx((a_inputs.prob_hi - a_inputs.prob_lo[0])/a_inputs.n_cell[0]),
      origin(a_inputs.prob_lo),
      center(a_inputs.center),
      radius(a_inputs.radius)
{
    CH_assert( IntVect(a_inputs.n_cell) > IntVect::Zero );
}


EBRestart::GeomParams::GeomParams( const Box&      a_domain,
                                   Real            a_dx,
                                   const RealVect& a_origin,
                                   const RealVect& a_center,
                                   Real            a_radius )
    : domain(a_domain),
      dx(a_dx),
      origin(a_origin),
      center(a_center),
      radius(a_radius)
{
}


/** Each time this is called, we move the embedded boundary a little. */
void
EBRestart::evolveGeomParams( EBRestart::GeomParams& a_geometry )
{
  a_geometry.center[0]+= a_geometry.dx/10.0;
  a_geometry.center[1]+= a_geometry.dx/20.0;
  a_geometry.radius += a_geometry.dx/50.0;
}


int
EBRestart::makeGeometry(EBIndexSpace&                a_ebIndexSpace,
                        const EBRestart::GeomParams& a_geomParams)
{
  int eekflag = 0;
  SphereIF outside(a_geomParams.radius, a_geomParams.center, false);

  RealVect vectDx = RealVect::Unit;
  vectDx *= a_geomParams.dx;

  GeometryShop workshop(outside,0,vectDx);

  workshop.m_phase = 0;

  ProblemDomain pdomain( a_geomParams.domain );
  a_ebIndexSpace.define(pdomain,
                        a_geomParams.origin,
                        a_geomParams.dx,
                        workshop);

  return eekflag;
}


/** Initializes the DisjointBoxLayouts, each step. */
void
EBRestart::makeHierarchy(Vector<DisjointBoxLayout>& a_dbl,
                         const ProblemDomain&       a_baseDomain,
                         const IntVectSet&          a_baseTags,
                         const Vector<int>&         a_refRatios,
                         const InputParams&         a_inputs)
{
  a_dbl.resize(a_inputs.nlevs);
  Real fillRatio     = 0.85;
  BRMeshRefine regridder(a_baseDomain, a_refRatios, fillRatio,
                         a_inputs.blocking_factor,
                         a_inputs.buffer_size, a_inputs.max_size);

  Vector<Vector<Box> > oldGrids(a_inputs.nlevs,1), newGrids(a_inputs.nlevs);
  oldGrids[0][0]=a_baseDomain.domainBox();
  for ( int l=1; l<a_inputs.nlevs; ++l )
  {
    oldGrids[l][0]=coarsen(oldGrids[l-1][0], a_refRatios[l-1]);
    // FIXME: is a_refRatios[l-1] correct??
  }

  regridder.regrid(newGrids, a_baseTags, 0, 1, oldGrids);
  // 0 is baselevel, 1 is toplevel.

  Vector<int> procs;
  for (int l=0; l<a_inputs.nlevs; l++)
  {
    newGrids[l].sort();
    LoadBalance(procs, newGrids[l]);
    a_dbl[l] = DisjointBoxLayout(newGrids[l], procs);
  }
}


/** Applies Chombo::coarsen(box,refratio) n times. */
Box
EBRestart::coarsenPow( const Box& a_box, int a_refratio, int n )
{
  if ( n==0 ) return a_box;
  else return coarsenPow( coarsen(a_box,a_refratio), a_refratio, n-1 );
}


/** Copied & pasted from one of the many other places this is defined under
 *  Chombo/lib/test.
*/
int
EBRestart::makeEBISL(EBISLayout&                a_ebisl,
                     const DisjointBoxLayout&   a_grids,
                     const Box&                 a_domain,
                     const EBIndexSpace*        a_ebIndexSpace,
                     const int&                 a_nghost)
{
  CH_assert(a_ebIndexSpace->isDefined());
  a_ebIndexSpace->fillEBISLayout(a_ebisl, a_grids, a_domain, a_nghost);
  return 0;
}


void
EBRestart::initEBGrids( Vector<DisjointBoxLayout>&            a_grids,
                        const Vector<LevelData<EBCellFAB>* >& a_ebvector,
                        const EBIndexSpace*                   a_ebIndexSpace,
                        int                                   a_ghosts,
                        const EBRestart::InputParams&         a_inputs,
                        const EBRestart::GeomParams&          a_geometry )
{
  //
  // Generate the grids and initialize a_ebvector with them.
  //
  int refRatio = a_inputs.refratio;
  int nlevs = a_inputs.nlevs;
  int ncomps = a_inputs.ncomps;

  IntVectSet tags( a_ebIndexSpace->irregCells(2) );
  EBRestart::makeHierarchy(a_grids, EBRestart::coarsenPow(a_geometry.domain,
                                                          refRatio, nlevs-1),
                           tags, Vector<int>(nlevs,refRatio), a_inputs);

  for ( int l=nlevs-1; l>=0; --l )
  {
    EBISLayout ebisl;
    EBRestart::makeEBISL(ebisl, a_grids[l],
                         EBRestart::coarsenPow(a_geometry.domain,
                                               refRatio,nlevs-1-l),
                         a_ebIndexSpace, a_ghosts );
    EBCellFactory ebCellFactory(ebisl);
    a_ebvector[l]->define(a_grids[l], ncomps, IntVect::Unit*a_ghosts,
                          ebCellFactory);
  }
}


std::string
EBRestart::filename( int a_step )
{
  std::ostringstream out;
  out << "step" << a_step << ".hdf5";
  return out.str();
}

//
// Generate field data.  We're not solving any real problem here;
// the only purpose of this is to come up with some numbers that couldn't
// happen by accident, so we can check if the hdf5 I/O is working.
//
// New data is a function of prev_checksum (a vector -- one element per
// level).
// Updates curr_checksum.
//
void
EBRestart::fillData( Vector<LevelData<EBCellFAB>* >& a_ebvector,
                     int                             a_nlevs,
                     int                             a_ncomps,
                     const EBRestart::CheckSumVect&  a_prev_checksums
             )
{
  static int g_i;
  g_i = 0;
  for (int lev=0; lev<a_nlevs; ++lev)
  {
    LevelData<EBCellFAB>& ld( *a_ebvector[lev] );
    DataIterator dit = ld.dataIterator();
    for ( dit.begin(); dit.ok(); ++dit )
    {
      EBCellFAB& ebcf( ld[dit] );

      //
      // Single-valued cells
      //
      BaseFab<Real>& singFab(ebcf.getSingleValuedFAB());
      for ( int i=0; i<singFab.box().numPts()*a_ncomps; ++i )
      {
        long x = (i+a_prev_checksums[lev].sum)
               % (a_prev_checksums[lev].len_irreg+100);
        singFab.dataPtr()[i] = x;
//      singFab.dataPtr()[i] = g_i++;
      }

      //
      // Multivalued cells
      //
      Box box( ld.disjointBoxLayout().get(dit) );
      box.grow( ld.ghostVect() );
      const IntVectSet& irregIVS( ebcf.getEBISBox().boundaryIVS(box) );
      const EBGraph& graph( ebcf.getEBISBox().getEBGraph() );
      Vector<Real> multidat( a_ncomps );
      const IntVectSet& multiIVS( ebcf.getMultiValuedFAB().getIVS() );
      for ( VoFIterator it(irregIVS,graph); it.ok(); ++it )
      {
        if ( multiIVS.contains( it().gridIndex() ) )
        {
          for ( int c=0;c<a_ncomps;++c )
          {
//          multidat[c] = 111000 + g_i++;
            multidat[c] = 111000 + g_i++ + 10*c + a_prev_checksums[lev].sum%71;
          }
          ebcf.assign( &multidat[0], it(), Interval(0,a_ncomps-1) );
        }
      }
    }
  }
}


/** The "checksum" is a function of all the field data at a given step.
 *  That field data, in turn, is a function of the checksum at the previous
 *  step.  That way I can be pretty sure that when this restart thing
 *  produces the expected hdf5 files, it's not just fooling me.
*/
void
EBRestart::updateChecksums( Vector<LevelData<EBCellFAB>*>&  a_ebvector,
                            int                             a_nlevs,
                            EBRestart::CheckSumVect&        a_checksums )
{
  for ( int lev=0; lev<a_nlevs; ++lev )
  {
    LevelData<EBCellFAB>& field( *a_ebvector[lev] );
    DataIterator dit = field.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      EBCellFAB& ebcf((*a_ebvector[lev])[dit]);
      BaseFab<Real>& fab(ebcf.getSingleValuedFAB());
      a_checksums[lev].len_reg += fab.box().numPts();
      a_checksums[lev].len_irreg +=
        ebcf.getEBISBox().getEBGraph().getIrregCells(fab.box()).numPts();
      for ( int i=0; i<fab.box().numPts(); ++i )
      {
        long x = long(fab.dataPtr()[i]);
        a_checksums[lev].sum += x*x;
      }
    }

    //
    // MPI-reduce the various members of a_checksums[lev].
    //
#ifdef CH_MPI
    Real all_proc_sum;
    Real proc_sum( a_checksums[lev].sum );
    if ( MPI_SUCCESS !=
        MPI_Allreduce(&proc_sum, &all_proc_sum, 1, MPI_CH_REAL,
                      MPI_SUM, Chombo_MPI::comm) )
    {
      MayDay::Error("Failure in MPI_Allreduce().");
    }
    a_checksums[lev].sum = long(all_proc_sum);
    // There's no MPI_CH_LONG.  Maybe I should make .sum Real then.

    //
    // Reduce a_checksums[lev].len_reg
    //
    Real all_proc_len_reg;
    Real proc_len_reg( a_checksums[lev].len_reg );
    if ( MPI_SUCCESS !=
        MPI_Allreduce(&proc_len_reg, &all_proc_len_reg, 1, MPI_CH_REAL,
                      MPI_SUM, Chombo_MPI::comm) )
    {
      MayDay::Error("Failure in MPI_Allreduce().");
    }
    a_checksums[lev].len_reg = long(all_proc_len_reg);

    //
    // Reduce a_checksums[lev].len_irreg
    //
    Real all_proc_len_irreg;
    Real proc_len_irreg( a_checksums[lev].len_irreg );
    if ( MPI_SUCCESS !=
        MPI_Allreduce(&proc_len_irreg, &all_proc_len_irreg, 1, MPI_CH_REAL,
                      MPI_SUM, Chombo_MPI::comm) )
    {
      MayDay::Error("Failure in MPI_Allreduce().");
    }
    a_checksums[lev].len_irreg = long(all_proc_len_irreg);
#endif
  }
}


/** So we can look at it in ChomboVis */
void
EBRestart::dumpOldStyle( std::string                          a_filename,
                         Vector<std::string>                  a_compNames,
                         int                                  a_ghosts,
                         const Box&                           a_domainCoar,
                         const EBIndexSpace*                  a_ebIndexSpace,
                         const Vector<LevelData<EBCellFAB>*>& a_ebvector,
                         const Vector<DisjointBoxLayout>&     a_grids,
                         int                                  a_nlevs,
                         int                                  a_ncomps,
                         Real                                 a_dx,
                         int                                  a_refRatio )
{
#ifdef CH_USE_HDF5
  // See if there are any multivalued cells:
  for ( int lev=0; lev<a_nlevs; ++lev )
  {
    writeEBHDF5( a_filename,
                 a_grids,
                 a_ebvector,
                 a_compNames,
                 ProblemDomain(a_domainCoar),
                 a_dx,
                 0.1,
                 0.0,
                 Vector<int>(a_nlevs,a_refRatio),
                 a_nlevs,
                 true,
                 Vector<Real>(a_ncomps,1.234),
                 a_ghosts*IntVect::Unit );
  }
#endif // CH_USE_HDF5
}
