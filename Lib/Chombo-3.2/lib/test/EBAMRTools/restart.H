#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RESTART_H_
#define _RESTART_H_

#include "UsingNamespace.H"

namespace EBRestart
{

struct InputParams;
struct GeomParams;

void
dumpOldStyle( std::string                          a_filename,
              Vector<std::string>                  a_compNames,
              int                                  a_ghosts,
              const Box&                           a_domainCoar,
              const EBIndexSpace*                  a_ebIndexSpace,
              const Vector<LevelData<EBCellFAB>*>& a_ebvector,
              const Vector<DisjointBoxLayout>&     a_grids,
              int                                  a_nlevs,
              int                                  a_ncomps,
              Real                                 a_dx,
              int                                  a_refRatio );

void
dumpHDF( std::string                            filename,
         Vector<std::string>                    compNames,
         int                                    ghosts,
         const Box&                             domainCoar,
         const EBIndexSpace*                    ebIndexSpace,
         const Vector<LevelData<EBCellFAB>*>&   ebvector,
         const Vector<DisjointBoxLayout>&       a_grids,
         int                                    nlevs,
         int                                    ncomps,
         Real                                   dx,
         int                                    refRatio );

void
loadHDF( int                                    step,
         int                                    ghosts,
         const EBIndexSpace*                    ebIndexSpace,
         const Vector<LevelData<EBCellFAB>* >&  ebvector,
         Vector<DisjointBoxLayout>&             grids,
         const GeomParams&                      geometry,
         const InputParams&                     inputs );


struct InputParams
{
  InputParams( int argc, char**argv, const char* infilename );

  Vector<int>  n_cell;
  Vector<Real> prob_lo;
  Real         prob_hi;
  int          nsteps;
  int          nlevs;
  int          ncomps;
  int          refratio;
  int          blocking_factor;
  int          buffer_size;
  int          max_size;
  int          ghosts;
  int          first_step;
  int          test_mode;
  int          test_pass; // Equals 1 or 2.  Don't set in inputs file!!
  int          verbose;

  //
  // Note that radius and center are members of GeomParams too.  Always use
  // those!  The only reason we have radius and center here too is that they're
  // mentioned in the inputs file, and we want only this class to deal with
  // ParmParse.  So, right at the get-go, the program loads all the inputs into
  // this class, and then GeomParams's constructor takes an InputParams arg.
  // During the course of the program, radius and center change -- but only in
  // GeomParams, those being the authoritative versions after all.
  //
private:
  Real         radius;
  RealVect     center;
  friend class GeomParams;
};


struct GeomParams
{
  GeomParams( const InputParams& inputs );
  GeomParams( const Box&      _domain,
              Real             _dx,
              const RealVect& _origin,
              const RealVect& _center,
              Real            _radius );
  Box      domain;
  Real     dx;
  RealVect origin;
  RealVect center;
  Real     radius;
};


void evolveGeomParams( GeomParams& );


void makeHierarchy(Vector<DisjointBoxLayout>&   dbl,
                   const ProblemDomain&         baseDomain,
                   const IntVectSet&            baseTags,
                   const Vector<int>&           refRatios,
                   const InputParams&           inputs);

// Define a sphere EBIS.
int makeGeometry(EBIndexSpace&      ebIndexSpace,
                 const GeomParams&  geomParams);

void makeLevelSet(Vector<LevelData<FArrayBox>* >&   levelSet,
                  const                             GeomParams&);

// Make the corresponding layout
int makeEBISL(EBISLayout&               ebisl,
              const DisjointBoxLayout&  grids,
              const Box&                domain,
              const EBIndexSpace*       ebIndexSpace,
              const int&                nghost);

void dumpmemoryatexit();

// Apply coarsen n times.
Box coarsenPow( const Box& box, int refratio, int n );

// Generate hdf5 file name.
std::string filename( int step );


/** Stores step-specific, level-specific data we use to verify that loading
 *  was successful.
 *  The sum member isn't necessarily a sum; it could be some other function
 *  of all the data, indeed right now it's the sum of squares.
*/
struct CheckSum
{
  CheckSum() : sum(0), len_reg(0), len_irreg(0)
  {
  }

  long sum;
  int  len_reg;    // # of regular cells
  int  len_irreg;  // # of vofs in irregular cells;
};

// Overall checksum structure: by step, by level.
typedef std::vector<CheckSum> CheckSumVect;
typedef std::map<int, CheckSumVect> CheckSumMap;

void initEBGrids( Vector<DisjointBoxLayout>&            grids,
                  const Vector<LevelData<EBCellFAB>* >& ebvector,
                  const EBIndexSpace*                   ebIndexSpace,
                  int                                   ghosts,
                  const InputParams&                    inputs,
                  const GeomParams&                     geometry );


void fillData( Vector<LevelData<EBCellFAB>* >&  ebvector,
               int                              nlevs,
               int                              ncomps,
               const CheckSumVect&              prev_checksums
             );

void
updateChecksums( Vector<LevelData<EBCellFAB>*>&     ebvector,
                 int                                nlevs,
                 CheckSumVect&                      checksums );


void stepOnce(
  int                             step,
  RefCountedPtr<EBIndexSpace>&    ebIndexSpace,
  RefCountedPtr<EBIndexSpace>&    ebIndexSpace_old,
  Vector<LevelData<EBCellFAB>* >& ebvector,
  Vector<LevelData<EBCellFAB>* >& ebvector_old,
  Vector<DisjointBoxLayout>&      grids,
  int                             ghosts,
  GeomParams&                     geometry,
  const InputParams&              inputs,
  const CheckSumVect&             prev_checksums,
  CheckSumVect&                   curr_checksums);


void mainFunc( const InputParams& );

} // namespace EBRestart

#endif // include guard
