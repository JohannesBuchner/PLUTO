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
using std::endl;

#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "CH_HDF5.H"
#include "CH_Attach.H"
#include "parstream.H"
#include "BoxIterator.H"
#include "FABView.H"

#include "NewPoissonOp.H"
#include "AMRPoissonOp.H"
#include "BCFunc.H"
#include "BiCGStabSolver.H"
#include "CH_Timer.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testMultiGrid" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

///
// Parse the standard test options (-v -q) out of the command line.
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

//   u = x*x + y*y + z*z
//   du/dx = 2*x
//   du/dy = 2*y
//   du/dz = 2*z
//   Laplace(u) = 2*CH_SPACEDIM

//XX -- this is just plain bad
//XX#ifndef __USE_GNU
//XX#define __USE_GNU
//XX#endif
//XX#include <fenv.h>
//XX#undef  __USE_GNU
//#include <fenv.h>

static int nCells = 256;
static Box domain = Box(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(nCells-1,nCells-1,nCells-1)));
static Real dx = 1.6/nCells;
static int blockingFactor = 8;
static Real xshift = 0.0;

int
testMultiGrid();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
  //  AttachDebugger();
#endif
#if defined(FE_INVALID) && defined(FE_DIVBYZERO) && defined(FE_OVERFLOW)
#ifdef __USE_GNU
  //int except =  FE_INVALID | FE_UNDERFLOW ;
  int except =  FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW ;
  feenableexcept(except);
#endif
#endif

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int overallStatus = 0;
  int status = testMultiGrid();

  if ( status == 0 )
  {
    pout() << indent << pgmname << " passed." << endl ;
  }
  else
  {
    overallStatus = 1;
    pout() << indent << pgmname << " failed with return code " << status << endl ;
  }

  xshift = 0.2;
  blockingFactor = 4;
  CH_TIMER_REPORT();
  CH_TIMER_PRUNE(0.002);
  //CH_TIMER_PRUNE(0.005);
  status = testMultiGrid();

  if ( status == 0 )
  {
    pout() << indent << pgmname << " passed." << endl ;
  }
  else
  {
    overallStatus = 1;
    pout() << indent << pgmname << " failed with return code " << status << endl ;
  }


#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return overallStatus;
}

extern "C"
{
  void Parabola_neum(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    switch (*dir)
    {
    case 0:
      a_values[0]=2*(pos[0]-xshift);
      return;
    case 1:
      a_values[0]=2*pos[1];
    return;
    case 2:
      a_values[0]=2*pos[2];
      return;
    default:
      MayDay::Error("no such dimension");
    };
  }

  void Parabola_diri(Real* pos,
                     int* dir,
                     Side::LoHiSide* side,
                     Real* a_values)
  {
    a_values[0] = D_TERM((pos[0]-xshift)*(pos[0]-xshift),+pos[1]*pos[1],+pos[2]*pos[2]);
  }

  void DirParabolaBC(FArrayBox& a_state,
                     const Box& valid,
                     const ProblemDomain& a_domain,
                     Real a_dx,
                     bool a_homogeneous)
  {

    for (int i=0; i<CH_SPACEDIM; ++i)
      {
        DiriBC(a_state,
               valid,
               dx,
               a_homogeneous,
               Parabola_diri,
               i,
               Side::Lo);
        DiriBC(a_state,
               valid,
               dx,
               a_homogeneous,
               Parabola_diri,
               i,
               Side::Hi);
      }
  }

  void NeumParabolaBC(FArrayBox& a_state,
                      const Box& valid,
                      const ProblemDomain& a_domain,
                      Real a_dx,
                      bool a_homogeneous)
  {

    for (int i=0; i<CH_SPACEDIM; ++i)
      {
        NeumBC(a_state,
               valid,
               dx,
               a_homogeneous,
               Parabola_neum,
               i,
               Side::Lo);
        NeumBC(a_state,
               valid,
               dx,
               a_homogeneous,
               Parabola_neum,
               i,
               Side::Hi);
      }
  }
}

static BCValueFunc pointFunc = Parabola_diri;

void parabola(const Box& box, int comps, FArrayBox& t)
{
  RealVect pos;
  Side::LoHiSide side;
  int dir;
  int num = 1;
  ForAllXBNN(Real,t, box, 0, comps)
    {
      num=nR;
      D_TERM(pos[0]=dx*(iR+0.5);, pos[1]=dx*(jR+0.5);, pos[2]=dx*(kR+0.5));
      pointFunc(&(pos[0]), &dir, &side, &tR);
    }EndFor;
}

void makeGrids(DisjointBoxLayout& a_dbl, const Box& a_domain)
{

  BRMeshRefine br;
  Box domain = a_domain;
  domain.coarsen(blockingFactor);
  domain.refine(blockingFactor);
 CH_assert(domain == a_domain);
  domain.coarsen(blockingFactor);

  ProblemDomain junk(domain);
  IntVectSet pnd(domain);
  IntVectSet tags;
  for (BoxIterator bit(domain); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      if (D_TERM(true, && iv[1]< 2*iv[0] && iv[1]>iv[0]/2, && iv[2] < domain.bigEnd(2)/2))
        {
          tags|=iv;
        }
    }
  Vector<Box> boxes;
  br.makeBoxes(boxes, tags, pnd, junk, 32/blockingFactor, 1);
  Vector<int> procs;
  LoadBalance(procs, boxes);
  for (int i=0; i<boxes.size(); ++i) boxes[i].refine(blockingFactor);
  a_dbl.define(boxes, procs);
}

struct setvalue
{
  static Real val;
  static void setFunc(const Box& box,
                      int comps, FArrayBox& t)
  {
    t.setVal(val);
  }
};

Real setvalue::val = 0;

int
testMultiGrid()
{

  ProblemDomain regularDomain(domain);
  pout()<<"\n single grid MultiGrid solver \n";
  // single grid solver test
  {
    Box phiBox = domain;
    phiBox.grow(1);

    FArrayBox phi(phiBox, 1);
    FArrayBox rhs(domain, 1);
    FArrayBox error(domain, 1);
    FArrayBox phi_exact(domain, 1);
    FArrayBox residual(domain, 1);
    FArrayBox correction(phiBox, 1);

    phi.setVal(0.0);
    rhs.setVal(2*CH_SPACEDIM);
    parabola(domain, 1, phi_exact);

    RealVect pos(IntVect::Unit);
    pos*=dx;

    NewPoissonOpFactory opFactory;

    opFactory.define(pos, DirParabolaBC);
    MultiGrid<FArrayBox> solver;
    BiCGStabSolver<FArrayBox> bottomSolver;
    bottomSolver.m_verbosity=0;
    MGLevelOp<FArrayBox>* op = opFactory.MGnewOp(regularDomain,0);

    solver.define(opFactory, &bottomSolver, regularDomain);

    solver.m_numMG=1;
    int iter = 5;

    pout()<< "homogeneous solver mode : solver.oneCycle(correction, reisdual)\n"
         <<std::endl;

    op->residual(residual, phi, rhs);
    Real rnorm = residual.norm();
    op->axby(error, phi, phi_exact, 1, -1);
    Real norm = error.norm(0);
    pout()<<"initial residual "<<rnorm<<"  Error max norm = "<<norm<<std::endl;

    solver.init(correction, residual);
    for (int i=0; i<iter; ++i)
      {

        correction.setVal(0.0);
        solver.oneCycle(correction, residual);
        op->incr(phi, correction, 1);
        op->residual(residual, phi, rhs);
        rnorm = residual.norm();
        op->axby(error, phi, phi_exact, 1, -1);
        Real norm = error.norm(0);
        pout()<<indent<<"residual L2 norm "<< rnorm
              <<"  Error max norm = "<<norm<<std::endl;
      }

    pout()<< "\n\n inhomogeneous solver mode : solver.oneCycle(a_phi, a_rhs)\n"
          <<std::endl;

    solver.m_homogeneous = false;

    op->scale(phi, 0.0);
    op->residual(residual, phi, rhs);
    rnorm = residual.norm();
    pout()<<"initial residual "<<rnorm<<"\n";
    solver.init(phi, rhs);
    for (int i=0; i<iter; ++i)
      {
        solver.oneCycle(phi, rhs);
        op->residual(residual, phi, rhs);
        rnorm = residual.norm();
        op->axby(error, phi, phi_exact, 1, -1);
        Real norm = error.norm(0);
        pout()<<indent<<"   residual L2 norm "<< rnorm
              <<"  Error max norm = "<<norm<<std::endl;
      }

    delete op;

  }
  pout()<<"\n level solver \n";
  //  Level solve
  {
    DisjointBoxLayout  dbl;

    makeGrids(dbl, domain);

    dbl.close();

    DataIterator dit(dbl);
    LevelData<FArrayBox> phi(dbl, 1, IntVect::Unit);
    LevelData<FArrayBox> correction(dbl, 1, IntVect::Unit);
    LevelData<FArrayBox> phi_exact(dbl, 1);
    LevelData<FArrayBox> error(dbl, 1);
    LevelData<FArrayBox> rhs(dbl, 1);
    LevelData<FArrayBox> residual(dbl, 1);

    setvalue::val = 2*CH_SPACEDIM;
    rhs.apply(setvalue::setFunc);
    setvalue::val = 0;
    phi.apply(setvalue::setFunc);
    phi_exact.apply(parabola);

    RealVect pos(IntVect::Unit);
    pos*=dx;

    AMRPoissonOpFactory opFactory;

    opFactory.define(regularDomain, dbl, pos[0], DirParabolaBC, 1);
    AMRLevelOpFactory<LevelData<FArrayBox> >& castFact  =
      (AMRLevelOpFactory<LevelData<FArrayBox> >&)opFactory;
    MultiGrid<LevelData<FArrayBox> > solver;
    BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;

    MGLevelOp<LevelData<FArrayBox> >* op = castFact.MGnewOp(regularDomain,0);

    solver.define(castFact, &bottomSolver, regularDomain);

    bottomSolver.m_verbosity = 0;
    int iter = 3;
    op->scale(phi, 0.0);
    op->axby(error, phi, phi_exact, 1, -1);
    op->residual(residual, phi, rhs, false);
    Real rnorm = op->norm(residual, 2);
    Real enorm = op->norm(error, 0);

    pout()<< "homogeneous solver mode : solver.oneCycle(correction, residual)\n";

    pout()<<"\nInitial residual norm "<<rnorm<<" Error max norm "<<enorm<<"\n\n";
    solver.init(correction, residual);
    for (int i=0; i<iter; ++i)
      {
        op->scale(correction, 0.0);
        solver.oneCycle(correction, residual);
        op->incr(phi, correction, 1.0);
        op->axby(error, phi, phi_exact, 1, -1);
        op->residual(residual, phi, rhs, false);
        rnorm = op->norm(residual, 2);
        enorm = op->norm(error, 0);
        pout()<<indent<<"residual norm "<<rnorm<<"   Error max norm = "<<enorm<<std::endl;
      }

    delete op;
  }

  return 0;
}
