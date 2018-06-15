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
#include "parstream.H"
#include "BoxIterator.H"
#include "FABView.H"

#include "NewPoissonOp.H"
#include "AMRPoissonOp.H"
#include "BCFunc.H"
#include "BiCGStabSolver.H"
#include "RelaxSolver.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testBiCGStab" ;
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
// #define __USE_GNU
// #include <fenv.h>
// #undef  __USE_GNU

//static Box domain = Box(IntVect(D_DECL(-64,-64,-64)), IntVect(D_DECL(63,63,63)));
static Box domain = Box(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(63,63,63)));
static Real dx = 0.0125;
static int blockingFactor = 8;
static Real xshift = 0.0;

int
testBiCGStab();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  //  int except =  FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW |  FE_INVALID ;
//  feenableexcept(except);

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int overallStatus = 0;
  int status = testBiCGStab();

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
  status = testBiCGStab();

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
      a_values[0]=2.*(pos[0]-xshift);
      return;
    case 1:
      a_values[0]=2.*pos[1];
    return;
    case 2:
      a_values[0]=2.*pos[2];
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
                      const ProblemDomain& a_domain,
                      Real a_dx,
                      bool a_homogeneous)
  {
    Box valid = a_state.box();
    valid.grow(-1);
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
testBiCGStab()
{

  ProblemDomain regularDomain(domain);
  pout()<<"\n GSRB unigrid solver \n";
  // GSRB single grid solver test
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
    correction.setVal(0.0);
    rhs.setVal(2*CH_SPACEDIM);

    parabola(domain, 1, phi_exact);
    RealVect pos(IntVect::Unit);
    pos*=dx;

    NewPoissonOp op;

    op.define(pos, regularDomain, DirParabolaBC);

    for (int i=0; i<15; ++i)
      {
        op.residual(residual, phi, rhs);
        Real rnorm = residual.norm();
        op.preCond(correction, residual);
        op.incr(phi, correction, 1.0);

        op.axby(error, phi, phi_exact, 1, -1);
        Real norm = error.norm(0);
        pout()<<indent<<"Residual L2 norm "<<rnorm<<"  Error max norm = "
              <<norm<<std::endl;
      }
  }

  pout()<<"\n unigrid solver \n";

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

    NewPoissonOp op;

    op.define(pos, regularDomain, DirParabolaBC);
    BiCGStabSolver<FArrayBox> solver;

    solver.define(&op, true);

    int iter = 1;

    pout()<< "homogeneous solver mode : solver.solve(correction, reisdual)\n"
          << "solver.i_max= "<<solver.m_imax<<" iter = "<<iter<<std::endl;

    op.residual(residual, phi, rhs);
    Real rnorm = residual.norm();
    pout()<<"initial residual "<<rnorm<<"\n";
    for (int i=0; i<iter; ++i)
      {

        correction.setVal(0.0);
        solver.solve(correction, residual);
        op.incr(phi, correction, 1);
        op.residual(residual, phi, rhs);
        rnorm = residual.norm();
        op.axby(error, phi, phi_exact, 1, -1);
        Real norm = error.norm(0);
        pout()<<indent<<"residual L2 norm "<< rnorm
              <<"  Error max norm = "<<norm<<std::endl;
      }

    pout()<< "\n\n inhomogeneous solver mode : solver.solve(a_phi, a_rhs)\n"<<std::endl;
    solver.setHomogeneous(false);

    op.scale(phi, 0.0);
    op.residual(residual, phi, rhs);
    rnorm = residual.norm();
    pout()<<"initial residual "<<rnorm<<"\n";
    for (int i=0; i<iter; ++i)
      {
        solver.solve(phi, rhs);
        op.residual(residual, phi, rhs);
        rnorm = residual.norm();
        op.axby(error, phi, phi_exact, 1, -1);
        Real norm = error.norm(0);
        pout()<<indent<<"   residual L2 norm "<< rnorm
              <<"  Error max norm = "<<norm<<std::endl;
      }

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

    AMRPoissonOp amrop;

    amrop.define(dbl, pos[0], regularDomain, DirParabolaBC);

    BiCGStabSolver<LevelData<FArrayBox> > bsolver;
    RelaxSolver<LevelData<FArrayBox> > rsolver;

    bsolver.define(&amrop, true);
    rsolver.define(&amrop, true);

    int iter = 1;
    amrop.scale(phi, 0.0);
    amrop.axby(error, phi, phi_exact, 1, -1);
    amrop.residual(residual, phi, rhs, false);
    Real rnorm = amrop.norm(residual, 2);
    Real enorm = amrop.norm(error, 0);

    pout()<< "homogeneous solver mode : solver.solve(correction, residual)\n"
          << "solver.i_max= "<<bsolver.m_imax<<" iter = "<<iter<<std::endl;
    pout()<<"\nInitial residual norm "<<rnorm<<" Error max norm "<<enorm<<"\n\n";
    pout()<<indent2<<"BiCGStab\n";
    for (int i=0; i<iter; ++i)
      {
        amrop.scale(correction, 0.0);
        bsolver.solve(correction, residual);
        amrop.incr(phi, correction, 1.0);
        amrop.axby(error, phi, phi_exact, 1, -1);
        amrop.residual(residual, phi, rhs, false);
        rnorm = amrop.norm(residual, 2);
        enorm = amrop.norm(error, 0);
        pout()<<indent<<"residual norm "<<rnorm<<"   Error max norm = "<<enorm<<std::endl;
      }
    amrop.scale(phi, 0.0);
    amrop.residual(residual, phi, rhs, false);
    pout()<<indent2<<"RelaxSolver\n";
    for (int i=0; i<iter; ++i)
      {
        amrop.scale(correction, 0.0);
        rsolver.solve(correction, residual);
        amrop.incr(phi, correction, 1.0);
        amrop.axby(error, phi, phi_exact, 1, -1);
        amrop.residual(residual, phi, rhs, false);
        rnorm = amrop.norm(residual, 2);
        enorm = amrop.norm(error, 0);
        pout()<<indent<<"residual norm "<<rnorm<<"   Error max norm = "<<enorm<<std::endl;
      }


    pout()<< "\n\ninhomogeneous solver mode : solver.solve(phi, rhs)\n"
          << "solver.i_max= "<<bsolver.m_imax<<" iter = "<<iter<<std::endl;
    bsolver.setHomogeneous(false);
    rsolver.setHomogeneous(false);

    amrop.scale(phi, 0.0);
    amrop.residual(residual, phi, rhs, false);
    pout()<<indent2<<"BiCGStab\n";
    for (int i=0; i<iter; ++i)
      {
        bsolver.solve(phi, rhs);
        amrop.axby(error, phi, phi_exact, 1, -1);
        amrop.residual(residual, phi, rhs, false);
        rnorm = amrop.norm(residual, 2);
        enorm = amrop.norm(error, 0);
        pout()<<indent<<"residual norm "<<rnorm<<"   Error max norm = "<<enorm<<std::endl;
       }
    amrop.scale(phi, 0.0);
    amrop.residual(residual, phi, rhs, false);
    pout()<<indent2<<"RelaxSolver\n";
    for (int i=0; i<iter; ++i)
      {
        rsolver.solve(phi, rhs);
        amrop.axby(error, phi, phi_exact, 1, -1);
        amrop.residual(residual, phi, rhs, false);
        rnorm = amrop.norm(residual, 2);
        enorm = amrop.norm(error, 0);
        pout()<<indent<<"residual norm "<<rnorm<<"   Error max norm = "<<enorm<<std::endl;
      }

  }

  return 0;
}
