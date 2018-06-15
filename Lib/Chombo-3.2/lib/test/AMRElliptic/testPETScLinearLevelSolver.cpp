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

//#include "NewPoissonOp.H"
//#include "AMRPoissonOp.H"
//#include "BCFunc.H"
//#include "RelaxSolver.H"
#ifdef CH_USE_PETSC
#include "PetscSolver.H"
#endif
#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testPETScLinearLevelSolver" ;
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

int
testPETScLinearLevelSolver();

int
main(int argc ,char* argv[])
{
#ifdef CH_USE_PETSC
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc, &argv,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif 
#endif // end petsc conditional

  //  int except =  FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW |  FE_INVALID ;
//  feenableexcept(except);

  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int overallStatus = 0;
  int status = testPETScLinearLevelSolver();

  if ( status == 0 )
    {
      pout() << indent << pgmname << " passed." << endl ;
    }
  else
    {
      overallStatus = 1;
      pout() << indent << pgmname << " failed with return code " << status << endl ;
    }
  
#ifdef CH_USE_PETSC
  ierr = PetscFinalize(); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Finalize();
#endif // mpi conditional
#endif // petsc conditional
  
  return overallStatus;
}

void makeGrids(DisjointBoxLayout& a_dbl, const Box& a_domain)
{
  int blockingFactor = 8;  
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
      //if (D_TERM(true, && iv[1]< 2*iv[0] && iv[1]>iv[0]/2, && iv[2] < domain.bigEnd(2)/2))
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

const Real kx[] = {  2*M_PI, 3*M_PI, 4*M_PI };
const Real ky[] = {  M_PI,   2*M_PI, 5*M_PI };
const Real kz[] = {  3*M_PI, 4*M_PI, 7*M_PI };
const int nk = sizeof(kx)/sizeof(kx[0]);
 
int
testPETScLinearLevelSolver()
{
  int ncells = 128;
  Box domain = Box(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(ncells-1,ncells-1,ncells-1)));
  Real dx = 1./(Real)ncells;
  ProblemDomain regularDomain(domain); // defaults to non-periodic
  
  DisjointBoxLayout  dbl;
  makeGrids( dbl, domain );
  dbl.close();
  
  pout()<<"\n level solver \n";
  
  DataIterator dit(dbl);
  LevelData<FArrayBox> phi(dbl, 1, IntVect::Unit);
  LevelData<FArrayBox> exact(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox> error(dbl, 1, IntVect::Zero);
  LevelData<FArrayBox> rhs(dbl, 1, IntVect::Zero);
  
  // set RHS
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisRhs = rhs[dit];
      FArrayBox& exFab = exact[dit];
      
      BoxIterator bit(thisRhs.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc(iv);
          loc *= dx;
          loc += 0.5*dx*RealVect::Unit;
          
          for (int i=0; i<nk; i++)
            {
              Real fact = sqrt(D_TERM(kx[i]*kx[i],+ky[i]*ky[i],+kz[i]*kz[i]));
              RealVect x = loc;
              thisRhs(iv,0) = fact*(D_TERM(sin(kx[i]*x[0]),
                                           *sin(ky[i]*x[1]),
                                           *sin(kz[i]*x[2])));
              exFab(iv,0) = (1./fact)*(D_TERM(sin(kx[i]*x[0]),
                                              *sin(ky[i]*x[1]),
                                              *sin(kz[i]*x[2])));
            }
        }
    }

#ifdef CH_USE_PETSC  
  PetscSolverPoisson<LevelData<FArrayBox> > solver;
  solver.define( dx, true );
  solver.m_beta = -1.0; // solving m_alpha u + m_beta del^2 u = f

  pout()<< "homogeneous solver mode -- default parameters\n";
  for (int kk=0;kk<2;kk++)
    {
      // first solve will create the internal KSP object
      solver.solve( phi, rhs );
      
      // get error
      Real enorm_2=0.,phinorm_2=0.;
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& exFab = exact[dit];
          FArrayBox& phifab = phi[dit];
          FArrayBox& errorfab = error[dit];
          
          //errorfab.axby( exFab.box(), phifab, exFab, 1, -1 );
          BoxIterator bit(exFab.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              errorfab(iv,0) = exFab(iv,0) - phifab(iv,0);
              phinorm_2 += exFab(iv,0)*exFab(iv,0);
              enorm_2 += errorfab(iv,0)*errorfab(iv,0);
            }
        }
      pout()<<indent<<"Relative error 2 norm = "<< sqrt(enorm_2/phinorm_2) <<std::endl;
      if (kk==1)break;

      // Hacking into solver to set some PETSc parameters - not recommended
      // PETSC will read parameters from a resource file (.petscrc by default) and then the command line.
      // So the first solve here will use parameters in .petscrc if avaliable and comand line args if provided.
      // this second solver will override these parameters.  It is recommended that solver parameters are set 
      // at run time.
      {
        PC pc; PetscErrorCode ierr;
        KSP ksp = solver.getKSP();
        ierr = KSPSetType( ksp, KSPCG );                            CHKERRQ(ierr);
        ierr = KSPSetComputeSingularValues( ksp, PETSC_TRUE ); CHKERRQ(ierr);
        ierr = KSPGetPC( ksp, &pc );                                   CHKERRQ(ierr);
        ierr = PCSetType( pc, PCGAMG );                               CHKERRQ(ierr); 
        ierr = PCGAMGSetNSmooths(pc, 1);                              CHKERRQ(ierr); 
        ierr = KSPSetTolerances(ksp,1.e-9,PETSC_DEFAULT,PETSC_DEFAULT,100); CHKERRQ(ierr); 
        if (verbose)
          {
            ierr = KSPMonitorSet(ksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr); 
          }
      }

      pout()<< "homogeneous solver mode -- gamg solver parameters\n";
    }
#else  
  pout()<<indent<<"Warning: PETSc test not built with PETSc!" <<std::endl;
  return 0;
#endif

  return 0;
}
