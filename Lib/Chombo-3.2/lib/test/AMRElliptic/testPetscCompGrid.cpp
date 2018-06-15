#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_USE_PETSC

#include <petscksp.h>

#include "AMRPoissonOp.H"
#include "PetscCompGridPois.H"
#include "MultilevelLinearOp.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"

#include "UsingNamespace.H"

static char help[] = "Test PetscCompGridPois to create Jacobean for Laplace(u) with Diri BCs.\n\n";

static PetscInt s_nCells0 = 8;
static PetscReal s_error_thresh = 0.25;
static PetscInt s_blockingfactor = 4;
static PetscInt s_nesting = 4;
static PetscInt s_maxboxsz = 16;
static PetscInt s_refRatio = 2;
static PetscInt s_order = 2;
static PetscBool s_debug = PETSC_FALSE;
static PetscBool s_amrfas = PETSC_FALSE;
static PetscBool s_amrmg = PETSC_FALSE;
static PetscBool s_plot = PETSC_TRUE;
static bool s_amr_type_iserror = false;
static PetscBool s_corner_stencil = PETSC_TRUE;
static PetscBool s_matlab = PETSC_FALSE;

/**
 *  Multigrid Poisson solver:
 *    Laplace(u) = rhs
 *
 *  A numerical test on page 64 of
 *   `A MULTIGRID TUTORIAL by Briggs, Henson & McCormick'
 *   is duplicated here.
 *  
 *  exact u := (x^2-x^4)*(y^4-y^2)*(z^4-z^2)
 *
 */

// RHS
Real rhsFunc(RealVect& a_x)
{
  Real r = 0.0;
  const RealVect x2 = a_x*a_x;
  const RealVect a = x2*(1-x2);
  const RealVect b = 2-12*x2;
  // cross products of coordinates prevents the use of D_TERM
  if (SpaceDim==1)
    r = b[0];
  else if (SpaceDim==2)
    r = b[0]*a[1] + a[0]*b[1];
  else if (SpaceDim==3)
    r = b[0]*a[1]*a[2] + a[0]*b[1]*a[2] + a[0]*a[1]*b[2];
  else
    MayDay::Error("Invalid Dimension!");
  return -r; // lets keep soluion positive
}

// exact u
Real exactSolution(RealVect& a_x)
{
  RealVect e = a_x*a_x;
  e *= 1. - e;
  return e.product(); 
}

#undef __FUNCT__
#define __FUNCT__ "plotAll"
PetscErrorCode plotAll( Vector<LevelData<FArrayBox> *> &a_phi,
                        Vector<LevelData<FArrayBox> *> &a_rhs,
                        Vector<RefCountedPtr<LevelData<FArrayBox> > > &a_exact,
                        Real a_errNorm[2], string a_fname, Real a_cdx,
                        Vector<DisjointBoxLayout> &a_grids,
                        Vector<int> &a_refratios,
                        Vector<ProblemDomain> &a_domains,
                        PetscCompGridPois &a_petscop,
                        Vec a_x,
                        int a_sub_id = -1 )
{
  CH_TIME("plotAll");
  int nLev = a_phi.size();
  PetscErrorCode ierr;
  Vector<LevelData<FArrayBox>* > plotData(nLev, NULL);
  
  if ( a_x )
    {
      ierr = a_petscop.putPetscInChombo(a_x,a_phi); CHKERRQ(ierr);
    }

  for (int ilev=0;ilev<nLev;ilev++) 
    {      
      plotData[ilev] = new LevelData<FArrayBox>(a_grids[ilev],4*COMP_POIS_DOF,IntVect::Unit);
    }

  a_errNorm[0] = a_errNorm[1] = 0;
  Real dx = a_cdx;
  for (int ilev=0;ilev<nLev;ilev++,dx/=s_refRatio) 
    {
      Interval phiInterval(0,COMP_POIS_DOF-1);
      a_phi[ilev]->copyTo(phiInterval, *plotData[ilev], phiInterval);
      Interval rhsInterval(COMP_POIS_DOF,2*COMP_POIS_DOF-1);
      a_rhs[ilev]->copyTo(phiInterval, *plotData[ilev], rhsInterval);
      Interval exInterval(2*COMP_POIS_DOF,3*COMP_POIS_DOF-1);
      a_exact[ilev]->copyTo(phiInterval, *plotData[ilev], exInterval);
      // use phi for error
      const DisjointBoxLayout& dbl = a_grids[ilev];
      for (DataIterator dit(dbl); dit.ok(); ++dit)
        {
          FArrayBox& exactfab = (*a_exact[ilev])[dit];
          FArrayBox& phiFAB = (*a_phi[ilev])[dit];
          Box region = exactfab.box();
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              for (int i=0;i<COMP_POIS_DOF;i++)
                phiFAB(iv,i) = phiFAB(iv,i) - exactfab(iv,i);
            }
        }
      
      // zero error on covered
      if (ilev!=nLev-1) 
        {
          const DisjointBoxLayout& dbl = a_grids[ilev];
          // zero out fine cover
          DisjointBoxLayout dblCoarsenedFine;
          Copier copier;
          coarsen(dblCoarsenedFine, a_grids[ilev+1], a_refratios[ilev]); // coarsens entire grid
          copier.define(dblCoarsenedFine, dbl, IntVect::Zero);
          LevelDataOps<FArrayBox> ops;
          ops.copyToZero(*a_phi[ilev],copier);
        }

      // copy in
      Interval errInterval(3*COMP_POIS_DOF,4*COMP_POIS_DOF-1);
      a_phi[ilev]->copyTo(phiInterval, *plotData[ilev], errInterval);

      // get error norms
      for (DataIterator dit(dbl); dit.ok(); ++dit)
        {
          Box region = dbl[dit];
          FArrayBox& phifab = (*a_phi[ilev])[dit];
          Real mnorm = phifab.norm(region,0);
          if (mnorm>a_errNorm[0]) a_errNorm[0] = mnorm;
          mnorm = phifab.norm(region,1)*D_TERM(dx,*dx,*dx);
          a_errNorm[1] += mnorm;
        }
    }
  {
    double error;
    MPI_Allreduce( &a_errNorm[0], &error, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD );
    a_errNorm[0] = error;
    MPI_Allreduce( &a_errNorm[1], &error, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD );
    a_errNorm[1] = error;
  }

  pout() << "\t\t plot |error|_inf=" << a_errNorm[0] << endl;
  
  // plot
  if (true){  
    CH_TIME("plot");
    char suffix[30];
    if (a_sub_id>=0) sprintf(suffix, "%dd.%d.hdf5",SpaceDim,a_sub_id);
    else sprintf(suffix, "%dd.hdf5",SpaceDim);
    a_fname += suffix;
    Vector<string> varNames(4*COMP_POIS_DOF);
    int kk=0;
    for (int i=0; i<COMP_POIS_DOF; ++i,kk++) varNames[kk] = "phi ";
    for (int i=0; i<COMP_POIS_DOF; ++i,kk++) varNames[kk] = "rhs ";
    for (int i=0; i<COMP_POIS_DOF; ++i,kk++) varNames[kk] = "exa ";
    for (int i=0; i<COMP_POIS_DOF; ++i,kk++) varNames[kk] = "err ";
    kk=0;
    for (int i=0; i<COMP_POIS_DOF; ++i,kk++) varNames[kk][3] = '1' + i;
    for (int i=0; i<COMP_POIS_DOF; ++i,kk++) varNames[kk][3] = '1' + i;
    for (int i=0; i<COMP_POIS_DOF; ++i,kk++) varNames[kk][3] = '1' + i;
    for (int i=0; i<COMP_POIS_DOF; ++i,kk++) varNames[kk][3] = '1' + i;

    Real bogusVal = 1.0;
    WriteAMRHierarchyHDF5(a_fname,
                          a_grids,
                          plotData,
                          varNames,
                          a_domains[0].domainBox(),
                          a_cdx,
                          bogusVal,
                          bogusVal,
                          a_refratios,
                          nLev);
  }

  for (int ilev=0;ilev<nLev;ilev++) 
    {
      delete plotData[ilev];
    }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "setupGrids"
PetscErrorCode setupGrids( Vector<int> &a_refratios,  // out
                           Vector<ProblemDomain> &a_domains, // out
                           Vector<DisjointBoxLayout> &a_grids, // in/out
                           int a_nLevs,
                           Real &a_cdx, // out
                           Vector<RefCountedPtr<LevelData<FArrayBox> > > a_exact,
                           Vector<RefCountedPtr<LevelData<FArrayBox> > > a_phi,
                           Real a_max_norm_error,
                           Vector<IntVectSet> &a_tagVects,
                           bool &a_same // out
                           )
{
  CH_TIME("setupGrids");
  const int new_level=a_nLevs-1;
  bool isperiodic[3] = {false,false,false};
  PetscFunctionBeginUser;

  a_cdx = 1./s_nCells0; // out

  a_refratios.resize(a_nLevs-1);
  a_domains.resize(a_nLevs);
  a_grids.resize(a_nLevs);

  if (a_nLevs==1)
    {
      Vector<int> procAssign;
      Vector<Box> boxvector;
      Box levbox(IntVect::Zero,(s_nCells0-1)*IntVect::Unit), dombox=levbox;
      pout() << "setupGrids for level " << new_level+1 << "/" << a_nLevs << ". level domain " << levbox << endl;
      domainSplit(levbox, boxvector, s_maxboxsz, s_blockingfactor); // need blocking factor of 4 for C-F
      LoadBalance( procAssign, boxvector );
      a_grids[new_level].define(boxvector,procAssign);
      a_domains[new_level].define(dombox,isperiodic);
      a_same = false;
    }
  else
    {
      Vector<Vector<Box> > old_grids(a_nLevs-1);
      int new_finest_level;
      Real fillRatio = .8;
      BRMeshRefine refiner;
      Box curr_dom_box;
      static double thresh;

      // refratio
      a_refratios[new_level-1] = s_refRatio;
      // make new domain
      {
        Box dombox = a_domains[new_level-1].domainBox();
        dombox.refine(s_refRatio);
        a_domains[new_level].define(dombox,isperiodic);
      }

      refiner.define(a_domains[0],a_refratios,fillRatio,s_blockingfactor,s_nesting,s_maxboxsz);

      Real dx = a_cdx;
      if (a_nLevs==2)
        {
          // get max_grad on coarsest level - has the whole domain
          double max_ = 0.;
          const DisjointBoxLayout& dbl = a_grids[0];
          for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
            {
              FArrayBox& exactFAB = (*a_exact[0])[dit];
              FArrayBox& phiFAB = (*a_phi[0])[dit];
              Box region = dbl[dit];
              for (BoxIterator bit(region); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  Real grad = 0., tt;
                  if (s_amr_type_iserror)
                    {
                      tt = abs(phiFAB(iv,0)-exactFAB(iv,0));
                      if (tt>max_) max_ = tt;
                    }
                  else
                    {
                      for (int dir=0; dir<SpaceDim; ++dir)
                        {
                          IntVect liv(iv), riv(iv); liv.shift(dir,-1); riv.shift(dir,1); 
                          tt = (exactFAB(riv,0)-exactFAB(liv,0))/(2.*dx);
                          grad += tt*tt;
                        }
                      grad = sqrt(grad);
                      if (grad>max_) max_ = grad;
                    }
                }
            }
          
          MPI_Allreduce( &max_, &thresh, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD );
          // set thresh
          thresh /= 2.; 
        }
      
      // make tags
      curr_dom_box = a_domains[0].domainBox(); // used for hard wired refinement
      for (int ilev=0;ilev<new_level;ilev++)
        {
          const DisjointBoxLayout& dbl = a_grids[ilev];
          for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
            {
              FArrayBox& exactFAB = (*a_exact[ilev])[dit];
              FArrayBox& phiFAB = (*a_phi[ilev])[dit];
              Box region = dbl[dit];
              for (BoxIterator bit(region); bit.ok(); ++bit)
                {
                  const IntVect& iv = bit();
                  if (s_amr_type_iserror)
                    {
                      Real error = abs(phiFAB(iv,0)-exactFAB(iv,0));
                      if (error > thresh) a_tagVects[ilev] |= iv; // real AMR
                    }
                  else // grad
                    {
                      Real grad = 0., tt;
                      for (int dir=0; dir<SpaceDim; ++dir)
                        {
                          IntVect liv(iv), riv(iv); liv.shift(dir,-1); riv.shift(dir,1); 
                          tt = (exactFAB(riv,0)-exactFAB(liv,0))/(2.*dx);
                          grad += tt*tt;
                        }
                      grad = sqrt(grad);
                      if (grad > thresh) a_tagVects[ilev] |= iv; // real AMR
                    }
                }
            } // grid

          // hardwire refinement in fake place
          IntVect se = a_domains[ilev].domainBox().smallEnd();
          //se.shift(1,curr_dom_box.size(1)-1);
          a_tagVects[ilev] |= se; 

          curr_dom_box.refine(s_refRatio);
          dx /= s_refRatio;
          thresh *= s_error_thresh; 
        }
  
      // make old_grids, need to make all because coarse grids can change after refinement
      for (int ilev=0;ilev<new_level;ilev++)
        {
          const DisjointBoxLayout& dbl = a_grids[ilev];
          old_grids[ilev].resize(dbl.size());
          LayoutIterator lit = dbl.layoutIterator();
          int boxIndex = 0;
          for (lit.begin(); lit.ok(); ++lit, ++boxIndex) 
            {
              old_grids[ilev][boxIndex] = dbl[lit()];
            }
        }
      
      pout() << "setupGrids for level " << new_level+1 << "/" << a_nLevs << ". level domain " << a_domains[new_level].domainBox() << endl;

      // want to keep these
      Vector<IntVectSet> tagVects_copy(a_nLevs-1);
      for (int ilev=0;ilev<new_level;ilev++)
        {
          tagVects_copy[ilev] = a_tagVects[ilev];
        }
      Vector<Vector<Box> > new_grids;
      new_finest_level = refiner.regrid( new_grids, tagVects_copy, 
                                         0, new_level-1, 
                                         old_grids );

      if (new_finest_level == new_level-1) 
        {
          a_same = true;
        }
      else
        {
          // test to see if grids have changed, create new DBLs
          a_same = true;
          for (int ilev=1; ilev<=new_finest_level; ++ilev)
            {
              int numGridsNew = new_grids[ilev].size();
              Vector<int> procIDs(numGridsNew);
              LoadBalance(procIDs, new_grids[ilev]);
              const DisjointBoxLayout newDBL( new_grids[ilev], procIDs,
                                              a_domains[ilev] );
              const DisjointBoxLayout oldDBL = a_grids[ilev];
              a_same &= oldDBL.sameBoxes(newDBL);
              a_grids[ilev] = newDBL;
            }
        }
      if (a_same)
        {
          pout() << "setupGrids -- grids unchanged" << endl;
        }
    }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "go"
PetscErrorCode go(int nGrids, int &status)
{
  CH_TIME("go");
  PetscErrorCode ierr;
  Real errNorm[20][2];
  Real convergeRate[20-1][2];
  bool sameGrids;
  const Real log2r = 1.0/log(2.0);
  const Real targetConvergeRate = 2*0.9;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > phi;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > rhs;
  Vector<LevelData<FArrayBox> *> phi2;
  Vector<LevelData<FArrayBox> *> rhs2;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > exact;
  Vector<int> refratios;
  Vector<ProblemDomain> domains;
  Vector<DisjointBoxLayout> grids;
  Vector<IntVectSet> tagVects(nGrids-1);
  Real cdx,dx;
  PetscCompGridPois petscop(0.,-1.,s_order);
  RefCountedPtr<ConstDiriBC> bcfunc = RefCountedPtr<ConstDiriBC>(new ConstDiriBC(1,petscop.getGhostVect()));
  BCHolder bc(bcfunc);
  PetscFunctionBeginUser;

  petscop.setCornerStencil(s_corner_stencil);
  petscop.setMatlab(s_matlab);

  for (int nLev=1,iMaxLev=0; nLev<=nGrids; nLev++,iMaxLev++)
    {
      ierr = setupGrids( refratios, domains, grids, nLev, cdx, 
                         exact, phi,
                         iMaxLev==0 ? 0.: errNorm[iMaxLev-1][0],
                         tagVects, sameGrids ); CHKERRQ(ierr);

      if (sameGrids){ nGrids = iMaxLev; break;}

      // allocate vectors, set RHS and exact
      phi.resize(nLev); rhs.resize(nLev); exact.resize(nLev);
      phi2.resize(nLev); rhs2.resize(nLev); 
      dx = cdx; 
      for (int ilev=0;ilev<nLev;ilev++) 
        {
          phi[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[ilev],COMP_POIS_DOF,petscop.getGhostVect()));
          rhs[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[ilev],COMP_POIS_DOF,IntVect::Zero));
          exact[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids[ilev],COMP_POIS_DOF,IntVect::Unit));
          rhs2[ilev] = &(*rhs[ilev]);
          phi2[ilev] = &(*phi[ilev]);

          const DisjointBoxLayout &dbl = grids[ilev];
          for (DataIterator dit(dbl) ; dit.ok(); ++dit)
            {
              FArrayBox& rhsFAB = (*rhs[ilev])[dit];
              Box region = dbl[dit];
              for (BoxIterator bit(region); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv); loc *= dx; loc += 0.5*dx*RealVect::Unit;
                  for (int i=0;i<COMP_POIS_DOF;i++)
                    rhsFAB(iv,i) = rhsFunc(loc);
                }

              FArrayBox& exactfab = (*exact[ilev])[dit];
              FArrayBox& phifab = (*phi[ilev])[dit];
              Box eregion = exactfab.box();
              for (BoxIterator bit(eregion); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect loc(iv); loc *= dx; loc += 0.5*dx*RealVect::Unit;
                  for (int i=0;i<COMP_POIS_DOF;i++)
                    exactfab(iv,i) = exactSolution(loc);
                }
              phifab *= .0;
            }
          if (ilev!=nLev-1) dx /= refratios[ilev];
        }

      // form matrix
      petscop.define(domains[0],grids,refratios,bc,cdx*RealVect::Unit);
      petscop.setVerbose(1);

      ierr = petscop.createMatrix(); CHKERRQ(ierr);

      // solve
      Vec  x, b;      /* approx solution, RHS */
      {
        CH_TIME("PETSC-solve");
        Mat  A;         /* linear system matrix */
        KSP  ksp;       /* linear solver context */
        A = petscop.getMatrix();
        ierr = MatGetVecs(A,&x,&b); CHKERRQ(ierr);
        ierr = petscop.putChomboInPetsc(rhs2,b); CHKERRQ(ierr);
        ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
        ierr = KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
        ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
        ierr = KSPDestroy(&ksp); CHKERRQ(ierr); 
      }

      if (s_plot)
        {
          ierr = plotAll(phi2,rhs2,exact,errNorm[iMaxLev],"testPetscCompGrid.",
                         cdx,grids,refratios,domains,petscop,x,nLev); CHKERRQ(ierr);
        }
      ierr = VecDestroy(&x); CHKERRQ(ierr);
      ierr = VecDestroy(&b); CHKERRQ(ierr);
    } // 
#ifdef CH_USE_FAS
  // do AMRFAS
  if (s_amrfas)
    {
    int status;
    AMRFAS<LevelData<FArrayBox> > amrFASSolver;
    FASPoissonOpFactory opFactory(2);

    // solving poisson problem here
    opFactory.define( bc, 0., -1. );
    amrFASSolver.m_pre = 4;
    amrFASSolver.m_post = 4;
    amrFASSolver.m_verbosity = 5;
    s_fas_verbosity = 1; // this is in AMRFASI.H
    amrFASSolver.setAvoidNorms(false);
    amrFASSolver.setInternalCoarseningRate(2);
    amrFASSolver.setNumVcycles(1); 
    amrFASSolver.setCycleType(FAS_FULL);
    amrFASSolver.m_max_iter = 100;
    amrFASSolver.m_atol = 1.0e-9;
    amrFASSolver.m_rtol = 1.0e-9;
    amrFASSolver.m_stagnate = 1.e-2;
    amrFASSolver.setSmootherType(FAS_RICH);
    amrFASSolver.setSmoothingDampingFactor(CH_SPACEDIM==3 ? .7 : .8);
    amrFASSolver.setFMGProlOrderP(2);
    amrFASSolver.setProlOrderP(1);

    amrFASSolver.define( domains[0], cdx*RealVect::Unit, 
                         grids, refratios, opFactory );

    for (int ilev=0;ilev<nGrids;ilev++) 
      {
        const DisjointBoxLayout &dbl = grids[ilev];
        for (DataIterator dit(dbl) ; dit.ok(); ++dit)
          {
            FArrayBox& fab = (*phi[ilev])[dit];
            fab.setVal(0.0,dbl[dit],0,fab.nComp());
          }
      }
    
    amrFASSolver.solve( phi, rhs, &status );

    PetscCompGridPois petscop(0.,-1.,2);
    petscop.define(domains[0],grids,refratios,bc,cdx*RealVect::Unit);
    Real errNorms[2];
    ierr = plotAll( phi2,rhs2,exact,errNorms,"testPetscCompGridAMRFAS.",
                    cdx,grids,refratios,domains,petscop,PETSC_NULL,nGrids); CHKERRQ(ierr);
  }
#endif
   // do AMRMultigrid
  if (s_amrmg)
  {
    AMRMultiGrid<LevelData<FArrayBox> > amrMG;
    AMRPoissonOpFactory opFactory;
    BiCGStabSolver<LevelData<FArrayBox> > biCGStab;
    PetscCompGridPois petscop(0.,-1.,2);

    petscop.define(domains[0],grids,refratios,bc,cdx*RealVect::Unit);

    // solving poisson problem here
    opFactory.define( domains[0], grids, refratios, cdx, bc, 0., -1. );
    
    amrMG.define( domains[0],
                  opFactory,
                  &biCGStab,
                  nGrids);

    amrMG.setSolverParameters(4,4,100,1,100,1.0e-9,1.0e-9,1.0e-9);
    amrMG.m_verbosity = 3;

    amrMG.solve( phi2, rhs2, nGrids-1, 0, true, true );
    Real errNorms[2];
    ierr = plotAll( phi2,rhs2,exact,errNorms,"testPetscCompGridGMG.",
                    cdx,grids,refratios,domains,petscop,PETSC_NULL,nGrids); CHKERRQ(ierr);
  }
  
  pout() << "\nConvergence Rates :  Inf norm     1 norm\n";
  pout() <<   "----------------------------------------\n" << std::setprecision(3);
  if (s_plot)
    {
      for (int i=0; i<nGrids-1; i++)
        {
          pout() << "            ";
          for (int j=0;j<2;j++)
            {
              const Real ratio = errNorm[i][j]/errNorm[i+1][j];
              convergeRate[i][j] = log(ratio)*log2r;
              if (j==1)convergeRate[i][j] *= 2; // this is not quite right
              if (convergeRate[i][j] < targetConvergeRate && i>0 ) // first test misses because of griding
                {
                  status += 1;
                }
              pout() << "         " << convergeRate[i][j];
            }
          pout() << endl;
        }
      if (status==0)
        {
          pout() <<  "All tests passed!\n";
        }
      else
        {
          pout() << status << " tests failed!\n";
        }
    }

  PetscFunctionReturn(0);
}
#endif

int main(int argc, char* argv[])
{ 
  int status = 0;
#ifdef CH_USE_PETSC
  PetscInt nlev;
  PetscErrorCode ierr;
  char string[64]; 
  PetscBool set;
  PetscInitialize(&argc,&argv,(char*)0,help);

  PetscOptionsGetInt(PETSC_NULL,"-n",&s_nCells0,PETSC_NULL);
  PetscOptionsGetBool(PETSC_NULL,"-bc_debug",&s_debug,PETSC_NULL);
  PetscOptionsGetBool(PETSC_NULL,"-plot_fas",&s_amrfas,PETSC_NULL);
  PetscOptionsGetBool(PETSC_NULL,"-plot_mg",&s_amrmg,PETSC_NULL);
  PetscOptionsGetBool(PETSC_NULL,"-corner_stencil",&s_corner_stencil,PETSC_NULL);
  PetscOptionsGetBool(PETSC_NULL,"-matlab",&s_matlab,PETSC_NULL);
  PetscOptionsGetBool(PETSC_NULL,"-plot",&s_plot,PETSC_NULL);
  PetscOptionsGetReal(PETSC_NULL,"-error_thresh",&s_error_thresh,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL,"-blocking_factor",&s_blockingfactor,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL,"-order",&s_order,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL,"-nesting",&s_nesting,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL,"-max_box_size",&s_maxboxsz,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL,"-refinement_ratio",&s_refRatio,PETSC_NULL);
  PetscOptionsGetString(PETSC_NULL,"-amr_type",string,64,&set);
  if (set && strcmp(string,"error")==0) s_amr_type_iserror = true;
  else s_amr_type_iserror = false;

  nlev = 4;
  PetscOptionsGetInt(PETSC_NULL,"-nlevels",&nlev,PETSC_NULL);
  ierr = go(nlev,status); CHKERRQ(ierr);

  CH_TIMER_REPORT();

  ierr = PetscFinalize(); CHKERRQ(ierr);
#endif
  return status;
}
