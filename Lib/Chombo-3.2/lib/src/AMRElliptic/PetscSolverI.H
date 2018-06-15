#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCSOLVERI_H_
#define _PETSCSOLVERI_H_

#include "LevelData.H"
#include "FluxBox.H"
#include "BoxIterator.H"
#include "CH_Timer.H"
#include "memusage.H"
#include "NamespaceHeader.H"

// *******************************************************
// PetscSolver Implementation
// *******************************************************
template <class T>
PetscSolver<T>::PetscSolver()
  :m_homogeneous(false),
   m_mat(0), // m_xx, m_rr, m_bb;
   m_snes(0),
   m_ksp(0),
   m_defined(0),
   m_function(0),
   m_jacobian(0),
   m_null(false),
   m_nz_init_guess(false),
   m_gid0(0)
{
  m_dx = 1.;
  m_prestring[0] = '\0';
}
// *******************************************************
template <class T>
void PetscSolver<T>::destroy()
{
  if ( m_defined )
  {
#ifdef CH_USE_PETSC
    if (m_mat)
      {
        MatDestroy(&m_mat);
        m_mat = 0;
      }
    VecDestroy(&m_bb);
    VecDestroy(&m_xx);
    VecDestroy(&m_rr);
#endif
    m_defined = 0;
  }
#ifdef CH_USE_PETSC
  if ( m_ksp )
    {
      KSPDestroy( &m_ksp );
      m_ksp = 0;
    }
  if ( m_snes )
    {
      SNESDestroy( &m_snes );
      m_snes = 0;
    }
#endif
}
// *******************************************************
template <class T>
PetscSolver<T>::~PetscSolver()
{
  destroy();
}

// *******************************************************
template <class T>
void PetscSolver<T>::define( Real a_dx,
                             bool a_homogeneous )
{
  m_homogeneous = a_homogeneous; // not used!!!
  m_dx = a_dx;
  CH_assert(m_dx!=0.);
}

// *******************************************************
template <class T>
void PetscSolver<T>::define( LinearOp<T> *a_op,
                             bool a_homogeneous )
{
  define( a_op->dx(), a_homogeneous);
}

// *******************************************************
template <class T>
void PetscSolver<T>::setNull( bool n /*= true*/ )
{
  m_null = n; CH_assert(0);
}

// *******************************************************
template <class T>
void PetscSolver<T>::solve( T & a_phi,
                            const T & a_rhs )
{
  T& phi = a_phi;
  const T& rhs = a_rhs;
  solve_private( phi, rhs );
}

// *******************************************************
template <class T>
void PetscSolver<T>::solve_mfree( T & a_phi,
                                  const T & a_rhs, 
                                  LinearOp<T> *a_op )
{
  T& phi = a_phi;
  const T& rhs = a_rhs;
  LinearOp<T>* op = a_op;
  solve_mfree_private( phi, rhs , op);
}

// *******************************************************
//   create_mat_vec
//     Create 'm_mat', 'm_xx', ....  Constructs 'm_gids'.
//
template <class T>
int PetscSolver<T>::create_mat_vec( const T& a_phi )
{
  CH_TIME("PetscSolver::create_mat_vec");
#ifdef CH_USE_PETSC
  const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
  DataIterator dit( dbl );
  PetscErrorCode ierr;
  const PetscInt nc = a_phi.nComp();
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif

  if ( !m_mat )
    {
      // print_memory_line("Before AMG set up");
      m_defined = 2;
 
      IntVect idghosts = a_phi.ghostVect();
      if (idghosts == IntVect::Zero)
        {
          MayDay::Error("PetscSolver<T>::create_mat_vec: No ghost cells in input LevelData<>.");
        }
      m_bccode.define(dbl, 1, idghosts);

      // global ids with ghost cells
      m_gids.define(dbl, 1, idghosts);

      // get first (zero based) id on this processor      
      PetscInt data = 0;
      for ( dit = a_phi.dataIterator() ; dit.ok() ; ++dit )
        {
          const Box &box = dbl.get(dit()); // not ghosted
          data += box.numPts();
          BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
          gidsFab.setVal(-1); // flag for BC eventually
        }
      const PetscInt NN = nc*data;
#ifdef CH_MPI
      PetscInt result;
      MPI_Datatype mtype;
      PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
      MPI_Scan( &data, &result, 1, mtype, MPI_SUM, wcomm );
      m_gid0 = result - data;
#else
      m_gid0 = 0;
#endif
      PetscInt gid = m_gid0;
      for ( dit = a_phi.dataIterator() ; dit.ok() ; ++dit )
        {
          BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
          const Box& box = dbl.get(dit());
          BoxIterator bit(box);
          for (bit.begin(); bit.ok(); bit.next(), gid++ )
            {
              IntVect iv = bit();
              gidsFab(iv,0) = gid;
            }
        }
      m_gids.exchange();

      // create matrix
      PetscInt nnzrow = getNNZPerRow();
      PetscInt *d_nnz=PETSC_NULL, *o_nnz=PETSC_NULL;
      
      if ( supportNNZExact() )
        { 
          PetscInt nglobal, nn = NN;
#ifdef CH_MPI
          MPI_Comm wcomm = Chombo_MPI::comm;
          MPI_Datatype mtype;
          PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
          MPI_Allreduce(&nn,&nglobal,1,mtype,MPI_SUM,wcomm);
#else
          nglobal = NN;
#endif    
          ierr = PetscMalloc( (NN+1)*sizeof(PetscInt), &d_nnz ); CHKERRQ(ierr);
          ierr = PetscMalloc( (NN+1)*sizeof(PetscInt), &o_nnz ); CHKERRQ(ierr);
          for (PetscInt kk=0;kk<NN;kk++) d_nnz[kk] = o_nnz[kk] = 0;

          ierr = formMatrix( m_mat, &a_phi, nc*m_gid0, NN, d_nnz, o_nnz  ); CHKERRQ(ierr);
          CH_assert(!m_mat);
          // fix bounds
          for (PetscInt kk=0;kk<NN;kk++) 
            {
              if (d_nnz[kk] > NN) d_nnz[kk] = NN;
              if (o_nnz[kk] > nglobal-NN) o_nnz[kk] = nglobal-NN;
            }
          nnzrow = 0;
        }

      ierr = MatCreate(wcomm,&m_mat);CHKERRQ(ierr);
      ierr = MatSetSizes(m_mat,NN,NN,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
      ierr = MatSetBlockSize(m_mat,nc);CHKERRQ(ierr);
      ierr = MatSetType(m_mat,MATAIJ);CHKERRQ(ierr);
      ierr = MatSeqAIJSetPreallocation(m_mat,nnzrow, d_nnz);CHKERRQ(ierr);
      ierr = MatMPIAIJSetPreallocation(m_mat,nnzrow, d_nnz, nnzrow/2, o_nnz);CHKERRQ(ierr);
      ierr = MatSetOption(m_mat,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) ;CHKERRQ(ierr);
      ierr = MatSetFromOptions( m_mat ); CHKERRQ(ierr);
      
      if ( d_nnz )
        {
          ierr = PetscFree( d_nnz );  CHKERRQ(ierr);
          ierr = PetscFree( o_nnz );  CHKERRQ(ierr);
        }

      // create vectors
      ierr = VecCreate( wcomm, &m_bb ); CHKERRQ(ierr);
      ierr = VecSetFromOptions( m_bb ); CHKERRQ(ierr);
      ierr = VecSetSizes( m_bb, NN, PETSC_DECIDE ); CHKERRQ(ierr);
      ierr = VecDuplicate( m_bb, &m_rr ); CHKERRQ(ierr);
      ierr = VecDuplicate( m_bb, &m_xx ); CHKERRQ(ierr);
    }
#endif
  return 0;
}

// *******************************************************
//    setup_solver
//      - creates solver if needed.  forms matrix, sets up KSP
//
template <class T>
int PetscSolver<T>::setup_solver( const T& a_phi )
{
  CH_TIMERS("PetscSolver::setup_solver");
  CH_TIMER("solve-setup-1st", t1);
  CH_TIMER("solve-setup-rest", t2);
#ifdef CH_USE_PETSC
  const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
  DataIterator dit( dbl );
  PetscErrorCode ierr;
  KSP ksp;
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif

  if ( m_defined == 0 )
    {
      m_defined = 2;
      // print_memory_line("Before AMG set up");

      ierr = create_mat_vec( a_phi ); CHKERRQ(ierr);

      CH_START(t1);

      // Add values to A
      ierr = formMatrix( m_mat, &a_phi ); CHKERRQ(ierr);

      // create solvers
      PetscBool ism = PETSC_FALSE;
      PetscOptionsGetBool(m_prestring,"-ksp_monitor",&ism,PETSC_NULL);
      if ( m_function && !m_snes )
        {
          ierr = SNESCreate( wcomm, &m_snes );                                       CHKERRQ(ierr);
          // these casts get rid of the 'const'. These are used in SNES as temp vector!  
          ierr = SNESSetFunction( m_snes, m_rr, m_function, (void*)this );         CHKERRQ(ierr);
	  ierr = SNESSetJacobian( m_snes, m_mat, m_mat, m_jacobian, (void*)this ); CHKERRQ(ierr);
          //ierr = SNESSetApplicationContext( m_snes, (void*)this );                   CHKERRQ(ierr);
          ierr = SNESSetFromOptions( m_snes );                                       CHKERRQ(ierr);
          ierr = SNESGetKSP( m_snes, &ksp );                                         CHKERRQ(ierr);
          if (ism)
            {
              ierr = KSPMonitorSet(ksp,ksp_monitor_pout,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
            }
        }
      else if ( !m_function && !m_ksp )
        {
          // create the KSP so that we can set KSP parameters
          KSPCreate( wcomm, &m_ksp );
          if ( strlen(m_prestring) > 0 )
            {
              ierr = KSPSetOptionsPrefix( m_ksp, m_prestring );    CHKERRQ(ierr);
            }
          ierr = KSPSetFromOptions(m_ksp);CHKERRQ(ierr);
          if (ism)
            {
              ierr = KSPMonitorSet(m_ksp,ksp_monitor_pout,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
            }
#if PETSC_VERSION_GE(3,5,0)
	  ierr = KSPSetOperators(m_ksp,m_mat,m_mat);CHKERRQ(ierr);
#else
          ierr = KSPSetOperators(m_ksp,m_mat,m_mat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
#endif
          ierr = KSPSetInitialGuessNonzero(m_ksp, m_nz_init_guess ? PETSC_TRUE : PETSC_FALSE );CHKERRQ(ierr);
          ksp = m_ksp;
        }
      else CH_assert(0);

      { // coordinates
        PC pc; 
        PetscInt sz,ind,bs,n,m;
#if PETSC_VERSION_LT(3,4,0) & PETSC_VERSION_RELEASE
        const PCType type;
#else
        PCType type;
#endif  
        ierr = KSPGetPC( ksp, &pc );     CHKERRQ(ierr);
        ierr = PCGetType( pc, &type );    CHKERRQ(ierr);
        ierr = MatGetBlockSize( m_mat, &bs );               CHKERRQ( ierr );
        if ( strcmp(type,PCGAMG) == 0 && bs > 1 )
          {
            PetscReal    *coords;
            DataIterator dit(a_phi.disjointBoxLayout());
        
            ierr = MatGetLocalSize( m_mat, &m, &n );  CHKERRQ(ierr);
            sz = CH_SPACEDIM*(m/bs);
            ierr = PetscMalloc( (sz+1)*sizeof(PetscReal), &coords ); CHKERRQ(ierr);
            for ( dit = a_phi.dataIterator(), ind = 0 ; dit.ok() ; ++dit )
              {
                const Box &box = a_phi.getBoxes()[dit];
                BoxIterator bit(box);
                for (bit.begin(); bit.ok(); bit.next())
                  {
                    IntVect iv = bit(); // coordinate in any scaled, shifted, rotated frame.
                    for (PetscInt k=0; k<CH_SPACEDIM; k++) coords[ind++] = (PetscReal)iv[k];
                  }
              }
            CH_assert(ind==sz);
            ierr = PCSetCoordinates( pc, CH_SPACEDIM, sz/CH_SPACEDIM, coords ); CHKERRQ(ierr);
            ierr = PetscFree( coords );  CHKERRQ(ierr);
          }
      }

      CH_STOP(t1);
      // print_memory_line("After AMG set up");
    }
  else if ( m_defined == 1 )
    {
      m_defined = 2;
      // form A -- m_mat
      CH_START(t2);
      ierr = MatZeroEntries( m_mat );        CHKERRQ(ierr);
      ierr = formMatrix( m_mat, &a_phi );     CHKERRQ(ierr);
      if ( m_ksp )
        {
          ksp = m_ksp;
        }
      else
        {
          ierr = SNESGetKSP( m_snes, &ksp ); CHKERRQ(ierr);
        }
#if PETSC_VERSION_GE(3,5,0)
      ierr = KSPSetOperators(ksp,m_mat,m_mat); CHKERRQ(ierr);
#else
      ierr = KSPSetOperators(ksp,m_mat,m_mat,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
#endif
      CH_STOP(t2);
    }
#endif
  return 0;
}

// *******************************************************
template <class T>
int PetscSolver<T>::solve_private( T& a_phi,
                                   const T& a_rhs )
{
  CH_TIMERS("PetscSolver::solve_private");
  CH_TIMER("formRHS", t2);
  CH_TIMER("solve", t3);
  CH_TIMER("output", t4);
#ifdef CH_USE_PETSC
  const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
  DataIterator dit( dbl );
  PetscErrorCode ierr;
  const PetscInt nc = a_rhs.nComp();
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
  // setup if needed
  ierr = setup_solver( a_phi ); CHKERRQ(ierr);
  CH_START(t2);
  // form RHS -- m_bb
  Real idx2 = 1.e0/(this->m_dx*this->m_dx);
  Real addV;

  //rhsOperation(a_rhs); 
  // this is an interface in case some operations are needed on the rhs
  // the default does nothing

  // add X and B from Chombo to PETSc and add stuff for EB to B
  ierr = VecSet( m_xx, 0.); CHKERRQ(ierr);
  ierr = VecSet( m_bb, 0.); CHKERRQ(ierr);
  for ( dit = a_rhs.dataIterator() ; dit.ok() ; ++dit )
    {
      const BaseFab<Real> &rhsFab = getRegFab(a_rhs,dit);
      const BaseFab<Real> &xFab = getRegFab(a_phi,dit);
      const Box& box = dbl.get(dit());
      // const EBCellFAB& ebfab = a_rhs[dit()];
      // const EBISBox& ebbox   = ebfab.getEBISBox();
      //const Box &box = ebbox.getRegion(); // for Helmholtz
      const BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); bit.next())
        {
          IntVect iv = bit();
          PetscInt mm, ki = nc*gidsFab(iv,0);
          for (mm=0;mm<nc;mm++,ki++)
            {
              PetscScalar v = rhsFab(iv,mm);
              //if (ki==2)pout() << "    *************************** b1(2)=" << v << endl;
              //  v += addBCrhsValue(iv,a_phi,dit,idx2);
              addV = addBCrhsValue(iv,a_phi,dit,idx2);
              if (addV == BaseFabRealSetVal)
                {
                  v = 0;
                }
              else
                {
                  v += addV;
                }
              ierr = VecSetValues(m_bb,1,&ki,&v,INSERT_VALUES); CHKERRQ(ierr);
              //if (ki==2)pout() << "    *************************** b2(2)=" << v << endl;
              //if ( nonZeroInit )
              //{
              v = xFab(iv,mm);
              ierr = VecSetValues(m_xx,1,&ki,&v,INSERT_VALUES); CHKERRQ(ierr);
              //}
            }
        }
    }
  ierr = VecAssemblyBegin( m_bb ); CHKERRQ(ierr);
  ierr = VecAssemblyEnd( m_bb ); CHKERRQ(ierr);
  //if (nonZeroInit)
  //{
  ierr = VecAssemblyBegin( m_xx ); CHKERRQ(ierr);
  ierr = VecAssemblyEnd( m_xx ); CHKERRQ(ierr);
  //}
  CH_STOP(t2);
  // null space for periodic BCs
  if ( m_null )
  {
    MatNullSpace nullsp;
    ierr = MatNullSpaceCreate(wcomm,PETSC_TRUE,0,PETSC_NULL,&nullsp);CHKERRQ(ierr);
    CH_assert(m_ksp); // not used yet, needs fix for nonlinear
    ierr = KSPSetNullSpace( m_ksp, nullsp ); CHKERRQ(ierr);
  }
  // solve
  CH_START(t3);
  if ( m_snes )
    {
      ierr = SNESSolve( m_snes, m_bb, m_xx );     CHKERRQ(ierr);

    }
  else
    {
      ierr = KSPSolve( m_ksp, m_bb, m_xx );     CHKERRQ(ierr);
    }
  CH_STOP(t3);
  // put solution into output
  CH_START(t4);
  ierr = putPetscInChombo( a_phi, m_xx );     CHKERRQ(ierr);
  a_phi.exchange();
  CH_STOP(t4);
#endif
  return 0;
}
#ifdef CH_USE_PETSC
template <class T>
PetscErrorCode PetscSolver<T>::apply_mfree(Mat A, Vec x, Vec f)
{
  CH_TIME("PetscSolver::apply_mfree");
  //Whenever I write any PETSc code, I look forward to pulling classes back from the void.
  PetscFunctionBegin;
  void *ctx;
  MatShellGetContext(A, &ctx);
  PetscSolver<T> *tthis = (PetscSolver<T> *)ctx;
  tthis->putPetscInChombo( tthis->m_phi_mfree, x );
  tthis->m_op_mfree->applyOp(tthis->m_Lphi_mfree,tthis->m_phi_mfree,tthis->m_homogeneous);
  tthis->putChomboInPetsc(f, tthis->m_Lphi_mfree );
  PetscFunctionReturn(0);
}
#endif
template <class T>
int PetscSolver<T>::solve_mfree_private( T& a_phi,
                                         const T& a_rhs, 
                                         LinearOp<T> *a_op )
{
  CH_TIME("PetscSolver::solve_mfree_private");
#ifdef CH_USE_PETSC
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
#else
  MPI_Comm wcomm = PETSC_COMM_SELF;
#endif
  PetscErrorCode ierr;
  // setup if needed
  ierr = setup_solver( a_phi ); CHKERRQ(ierr);

  //create an operator matrix shell with same dimensions as m_mat
  PetscInt m, n, M, N;
  ierr = MatGetSize(m_mat, &M, &N);CHKERRQ(ierr); CH_assert(M == N);
  ierr = MatGetLocalSize(m_mat, &m, &n);CHKERRQ(ierr);
  Mat L; 
  ierr = MatCreateShell(wcomm,m,n,N,N,(void *)this,&L);CHKERRQ(ierr);
  ierr = MatShellSetOperation(L,MATOP_MULT,(void(*)(void))apply_mfree);
  m_op_mfree = a_op;
  //allocate space for a vector and a matrix-vector product in Chombo-land
  a_op->create( m_phi_mfree , a_phi);
  a_op->create( m_Lphi_mfree , a_rhs);
  
  // form RHS -- m_bb
  Real idx2 = 1.e0/(this->m_dx*this->m_dx);
  Real addV;

  // add X and B from Chombo to PETSc and add stuff for EB to B. This is just copy-pasted...
  ierr = VecSet( m_xx, 0.); CHKERRQ(ierr);
  ierr = VecSet( m_bb, 0.); CHKERRQ(ierr);
  const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout(); 
  const PetscInt nc = a_rhs.nComp();
  for (DataIterator dit = a_rhs.dataIterator() ; dit.ok() ; ++dit )
    {
      const BaseFab<Real> &rhsFab = getRegFab(a_rhs,dit);
      const BaseFab<Real> &xFab = getRegFab(a_phi,dit);
      const Box& box = dbl.get(dit());
      const BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); bit.next())
        {
          IntVect iv = bit();
          PetscInt mm, ki = nc*gidsFab(iv,0);
          for (mm=0;mm<nc;mm++,ki++)
            {
              PetscScalar v = rhsFab(iv,mm);
              addV = addBCrhsValue(iv,a_phi,dit,idx2);
              if (addV == BaseFabRealSetVal)
                {
                  v = 0;
                }
              else
                {
                  v += addV;
                }
              ierr = VecSetValues(m_bb,1,&ki,&v,INSERT_VALUES); CHKERRQ(ierr);
              v = xFab(iv,mm);
              ierr = VecSetValues(m_xx,1,&ki,&v,INSERT_VALUES); CHKERRQ(ierr);
             
            }
        }
    }
  ierr = VecAssemblyBegin( m_bb ); CHKERRQ(ierr);
  ierr = VecAssemblyEnd( m_bb ); CHKERRQ(ierr);
  ierr = VecAssemblyBegin( m_xx ); CHKERRQ(ierr);
  ierr = VecAssemblyEnd( m_xx ); CHKERRQ(ierr);
  //the ksp needed here is a bit different from m_ksp, becasue the preconditioner and main operator are different 
  KSP ksp;
  KSPCreate( wcomm, &ksp );
  if ( strlen(m_prestring) > 0 )
    {
      ierr = KSPSetOptionsPrefix( m_ksp, m_prestring );CHKERRQ(ierr);
    }
  ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
  PetscBool ism;
  PetscOptionsGetBool(m_prestring,"-ksp_monitor",&ism,PETSC_NULL);
  if (ism)
    {
      ierr = KSPMonitorSet(ksp,ksp_monitor_pout,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
    }
#if PETSC_VERSION_GE(3,5,0)
  ierr = KSPSetOperators(ksp,L,m_mat); CHKERRQ(ierr);
#else
  ierr = KSPSetOperators(ksp,L,m_mat,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
#endif

  ierr = KSPSetInitialGuessNonzero(ksp, m_nz_init_guess ? PETSC_TRUE : PETSC_FALSE ); CHKERRQ(ierr);

  if (0){ // coordinates
    PC pc; PetscInt sz,ind,bs,n,m;
#if PETSC_VERSION_LT(3,4,0) & PETSC_VERSION_RELEASE
    const PCType type;
#else
    PCType type;
#endif  
    ierr = KSPGetPC( ksp, &pc );     CHKERRQ(ierr);
    ierr = PCGetType( pc, &type );    CHKERRQ(ierr);
    ierr = MatGetBlockSize( m_mat, &bs );               CHKERRQ( ierr );
    if ( strcmp(type,PCGAMG) == 0 && bs > 1 )
      {
        PetscReal    *coords;
        DataIterator dit(a_phi.disjointBoxLayout());
        
        ierr = MatGetLocalSize( m_mat, &m, &n );  CHKERRQ(ierr);
        sz = CH_SPACEDIM*(m/bs);
        ierr = PetscMalloc( (sz+1)*sizeof(PetscReal), &coords ); CHKERRQ(ierr);
        for ( dit = a_phi.dataIterator(), ind = 0 ; dit.ok() ; ++dit )
          {
            const Box &box = a_phi.getBoxes()[dit];
            BoxIterator bit(box);
            for (bit.begin(); bit.ok(); bit.next())
              {
                IntVect iv = bit(); // coordinate in any scaled, shifted, rotated frame.
                for (PetscInt k=0; k<CH_SPACEDIM; k++) coords[ind++] = (PetscReal)iv[k];
              }
          }
        CH_assert(ind==sz);
        ierr = PCSetCoordinates( pc, CH_SPACEDIM, sz/CH_SPACEDIM, coords ); CHKERRQ(ierr);
        ierr = PetscFree( coords );  CHKERRQ(ierr);
      }
  }

  //carry out the solve
  ierr = KSPSolve( ksp, m_bb, m_xx );CHKERRQ(ierr);

  //put solution into output
  ierr = putPetscInChombo( a_phi, m_xx );CHKERRQ(ierr);
  a_phi.exchange();

  //clean up 
  a_op->clear( m_phi_mfree);
  a_op->clear( m_Lphi_mfree);
  KSPDestroy( &ksp );
#endif
  return 0;
}

#ifdef CH_USE_PETSC
// *******************************************************
template <class T>
PetscErrorCode PetscSolver<T>::putPetscInChombo( T& a_phi, const Vec xx )
{
  PetscErrorCode ierr;
  PetscScalar *arr;
  const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
  DataIterator dit( dbl );
  const PetscInt nc = a_phi.nComp();

  ierr = VecGetArray(xx,&arr);  CHKERRQ(ierr);
  for ( dit = a_phi.dataIterator() ; dit.ok() ; ++dit )
    {
      BaseFab<Real> &phiFab = getRegFab(a_phi,dit);
      const Box& box = dbl.get(dit());
      const BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); bit.next())
        {
          IntVect iv = bit();
          PetscInt mm, ki = nc*gidsFab(iv,0);
          for (mm=0;mm<nc;mm++,ki++)
            {
              phiFab(iv,mm) = arr[ki - nc*m_gid0];
            }
        }
    }
  ierr =  VecRestoreArray(xx,&arr); CHKERRQ(ierr);
  return 0;
}

// *******************************************************
template <class T>
PetscErrorCode PetscSolver<T>::putChomboInPetsc( Vec out, const T& a_phi )
{
  PetscErrorCode ierr;
  const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
  DataIterator dit( dbl );
  const PetscInt nc = a_phi.nComp();
  Real idx2 = 1.e0/(this->m_dx*this->m_dx);

  // add BC stuff to RHS (EB)
  ierr = VecSet( out, 0.);  CHKERRQ(ierr);
  for ( dit = a_phi.dataIterator() ; dit.ok() ; ++dit )
    {
      const BaseFab<Real> &phiFab = getRegFab(a_phi,dit);
      const Box& box = dbl.get(dit());
      const BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); bit.next())
        {
          IntVect iv = bit();
          PetscInt mm, ki = nc*gidsFab(iv,0);
          for (mm=0;mm<nc;mm++,ki++)
            {
              PetscScalar v = phiFab(iv,mm);
              v += addBCrhsValue(iv,a_phi,dit,idx2);
              if (isnan(v)) v = 0;
              ierr = VecSetValues(out,1,&ki,&v,INSERT_VALUES);
            }
        }
    }
  ierr = VecAssemblyBegin( out );  CHKERRQ(ierr);
  ierr = VecAssemblyEnd( out );  CHKERRQ(ierr);

  return 0;
}
#endif

// *******************************************************
//  computes m_bb - A m_xx, utility method (not used)
//
template <class T>
Real PetscSolver<T>::computeResidual()
{
  Vec tempVec;
  LevelData<FArrayBox > Residual;
  PetscScalar alpha = -1;
  PetscErrorCode ierr;

  ierr = VecDuplicate(m_xx, &tempVec);CHKERRQ(ierr);
  ierr = MatMult(m_mat, m_xx, tempVec);CHKERRQ(ierr);
  ierr = VecAXPY(tempVec,alpha,m_bb);CHKERRQ(ierr);

  IntVect idghosts = m_gids.ghostVect();
  const DisjointBoxLayout &dbl = m_gids.disjointBoxLayout();
  Residual.define( dbl, 0, idghosts );

  ierr = putPetscInChombo( Residual, tempVec );CHKERRQ(ierr);
  VecDestroy( &tempVec );
  Residual.exchange();

  //viewBFR(&Residual );

  return 0;
}

// *******************************************************
//  computes |phi|_inf
//
template <class T>
Real PetscSolver<T>::normInfinity( const T& a_phi ) const
{
  PetscErrorCode ierr;
  PetscReal      norm;
  
  ierr = putPetscInChombo( a_phi, m_rr );     //CHKERRQ(ierr);
  ierr = VecNorm(m_rr,NORM_INFINITY,&norm);CHKERRQ(ierr);
  //viewBFR(&Residual );
  return norm;
}

// *******************************************************
// PetscSolver<T>::applyOp
//  apply the matrix - use for debugging
//
template <class T>
int PetscSolver<T>::applyOp( T & a_out, const T & a_phi )
{
  //const int nc = a_phi.nComp();
  const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
  DataIterator dit( dbl );
  PetscErrorCode ierr;

  // setup if needed
  ierr = create_mat_vec( a_phi );

  // add BC stuff to RHS (EB)
  ierr = putChomboInPetsc( m_bb, a_phi );  CHKERRQ(ierr);

  // apply op
  ierr = MatMult( m_mat, m_bb, m_xx );  CHKERRQ(ierr);

  // get data back to Chombo
  ierr = putPetscInChombo( a_out, m_xx );     CHKERRQ(ierr);

  return ierr;
}
#include "NamespaceFooter.H"
#endif
