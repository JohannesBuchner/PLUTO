#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LevelData.H"
#include "FluxBox.H"
#include "BoxIterator.H"

#include "PetscSolver.H"
#include "ViscousTensorOp.H"

#include "NamespaceHeader.H"

// *******************************************************
// PetscSolverPoisson Implementation
// *******************************************************
template< >
PetscSolverPoisson<LevelData<FArrayBox> >::PetscSolverPoisson()
  :PetscSolverFAB<LevelData<FArrayBox> >(),
   m_alpha(0.),
   m_beta(1.0)
{
}

// *******************************************************
template< >
void PetscSolverPoisson<LevelData<FArrayBox> >::define( Real a_dx,
                                                        bool a_homogeneous)
{
  PetscSolver<LevelData<FArrayBox> >::define(a_dx,a_homogeneous);
}

#ifdef CH_USE_PETSC
// *******************************************************
//   PetscSolverPoisson<LevelData<FArrayBox> >::formMatrix
//
template< >
int PetscSolverPoisson<LevelData<FArrayBox> >::formMatrix( Mat a_mat, const LevelData<FArrayBox> *a_dummy,
                                                           PetscInt dummy_my0, PetscInt dummy_nloc,
                                                           PetscInt *dummy_d, PetscInt *dummy_o )
{
  CH_assert(this->m_dx!=0.0);
  int ierr, nc = 1;
  Real idx2 = 1.e0/(this->m_dx*this->m_dx) * this->m_beta;
  Real vo = 1.e0 * idx2;
  Real vd = -CH_SPACEDIM * 2.e0 * idx2 + this->m_alpha;

  DisjointBoxLayout dbl = m_gids.disjointBoxLayout();
  for ( DataIterator dit = m_gids.dataIterator() ; dit.ok() ; ++dit )
    {
      const Box &box = dbl[dit];
      const BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); bit.next())
        {
          IntVect iv = bit();
          PetscInt ki = nc*gidsFab(iv,0);
          PetscInt i = ki;
          // diag
          ierr = MatSetValues(a_mat,1,&i,1,&i,&vd,ADD_VALUES); CHKERRQ(ierr);

          for ( int dim_dir = 0 ; dim_dir < CH_SPACEDIM ; dim_dir++ )
            {
              for ( int i_dir = -1 ; i_dir <= 1 ; i_dir += 2 )
                {
                  IntVect shiftiv = IntVect::Zero ;
                  shiftiv[dim_dir] = i_dir;
                  IntVect jv = iv + shiftiv;
                  PetscInt j = nc*gidsFab(jv,0);
                  if ( j >= 0 )
                    {
                      ierr = MatSetValues(a_mat,1,&i,1,&j,&vo,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else
                    {
                      vo *= -1.;
                      ierr = MatSetValues(a_mat,1,&i,1,&i,&vo,ADD_VALUES); CHKERRQ(ierr);
                      vo *= -1.;
                    }
                }
            }
        }
    }
  ierr = MatAssemblyBegin(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  //     PetscViewer viewer;
  //     PetscViewerASCIIOpen( wcomm, "A.m", &viewer);
  //     PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  //     MatView( a_mat, viewer );
  //     PetscViewerDestroy(viewer);
  //     MPI_Barrier(wcomm);
  //     exit(12);
  return 0;
}
#endif
// *******************************************************
// PetscSolverViscousTensor Implementation
// *******************************************************
template < >
PetscSolverViscousTensor<LevelData<FArrayBox> >::PetscSolverViscousTensor()
  :PetscSolverFAB<LevelData<FArrayBox> >(),
   m_alpha(1.),
   m_beta(1.),
   m_dxCrse(0.),
   m_a(0),
   m_eta(0),
   m_lamb(0)
{
}

// *******************************************************
template < >
void PetscSolverViscousTensor<LevelData<FArrayBox> >::define( LinearOp<LevelData<FArrayBox> >* a_operator,
                                                              bool a_homogeneous )
{
  PetscSolver<LevelData<FArrayBox> >::define( a_operator, a_homogeneous );
  ViscousTensorOp *op = dynamic_cast<ViscousTensorOp*>( a_operator );  CH_assert(op);
  m_dxCrse = op->dxCrse();

  // as a start, grab eta, lamdba, and 
  define( op->getAlpha(), op->getBeta(), 
          &(*op->getACoef()),
          &(*op->getEta()),
          &(*op->getLambda()));
}

#ifdef CH_USE_PETSC
// *******************************************************
//  PetscSolverViscousTensor<LevelData<FArrayBox> >::formMatrix
//
template < >
int PetscSolverViscousTensor<LevelData<FArrayBox> >::formMatrix( Mat a_mat, 
                                                                 const LevelData<FArrayBox> *a_dummy,
                                                                 PetscInt dummy_my0, PetscInt dummy_nloc,
                                                                 PetscInt *dummy_d, PetscInt *dummy_o )
{
  CH_assert(m_dx>0.0);
  CH_assert(CH_SPACEDIM==2);

  int  ierr,nc=2;
  Real hinv2=m_beta/(m_dx*m_dx), dxc=m_dxCrse, foff=(m_dx - dxc)/(3.*m_dx + dxc), fdia=2.*(dxc - m_dx)/(dxc + m_dx);

  if (!m_eta)
    {
      MayDay::Error("m_eta not set in PetscSolverViscousTensor -- m_eta must be set.");
    }
  DisjointBoxLayout dbl = m_gids.disjointBoxLayout();
  for ( DataIterator dit = m_gids.dataIterator() ; dit.ok() ; ++dit )
    {
      const Box &box = dbl[dit];
      const BaseFab<PetscInt> &gidsFab = this->m_gids[dit]; 
      const FluxBox &etaFab = (*m_eta)[dit()];
      const FArrayBox &eta_x = etaFab[0];
      const FArrayBox &eta_y = etaFab[1];
      //
      BoxIterator bit(box);
      for (bit.begin(); bit.ok(); bit.next())
        {
          Real eta_x0,eta_x1,lam_x0,lam_x1,eta_y0,eta_y1,lam_y0,lam_y1;
          IntVect iv = bit();
          PetscInt ki=nc*gidsFab(iv,0),jj;
          PetscInt kiwst,kiest,kjsth,kjnth,kIwstJnth,kIestJnth,kIwstJsth,kIestJsth;
          // get IV for stencil
          {
            IntVect ilft = IntVect::Zero, iest = IntVect::Zero, jsth = IntVect::Zero, jnth = IntVect::Zero, tempi = IntVect::Zero;
            ilft[0] = -1; iest[0] = 1; jsth[1] = -1; jnth[1] = 1;
            ilft += iv; iest += iv; jsth += iv; jnth += iv;
            // get PETSc indices for start of block -- PETSc will ignore egative indices so this is sortof D. BC
            kiwst  = nc*gidsFab(ilft,0);
            kiest = nc*gidsFab(iest,0);
            kjsth    = nc*gidsFab(jsth,0);
            kjnth    = nc*gidsFab(jnth,0);
            tempi[0] = -1; tempi[1] = 1; tempi += iv;
            kIwstJnth = nc*gidsFab(tempi,0);
            tempi[0] = 1; tempi[1] = 1; tempi += iv;
            kIestJnth = nc*gidsFab(tempi,0);
            tempi[0] = -1; tempi[1] = -1; tempi += iv;
            kIwstJsth = nc*gidsFab(tempi,0);
            tempi[0] = 1; tempi[1] = -1; tempi += iv;
            kIestJsth = nc*gidsFab(tempi,0);

            // get stencil coeficients
            eta_x0 = hinv2*eta_x(iv,0);     eta_y0 = hinv2*eta_y(iv,0);
            eta_x1 = hinv2*eta_x(iest,0);   eta_y1 = hinv2*eta_y(jnth,0);
            if ( m_lamb )
            {
              const FluxBox &lamFab = (*m_lamb)[dit()];
              const FArrayBox &lam_x = lamFab[0];
              const FArrayBox &lam_y = lamFab[1];
              lam_x0 = hinv2*lam_x(iv,0);     lam_y0 = hinv2*lam_y(iv,0);
              lam_x1 = hinv2*lam_x(iest,0); lam_y1 = hinv2*lam_y(jnth,0);
            }
            else
            {
              lam_x0 = 2.*eta_x0;
              lam_y0 = 2.*eta_y0;
              lam_x1 = 2.*eta_x1;
              lam_y1 = 2.*eta_y1;
            }
          }
          // loop over two equations, should probably just hard wire this
          for (int kk=0,vv=1;kk<2;kk++,vv--)
            {
              const PetscInt ii = ki + kk; // the row (cell) we are working on
              Real vd,tt;

              // add one eta everywhere -- same for u and v
              vd = eta_x0;
              if ( kiwst>=0 )
              {
                jj = kiwst + kk;         ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
              }
              else
                {
                  jj = kiest + kk; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  tt = fdia*vd;                  ierr = MatSetValues(a_mat,1,&ii,1,&ii,&tt,ADD_VALUES); CHKERRQ(ierr);
              }
              vd = eta_x1;
              if ( kiest>=0 )
              {
                jj = kiest + kk;        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
              }
              else
              {
                jj = kiwst + kk; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                tt = fdia*vd;                  ierr = MatSetValues(a_mat,1,&ii,1,&ii,&tt,ADD_VALUES); CHKERRQ(ierr);
              }
              vd = eta_y0;
              if ( kjsth>=0 )
              {
                jj = kjsth + kk;        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
              }
              else
              {
                jj = kjnth + kk; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                tt = fdia*vd;                  ierr = MatSetValues(a_mat,1,&ii,1,&ii,&tt,ADD_VALUES); CHKERRQ(ierr);
              }
              vd = eta_y1;
              if ( kjnth>=0 )
              {
                jj = kjnth + kk;        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
              }
              else
              {
                jj = kjsth + kk; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                tt = fdia*vd;                  ierr = MatSetValues(a_mat,1,&ii,1,&ii,&tt,ADD_VALUES); CHKERRQ(ierr);
              }
              jj = ii; vd = -(eta_x0 + eta_y0 + eta_x1 + eta_y1);
              ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);

              // add extra eta and the lambda for each direction
              if (kk==0)
              {
                  vd = eta_x0 + lam_x0;
                  if ( kiwst>=0 )
                  {
                    jj = kiwst + kk;        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                  }
                  else
                  {
                    jj = kiest + kk; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    tt = fdia*vd;                  ierr = MatSetValues(a_mat,1,&ii,1,&ii,&tt,ADD_VALUES); CHKERRQ(ierr);
                  }
                  //
                  vd = eta_x1 + lam_x1;
                  if ( kiest>=0 )
                  {
                    jj = kiest + kk;        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                  }
                  else
                  {
                    jj = kiwst + kk; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    tt = fdia*vd;                  ierr = MatSetValues(a_mat,1,&ii,1,&ii,&tt,ADD_VALUES); CHKERRQ(ierr);
                  }
                  jj = ii; vd = -(eta_x0 + eta_x1 + lam_x0 + lam_x1);
                  ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
              }
              else {
                vd = eta_y0 + lam_y0;
                if ( kjsth>=0 )
                  {
                    jj = kjsth + kk;     ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                  }
                else
                  {
                    jj = kjnth + kk; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    tt = fdia*vd;                  ierr = MatSetValues(a_mat,1,&ii,1,&ii,&tt,ADD_VALUES); CHKERRQ(ierr);
                  }
                //
                vd = eta_y1 + lam_y1;
                if ( kjnth>=0 )
                  {
                    jj = kjnth + kk;     ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                  }
                else
                  {
                    jj = kjsth + kk; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    tt = fdia*vd;                  ierr = MatSetValues(a_mat,1,&ii,1,&ii,&tt,ADD_VALUES); CHKERRQ(ierr);
                  }
                jj = ii; vd = -(eta_y0 + eta_y1 + lam_y0 + lam_y1);
                ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
              }

              // u -- v coupling terms -- 8 point stencil
              // S
              if (kk==0) vd = 0.25*(-lam_x1 + lam_x0);
              else vd = 0.25*(-eta_x1 + eta_x0);
              if ( kjsth>=0 )
                {
                  jj = kjsth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                }
              else
                {
                  jj = kjnth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = ki + vv;    tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                }

              // S.W.
              if (kk==0) vd = 0.25*(eta_y0 + lam_x0);
              else vd = 0.25*(eta_x0 + lam_y0);
              if ( kIwstJsth >= 0 && kiwst >= 0 && kjsth >= 0 )
                {
                  jj = kIwstJsth + vv;
                  ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                }
              else if ( kjsth >= 0 || kiwst >= 0 )
                {
                  if ( kIwstJsth >= 0 ) // reentrant corner: treat lam/eta like edge and eta/lamb like interior
                    {
                      if ( kiwst >= 0 )
                        {
                          if (kk==0) tt = 0.25*lam_x0;
                          else tt = 0.25*eta_x0;
                          jj = kIwstJsth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          // change 'vd' for "face" distribution
                          if (kk==0) vd = 0.25*eta_y0;
                          else vd = 0.25*lam_y0;
                          jj = kiwst + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          jj = kIwstJnth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                        }
                      else
                        {
                          if (kk==0) tt = 0.25*eta_y0;
                          else tt = 0.25*lam_y0;
                          jj = kIwstJsth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          // change 'vd' for "face" distribution
                          if (kk==0) vd = 0.25*lam_x0;
                          else vd = 0.25*eta_x0;
                          jj = kjsth + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          jj = kIestJsth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                        }
                    }
                  else if ( kiwst >= 0 && kjsth >= 0 ) // inner corner corner
                    {
                      //vd *= .5;
                      if (kk==0) vd = 0.25*eta_y0;
                      else vd = 0.25*lam_y0;
                      jj = kjsth + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIestJsth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      if (kk==0) vd = 0.25*lam_x0;
                      else vd = 0.25*eta_x0;
                      jj = kiwst + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIwstJnth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else if ( kjsth >= 0 ) // east face
                    {
                      jj = kjsth + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIestJsth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else if ( kiwst >= 0 )  // north face
                    {
                      jj = kiwst + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIwstJnth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                }
              else // ouside corner
                {
                  jj = kIestJnth + vv; tt = foff*foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = kiest + vv;     tt = fdia*foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = kjnth + vv;                        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = ki + vv;        tt = fdia*fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                }
              // S.E.
              if (kk==0) vd = 0.25*(-lam_x1 - eta_y0);
              else vd = 0.25*(-lam_y0 - eta_x1);
              if ( kIestJsth >= 0 && kiest >= 0 && kjsth >= 0 )
                {
                  jj = kIestJsth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                }
              else if ( kiest >= 0 || kjsth >= 0 )
                {
                  if ( kIestJsth >= 0 ) // reentrant corner: treat lam/eta like edge and eta/lamb like interior
                    {
                      if ( kiest >= 0 )
                        {
                          if (kk==0) tt = -0.25*lam_x1;
                          else tt = -0.25*eta_x1;
                          jj = kIestJsth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          // change 'vd' for "face" distribution
                          if (kk==0) vd = -0.25*eta_y0;
                          else vd = -0.25*lam_y0;
                          jj = kiest + vv;  tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          jj = kIestJnth + vv;  tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                        }
                      else
                        {
                          if (kk==0) tt = -0.25*eta_y0;
                          else tt = -0.25*lam_y0;
                          jj = kIestJsth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          // change 'vd' for "face" distribution
                          if (kk==0) vd = -0.25*lam_x1;
                          else vd = -0.25*eta_x1;
                          jj = kjsth + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          jj = kIwstJsth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                        }
                    }
                  else if ( kiest >= 0 && kjsth >= 0 )  // corner corner
                    {
                      if (kk==0) vd = -0.25*eta_y0;
                      else vd = -0.25*lam_y0;
                
                      jj = kjsth + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIwstJsth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      if (kk==0) vd = -0.25*lam_x1;
                      else vd = -0.25*eta_x1;
                      jj = kiest + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIestJnth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else if ( kiest >= 0 )
                    {
                      jj = kiest + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIestJnth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else if ( kjsth >= 0 )
                    {
                      jj = kjsth + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIwstJsth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                }
              else
                {
                  jj = kIwstJnth + vv; tt = foff*foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = kiwst + vv;     tt = fdia*foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = kjnth + vv;                        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = ki + vv;        tt = fdia*fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                }
              // W
              if (kk==0) vd = 0.25*(-eta_y1 + eta_y0);
              else vd = 0.25*(-lam_y1 + lam_y0);
              if ( kiwst>=0 )
                {
                  jj = kiwst + vv;   ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                }
              else
                {
                  jj = kiest + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = ki + vv;    tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                }
              // E
              if (kk==0) vd = 0.25*(eta_y1 - eta_y0);
              else vd = 0.25*(lam_y1 - lam_y0);
        
              if ( kiest>=0 )
                {
                  jj = kiest + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                }
              else
                {
                  jj = kiwst + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = ki + vv;    tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                }
              // N.W.
              if (kk==0) vd = 0.25*(-lam_x0 - eta_y1);
              else vd = 0.25*(-lam_y1 - eta_x0);
              if ( kIwstJnth >= 0 && kjnth >= 0 && kiwst >= 0 )
                {
                  jj = kIwstJnth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                }
              else if ( kjnth >= 0 || kiwst >= 0 )
                {
                  if ( kIwstJnth >= 0 ) // reentrant corner: treat lam/eta like edge and eta/lamb like interior
                    {
                
                       if ( kiwst >= 0 )
                        {
                          if (kk==0) tt = -0.25*lam_x0;
                          else tt = -0.25*eta_x0;
                          jj = kIwstJnth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          // change 'vd' for "face" distribution
                          if (kk==0) vd = -0.25*eta_y1;
                          else vd = -0.25*lam_y1;
                          jj = kiwst + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          jj = kIwstJsth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                        }
                      else
                        {
                          if (kk==0) tt = -0.25*eta_y1;
                          else tt = -0.25*lam_y1;
                          jj = kIwstJnth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          // change 'vd' for "face" distribution
                          if (kk==0) vd = -0.25*lam_x0;
                          else vd = -0.25*eta_x0;
                          jj = kjnth + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          jj = kIestJnth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                        }

                    }
                  else if ( kiwst >= 0 && kjnth >= 0 ) // corner corner
                    {
                      if (kk==0) vd = -0.25*eta_y1;
                      else vd = -0.25*lam_y1;
                      jj = kjnth + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIestJnth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      if (kk==0) vd = -0.25*lam_x0;
                      else vd = -0.25*eta_x0;
                      jj = kiwst + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIwstJsth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else if ( kjnth >= 0 )
                    {
                      jj = kjnth + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIestJnth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else if ( kiwst >= 0 )
                    {
                      jj = kiwst + vv;         tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIwstJsth + vv;     tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                }
              else
                {
                  jj = kIestJsth + vv; tt = foff*foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = kiest + vv;     tt = fdia*foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = kjsth + vv;                        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = ki + vv;        tt = fdia*fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                }
              // N
              if (kk==0) vd = 0.25*(lam_x1 - lam_x0);
              else vd = 0.25*(eta_x1 - eta_x0);
              if ( kjnth>=0 )
                {
                  jj = kjnth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                }
              else
                {
                  jj = kjsth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = ki + vv;    tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                }
              // N.E.
              if (kk==0) vd = 0.25*(eta_y1 + lam_x1);
              else vd = 0.25*(eta_x1 + lam_y1);
              if ( kIestJnth >= 0 && kiest >= 0 && kjnth >= 0 )
                {
                  jj = kIestJnth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
                }
              else if ( kiest >= 0 || kjnth >= 0 )
                {
                  if ( kIestJnth >= 0 ) // reentrant corner: treat lam/eta like edge and eta/lamb like interior
                    {
                      if ( kiest >= 0 )
                        {
                          if (kk==0) tt = 0.25*lam_x1;
                          else tt = 0.25*eta_x1;
                          jj = kIestJnth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          // change 'vd' for "face" distribution
                          if (kk==0) vd = 0.25*eta_y1;
                          else vd = 0.25*lam_y1;
                          jj = kiest + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          jj = kIestJsth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                        }
                      else
                        {
                          if (kk==0) tt = 0.25*eta_y1;
                          else tt = 0.25*lam_y1;
                          jj = kIestJnth + vv; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          // change 'vd' for "face" distribution
                          if (kk==0) vd = 0.25*lam_x1;
                          else vd = 0.25*eta_x1;
                          jj = kjnth + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                          jj = kIwstJnth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                        }
                    }
                  else if ( kiest >= 0 && kjnth >= 0 ) // corner corner
                    {
                      if (kk==0) vd = 0.25*eta_y1;
                      else vd = 0.25*lam_y1;
                      jj = kjnth + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIwstJnth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      if (kk==0) vd = 0.25*lam_x1;
                      else vd = 0.25*eta_x1;
                      jj = kiest + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIestJsth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else if ( kiest >= 0 )
                    {
                      jj = kiest + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIestJsth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                  else if ( kjnth >= 0 )
                    {
                      jj = kjnth + vv;     tt = fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                      jj = kIwstJnth + vv; tt = foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                    }
                }
              else
                {
                  jj = kIwstJsth + vv; tt = foff*foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = kiwst + vv;     tt = fdia*foff*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = kjsth + vv;                        ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                  jj = ki + vv;        tt = fdia*fdia*vd; ierr = MatSetValues(a_mat,1,&ii,1,&jj,&tt,ADD_VALUES); CHKERRQ(ierr);
                }
              // diagonal alpha term
              jj = ii; vd = m_alpha * ((m_a==0) ? 1. : (*m_a)[dit()](iv,0));
              ierr = MatSetValues(a_mat,1,&ii,1,&jj,&vd,ADD_VALUES); CHKERRQ(ierr);
            } // kk=1:2 loop
        } // bit()
    } // dit ()
  ierr = MatAssemblyBegin(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(a_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  //     PetscViewer viewer;
  //     PetscViewerASCIIOpen( wcomm, "A.m", &viewer);
  //     PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  //     MatView( a_mat, viewer );
  //     PetscViewerDestroy(&viewer);
  //     MPI_Barrier(wcomm);
  //     exit(12);
  return 0;
}
#endif

#include "NamespaceFooter.H"
