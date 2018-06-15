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

#include "PetscCompGridVTO.H"
#include "FluxBox.H"
#include "NamespaceHeader.H"

// derived ViscousTensor Operator class

void
PetscCompGridVTO::clean()
{
  PetscCompGrid::clean();
}

void 
PetscCompGridVTO::createOpStencil(IntVect a_iv, int a_ilev,const DataIndex &a_di, StencilTensor &a_sten)
{
  CH_TIME("PetscCompGridVTO::createOpStencil");
  Real dx=m_dxs[a_ilev][0],idx2=1./(dx*dx);
  Real eta_x0,eta_x1,lam_x0,lam_x1,eta_y0,eta_y1,lam_y0,lam_y1;
  IntVect ivwst=IntVect::Zero,ivest=IntVect::Zero,ivsth=IntVect::Zero,ivnth=IntVect::Zero;
  D_TERM(,,);
  //IntVect ivSW(-1,-1),ivSE(1,-1),ivNE(1,1),ivNW(-1,1);
  IntVect ivSW(IntVect::Unit),ivSE(IntVect::Unit),ivNE(IntVect::Unit),ivNW(IntVect::Unit);
  ivSW[1] = ivSW[0] = -1; // >= 2D
  ivSE[1] = -1; // >= 2D
  ivNW[0] = -1;

  ivwst[0] = -1; ivest[0] = 1; ivsth[1] = -1; ivnth[1] = 1;
  {    // get stencil coeficients
    Real beta_hinv2=m_beta*idx2;
    const FluxBox &etaFab = (*m_eta[a_ilev])[a_di];
    const FArrayBox &eta_x = etaFab[0];
    const FArrayBox &eta_y = etaFab[1];
    const FluxBox &lamFab = (*m_lamb[a_ilev])[a_di];
    const FArrayBox &lam_x = lamFab[0];
    const FArrayBox &lam_y = lamFab[1];
    eta_x0 = beta_hinv2*eta_x(a_iv,0);           eta_y0 = beta_hinv2*eta_y(a_iv,0);
    eta_x1 = beta_hinv2*eta_x(ivest+a_iv,0);     eta_y1 = beta_hinv2*eta_y(ivnth+a_iv,0);
    lam_x0 = beta_hinv2*lam_x(a_iv,0);           lam_y0 = beta_hinv2*lam_y(a_iv,0);
    lam_x1 = beta_hinv2*lam_x(ivest+a_iv,0);     lam_y1 = beta_hinv2*lam_y(ivnth+a_iv,0);
  }

  // loop over two equations, should probably just hard wire this
  for (int kk=0,vv=1;kk<2;kk++,vv--)
    {
      Real vd;
      // add one eta everywhere -- same for u and v
      vd = eta_x0;
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivwst+a_iv,a_ilev)];  v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
      vd = eta_x1;
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivest+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
      vd = eta_y0;
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivsth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
      vd = eta_y1;
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivnth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
      vd = -(eta_x0 + eta_y0 + eta_x1 + eta_y1);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }

      // add extra eta and the lambda for each direction
      if (kk==0)
        {
          vd = eta_x0 + lam_x0;
          { 
            StencilTensorValue &v1 = a_sten[IndexML(ivwst+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
            v1.addValue(kk,kk,vd);
          }
          vd = eta_x1 + lam_x1;
          { 
            StencilTensorValue &v1 = a_sten[IndexML(ivest+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
            v1.addValue(kk,kk,vd);
          }
          vd = -(eta_x0 + eta_x1 + lam_x0 + lam_x1);
          { 
            StencilTensorValue &v1 = a_sten[IndexML(a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
            v1.addValue(kk,kk,vd);
          }
        }
      else {
        vd = eta_y0 + lam_y0;
        { 
          StencilTensorValue &v1 = a_sten[IndexML(ivsth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
          v1.addValue(kk,kk,vd);
        }
        vd = eta_y1 + lam_y1;
        { 
          StencilTensorValue &v1 = a_sten[IndexML(ivnth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
          v1.addValue(kk,kk,vd);
        }
        vd = -(eta_y0 + eta_y1 + lam_y0 + lam_y1);
        { 
          StencilTensorValue &v1 = a_sten[IndexML(a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
          v1.addValue(kk,kk,vd);
        }
      }

      // u -- v coupling terms -- 8 point stencil
      // S
      if (kk==0) vd = 0.25*(-lam_x1 + lam_x0);
      else vd = 0.25*(-eta_x1 + eta_x0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivsth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // S.W.
      if (kk==0) vd = 0.25*(eta_y0 + lam_x0);
      else vd = 0.25*(eta_x0 + lam_y0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivSW+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // S.E.
      if (kk==0) vd = 0.25*(-lam_x1 - eta_y0);
      else vd = 0.25*(-lam_y0 - eta_x1);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivSE+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // W
      if (kk==0) vd = 0.25*(-eta_y1 + eta_y0);
      else vd = 0.25*(-lam_y1 + lam_y0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivwst+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // E
      if (kk==0) vd = 0.25*(eta_y1 - eta_y0);
      else vd = 0.25*(lam_y1 - lam_y0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivest+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // N.W.
      if (kk==0) vd = 0.25*(-lam_x0 - eta_y1);
      else vd = 0.25*(-lam_y1 - eta_x0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivNW+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // N
      if (kk==0) vd = 0.25*(lam_x1 - lam_x0);
      else vd = 0.25*(eta_x1 - eta_x0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivnth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // N.E.
      if (kk==0) vd = 0.25*(eta_y1 + lam_x1);
      else vd = 0.25*(eta_x1 + lam_y1);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivNE+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // diagonal alpha term
      vd = m_alpha * (!(m_a[a_ilev]) ? 1. : (*m_a[a_ilev])[a_di](a_iv,0));
      { 
        StencilTensorValue &v1 = a_sten[IndexML(a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
    } // kk=1:2 loop
  
  // debug
  if (0)
    {
      Real summ=0.;
      StencilTensor::const_iterator end = a_sten.end(); 
      for ( StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it) {
        for (int i=0; i<CH_SPACEDIM; ++i) 
          {
            for (int j=0; j<CH_SPACEDIM; ++j) 
              {
                summ += it->second.value(i,j);
              }
          }
      }
      
      if (abs(summ>1.e-10))
        {
          for (StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it) {
            pout() << it->first << "  ERROR: \n";
            for (int i=0; i<CH_SPACEDIM; ++i) 
              {
                pout() << "\t\t";
                for (int j=0; j<CH_SPACEDIM; ++j) 
                  {
                    pout() << it->second.value(i,j) << ", ";
                  }
                pout() << endl;
              }
          }
          pout() << "\t ERROR summ: " << summ << endl; 
        }
    }
}

#include "NamespaceFooter.H"

#endif
