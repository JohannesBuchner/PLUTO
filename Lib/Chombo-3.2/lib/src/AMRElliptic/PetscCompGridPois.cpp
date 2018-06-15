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

#include "PetscCompGridPois.H"
#include "NamespaceHeader.H"

// derived Poisson class

void
PetscCompGridPois::clean()
{
  PetscCompGrid::clean();
}

void 
PetscCompGridPois::createOpStencil( IntVect a_iv, int a_ilev,const DataIndex &a_di_dummy, 
                                    StencilTensor &a_sten)
{
  CH_TIME("PetscCompGridPois::createOpStencil");
  Real dx=m_dxs[a_ilev][0],idx2=1./(dx*dx);
  
  if (!isCornerStencil() || SpaceDim==1)
    {
      if (m_order==2)
        {
          StencilTensorValue &v0 = a_sten[IndexML(a_iv,a_ilev)]; v0.define(COMP_POIS_DOF);
          for (int i=0; i<COMP_POIS_DOF; ++i) v0.setValue(i,i,m_alpha - m_beta*2.*SpaceDim*idx2);
          for (int dir=0; dir<CH_SPACEDIM; ++dir)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int isign = sign(sit());
                  IntVect jiv(a_iv); jiv.shift(dir,isign);
                  StencilTensorValue &v1 = a_sten[IndexML(jiv,a_ilev)]; v1.define(COMP_POIS_DOF);
                  for (int i=0; i<COMP_POIS_DOF; ++i) v1.setValue(i,i,m_beta*idx2);
                }
            }
        }
      else
        {         
          StencilTensorValue &v0 = a_sten[IndexML(a_iv,a_ilev)]; v0.define(COMP_POIS_DOF);
          for (int i=0; i<COMP_POIS_DOF; ++i) v0.setValue(i,i,m_alpha - m_beta*(30./12.)*SpaceDim*idx2);
          for (int dir=0; dir<CH_SPACEDIM; ++dir)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  int isign = sign(sit());
                  IntVect jiv(a_iv); jiv.shift(dir,isign);
                  StencilTensorValue &v1 = a_sten[IndexML(jiv,a_ilev)]; v1.define(COMP_POIS_DOF);
                  for (int i=0; i<COMP_POIS_DOF; ++i) v1.setValue(i,i,m_beta*(16./12.)*idx2);
                  jiv.shift(dir,isign);
                  StencilTensorValue &v2 = a_sten[IndexML(jiv,a_ilev)]; v2.define(COMP_POIS_DOF);
                  for (int i=0; i<COMP_POIS_DOF; ++i) v2.setValue(i,i,-m_beta*(1./12.)*idx2);
                }
            }
        }
    }
  else
    {
      if (m_order!=2)
        {
          MayDay::Error("PetscCompGridPois::createOpStencil - only 2nd order implemented with corner stencils");
        }
      Real stenvals[2][4] = {{20./6., 4./6., 1./6., 0}, {64./15., 7./15, 1./10., 1./30.}};

      StencilTensorValue &v0 = a_sten[IndexML(a_iv,a_ilev)]; v0.define(COMP_POIS_DOF);
      for (int i=0; i<COMP_POIS_DOF; ++i) v0.setValue(i,i,m_alpha - stenvals[CH_SPACEDIM-2][0]*m_beta*idx2);
      for (int dir1=0; dir1<CH_SPACEDIM; ++dir1)
        for (SideIterator sit1; sit1.ok(); ++sit1)
          {
            int isign1 = sign(sit1());
            IntVect jiv1(a_iv); jiv1.shift(dir1,isign1);
            StencilTensorValue &v1 = a_sten[IndexML(jiv1,a_ilev)]; v1.define(COMP_POIS_DOF);
            for (int i=0; i<COMP_POIS_DOF; ++i) v1.setValue(i,i,stenvals[CH_SPACEDIM-2][1]*m_beta*idx2);

            for (int dir2=dir1+1; dir2<CH_SPACEDIM; ++dir2)
              for (SideIterator sit2; sit2.ok(); ++sit2)
                {
                  int isign2 = sign(sit2());
                  IntVect jiv2(jiv1); jiv2.shift(dir2,isign2);
                  StencilTensorValue &v2 = a_sten[IndexML(jiv2,a_ilev)]; v2.define(COMP_POIS_DOF);
                  for (int i=0; i<COMP_POIS_DOF; ++i) v2.setValue(i,i,stenvals[CH_SPACEDIM-2][2]*m_beta*idx2);

                  for (int dir3=dir2+1; dir3<CH_SPACEDIM; ++dir3)
                    for (SideIterator sit3; sit3.ok(); ++sit3)
                      {
                        int isign3 = sign(sit3());
                        IntVect jiv3(jiv2); jiv3.shift(dir3,isign3);
                        StencilTensorValue &v3 = a_sten[IndexML(jiv3,a_ilev)]; v3.define(COMP_POIS_DOF);
                        for (int i=0; i<COMP_POIS_DOF; ++i) v3.setValue(i,i,stenvals[CH_SPACEDIM-2][3]*m_beta*idx2);
                      }
                }
          }
    }
}

#include "NamespaceFooter.H"

#endif
