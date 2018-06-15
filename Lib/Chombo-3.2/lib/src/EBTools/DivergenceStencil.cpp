#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DivergenceStencil.H"
#include "DivergenceStencilF_F.H"
#include "NamespaceHeader.H"

DivergenceStencil::
DivergenceStencil(const EBCellFAB       & a_dvrgFluxFAB,
                  const EBFluxFAB       & a_faceFluxFAB,
                  const BaseIVFAB<Real> & a_bndyFluxFAB,
                  const Box             & a_grid,
                  const EBISBox         & a_ebisBox,
                  const RealVect        & a_dx,
                  bool a_useEBFlux)
{
  m_grid = a_grid;
  m_dx   = a_dx;
  m_useEBFlux = a_useEBFlux;

  IntVectSet ivs = a_ebisBox.getIrregIVS(a_grid);
  VoFIterator vofit(ivs, a_ebisBox.getEBGraph());
  Vector<VolIndex> vofs = vofit.getVector();
  Vector<RefCountedPtr<BaseIndex  > > vofptr(vofs.size());

  //get this of vofs
  for (int ivof = 0; ivof < vofs.size(); ivof++)
    {
      vofptr[ivof] = RefCountedPtr<BaseIndex  >(new VolIndex(vofs[ivof]));
    }

  //get the stencil for the bndry faces because that is easy
  Vector<RefCountedPtr<BaseStencil> > bndryStencils;
  if (m_useEBFlux)
    bndryStencils.resize(vofs.size());
  Vector<RefCountedPtr<BaseStencil> > faceStencils[CH_SPACEDIM];
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      faceStencils[idir].resize(vofs.size());
    }
  for (int ivof = 0; ivof < vofs.size(); ivof++)
    {
      VoFStencil* vofstenptr = new VoFStencil();
      Real bndryArea  = a_ebisBox.bndryArea(vofs[ivof]);
      RealVect normal = a_ebisBox.normal(   vofs[ivof]);
      for (int idir=0; idir < SpaceDim; idir++)
        {
          vofstenptr->add(vofs[ivof], -normal[idir]*bndryArea/a_dx[idir]);
        }
      if (m_useEBFlux)
        bndryStencils[ivof] = RefCountedPtr<BaseStencil>(vofstenptr);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FaceStencil* facestenptr = new FaceStencil();
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex> faces = a_ebisBox.getFaces(vofs[ivof], idir, sit());
              for (int iface = 0; iface < faces.size(); iface++)
                {
                  int isign = sign(sit());
                  Real areaFrac = a_ebisBox.areaFrac(faces[iface]);
                  facestenptr->add(faces[iface], isign*areaFrac/a_dx[idir]);
                }
            }
          faceStencils[idir][ivof] = RefCountedPtr<BaseStencil>(facestenptr);
        }
      delete vofstenptr;
    }
  if (a_useEBFlux)
    {
      m_bdryStencil = RefCountedPtr< AggStencil<BaseIVFAB<Real>, EBCellFAB> >
        (new AggStencil<BaseIVFAB<Real>, EBCellFAB>(vofptr,  bndryStencils,     a_bndyFluxFAB,       a_dvrgFluxFAB));
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_openStencil[idir] = RefCountedPtr< AggStencil<EBFaceFAB, EBCellFAB> >
        (new AggStencil<EBFaceFAB,   EBCellFAB>(vofptr, faceStencils[idir], a_faceFluxFAB[idir], a_dvrgFluxFAB));
    }
}

/**********/

void
DivergenceStencil::
divergence(EBCellFAB             & a_divF,
           const EBFluxFAB       & a_faceFlux,
           const BaseIVFAB<Real> & a_bndryFlux,
           const int             & a_destVar,
           bool a_incrementOnly)
{
  //call fortran assuming regular everywhere
  if (!a_incrementOnly)
    a_divF.setVal(0.);

  int ncomp = a_divF.nComp();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const BaseFab<Real>& flux = a_faceFlux[idir].getSingleValuedFAB();
      m_openStencil[idir]->cache(a_divF);

      FORT_DIVERGESTENF(CHF_BOX(m_grid),
                        CHF_FRA(a_divF.getSingleValuedFAB()),
                        CHF_CONST_FRA(flux),
                        CHF_CONST_INT(idir),
                        CHF_CONST_INT(ncomp),
                        CHF_CONST_REAL(m_dx[idir]));

      m_openStencil[idir]->uncache(a_divF);
    }
  //now do irregular calcs
  for (int icomp = 0; icomp < ncomp; icomp++)
    {

      if (m_useEBFlux)
        {
          m_bdryStencil->apply(a_divF, a_bndryFlux, icomp, true);
        }
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_openStencil[idir]->apply(a_divF, a_faceFlux[idir], icomp, true);
        }
    }

}
#include "NamespaceFooter.H"

