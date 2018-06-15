#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IFInspector.H"
#include "NormalDerivative.H"

#include "NamespaceHeader.H"

////////////////////

void IFInspector::fillValues(EBCellFAB     & a_ebFab,
                             const BaseIF  & a_imIF,
                             const RealVect& a_origin,
                             const RealVect& a_dx)
{
  //a particular box in the layout
  const Box& box = a_ebFab.box();

  //a particular box, which is linked to the results of computeVoFInternals
  const EBISBox& ebisBox = a_ebFab.getEBISBox();

  //This is all cells in the calculation
  IntVectSet ivs(box);
  for (VoFIterator vofit(ivs,ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect thisIv = vof.gridIndex();

      IndexTM<Real,GLOBALDIM> coord;

      for (int i=0;i<SpaceDim;i++)
        {
          coord[i]=a_origin[i]+thisIv[i]*a_dx[i];
        }
      a_ebFab(vof,0)=a_imIF.value(coord);
    }
}

////////////////////

void IFInspector::fillValues(LevelData<EBCellFAB>& a_levelFab,
                             const BaseIF& a_imIF,
                             const RealVect& a_origin,
                             const RealVect& a_dx)
{
  for (DataIterator dit = a_levelFab.dataIterator();dit.ok();++dit)
    {
      EBCellFAB & ebFab=a_levelFab[dit()];
      fillValues(ebFab,a_imIF,a_origin,a_dx);
    }
}

////////////////////

void IFInspector::fillNormal(EBCellFAB     & a_ebFab,
                             const BaseIF  & a_imIF,
                             const RealVect& a_origin,
                             const RealVect& a_dx)
{
  //a particular box in the layout
  const Box& box = a_ebFab.box();

  //a particular box, which is linked to the results of computeVoFInternals
  const EBISBox& ebisBox = a_ebFab.getEBISBox();

  //This is all cells in the calculation
  IntVectSet ivs(box);
  for (VoFIterator vofit(ivs,ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect thisIv = vof.gridIndex();

      IndexTM<Real,GLOBALDIM> coord;

      for (int i=0;i<SpaceDim;i++)  //todo: think out the third coordinate here and below
        {
          coord[i]=a_origin[i]+thisIv[i]*a_dx[i];
        }

      NormalDerivative<GLOBALDIM> normalDerivative;

      for (int idir=0; idir < GLOBALDIM; ++idir)
        {
          IndexTM<int,GLOBALDIM> zero = IndexTM<int,GLOBALDIM>::Zero;
          a_ebFab(vof,idir) = normalDerivative.evaluate(zero,
                                                        idir,
                                                        coord,
                                                        a_imIF);
        }
    }
}

////////////////////

void IFInspector::fillNormal(LevelData<EBCellFAB> & a_levelFab,
                             const BaseIF         & a_imIF,
                             const RealVect       & a_origin,
                             const RealVect       & a_dx)
{
  for (DataIterator dit = a_levelFab.dataIterator();dit.ok();++dit)
    {
      EBCellFAB & ebFab=a_levelFab[dit()];
      fillNormal(ebFab,a_imIF,a_origin,a_dx);
    }
}

////////////////////

void IFInspector::fillGradNormal(EBCellFAB     & a_ebFab,
                                 const BaseIF  & a_imIF,
                                 const RealVect& a_origin,
                                 const RealVect& a_dx)
{
  //a particular box in the layout
  const Box& box = a_ebFab.box();

  //a particular box, which is linked to the results of computeVoFInternals
  const EBISBox& ebisBox = a_ebFab.getEBISBox();

  //This is all cells in the calculation
  IntVectSet ivs(box);
  for (VoFIterator vofit(ivs,ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect thisIv = vof.gridIndex();

      IndexTM<Real,GLOBALDIM> coord;

      for (int i=0;i<SpaceDim;i++)
        {
          coord[i]=a_origin[i]+thisIv[i]*a_dx[i];
        }

      NormalDerivative<GLOBALDIM> normalDerivative;
      int k = 0;
      for (int idir = 0; idir < GLOBALDIM ; ++idir)
        {
          for (int jdir = 0; jdir < GLOBALDIM ; ++jdir)
            {
              IndexTM<int,GLOBALDIM> deriv = BASISV_TM<int,GLOBALDIM>(jdir);
              a_ebFab(vof,k++) = normalDerivative.evaluate(deriv,
                                                           idir,
                                                           coord,
                                                           a_imIF);
            }
        }
    }
}

////////////////////

void IFInspector::fillGradNormal(LevelData<EBCellFAB>& a_levelFab,
                                 const BaseIF& a_imIF,
                                 const RealVect& a_origin,
                                 const RealVect& a_dx)
{
  for (DataIterator dit = a_levelFab.dataIterator();dit.ok();++dit)
    {
      EBCellFAB & ebFab= a_levelFab[dit()];
      fillGradNormal(ebFab,a_imIF,a_origin,a_dx);
    }
}

#include "NamespaceFooter.H"
