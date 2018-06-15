#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBImplicitIntegrationStrategy.H"
#include "LevelTGAF_F.H" // For fortran stuff

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
EBImplicitIntegrationStrategy::
EBImplicitIntegrationStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
EBImplicitIntegrationStrategy::
~EBImplicitIntegrationStrategy()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
EBImplicitIntegrationStrategy::
setSourceGhostCells(LevelData<EBCellFAB>& a_src)
{
  int ncomp = a_src.nComp();
  const DisjointBoxLayout& dbl = this->disjointBoxLayout();
  const EBLevelGrid& grid = this->grid();
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
  {
    const Box& box   = dbl.get(dit());
    const Box& srcBox = a_src[dit()].box();
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        int iside = sign(sit());
        Box bc_box = adjCellBox(box, idir, sit(), 1);

        for (int jdir = 0; jdir < SpaceDim; jdir++)
        {
          //want corners too
          if (jdir != idir)
          {
            bc_box.grow(jdir, 1);
          }
        }

        //if fails might not have a ghost cell.
        bc_box &= grid.getDomain().domainBox();

        CH_assert(srcBox.contains(bc_box));
        if (box.size(idir) >= 4)
        {
          FORT_HORESGHOSTBC(CHF_FRA(a_src[dit()].getSingleValuedFAB()),
              CHF_BOX(bc_box),
              CHF_CONST_INT(idir),
              CHF_CONST_INT(iside),
              CHF_CONST_INT(ncomp));
        }
        else
        {
          // valid region not wide enough to apply HOExtrap -- drop
          // to linear extrap
          FORT_RESGHOSTBC(CHF_FRA(a_src[dit()].getSingleValuedFAB()),
              CHF_BOX(bc_box),
              CHF_CONST_INT(idir),
              CHF_CONST_INT(iside),
              CHF_CONST_INT(ncomp));
        }
        IntVectSet ivs = grid.getEBISL()[dit()].getIrregIVS(bc_box);
        for (VoFIterator vofit(ivs, grid.getEBISL()[dit()].getEBGraph()); vofit.ok(); ++vofit)
        {
          for (int icomp = 0; icomp < ncomp; icomp++)
          {
            Real valNeigh = 0;
            Vector<FaceIndex> faces = grid.getEBISL()[dit()].getFaces(vofit(), idir, flip(sit()));
            for (int iface = 0; iface < faces.size(); iface++)
            {
              VolIndex vofNeigh = faces[iface].getVoF(flip(sit()));
              valNeigh += a_src[dit()](vofNeigh, icomp);
            }
            if (faces.size() > 1) valNeigh /= faces.size();
            a_src[dit()](vofit(), icomp) = valNeigh;
          }
        }
      }
    }
  }
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
