#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBDebugOut.H"
#include "EBCellFAB.H"
#include "LayoutIterator.H"
#include "parstream.H"
#include "FaceIndex.H"
#include "Vector.H"
#include "VolIndex.H"
#include "LoHiSide.H"
#include "VoFIterator.H"
#include "EBFluxFAB.H"
#include "DebugOut.H"
#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include "EBFaceFAB.H"
#include <iomanip>
#include "NamespaceHeader.H"
IntVect EBDebugPoint::s_ivd = IntVect::Zero;

void printMaxMinEBFF(EBFluxFAB* a_data)
{


  for (int idir = 0; idir < SpaceDim; idir ++)
    {
      for (int comp = 0; comp < a_data->nComp(); comp++)
        {
          Real maxVal = -1.e99;
          Real minVal =  1.e99;
          FaceIndex fmax, fmin;

          const EBFaceFAB& dataEBFAB = (*a_data)[idir];
          const Box&        dblBox   = a_data->getRegion();
          const IntVectSet ivsBox(dblBox);
          const EBISBox& ebisBox = dataEBFAB.getEBISBox();

          for (FaceIterator fit(ivsBox, ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
               fit.ok(); ++fit)
            {
              const FaceIndex& face = fit();
              const Real& val = dataEBFAB(face,comp);
              if (val > maxVal)
                {
                  maxVal = val;
                  fmax = face;
                }
              if (val < minVal)
                {
                  minVal = val;
                  fmin = face;
                }
            }
          pout() << " "
                 << setprecision(4)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 <<  "comp = " <<  comp
                 << ", dir = " <<  idir
                 << ", max = " << maxVal
                 << ", min = " << minVal
                 << ", imaxlo =" << fmax.gridIndex(Side::Lo)
                 << ", imaxhi =" << fmax.gridIndex(Side::Hi)
                 << ", iminlo =" << fmin.gridIndex(Side::Lo)
                 << ", iminhi =" << fmin.gridIndex(Side::Hi) << endl;
        }
    }
}
void printMaxMinLDFlux(LevelData<EBFluxFAB>* a_data)
{


  const DisjointBoxLayout& dbl = a_data->disjointBoxLayout();
  for (int idir = 0; idir < SpaceDim; idir ++)
    {
      for (int comp = 0; comp < a_data->nComp(); comp++)
        {
          Real maxVal = -1.e99;
          Real minVal =  1.e99;
          FaceIndex fmax, fmin;

          for (DataIterator dit=a_data->dataIterator();dit.ok();++dit)
            {
              const EBFaceFAB& dataEBFAB = (*a_data)[dit()][idir];
              const Box&        dblBox   = dbl.get(dit());
              const IntVectSet ivsBox(dblBox);
              const EBISBox& ebisBox = dataEBFAB.getEBISBox();
              for (FaceIterator fit(ivsBox, ebisBox.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
                   fit.ok(); ++fit)
                {
                  const FaceIndex& face = fit();
                  const Real& val = dataEBFAB(face,comp);
                  if (val > maxVal)
                    {
                      maxVal = val;
                      fmax = face;
                    }
                  if (val < minVal)
                    {
                      minVal = val;
                      fmin = face;
                    }
                }
            }
          pout() << " "
                 << setprecision(4)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 <<  "comp = " <<  comp
                 << ", dir = " <<  idir
                 << ", max = " << maxVal
                 << ", min = " << minVal
                 << ", imaxlo =" << fmax.gridIndex(Side::Lo)
                 << ", imaxhi =" << fmax.gridIndex(Side::Hi)
                 << ", iminlo =" << fmin.gridIndex(Side::Lo)
                 << ", iminhi =" << fmin.gridIndex(Side::Hi) << endl;
        }
    }
}
void printMaxMinLDCell(LevelData<EBCellFAB>* a_data)
{
  const DisjointBoxLayout& dbl = a_data->disjointBoxLayout();
  for (int comp = 0; comp < a_data->nComp(); comp++)
    {
      Real maxVal = -1.e99;
      Real minVal =  1.e99;
      IntVect bogusIV = -999*IntVect::Unit;
      VolIndex vofmax(bogusIV, -999);
      VolIndex vofmin(bogusIV, -999);
      for (DataIterator dit=a_data->dataIterator();dit.ok();++dit)
        {
          const EBCellFAB& dataEBFAB = (*a_data)[dit()];
          const Box&        dblBox   = dbl.get(dit());
          const IntVectSet ivsBox(dblBox);
          const EBISBox& ebisBox = dataEBFAB.getEBISBox();
          for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              const Real& val = dataEBFAB(vof,comp);
              if (val > maxVal)
                {
                  maxVal = val;
                  vofmax = vof;
                }
              if (val < minVal)
                {
                  minVal = val;
                  vofmin = vof;
                }
            }
        }
      pout() << " "
             << setprecision(4)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific)
             <<  "comp = " <<  comp
             << ", max = " << maxVal
             << ", min = " << minVal
             << ", imax =" << vofmax.gridIndex()
             << ", imin =" << vofmin.gridIndex() << endl;
    }
}
void printMaxMinEBCF(EBCellFAB* a_data)
{
  for (int comp = 0; comp < a_data->nComp(); comp++)
    {
      Real maxVal = -1.e99;
      Real minVal =  1.e99;
      IntVect bogusIV = -999*IntVect::Unit;
      VolIndex vofmax(bogusIV, -999);
      VolIndex vofmin(bogusIV, -999);
      const EBCellFAB& dataEBFAB = (*a_data);
      const Box&        dblBox   = a_data->getRegion();
      const IntVectSet ivsBox(dblBox);
      const EBISBox& ebisBox = dataEBFAB.getEBISBox();
      for (VoFIterator vofit(ivsBox, ebisBox.getEBGraph());vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const Real& val = dataEBFAB(vof,comp);
          if (val > maxVal)
            {
              maxVal = val;
              vofmax = vof;
            }
          if (val < minVal)
            {
              minVal = val;
              vofmin = vof;
            }
        }
      pout() << " "
             << setprecision(4)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific)
             <<  "comp = " <<  comp
             << ", max = " << maxVal
             << ", min = " << minVal
             << ", imax =" << vofmax.gridIndex()
             << ", imin =" << vofmin.gridIndex() << endl;
    }
}
void printPointEBCF(EBCellFAB* a_dat)
{

  IntVect iv = EBDebugPoint::s_ivd;
  VolIndex vof(iv, 0);
  if (a_dat->box().contains(iv))
    {
      pout() << "data for iv " << iv << ":";
      for (int ivar = 0; ivar < a_dat->nComp(); ivar++)
        {
          pout()               << setprecision(4)
                               << setiosflags(ios::showpoint)
                               << setiosflags(ios::scientific)
                               << " " << (*a_dat)(vof, ivar) << " ";
        }
      pout() << endl;
    }
}

void printLocalEBCF(EBCellFAB* a_dat)
{

  IntVect ivd = EBDebugPoint::s_ivd;
  Box bdeb(ivd, ivd);
  bdeb.grow(1);
  for (BoxIterator bit(bdeb); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if ((*a_dat).box().contains(iv))
        {
          Vector<VolIndex> vofs = a_dat->getEBISBox().getVoFs(iv);
          for (int ivof = 0; ivof < vofs.size(); ivof++)
            {
              VolIndex& vof = vofs[ivof];
              pout() << "data for iv " << iv << ":";
              for (int ivar = 0; ivar < a_dat->nComp(); ivar++)
                {
                  pout()                  << setprecision(4)
                                          << setiosflags(ios::showpoint)
                                          << setiosflags(ios::scientific)
                                          << " " << (*a_dat)(vof, ivar) << " ";
                }
              pout() << endl;
            }
        }
    }
}
void printPointLDCell(LevelData<EBCellFAB>* a_dat)
{
  int ibox= 0;
  for (DataIterator dit = a_dat->dataIterator(); dit.ok(); ++dit)
    {
      pout() << "box number "  << ibox << endl;
      printPointEBCF(&((*a_dat)[dit()]));
      ibox++;
    }
}
void printLocalLDCell(LevelData<EBCellFAB>* a_dat)
{
  int ibox= 0;
  for (DataIterator dit = a_dat->dataIterator(); dit.ok(); ++dit)
    {
      pout() << "box number "  << ibox << endl;
      printLocalEBCF(&((*a_dat)[dit()]));
      ibox++;
    }
}

void printPointEBFF(EBFluxFAB* a_dat)
{
  IntVectSet ivs(EBDebugPoint::s_ivd);
  if (a_dat->getRegion().contains(EBDebugPoint::s_ivd))
    {
      const EBGraph& ebgraph = (*a_dat).getEBISBox().getEBGraph();
      pout() << "flux data for iv " << EBDebugPoint::s_ivd << ":" << endl;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          for (FaceIterator faceit(ivs, ebgraph, idir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
            {
              pout() << "face" << faceit() << ": ";
              for (int ivar = 0; ivar < a_dat->nComp(); ivar++)
                {
                  pout() << " " << (*a_dat)[idir](faceit(), ivar) << " ";
                }
              pout() << endl;
            }
        }
    }
}

void printPointEBFace(EBFaceFAB* a_dat)
{
  IntVectSet ivs(EBDebugPoint::s_ivd);
  if (a_dat->getRegion().contains(EBDebugPoint::s_ivd))
    {
      const EBGraph& ebgraph = (*a_dat).getEBISBox().getEBGraph();
      pout() << "flux data for iv " << EBDebugPoint::s_ivd << ":" << endl;
      int idir = a_dat->direction();
      for (FaceIterator faceit(ivs, ebgraph, idir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
        {
          pout() << "face" << faceit() << ": ";
          for (int ivar = 0; ivar < a_dat->nComp(); ivar++)
            {
              pout() << " " << (*a_dat)(faceit(), ivar) << " ";
            }
          pout() << endl;
        }
    }
}

void printPointLDFlux(LevelData<EBFluxFAB>* a_dat)
{
  for (DataIterator dit = a_dat->dataIterator(); dit.ok(); ++dit)
    {
      printPointEBFF(&((*a_dat)[dit()]));
    }
}

void dumpEBFaceIVS(const EBFaceFAB* a_fab, const IntVectSet* a_ivs, Real a_thresh)
{
  int idir = a_fab->direction();
  const EBGraph& ebgraph = a_fab->getEBISBox().getEBGraph();
  const int ncomp = a_fab->nComp();
  for (FaceIterator faceit(*a_ivs, ebgraph, idir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      bool writeThis = true;
      //for (int ivar = 0; ivar < ncomp; ivar++)
      //  {
      //    Real data = (*a_fab)(face, ivar);
      //    if (Abs(data) >= a_thresh)
      //      {
      //        writeThis = true;
      //      }
      //  }
      if (writeThis)
        {
          const VolIndex& voflo = face.getVoF(Side::Lo);
          const VolIndex& vofhi = face.getVoF(Side::Hi);
          pout() << "face= (("
                 << voflo.gridIndex() << ", " << voflo.cellIndex() << "), ("
                 << vofhi.gridIndex() << ", " << vofhi.cellIndex() << "))" ;

          pout() << ";      data=";
          for (int ivar = 0; ivar < ncomp; ivar++)
            {
              pout() << " "
                     << setprecision(8)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << (*a_fab)(face, ivar);
            }
          pout() << endl;
        }
    }
}
void dumpEBFABIVS(const EBCellFAB* a_fab, const IntVectSet* a_ivs, Real a_thresh)
{
  const EBGraph& ebgraph = a_fab->getEBISBox().getEBGraph();
  const int ncomp = a_fab->nComp();
  for (VoFIterator vofit(*a_ivs, ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      bool writeThis = true;
//      for (int ivar = 0; ivar < ncomp; ivar++)
//        {
//          Real data = (*a_fab)(vof, ivar);
////          if (Abs(data) >= a_thresh)
////            {
////              writeThis = true;
////            }
//        }
      if (writeThis)
        {
          pout() << "vof= " << vof.gridIndex() << ", " << vof.cellIndex();

          pout() << ";      data=";
          for (int ivar = 0; ivar < ncomp; ivar++)
            {
              pout() << " "
                     << setprecision(8)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << (*a_fab)(vof, ivar);
            }
          pout() << endl;
        }
    }
}
void dumpEBFAB(const EBCellFAB* a_fab)
{
  const EBCellFAB& fab = *a_fab;
  Box box = fab.getRegion();
  IntVectSet ivs(box);
  pout() << "valid and ghost data in ebcellfab" << endl;
  dumpEBFABIVS(a_fab, &ivs);
}
void dumpEBFace(const EBFaceFAB* a_fab)
{
  dumpEBFaceThresh(a_fab, 0.0);
}

void dumpEBFlux(const EBFluxFAB* a_fab)
{
  dumpEBFluxThresh(a_fab, 0.0);
}

void dumpEBFaceThresh(const EBFaceFAB* a_fab, Real a_thresh)
{
  const EBFaceFAB& fab = *a_fab;
  Box box = fab.getCellRegion();
  IntVectSet ivs(box);
  pout() << "valid and ghost data in ebcellfab" << endl;
  dumpEBFaceIVS(a_fab, &ivs, a_thresh);
}

void dumpEBFluxThresh(const EBFluxFAB* a_fab, Real a_thresh)
{
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      pout() << "face direction = " << idir << endl;
      dumpEBFaceThresh(&((*a_fab)[idir]), a_thresh);
    }
}

void dumpEBLevelFluxThresh(const LevelData<EBFluxFAB>* a_level, Real a_thresh)
{
  const LevelData<EBFluxFAB>& lev = *a_level;
  pout() << "eblevelflux:" << endl;
  for (DataIterator dit = lev.dataIterator(); dit.ok(); ++dit)
    {
      pout() << "box = " << lev.disjointBoxLayout().get(dit()) << endl;
      dumpEBFluxThresh(&(lev[dit()]), a_thresh);
    }
}
void dumpEBLevelFlux(const LevelData<EBFluxFAB>* a_fab)
{
  dumpEBLevelFluxThresh(a_fab, 0.0);
}

void dumpEBLevelIrreg(const LevelData<EBCellFAB>* a_level)
{
  dumpEBLevelIrregThresh(a_level, 0.0);
}

void dumpEBLevel(const LevelData<EBCellFAB>* a_level)
{
  dumpEBLevelThresh(a_level, 0.0);
}

void dumpEBLevelIrregThresh(const LevelData<EBCellFAB>* a_level, Real a_thresh)
{
  const LevelData<EBCellFAB>& lev = *a_level;
  pout() << "valid irregular data in eblevel" << endl;
  for (DataIterator dit = lev.dataIterator(); dit.ok(); ++dit)
    {
      Box box = lev.disjointBoxLayout().get(dit());
      IntVectSet ivs = lev[dit()].getEBISBox().getIrregIVS(box);
      dumpEBFABIVS(&(lev[dit()]), &ivs, a_thresh);
    }
}

void dumpEBLevelAll(const LevelData<EBCellFAB>* a_level)
{
  const LevelData<EBCellFAB>& lev = *a_level;
  pout() << "valid irregular data in eblevel" << endl;
  for (DataIterator dit = lev.dataIterator(); dit.ok(); ++dit)
    {
      Box box = lev[dit()].box();
      box &= lev[dit()].getEBISBox().getDomain();
      IntVectSet ivs(box);
      dumpEBFABIVS(&(lev[dit()]), &ivs, 0.0);
    }
}
void dumpEBLevelThresh(const LevelData<EBCellFAB>* a_level, Real a_thresh)
{
  const LevelData<EBCellFAB>& lev = *a_level;
  pout() << "valid data in eblevel" << endl;
  for (DataIterator dit = lev.dataIterator(); dit.ok(); ++dit)
    {
      Box box = lev.disjointBoxLayout().get(dit());
      IntVectSet ivs(box);
      dumpEBFABIVS(&(lev[dit()]), &ivs, a_thresh);
    }
}

void dumpEBAMR(const Vector<LevelData<EBCellFAB>* >* a_amr)
{
  for (int ilev = 0; ilev < a_amr->size(); ilev++)
    {
      pout() << "ilev = " << ilev << " : " ;
      dumpEBLevel((*a_amr)[ilev]);
    }
}

void dumpEBAMRIrreg(const Vector<LevelData<EBCellFAB>* >* a_amr)
{
  for (int ilev = 0; ilev < a_amr->size(); ilev++)
    {
      pout() << "ilev = " << ilev << " : " ;
      dumpEBLevelIrreg((*a_amr)[ilev]);
    }
}


void dumpEBAMRThresh(const Vector<LevelData<EBCellFAB>* >* a_amr,     Real  a_thresh)
{
  for (int ilev = 0; ilev < a_amr->size(); ilev++)
    {
      pout() << "ilev = " << ilev << " : " ;
      dumpEBLevelThresh((*a_amr)[ilev], a_thresh);
    }
}

void dumpEBAMRIrregThresh(const Vector<LevelData<EBCellFAB>* >* a_amr, Real  a_thresh)
{
  for (int ilev = 0; ilev < a_amr->size(); ilev++)
    {
      pout() << "ilev = " << ilev << " : " ;
      dumpEBLevelIrregThresh((*a_amr)[ilev], a_thresh);
    }
}

void dumpEBFABIrreg(const EBCellFAB* a_fab)
{
  const EBCellFAB& fab = *a_fab;
  const EBGraph& ebgraph = fab.getEBISBox().getEBGraph();
  Box box = fab.getRegion();
  IntVectSet ivs= ebgraph.getIrregCells(box);

  pout() << "valid and ghost irregular cell data in ebcellfab" << endl;
  dumpEBFABIVS(a_fab, &ivs);
}

void dumpEBFABIrregGeometry(const EBCellFAB* a_fab)
{
  const EBCellFAB& fab = *a_fab;
  const EBISBox& ebisbox = fab.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  Box box = fab.getRegion();
  IntVectSet ivs= ebgraph.getIrregCells(box);

  pout() << "valid and ghost irregular geometry data in ebcellfab" << endl;
  pout() << setprecision(8)
         << setiosflags(ios::showpoint)
         << setiosflags(ios::scientific);

  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      pout() << "vof= " << vof.gridIndex() << ", " << vof.cellIndex();

      pout() << setiosflags(ios::showpos);
      pout() << "; \t volFrac=" << ebisbox.volFrac(vof);
      pout() << "; \t cntroid=" << ebisbox.centroid(vof);
      pout() << "; \t bdrArea=" << ebisbox.bndryArea(vof);
      pout() << "; \t bdrCntr=" << ebisbox.bndryCentroid(vof);
      pout() << "; \t bdrNrml=" << ebisbox.normal(vof);
      pout() << endl;
      pout() << resetiosflags(ios::showpos);

      int nFace = ebisbox.numFacePhase(vof);
      if (nFace>1)
        {
          for (int iface=0; iface<nFace; iface++)
            {
              pout() << "   face= " << iface;
              pout() << setiosflags(ios::showpos);
              pout() << "; \t bdrArea=" << ebisbox.bndryArea(vof, iface);
              pout() << "; \t bdrCntr=" << ebisbox.bndryCentroid(vof, iface);
              pout() << "; \t bdrNrml=" << ebisbox.normal(vof, iface);
              pout() << endl;
              pout() << resetiosflags(ios::showpos);
            }
        }
    }
  // pout() << resetiosflags(ios::showpoint)
  //        << resetiosflags(ios::scientific);
}

void
getMaxEBFAB(const EBCellFAB*  ldptr)
{

  const EBCellFAB& fab = *ldptr;
  for (int icomp = 0; icomp < fab.nComp(); icomp++)
    {
      Real maxval = -1.0e16;
      Real minval =  1.0e16;
      VolIndex vofmax, vofmin;
      pout() << "c = " << icomp << ", ";

      const EBGraph&   ebg = fab.getEBISBox().getEBGraph();
      const Box& grid = fab.box();
      IntVectSet ivs(grid);
      for (VoFIterator vofit(ivs, ebg); vofit.ok(); ++vofit)
        {
          Real datval = fab(vofit(), icomp);
          if (datval > maxval)
            {
              maxval = datval;
              vofmax = vofit();
            }
          if (datval < minval)
            {
              minval = datval;
              vofmin = vofit();
            }
        }
      pout() << "max=" <<  maxval << " at " << vofmax << ", ";
      pout() << "min=" <<  minval << " at " << vofmin << endl;
    }
}

void
getMaxEBLevel(const LevelData<EBCellFAB>*  ldptr)
{

  const LevelData<EBCellFAB>& ld = *ldptr;
  for (int icomp = 0; icomp < ld.nComp(); icomp++)
    {
      Real maxval = -1.0e16;
      Real minval =  1.0e16;
      VolIndex vofmax, vofmin;
      pout() << "c = " << icomp << ", ";
      const DisjointBoxLayout& dbl = ld.disjointBoxLayout();
      for (DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
        {
          const EBCellFAB& fab = ld[dit()];
          const EBGraph&   ebg = fab.getEBISBox().getEBGraph();
          const Box& grid = dbl.get(dit());
          IntVectSet ivs(grid);
          for (VoFIterator vofit(ivs, ebg); vofit.ok(); ++vofit)
            {
              Real datval = ld[dit()](vofit(), icomp);
              if (datval > maxval)
                {
                  maxval = datval;
                  vofmax = vofit();
                }
              if (datval < minval)
                {
                  minval = datval;
                  vofmin = vofit();
                }
            }
        }
      pout() << "max=" <<  maxval << " at " << vofmax << ", ";
      pout() << "min=" <<  minval << " at " << vofmin << endl;
    }
}
void
dumpLDEBCF(const LevelData<EBCellFAB>*  ldptr)
{
  const LevelData<EBCellFAB>& ld = *ldptr;
  for (DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
    {
      dumpEBFAB(&(ld[dit()]));
    }
}

void
dumpEBLDDBL(const LevelData<EBCellFAB>*  memLDF_Ptr)
{
  const DisjointBoxLayout& dbl = memLDF_Ptr->disjointBoxLayout();
  dumpDBL(&dbl);
}

void dumpLDBIVF(const LayoutData< BaseIVFAB<Real> >* a_ldptr)
{
  const LayoutData< BaseIVFAB<Real> >& ld = *a_ldptr;
  for (DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
    {
      dumpIVFAB(&ld[dit()]);
    }
}

void dumpIVFAB(const BaseIVFAB<Real>* a_vectPtr)
{
  const BaseIVFAB<Real>& ivfab = *a_vectPtr;
  const int ncomp = ivfab.nComp();
  const EBGraph& ebgraph = ivfab.getEBGraph();
  const IntVectSet& ivs = ivfab.getIVS();

  pout() << "data in base ivfab" << endl;
  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      pout() << "vof= " << vof.gridIndex() << ", " << vof.cellIndex();
      //RealVect centroid = ebgraph.centroid(vof);
      //pout() << ";      cent=";
      //for (int idir = 0; idir < SpaceDim; idir++)
      //  {
      //    pout() << " "
      //           << setprecision(8)
      //           << setiosflags(ios::showpoint)
      //           << setiosflags(ios::scientific)
      //           << centroid[idir];
      //  }
      pout() << ";      data=";
      for (int ivar = 0; ivar < ncomp; ivar++)
        {
          pout() << " "
                 << setprecision(8)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << ivfab(vof, ivar);
        }
      pout() << endl;
    }
}
void dumpIFFAB(const BaseIFFAB<Real>* a_vectPtr)
{
  const BaseIFFAB<Real>& iffab = *a_vectPtr;
  const int ncomp = iffab.nComp();
  const int dir = iffab.direction();
  const EBGraph& ebgraph = iffab.getEBGraph();
  const IntVectSet& ivs = iffab.getIVS();
  FaceIterator faceit(ivs, ebgraph, dir,
                      FaceStop::SurroundingWithBoundary);
  pout() << "data in base iffab" << endl;
  for (faceit.reset();faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();
      pout() << "face= ";
      for (SideIterator sit; sit.ok(); ++sit)
        {
          pout() << " " << face.gridIndex(sit());
        }
      for (int ivar = 0; ivar < ncomp; ivar++)
        {
          pout() << " "
                 << setprecision(8)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << iffab(face, ivar);
        }
      pout() << endl;
    }
}
void dumpVVoF(const Vector<VolIndex>* a_vectPtr)
{
  if (a_vectPtr == NULL) return;
  Vector<VolIndex> vect = *a_vectPtr;
  pout() << "vector of vofs contains" << endl;
  for (int iveco = 0; iveco < vect.size(); iveco++)
    {
      const VolIndex& vof = vect[iveco];
      pout() <<  vof.gridIndex() << "   "  << vof.cellIndex() << "  ";;
    }
  pout() << endl;
}

void dumpVFace(const Vector<FaceIndex>* a_vectPtr)
{
  if (a_vectPtr == NULL) return;
  Vector<FaceIndex> vect = *a_vectPtr;
  pout() << "vector of faces contains" << endl;
  for (int iveco = 0; iveco < vect.size(); iveco++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          VolIndex vof = vect[iveco].getVoF(sit());
          pout() <<  vof.gridIndex() << "   " << vof.cellIndex() << "  ";;
        }
      if (vect[iveco].isBoundary())
        {
          pout() <<  "on boundary";
        }
      pout() << endl;
    }
}

void dumpFace(const FaceIndex* a_vectPtr)
{
  if (a_vectPtr == NULL) return;
  const FaceIndex& face = *a_vectPtr;
  pout() << "faces contains" << endl;
  for (SideIterator sit; sit.ok(); ++sit)
    {
      VolIndex vof = face.getVoF(sit());
      pout() <<  vof.gridIndex() << "   " << vof.cellIndex() << "  ";
    }
  pout() << endl;
  if (face.isBoundary())
    {
      pout() <<  " face is on on boundary";
      pout() << endl;
    }
}

void dumpFaceSten(const FaceStencil* a_vectPtr)
{
  if (a_vectPtr == NULL) return;
  const FaceStencil& sten = *a_vectPtr;
  pout() << "stencil contains:" << endl;
  for (int isten = 0; isten <  sten.size(); isten++)
    {
      const FaceIndex& face = sten.face(isten);
      const Real& weight = sten.weight(isten);
      pout() << "face: ";
      for (SideIterator sit; sit.ok(); ++sit)
        {
          VolIndex vof = face.getVoF(sit());
          pout() <<  vof.gridIndex() << "   " << vof.cellIndex() << "  ";;
        }
      pout() << "weight: " << weight;
      pout() << endl;
    }
}

void dumpVoFSten(const VoFStencil* a_vectPtr)
{
  if (a_vectPtr == NULL) return;
  const VoFStencil& sten = *a_vectPtr;
  pout() << "stencil contains:" << endl;
  for (int isten = 0; isten <  sten.size(); isten++)
    {
      const VolIndex& vof = sten.vof(isten);
      const Real& weight = sten.weight(isten);
      pout() << "vof: (";
      pout() <<  vof.gridIndex() << ", " << vof.cellIndex() << ") ";;
      pout() << "var: " << sten.variable(isten) << " ";
      pout() << "weight: " << weight;
      pout() << endl;
    }
}
#include "NamespaceFooter.H"
