#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stdio.h>

#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "IndexTM.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "PolyGeom.H"
#include "EBAMRIO.H"

#include "EBMenagerieUtils.H"

#include "UsingNamespace.H"

void makeLayout(DisjointBoxLayout& a_dbl,
                const Box&         a_domain)
{
  ParmParse pp;

  int maxsize;
  Vector<Box> vbox(1,a_domain);

  pp.get("maxboxsize",maxsize);

  domainSplit(a_domain,vbox,maxsize);

  Vector<int> procAssign;
  LoadBalance(procAssign,vbox);

  a_dbl.define(vbox,procAssign);
}

void makeEBISL(EBISLayout&              a_ebisl,
               const DisjointBoxLayout& a_grids,
               const Box&               a_domain,
               const int&               a_nghost)
{
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());

  ebisPtr->fillEBISLayout(a_ebisl,a_grids,a_domain,a_nghost);
}

void fillData(LevelData<EBCellFAB>& a_level,
              const RealVect&       a_origin,
              const Real&           a_dx,
              const BaseIF&         a_implicit)
{
  // Set the data in each VoF to the value of the implicit function at the
  // cell center of the cell containing the VoF
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
  {
    EBCellFAB& ebcell = a_level[dit()];
    const EBISBox& ebisbox = ebcell.getEBISBox();
    FArrayBox& data = ebcell.getFArrayBox();
    const Box& box = ebcell.box();

    for (BoxIterator bit(box); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      RealVect coord(a_origin);
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        coord[idir] += a_dx * (iv[idir] + 0.5);
      }

      Real value = a_implicit.value(coord);

      Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
      int size = vofs.size();

      if (size == 0)
      {
        data(iv,0) = value;
      }
      else
      {
        for (int i = 0; i < size; i++)
        {
          ebcell(vofs[i],0) = value;
        }
      }
    }
  }
}

void fillData(LevelData<EBCellFAB>& a_level,
              const RealVect&       a_origin,
              const RealVect&       a_dx,
              const BaseIF&         a_implicit)
{
  // Set the data in each VoF to the value of the implicit function at the
  // cell center of the cell containing the VoF
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB& ebcell = a_level[dit()];
      const EBISBox& ebisbox = ebcell.getEBISBox();
      FArrayBox& data = ebcell.getFArrayBox();
      const Box& box = ebcell.box();

      for (BoxIterator bit(box); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();

          RealVect coord(a_origin);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              coord[idir] += a_dx[idir] * (iv[idir] + 0.5);
            }

          Real value = a_implicit.value(coord);

          Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
          int size = vofs.size();

          if (size == 0)
            {
              data(iv,0) = value;
            }
          else
            {
              for (int i = 0; i < size; i++)
                {
                  ebcell(vofs[i],0) = value;
                }
            }
        }
    }
}

void fillData(LevelData<EBCellFAB>& a_data,
              const RealVect&       a_origin,
              const RealVect&       a_dx,
              const BaseIF&         a_implicit,
              const IntVect&        a_ghost)
{
  // Set the data in each VoF to the value of the implicit function at the
  // cell center of the cell containing the VoF
  const DisjointBoxLayout& dbl = a_data.disjointBoxLayout();
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      const Box& dblBox = dbl.get(dit());
      EBCellFAB& ebcell = a_data[dit()];
      const EBISBox& ebisbox = ebcell.getEBISBox();
      FArrayBox& data = ebcell.getFArrayBox();
      Box ghostedBox = dblBox;
      ghostedBox.grow(a_ghost);

      for (BoxIterator bit(ghostedBox); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();

          RealVect coord(a_origin);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              coord[idir] += a_dx[idir] * (iv[idir] + 0.5);
            }

          Real value = a_implicit.value(coord);

          if (!dblBox.contains(iv))
            {//do this for domain ghosts (interior ghosts get over-written by an exchange)
              data(iv,0) = value;
            }
          else
            {
              Vector<VolIndex> vofs = ebisbox.getVoFs(iv);
              int size = vofs.size();
              if (size == 0)
                {
                  data(iv,0) = value;
                }
              else
                {
                  for (int i = 0; i < size; i++)
                    {
                      ebcell(vofs[i],0) = value;
                    }
                }
            }
        }
    }
  a_data.exchange();
}

void makeLayouts(Vector<DisjointBoxLayout>& a_dbl,
                 const Vector<int>&         a_refRatio,
                 const Vector<RealVect>&    a_dx,
                 const RealVect&            a_origin,
                 const Vector<Box>&         a_domain)
{
  pout() << "Begin making layouts " << endl;
  ParmParse pp;

  int maxBoxSize;
  pp.get("maxboxsize",maxBoxSize);

  Real fillRatio;
  pp.get("fill_ratio",fillRatio);

  int blockFactor;
  pp.get("block_factor",blockFactor);

  int bufferSize;
  pp.get("buffer_size",bufferSize);

  int numLevels = a_domain.size();
  a_dbl.resize(numLevels);

//   for (int ilev = 0;ilev<numLevels;ilev++)
//     {
//       Vector<Box> vbox(1,a_domain[ilev]);

//       domainSplit(a_domain[ilev],vbox,maxBoxSize);

//       Vector<int> procAssign;
//       LoadBalance(procAssign,vbox);

//       a_dbl[ilev].define(vbox,procAssign);
//     }

  ///////

  // first, construct tags
  Vector<IntVectSet> tagsVect(numLevels);
  Vector<Vector<Box> > boxes;
  Vector<Vector<Box> > vectBoxes(numLevels);

  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      domainSplit(a_domain[ilev], vectBoxes[ilev],
                  maxBoxSize);
    }

  if (numLevels != 1)
    {
      // do tagging
      tagCells(tagsVect,a_refRatio,a_dx,a_origin,a_domain);

      BRMeshRefine meshrefine(a_domain[0], a_refRatio,
                              fillRatio,   blockFactor,
                              bufferSize,  maxBoxSize);

      int top_level = numLevels - 1;

      // for now, just assume that all levels get regridded at the same time
      int lbase = 0;
      int new_finest_level = meshrefine.regrid(boxes,
                                               tagsVect,
                                               lbase,
                                               top_level,
                                               vectBoxes);

      CH_assert(new_finest_level == numLevels-1);

      for (int ilev = 0; ilev < numLevels; ilev++)
        {
          Vector<int> procAssign;
          LoadBalance(procAssign,boxes[ilev]);

          a_dbl[ilev].define(boxes[ilev],procAssign);
        }
    }
  else
    {
      Vector<int> procAssign;
      LoadBalance(procAssign,vectBoxes[0]);

      a_dbl[0].define(vectBoxes[0],procAssign);
    }
}

void tagCells(Vector<IntVectSet>&     a_tags,
              const Vector<int>&      a_refRatio,
              const Vector<RealVect>& a_dx,
              const RealVect&         a_origin,
              const Vector<Box>&      a_domain)
{
  int numLevels = a_tags.size();

  for (int lev=0; lev<numLevels-1; lev++)//only tag on coarser levels
    {
      IntVectSet& levelTags = a_tags[lev];

      tagCellsLevel(levelTags, a_domain[lev],a_dx[lev],a_origin);

      //fix the no-tag case (refine based on coarser)
      const int numPts = levelTags.numPts();
      if (numPts == 0)
        {
          if (lev==0)
            {//tag based on an arbitrary box
              Box tagBox = a_domain[lev];
              tagBox.coarsen(a_refRatio[lev]);
              levelTags |= tagBox;
              MayDay::Warning("no tags for coarsest level::arbitrary refinement specified");
            }
          else
            {//tag based on coarser tags
              MayDay::Warning("tagging level based on coarser tags");
              IntVectSet refTags = a_tags[lev-1];
              refTags.refine(a_refRatio[lev]);
              levelTags |= refTags;
            }
        }
    }
}

void tagCellsLevel(IntVectSet&     a_tags,
                   const Box&      a_domain,
                   const RealVect& a_dx,
                   const RealVect& a_origin)
{
  // parse input file
  ParmParse pp;

  //tag boxes based on low and high points in physical coordinates specified in the input file
  int numBox = 2;
  Vector<const char*> nameLo(numBox);
  Vector<const char*> nameHi(numBox);
  nameLo[0] = "tag_box1_lo";
  nameHi[0] = "tag_box1_hi";
  nameLo[1] = "tag_box2_lo";
  nameHi[1] = "tag_box2_hi";

  for (int iBox = 0; iBox < numBox; iBox++)
    {
      Vector<Real> lo(SpaceDim);
      Vector<Real> hi(SpaceDim);

      pp.getarr(nameLo[iBox],lo,0,SpaceDim);
      pp.getarr(nameHi[iBox],hi,0,SpaceDim);

      //figure out the lo and hi iv's of the tag box
      IntVect ivLo,ivHi;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Real xivLo = (lo[idir] - a_origin[idir])/a_dx[idir];
          Real xivHi = (hi[idir] - a_origin[idir])/a_dx[idir];
          ivLo[idir] = int(xivLo - 0.5);
          ivHi[idir] = int(xivHi - 0.5);
        }

      //make the box and add it to a_tags
      Box tagBox = Box(ivLo,ivHi);
      tagBox &= a_domain;
      a_tags |= tagBox;
    }
}

void makeEBISL(Vector<EBISLayout>&              a_ebisl,
               const Vector<DisjointBoxLayout>& a_grids,
               const Vector<Box>&               a_domain,
               const int&                       a_nghost)
{
  pout() << "Begin making EBISLs " << endl;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());

  int numLevels = a_domain.size();
  a_ebisl.resize(numLevels);
  for (int ilev = 0;ilev<numLevels;ilev++)
    {
      ebisPtr->fillEBISLayout(a_ebisl[ilev],a_grids[ilev],a_domain[ilev],a_nghost);
    }
}


void createEBDistributionFiles()
{
#ifdef CH_USE_HDF5
  pout() << "Begin createEBDistributionFiles " << endl;

  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  CH_assert(ebisPtr->isDefined());

  ProblemDomain domain;
  Vector<LevelData<EBCellFAB>*> allBoxes, IrregBoxes;
  int numLevels = ebisPtr->numLevels();

  domain = ebisPtr->getBox(numLevels-1);
  Real dx = ebisPtr->dx(numLevels-1);
  Vector<int> refRatio(numLevels, 2);
  Vector<DisjointBoxLayout> grids;
  Vector<std::string> names(1, std::string("procID"));

  for (int i=numLevels-1; i>=0; i--)
    {
      DisjointBoxLayout dbl = ebisPtr->levelGrids(i);
      grids.push_back(dbl);
      EBISLayout ebisl;
      ebisPtr->fillEBISLayout(ebisl, dbl, ebisPtr->getBox(i), 1);
      EBCellFactory factory(ebisl);
      allBoxes.push_back(new LevelData<EBCellFAB>(dbl, 1, IntVect::Zero, factory));
      LevelData<EBCellFAB>& ld = *(allBoxes.back());
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          int p = procID();
          ld[dit].setVal(p);
        }
    }

  writeEBHDF5(std::string("EBIndexSpace.hdf5"),
              grids,
              allBoxes,
              names,
              domain.domainBox(),
              dx,
              0,
              0,
              refRatio,
              numLevels,
              false,
              Vector<Real>(1,0));

  for (int i=0; i<numLevels; ++i)
    {
      delete allBoxes[i];
    }
#endif
}
