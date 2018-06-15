#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeSetOperations.cpp
// petermc, 12 June 2003

#include "NodeSetOperations.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "HalveIntF_F.H"
#include "BaseFabIntPlusF_F.H"
#include "MaskValueF_F.H"
#include "NamespaceHeader.H"
using std::cerr;
using std::endl;

// ---------------------------------------------------------
void
interiorNodes(IntVectSet& a_ivs,
              const Box& a_base_domain,
              const DisjointBoxLayout& a_boxLayout,
              const Box& a_box)
{
  ProblemDomain probdomain(a_base_domain);
  interiorNodes(a_ivs, probdomain, a_boxLayout, a_box);
}


// ---------------------------------------------------------
void
interiorNodes(IntVectSet& a_ivs,
              const ProblemDomain& a_base_domain,
              const DisjointBoxLayout& a_boxLayout,
              const Box& a_box)
{
  // define bgrow as the cells of a_box plus one ghost layer.
  // IntVectSet bgrow(grow(a_box, 1));
  Box bgrow(a_box);
  bgrow.grow(1);
  // intersect bgrow with the base domain.
  // (no effect in directions where base domain is periodic)
  bgrow &= a_base_domain;
  // shell is bgrow with all boxes at this level removed.
  IntVectSet shell(bgrow);

  Box periodicTestBox(a_base_domain.domainBox());
  if (a_base_domain.isPeriodic())
    for (int idir = 0; idir < SpaceDim; idir++)
      if (a_base_domain.isPeriodic(idir))
        periodicTestBox.grow(idir, -1);

  IntVect shiftMult = a_base_domain.domainBox().size();

  // (use LayoutIterator to get _all_ boxes, even those in other procs)
  for (LayoutIterator lit = a_boxLayout.layoutIterator(); lit.ok(); ++lit)
    {
      Box shiftedBox = a_boxLayout.get(lit);
      shell -= shiftedBox;

      // only do this IF we're periodic _and_ both boxes
      // adjoin the domain box boundary somewhere
      if (a_base_domain.isPeriodic() &&
          !periodicTestBox.contains(bgrow) &&
          !periodicTestBox.contains(shiftedBox))
        {
          ShiftIterator shiftIt = a_base_domain.shiftIterator();
          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
            {
              IntVect shiftVect(shiftMult*shiftIt());
              shiftedBox.shift(shiftVect);
              shell -= shiftedBox;
              shiftedBox.shift(-shiftVect);
            }
        }
    }
  // define interiorPlus as the cells removed from bgrow to yield shell.
  // const IntVectSet interiorPlus(bgrow - shell);
  IntVectSet interiorPlus(bgrow);
  interiorPlus -= shell;
  // initialize interior as interiorPlus.
  a_ivs.define(interiorPlus);
  Box offsetBox(IntVect::Zero, IntVect::Unit);
  // intersect a_ivs with interiorPlus with offsets
  BoxIterator offset(offsetBox);
  for (offset.begin(); offset.ok(); ++offset)
    {
      IntVectSet interiorPlus_off(interiorPlus);
      interiorPlus_off.shift(offset());
      a_ivs &= interiorPlus_off;
    }
}


// ---------------------------------------------------------
void
interiorBoundaryNodes(LayoutData< Vector<IntVectSet> >& a_IVSV,
                      const DisjointBoxLayout& a_boxes,
                      const Box& a_domain)
{
  ProblemDomain probdomain(a_domain);
  interiorBoundaryNodes(a_IVSV, a_boxes, probdomain);
}


// ---------------------------------------------------------
void
interiorBoundaryNodes(LayoutData< Vector<IntVectSet> >& a_IVSV,
                      const DisjointBoxLayout& a_boxes,
                      const ProblemDomain& a_domain)
{
  { // Revise DenseIntVectSet threshold upward if necessary.
    IntVectSet sampleIVS;
    int maxd = IntVectSet::s_maxDense;
    for (LayoutIterator lit = a_boxes.layoutIterator(); lit.ok(); ++lit)
      {
        const Box& bxCells = a_boxes.get(lit);
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            Box face = grow(surroundingNodes(bdryLo(bxCells, idir, 1)),
                            IntVect::Unit - BASISV(idir));
            long npts = face.numPts();
            if (npts > maxd)
              {
                maxd = npts;
                sampleIVS.setMaxDense(maxd);
              }
          }
      }
  }

  // boxesNodes contains the indices of the nodes of the boxes in a_boxes,
  // but they're represented as CELLs.
  // DataIterator dit = a_boxes.dataIterator();
  // BoxLayout boxesNodes;
  // boxesNodes.deepCopy(a_boxes);
  // for (dit.begin(); dit.ok(); ++dit)
  // {
  // Box bx = a_boxes[dit];
  // for (int idir = 0; idir < SpaceDim; idir++)
  // bx.growHi(idir, 1);
  // boxesNodes[dit] = bx;
  // }
  // boxesNodes.close();

  // periodicTestBox is the domain box with one layer of cells removed
  // from each end in periodic directions.
  Box periodicTestBox(a_domain.domainBox());
  if (a_domain.isPeriodic())
    for (int idir = 0; idir < SpaceDim; idir++)
      if (a_domain.isPeriodic(idir))
        periodicTestBox.grow(idir, -1);

  // We use this to move things around in the periodic case.
  IntVect shiftMult = a_domain.domainBox().size();

  a_IVSV.define(a_boxes);
  for (DataIterator dit = a_boxes.dataIterator(); dit.ok(); ++dit)
    {
      Vector<IntVectSet>& IVSV = a_IVSV[dit];

      const Box& bxCells = a_boxes.get(dit);

      // bxNodes is CELL-centered, contains indices of nodes of bxCells
      // const Box& bxNodes = boxesNodes.get(dit);

      // on each side, find interior boundary nodes of bxCells
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide mySide = sit();

              Box faceboxOriginal = adjCellBox(bxCells, idir, mySide, 1);

              Box faceboxIn = faceboxOriginal & a_domain;
              // empty faceboxIn means physical boundary
              IntVectSet thisIVS(faceboxIn);
              Box facegrow(faceboxIn);
              if (! faceboxIn.isEmpty())
                {
                  // facegrow (CELL-centered) = faceboxIn expanded by:
                  // one layer in idir, side INTO the box;
                  // one layer in all non-idir directions;
                  // and then truncated to physical domain.
                  facegrow.growDir(idir, flip(mySide), 1);
                  for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                    if (idirOther != idir) facegrow.grow(idirOther, 1);

                  facegrow &= a_domain;

                  // box and facegrow for idir==0, mySide==Side::Hi
                  //
                  //              +---+
                  //      +-------+-+ |
                  //      |       | | |
                  //      |       | | |
                  //      |       | | |
                  //      +-------+-+ |
                  //              +---+
                  //

                  // shell is facegrow with all boxes on this level removed.
                  IntVectSet shell(facegrow);
                  // (use LayoutIterator to get _all_ boxes, even those in other procs)

                  for (LayoutIterator lit = a_boxes.layoutIterator(); lit.ok(); ++lit)
                    {
                      Box shiftedBox = a_boxes.get(lit);
                      shell -= shiftedBox;
                      // added by petermc, 4 Feb 2003
                      // only do this IF we're periodic _and_ both boxes
                      // adjoin the domain box boundary somewhere
                      if (a_domain.isPeriodic() &&
                          !periodicTestBox.contains(facegrow) &&
                          !periodicTestBox.contains(shiftedBox))
                        {
                          ShiftIterator shiftIt = a_domain.shiftIterator();
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                            {
                              IntVect shiftVect(shiftMult*shiftIt());
                              shiftedBox.shift(shiftVect);
                              shell -= shiftedBox;
                              shiftedBox.shift(-shiftVect);
                            }
                        }
                    }

                  // define interiorPlus as the cells removed from facegrow
                  // to yield shell.
                  // const IntVectSet interiorPlus(facegrow - shell);
                  IntVectSet interiorPlus(facegrow);
                  interiorPlus -= shell;
                  // initialize interior as interiorPlus.
                  thisIVS = interiorPlus;

                  Box offsetBox(IntVect::Zero, IntVect::Unit);
                  // intersect thisIVS with interiorPlus with offsets
                  for (BoxIterator offset(offsetBox); offset.ok(); ++offset)
                    {
                      IntVectSet interiorPlus_off(interiorPlus);
                      interiorPlus_off.shift(offset());
                      thisIVS &= interiorPlus_off;
                    }
                }

              if (! thisIVS.isEmpty() )
                {
                  // IntVectSet IVS(bxNodes);
                  // petermc, 16 July 2003:  bxNodes is too big and forces
                  // a TreeIntVectSet.
                  IntVectSet IVS(facegrow);
                  IVS &= thisIVS;
                  // remove nodes that are already present in other components
                  for (int lcomp = 0; lcomp < IVSV.size(); lcomp++)
                    IVS -= IVSV[lcomp];
                  // add the new nodes as a new component of IVSV
                  if (! IVS.isEmpty() ) IVSV.push_back(IVS);
                }
            }
        }
    }
}


// ---------------------------------------------------------
void
interiorBoundaryNodes(LayoutData< Vector<IntVectSet> >& a_IVSV,
                      const DisjointBoxLayout& a_dest,
                      const DisjointBoxLayout& a_src,
                      const Box& a_domain)
{
  ProblemDomain probdomain(a_domain);
  interiorBoundaryNodes(a_IVSV, a_dest, a_src, probdomain);
}


// ---------------------------------------------------------
void
interiorBoundaryNodes(LayoutData< Vector<IntVectSet> >& a_IVSV,
                      const DisjointBoxLayout& a_dest,
                      const DisjointBoxLayout& a_src,
                      const ProblemDomain& a_domain)
{
  { // Revise DenseIntVectSet threshold upward if necessary.
    IntVectSet sampleIVS;
    int maxd = IntVectSet::s_maxDense;
    for (LayoutIterator lit = a_src.layoutIterator(); lit.ok(); ++lit)
      {
        const Box& bxCells = a_src.get(lit);
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            // Box face = grow(surroundingNodes(bdryLo(bxCells, idir, 1)),
            // IntVect::Unit - BASISV(idir));
            // NEW, 7 Dec 2005
            Box face = surroundingNodes(bdryLo(bxCells, idir, 1));
            long npts = face.numPts();
            if (npts > maxd)
              {
                maxd = npts;
                sampleIVS.setMaxDense(maxd);
              }
          }
      }
  }

  // destNodeIndices contains the indices of the NODEs of the boxes in a_dest,
  // but they're represented as CELLs because that's what IntVectSet uses.
  BoxLayout destNodeIndices;
  destNodeIndices.deepCopy(a_dest);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      destNodeIndices.growSide(idir, 1, Side::Hi);
    }
  destNodeIndices.close();

  Box periodicTestBox(a_domain.domainBox());
  if (a_domain.isPeriodic())
    for (int idirPeriodic = 0; idirPeriodic < SpaceDim; idirPeriodic++)
      if (a_domain.isPeriodic(idirPeriodic))
        periodicTestBox.grow(idirPeriodic, -1);

  IntVect shiftMult = a_domain.domainBox().size();

  // destSrcIntersec[dit] is true iff the current source box
  // shares nodes with a_dest[dit].
  LayoutData<bool> destSrcIntersect(a_dest);

  // cout << "interiorBoundaryNodes:  from " << a_src.size()
  //      << " boxes to " << a_dest.size() << " boxes"
  //      << " on " << a_domain.longside() << endl;
  // Define a_IVSV to have a Vector<IntVectSet> for each box of a_dest.
  a_IVSV.define(a_dest);
  for (LayoutIterator litSrc = a_src.layoutIterator(); litSrc.ok(); ++litSrc)
    {
      const Box& srcbox = a_src.get(litSrc);

      // srcNodesBox is CELL-centered, contains indices of NODEs of srcbox
      Box srcNodesBox = srcbox;
      for (int idir = 0; idir < SpaceDim; idir++)
        srcNodesBox.growHi(idir, 1);

      for (DataIterator dit = a_dest.dataIterator(); dit.ok(); ++dit)
        {
          destSrcIntersect[dit] = srcNodesBox.intersects(destNodeIndices.get(dit()));
        }

      // bool destSrcIntersect[dit]:
      // Do destNodeIndicesBox[dit] and srcNodesBox intersect?

      // On each face, find interior boundary nodes of srcbox.
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide mySide = sit();

              Box faceboxOriginal = adjCellBox(srcbox, idir, mySide, 1);

              Box faceboxIn = faceboxOriginal & a_domain;
              // empty faceboxIn means physical boundary
              IntVectSet srcIVS(faceboxIn);
              Box facegrow(faceboxIn);
              if (! faceboxIn.isEmpty())
                {
                  // facegrow = faceboxIn expanded by:
                  // one layer in idir, side INTO the box;
                  // one layer in all non-idir directions;
                  // and then truncated to physical domain.
                  facegrow.growDir(idir, flip(mySide), 1);
                  for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                    if (idirOther != idir) facegrow.grow(idirOther, 1);

                  facegrow &= a_domain;

                  // box and facegrow for idir==0, mySide==Side::Hi
                  //
                  //              +---+
                  //      +-------+-+ |
                  //      |       | | |
                  //      |       | | |
                  //      |       | | |
                  //      +-------+-+ |
                  //              +---+
                  //

                  // shell is facegrow with all boxes on this level removed.
                  IntVectSet shell(facegrow);
                  // (use LayoutIterator to get _all_ boxes, even those in other procs)
                  for (LayoutIterator litSrc2 = a_src.layoutIterator(); litSrc2.ok(); ++litSrc2)
                    {
                      Box shiftedBox = a_src.get(litSrc2);
                      shell -= shiftedBox;
                      // added by petermc, 4 Feb 2003
                      // only do this IF we're periodic _and_ both boxes
                      // adjoin the domain box boundary somewhere
                      if (a_domain.isPeriodic() &&
                          !periodicTestBox.contains(faceboxIn) &&
                          !periodicTestBox.contains(shiftedBox))
                        {
                          ShiftIterator shiftIt = a_domain.shiftIterator();
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                            {
                              IntVect shiftVect(shiftMult*shiftIt());
                              shiftedBox.shift(shiftVect);
                              shell -= shiftedBox;
                              shiftedBox.shift(-shiftVect);
                            }
                        }

                    }

                  // define interiorPlus as the cells removed from facegrow
                  // to yield shell.
                  // const IntVectSet interiorPlus(facegrow - shell);
                  IntVectSet interiorPlus(facegrow);
                  interiorPlus -= shell;
                  // initialize interior as interiorPlus.
                  srcIVS = interiorPlus;

                  Box offsetBox(IntVect::Zero, IntVect::Unit);
                  // intersect srcIVS with interiorPlus with offsets
                  for (BoxIterator offset(offsetBox); offset.ok(); ++offset)
                    {
                      IntVectSet interiorPlus_off(interiorPlus);
                      interiorPlus_off.shift(offset());
                      srcIVS &= interiorPlus_off;
                    }
                }

              if (! srcIVS.isEmpty() )
                {
                  for (DataIterator dit = a_dest.dataIterator(); dit.ok(); ++dit)
                    {
                      if (destSrcIntersect[dit])
                        // get nodes of destination box
                        {
                          Box destNodeIndicesBox(destNodeIndices.get(dit));
                          destNodeIndicesBox &= srcNodesBox;
                          if (! destNodeIndicesBox.isEmpty() )
                            {
                              // petermc, 16 Jul 2003:  destNodeIndicesBox was
                              // too big and forced TreeIntVectSet.
                              // Replace these two lines with reverse order.
                              // IntVectSet IVS(destNodeIndicesBox);
                              // IVS &= srcIVS;
                              IntVectSet IVS(srcIVS);
                              IVS &= destNodeIndicesBox;
                              // IVS: nodes of destination box that are also
                              // interior boundary nodes of source box
                              Vector<IntVectSet>& IVSV = a_IVSV[dit];
                              // remove nodes that are already present in other components
                              for (int lcomp = 0; lcomp < IVSV.size(); lcomp++)
                                IVS -= IVSV[lcomp];
                              // add the new nodes as a new component
                              if (! IVS.isEmpty() ) IVSV.push_back(IVS);
                            }
                        }
                    }
                }
            }
        }
    }
}


// ---------------------------------------------------------
void
exteriorBoundaryNodes(LayoutData< Vector<IntVectSet> >& a_exterior,
                      const LayoutData< Vector<IntVectSet> >& a_interior,
                      const DisjointBoxLayout& a_boxes)
{
  { // Revise DenseIntVectSet threshold upward if necessary.
    IntVectSet sampleIVS;
    int maxd = IntVectSet::s_maxDense;
    for (LayoutIterator lit = a_boxes.layoutIterator(); lit.ok(); ++lit)
      {
        const Box& bxCells = a_boxes.get(lit);
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            Box faceNodes = surroundingNodes(bxCells);
            faceNodes.setBig(idir, faceNodes.smallEnd(idir));
            long npts = faceNodes.numPts();
            if (npts > maxd)
              {
                maxd = npts;
                sampleIVS.setMaxDense(maxd);
              }
          }
      }
  }

  a_exterior.define(a_boxes);
  for (DataIterator dit = a_boxes.dataIterator(); dit.ok(); ++dit)
    {
      Box box = a_boxes.get(dit);
      Vector<IntVectSet> ibn = a_interior[dit];
      Vector<IntVectSet>& ebn = a_exterior[dit];
      // on each side, find exterior boundary nodes
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide mySide = sit();

              // faceNodes will contain the cells with indices
              // that are the same as the nodes on this face.
              Box faceNodes = adjCellBox(box, idir, mySide, 1);
              if (mySide == Side::Lo)
                {
                  // On low side, shift up by unit vector in direction idir,
                  // in order to get indices of the correct nodes.
                  faceNodes.shift(idir, 1);
                  // (On high side, the node indices are the same as
                  // the cell indices.)
                }
              // Expand by one in all directions other than idir.
              for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                if (idirOther != idir) faceNodes.growHi(idirOther, 1);

              // IVS will be the face nodes without interior boundary nodes.
              // What's left will be the exterior boundary nodes.
              IntVectSet IVS(faceNodes);
              // remove the interior boundary nodes from IVS.
              for (int lcomp = 0; lcomp < ibn.size(); lcomp++)
                IVS -= ibn[lcomp];
              // add remaining nodes as a new component of
              // exterior boundary nodes.
              if (! IVS.isEmpty() ) ebn.push_back(IVS);
            }
        }
    }
}


// ---------------------------------------------------------
void
zeroBoundaryNodes(BoxLayoutData<NodeFArrayBox>& a_dest,
                  const LayoutData< Vector<IntVectSet> >& a_IVSV)
{
  int numComps = a_dest.nComp();
  for (DataIterator dit = a_dest.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& destFab = a_dest[dit].getFab();
      Vector<IntVectSet> IVS = a_IVSV[dit];
      for (int lcomp = 0; lcomp < IVS.size(); lcomp++)
        for (IVSIterator it(IVS[lcomp]); it.ok(); ++it)
          {
            IntVect iv = it();
            for (int ivar = 0; ivar < numComps; ivar++)
              destFab(iv, ivar) = 0.;
          }
    }
}


// ---------------------------------------------------------
void
copyInteriorNodes(LevelData<NodeFArrayBox>& a_dest,
                  const LevelData<NodeFArrayBox>& a_src,
                  const LayoutData< Vector<IntVectSet> >& a_IVSV)
{
  const DisjointBoxLayout& destBoxes = a_dest.getBoxes();
  const DisjointBoxLayout& srcBoxes = a_src.getBoxes();
  const DisjointBoxLayout& ivsBoxes = (DisjointBoxLayout&) a_IVSV.boxLayout();
  CH_assert(destBoxes.numBoxes(procID()) == ivsBoxes.numBoxes(procID()));

  // (1) define psi on destBoxes
  // (2) copy a_src to psi;
  // (3) for each box in destBoxes
  //         for each interior node of srcBoxes in it,
  //             copy psi to a_dest there.

  // (1) define psi on destBoxes
  int numComps = a_src.nComp();
  BoxLayoutData<NodeFArrayBox> psi(destBoxes, numComps);

  // (2) copy a_src to psi;
  // petermc, 13 May 2003:  change Copier ghosts from Unit to Zero.
  Copier myCopier(srcBoxes, destBoxes, IntVect::Zero);
  a_src.copyTo(a_src.interval(), psi, psi.interval(), myCopier);

  // (3) for each box in destBoxes
  for (DataIterator dit = destBoxes.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& psiFab = psi[dit].getFab();
      FArrayBox& destFab = a_dest[dit].getFab();

      // copy psi to a_dest there.

      // destNodeIndices:  indices of _all_ nodes in destBoxes[dit]
      Box destNodeIndices = destBoxes.get(dit);
      for (int dir = 0; dir < SpaceDim; dir++)
        destNodeIndices.growHi(dir, 1);

      // (1) interior of each source box

      for (LayoutIterator lit = srcBoxes.layoutIterator(); lit.ok(); ++lit)
        {
          const Box thisBox = srcBoxes.get(lit);
          Box inner(thisBox);
          for (int dir = 0; dir < SpaceDim; dir++)
            inner.growLo(dir, -1);

          inner &= destNodeIndices;

          // destFab.copy(psiFab, inner);  // incompatible box types
          for (BoxIterator bit(inner); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              for (int ivar = 0; ivar < numComps; ivar++)
                destFab(iv, ivar) = psiFab(iv, ivar);
            }
        }

      // (2) interior boundary points of source boxes

      Vector<IntVectSet> IVS = a_IVSV[dit];
      for (int lcomp = 0; lcomp < IVS.size(); lcomp++)
        for (IVSIterator it(IVS[lcomp]); it.ok(); ++it)
          {
            IntVect iv = it();
            for (int ivar = 0; ivar < numComps; ivar++)
              destFab(iv, ivar) = psiFab(iv, ivar);
          }
    }
}


// ---------------------------------------------------------
void
fullIntVectSets(LayoutData< BitSet >& a_IVSVfull,
                const LayoutData< Vector<IntVectSet> >& a_IVSV)
{
  const BoxLayout& grids = a_IVSV.boxLayout();
  a_IVSVfull.define(grids);
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      const Vector<IntVectSet>& IVSvec = a_IVSV[dit];
      int ncomps = IVSvec.size();
      BitSet& fullvec = a_IVSVfull[dit];
      fullvec = BitSet(ncomps, false);
      for (int lcomp = 0; lcomp < ncomps; lcomp++)
        {
          const IntVectSet& IVS = IVSvec[lcomp];
          CH_assert(!IVS.isEmpty()); // shouldn't be empty, should it?
          Box container(IVS.minBox());
          if (container.numPts() == IVS.numPts()) fullvec.setTrue(lcomp);
        }
    }
}


// ---------------------------------------------------------
void
interiorBoundaryNodes(LayoutData< Vector<Box> >& a_IVSV,
                      const DisjointBoxLayout& a_boxes,
                      const Box& a_domain)
{
  ProblemDomain probdomain(a_domain);
  interiorBoundaryNodes(a_IVSV, a_boxes, probdomain);
}


// ---------------------------------------------------------
void
interiorBoundaryNodes(LayoutData< Vector<Box> >& a_IVSV,
                      const DisjointBoxLayout& a_boxes,
                      const ProblemDomain& a_domain)
{
  // boxesNodes contains the indices of the nodes of the boxes in a_boxes,
  // but they're represented as CELLs.

  // periodicTestBox is the domain box with one layer of cells removed
  // from each end in periodic directions.
  Box periodicTestBox(a_domain.domainBox());
  if (a_domain.isPeriodic())
    for (int idir = 0; idir < SpaceDim; idir++)
      if (a_domain.isPeriodic(idir))
        periodicTestBox.grow(idir, -1);

  // We use this to move things around in the periodic case.
  IntVect shiftMult = a_domain.domainBox().size();

  a_IVSV.define(a_boxes);
  for (DataIterator dit = a_boxes.dataIterator(); dit.ok(); ++dit)
    {
      Vector<Box>& IVSV = a_IVSV[dit];

      const Box& bxCells = a_boxes.get(dit);

      // bxNodes is CELL-centered, contains indices of nodes of bxCells
      // const Box& bxNodes = boxesNodes.get(dit);

      // on each side, find interior boundary nodes of bxCells
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide mySide = sit();

              Box faceboxOriginal = adjCellBox(bxCells, idir, mySide, 1);

              Box faceboxIn = faceboxOriginal & a_domain;
              // empty faceboxIn means physical boundary
              // IntVectSet thisIVS(faceboxIn);
              Vector<Box> thisIVS;
              Box facegrow(faceboxIn);
              if (! faceboxIn.isEmpty() )
                {
                  thisIVS.push_back(faceboxIn);
                  // facegrow (CELL-centered) = faceboxIn expanded by:
                  // one layer in idir, side INTO the box;
                  // one layer in all non-idir directions;
                  // and then truncated to physical domain.
                  facegrow.growDir(idir, flip(mySide), 1);
                  for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                    if (idirOther != idir) facegrow.grow(idirOther, 1);

                  facegrow &= a_domain;

                  // box and facegrow for idir==0, mySide==Side::Hi
                  //
                  //              +---+
                  //      +-------+-+ |
                  //      |       | | |
                  //      |       | | |
                  //      |       | | |
                  //      +-------+-+ |
                  //              +---+
                  //

                  // shell is facegrow with all boxes on this level removed.
                  // IntVectSet shell(facegrow);
                  Vector<Box> shell;
                  shell.push_back(facegrow);
                  // (use LayoutIterator to get _all_ boxes, even those in other procs)
                  for (LayoutIterator lit = a_boxes.layoutIterator(); lit.ok(); ++lit)
                    {
                      Box shiftedBox = a_boxes.get(lit);
                      // shell -= shiftedBox;
                      removeBoxFromBoxes(shell, shiftedBox); // NEW
                      // added by petermc, 4 Feb 2003
                      // only do this IF we're periodic _and_ both boxes
                      // adjoin the domain box boundary somewhere
                      if (a_domain.isPeriodic() &&
                          !periodicTestBox.contains(facegrow) &&
                          !periodicTestBox.contains(shiftedBox))
                        {
                          ShiftIterator shiftIt = a_domain.shiftIterator();
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                            {
                              IntVect shiftVect(shiftMult*shiftIt());
                              shiftedBox.shift(shiftVect);
                              // shell -= shiftedBox;
                              removeBoxFromBoxes(shell, shiftedBox); // NEW
                              shiftedBox.shift(-shiftVect);
                            }
                        }
                    }

                  // define interiorPlus as the cells removed from facegrow
                  // to yield shell.
                  // const IntVectSet interiorPlus(facegrow - shell);
                  // IntVectSet interiorPlus(facegrow);
                  Vector<Box> interiorPlus;
                  interiorPlus.push_back(facegrow);
                  // interiorPlus -= shell;
                  removeBoxesFromBoxes(interiorPlus, shell); // NEW
                  // initialize interior as interiorPlus.
                  thisIVS = interiorPlus;

                  Box offsetBox(IntVect::Zero, IntVect::Unit);
                  // intersect thisIVS with interiorPlus with offsets
                  for (BoxIterator offset(offsetBox); offset.ok(); ++offset)
                    {
                      // IntVectSet interiorPlus_off(interiorPlus);
                      Vector<Box> interiorPlus_off;
                      interiorPlus_off.append(interiorPlus);
                      // interiorPlus_off.shift(offset());
                      shiftBoxes(interiorPlus_off, offset()); // NEW
                      // thisIVS &= interiorPlus_off;
                      intersectBoxes(thisIVS, interiorPlus_off); // NEW
                    }
                }

              if (thisIVS.size() > 0)
                {
                  // IntVectSet IVS(bxNodes);
                  // petermc, 16 July 2003:  bxNodes is too big and forces
                  // a TreeIntVectSet.
                  // IntVectSet IVS(facegrow);
                  Vector<Box> IVS;
                  if (! facegrow.isEmpty() )
                    IVS.push_back(facegrow);
                  // IVS &= thisIVS;
                  intersectBoxes(IVS, thisIVS);
                  // To IVSV, append IVS, making sure there is no overlap.
                  appendBoxes(IVSV, IVS); // NEW
                  // remove nodes that are already present in other components
                  // for (int lcomp = 0; lcomp < IVSV.size(); lcomp++)
                    // IVS -= IVSV[lcomp];
                  // From IVS, remove nodes that are already present in IVSV.
                  // add the new nodes as a new component of IVSV
                  // if (IVS.size() > 0) IVSV.append(IVS);
                }
            }
        }
      cellsToNodes(IVSV);
    }
}


// ---------------------------------------------------------
void
interiorBoundaryNodes(LayoutData< Vector<Box> >& a_IVSV,
                      const DisjointBoxLayout& a_dest,
                      const DisjointBoxLayout& a_src,
                      const Box& a_domain)
{
  ProblemDomain probdomain(a_domain);
  interiorBoundaryNodes(a_IVSV, a_dest, a_src, probdomain);
}


// ---------------------------------------------------------
void
interiorBoundaryNodes(LayoutData< Vector<Box> >& a_IVSV,
                      const DisjointBoxLayout& a_dest,
                      const DisjointBoxLayout& a_src,
                      const ProblemDomain& a_domain)
{
  DataIterator dit = a_dest.dataIterator();

  // destNodes contains the indices of the NODEs of the boxes in a_dest.
  BoxLayout destNodes;
  getSurroundingNodes(destNodes, a_dest);

  Box periodicTestBox(a_domain.domainBox());
  if (a_domain.isPeriodic())
    for (int idirPeriodic = 0; idirPeriodic < SpaceDim; idirPeriodic++)
      if (a_domain.isPeriodic(idirPeriodic))
        periodicTestBox.grow(idirPeriodic, -1);

  IntVect shiftMult = a_domain.domainBox().size();

  // destSrcIntersec[dit] is true iff the current source box
  // shares nodes with a_dest[dit].
  LayoutData<bool> destSrcIntersect(a_dest);

  // cout << "interiorBoundaryNodes:  from " << a_src.size()
  //      << " boxes to " << a_dest.size() << " boxes"
  //      << " on " << a_domain.longside() << endl;
  // Define a_IVSV to have a Vector<IntVectSet> for each box of a_dest.
  a_IVSV.define(a_dest);
  for (LayoutIterator lit = a_src.layoutIterator(); lit.ok(); ++lit)
    {
      const Box& srcbox = a_src.get(lit);

      // srcNodesBox is node-centered, contains indices of nodes of srcbox
      Box srcNodesBox = surroundingNodes(srcbox);

      for (dit.begin(); dit.ok(); ++dit)
        {
          Box destNodesBox(destNodes.get(dit));
          destNodesBox &= srcNodesBox;
          destSrcIntersect[dit] = (! destNodesBox.isEmpty() );
        }

      // on each side, find interior boundary nodes of srcbox
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide mySide = sit();

              Box faceboxOriginal = adjCellBox(srcbox, idir, mySide, 1);

              Box faceboxIn = faceboxOriginal & a_domain;
              // empty faceboxIn means physical boundary
              // IntVectSet srcIVS(faceboxIn);
              Vector<Box> srcIVS;
              srcIVS.push_back(faceboxIn);
              Box facegrow(faceboxIn);
              if (! faceboxIn.isEmpty())
                {
                  // facegrow = faceboxIn expanded by:
                  // one layer in idir, side INTO the box;
                  // one layer in all non-idir directions;
                  // and then truncated to physical domain.
                  facegrow.growDir(idir, flip(mySide), 1);
                  for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                    if (idirOther != idir) facegrow.grow(idirOther, 1);

                  facegrow &= a_domain;

                  // box and facegrow for idir==0, mySide==Side::Hi
                  //
                  //              +---+
                  //      +-------+-+ |
                  //      |       | | |
                  //      |       | | |
                  //      |       | | |
                  //      +-------+-+ |
                  //              +---+
                  //

                  // shell is facegrow with all boxes on this level removed.
                  // IntVectSet shell(facegrow);
                  Vector<Box> shell;
                  shell.push_back(facegrow);
                  // (use LayoutIterator to get _all_ boxes, even those in other procs)
                  for (LayoutIterator lit2 = a_src.layoutIterator(); lit2.ok(); ++lit2)
                    {
                      Box shiftedBox = a_src.get(lit2);
                      // shell -= shiftedBox;
                      removeBoxFromBoxes(shell, shiftedBox); // NEW
                      // added by petermc, 4 Feb 2003
                      // only do this IF we're periodic _and_ both boxes
                      // adjoin the domain box boundary somewhere
                      if (a_domain.isPeriodic() &&
                          !periodicTestBox.contains(faceboxIn) &&
                          !periodicTestBox.contains(shiftedBox))
                        {
                          ShiftIterator shiftIt = a_domain.shiftIterator();
                          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                            {
                              IntVect shiftVect(shiftMult*shiftIt());
                              shiftedBox.shift(shiftVect);
                              // shell -= shiftedBox;
                              removeBoxFromBoxes(shell, shiftedBox); // NEW
                              shiftedBox.shift(-shiftVect);
                            }
                        }

                    }

                  // define interiorPlus as the cells removed from facegrow
                  // to yield shell.
                  // const IntVectSet interiorPlus(facegrow - shell);
                  // IntVectSet interiorPlus(facegrow);
                  Vector<Box> interiorPlus;
                  interiorPlus.push_back(facegrow);
                  // interiorPlus -= shell;
                  removeBoxesFromBoxes(interiorPlus, shell); // NEW
                  // initialize interior as interiorPlus.
                  srcIVS = interiorPlus;

                  Box offsetBox(IntVect::Zero, IntVect::Unit);
                  // intersect srcIVS with interiorPlus with offsets
                  for (BoxIterator offset(offsetBox); offset.ok(); ++offset)
                    {
                      // IntVectSet interiorPlus_off(interiorPlus);
                      Vector<Box> interiorPlus_off;
                      interiorPlus_off.append(interiorPlus);
                      // interiorPlus_off.shift(offset());
                      shiftBoxes(interiorPlus_off, offset()); // NEW
                      // srcIVS &= interiorPlus_off;
                      intersectBoxes(srcIVS, interiorPlus_off); // NEW
                    }
                }

              if (srcIVS.size() > 0)
                {
                  for (dit.begin(); dit.ok(); ++dit)
                    {
                      if (destSrcIntersect[dit])
                        // get nodes of destination box
                        {
                          Box destNodesBox(destNodes.get(dit));
                          destNodesBox &= srcNodesBox;
                          if (! destNodesBox.isEmpty() )
                            {
                              // petermc, 16 Jul 2003:  destNodesBox was
                              // too big and forced TreeIntVectSet.
                              // Replace these two lines with reverse order.
                              // IntVectSet IVS(destNodesBox);
                              // IVS &= srcIVS;
                              // IntVectSet IVS(srcIVS);
                              Vector<Box> IVS = srcIVS;
                              // petermc, 16 nov 2005:  convert NODE to CELL
                              Box destNodesBoxShift(destNodesBox);
                              destNodesBoxShift.shiftHalf(IntVect::Unit);
                              intersectBoxes(IVS, destNodesBoxShift);
                              // IVS: nodes of destination box that are also
                              // interior boundary nodes of source box
                              Vector<Box>& IVSV = a_IVSV[dit];
                              // To IVSV, append IVS, making sure there is no overlap.
                              appendBoxes(IVSV, IVS);
                              // remove nodes that are already present in other components
                              // for (int lcomp = 0; lcomp < IVSV.size(); lcomp++)
                              // removeBoxesFromBoxes(IVS, IVSV[lcomp]); // NEW
                                // IVS -= IVSV[lcomp];
                              // add the new nodes as a new component
                              // if (! IVS.isEmpty() ) IVSV.push_back(IVS);
                              // if (IVS.size() > 0) IVSV.append(IVS);
                            }
                        }
                    }
                }
            }
        }
    }

  for (dit.begin(); dit.ok(); ++dit)
    cellsToNodes(a_IVSV[dit]);
}


// ---------------------------------------------------------
void
exteriorBoundaryNodes(LayoutData< Vector<Box> >& a_exterior,
                      const LayoutData< Vector<Box> >& a_interior,
                      const DisjointBoxLayout& a_boxes)
{
  a_exterior.define(a_boxes);
  for (DataIterator dit = a_boxes.dataIterator(); dit.ok(); ++dit)
    {
      Box box = a_boxes.get(dit);
      // Vector<IntVectSet> ibn = a_interior[dit];
      // Vector<Box> ibn = a_interior[dit];
      Vector<Box> ibn;
      ibn.append(a_interior[dit]);
      nodesToCells(ibn);
      // Vector<IntVectSet>& ebn = a_exterior[dit];
      Vector<Box>& ebn = a_exterior[dit];
      // on each side, find exterior boundary nodes
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              Side::LoHiSide mySide = sit();

              // faceNodes will contain the cells with indices
              // that are the same as the nodes on this face.
              Box faceNodes = adjCellBox(box, idir, mySide, 1);
              // Expand by one in all directions other than idir.
              for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                if (idirOther != idir) faceNodes.growHi(idirOther, 1);

              // IVS will be the face nodes without interior boundary nodes.
              // What's left will be the exterior boundary nodes.
              // IntVectSet IVS(faceNodes);
              Vector<Box> IVS;
              IVS.push_back(faceNodes);
              // remove the interior boundary nodes from IVS.
              removeBoxesFromBoxes(IVS, ibn); // NEW
              // for (int lcomp = 0; lcomp < ibn.size(); lcomp++)
              // IVS -= ibn[lcomp];
              // add remaining nodes as a new component of
              // exterior boundary nodes.
              // if (! IVS.isEmpty() ) ebn.push_back(IVS);
              appendBoxes(ebn, IVS); // NEW
            }
        }
      cellsToNodes(ebn);
    }
}


// ---------------------------------------------------------
void appendBoxes(Vector<Box>&         a_boxes,
                 const Vector<Box>&   a_new)
{
  // We assume the boxes in a_boxes are disjoint in input;
  // they'll also be disjoint on output.

  // If nothing is to be added, then a_boxes is unchanged.
  if (a_new.size() == 0) return;
  // newTotally will hold points of a_new that are NOT in any of a_boxes.
  Vector<Box> newTotally;
  for (int i = 0; i < a_new.size(); i++)
    {
      Vector<Box> revisedBoxes;
      // revisedBoxes := a_new[i] - a_boxes
      removeBoxesFromBox(revisedBoxes, a_boxes, a_new[i]);
      if (revisedBoxes.size() > 0)
        newTotally.append(revisedBoxes);
    }
  if (newTotally.size() > 0)
    a_boxes.append(newTotally);
}


// ---------------------------------------------------------
void removeBoxFromBoxes(Vector<Box>&   a_boxes,
                        const Box&     a_remove)
{
  Vector<Box> newBoxes;
  for (int i = 0; i < a_boxes.size(); i++)
    {
      Vector<Box> revisedBoxes;
      // revisedBoxes := a_boxes[i] - a_remove
      removeBoxFromBox(revisedBoxes, a_remove, a_boxes[i]);
      if (revisedBoxes.size() > 0)
        newBoxes.append(revisedBoxes);
    }
  a_boxes = newBoxes;
}


// ---------------------------------------------------------
void removeBoxesFromBox(Vector<Box>&         a_boxes,
                        const Vector<Box>&   a_remove,
                        const Box&           a_base)
{
  a_boxes.clear();
  a_boxes.push_back(a_base);
  removeBoxesFromBoxes(a_boxes, a_remove);
}


// ---------------------------------------------------------
void removeBoxFromBox(Vector<Box>&   a_boxes,
                      const Box&     a_remove,
                      const Box&     a_base)
{
  a_boxes.clear();
  if (a_base.isEmpty())
    { // started with empty box; end with empty list
      return;
    }

  IndexType ix = a_base.ixType();
  CH_assert( a_remove.ixType() == ix );
  Box bxRemove(a_remove);
  bxRemove &= a_base;
  if (bxRemove.isEmpty())
    { // nothing to remove
      a_boxes.push_back(a_base);
      return;
    }

  // a_boxes := a_base - bxRemove
  int specialDir = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    if ( (bxRemove.smallEnd(idir) == a_base.smallEnd(idir)) ||
         (bxRemove.bigEnd(idir) == a_base.bigEnd(idir)) )
      {
        specialDir = idir;
        break;
      }

  // specialDir  0  1  2
  // dir1        1  0  0
  // dir2        2  2  1
  int dir1 = (specialDir == 0) ? 1 : 0;
  int dir2 = (specialDir == 2) ? 1 : 2;

  int baseSpecLo = a_base.smallEnd(specialDir);
  int baseSpecHi = a_base.bigEnd(specialDir);
  int base1Lo = a_base.smallEnd(dir1);
  int base1Hi = a_base.bigEnd(dir1);
  int base2Lo = a_base.smallEnd(dir2);
  int base2Hi = a_base.bigEnd(dir2);

  int removeSpecLo = bxRemove.smallEnd(specialDir);
  int removeSpecHi = bxRemove.bigEnd(specialDir);
  int remove1Lo = bxRemove.smallEnd(dir1);
  int remove1Hi = bxRemove.bigEnd(dir1);
  int remove2Lo = bxRemove.smallEnd(dir2);
  int remove2Hi = bxRemove.bigEnd(dir2);

  // [LO, LO, LO] : [lo-1, HI, HI]
  // [hi+1, LO, LO] : [HI, HI, HI]
  // [lo, lo, hi+1] : [hi, HI, HI]
  // [lo, hi+1, LO] : [hi, HI, hi]
  // [lo, LO, LO] : [hi, hi, lo-1]
  // [lo, LO, lo] : [hi, lo-1, HI]

  // [LO, LO, LO] : [lo-1, HI, HI]
  if (baseSpecLo <= removeSpecLo-1)
    a_boxes.push_back(Box(baseSpecLo * BASISV(specialDir) +
                          base1Lo * BASISV(dir1) +
                          base2Lo * BASISV(dir2),
                          (removeSpecLo-1) * BASISV(specialDir) +
                          base1Hi * BASISV(dir1) +
                          base2Hi * BASISV(dir2),
                          ix));

  // [hi+1, LO, LO] : [HI, HI, HI]
  if (removeSpecHi+1 <= baseSpecHi)
    a_boxes.push_back(Box((removeSpecHi+1) * BASISV(specialDir) +
                          base1Lo * BASISV(dir1) +
                          base2Lo * BASISV(dir2),
                          baseSpecHi * BASISV(specialDir) +
                          base1Hi * BASISV(dir1) +
                          base2Hi * BASISV(dir2),
                          ix));

  // [lo, lo, hi+1] : [hi, HI, HI]
  if (remove2Hi+1 <= base2Hi)
    a_boxes.push_back(Box(removeSpecLo * BASISV(specialDir) +
                          remove1Lo * BASISV(dir1) +
                          (remove2Hi+1) * BASISV(dir2),
                          removeSpecHi * BASISV(specialDir) +
                          base1Hi * BASISV(dir1) +
                          base2Hi * BASISV(dir2),
                          ix));

  // [lo, hi+1, LO] : [hi, HI, hi]
  if (remove1Hi+1 <= base1Hi)
    a_boxes.push_back(Box(removeSpecLo * BASISV(specialDir) +
                          (remove1Hi+1) * BASISV(dir1) +
                          base2Lo * BASISV(dir2),
                          removeSpecHi * BASISV(specialDir) +
                          base1Hi * BASISV(dir1) +
                          remove2Hi * BASISV(dir2),
                          ix));

  // [lo, LO, LO] : [hi, hi, lo-1]
  if (base2Lo <= remove2Lo-1)
    a_boxes.push_back(Box(removeSpecLo * BASISV(specialDir) +
                          base1Lo * BASISV(dir1) +
                          base2Lo * BASISV(dir2),
                          removeSpecHi * BASISV(specialDir) +
                          remove1Hi * BASISV(dir1) +
                          (remove2Lo-1) * BASISV(dir2),
                          ix));

  // [lo, LO, lo] : [hi, lo-1, HI]
  if (base1Lo <= remove1Lo-1)
    a_boxes.push_back(Box(removeSpecLo * BASISV(specialDir) +
                          base1Lo * BASISV(dir1) +
                          remove2Lo * BASISV(dir2),
                          removeSpecHi * BASISV(specialDir) +
                          (remove1Lo-1) * BASISV(dir1) +
                          base2Hi * BASISV(dir2),
                          ix));
}


// ---------------------------------------------------------
void removeBoxesFromBoxes(Vector<Box>&         a_boxes,
                          const Vector<Box>&   a_remove)
{
  for (int j = 0; j < a_remove.size(); j++)
    {
      // a_boxes -= a_remove[j]
      removeBoxFromBoxes(a_boxes, a_remove[j]);
    }
}


// ---------------------------------------------------------
void intersectBoxes(Vector<Box>&         a_boxes,
                    const Vector<Box>&   a_new)
{
  // a_boxes := a_boxes & a_new

  Vector<Box> revisedBoxes;
  for (int i = 0; i < a_boxes.size(); i++)
    {
      const Box& bxOrig = a_boxes[i];
      for (int j = 0; j < a_new.size(); j++)
        {
          Box bxNew(a_new[j]);
          bxNew &= bxOrig;
          if (! bxNew.isEmpty() )
            revisedBoxes.push_back(bxNew);
        }
    }
  a_boxes = revisedBoxes;
}


// ---------------------------------------------------------
void intersectBoxes(Vector<Box>&   a_boxes,
                    const Box&     a_new)
{
  // Replace a_boxes with a_boxes & a_new.
  // a_boxes[i] &= a_new
  // But it's slightly more complicated than this,
  // because some boxes may become empty,
  // and the returned Vector should not include any empty boxes.
  Vector<Box> revisedBoxes;
  for (int i = 0; i < a_boxes.size(); i++)
    {
      Box revisedBox = a_boxes[i];
      revisedBox &= a_new;
      if (! revisedBox.isEmpty() )
        revisedBoxes.push_back(revisedBox);
    }
  a_boxes = revisedBoxes;
}


// ---------------------------------------------------------
void shiftBoxes(Vector<Box>&     a_boxes,
                const IntVect&   a_offset)
{
  for (int i = 0; i < a_boxes.size(); i++)
    {
      a_boxes[i].shift(a_offset);
    }
}


// ---------------------------------------------------------
void cellsToNodes(Vector<Box>&   a_boxes)
{
  IntVect shift = -IntVect::Unit;
  for (int i = 0; i < a_boxes.size(); i++)
    {
      CH_assert( a_boxes[i].ixType() == IndexType::TheCellType() );
      a_boxes[i].shiftHalf(shift);
    }
}


// ---------------------------------------------------------
void nodesToCells(Vector<Box>&   a_boxes)
{
  IntVect shift = IntVect::Unit;
  for (int i = 0; i < a_boxes.size(); i++)
    {
      CH_assert( a_boxes[i].ixType() == IndexType::TheNodeType() );
      a_boxes[i].shiftHalf(shift);
    }
}


// ---------------------------------------------------------
void
copyInteriorNodes(LevelData<NodeFArrayBox>& a_dest,
                  const LevelData<NodeFArrayBox>& a_src,
                  const LayoutData< Vector<Box> >& a_IVSV)
{
  const DisjointBoxLayout& destBoxes = a_dest.getBoxes();
  const DisjointBoxLayout& srcBoxes = a_src.getBoxes();
  CH_assert(destBoxes.numBoxes(procID()) == a_IVSV.boxLayout().numBoxes(procID()));

  // (1) define psi on destBoxes
  // (2) copy a_src to psi;
  // (3) for each box in destBoxes
  //         for each interior node of srcBoxes in it,
  //             copy psi to a_dest there.

  // (1) define psi on destBoxes
  int numComps = a_src.nComp();
  BoxLayoutData<NodeFArrayBox> psi(destBoxes, numComps);

  // (2) copy a_src to psi;
  // petermc, 13 May 2003:  change Copier ghosts from Unit to Zero.
  Copier myCopier(srcBoxes, destBoxes, IntVect::Zero);
  a_src.copyTo(a_src.interval(), psi, psi.interval(), myCopier);

  // (3) for each box in destBoxes
  for (DataIterator dit = destBoxes.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& psiFab = psi[dit].getFab();
      FArrayBox& destFab = a_dest[dit].getFab();

      // copy psi to a_dest there.

      // destNodes:  indices of _all_ nodes in destBoxes[dit]
      Box destNodes = destBoxes.get(dit);
      for (int dir = 0; dir < SpaceDim; dir++)
        destNodes.growHi(dir, 1);

      // (1) interior of each source box

      for (LayoutIterator lit = srcBoxes.layoutIterator(); lit.ok(); ++lit)
        {
          const Box thisBox = srcBoxes.get(lit);
          Box inner(thisBox);
          for (int dir = 0; dir < SpaceDim; dir++)
            inner.growLo(dir, -1);

          inner &= destNodes;

          if (! inner.isEmpty() )
            {
              // destFab.copy(psiFab, inner);  // incompatible box types
              for (BoxIterator bit(inner); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  for (int ivar = 0; ivar < numComps; ivar++)
                    destFab(iv, ivar) = psiFab(iv, ivar);
                }
            }
        }

      // (2) interior boundary points of source boxes

      Vector<Box> IVS = a_IVSV[dit];
      for (int i = 0; i < IVS.size(); i++)
        destFab.copy(psiFab, IVS[i]);
      //      for (int lcomp = 0; lcomp < IVS.size(); lcomp++)
      //        for (IVSIterator it(IVS[lcomp]); it.ok(); ++it)
      //          {
      //            IntVect iv = it();
      //            for (int ivar = 0; ivar < numComps; ivar++)
      //              destFab(iv, ivar) = psiFab(iv, ivar);
      //          }
    }
}


// ---------------------------------------------------------
void
zeroBoundaryNodes(BoxLayoutData<NodeFArrayBox>& a_dest,
                  const LayoutData< Vector<Box> >& a_IVSV)
{
  int numComps = a_dest.nComp();
  for (DataIterator dit = a_dest.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& destFab = a_dest[dit].getFab();
      Vector<Box> IVS = a_IVSV[dit];
      for (int lcomp = 0; lcomp < IVS.size(); lcomp++)
        destFab.setVal(0., IVS[lcomp], 0, numComps);

      //        for (IVSIterator it(IVS[lcomp]); it.ok(); ++it)
      //          {
      //            IntVect iv = it();
      //            for (int ivar = 0; ivar < numComps; ivar++)
      //              destFab(iv, ivar) = 0.;
      //          }
    }
}


// ---------------------------------------------------------
void
exteriorAndInteriorNodes(Vector< LayoutData< Vector<Box> >* >& a_exterior,
                         Vector< LayoutData< Vector<Box> >* >& a_intFinerCoarsened,
                         const Vector<DisjointBoxLayout>& a_layouts,
                         const Vector<ProblemDomain>& a_domain,
                         const Vector<int>& a_nRefFine)
{
  int numLevels = a_layouts.size();
  a_exterior.resize(numLevels);
  a_intFinerCoarsened.resize(numLevels-1);
  for (int lev = 0; lev < numLevels; lev++)
    {
      LayoutData< Vector<Box> > interior;
      interiorBoundaryNodes(interior, a_layouts[lev], a_domain[lev]);
      a_exterior[lev] = new LayoutData< Vector<Box> >;
      exteriorBoundaryNodes(*(a_exterior[lev]), interior, a_layouts[lev]);

      if (lev < numLevels-1)
        {
          DisjointBoxLayout coarsenedFinerGrids;
          coarsen(coarsenedFinerGrids, a_layouts[lev+1], a_nRefFine[lev]);
          LayoutData< Vector<Box> >* IVSVintFinerCoarsened =
            new LayoutData< Vector<Box> >;
          interiorBoundaryNodes(*IVSVintFinerCoarsened,
                                a_layouts[lev],
                                coarsenedFinerGrids,
                                a_domain[lev]);
          a_intFinerCoarsened[lev] = IVSVintFinerCoarsened;
        }
    }
}


// ---------------------------------------------------------
void
exteriorAndInteriorNodes(Vector< LayoutData< Vector<IntVectSet> >* >& a_exterior,
                         Vector< LayoutData< Vector<IntVectSet> >* >& a_intFinerCoarsened,
                         const Vector<DisjointBoxLayout>& a_layouts,
                         const Vector<ProblemDomain>& a_domain,
                         const Vector<int>& a_nRefFine)
{
  int numLevels = a_layouts.size();
  a_exterior.resize(numLevels);
  a_intFinerCoarsened.resize(numLevels-1);
  for (int lev = 0; lev < numLevels; lev++)
    {
      LayoutData< Vector<IntVectSet> > interior;
      interiorBoundaryNodes(interior, a_layouts[lev], a_domain[lev]);
      a_exterior[lev] = new LayoutData< Vector<IntVectSet> >;
      exteriorBoundaryNodes(*(a_exterior[lev]), interior, a_layouts[lev]);

      if (lev < numLevels-1)
        {
          DisjointBoxLayout coarsenedFinerGrids;
          coarsen(coarsenedFinerGrids, a_layouts[lev+1], a_nRefFine[lev]);
          LayoutData< Vector<IntVectSet> >* IVSVintFinerCoarsened =
            new LayoutData< Vector<IntVectSet> >;
          interiorBoundaryNodes(*IVSVintFinerCoarsened,
                                a_layouts[lev],
                                coarsenedFinerGrids,
                                a_domain[lev]);
          a_intFinerCoarsened[lev] = IVSVintFinerCoarsened;
        }
    }
}


// ---------------------------------------------------------
void
getMaskInteriorNodes(LevelData<NodeFArrayBox>& a_mask,
                     const DisjointBoxLayout& a_dest,
                     const DisjointBoxLayout* a_srcPtr,
                     const Box& a_domain,
                     int a_onoff)
{
  ProblemDomain probdomain(a_domain);
  getMaskInteriorNodes(a_mask, a_dest, a_srcPtr, probdomain, a_onoff);
}


// ---------------------------------------------------------
void
getMaskInteriorNodes(LevelData<NodeFArrayBox>& a_mask,
                     const DisjointBoxLayout& a_dest,
                     const DisjointBoxLayout* a_srcPtr,
                     const ProblemDomain& a_domain,
                     int a_onoff)
{
  // destNodes contains the NODEs of the boxes in a_dest.
  BoxLayout destNodes;
  getSurroundingNodes(destNodes, a_dest);

  // srcNodes contains the NODEs of the boxes in *a_srcPtr.
  BoxLayout srcNodes;
  if (a_srcPtr == NULL)
    srcNodes = destNodes;
  else
    getSurroundingNodes(srcNodes, *a_srcPtr);

  // srcCount gives the number of corners that are in boxes.
  IntVect all2 = 2 * IntVect::Unit;
  int cmax = all2.product();
  BoxLayoutData< BaseFab<int> > srcCount(srcNodes, 1);
  for (DataIterator ditSrc = srcCount.dataIterator(); ditSrc.ok(); ++ditSrc)
    {
      BaseFab<int>& srcFab = srcCount[ditSrc];
      srcFab.setVal(cmax);
      // Halve srcFab on each face.
      halveIntFaces(srcFab, srcFab.box());
    }

  a_mask.define(a_dest, 1);
  // BoxLayoutData< BaseFab<int> > destCount(destNodes, 1);
  { // scope of intCount
    LayoutData < Vector< RefCountedPtr< BaseFab<int> > > > intCount;
    const Interval intvl0(0, 0);

    // pdNodes is NODE-centered version of a_domain.
    bool* periodicities = new bool[SpaceDim];
    for (int idir = 0; idir < SpaceDim; idir++)
      periodicities[idir] = a_domain.isPeriodic(idir);
    Box domainNodes = surroundingNodes(a_domain.domainBox());
    ProblemDomain pdNodes(domainNodes, periodicities);
    delete [] periodicities;

    srcCount.generalCopyTo(destNodes, intCount, intvl0, pdNodes);
    for (DataIterator ditDst = destNodes.dataIterator(); ditDst.ok(); ++ditDst)
      {
        Vector< RefCountedPtr< BaseFab<int> > >& intFab = intCount[ditDst];
        Box bxNodes = destNodes.get(ditDst);
        BaseFab<int> sumFab(bxNodes, 1);
        // plusReduce(sumFab, intFab);
        sumFab.setVal(0);
        for (int ind = 0; ind < intFab.size(); ind++)
          {
            const BaseFab<int>& piece = *(intFab[ind]);
            // Add over domain of intersection of sumFab and piece.
            // sumFab.plus(piece, 0, 0, 1);
            // piece.box() should be contained in sumFab.box().
            FORT_BASEFABINTPLUS(CHF_FIA(sumFab),
                                CHF_CONST_FIA(piece),
                                CHF_BOX(piece.box()));
          }

        // Set mask to 0 where sumFab == cmax, and 1 elsewhere.
        FArrayBox& maskFab = a_mask[ditDst].getFab();
        FORT_MASKVALUE(CHF_FRA1(maskFab, 0),
                       CHF_CONST_FIA1(sumFab, 0),
                       CHF_BOX(bxNodes),
                       CHF_CONST_INT(cmax),
                       CHF_CONST_INT(a_onoff));
      }
  }
}


// ---------------------------------------------------------
void
getMaskValidNodes(LevelData<NodeFArrayBox>& a_mask,
                  const DisjointBoxLayout& a_layout,
                  const DisjointBoxLayout* a_finerCoarsenedPtr,
                  const Box& a_domain)
{
  ProblemDomain probdomain(a_domain);
  getMaskValidNodes(a_mask, a_layout, a_finerCoarsenedPtr, probdomain);
}


// ---------------------------------------------------------
void
getMaskValidNodes(LevelData<NodeFArrayBox>& a_mask,
                  const DisjointBoxLayout& a_layout,
                  const DisjointBoxLayout* a_finerCoarsenedPtr,
                  const ProblemDomain& a_domain)
{
  // Set a_mask to 1 on interior nodes of a_layout.
  getMaskInteriorNodes(a_mask,
                       a_layout, NULL, a_domain, 1);
  if (a_finerCoarsenedPtr != NULL)
    {
      // Set maskFinerCoarsened to 0 on interior nodes of *a_finerCoarsenedPtr.
      LevelData<NodeFArrayBox> maskFinerCoarsened;
      getMaskInteriorNodes(maskFinerCoarsened,
                           a_layout, a_finerCoarsenedPtr, a_domain, 0);
      for (DataIterator dit = a_layout.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& maskFab = a_mask[dit].getFab();
          const FArrayBox& maskFiner = maskFinerCoarsened[dit].getFab();
          maskFab *= maskFiner;
        }
    }
}


// ---------------------------------------------------------
void
getMaskValidNodes(Vector< LevelData<NodeFArrayBox>* >& a_masks,
                  const Vector<DisjointBoxLayout>& a_layouts,
                  const Vector<ProblemDomain>& a_domain,
                  const Vector<int>& a_nRefFine)
{
  int numLevels = a_layouts.size();
  a_masks.resize(numLevels);
  for (int lev = 0; lev < numLevels; lev++)
    {
      a_masks[lev] = new LevelData<NodeFArrayBox>;
      DisjointBoxLayout* coarsenedFinerGridsPtr;
      if (lev < numLevels-1)
        { // take finer level into account
          coarsenedFinerGridsPtr = new DisjointBoxLayout;
          coarsen(*coarsenedFinerGridsPtr, a_layouts[lev+1], a_nRefFine[lev]);
        }
      else
        { // no finer level
          coarsenedFinerGridsPtr = NULL;
        }

      getMaskValidNodes(*(a_masks[lev]),
                        a_layouts[lev],
                        coarsenedFinerGridsPtr,
                        a_domain[lev]);

      if (coarsenedFinerGridsPtr != NULL)
        delete coarsenedFinerGridsPtr;
    }
}


// ---------------------------------------------------------
void getSurroundingNodes(BoxLayout&         a_gridsNodes,
                         const BoxLayout&   a_gridsCells)
{
  a_gridsNodes.deepCopy(a_gridsCells);
  a_gridsNodes.surroundingNodes();
  a_gridsNodes.close();
}


// ---------------------------------------------------------
void getEnclosedCells(BoxLayout&         a_gridsCells,
                      const BoxLayout&   a_gridsNodes)
{
  a_gridsCells.deepCopy(a_gridsNodes);
  a_gridsCells.enclosedCells();
  a_gridsCells.close();
}


// ---------------------------------------------------------
void halveIntFaces(BaseFab<int>& a_intFab,
                   const Box& a_bx)
{
  const Box& dataBox = a_intFab.box();
  SideIterator sit;
  for (sit.begin(); sit.ok(); sit.next())
    {
      Side::LoHiSide mySide = sit();
      IntVect corner = a_bx.sideEnd(mySide);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          // Set face to mySide-idir-face of a_bx
          Box face(a_bx);
          face.setRange(idir, corner[idir]);

          face &= dataBox;
          if ( !face.isEmpty() )
            {
              FORT_HALVEINT(CHF_FIA(a_intFab),
                            CHF_BOX(face));
            }
        }
    }
}

#include "NamespaceFooter.H"
