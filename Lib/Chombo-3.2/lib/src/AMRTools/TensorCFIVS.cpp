#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TensorCFIVS.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "NamespaceHeader.H"

bool
TensorCFIVS::isDefined() const
{
  return(isdefined);
}

const IntVectSet&
TensorCFIVS::getFineIVS() const
{
  CH_assert(isdefined);
  return fiinterp_ivs;
}

void
TensorCFIVS::setDefaultValues()
{
  isdefined = false;
  fiinterp_ivs.define();
}
TensorCFIVS::TensorCFIVS()
{
  setDefaultValues();
}

TensorCFIVS::~TensorCFIVS()
{
  setDefaultValues();
}

TensorCFIVS::TensorCFIVS(
             const Box& Domain,
             const Box& box_in,
             const DisjointBoxLayout& Levboxes,
             int Direction, Side::LoHiSide hiorlo)
{
  setDefaultValues();
  ProblemDomain probdomain(Domain);
  define(probdomain,  box_in, Levboxes,
         Direction, hiorlo);

}


TensorCFIVS::TensorCFIVS(
             const ProblemDomain& Domain,
             const Box& box_in,
             const DisjointBoxLayout& Levboxes,
             int Direction, Side::LoHiSide hiorlo)
{
  setDefaultValues();
  define(Domain,  box_in, Levboxes,
         Direction, hiorlo);

}

void
TensorCFIVS::define(
              const Box& Domain,
              const Box& box_in,
              const DisjointBoxLayout& fine_boxes,
              int direction,
              Side::LoHiSide hiorlo)
{
  ProblemDomain probdomain(Domain);
  define(probdomain, box_in, fine_boxes, direction, hiorlo);
}



void
TensorCFIVS::define(
              const ProblemDomain& Domain,
              const Box& box_in,
              const DisjointBoxLayout& fine_boxes,
              int direction,
              Side::LoHiSide hiorlo)
{
  isdefined = true;
  CH_assert(direction >= 0);
  CH_assert(direction < SpaceDim);
  CH_assert(!Domain.isEmpty());
  CH_assert(Domain.contains(box_in));
  CH_assert(fine_boxes.checkPeriodic(Domain));
  CH_assert((hiorlo == Side::Lo) || (hiorlo == Side::Hi));

  //shift direction
  int hilo = sign(hiorlo);

  //create fine stencil
  Box finebox = box_in;
  Box edgebox;
  CH_assert((hilo ==1) || (hilo == -1));
  if (hilo == -1)
    {
      edgebox = adjCellLo(finebox,direction,1);
    }
  else
    {
      edgebox = adjCellHi(finebox,direction,1);
    }

  // for the tensor solver, we will need to have
  // corner points if they are adjacent to a
  // neigboring fine-level grid
  // we do this by growing the edgebox in the tangential
  // directions.  we do _this_ by growing in all directions,
  // then shrinking in the normal direction
  edgebox.grow(1);
  edgebox.grow(direction, -1);

  edgebox &= Domain;
  if (!edgebox.isEmpty())
    {
      fiinterp_ivs.define(edgebox);

      // also keep an associated ivs of interior cells
      // which correspond to ghost cells; will subtract
      // interior cells from this associated IVS.  at the
      // end of the iteration, the only cells left in this
      // associated IVS will be IV's corresponding to
      // cells in fineinterp_ivs which do not have
      // fine interior cells (i.e. corner cells which
      // we will not want to interpolate.  then we can
      // subtract the associated ivs from fiinterp_ivs
      // to remove these cells from the IVS, leaving
      // only ghost cells which are legitimate CF interface
      // cells
      Box interiorCellsBox(edgebox);
      if (hilo == -1)
        {
         interiorCellsBox.shift(direction, 1);
        }
      else
        {
          interiorCellsBox.shift(direction, -1);
        }
      IntVectSet nonInteriorCells(interiorCellsBox);

      LayoutIterator lit = fine_boxes.layoutIterator();
      Box periodicTestBox(Domain.domainBox());
      if (Domain.isPeriodic())
        {
          for (int idir=0; idir<SpaceDim; idir++)
            {
              if (Domain.isPeriodic(idir))
                periodicTestBox.grow(idir,-1);
            }
        }
      for (lit.reset(); lit.ok(); ++lit)
        {
          fiinterp_ivs -= fine_boxes[lit()];
          nonInteriorCells -= fine_boxes[lit()];

          // only do this IF we're periodic _and_ both boxes
          // adjoin the domain box boundary somewhere
          if (Domain.isPeriodic()
              && (!periodicTestBox.contains(edgebox)
                  || !periodicTestBox.contains(interiorCellsBox) )
              && !periodicTestBox.contains(fine_boxes[lit()]))
            {
              ShiftIterator shiftIt = Domain.shiftIterator();
              IntVect shiftMult = Domain.domainBox().size();
              Box shiftedBox = fine_boxes[lit()];
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect(shiftMult*shiftIt());
                  shiftedBox.shift(shiftVect);
                  fiinterp_ivs -= shiftedBox;
                  nonInteriorCells -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }
      // finally subtract off intvects remaining in nonIteriorCells
      // since they don't represent cf-interface cells
      IntVect shiftVect(hilo*BASISV(direction));
      nonInteriorCells.shift(shiftVect);
      fiinterp_ivs -= nonInteriorCells;
    }
}
#include "NamespaceFooter.H"
