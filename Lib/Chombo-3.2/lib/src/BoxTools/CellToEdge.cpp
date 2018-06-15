#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// CellToEdge.cpp
// Dan Martin, Fri, Jan 14, 2000

#include <cassert>
#include "DataIterator.H"
#include "CellToEdge.H"
#include "CellToEdgeF_F.H"

#include "NamespaceHeader.H"

/// -----------------------------------------------------------
void CellToEdge(const LevelData<FArrayBox>& a_cellData,
                LevelData<FluxBox>& a_edgeData)

{
  // this is just a wrapper around single-gridBox version
  DataIterator dit = a_cellData.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      CellToEdge(a_cellData[dit()], a_edgeData[dit()]);
    }
}

///
void CellToEdge(const FArrayBox& a_cellData,
                FluxBox& a_edgeData)
{

  // loop over components -- assumption is that in cell-centered
  // data, direction changes faster than component
  int cellcomp;
  for (int comp = 0; comp < a_edgeData.nComp(); comp++)
    {
      // loop over directions
      for (int dir = 0; dir < SpaceDim; dir++)
        {
          // define faces over which we can do averaging
          Box edgeBox = surroundingNodes(a_cellData.box(), dir);
          edgeBox.grow(dir,-1);
          edgeBox &= a_edgeData[dir].box();

          if (a_cellData.nComp() == a_edgeData.nComp())
            {
              // straightforward cell->edge averaging
              cellcomp = comp;
            }
          else
            {
              // each cell comp represents a different spatial direction
              cellcomp = SpaceDim*comp + dir;
            }
          FORT_CELLTOEDGE(CHF_CONST_FRA1(a_cellData, cellcomp),
                          CHF_FRA1(a_edgeData[dir], comp),
                          CHF_BOX(edgeBox),
                          CHF_CONST_INT(dir));
        }
    }
}

// ------------------------------------------------------------
void CellToEdge(const FArrayBox& a_cellData, FArrayBox& a_edgeData,
                const int a_dir)

{

  Box edgeBox = surroundingNodes(a_cellData.box(), a_dir);
  edgeBox.grow(a_dir,-1);
  edgeBox &= a_edgeData.box();

  int cellComp;
  for (int comp = 0; comp < a_edgeData.nComp(); comp++)
    {
      if (a_cellData.nComp() == a_edgeData.nComp())
        {
          // straightforward cell->edge averaging
          cellComp = comp;
        }
      else
        {
          // each cell comp represents a different spatial direction
          cellComp = SpaceDim*comp + a_dir;
        }

      FORT_CELLTOEDGE(CHF_CONST_FRA1(a_cellData, cellComp),
                      CHF_FRA1(a_edgeData, comp),
                      CHF_BOX(edgeBox),
                      CHF_CONST_INT(a_dir));
    }

}

// ------------------------------------------------------------
void CellToEdge(const FArrayBox& a_cellData, const int a_cellComp,
                FArrayBox& a_edgeData, const int a_edgeComp,
                const int a_dir)

{

  Box edgeBox = surroundingNodes(a_cellData.box(),a_dir);
  edgeBox.grow(a_dir, -1);
  edgeBox &= a_edgeData.box();

  FORT_CELLTOEDGE(CHF_CONST_FRA1(a_cellData,a_cellComp),
                  CHF_FRA1(a_edgeData,a_edgeComp),
                  CHF_BOX(edgeBox),
                  CHF_CONST_INT(a_dir));
}

#include "NamespaceFooter.H"
