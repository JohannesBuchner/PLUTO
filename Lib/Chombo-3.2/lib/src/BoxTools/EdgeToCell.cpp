#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// EdgeToCell.cpp
// Dan Martin, Fri, Jan 14, 2000

#include <cassert>

#include "DataIterator.H"
#include "EdgeToCell.H"
#include "EdgeToCellF_F.H"

#include "NamespaceHeader.H"

// ----------------------------------------------------------------
void EdgeToCell(const LevelData<FluxBox>& a_edgeData,
                LevelData<FArrayBox>& a_cellData)
{
  // this is just a wrapper around the single-grid version
  DataIterator dit = a_edgeData.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      EdgeToCell(a_edgeData[dit()], a_cellData[dit()]);
    }
}

// ----------------------------------------------------------------
void EdgeToCell(const FluxBox& a_edgeData,
                FArrayBox& a_cellData)
{
  // loop over components -- assumption is that in cell-centered
  // data, direction changes faster than component.
  for (int comp = 0; comp<a_edgeData.nComp(); comp++)
    {
      // loop over directions
      for (int dir = 0; dir < SpaceDim; dir++)
        {
          Box cellBox = a_cellData.box();
          cellBox &= a_edgeData.box();
          int cellcomp = SpaceDim*comp + dir;
          FORT_EDGETOCELL(CHF_CONST_FRA1(a_edgeData[dir],comp),
                          CHF_FRA1(a_cellData, cellcomp),
                          CHF_BOX(cellBox),
                          CHF_CONST_INT(dir));
        } // end loop over directions
    } // end loop over components
}

void EdgeToCell(const FluxBox& a_edgeData, const int a_edgeComp,
                FArrayBox& a_cellData, const int a_cellComp,
                const int a_dir)
{
   Box cellBox = a_cellData.box();
   cellBox &= a_edgeData.box();

   FORT_EDGETOCELL(CHF_CONST_FRA1(a_edgeData[a_dir],a_edgeComp),
                   CHF_FRA1(a_cellData, a_cellComp),
                   CHF_BOX(cellBox),
                   CHF_CONST_INT(a_dir));
}

void EdgeToCell(const FluxBox& a_edgeData, const int a_edgeComp,
                FArrayBox& a_cellData, const int a_cellComp,
                const Box& a_cellBox, const int a_dir)
{

  FORT_EDGETOCELL(CHF_CONST_FRA1(a_edgeData[a_dir],a_edgeComp),
                  CHF_FRA1(a_cellData, a_cellComp),
                  CHF_BOX(a_cellBox),
                  CHF_CONST_INT(a_dir));
}

// max functions

// ----------------------------------------------------------------
void EdgeToCellMax(const LevelData<FluxBox>& a_edgeData,
                   LevelData<FArrayBox>& a_cellData)

{
  // this is just a wrapper around the single-grid version
  DataIterator dit = a_edgeData.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      EdgeToCellMax(a_edgeData[dit()], a_cellData[dit()]);
    }
}

// ----------------------------------------------------------------
void EdgeToCellMax(const FluxBox& a_edgeData,
                   FArrayBox& a_cellData)
{

  // loop over components -- assumption is that in cell-centered
  // data, direction changes faster than component.
  for (int comp = 0; comp<a_edgeData.nComp(); comp++)
    {
      // loop over directions
      for (int dir = 0; dir < SpaceDim; dir++)
        {
          Box cellBox = a_cellData.box();
          cellBox &= a_edgeData.box();
          int cellcomp = SpaceDim*comp + dir;
          FORT_EDGETOCELL(CHF_CONST_FRA1(a_edgeData[dir],comp),
                          CHF_FRA1(a_cellData, cellcomp),
                          CHF_BOX(cellBox),
                          CHF_CONST_INT(dir));
        } // end loop over directions
    } // end loop over components
}

void EdgeToCellMax(const FluxBox& a_edgeData, const int a_edgeComp,
                   FArrayBox& a_cellData, const int a_cellComp,
                   const int a_dir)
{

  const Box& cellBox = a_cellData.box();
  FORT_EDGETOCELLMAX(CHF_CONST_FRA1(a_edgeData[a_dir],a_edgeComp),
                     CHF_FRA1(a_cellData, a_cellComp),
                     CHF_BOX(cellBox),
                     CHF_CONST_INT(a_dir));
}

void EdgeToCellMax(const FluxBox& a_edgeData, const int a_edgeComp,
                   FArrayBox& a_cellData, const int a_cellComp,
                   const Box& a_cellBox, const int a_dir)
{

  FORT_EDGETOCELLMAX(CHF_CONST_FRA1(a_edgeData[a_dir],a_edgeComp),
                     CHF_FRA1(a_cellData, a_cellComp),
                     CHF_BOX(a_cellBox),
                     CHF_CONST_INT(a_dir));
}

#include "NamespaceFooter.H"
