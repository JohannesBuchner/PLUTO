#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// EdgeToCell.H
// Dan Martin, Fri, Jan 14, 2000

#ifndef _EDGETOCELL_H_
#define _EDGETOCELL_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

/// averages cell-centered cellData to edge-centered edgeData
void EdgeToCell(const LevelData<FluxBox>& a_edgeData,
                LevelData<FArrayBox>& a_cellData);

/// averages cellData->edgeData (all directions) for one FAB
/** averages cell-centered cellData to SpaceDim-dimensioned
    EdgeData -- assumes all boundary conditions already set.
    Does this for all grids. */
void EdgeToCell(const FluxBox& a_edgeData,
                FArrayBox& a_cellData);

/// averages cellData->edgeData in 1 direction, for one component in a FAB
void EdgeToCell(const FluxBox& a_edgeData, const int a_edgeComp,
                FArrayBox& a_cellData, const int a_cellComp,
                const int a_dir);

/// averages cellData->edgeData in 1 direction, for one component in a FAB, on a given cell box
void EdgeToCell(const FluxBox& a_edgeData, const int a_edgeComp,
                FArrayBox& a_cellData, const int a_cellComp,
                const Box& a_cellBox, const int a_dir);

// max functions
///  cell-centered cellData from edge-centered edgeData using max
/** instead of averaging faces to cells, this function puts the max of
    the two opposing face-centered values into each cell.
*/
void EdgeToCellMax(const LevelData<FluxBox>& a_edgeData,
                   LevelData<FArrayBox>& a_cellData);

///  edgeData->cellData using max (all directions) for one FAB
/** fills cell-centered cellData from EdgeData using max
    -- assumes all boundary conditions already set.
    Does this for all grids. */
void EdgeToCellMax(const FluxBox& a_edgeData,
                   FArrayBox& a_cellData);

/// fills cellData using edgeData in 1 direction using max, for one component in a FAB
void EdgeToCellMax(const FluxBox& a_edgeData, const int a_edgeComp,
                   FArrayBox& a_cellData, const int a_cellComp,
                   const int a_dir);

/// fills cellData->edgeData in 1 direction, for one component in a FAB, on a given cell box using max
void EdgeToCellMax(const FluxBox& a_edgeData, const int a_edgeComp,
                   FArrayBox& a_cellData, const int a_cellComp,
                   const Box& a_cellBox, const int a_dir);

#include "NamespaceFooter.H"

#endif
