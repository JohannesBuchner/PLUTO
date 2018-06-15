#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CELLTOEDGE_H_
#define _CELLTOEDGE_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

///
void CellToEdge(const LevelData<FArrayBox>& a_cellData,
                LevelData<FluxBox>& a_edgeData);

///
void CellToEdge(const FArrayBox& a_cellData,
                FluxBox& a_edgeData);

/// this one does a single direction
void CellToEdge(const FArrayBox& a_cellData,
                FArrayBox& a_edgeData, const int a_dir);

/// this one does a single direction and component
void CellToEdge(const FArrayBox& a_cellData, const int a_cellComp,
                FArrayBox& a_edgeData, const int a_edgeComp,
                const int a_dir);

#include "NamespaceFooter.H"

#endif
