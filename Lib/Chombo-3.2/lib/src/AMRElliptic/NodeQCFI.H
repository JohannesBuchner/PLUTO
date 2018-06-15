#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// NodeQCFI.H
// petermc, Mon, Apr 23, 2001
// petermc, 1 Nov 2001, redid this so that it uses NodeQCFI2
#ifndef _NODEQCFI_H_
#define _NODEQCFI_H_

#include <cmath>
#include <cstdlib>
#include "BaseFab.H"
#include "NodeFArrayBox.H"
#include "NodeCFIVS.H"
#include "LevelData.H"
#include "NodeQuadCFInterp2.H"
#include "NodeBCFunc.H"
#include "NamespaceHeader.H"

///  Class to interpolate quadratically at coarse/fine interface.
class NodeQCFI
/**

*/
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.  User must subsequently call define().
   */
  NodeQCFI();

  ///
  /** Constructor calls setDefaultValues() and then
      calls define() with the same arguments.
  */
  NodeQCFI(const DisjointBoxLayout& a_grids,
           Real a_dx,
           const ProblemDomain& a_domain,
           const LayoutData<NodeCFIVS>* const a_loCFIVS,
           const LayoutData<NodeCFIVS>* const a_hiCFIVS,
           int a_refToCoarse,
           NodeBCFunc               a_bc,
           int a_interpolationDegree = 2,
           int a_ncomp = 1,
           bool a_verbose = false);

  ///
  /** Constructor calls setDefaultValues() and then
      calls define() with the same arguments.
  */
  NodeQCFI(const DisjointBoxLayout& a_grids,
           Real a_dx,
           const Box& a_domain,
           const LayoutData<NodeCFIVS>* const a_loCFIVS,
           const LayoutData<NodeCFIVS>* const a_hiCFIVS,
           int a_refToCoarse,
           NodeBCFunc               a_bc,
           int a_interpolationDegree = 2,
           int a_ncomp = 1,
           bool a_verbose = false);

  ///
  /** Destructor.
   */
  ~NodeQCFI();

  ///
  /** Full define function.  Makes all coarse-fine
      information and sets internal variables.  The current level
      is taken to be the fine level.
  */
  void define(/// CELL-centered grids at this level
              const DisjointBoxLayout& a_grids,
              /// mesh spacing at this level
              Real a_dx,
              /// CELL-centered physical domain at this level
              const ProblemDomain& a_domain,
              /// pointer to object storing coarse/fine interface nodes
              const LayoutData<NodeCFIVS>* const a_loCFIVS,
              /// pointer to object storing coarse/fine interface nodes
              const LayoutData<NodeCFIVS>* const a_hiCFIVS,
              /// refinement ratio between this and next coarser level
              int a_refToCoarse,
              NodeBCFunc               a_bc,
              /// degree of interpolation; 1 for (bi)linear, 2 for (bi)quadratic
              int a_interpolationDegree = 2,
              /// number of components of data
              int a_ncomp = 1,
              /// verbose output flag
              bool a_verbose = false);

  ///
  /** Full define function.  Makes all coarse-fine
      information and sets internal variables.  The current level
      is taken to be the fine level.
  */
  void define(/// CELL-centered grids at this level
              const DisjointBoxLayout& a_grids,
              /// mesh spacing at this level
              Real a_dx,
              /// CELL-centered physical domain at this level
              const Box& a_domain,
              /// pointer to object storing coarse/fine interface nodes
              const LayoutData<NodeCFIVS>* const a_loCFIVS,
              /// pointer to object storing coarse/fine interface nodes
              const LayoutData<NodeCFIVS>* const a_hiCFIVS,
              /// refinement ratio between this and next coarser level
              int a_refToCoarse,
              NodeBCFunc               a_bc,
              /// degree of interpolation, 1 for (bi)linear, 2 for (bi)quadratic
              int a_interpolationDegree = 2,
              /// number of components of data
              int a_ncomp = 1,
              /// verbose output flag
              bool a_verbose = false);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object was created with the defining
      constructor or if define() has been called.
  */
  bool isDefined() const;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Set whether to give output.  Default is <tt>false</tt>.
   */
  void setVerbose( bool a_verbose );

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Coarse / Fine (inhomogeneous) interpolation operator.
      Fill the nodes of <i>a_phi</i> on the coarse/fine interface with
      interpolated data from <i>a_phiCoarse</i>.
  */
  void coarseFineInterp(/// data at this level
                        LevelData<NodeFArrayBox>& a_phiFine,
                        /// data at the next coarser level
                        const LevelData<NodeFArrayBox>& a_phiCoarse,
                        /// whether to apply inhomogeneous physical boundary conditions?
                        bool a_inhomogeneous);

  /*@}*/

protected:

  NodeBCFunc               m_bc;
  /** CELL-centered grids at the current level (the finer level)
   */
  DisjointBoxLayout m_grids;

  /** number of components of data, needed for setting size of work array
   */
  int m_ncomp;

  /** refinement ratio between this and the next coarser level
   */
  int m_refToCoarse;

  /** the number of coarsenings to be done:
      this is log2(<i>m_refToCoarse</i>)
  */
  int m_coarsenings;

  /** has full define function been called?
   */
  bool m_isDefined;

  /** mesh spacing at this (fine) level
   */
  Real m_dx;

  /** mesh spacing at coarsest level refined by 2
   */
  Real m_dxPenultimate;

  /** CELL-centered physical domain of coarsest level refined by 2
   */
  ProblemDomain m_domainPenultimate;

  /** data at intermediate coarsening levels;
      Vector length <i>m_coarsenings</i>-1.
  */
  Vector< LevelData<NodeFArrayBox>* > m_inter;

  /** Vector of interpolating objects, of length <i>m_coarsenings</i>.
      <i>m_qcfi2</i>[i] averages refined versions of the SAME grids if i > 1,
      and in general DIFFERENT grids if i == 0.
  */
  Vector< NodeQuadCFInterp2* > m_qcfi2;

  /** pointers to objects storing coarse/fine interface nodes
      between levels of refinement, including intermediate ones
  */
  Vector< LayoutData<NodeCFIVS>* > m_loCFIVScoarser;

  /** pointers to objects storing coarse/fine interface nodes
      between levels of refinement, including intermediate ones
  */
  Vector< LayoutData<NodeCFIVS>* > m_hiCFIVScoarser;

  /** if <tt>true</tt>, print out extra information
   */
  bool m_verbose ;

private:
  //internally useful functions

  void clearMemory();

  void setDefaultValues();
};

#include "NamespaceFooter.H"
#endif
