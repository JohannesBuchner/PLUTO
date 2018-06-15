#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves fri, sept 14, 2001
// mfbarad made it eb aware in June, 2005
// and on and on...

#ifndef _EBQUADCFINTERP_H_
#define _EBQUADCFINTERP_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "QuadCFInterp.H"
#include "BaseIVFAB.H"
#include "CornerCopier.H"
#include "VoFIterator.H"
#include "EBCFData.H"
#include "NamespaceHeader.H"

///
/**
   Fills one layer of fine level ghost cells
   by piecewise linear interpolation
   from the coarse level.
 */
class EBQuadCFInterp: public QuadCFInterp
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBQuadCFInterp();

  ///
  virtual ~EBQuadCFInterp();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.  \\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF. \\
     Set doEBCFCrossing=false if you are sure that tall
     irregular cells are on the finest grid.
  */
  EBQuadCFInterp(const DisjointBoxLayout& a_dblFine,
                 const DisjointBoxLayout& a_dblCoar,
                 const EBISLayout&        a_ebislFine,
                 const EBISLayout&        a_ebislCoar,
                 const ProblemDomain&     a_domainCoar,
                 const int&               a_nref,
                 const int&               a_nvar,
                 const LayoutData<IntVectSet>&  a_cfivs,
                 const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance(),
                 bool a_doEBCFCrossing = true);

  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF. \\
     Set doEBCFCrossing=false if you are sure that tall
     irregular cells are on the finest grid.
  */
  void define(const DisjointBoxLayout& a_dblFine,
              const DisjointBoxLayout& a_dblCoar,
              const EBISLayout&        a_ebislFine,
              const EBISLayout&        a_ebislCoar,
              const ProblemDomain&     a_domainCoar,
              const int&               a_nref,
              const int&               a_nvar,
              const LayoutData<IntVectSet>&  a_cfivs,
              const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance(),
              bool a_doEBCFCrossing = true);


  ///
  /**
     Interpolate the fine data from the coarse data
     over the intersection
     of the fine layout with the refined coarse layout. \\
     {\bf Arguments:}\\
     coarDataOld: The data over the coarse layout.\\
     fineData: The data over the fine layout.  \\
     Fine and coarse data must
     have the same number of variables.\\
     variables:  The variables to interpolate.  Those not
     in this range will be left alone.   This range of variables
     must be in both the coarse and fine data.
     {\bf This routine assumes that coarse data that lives under finer
     grids has been coarsened from the finer data to O(h^3).  If embedded boundaries
     cross coarse fine boundaries and this is not the case, O(1) errors will happen}
  */
  void
  interpolate(LevelData<EBCellFAB>&       a_fineData,
              const LevelData<EBCellFAB>& a_coarData,
              const Interval&             a_variables,
              bool a_doOnlyRegularInterp = false);

  bool isDefined() const;

  ///
  /**
     Use the stencils to actually do EB aware interpolation
  */
  void interpEBCFCrossing(LevelData<EBCellFAB>&       a_fineData,
                          const LevelData<EBCellFAB>& a_coarData,
                          const Interval&             a_variables);

  ///
  /**
     Use the corner-stencils to actually do EB aware interpolation
     Use coarse value when stencil is empty.
  */
  void interpEBCFCorners(LevelData<EBCellFAB>&       a_fineData,
                         const LevelData<EBCellFAB>& a_coarData,
                         const Interval&             a_variables);

  //need this because ebtensorcfi needs it
  RefCountedPtr<EBCFData> getEBCFData() const
  {return m_ebcfdata;}

protected:

  ///
  /**
     Builds the IntVectSets for where we need to do
     EB aware CF interpolation.
     Returns a bool that indicates if we need to do EB aware CF interp
  */
  bool getEBCFIVS(const LayoutData<IntVectSet>& a_cfivs);

  ///
  /**
     Builds a stencil for each vof that needs an
     EB aware CF interpolation
  */
  void buildEBCFStencils();

  ///
  /**
     Builds a stencil for each corner-vof that needs an
     EB aware CF interpolation
  */
  void buildEBCFCornerStencils(const LayoutData<IntVectSet>& a_cfivs);



  RefCountedPtr<EBCFData> m_ebcfdata;


  LevelData<EBCellFAB> m_ebBufferCoarsenedFine;

  LayoutData<BaseIVFAB<VoFStencil> > m_coarStencilLo[SpaceDim];
  LayoutData<BaseIVFAB<VoFStencil> > m_coarStencilHi[SpaceDim];

  //These are the stencils that interpolate with phistar to get ghost val
  //phiStar is indexed at the same point as the ghost phi
  //(UNLIKE in regular, non eb quadcfinterp) but lives between
  //coarse cell centers (just like in non-eb quadcfinterp)
  LayoutData<BaseIVFAB<VoFStencil> > m_fineStencilLo[SpaceDim];
  LayoutData<BaseIVFAB<VoFStencil> > m_fineStencilHi[SpaceDim];
  LayoutData<BaseIVFAB<VoFStencil> > m_stencilCorners;
  LayoutData<BaseIVFAB<VoFStencil> > m_stencilEdges;

  //fine index space fine layout
  LayoutData<IntVectSet> m_ebcfivsLo[SpaceDim];
  LayoutData<IntVectSet> m_ebcfivsHi[SpaceDim];

  // Copier object to handle corner cells
  CornerCopier m_cornerCopier;

  bool m_doEBCFCrossing;
  int m_refRat;
  static IntVect s_ivDebFine;
  static IntVect s_ivDebCoar;
  ProblemDomain m_domainFine, m_domainCoar;
private:

  //disallowed for all the usual reasons
  EBQuadCFInterp(const EBQuadCFInterp& ebcin)
  {
    MayDay::Error("ebqcfi 2 invalid operator");
  }
  void operator=(const EBQuadCFInterp& fabin)
  {
    MayDay::Error("ebqcfi 3 invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
