#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves fri aug 17, 2001

#ifndef _EBMGINTERP_H_
#define _EBMGINTERP_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "EBIndexSpace.H"
#include "EBStencil.H"
#include "EBLevelGrid.H"
#include "EBConstantCFInterp.H"
#include "NamespaceHeader.H"

/// Piecewise constant interpolation.
/**
   Finedata += coarse data.
  */
class EBMGInterp
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBMGInterp();

  ///
  ~EBMGInterp();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  EBMGInterp(const DisjointBoxLayout& dblFine,
             const DisjointBoxLayout& dblCoar,
             const EBISLayout& ebislFine,
             const EBISLayout& ebislCoar,
             const ProblemDomain& domainCoar,
             const int& nref,
             const int& nvar,
             const EBIndexSpace* ebisPtr,
             const IntVect& ghostCellsPhi,
             const bool& layoutChanged = true,
             const bool& doLinear  = false);


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
     at each VoF.
  */
  void define(const DisjointBoxLayout& dblFine,
              const DisjointBoxLayout& dblCoar,
              const EBISLayout& ebislFine,
              const EBISLayout& ebislCoar,
              const ProblemDomain& domainCoar,
              const int& nref,
              const int& nvar,
              const EBIndexSpace* ebisPtr,
              const IntVect& ghostCellsPhi,
              const bool& layoutChanged = true,
              const bool& doLinear  = false);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Piecewise constant interpolation function
     FineData += coarseData
  */
  void
  pwcInterp(LevelData<EBCellFAB>&       a_fineData,
            const LevelData<EBCellFAB>& a_coarseData,
            const Interval&             a_variables);


  ///
  /**
     Piecewise linear interpolation function
     FineData += coarseData + slope*dist
     only works if you turned on useLinear in define function.
     Ghost cells for the coarse data must be filled coming into the function.
  */
  void
  pwlInterp(LevelData<EBCellFAB>&       a_fineData,
            const LevelData<EBCellFAB>& a_coarseData,
            const Interval&             a_variables);

  ///
  /**
     Piecewise linear interpolation function
     FineData += coarseData + slope*dist
     only works if you turned on useLinear in define function.
     no copies here
     Ghost cells for the coarse data must be filled coming into the function.
  */
  void
  pwlInterpMG(LevelData<EBCellFAB>&       a_fineData,
              const LevelData<EBCellFAB>& a_coarseData,
              const Interval&             a_variables);

  ///
  /**
     Piecewise constant interpolation function for unchanged layouts (no agglomeration)
     there are no copies here
     FineData += coarseData
  */
  void
  pwcInterpMG(LevelData<EBCellFAB>&       a_fineData,
              const LevelData<EBCellFAB>& a_coarseData,
              const Interval&             a_variables);


protected:
  void
  defineLinearStencil(const DataIndex       & a_dit,
                      const Vector<VolIndex>& a_allFineVoFs,
                      const EBISLayout      & a_ebislStenFine,
                      const EBISLayout      & a_ebislStenCoar,
                      const Box             & a_boxFine,
                      const Box             & a_boxCoar);

  void
  defineConstantStencil(const DataIndex       & a_dit,
                        const Vector<VolIndex>& a_allFineVoFs,
                        const EBISLayout      & a_ebislStenFine,
                        const EBISLayout      & a_ebislStenCoar,
                        const Box             & a_boxFine,
                        const Box             & a_boxCoar);

  void fillGhostCellsPWC(LevelData<EBCellFAB>& a_data,
                         const EBISLayout&     a_ebisl,
                         const ProblemDomain&  a_dom);
  void
  pwcInterpFAB(EBCellFAB&       a_refCoar,
               const Box&       a_coarBox,
               const EBCellFAB& a_fine,
               const DataIndex& a_datInd,
               const Interval&  a_variables) const;


  void
  pwlInterpFAB(EBCellFAB&       a_refCoar,
               const Box&       a_coarBox,
               const EBCellFAB& a_fine,
               const DataIndex& a_datInd,
               const Interval&  a_variables) const;

  void setDefaultValues();

  void defineStencils();

  bool m_isDefined;
  IntVect m_ghost;
  DisjointBoxLayout m_coarGrids;
  DisjointBoxLayout m_fineGrids;
  ProblemDomain m_coarDomain;
  ProblemDomain m_fineDomain;
  EBISLayout m_coarEBISL;
  EBISLayout m_fineEBISL;

  int m_refRat;
  int m_nComp;

  bool m_doLinear;
  bool m_layoutChanged;
  bool m_coarsenable;

  //LevelData<EBCellFAB> m_refinedCoarseData;
  //if coarsenable, these are the coarsened fine
  // otherwise, refined coarse.
  //actual buffer is no longer member data
  //not defined if !m_layoutChanged
  DisjointBoxLayout m_buffGrids;
  EBISLayout        m_buffEBISL;

  //the following copiers only get defined where needed.
  Copier m_copierFtoRC;//this goes from m_fineGrids to m_refinedCoarseGrids
  Copier m_copierRCtoF;//this goes from m_refinedCoarseGrids to m_fineGrids


  //stencils for prolongation
  LayoutData<RefCountedPtr<EBStencil> >  m_interpEBStencil;

  //stencils for prolongation
  LayoutData<RefCountedPtr<EBStencil> >  m_linearEBStencil;

private:
    //disallowed for all the usual reasons
  EBMGInterp(const EBMGInterp& ebcin)
  {
    MayDay::Error("ebmgi 2 invalid operator");
  }
  void operator=(const EBMGInterp& fabin)
  {
    MayDay::Error("ebmgi 3 invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
