#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#ifndef _EBISLAYOUT_H_
#define _EBISLAYOUT_H_

#include "REAL.H"
#include "IntVect.H"
#include "IntVectSet.H"
#include "VolIndex.H"
#include "FaceIndex.H"
#include "IntVectSet.H"
#include "EBISBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "EBGraph.H"

#include "NamespaceHeader.H"

class EBIndexSpace;
class EBISLayoutImplem;

///
/**
   EBISLayout represents the geometric information
   on a union of rectangles.    This is
   another ref-counted pointer object so copying
   it is cheap but has the ref-counted pointer semantic.
*/
class EBISLayout
{
public:
  ///
  /**
     Access the {\tt EBISBox} associated with the input DataIndex.
     Only constant access is permitted.
  */

  const EBISBox& operator[](const DataIndex& a_index) const;

  ///
  /**
     Sets the maximum level of refinement that this {\tt EBISLayout}
     will have to perform.  Creates and holds new {\tt EBISLayout}s
     at intermediate levels of refinement.
     Default is one (no refinement done).
  */
  void setMaxRefinementRatio(const int& a_maxRefine, const EBIndexSpace* const a_ebisPtr);

  ///
  /**
     Sets the maximum level of coarsening that this {\tt EBISLayout}
     will have to perform.  Creates and holds new {\tt EBISLayout}s
     at intermediate levels of refinement.
     Default is one (no coarsening done).
  */
  void setMaxCoarseningRatio(const int& a_maxCoarsen, const EBIndexSpace* const a_ebisPtr);

  ///
  /**
     Get the nghost used to create this EBISLayout
  */
  int getGhost() const
  {
    return m_nghost;
  }

  ///
  EBISLayout();

  ///
  /**
     Returns the ghosted layout associated with the ebisl.
   */
  const BoxLayout& getGrownLayout() const;

  ///
  /**
     Returns the disjoint layout that was sent in.
   */
  const DisjointBoxLayout& getDisjointLayout() const;

  ///
  ~EBISLayout();

  ///
  /**
     Returns the index of the VoF corresponding to coarsening
     the input VoF by the input ratio.  It is an error if
     the ratio is greater than the maximum coarsening ratio
     or if the vof does not exist at the input data index.
   */
  VolIndex coarsen(const VolIndex& a_vof,
                   const int& a_ratio,
                   const DataIndex& a_datInd) const;

  ///
  /**
     Returns the indices of the VoFs corresponding to refining
     the input VoF by the input ratio.  It is an error if
     the ratio is greater than the maximum refinement ratio
     or if the vof does not exist at the input data index.
  */
  Vector<VolIndex> refine(const VolIndex& a_vof,
                          const int& a_ratio,
                          const DataIndex& a_datInd) const;

  ///
  /**
     Returns the indices of the Faces corresponding to refining
     the input Face by the input ratio.  It is an error if
     the ratio is greater than the maximum refinement ratio
     or if the face does not exist at the input data index.
  */
  Vector<FaceIndex> refine(const FaceIndex& a_face,
                           const int& a_ratio,
                           const DataIndex& a_datInd) const;

  ///
  int getMaxCoarseningRatio() const;

  ///
  int getMaxRefinementRatio() const;

  ///
  void define(const ProblemDomain& a_domain,
              const DisjointBoxLayout& a_grids,
              const int& a_nghost,
              const LevelData<EBGraph>& a_graph,
              const LevelData<EBData> & a_data);

  bool isDefined() const;

  //multifluid things.  tamper not with these unless you REALLY know what you are doing.
  const EBIndexSpace* getEBIS() const;
  void setEBIS(const EBIndexSpace* const a_ebisPtr);
private:
  ///
  RefCountedPtr<EBISLayoutImplem> m_implem;
  int                             m_nghost;
};

///
/**
   EBISLayoutImplem represents the geometric information
   on a union of rectangles.  This is not really part
   of the public interface.  Users should use EBISLayout.
*/
class EBISLayoutImplem
{
public:

  ///
  const EBISBox& operator[](const DataIndex&) const;

  ///
  EBISLayoutImplem();

  ///
  ~EBISLayoutImplem();

  ///
  const BoxLayout& getGrownLayout() const
  {
    return m_blGhostDom;
  }

  const DisjointBoxLayout& getDisjointLayout() const
  {
    return m_dblInputDom;
  }

  ///
  /**
     Returns the index of the VoF corresponding to coarsening
     the input VoF by the input ratio.  It is an error if
     the ratio is greater than the maximum coarsening ratio
     or if the vof does not exist at the input data index.
   */
  VolIndex coarsen(const VolIndex& a_vof,
                   const int& a_ratio,
                   const DataIndex& a_datInd) const;

  ///
  /**
     Returns the index of the face corresponding to coarsening
     the input face by the input ratio.  It is an error if
     the ratio is greater than the maximum coarsening ratio
     or if the face does not exist at the input data index.
   */
  FaceIndex coarsen(const FaceIndex& a_face,
                    const int&       a_ratio,
                    const DataIndex& a_datInd) const;

  ///
  /**
     Returns the indices of the VoFs corresponding to refining
     the input VoF by the input ratio.  It is an error if
     the ratio is greater than the maximum refinement ratio
     or if the vof does not exist at the input data index.
  */
  Vector<VolIndex> refine(const VolIndex&  a_vof,
                          const int&       a_ratio,
                          const DataIndex& a_datInd) const;

  ///
  /**
     Returns the indices of the Faces corresponding to refining
     the input Face by the input ratio.  It is an error if
     the ratio is greater than the maximum refinement ratio
     or if the face does not exist at the input data index.
  */
  Vector<FaceIndex> refine(const FaceIndex&  a_face,
                           const int&       a_ratio,
                           const DataIndex& a_datInd) const;

  ///
  /**
     Sets the maximum level of refinement that this {\tt EBISLayout}
     will have to perform.  Creates and holds new {\tt EBISLayout}s
     at intermediate levels of refinement.
     Default is one (no refinement done).
  */
  void setMaxRefinementRatio(const int& a_maxRefine, const EBIndexSpace* const a_ebisPtr);

  ///
  /**
     Sets the maximum level of coarsening that this {\tt EBISLayout}
     will have to perform.  Creates and holds new {\tt EBISLayout}s
     at intermediate levels of refinement.
     Default is one (no coarsening done).
  */
  void setMaxCoarseningRatio(const int& a_maxCoarsen, const EBIndexSpace* const a_ebisPtr);

  ///
  int getMaxCoarseningRatio() const;

  ///
  int getMaxRefinementRatio() const;

  ///
  void define(const ProblemDomain& a_domain,
              const DisjointBoxLayout& a_grids,
              const int& a_nghost,
              const LevelData<EBGraph>& a_graph,
              const LevelData<EBData> & a_data);

  ///
  static void setVerbose(bool a_verbose);

  bool isDefined() const {return m_defined;}

  const EBIndexSpace* getEBIS() const
  {
    return m_ebisPtr;
  }
  void setEBIS(const EBIndexSpace* const a_ebisPtr)
  {
    m_ebisPtr = a_ebisPtr;
  }
private:

  const EBIndexSpace* m_ebisPtr;
  ///
  BoxLayout m_blGhostDom;

  ///
  DisjointBoxLayout m_dblInputDom;

  ///
  Vector<EBISLayout> m_fineLevels;

  ///
  Vector<EBISLayout> m_coarLevels;

  ///
  int m_nghost;

  ///
  int m_maxCoarseningRatio;

  ///
  int m_maxRefinementRatio;

  ///
  ProblemDomain m_domain;

  ///
  LayoutData<EBISBox> m_ebisBoxes;

  bool m_defined;

private:
  static bool s_verbose;
  void operator=(const EBISLayoutImplem& ebiin)
  {;}
  EBISLayoutImplem(const EBISLayoutImplem& ebiin)
  {;}

};

#include "NamespaceFooter.H"
#endif
