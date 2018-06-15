#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCCOMPGRIDPOIS_H_
#define _PETSCCOMPGRIDPOIS_H_
#ifdef CH_USE_PETSC

#include "PetscCompGrid.H"
#include "NamespaceHeader.H"

//! \class PetscCompGridPois
//! This base class organizes the construction of a PETSc matrix, 
//! and solve with an AMR hierarchy
#define COMP_POIS_DOF 1
class PetscCompGridPois : public PetscCompGrid
{
public:
  //! Base class constructor. Called by all subclass constructors.
  PetscCompGridPois(Real a_al=0., Real a_be=-1.0, int a_order=2) : 
    PetscCompGrid(COMP_POIS_DOF), m_alpha(a_al), m_beta(a_be), m_order(a_order), m_cornerStencil(false)
  {
  }
  virtual void clean();

  //! Destructor.
  virtual ~PetscCompGridPois()
  {   
    clean();
  }
  bool isCornerStencil() const {return m_cornerStencil;}
  void setCornerStencil(bool a_b=true) {m_cornerStencil = a_b;}
  virtual IntVect getGhostVect()const
  {
    return (m_order==2) ? IntVect::Unit : 2*IntVect::Unit;
  }
  Real getAlpha()const{return m_alpha;}
  Real getBeta()const{return m_beta;}
  void setAlpha(Real a_b) {m_alpha = a_b;}
  void setBeta(Real a_b) {m_beta = a_b;}
  int getOrder()const{return m_order;}
  void setOrder(int a_b) {m_order = a_b;}
protected:
  virtual void createOpStencil(IntVect,int,const DataIndex&,StencilTensor &);
  Real              m_alpha;
  Real              m_beta;
  int               m_order;
  bool              m_cornerStencil;
};

#include "NamespaceFooter.H"
#endif
#endif
