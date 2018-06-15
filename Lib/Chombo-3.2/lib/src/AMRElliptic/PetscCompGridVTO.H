#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCCOMPGRIDVTO_H_
#define _PETSCCOMPGRIDVTO_H_
#ifdef CH_USE_PETSC
#include "PetscCompGrid.H"
#include "NamespaceHeader.H"

//! \class PetscCompGridVTO
//! This base class organizes the construction of a PETSc matrix, 
//! and solve with an AMR 
class PetscCompGridVTO : public PetscCompGrid
{
public:
  //! Base class constructor. Called by all subclass constructors.
  PetscCompGridVTO(Real a_al=0., Real a_beta=1.0) : 
    PetscCompGrid(CH_SPACEDIM), m_alpha(a_al), m_beta(a_beta)
  {
  }
  virtual void clean();

  void defineCoefs( Real a_alpha, Real a_beta, 
                    Vector<RefCountedPtr<LevelData<FluxBox> > > &a_mu,
                    Vector<RefCountedPtr<LevelData<FluxBox> > > &a_lambda,
                    Vector<RefCountedPtr<LevelData<FArrayBox> > > &a_a)
  {
    m_alpha = a_alpha;
    m_beta = a_beta;
    m_a = a_a;
    m_eta = a_mu;
    m_lamb = a_lambda;
  }

  //! Destructor.
  virtual ~PetscCompGridVTO()
  {   
    clean();
  }
  virtual IntVect getGhostVect()const
  {
    return IntVect::Unit;
  }
  Real getAlpha()const{return m_alpha;}
  Real getBeta()const{return m_beta;}
  void setAlpha(Real a_b) {m_alpha = a_b;}
  void setBeta(Real a_b) {m_beta = a_b;}
protected:
  virtual void createOpStencil(IntVect,int,const DataIndex&,StencilTensor&);
  Real              m_alpha;
  Real              m_beta;
  //coefficient storage. 
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_eta;
  Vector<RefCountedPtr<LevelData<FluxBox> > > m_lamb;
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_a;
};

#include "NamespaceFooter.H"
#endif
#endif
