#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PETSCCOMPGRID_H_
#define _PETSCCOMPGRID_H_
#ifdef CH_USE_PETSC
#include "Stencil.H"
#include "BCFunc.H"
#include "DisjointBoxLayout.H"
#include "BoxIterator.H"
#include "FourthOrderInterpStencil.H"
#include <petsc-private/pcimpl.h>

#include "NamespaceHeader.H"

typedef enum {GHOST=-1,FINE_COVERED=-2,DIRI_BC=-3,
              NEUM_BC=-4,ANY_BC=-5,UNKNOWN=-6} GID_type; 

std::ostream& operator<< (std::ostream& os, GID_type);

//! \class PetscCompGrid
//! This base class organizes the construction of a PETSc matrix, with an AMR hierarchy
class PetscCompGrid
{
public: 

  //! Base class constructor. Called by all subclass constructors.
  PetscCompGrid(int a_dof) : m_gid0(0),m_mat(0),m_CFStencilRad(2),m_writeMatlab(false),
			     m_verbose(0),m_dof(a_dof)
  {
#if defined(PETSC_USE_LOG)
    PetscLogEventRegister("createMatrix",PETSC_VIEWER_CLASSID,&m_event0);
    PetscLogEventRegister("PetscCompGrid",PETSC_VIEWER_CLASSID,&m_event1);
#endif
  }
  
  virtual void define( const ProblemDomain &a_cdomain,
                       Vector<DisjointBoxLayout> &a_grids, 
                       Vector<int> &a_refratios, 
                       BCHolder a_bc,
                       const RealVect &a_cdx,
                       int a_numLevels=-1, int a_ibase=0);
  
  //! Destructor.
  virtual ~PetscCompGrid();
  virtual void clean();

  Mat getMatrix() const { return m_mat; }
  void setMatlab(bool b=true) {m_writeMatlab = b;}
  void setVerbose(int a_v){ m_verbose=a_v;}

  // main Mat creation routines
  PetscErrorCode createMatrix();
  // solver uses these, not used in base class, used in app code
  PetscErrorCode putChomboInPetsc( const Vector<LevelData<FArrayBox> * > &rhs, Vec b )const;
  PetscErrorCode putPetscInChombo( Vec b, Vector<LevelData<FArrayBox> * > &rhs )const;
  // derived type needs to tell me how many ghosts it needs
  virtual IntVect getGhostVect()const = 0;
  virtual void addExtraCovered(GID_type,int,const DataIndex&,BaseFab<PetscInt>&) {} // default no op
protected:
  // core stencil ops - base class has basic implementations except for equations
  virtual void createOpStencil(IntVect,int,const DataIndex&,StencilTensor &) = 0;
  virtual void applyBCs(IntVect,int,const DataIndex&,Box,StencilTensor &);
  virtual void InterpToFine(IntVect,int,const DataIndex&,StencilTensor &);
  virtual void InterpToCoarse(IntVect,int,const DataIndex&,StencilTensor &);
  // helper for FourthOrderInterpStencil
  IntVect getCFStencil(const ProblemDomain &a_cdom, const IntVect a_ivc);
  // at end add stencil to matrix - done
  PetscErrorCode AddStencilToMat(IntVect,int,const DataIndex&,StencilTensor &, Mat);
  // utils
  void NodeDefine(StenScalarNode &a_node, IntVect a_iv, int a_lev, Real a_val)
  {
    a_node.first.setIV(a_iv);
    a_node.first.setLevel(a_lev);
    a_node.second.setValue(a_val);
  }
  // data
  Vector<ProblemDomain> m_domains;
  Vector<DisjointBoxLayout> m_grids;
  Vector<int> m_refRatios;
  Vector<RealVect> m_dxs;
public:
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt> > > > m_GIDs;
protected:
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt> > > > m_crsSupportGIDs;
  Vector<RefCountedPtr<LevelData<BaseFab<PetscInt> > > > m_fineCoverGIDs;
public:
  int m_gid0;
protected:
  Mat m_mat;
  BCHolder m_bc;
  // all possible stencils, on box(-m_CFStencilRad:+m_CFStencilRad)^SpaceDim
  int m_CFStencilRad;
  BaseFab<FourthOrderInterpStencil*> m_FCStencils;
  bool m_writeMatlab;
  int m_verbose;
  const int m_dof;
#if defined(PETSC_USE_LOG)
  PetscLogEvent m_event0;
  PetscLogEvent m_event1;
#endif
private:
  // Disallowed for all the usual reasons
  void operator=(const PetscCompGrid& a_input)
  {
    MayDay::Error("invalid operator");
  }
  // Disallowed for all the usual reasons
  PetscCompGrid(const PetscCompGrid& a_input):m_dof(0)
  {
    MayDay::Error("invalid operator");
  }
};

class CompBC : public BCFunction
{
public:
  CompBC(int a_nSource, IntVect a_nGhosts);
  ~CompBC();
  void define(int a_nSource, IntVect a_nGhosts);

  virtual void createCoefs() = 0;
  virtual void operator()( FArrayBox&           a_state,
                           const Box&           a_valid,
                           const ProblemDomain& a_domain,
                           Real                 a_dx,
                           bool                 a_homogeneous) = 0;
 
  IntVect nGhosts()const{return m_nGhosts;}
  int nSources()const{return m_nSources;}
  PetscReal getCoef(int a_iSrc, int a_iGhost=0);
protected:
  PetscReal *m_Rcoefs;
  IntVect    m_nGhosts;
  int        m_nSources; // degreee is m_nSources-1
  bool       m_isDefined;
};

class ConstDiriBC : public CompBC
{
public:
  ConstDiriBC(int a_nSource=1, IntVect a_nGhosts=IntVect::Unit) : CompBC(a_nSource,a_nGhosts) {}
  virtual void createCoefs();
  virtual void operator()( FArrayBox&           a_state,
                           const Box&           a_valid,
                           const ProblemDomain& a_domain,
                           Real                 a_dx,
                           bool                 a_homogeneous);
};

#include "NamespaceFooter.H"

#endif
#endif
