#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BCFUNC_H_
#define _BCFUNC_H_

#include "IntVect.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "RefCountedPtr.H"
#include "DataIndex.H"
#include "NamespaceHeader.H"

///
/**
   function interface for ghost cell boundary conditions
   of AMRPoissonOp.   If you are using Neumann or Dirichlet
   boundary conditions,   it is easiest to use the functions
   provided.
 */
typedef void(*BCFunc)(FArrayBox&           a_state,
                      const Box&           a_valid,
                      const ProblemDomain& a_domain,
                      Real                 a_dx,
                      bool                 a_homogeneous);

/// For use with pure periodic BC
void doNothingBC(FArrayBox&           a_state,
                 const Box&           a_valid,
                 const ProblemDomain& a_domain,
                 Real                 a_dx,
                 bool                 a_homogeneous);

//! \class BCFunction
//! This base class represents a boundary condition for an elliptic or parabolic
//! partial differential equation.
class BCFunction
{
public:

  //! Base class constructor. Called by all subclass constructors.
  BCFunction()
  {
  }

  //! Destructor.
  virtual ~BCFunction()
  {
  }

  //! Computes values of a solution, \f$\phi\f$, on ghost cells. These ghost values
  //! impose the boundary condition represented by the BCFunction object.
  //! \param a_state The values of \f$\phi\f$ on the box given by \a a_valid.
  //! \param a_valid The box on which the boundary condition is to be imposed.
  //! \param a_domain The problem domain on which this boundary condition is imposed.
  //! \param a_dx The grid spacing.
  //! \param a_homogeneous If set to true, ghost values are computed for a homogeneous
  //!                      boundary condition. This is useful for multigrid solves.
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous) = 0;

  //! Computes the values of \f$\phi\f$ on ghost cells specifying a data index.
  //! By default, this calls the version of the method without a DataIndex.
  //! \param a_state The values of \f$\phi\f$ on the box given by \a a_valid.
  //! \param a_valid The box on which the boundary condition is to be imposed.
  //! \param a_domain The problem domain on which this boundary condition is imposed.
  //! \param a_dx The grid spacing.
  //! \param a_index A DataIndex that can be used in the calculation.
  //! \param a_homogeneous If set to true, ghost values are computed for a homogeneous
  //!                      boundary condition. This is useful for multigrid solves.
  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          const DataIndex&     a_index,
                          bool                 a_homogeneous)
  {
    operator()(a_state, a_valid, a_domain, a_dx, a_homogeneous);
  }

  //! Sets the time associated with this boundary condition. By default, this
  //! method does nothing, so this aspect of BCFunction can be ignored for
  //! time-independent boundary conditions.
  //! \param a_time The time to be passed to the boundary condition.
  virtual void setTime(const Real& a_time)
  {
  }

  /// Fill the ghost cells for a single level
  void fillGhostCells(const LevelData<FArrayBox>& phi, const Real dx,
                      const bool homogeneous)
  {
    const DisjointBoxLayout& dbl = phi.disjointBoxLayout();
    LevelData<FArrayBox>& ph = (LevelData<FArrayBox>&) phi;
    for (DataIterator dit = ph.dataIterator(); dit.ok(); ++dit)
      operator()(ph[dit], dbl[dit], dbl.physDomain(), dx, homogeneous);
  }

  /// Fill the ghost cells for a Hierarchy of levels.
  void fillGhostCells(const Vector<LevelData<FArrayBox>*>& phi,
                      const Real dx0, const Vector<int>& refV,
                      const bool homogeneous)
  {
    Real dx = dx0;
    for (int i=0; i<phi.size(); i++)
      {
        fillGhostCells(*phi[i], dx, homogeneous);
        if (i<phi.size()-1)
          dx /= refV[i];
      }
  }

private:

  // Copy constructor and assignment operator are disallowed.
  BCFunction(const BCFunction&);
  BCFunction& operator=(const BCFunction&);

};

//! \class BCHolder
//! This class is a catch-all that can use a function pointer or a BCFunction
//! object to impose boundary conditions.
class BCHolder
{
public:

  //! Default constructor. Seems like this should be disallowed, since BCHolder
  //! isn't sub-classable, and there's no way to set the member pointer after
  //! creation.
  BCHolder():m_funcptr(NULL)
  {
  }

  //! Creates a BCHolder using a function pointer.
  BCHolder(BCFunc funcptr):m_funcptr(funcptr)
  {
  }

  //! Creates a BCHolder using a BCFunction instance.
  BCHolder(RefCountedPtr<BCFunction> refptr):m_funcptr(NULL),m_bc(refptr)
  {
  }

  //! Imposes a boundary condition on a solution whose state is described
  //! by the given FArrayBox. The signature for this method matches its
  //! counterpart in BCFunction.
  void operator()(FArrayBox&           a_state,
                  const Box&           a_valid,
                  const ProblemDomain& a_domain,
                  Real                 a_dx,
                  bool                 a_homogeneous,
                  const DataIndex      a_index = DataIndex())
  {
    if (m_funcptr != NULL)
    {
      m_funcptr(a_state, a_valid, a_domain, a_dx, a_homogeneous);
    }
    else
    {
      m_bc->operator()(a_state, a_valid, a_domain, a_dx, a_index, a_homogeneous);
    }
  }

  //! Sets the time associated with the boundary condition, a la BCFunction::setTime.
  //! If this BCHolder holds a function pointer, this call does nothing.
  void setTime(Real a_time)
  {
    if (!m_bc.isNull())
      m_bc->setTime(a_time);
  }

  /// Provides access to the bc function
  RefCountedPtr<BCFunction> getBCFunction()
  {
      return RefCountedPtr<BCFunction>(m_bc);
  }

protected:
  BCFunc m_funcptr;
  RefCountedPtr<BCFunction> m_bc;
};

///
/**
   given
   pos [x,y,z] position on center of cell edge
   int dir direction, x being 0
   int side -1 for low, +1 = high,
   fill in the a_values array
*/

typedef void(*BCValueFunc)(Real*           a_pos,
                           int*            a_dir,
                           Side::LoHiSide* a_side,
                           Real*           a_value);

class BCValueFunction
{
public:
  virtual ~BCValueFunction()
  {
  }

  virtual void operator()(Real*           a_pos,
                          int*            a_dir,
                          Side::LoHiSide* a_side,
                          Real*           a_value) = 0;
};

class BCValueHolder
{
public:
  BCValueHolder():m_funcptr(NULL)
  {
  }

  BCValueHolder(BCValueFunc funcptr):m_funcptr(funcptr)
  {
  }

  BCValueHolder(RefCountedPtr<BCValueFunction> refptr):m_funcptr(NULL),m_bc(refptr)
  {
  }

  virtual ~BCValueHolder()
  {
  }

  virtual void operator()(Real*           a_pos,
                          int*            a_dir,
                          Side::LoHiSide* a_side,
                          Real*           a_value)
  {
    if (m_funcptr != NULL)
    {
      m_funcptr(a_pos, a_dir, a_side, a_value);
    }
    else
    {
      m_bc->operator()(a_pos, a_dir, a_side, a_value);
    }
  }

protected:
  BCValueFunc m_funcptr;
  RefCountedPtr<BCValueFunction> m_bc;
};

// Class used by ConstDiriNeumBC to define BCs which can be passed anywhere a
// BCFunction/BCHolder is needed
class ConstBCFunction: public BCFunction
{
public:
  ConstBCFunction(const IntVect&  a_loSideType,
                  const RealVect& a_loSideValue,
                  const IntVect&  a_hiSideType,
                  const RealVect& a_hiSideValue);

  ~ConstBCFunction();

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous);

protected:
  IntVect  m_loSideType;
  RealVect m_loSideValue;

  IntVect  m_hiSideType;
  RealVect m_hiSideValue;
};

///
/**
   A helper function to produce the needed object for constant
   Dirichlet/Neumann on all the faces.  The return value can be passed to
   anything expecting a BCFunction/BCHolder.

   "a_loSideType/a_hiSideType" specify the type of boundary condition in a
   given direction by having the enter corresponding to the direction set to
   0 for Neumann or 1 for Dirichlet (on the low or high side, respectively).
   "a_loSideValue/a_hiSideValue" specify the (constant) value for boundary
   condition specified above.

   For example, in 2D if "a_loSideType" = (1,0), "a_hiSideType" = (1,1),
   "a_loSideValue" = (0.0,1.0) and "a_hiSideValue" = (0.0,0.0) then the
   boundary conditions are:

          Low  side x:  Dirichlet = 0.0
          Low  side y:  Neumann   = 1.0
          High side x:  Dirichlet = 0.0
          High side y:  Dirichlet = 0.0
 */
RefCountedPtr<BCFunction> ConstDiriNeumBC(const IntVect&  a_loSideType,
                                          const RealVect& a_loSideValue,
                                          const IntVect&  a_hiSideType,
                                          const RealVect& a_hiSideValue);

///
/**
   Neumann bc for a particular side, specified component interval
   For use in AMRPoissonOp.
 */
void NeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            const BCValueHolder&   a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            Interval&       a_interval);

///
/**
   Neumann bc for a particular side, all components
   For use in AMRPoissonOp.
 */
void NeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            const BCValueHolder&   a_value,
            int             a_dir,
            Side::LoHiSide  a_side);

///
/**
   Neumann bcs for all sides
   For use in AMRPoissonOp.
 */
void NeumBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value);

///
/**
   Dirichlet boundary conditions for a side, specified component interval
   For use in AMRPoissonOp.
 */
void DiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            Interval&       a_interval,
            int             a_order = 1);

///
/**
   Dirichlet boundary conditions for a side, all components.
   For use in AMRPoissonOp.
 */
void DiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_dir,
            Side::LoHiSide  a_side,
            int             a_order = 1);

///
/**
   Dirichlet boundary conditions for one side.
   For use in AMRPoissonOp.
 */
void DiriBC(FArrayBox&      a_state,
            const Box&      a_valid,
            Real            a_dx,
            bool            a_homogeneous,
            BCValueHolder   a_value,
            int             a_order = 1);

///
/**
   No slip vector bc (zero all comps).
   need a_state.ncomp == spacedim
   For use in ResistivityOp, for example.
 */
void NoSlipVectorBC(FArrayBox&     a_state,
                    const Box&     a_valid,
                    Real           a_dx,
                    int            a_dir,
                    Side::LoHiSide a_side,
                    int            a_order = 2);

///
/**
   0 normal comp, reflective for all other comps
   need a_state.ncomp == spacedim
   For use in ResistivityOp, for example.
 */
void ReflectiveVectorBC(FArrayBox&     a_state,
                        const Box&     a_valid,
                        Real           a_dx,
                        int            a_dir,
                        Side::LoHiSide a_side,
                        int            a_order = 2);

///
/**
   Extrapolation boundary conditions for a side, specified component interval
   For use in AMRPoissonOp.
 */
void ExtrapolateBC(FArrayBox&      a_state,
                   const Box&      a_valid,
                   Real            a_dx,
                   int             a_dir,
                   Side::LoHiSide  a_side,
                   Interval&       a_interval,
                   int             a_order = 1);

///
/**
   Extrapolation boundary conditions for a side, all components.
   For use in AMRPoissonOp.
 */
void ExtrapolateBC(FArrayBox&      a_state,
                   const Box&      a_valid,
                   Real            a_dx,
                   int             a_dir,
                   Side::LoHiSide  a_side,
                   int             a_order = 1);

///
/**
   Extrapolation boundary conditions for one side.
   For use in AMRPoissonOp.
 */
void ExtrapolateBC(FArrayBox&      a_state,
                   const Box&      a_valid,
                   Real            a_dx,
                   int             a_order = 1);

#include "NamespaceFooter.H"

#endif
