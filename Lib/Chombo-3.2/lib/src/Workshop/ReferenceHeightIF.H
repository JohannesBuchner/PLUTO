#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REFERENCEHEIGHTIF_H_
#define _REFERENCEHEIGHTIF_H_

#include "Notation.H"
#include "RealVect.H"
#include "IndexTM.H"
#include "BaseIF.H"
#include "PolynomialIF.H"

#include "NamespaceHeader.H"

///
/**
   This adapter wraps any BaseIF subclass and gives it the
   the concept of a reference height. The adapter passes through
   all function calls to the adapted object, except value(RealVect)..
   The reference height is used to evaluate the value function when
   SpaceDim=GLOBALDIM-1,

 */
class ReferenceHeightIF : public BaseIF
{
public:

  /**
   Construct the ReferenceHeightIF using another IF plus
   a reference height. The delegate IF will be copied.
  */
  ReferenceHeightIF(const  BaseIF                 & a_implicitFunction,
                    const  Real                   & a_referenceHeight,
                    const IndexTM<Real,GLOBALDIM> & a_origin)
    :
    m_implicitFunction(a_implicitFunction.newImplicitFunction()),
    m_referenceSurface(NULL),
    m_referenceHeight(a_referenceHeight),
    m_origin(a_origin)
  {
  }

  ReferenceHeightIF(const  BaseIF                 & a_implicitFunction,
                    const  PolynomialIF           & a_referenceSurface,
                    const  Real                   & a_referenceHeight,
                    const IndexTM<Real,GLOBALDIM> & a_origin)
    :
    m_implicitFunction(a_implicitFunction.newImplicitFunction()),
    m_referenceSurface(static_cast<PolynomialIF*>(a_referenceSurface.newImplicitFunction())),
    m_referenceHeight(a_referenceHeight),
    m_origin(a_origin)
  {
  }

  /// Destructor cleans up the member implicit function
  virtual ~ReferenceHeightIF()
  {
    delete m_implicitFunction;
  }

  Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
             const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    Real ifVal = m_implicitFunction->value(a_partialDerivative,a_point);

    //There is a constructor that doesn't create a reference surface.
   if (m_referenceSurface != NULL)
      {
        Real refSurfVal = m_referenceSurface->value(a_partialDerivative, a_point);
        ifVal          -= refSurfVal;
      }

    return (ifVal);
  }
  ///
  /**
     Return the value of the function at a_point. value(RealVect)
     to value(IndexTM). However, in the case where GLOBALDIM == 3 and
     SpaceDim == 2 the missing dimension is filled in with the reference height.
     Subclasses need not implement this function, but if

  */
  virtual Real value(const RealVect& a_point) const
  {
    IndexTM<Real,GLOBALDIM> pt;

    for (int idir = 0; idir < SpaceDim; ++idir)
      {
        pt[idir] = a_point[idir];
      }

    if (GLOBALDIM == 3 && SpaceDim == 2)
    {
      pt[SpaceDim] = m_referenceHeight + m_origin[SpaceDim];
    }

    Real ifVal = value(pt);

    if (m_referenceSurface != NULL)
      {
        Real refSurfVal = m_referenceSurface->value(a_point);
        ifVal          -= refSurfVal;
      }

    else
      {
        MayDay::Abort("No reference surface defined");
      }

    return (ifVal);
  }

   ///
  /**
     Return the value of the function at a_point (of type INdexTM).
  */
  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const
  {
    return m_implicitFunction->value(a_point);
  };

  ///
  /**
     Return a newly allocated derived class.  The responsibility
     for deleting the memory is left to the calling function.
  */
  virtual ReferenceHeightIF* newImplicitFunction() const
  {
    if (m_referenceSurface != NULL)
      {
        return new ReferenceHeightIF(*m_implicitFunction, *m_referenceSurface, m_referenceHeight, m_origin);
      }
    else
      {
        return new ReferenceHeightIF(*m_implicitFunction, m_referenceHeight,m_origin);
      }
  }
  virtual void print(ostream& a_out) const
  {
    m_implicitFunction->print(a_out);
  };

  ///
  /**
    Return the reference height for this IF
  */
  Real getReferenceHeight() const
  {
    return m_referenceHeight;
  }

  ///
  /**
    Return the physical origin of the domain
  */
  IndexTM<Real,GLOBALDIM> getOrigin() const
  {
    return m_origin;
  }

  ///
  /**
    Return the reference surface for this IF
  */
  void getReferenceSurfacePolynomial(Vector<PolyTerm>& a_polynomial,
                                     bool            & a_inside) const
  {
    return m_referenceSurface->GetParams(a_polynomial,a_inside);
  }

  ///
  /**
     Evaluate the reference surface at a point
  */
  Real evaluateReferenceSurfacePolynomial(const RealVect & a_point) const
  {
    return m_referenceSurface->value(a_point);
  }

  ///
  /**
     Pass this call onto the IFs contained in this IF class.
  */
  virtual void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                                const RealVect          & a_dx)
  {
    m_implicitFunction->boxLayoutChanged(a_newBoxLayout,a_dx);

    if (m_referenceSurface != NULL)
    {
      m_referenceSurface->boxLayoutChanged(a_newBoxLayout,a_dx);
    }
  }

///
/**
return the pointer to the implicit function m_implicitFunction
*/

  BaseIF* getPointer2IF()const
  {
    return m_implicitFunction;
  }

//return a pointer to a new ReferenceHeight that has a_newIF as m_implicitFunction and NULL as m_referenceSurface
  ReferenceHeightIF* newChangedIF(const BaseIF* a_newIF )const
  {
    return new ReferenceHeightIF( *a_newIF, m_referenceHeight, m_origin);
  }

  bool hasReferenceSurface() const
  {
    if (m_referenceSurface==NULL)
      {
        return false;
      }
    return true;
  }

private:

  // Default constructor is disabled in order to prevent initialization
  // without reference height. Base classes should all use the constructor
  // version that takes a reference height
  ReferenceHeightIF()
  {
  }

  BaseIF*                 m_implicitFunction;
  PolynomialIF*           m_referenceSurface;
  Real                    m_referenceHeight;
  IndexTM<Real,GLOBALDIM> m_origin;
};

#include "NamespaceFooter.H"
#endif
