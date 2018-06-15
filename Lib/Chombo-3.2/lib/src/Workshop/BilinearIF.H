#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BILINEARIF_H_
#define _BILINEARIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "IndexTM.H"
#include "FArrayBox.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
   This implicit function specifies a general bilinear function.
*/
class BilinearIF: public BaseIF
{
public:
  ///
  /**
     Constructor specifying the value to interpolate with a bilinear function (a_etaCorner)
     and the properties of the grid (a_origin and a_dx)
  */

BilinearIF(LevelData<NodeFArrayBox>     * a_etaCorner,
           const IndexTM<Real,GLOBALDIM>& a_origin3D,
           const IndexTM<Real,GLOBALDIM>& a_dx3D,
           DisjointBoxLayout            * a_grids);

  /// Destructor
virtual ~BilinearIF();

//copy constructor

BilinearIF(const BilinearIF& a_inputIF);

///
/**
   Put the grid index of the lower left node of the box which contain a_point in a_index
*/

virtual void findIndex(const IndexTM<Real,GLOBALDIM>& a_point,
                       IntVect                      & a_index) const;

  ///
  /**
     Return the value of the derivative at a_point.
  */
  virtual Real value(const IndexTM<int,GLOBALDIM>  & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM> & a_point            ) const;

  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual LevelData<NodeFArrayBox>* getEtaCorner() const;

  virtual DisjointBoxLayout* getGrid() const;

  virtual void getPolynomial(Vector<PolyTerm> & a_polynomial,
                             IntVect          & a_index);
private:

 BilinearIF()
  {
    MayDay::Abort("BilinearIF uses strong construction");
  }

  void operator=(const BilinearIF& a_inputIF)
  {
    MayDay::Abort("BilinearIF doesn't allow assignment");
  }

  //Node values of eta used for the bilinear interpolation
  LevelData<NodeFArrayBox>* m_etaCorner;

  //grids info to link indexes and physical coordinates
  IndexTM<Real,GLOBALDIM> m_origin3D;
  IndexTM<Real,GLOBALDIM> m_dx3D;
  DisjointBoxLayout* m_grids;



};

#include "NamespaceFooter.H"
#endif
