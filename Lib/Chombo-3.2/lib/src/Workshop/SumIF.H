#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SUMIF_H_
#define _SUMIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function is the sum of two (of more) implicit
    functions with which it is constructed.
 */
class SumIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying two implicit functions as inputs.
   */
  SumIF(const BaseIF& a_impFunc1,
        const BaseIF& a_impFunc2,
        const bool&   a_sign);

  ///
  /**
      Constructor specifying any number of implicit functions as inputs.
   */
  SumIF(const Vector<BaseIF*>& a_impFuncs);

  /// Copy constructor
  SumIF(const SumIF& a_inputIF);

  /// Destructor
  virtual ~SumIF();

 void setSign(bool a_sign);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual bool getSign() const;

  virtual BaseIF* getImplicitFunction(int a_num);

  virtual BaseIF* newImplicitFunction() const;

  ///
  /**
     Pass this call onto the IFs contained in this IF class.
  */
  virtual void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                                const RealVect          & a_dx)
  {
    for (int i = 0; i < m_numFuncs; i++)
    {
      m_impFuncs[i]->boxLayoutChanged(a_newBoxLayout,a_dx);
    }
  }

protected:
  int             m_numFuncs; // number of implicit functions
  Vector<BaseIF*> m_impFuncs; // implicit functions
  bool            m_sign;    //  if true -> sum of all IF //if false the first IF minus the second one (never happen when m_numFuncs > 2)

private:
  SumIF()
  {
    MayDay::Abort("SumIF uses strong construction");
  }

  void operator=(const SumIF& a_inputIF)
  {
    MayDay::Abort("SumIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
