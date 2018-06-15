#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SumIF.H"

#include "NamespaceHeader.H"

SumIF::SumIF(const BaseIF& a_impFunc1,
             const BaseIF& a_impFunc2,
             const bool&   a_sign)
{
  // Number of implicit function in sum
  m_numFuncs = 2;

  // Vector of implicit function pointers
  m_impFuncs.resize(m_numFuncs);

  // Make copies of the implicit functions
  m_impFuncs[0] = a_impFunc1.newImplicitFunction();
  m_impFuncs[1] = a_impFunc2.newImplicitFunction();

  //sum or difference
  m_sign = a_sign;
}

SumIF::SumIF(const Vector<BaseIF *>& a_impFuncs)
{
  // Number of implicit function in sum
  m_numFuncs = a_impFuncs.size();

  // Vector of implicit function pointers
  m_impFuncs.resize(m_numFuncs);

  // Make copies of the implicit functions

  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
  {
    if (a_impFuncs[ifunc] == NULL)
    {
      m_impFuncs[ifunc] = NULL;
    }
    else
    {
      m_impFuncs[ifunc] = a_impFuncs[ifunc]->newImplicitFunction();
    }
  }
  //sign to make the sum
  m_sign = true;
}

SumIF::SumIF(const SumIF& a_inputIF)
{
  // Number of implicit function in sum
  m_numFuncs = a_inputIF.m_impFuncs.size();

  // Vector of implicit function pointers
  m_impFuncs.resize(m_numFuncs);

  //get the sign
  m_sign =  a_inputIF.getSign();

  // Make copies of the implicit functions
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
  {
    if (a_inputIF.m_impFuncs[ifunc] == NULL)
    {
      m_impFuncs[ifunc] = NULL;
    }
    else
    {
      m_impFuncs[ifunc] = a_inputIF.m_impFuncs[ifunc]->newImplicitFunction();
    }
  }
}

SumIF::~SumIF()
{
  // Delete all the copies
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
  {
    if (m_impFuncs[ifunc] != NULL)
    {
      delete m_impFuncs[ifunc];
    }
  }
}

void SumIF::setSign(bool a_sign)
{
  m_sign = a_sign;
}

Real SumIF::value(const RealVect& a_point) const

{

  // sum of the implicit functions values
  Real retval;
  retval = 0.0;

  if (m_sign)
    {
      // sum values and return it
      if (m_numFuncs > 0)
        {
          retval = m_impFuncs[0]->value(a_point);

          for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
            {
              retval=retval + m_impFuncs[ifunc]->value(a_point);
            }
        }
    }
  else
    {
//do the difference
      if (m_numFuncs != 2)
        {
          MayDay::Abort("Cannot make the difference for m_numFuncs != 2, m_sign should be equal to true");
        }
      retval =  m_impFuncs[0]->value(a_point) - m_impFuncs[1]->value(a_point);
    }

  return retval;
}

Real SumIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{

  Real retval;
  retval = 0.0;

  if (m_sign)
    {
      if (m_numFuncs > 0)
        {
          retval = m_impFuncs[0]->value(a_point);

          for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
            {
              retval=retval + m_impFuncs[ifunc]->value(a_point);
            }
        }
    }

  else
    {
      //do the difference
      if (m_numFuncs != 2)
        {
          MayDay::Abort("Cannot make the difference for m_numFuncs != 2, m_sign should be equal to true");
        }
      retval =  m_impFuncs[0]->value(a_point) - m_impFuncs[1]->value(a_point);
    }

  return retval;
}

Real SumIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                           const IndexTM<Real,GLOBALDIM>& a_point) const
{

  Real retval;
  retval = 0.0;

  if (m_sign)
    {
      if (m_numFuncs > 0)
        {
          retval = m_impFuncs[0]->value(a_partialDerivative,a_point);

          for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
            {
              retval=retval + m_impFuncs[ifunc]->value(a_partialDerivative,a_point);
            }
        }
    }
  else
    {
      //do the difference
      if (m_numFuncs != 2)
        {
          MayDay::Abort("Cannot make the difference for m_numFuncs != 2, m_sign should be equal to true");
        }
      retval =  m_impFuncs[0]->value(a_partialDerivative,a_point) - m_impFuncs[1]->value(a_partialDerivative,a_point);
    }

  return retval;
}

bool SumIF::getSign() const
{
  return m_sign;
}

BaseIF* SumIF::getImplicitFunction(int a_num)
{
  return m_impFuncs[a_num];
}


BaseIF* SumIF::newImplicitFunction() const
{
  SumIF* sumPtr = new SumIF(m_impFuncs);
  sumPtr->setSign(m_sign);

  return static_cast<BaseIF*>(sumPtr);
}



#include "NamespaceFooter.H"
