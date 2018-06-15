#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IntVect.H"

#include "PolynomialIF.H"

#include "BilinearIF.H"
#include "NamespaceHeader.H"



BilinearIF::BilinearIF(LevelData<NodeFArrayBox>       * a_etaCorner,
                       const IndexTM<Real,GLOBALDIM>  & a_origin,
                       const IndexTM<Real,GLOBALDIM>  & a_dx,
                       DisjointBoxLayout              * a_grids)
  :  m_etaCorner(a_etaCorner),
     m_origin3D(a_origin),
     m_dx3D(a_dx),
     m_grids(a_grids)
{
}

BilinearIF::BilinearIF(const BilinearIF& a_inputIF)
{
  m_etaCorner = a_inputIF.getEtaCorner();
  m_origin3D = a_inputIF.m_origin3D;
  m_dx3D = a_inputIF.m_dx3D;
  m_grids = a_inputIF.getGrid();
}

BilinearIF::~BilinearIF()
{
}

LevelData<NodeFArrayBox>* BilinearIF::getEtaCorner() const
{
  return m_etaCorner;
}

DisjointBoxLayout* BilinearIF::getGrid() const
{
  return m_grids;
}

void BilinearIF::findIndex(const IndexTM<Real,GLOBALDIM>& a_point ,IntVect& a_index) const
{
  const int indicex = floor((a_point[0] - m_origin3D[0])/m_dx3D[0]);
  const int indicey = floor((a_point[1] - m_origin3D[1])/m_dx3D[1]);
  a_index.setVal(0,indicex);
  a_index.setVal(1,indicey);
}


Real BilinearIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivativeOp,
                       const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Real retval = LARGEREALVAL;
  //eta is the first component in FArrayBox
  int etacomp = 0;

  //find the coordinates of lower left corner of the box containing a_point
  IntVect index = IntVect::Zero;
  findIndex(a_point,index);

  //find values to interpolate:hil,hir,hsl,hsr
  Real hCorner[4];
  //iterate through LevelData of m_etaCorner
  for (DataIterator dit = m_etaCorner->dataIterator();dit.ok();++dit)
    {
      NodeFArrayBox& nodeFAB = (*m_etaCorner)[dit()];

      Box boxIterNode = (*m_grids)[dit()];
      boxIterNode.surroundingNodes();

      if (boxIterNode.contains(index))
        {
          hCorner[0] = nodeFAB(index,                       etacomp);
          hCorner[1] = nodeFAB(index + BASISV(0),           etacomp);
          hCorner[3] = nodeFAB(index + BASISV(1)+ BASISV(0),etacomp);
          hCorner[2] = nodeFAB(index + BASISV(1),           etacomp);
          break;
        }
    }

  //physical coordinates of lower left corner
  RealVect corner = RealVect::Zero;
  for (int idir = 0; idir < SpaceDim ; ++idir)
    {
      corner[idir] = index[idir] * m_dx3D[idir] + m_origin3D[idir];
    }

  //coeficients of bilinear interpolant
  const Real b = (hCorner[1]-hCorner[0])/m_dx3D[0];
  const Real c = (hCorner[2]-hCorner[0])/m_dx3D[1];
  const Real a = (hCorner[3]+hCorner[0]-hCorner[1]-hCorner[2])/(m_dx3D[0]*m_dx3D[1]);


  if ((a_partialDerivativeOp[0]==0)&&(a_partialDerivativeOp[1]==0))
    {
      retval = a*(a_point[0]-corner[0])*(a_point[1]-corner[1])+b*(a_point[0]-corner[0])+c*(a_point[1]-corner[1])+hCorner[0];
    }
  else if ((a_partialDerivativeOp[0]==1)&&(a_partialDerivativeOp[1]==0))
    {
      retval =  a*(a_point[1]-corner[1])+b;
    }
  else if ((a_partialDerivativeOp[0]==0)&&(a_partialDerivativeOp[1]==1))
    {
      retval = a*(a_point[0]-corner[0])+c;
    }
  else if ((a_partialDerivativeOp[0]==1)&&(a_partialDerivativeOp[1]==1))
    {
      retval = a;
    }
  else
    {
      retval = 0.0;
    }
  return retval;
}

Real BilinearIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  // Because we are intrested in the (0,0) derivative
  int ival[3];ival[0]=0;ival[1]=0;ival[2]=0;
  const IndexTM<int,GLOBALDIM> partialDerivativeOp(ival);

  return value(partialDerivativeOp,a_point);

}

Real BilinearIF::value(const RealVect& a_point) const
{
  // Because we are intrested in the (0,0) derivative
   int ival[3];
   ival[0] = 0;
   ival[1] = 0;
   ival[2] = 0;
   const IndexTM<int,GLOBALDIM> partialDerivativeOp(ival);

 // Because we are intrested in the a_point
   Real jval[3];
   jval[0] = a_point[0];
   jval[1] = a_point[1];
   jval[2] = m_dx3D[2];
  const IndexTM<Real,GLOBALDIM> point(jval);

  return value(partialDerivativeOp,point);
}


void BilinearIF::getPolynomial(Vector<PolyTerm> & a_polynomial,
                               IntVect          & a_index)
{
  //eta is the first component in FArrayBox
  int etacomp = 0;
  //find values to interpolate:hil,hir,hsl,hsr
  Real hCorner[4];

  for (DataIterator dit = m_etaCorner->dataIterator();dit.ok();++dit)
    {
      NodeFArrayBox& nodeFAB = (*m_etaCorner)[dit()];

      Box boxIterNode = (*m_grids)[dit()];
      boxIterNode.surroundingNodes();

      if (boxIterNode.contains(a_index))
        {
          hCorner[0] = nodeFAB(a_index,                       etacomp);
          hCorner[1] = nodeFAB(a_index + BASISV(0),           etacomp);
          hCorner[3] = nodeFAB(a_index + BASISV(1)+ BASISV(0),etacomp);
          hCorner[2] = nodeFAB(a_index + BASISV(1),           etacomp);
          break;
        }
    }

  //coeficients of bilinear interpolant
  const Real b = (hCorner[1]-hCorner[0])/m_dx3D[0];
  const Real c = (hCorner[2]-hCorner[0])/m_dx3D[1];
  const Real a = (hCorner[3]+hCorner[0]-hCorner[1]-hCorner[2])/(m_dx3D[0]*m_dx3D[1]);

  //physical coordinates of lower left corner
  RealVect corner = RealVect::Zero;
  for (int idir = 0; idir < SpaceDim ; ++idir)
    {
      corner[idir] = a_index[idir] * m_dx3D[idir] + m_origin3D[idir];
    }

  Real polyCoeff[4];
  polyCoeff[3] = a;
  polyCoeff[2] = b - a*corner[1];
  polyCoeff[1] = c - a*corner[0];
  polyCoeff[0] = a*corner[0]*corner[1] - b*corner[0] - c*corner[1] + hCorner[0];
  IntVect polyPower[4];
  polyPower[3][0] = 1;
  polyPower[3][1] = 1;
  polyPower[2][0] = 1;
  polyPower[2][1] = 0;
  polyPower[1][0] = 0;
  polyPower[1][1] = 1;
  polyPower[0][0] = 0;
  polyPower[0][1] = 0;

  a_polynomial.resize(4);

  a_polynomial[3].coef = polyCoeff[3];
  a_polynomial[3].powers = polyPower[3];

  a_polynomial[2].coef = polyCoeff[2];
  a_polynomial[2].powers = polyPower[2];

  a_polynomial[1].coef = polyCoeff[1];
  a_polynomial[1].powers = polyPower[1];

  a_polynomial[0].coef = polyCoeff[0];
  a_polynomial[0].powers = polyPower[0];
}


BaseIF* BilinearIF::newImplicitFunction() const
{
  BilinearIF* bilinearPtr = new BilinearIF( m_etaCorner,
                                            m_origin3D,
                                            m_dx3D,
                                            m_grids);

  return static_cast<BaseIF*>(bilinearPtr);
}

#include "NamespaceFooter.H"
