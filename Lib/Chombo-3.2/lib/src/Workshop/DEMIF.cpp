#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
using std::ifstream;

#include "DEMIF.H"

#include "NamespaceHeader.H"

DEMIF::DEMIF(const IntVect&     a_ncell,
             const int&         a_interpType,
             const RealVect&    a_dx,
             const std::string& a_demFile,
             const Real&        a_bottomBuffer,
             const Real&        a_truncElev,
             const Real&        a_highGround,
             const Real&        a_verticalScale)
{
  m_ncell         = a_ncell;
  m_dx            = a_dx;
  m_interpType    = a_interpType;
  CH_assert(m_interpType!=2);//quadratic not supported anymore
  if (m_interpType==1)
    {
      m_doCubic = false;
    }
  else
    {
      m_doCubic = true;
    }

  m_bottomBuffer  = a_bottomBuffer;
  m_truncElev     = a_truncElev;
  m_highGround    = a_highGround;
  m_verticalScale = a_verticalScale;

  // Read in the header info
  bool justhead = true;
  readDEM(justhead,
          a_demFile);

  //allocate memory for DEM matrix
  m_DEM = new Real *[m_ncols];
  for (int i = 0; i < m_ncols;i++)
    {
      m_DEM[i] = new Real[m_nrows];
    }

  // Read in the whole file
  justhead = false;
  readDEM(justhead,
          a_demFile);

  //cache constants based on m_cellsize
  cacheConstants();

  //fill in stray "wet" cells, change nodata value to land
  // and convert DEM to meters.
  fixDEM();
}

bool DEMIF::cacheConstants()
{
  m_hx     = m_cellsize;
  m_hy     = m_cellsize;
  m_hx2    = pow(m_cellsize,2);
  m_hx3    = pow(m_cellsize,3);
  m_hy2    = pow(m_cellsize,2);
  m_hy3    = pow(m_cellsize,3);
  bool done = true;
  return done;
}

DEMIF::DEMIF(const DEMIF& a_inputIF)
{
  m_ncols         = a_inputIF.m_ncols;
  m_nrows         = a_inputIF.m_nrows;
  m_NODATA        = a_inputIF.m_NODATA;
  m_xllcorner     = a_inputIF.m_xllcorner;
  m_yllcorner     = a_inputIF.m_yllcorner;
  m_cellsize      = a_inputIF.m_cellsize;
  m_cellvalue     = a_inputIF.m_cellvalue;
  m_dx            = a_inputIF.m_dx;
  m_ncell         = a_inputIF.m_ncell;
  m_interpType    = a_inputIF.m_interpType;
  m_highGround    = a_inputIF.m_highGround;
  m_minDEM        = a_inputIF.m_minDEM;
  m_maxDEM        = a_inputIF.m_maxDEM;
  m_bottomBuffer  = a_inputIF.m_bottomBuffer;
  m_truncElev     = a_inputIF.m_truncElev;
  m_verticalScale = a_inputIF.m_verticalScale;

  if (m_interpType==1)
    {
      m_doCubic = false;
    }
  else
    {
      m_doCubic = true;
    }

  //allocate and copy memory for DEM matrix
  m_DEM = new Real *[m_ncols];

  for (int i = 0; i < m_ncols;i++)
    {
      m_DEM[i] = new Real[m_nrows];
      for (int j = 0; j < m_nrows;j++)
        {
          m_DEM[i][j] = a_inputIF.m_DEM[i][j];
        }
    }

  //cache constants based on m_cellsize
  cacheConstants();
}

DEMIF::~DEMIF()
{
  for (int i = 0; i < m_ncols;i++)
    {
      delete [] m_DEM[i];
    }
  delete [] m_DEM;
}

//Reads in a Digital Elevation Model file
bool DEMIF::readDEM(const bool&        a_justhead,
                    const std::string& a_demFile)
{
  bool fileread;
  //pout()<<"Reading DEM Header...\n";
  ifstream DEM_file;
  DEM_file.open(a_demFile.c_str()); //the input file
  if (DEM_file.bad())
    {
      MayDay::Abort("Bad DEM_file in readDEM");
    }
  char astring[1024];
  m_ncols = -99;

  DEM_file>>astring;
  DEM_file>>m_ncols;

  DEM_file>>astring;
  DEM_file>>m_nrows;

  DEM_file>>astring;
  DEM_file>>m_xllcorner;

  DEM_file>>astring;
  DEM_file>>m_yllcorner;

  DEM_file>>astring;
  DEM_file>>m_cellsize;

  DEM_file>>astring;
  DEM_file>>m_cellvalue;

  DEM_file>>astring;
  DEM_file>>m_NODATA;

  if (SpaceDim==2)
    {
      if (m_nrows!=1)
        {//2D column data...switch m_nrows and m_ncols
          int temp = m_nrows;
          m_nrows = m_ncols;
          m_ncols = temp;
        }
      else if (m_ncols==1)
        {
          MayDay::Abort("DEMIF::ncols and nrows are 1...need a bigger dataset");
        }
      CH_assert(m_nrows==1);
    }
  else if (SpaceDim==3 && m_nrows==1)
    {
      MayDay::Abort("DEMIF::nrows must not equal 1 for 3D");
    }

  //make sure the grid is big enough
  if (m_ncols==-99)
    {
      MayDay::Abort("DEMIF::Something wrong with DEM file");
    }

  if (!a_justhead)
    {
      m_minDEM =  1e20;
      m_maxDEM = -1e20;
      for (int j=m_nrows-1;j>=0;j--)
        {
          for (int i=0;i<m_ncols;i++)
            {
              //load in the data
              DEM_file>>m_DEM[i][j];

              //vertical shift and scale
              if (rint(m_DEM[i][j]) != m_NODATA)
                {
                  Real val = (m_bottomBuffer + m_DEM[i][j]*m_cellvalue)*m_verticalScale;
                  m_DEM[i][j] = val;
                  if (val < m_minDEM)
                    {
                      m_minDEM = val; //this is the highest spot so far...
                    }
                  if (val > m_maxDEM)
                    {
                      m_maxDEM = val; //this is the deepest spot so far...
                    }
                }
            }
        }
      pout()<<"  DEM min = "<<m_minDEM<<endl;
      pout()<<"  DEM max = "<<m_maxDEM<<endl;
    }
  DEM_file.close();

//   pout()<<"Warning: Setting (xllcorner,yllcorner) of DEM to (0.,0.)"<<endl;
//   m_xllcorner=0.;
//   m_yllcorner=0.;

  fileread = true;
  return fileread;
}

// change nodata value to land
// and convert DEM to meters.
bool DEMIF::fixDEM()
{
  bool DEMfixed = true;

  pout()<<"Changing fixing NODATA cells in DEM\n";

  int kn = 0;
  int kt = 0;
  for (int i=0;i<m_ncols;i++)
    {
      for (int j=0;j<m_nrows;j++)
        {
          if (rint(m_DEM[i][j]) == m_NODATA)
            {
              m_DEM[i][j] = m_highGround;
              kn++;
            }
          if (m_DEM[i][j] > m_truncElev)
            {
              m_DEM[i][j] = m_highGround;
              kt++;
            }
        }
    }
  pout()<<"Filled: "<<kn<<" no-data cells."<<endl;
  pout()<<"Filled: "<<kt<<" shallow cells."<<endl;
  return DEMfixed;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Real DEMIF::value(const RealVect& a_point) const
{

  Real retval =0.;

  const RealVect& x = a_point;

  //where xllcorner is the x-lower-left corner of the cell-centered DEM
  //find nearest cell in DEM (recall that the DEM values are cell centered)
  Real xi,yj;
  xi = (x[0]-m_xllcorner)/m_cellsize - 0.5;
  yj = (x[1]-m_yllcorner)/m_cellsize - 0.5;
  Real xiFloor = floor(xi);
  Real yjFloor = floor(yj);

  int i,j;
  // centered on i+1/2
  i = (int)xiFloor;
  j = (int)yjFloor;

#if CH_SPACEDIM==2
  j = 0;
#endif

  int sz;//this is the size of the interp stencil
  if (m_doCubic)
    {//cubic
     sz = 1;
    }
  else
    {//linear
     sz = 0;
    }

  //the following checks to make sure we are inside of the DEM database (otherwise return a NODATA value)
  if (i-sz<0)
    {
      retval=m_highGround;
    }
  else if (i+sz>m_ncols-2)
    {
      retval=m_highGround;
    }
#if CH_SPACEDIM ==2
  else if (j != 0)
    {
      MayDay::Abort("j not equal to 0 for 2d in DEMIF");
    }
#elif CH_SPACEDIM == 3
  else if (j-sz<0)
    {
      retval=m_highGround;
    }
  else if (j+sz>m_nrows-2)
    {
      retval=m_highGround;
    }
#endif
  else
    {
      if (m_interpType==1)
        { //bilinear interpolation
          RealVect a   = RealVect::Zero;
          Real xLo = m_xllcorner + (i+0.5)*m_cellsize; //this is the x position of the low side of this dem cell (i,j)
          Real yLo = m_yllcorner + (j+0.5)*m_cellsize; //this is the y position of the low side of this dem cell (i,j)

          //a is the (bi)linear coefficient
          a[0] = (x[0]-xLo)/m_cellsize;
          a[1] = (x[1]-yLo)/m_cellsize;

          const Real& fLoLo = m_DEM[i  ][j  ];
          const Real& fHiLo = m_DEM[i+1][j  ];
          Real rLo = (1.0-a[0])*fLoLo + a[0]*fHiLo;
          retval = rLo;

#if CH_SPACEDIM==3
          const Real& fLoHi = m_DEM[i  ][j+1];
          const Real& fHiHi = m_DEM[i+1][j+1];
          Real rHi = (1.0-a[0])*fLoHi + a[0]*fHiHi;
          retval   = (1.0-a[1])*rLo   + a[1]*rHi;
#endif
        }
      else if (m_interpType==3)
        { //bicubic interpolation type1 (this performs cubic interpolation in 1D, matching all values)
          //now let's do some bicubic interpolation...
#if CH_SPACEDIM==3
          const Real& f1 = m_DEM[i-1][j-1];
          const Real& f2 = m_DEM[i  ][j-1];
          const Real& f3 = m_DEM[i+1][j-1];
          const Real& f4 = m_DEM[i+2][j-1];
#endif
          const Real& f5 = m_DEM[i-1][j  ];
          const Real& f6 = m_DEM[i  ][j  ];
          const Real& f7 = m_DEM[i+1][j  ];
          const Real& f8 = m_DEM[i+2][j  ];
#if CH_SPACEDIM==3
          const Real& f9 = m_DEM[i-1][j+1];
          const Real& f10= m_DEM[i  ][j+1];
          const Real& f11= m_DEM[i+1][j+1];
          const Real& f12= m_DEM[i+2][j+1];
          const Real& f13= m_DEM[i-1][j+2];
          const Real& f14= m_DEM[i  ][j+2];
          const Real& f15= m_DEM[i+1][j+2];
          const Real& f16= m_DEM[i+2][j+2];
#endif
          Real xd=x[0]-m_xllcorner-m_hx*(float(i+1)); //this is the x-distance from (i+1/2,j+1/2)
          Real x2 = pow(xd,2);
          Real x3 = pow(xd,3);
#if CH_SPACEDIM==3
          Real yd=x[1]-m_yllcorner-m_hy*(float(j+1)); //this is the y-distance from (i+1/2,j+1/2)
          Real y2 = pow(yd,2);
          Real y3 = pow(yd,3);
#endif

#if CH_SPACEDIM==3
          //cubic type1 in x direction: for row j-1
          Real r1 = (-3.*m_hx3*( f1 - 9.*(f2 + f3) + f4)    +
                     2.*m_hx2*( f1 -27.*(f2 - f3) - f4)*xd  +
                     12.*m_hx*( f1 -    (f2 + f3) + f4)*x2 +
                     8.*(-f1 + 3.*(f2 - f3) + f4)*x3)/(48.*m_hx3);
#endif
          //cubic type1 in x direction: for row j
          Real r2 = (-3.*m_hx3*( f5 - 9.*(f6 + f7) + f8)    +
                     2.*m_hx2*( f5 -27.*(f6 - f7) - f8)*xd  +
                     12.*m_hx*( f5 -    (f6 + f7) + f8)*x2 +
                     8.*(-f5 + 3.*(f6 - f7) + f8)*x3)/(48.*m_hx3);
          retval = r2;
#if CH_SPACEDIM==3
          //cubic type1 in x direction: for row j+1
          Real r3 = (-3.*m_hx3*( f9 - 9.*(f10 + f11) + f12)    +
                     2.*m_hx2*( f9 -27.*(f10 - f11) - f12)*xd  +
                     12.*m_hx*( f9 -    (f10 + f11) + f12)*x2 +
                     8.*(-f9 + 3.*(f10 - f11) + f12)*x3)/(48.*m_hx3);

          //cubic type1 in x direction: for row j+2
          Real r4 = (-3.*m_hx3*( f13 - 9.*(f14 + f15) + f16)    +
                     2.*m_hx2*( f13 -27.*(f14 - f15) - f16)*xd  +
                     12.*m_hx*( f13 -    (f14 + f15) + f16)*x2 +
                     8.*(-f13 + 3.*(f14 - f15) + f16)*x3)/(48.*m_hx3);

          //cubic type1 in y direction!!!
          retval  = (-3.*m_hy3*( r1 - 9.*(r2 + r3) + r4)    +
                     2.*m_hy2*( r1 -27.*(r2 - r3) - r4)*yd  +
                     12.*m_hy*( r1 -    (r2 + r3) + r4)*y2 +
                     8.*(-r1 + 3.*(r2 - r3) + r4)*y3)/(48.*m_hy3);
#endif
        }
      else if (m_interpType==4)
        { //bicubic interpolation type2 (this performs cubic interpolation in 1D, matching values and derivatives)

          //now let's do some bicubic interpolation...
#if CH_SPACEDIM==3
          const Real& f1 = m_DEM[i-1][j-1];
          const Real& f2 = m_DEM[i  ][j-1];
          const Real& f3 = m_DEM[i+1][j-1];
          const Real& f4 = m_DEM[i+2][j-1];
#endif
          const Real& f5 = m_DEM[i-1][j  ];
          const Real& f6 = m_DEM[i  ][j  ];
          const Real& f7 = m_DEM[i+1][j  ];
          const Real& f8 = m_DEM[i+2][j  ];
#if CH_SPACEDIM==3
          const Real& f9 = m_DEM[i-1][j+1];
          const Real& f10= m_DEM[i  ][j+1];
          const Real& f11= m_DEM[i+1][j+1];
          const Real& f12= m_DEM[i+2][j+1];
          const Real& f13= m_DEM[i-1][j+2];
          const Real& f14= m_DEM[i  ][j+2];
          const Real& f15= m_DEM[i+1][j+2];
          const Real& f16= m_DEM[i+2][j+2];
#endif
          Real xd=x[0]-m_xllcorner-m_hx*(float(i+1)); //this is the x-distance from (i+1/2,j+1/2)
          Real x2 = pow(xd,2);
          Real x3 = pow(xd,3);
#if CH_SPACEDIM==3
          Real yd=x[1]-m_yllcorner-m_hy*(float(j+1)); //this is the y-distance from (i+1/2,j+1/2)
          Real y2 = pow(yd,2);
          Real y3 = pow(yd,3);

          //cubic type2 in x direction: for row j-1
          Real r1 = (-f1 + 9.*(f2 + f3) - f4)   /    16.  +
            ( f1 -11.*(f2 - f3) - f4)*xd /(m_hx *8.) +
            ( f1 -    (f2 + f3) + f4)*x2/(m_hx2*4.) +
            (-f1 + 3.*(f2 - f3) + f4)*x3/(m_hx3*2.);
#endif
          //cubic type2 in x direction: for row j
          Real r2 = (-f5 + 9.*(f6 + f7) - f8)   /    16.  +
            ( f5 -11.*(f6 - f7) - f8)*xd /(m_hx *8.) +
            ( f5 -    (f6 + f7) + f8)*x2/(m_hx2*4.) +
            (-f5 + 3.*(f6 - f7) + f8)*x3/(m_hx3*2.);
          retval = r2;
#if CH_SPACEDIM==3
          //cubic type2 in x direction: for row j+1
          Real r3 = (-f9 + 9.*(f10 + f11) - f12)   /    16.  +
            ( f9 -11.*(f10 - f11) - f12)*xd /(m_hx *8.) +
            ( f9 -    (f10 + f11) + f12)*x2/(m_hx2*4.) +
            (-f9 + 3.*(f10 - f11) + f12)*x3/(m_hx3*2.);

          //cubic type2 in x direction: for row j+2
          Real r4 = (-f13 + 9.*(f14 + f15) - f16)   /    16.  +
            ( f13 -11.*(f14 - f15) - f16)*xd /(m_hx *8.) +
            ( f13 -    (f14 + f15) + f16)*x2/(m_hx2*4.) +
            (-f13 + 3.*(f14 - f15) + f16)*x3/(m_hx3*2.);

          //cubic type2 in y direction!!!
          retval  = (-r1 + 9.*(r2 + r3) - r4)   /    16.  +
            ( r1 -11.*(r2 - r3) - r4)*yd /(m_hy *8.) +
            ( r1 -    (r2 + r3) + r4)*y2/(m_hy2*4.) +
            (-r1 + 3.*(r2 - r3) + r4)*y3/(m_hy3*2.);
#endif
        }
      else
        {
          MayDay::Abort("Bad InterpType in DEMIF");
        }
    }

  retval -= a_point[SpaceDim-1];

  return retval;
}

BaseIF* DEMIF::newImplicitFunction() const
{
  DEMIF* demPtr = new DEMIF(*this);

  return static_cast<BaseIF*>(demPtr);
}

#include "NamespaceFooter.H"
