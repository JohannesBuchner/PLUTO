#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

// Flag everything as not defined or set
PatchPluto::PatchPluto()
{
  m_isDefined = false;
  m_isBoundarySet = false;
  m_isRiemannSet = false;
  m_isCurrentTimeSet = false;
  m_isCurrentBoxSet = false;
}

PatchPluto::~PatchPluto()
{
}

// Define this object and the boundary condition object
void PatchPluto::define(ProblemDomain& a_domain,
                          const Real&    a_dx,
                          const int&     a_level,
                          const int&     a_numGhost)
{

  // Store the domain and grid spacing
  m_domain = a_domain;
  m_dx = a_dx;
  m_level = a_level;
  m_numGhost = a_numGhost;

  m_isDefined = true;
}

// Factory method - this object is its own factory.  It returns a pointer
// to new PatchPluto object with its boundary condtions defined.
 PatchPluto* PatchPluto::new_patchPluto() const
{
   // Make the new object
   PatchPluto* retval =
     static_cast<PatchPluto*>(new PatchPluto());

   // Set the boundary and Riemann solver values
   retval->setBoundary(left_bound_side,
                       right_bound_side);
   retval->setRiemann(rsolver);

   // Return the new object
   return retval;
}

// Set the current physical time - used for time dependent boundary conditions
void PatchPluto::setCurrentTime(const Real& a_currentTime)
{
  m_currentTime = a_currentTime;
  m_isCurrentTimeSet = true;
}

// Set the current box for the advanceStep call
void PatchPluto::setCurrentBox(const Box& a_currentBox)
{
  m_currentBox = a_currentBox;
  m_isCurrentBoxSet = true;
}

// Set the values corresponding to boundary conditions
void PatchPluto::setBoundary(const int *leftbound,
                             const int *rightbound)
{

 for (int i = 0; i < 3; i++)
    {
      left_bound_side[i] =  leftbound[i];
     right_bound_side[i] = rightbound[i];
    }

 m_isBoundarySet = true;

}

// Set the riemann solver to be used
void PatchPluto::setRiemann(Riemann_Solver *solver)
{

   rsolver = solver; 
   m_isRiemannSet = true;

}

// Set the grid structure to be passed to updateState
// (or advanceSolution and SWEEP in the next future)
void PatchPluto::setGrid(const Box&    a_box, Grid *grid, FArrayBox& a_dV)
{
  int    i, j, k, idim;
  int    np_int, np_tot, np_int_glob, np_tot_glob;
  double stretch, scrh, scrh1, scrh2, scrh3, dxmin[3];
  double ***dV[CHOMBO_NDV];

  for (int dir = 0; dir < SpaceDim; ++dir){  
    grid->level            = m_level;
    grid->np_tot[dir]      = a_box.size()[dir]+2*m_numGhost;
    grid->np_int[dir]      = a_box.size()[dir];
    grid->np_tot_glob[dir] = m_domain.domainBox().size()[dir]+2*m_numGhost;
    grid->np_int_glob[dir] = m_domain.domainBox().size()[dir];
    grid->nghost[dir]      = m_numGhost;
    grid->beg[dir]         = a_box.loVect()[dir]+m_numGhost;
    grid->end[dir]         = a_box.hiVect()[dir]+m_numGhost;
    grid->gbeg[dir]        = m_domain.domainBox().loVect()[dir]+m_numGhost;
    grid->gend[dir]        = m_domain.domainBox().hiVect()[dir]+m_numGhost;
    grid->lbeg[dir]        = m_numGhost;
    grid->lend[dir]        = grid->np_int[dir] + m_numGhost - 1;
   
    grid->lbound[dir] = 0;
    grid->rbound[dir] = 0;
   // Set flags to compute boundary conditions
   // is_gbeg (left boundary)
    Box tmp = a_box;
    tmp.shift(dir,-1);
    tmp &= m_domain;
    tmp.shift(dir,1);
    if (tmp != a_box) grid->lbound[dir] = 1;
   // is_gend (right_boundary)
    tmp = a_box;
    tmp.shift(dir,1);
    tmp &= m_domain;
    tmp.shift(dir,-1);
    if (tmp != a_box) grid->rbound[dir] = 1;

  }

  for (int dir = 0; dir < 3; ++dir){
    if (dir >= SpaceDim) {
      grid->np_int[dir] = grid->np_int_glob[dir] = 1;
      grid->np_tot[dir] = grid->np_tot_glob[dir] = 1;
      grid->nghost[dir] = 0;
      grid->beg[dir]    = grid->end[dir] = 0;
      grid->gbeg[dir]   = grid->gend[dir] = 0;
      grid->lbeg[dir]   = grid->lend[dir] = 0;
      grid->lbound[dir] = grid->rbound[dir] = 0;
    }
    np_tot_glob = grid->np_tot_glob[dir];
    np_int_glob = grid->np_int_glob[dir];
    np_tot = grid->np_tot[dir];
    np_int = grid->np_int[dir];

    grid->x[dir]  = ARRAY_1D(np_tot, double);
    grid->xr[dir] = ARRAY_1D(np_tot, double);
    grid->xl[dir] = ARRAY_1D(np_tot, double);
    grid->dx[dir] = ARRAY_1D(np_tot, double);

    stretch = 1.;
    if (dir == JDIR) stretch = g_x2stretch;
    if (dir == KDIR) stretch = g_x3stretch;

    for (i = 0; i < np_tot; i++){
      #if (CHOMBO_LOGR == YES)
      if (dir == IDIR) {
        double argl = Real(i + grid->beg[dir] - 2*grid->nghost[dir]);
        double argr = Real(i + grid->beg[dir] - 2*grid->nghost[dir] + 1);

        grid->xl[dir][i] = g_domBeg[IDIR]*exp(argl*m_dx);
        grid->xr[dir][i] = g_domBeg[IDIR]*exp(argr*m_dx);  
        grid->x[dir][i]  = 0.5*(grid->xr[dir][i] + grid->xl[dir][i]);
        grid->dx[dir][i] = grid->xr[dir][i] - grid->xl[dir][i];
      } else {
      #endif
        double argl = Real(i + grid->beg[dir] - 2*grid->nghost[dir]);
        grid->xl[dir][i] = argl*m_dx*stretch + g_domBeg[dir];
        grid->x[dir][i]  = grid->xl[dir][i] + 0.5*m_dx*stretch;
        grid->xr[dir][i] = grid->xl[dir][i] + m_dx*stretch;
        grid->dx[dir][i] = m_dx*stretch; 
      #if (CHOMBO_LOGR == YES)
      }
      #endif
    }
  }

  CH_assert(m_isBoundarySet);
  grid->lbound[IDIR] *=   left_bound_side[IDIR];
  grid->rbound[IDIR] *=  right_bound_side[IDIR];
  grid->lbound[JDIR] *=   left_bound_side[JDIR];
  grid->rbound[JDIR] *=  right_bound_side[JDIR];
  grid->lbound[KDIR] *=   left_bound_side[KDIR];
  grid->rbound[KDIR] *=  right_bound_side[KDIR];
  MakeGeometry(grid);

// compute cell volumes (dV/m_dx^3) and cylindrical radius

  #if GEOMETRY != CARTESIAN
  double mdx_dim = D_EXPAND(m_dx, *m_dx, *m_dx);

  NX1_TOT = grid->np_tot[IDIR];
  NX2_TOT = grid->np_tot[JDIR];
  NX3_TOT = grid->np_tot[KDIR];
  for (i = 0; i < CHOMBO_NDV; i++) {
    dV[i] = ArrayMap(NX3_TOT, NX2_TOT, NX1_TOT, a_dV.dataPtr(i));
  }

  for (k = 0; k < NX3_TOT; k++) {
  for (j = 0; j < NX2_TOT; j++) {
  for (i = 0; i < NX1_TOT; i++) {
    dV[0][k][j][i] = grid->dV[k][j][i]/mdx_dim;
    #if CHOMBO_CONS_AM == YES
    #if GEOMETRY == POLAR || GEOMETRY == CYLINDRICAL
    dV[1][k][j][i] = grid->x[IDIR][i];
    #else
    dV[1][k][j][i] = D_EXPAND(grid->x[IDIR][i], *= sin(grid->x[JDIR][j]);
    #endif
    #endif
  }}}


  for (i = 0; i < CHOMBO_NDV; i++) FreeArrayMap(dV[i]);
#endif /* != CARTESIAN */

// compute dl_min

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)

  grid->dl_min[IDIR] = m_dx;
  grid->dl_min[JDIR] = m_dx*g_x2stretch;
  grid->dl_min[KDIR] = m_dx*g_x3stretch;

#else

  IBEG = grid->lbeg[IDIR]; IEND = grid->lend[IDIR];
  JBEG = grid->lbeg[JDIR]; JEND = grid->lend[JDIR];
  KBEG = grid->lbeg[KDIR]; KEND = grid->lend[KDIR];

  for (idim = 0; idim < 3; idim++)  dxmin[idim] = 1.e30;

  for (i = IBEG; i <= IEND; i++) {
  for (j = JBEG; j <= JEND; j++) {
  for (k = KBEG; k <= KEND; k++) {

  #if (TIME_STEPPING == RK2)
    D_EXPAND(scrh1 = Length_1(i, j, k, grid); ,
             scrh2 = Length_2(i, j, k, grid); ,
             scrh3 = Length_3(i, j, k, grid); )

    scrh = D_EXPAND(1./scrh1, +1./scrh2, +1./scrh3);
    scrh = (double)DIMENSIONS/scrh;
    D_EXPAND(dxmin[IDIR] = MIN (dxmin[IDIR], scrh); ,
             dxmin[JDIR] = dxmin[IDIR]; ,
             dxmin[KDIR] = dxmin[IDIR]; )
  #else
    scrh = Length_1(i, j, k, grid);
    dxmin[IDIR] = MIN (dxmin[IDIR], scrh);

    scrh = Length_2(i, j, k, grid);
    dxmin[JDIR] = MIN (dxmin[JDIR], scrh);

    scrh = Length_3(i, j, k, grid);
    dxmin[KDIR] = MIN (dxmin[KDIR], scrh);
  #endif

  }}}

  grid->dl_min[IDIR] = dxmin[IDIR];
  grid->dl_min[JDIR] = dxmin[JDIR];
  grid->dl_min[KDIR] = dxmin[KDIR];

 #endif

}

// Loop over the patches of a level to assign initial conditions

void PatchPluto::initialize(LevelData<FArrayBox>& a_U)
{

 CH_assert(m_isDefined);

 DataIterator dit = a_U.boxLayout().dataIterator();

 // Iterator of all grids in this level
 for (dit.begin(); dit.ok(); ++dit)
 {
   // Storage for current grid
   FArrayBox& U = a_U[dit()];
   // Set up initial condition in this patch
   startup(U);

 }
   
}

// Number of conserved variables
int PatchPluto::numConserved()
{
   return NVAR;
}

// Generate default names for the conserved variables, "variable#"
Vector<string> PatchPluto::ConsStateNames()
{
  Vector<string> retval;
  int nv;
  char varNameChar[80];

  for (nv = 0; nv < NVAR; nv++){
    if (nv == RHO) retval.push_back("Density");
    if (nv == MX1) retval.push_back("X-momentum");
    if (nv == MX2) retval.push_back("Y-momentum");
    if (nv == MX3) retval.push_back("Z-momentum");
#if PHYSICS == MHD || PHYSICS == RMHD
    if (nv == BX1) retval.push_back("X-magnfield");
    if (nv == BX2) retval.push_back("Y-magnfield");
    if (nv == BX3) retval.push_back("Z-magnfield");
#endif
#if EOS != ISOTHERMAL
    if (nv == ENG) retval.push_back("energy-density");
#endif

#if NTRACER > 0   
    if (nv >= TRC && nv < (TRC + NTRACER)){
      sprintf(varNameChar,"rho*tracer%d",nv+1-TRC);
      retval.push_back(string(varNameChar));
    }
#endif
#if ENTROPY_SWITCH
    if (nv == ENTR) retval.push_back("rho*entropy");
#endif

#ifdef GLM_MHD 
    if (nv == PSI_GLM) retval.push_back("psi_glm");
#endif

#if COOLING == SNEq
    if (nv == X_HI) retval.push_back("rho*x_hi");
#endif

#if COOLING == H2_COOL
    if (nv >= X_HI && nv < NFLX+NIONS){
      sprintf(varNameChar,"rho*fX%d",nv-NFLX);
      retval.push_back(string(varNameChar));
    }
#endif

#if COOLING == MINEq
    if (nv >= X_HI && nv < NFLX+NIONS){
      sprintf(varNameChar,"rho*fX%d",nv-NFLX);
      retval.push_back(string(varNameChar));
    }
#endif
  }
  return retval; 
}

// Generate default names for the primitive variables, "variable#"
Vector<string> PatchPluto::PrimStateNames()
{
  Vector<string> retval;
  int nv;
  static Output output;
  
  if (output.var_name == NULL){
    for (nv = 0; nv < NVAR; nv++){
      output.var_name = ARRAY_2D(64,128,char);
    }
    SetDefaultVarNames(&output);
  }

  for (nv = 0; nv < NVAR; nv++){
    retval.push_back(output.var_name[nv]);
  }
  return retval; 
}

// Number of flux variables
int PatchPluto::numFluxes()
{
  // In some computations there may be more fluxes than conserved variables
  return NVAR;
}

// Number of primitive variables
int PatchPluto::numPrimitives()
{
  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return NVAR;
}

// Component index within the primitive variables of the pressure
int PatchPluto::pressureIndex()
{
  #if EOS != ISOTHERMAL
   return PRS;
  #else 
   return -1;
  #endif

}

// Return true if everything is defined and setup
bool PatchPluto::isDefined() const
{
  return m_isDefined        &&
         m_isBoundarySet    &&
         m_isRiemannSet     &&
         m_isCurrentTimeSet &&
         m_isCurrentBoxSet;
}
