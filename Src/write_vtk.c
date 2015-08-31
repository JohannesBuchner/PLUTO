/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Write data in VTK format.

  Collection of basic functions to write VTK files using the
  simple legacy format.
  Files can be written in serial or parallel mode and consists of the
  basic five parts :

  -# File version and identifier
  -# Header consisting of a string
  -# File format
  -# Dataset structure: describes the gwometry and topology of the 
     dataset. 
  -# Dataset attributes. This section is used to write the actual binary
     data as a vector or scalar data.

  The mesh topology and the variable datasets are written using single 
  precision (4 bytes) binary format.
  VTK file are usually written following big endian order. 
  Therefore, we swap endianity only if local architecture has little 
  endian ordering.

  The WriteVTK_Header() function provides the basic functionality for 
  steps 1, 2, 3 and 4. Only processor 0 does the actual writing.
  For cartesian/cylindrical geometries the default grid topology is 
  "RECTILINEAR_GRIDS" whereas for polar/spherical we employ 
  "STRUCTURED_GRID" to provide a convenient mapping to a cartesian mesh.\n
  <b> Note for 2D datasets</b>: in order to produce a strictly 2D 
  dataset we always set the third coordinate (x3) to zero. 
  For this reason, in 2D spherical cordinates we swap the role of the  
  "y" and "z" coordinates.
            
  The WriteVTK_Vector() is fully parallel and is used to write data with 
  the vector attribute (by default these include velocity and magnetic 
  fields).
   
  The WriteVTK_Scalar() is fully parallel and is used to write data with 
  the scalar attribute (by default these include density, pressure, tracer
  and user-defined variables).

  \b Reference

  http://www.vtk.org/VTK/img/file-formats.pdf

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 17, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define RECTILINEAR_GRID    14
#define STRUCTURED_GRID     35

#ifndef VTK_FORMAT
  #if GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL
    #define VTK_FORMAT  RECTILINEAR_GRID
  #else
    #define VTK_FORMAT  STRUCTURED_GRID
  #endif
#endif

#ifndef VTK_TIME_INFO       /* Enable this keyword only if you're using */
 #define VTK_TIME_INFO  NO  /* VisIt to display data results            */
#endif

/* ---------------------------------------------------------
    The following macros are specific to this file only 
    and are used to ease up serial/parallel implementation
    for writing strings and real arrays 
   --------------------------------------------------------- */   
    
#ifdef PARALLEL
 #define VTK_HEADER_WRITE_STRING(header) \
         AL_Write_header (header, strlen(header), MPI_CHAR, SZ_Float_Vect);
 #define VTK_HEADER_WRITE_FLTARR(arr, nelem) \
         AL_Write_header (arr, nelem, MPI_FLOAT, SZ_Float_Vect);
 #define VTK_HEADER_WRITE_DBLARR(arr, nelem) \
         AL_Write_header (arr, nelem, MPI_DOUBLE, SZ_Float_Vect);
#else
 #define VTK_HEADER_WRITE_STRING(header) \
         fprintf (fvtk,header);
 #define VTK_HEADER_WRITE_FLTARR(arr,nelem) \
         fwrite(arr, sizeof(float), nelem, fvtk);
 #define VTK_HEADER_WRITE_DBLARR(arr,nelem) \
         fwrite(arr, sizeof(double), nelem, fvtk);
#endif

/* ********************************************************************* */
void WriteVTK_Header (FILE *fvtk, Grid *grid)
/*!
 * Write VTK header in parallel or serial mode.
 * In parallel mode only processor 0 does the actual writing 
 * (see al_io.c/AL_Write_header).
 *
 *
 * \param [in]  fvtk  pointer to file
 * \param [in]  grid  pointer to an array of Grid structures
 *
 * \todo  Write the grid using several processors. 
 *********************************************************************** */
{
  int    i, j, k;
  int    nx1, nx2, nx3;
  char   header[1024];
  float  x1, x2, x3;
  static float  ***node_coord, *xnode, *ynode, *znode;

/* ------------------------------------------------
          get global dimensions
   ------------------------------------------------ */

  nx1 = grid[IDIR].gend + 1 - grid[IDIR].nghost;
  nx2 = grid[JDIR].gend + 1 - grid[JDIR].nghost;
  nx3 = grid[KDIR].gend + 1 - grid[KDIR].nghost;

/* -------------------------------------------------------------
     Allocate memory and define node coordinates only once.
   ------------------------------------------------------------- */

  if (node_coord == NULL){
    node_coord = ARRAY_3D(nx2 + JOFFSET, nx1 + IOFFSET, 3, float);

    #if VTK_FORMAT == RECTILINEAR_GRID
     xnode = ARRAY_1D(nx1 + IOFFSET, float);
     ynode = ARRAY_1D(nx2 + JOFFSET, float);
     znode = ARRAY_1D(nx3 + KOFFSET, float);

     for (i = 0; i < nx1 + IOFFSET; i++){
       x1 = (float)(grid[IDIR].xl_glob[i+IBEG]);
       if (IsLittleEndian()) SWAP_VAR(x1);
       xnode[i] = x1;
     }
     for (j = 0; j < nx2 + JOFFSET; j++){
       x2 = (float)(grid[JDIR].xl_glob[j+JBEG]);
       if (IsLittleEndian()) SWAP_VAR(x2);
       ynode[j] = x2;
     }
     for (k = 0; k < nx3 + KOFFSET; k++){
       x3 = (float)(grid[KDIR].xl_glob[k+KBEG]);
       if (IsLittleEndian()) SWAP_VAR(x3);
       #if DIMENSIONS == 2
        znode[k] = 0.0;
       #else
        znode[k] = x3;
       #endif
     }
    #endif
  }

/* ----------------------------------------------------------
    Part I, II, III: 
    Write file header on string "header" 
   ---------------------------------------------------------- */

  sprintf(header,"# vtk DataFile Version 2.0\n"); 
  sprintf(header+strlen(header),"PLUTO %s VTK Data\n",PLUTO_VERSION);
  sprintf(header+strlen(header),"BINARY\n");
  #if VTK_FORMAT == RECTILINEAR_GRID
   sprintf(header+strlen(header),"DATASET %s\n","RECTILINEAR_GRID");
  #elif VTK_FORMAT == STRUCTURED_GRID
   sprintf(header+strlen(header),"DATASET %s\n","STRUCTURED_GRID");
  #endif

  VTK_HEADER_WRITE_STRING(header);

  /* -- generate time info -- */

  #if VTK_TIME_INFO == YES
   sprintf (header,"FIELD FieldData 1\n");
   sprintf (header+strlen(header),"TIME 1 1 double\n");
   double tt=g_time;
   if (IsLittleEndian()) SWAP_VAR(tt);
   VTK_HEADER_WRITE_STRING(header);
   VTK_HEADER_WRITE_DBLARR(&tt, 1);
   VTK_HEADER_WRITE_STRING("\n");

  #endif /* VTK_TIME_INFO */

  sprintf(header,"DIMENSIONS %d %d %d\n",
                  nx1 + IOFFSET, nx2 + JOFFSET, nx3 + KOFFSET);
  VTK_HEADER_WRITE_STRING(header);
  
#if VTK_FORMAT == RECTILINEAR_GRID

  /* -- reset header string and keep going -- */

   sprintf(header,"X_COORDINATES %d float\n", nx1 + IOFFSET);
   VTK_HEADER_WRITE_STRING(header);
   VTK_HEADER_WRITE_FLTARR(xnode, nx1 + IOFFSET);

   sprintf(header,"\nY_COORDINATES %d float\n", nx2 + JOFFSET);
   VTK_HEADER_WRITE_STRING(header);
   VTK_HEADER_WRITE_FLTARR(ynode, nx2 + JOFFSET);

   sprintf(header,"\nZ_COORDINATES %d float\n", nx3 + KOFFSET);
   VTK_HEADER_WRITE_STRING(header);
   VTK_HEADER_WRITE_FLTARR(znode, nx3 + KOFFSET);

   sprintf (header,"\nCELL_DATA %d\n", nx1*nx2*nx3);
   VTK_HEADER_WRITE_STRING (header);

#elif VTK_FORMAT == STRUCTURED_GRID

  sprintf(header,"POINTS %d float\n", (nx1+IOFFSET)*(nx2+JOFFSET)*(nx3+KOFFSET));
  VTK_HEADER_WRITE_STRING(header);

/* ---------------------------------------------------------------
    Part IV: (structured) grid information 
   --------------------------------------------------------------- */

  x1 = x2 = x3 = 0.0;
  for (k = 0; k < nx3 + KOFFSET; k++){
    for (j = 0; j < nx2 + JOFFSET; j++){ 
    for (i = 0; i < nx1 + IOFFSET; i++){
      D_EXPAND(x1 = grid[IDIR].xl_glob[IBEG + i];  ,
               x2 = grid[JDIR].xl_glob[JBEG + j];  ,
               x3 = grid[KDIR].xl_glob[KBEG + k];)
       
      #if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)
       node_coord[j][i][0] = x1;
       node_coord[j][i][1] = x2;
       node_coord[j][i][2] = x3;
      #elif GEOMETRY == POLAR
       node_coord[j][i][0] = x1*cos(x2);
       node_coord[j][i][1] = x1*sin(x2);
       node_coord[j][i][2] = x3;
      #elif GEOMETRY == SPHERICAL
       #if DIMENSIONS == 2
        node_coord[j][i][0] = x1*sin(x2);
        node_coord[j][i][1] = x1*cos(x2);
        node_coord[j][i][2] = 0.0;
       #elif DIMENSIONS == 3
        node_coord[j][i][0] = x1*sin(x2)*cos(x3);
        node_coord[j][i][1] = x1*sin(x2)*sin(x3);
        node_coord[j][i][2] = x1*cos(x2);
       #endif
      #endif
      
      if (IsLittleEndian()){
        SWAP_VAR(node_coord[j][i][0]);
        SWAP_VAR(node_coord[j][i][1]);
        SWAP_VAR(node_coord[j][i][2]);
      }
    }}
    VTK_HEADER_WRITE_FLTARR(node_coord[0][0],3*(nx1+IOFFSET)*(nx2+JOFFSET));
  }

  sprintf (header,"\nCELL_DATA %d\n", nx1*nx2*nx3);
  VTK_HEADER_WRITE_STRING(header);

#endif
}
#undef STRUCTERED_GRID    
#undef RECTILINEAR_GRID  

/* ********************************************************************* */
void WriteVTK_Vector (FILE *fvtk, Data_Arr V, double unit,
                      char *var_name, Grid *grid)
/*!
 * Write VTK vector field data.
 * This is enabled only when VTK_VECTOR_DUMP is set to \c YES.
 * For generality purposes, vectors are written always with 3 
 * components, even when there're only 2 being used.
 *
 * The following Maple script has been used to find vector 
 * components from cyl/sph to cartesian:
 *
 * \code
   restart;
   with(linalg);
   Acyl := matrix (3,3,[ cos(phi), sin(phi),  0,
                      -sin(phi), cos(phi),  0,
                        0     ,    0    , 1]);
   Asph := matrix (3,3,[ sin(theta)*cos(phi), sin(theta)*sin(phi),  cos(theta),
                   cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta),
                   -sin(phi)          , cos(phi)           , 0]);
   Bcyl := simplify(inverse(Acyl));
   Bsph := simplify(inverse(Asph));
 * \endcode
 *
 * \param [in]  fvtk    pointer to file
 * \param [in]  V       a 4D array [nv][k][j][i] containing the vector
 *                      components (nv) defined at cell centers (k,j,i).
 *                      The index nv = 0 marks the vector first component.
 * \param [in] unit     the corresponding cgs unit (if specified, 1 otherwise)
 * \param [in] var_name the variable name appearing in the VTK file
 * \param [in]    grid  pointer to an array of Grid structures
 *********************************************************************** */
{
  int i,j,k;
  int vel_field, mag_field;
  char header[512];
  static Float_Vect ***vect3D;
  double v[3], x1, x2, x3;

  if (vect3D == NULL){
    vect3D = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, Float_Vect);
  }

/* --------------------------------------------------------
               Write VTK vector fields 
   -------------------------------------------------------- */

  v[0] = v[1] = v[2] = 0.0;
  x1 = x2 = x3 = 0.0;

  vel_field = (strcmp(var_name,"vx1") == 0);
  mag_field = (strcmp(var_name,"bx1") == 0);
  if (vel_field || mag_field) { 
    DOM_LOOP(k,j,i){ 
      D_EXPAND(v[0] = V[0][k][j][i]; x1 = grid[IDIR].x[i]; ,
               v[1] = V[1][k][j][i]; x2 = grid[JDIR].x[j]; ,
               v[2] = V[2][k][j][i]; x3 = grid[KDIR].x[k];)
   
      VectorCartesianComponents(v, x1, x2, x3);
      vect3D[k][j][i].v1 = (float)v[0]*unit;
      vect3D[k][j][i].v2 = (float)v[1]*unit;
      vect3D[k][j][i].v3 = (float)v[2]*unit;

      if (IsLittleEndian()){
        SWAP_VAR(vect3D[k][j][i].v1);
        SWAP_VAR(vect3D[k][j][i].v2);
        SWAP_VAR(vect3D[k][j][i].v3);
      }

    } /* endfor DOM_LOOP(k,j,i) */

    if (vel_field)
      sprintf (header,"\nVECTORS %dD_Velocity_Field float\n", DIMENSIONS);
    else
      sprintf (header,"\nVECTORS %dD_Magnetic_Field float\n", DIMENSIONS);

    VTK_HEADER_WRITE_STRING(header);
    WriteBinaryArray (vect3D[0][0], sizeof(Float_Vect), SZ_Float_Vect, fvtk, -1);
  }
}

/* ********************************************************************* */
void WriteVTK_Scalar (FILE *fvtk, double ***V, double unit,
                      char *var_name, Grid *grid)
/*!
 * Write VTK scalar field.
 * 
 * \param [in]   fvtk       pointer to file (handle)
 * \param [in]   V          pointer to 3D data array
 * \param [in] unit     the corresponding cgs unit (if specified, 1 otherwise)
 * \param [in]   var_name   the variable name appearing in the VTK file
 * \param [in]   grid       pointer to an array of Grid structures
 *********************************************************************** */
{
  int i,j,k;
  char header[512];
  float ***Vflt;

  sprintf (header,"\nSCALARS %s float\n", var_name);
  sprintf (header,"%sLOOKUP_TABLE default\n",header);

  #ifdef PARALLEL
   MPI_Barrier (MPI_COMM_WORLD);
   AL_Write_header (header, strlen(header), MPI_CHAR, SZ_float);
  #else
   fprintf (fvtk, "%s",header);
  #endif

  Vflt = Convert_dbl2flt(V, unit, IsLittleEndian());
  WriteBinaryArray (Vflt[0][0], sizeof(float), SZ_float, fvtk, -1);
}
