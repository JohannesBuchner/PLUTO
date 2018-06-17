; ----------------------------------------------------------------------
;  PLOAD.PRO: main driver to read different data formats (flt, dbl, 
;             flt.h5, dbl.h5, hdf5 written by the PLUTO code.
;
;  LAST MODIFIED:  Apr 25, 2015
; ----------------------------------------------------------------------

COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
                    dx1,dx2,dx3, geometry, xpos, ypos,$
                    AMRLevel, AMRBoxes;  ** Chombo data structure **
                                      ;  ** loaded when HDF5LOAD is called **

COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
                             Ax1, Ax2, Ax3, $
                             Bx1, Bx2, Bx3, $
                             Bx1s, Bx2s, Bx3s,$
                             Ex1, Ex2, Ex3,qg,$
                         ; ----------------------------------------
                             v1, v2, v3, $   ; Kept for backward
                             b1, b2, b3, $   ; compatibility with 
                             A1, A2, A3, $   ; PLUTO 3 data name
                             b1s, b2s, b3s, $ ;
                             pr,            $ ;
                         ; -----------------------------------------
                  prs, phi_glm, psi_glm, $
                  tr1, tr2, tr3, tr4, $
                  x_HI, x_HII, x_H2, x_heI, x_heII,$
                  x_cI, x_cII, x_cIII, x_cIV, x_cV,$
                  x_nI, x_nII, x_nIII, x_nIV, x_nV,$
                  x_oI, x_oII, x_oIII, x_oIV, x_oV,$
                  x_neI, x_neII, x_neIII, x_neIV, x_neV,$
                  x_sI, x_sII, x_sIII, x_sIV, x_sV, $
                  rho_d, vx1_d, vx2_d, vx3_d, $  ; Dust
                  particles, particles_reduced, $ ; Particles
                  vars, vname
COMMON PLUTO_RUN,  t, dt, nlast, first_call

; *************************************************************************
;
;   Match a variable name to its storage area, e.g.,
;   variable named "rho" should be assigned to the
;   rho variable contained in the PLUTO_VAR common block
;
; *************************************************************************
PRO MATCH_VARNAME, vpt, name, silent=silent

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

 match  = 0
 match3 = 0

 ; ----------------------------------------------------------------
 ; Backward compatibility for vector names (PLUTO 3.xx or earlier)
 ; ----------------------------------------------------------------

 CASE name OF
   "v1":  BEGIN v1  = vpt & match3 = 1 & END
   "v2":  BEGIN v2  = vpt & match3 = 1 & END
   "v3":  BEGIN v3  = vpt & match3 = 1 & END
   "b1":  BEGIN b1  = vpt & match3 = 1 & END
   "b2":  BEGIN b2  = vpt & match3 = 1 & END
   "b3":  BEGIN b3  = vpt & match3 = 1 & END
   "b1s": BEGIN b1s = vpt & match3 = 1 & END
   "b2s": BEGIN b2s = vpt & match3 = 1 & END
   "b3s": BEGIN b3s = vpt & match3 = 1 & END
   "A1":  BEGIN A1  = vpt & match3 = 1 & END
   "A2":  BEGIN A2  = vpt & match3 = 1 & END
   "A3":  BEGIN A3  = vpt & match3 = 1 & END
   "pr":  BEGIN pr  = vpt & match3 = 1 & END
 ELSE:
 ENDCASE

 ; --------------------------------------------
 ;    PLUTO 4 standard variable names
 ; --------------------------------------------

 CASE name OF
   "rho":  BEGIN rho  = vpt & match = 1 & END
   "vx1":  BEGIN vx1  = vpt & match = 1 & END
   "vx2":  BEGIN vx2  = vpt & match = 1 & END
   "vx3":  BEGIN vx3  = vpt & match = 1 & END
   "Ax1":  BEGIN Ax1  = vpt & match = 1 & END
   "Ax2":  BEGIN Ax2  = vpt & match = 1 & END
   "Ax3":  BEGIN Ax3  = vpt & match = 1 & END
   "Bx1":  BEGIN Bx1  = vpt & match = 1 & END
   "Bx2":  BEGIN Bx2  = vpt & match = 1 & END
   "Bx3":  BEGIN Bx3  = vpt & match = 1 & END
   "bx1":  BEGIN Bx1  = vpt & match = 1 & END
   "bx2":  BEGIN Bx2  = vpt & match = 1 & END
   "bx3":  BEGIN Bx3  = vpt & match = 1 & END
   "Ex1":  BEGIN Ex1  = vpt & match = 1 & END
   "Ex2":  BEGIN Ex2  = vpt & match = 1 & END
   "Ex3":  BEGIN Ex3  = vpt & match = 1 & END
   "qg":   BEGIN qg   = vpt & match = 1 & END
   "Bx1s": BEGIN Bx1s = vpt & match = 1 & END
   "Bx2s": BEGIN Bx2s = vpt & match = 1 & END
   "Bx3s": BEGIN Bx3s = vpt & match = 1 & END
   "psi_glm": BEGIN psi_glm = vpt & match = 1 & END
   "phi_glm": BEGIN phi_glm = vpt & match = 1 & END

   "prs": BEGIN prs = vpt & match = 1 & END
   "tr1": BEGIN tr1 = vpt & match = 1 & END
   "tr2": BEGIN tr2 = vpt & match = 1 & END
   "tr3": BEGIN tr3 = vpt & match = 1 & END
   "tr4": BEGIN tr4 = vpt & match = 1 & END
 
  ; -- Ion fractions --

   "X_HI":   BEGIN x_hI   = vpt & match = 1 & END
   "X_HII":  BEGIN x_hII  = vpt & match = 1 & END
   "X_H2":   BEGIN x_h2   = vpt & match = 1 & END
   "X_HeI":  BEGIN x_heI  = vpt & match = 1 & END
   "X_HeII": BEGIN x_heII = vpt & match = 1 & END
   "X_CI":   BEGIN x_cI   = vpt & match = 1 & END
   "X_CII":  BEGIN x_cII  = vpt & match = 1 & END
   "X_CIII": BEGIN x_cIII = vpt & match = 1 & END
   "X_CIV":  BEGIN x_cIV  = vpt & match = 1 & END
   "X_CV":   BEGIN x_cV   = vpt & match = 1 & END

   "X_NI":   BEGIN x_nI   = vpt & match = 1 & END
   "X_NII":  BEGIN x_nII  = vpt & match = 1 & END
   "X_NIII": BEGIN x_nIII = vpt & match = 1 & END
   "X_NIV":  BEGIN x_nIV  = vpt & match = 1 & END
   "X_NV":   BEGIN x_nV   = vpt & match = 1 & END

   "X_OI":   BEGIN x_oI   = vpt & match = 1 & END
   "X_OII":  BEGIN x_oII  = vpt & match = 1 & END
   "X_OIII": BEGIN x_oIII = vpt & match = 1 & END
   "X_OIV":  BEGIN x_oIV  = vpt & match = 1 & END
   "X_OV":   BEGIN x_oV   = vpt & match = 1 & END

   "X_NeI":   BEGIN x_neI   = vpt & match = 1 & END
   "X_NeII":  BEGIN x_neII  = vpt & match = 1 & END
   "X_NeIII": BEGIN x_neIII = vpt & match = 1 & END
   "X_NeIV":  BEGIN x_neIV  = vpt & match = 1 & END
   "X_NeV":   BEGIN x_neV   = vpt & match = 1 & END

   "X_SI":   BEGIN x_sI   = vpt & match = 1 & END
   "X_SII":  BEGIN x_sII  = vpt & match = 1 & END
   "X_SIII": BEGIN x_sIII = vpt & match = 1 & END
   "X_SIV":  BEGIN x_sIV  = vpt & match = 1 & END
   "X_SV":   BEGIN x_sV   = vpt & match = 1 & END

  ; -- Dust --

   "rho_d":  BEGIN rho_d  = vpt & match = 1 & END
   "vx1_d":  BEGIN vx1_d  = vpt & match = 1 & END
   "vx2_d":  BEGIN vx2_d  = vpt & match = 1 & END
   "vx3_d":  BEGIN vx3_d  = vpt & match = 1 & END

 ELSE:
 ENDCASE

 ; ---- Do dome printing now ---- 

 IF (match eq 1) THEN BEGIN
   IF (NOT KEYWORD_SET(SILENT)) THEN PRINT,"> Reading ",name
 ENDIF

 IF (match3 eq 1) THEN BEGIN
   IF (NOT KEYWORD_SET(SILENT)) THEN PRINT,"> Reading ",name, " (PLUTO-3 Data)"
 ENDIF

; -- Now check user-defined variables --

  CATCH, err_status                ; Error handler just in case
  IF (err_status NE 0) THEN BEGIN  ; MATCH_USERDEF_VARNAME does not exist
    vpt = 0; free memory
    RETURN
  ENDIF

  MATCH_USERDEF_VARNAME,vpt, name, silent=silent
  
 vpt = 0; free memory

END

;+
; NAME: PRINT_MEMORY_USAGE
;
; PURPOSE: print memory usage in megabytes
;
;-

PRO PRINT_MEMORY_USAGE, SILENT=SILENT

  IF (KEYWORD_SET (SILENT)) THEN RETURN
  
  PRINT, "> Used memory:", $
    STRCOMPRESS(STRING(MEMORY(/CURRENT,/L64)/1024.d0/1024.d0))," Mb"
END

;************************************************************
;
; NAME:       READ_PLUTO_GRID
;
; AUTHOR:     A. Mignone
;
; LAST MODIFIED:      Aug 26, 2012
;
; SYNOPSIS:   read_pluto_grid
;
; DESCRIPTION:  load grid and geometry information;
;************************************************************
PRO READ_PLUTO_GRID, dir, silent=silent

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    read grid files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 OPENR, U, dir+'grid.out', /GET_LUN, ERROR=err
 IF (NOT (err EQ 0)) THEN BEGIN
   print,"! cannot find "+DIR+"grid.out"
   RETURN
 ENDIF

; -----------------------------------------------------------------
;  read the very first character: 
;  if it matches "#" then we use the PLUTO 4.0 grid data file
;  otherwise we stick to the old format
; -----------------------------------------------------------------

 s = ''
 READF, U, s
 id = STRMID(s,0,1)
 IF (id EQ '#') THEN BEGIN; reading PLUTO 4 Grid file

   IF (NOT KEYWORD_SET(SILENT)) THEN PRINT,"> Using PLUTO 4 Grid file"
   WHILE (id EQ '#') DO BEGIN
     POINT_LUN, -U, pos
     READF, U, s 
     id = STRMID(s, 0, 1)
   ENDWHILE
   POINT_LUN, U, pos
   
   j = LONG64(0) & xL = DOUBLE(0.0) & xR = DOUBLE(0.0)

   READF, U, fn1
   NX1 = LONG64(fn1) &  x1  = DBLARR(NX1) &  dx1 = DBLARR(NX1)
   FOR i = 0L, NX1 - 1L DO BEGIN
     READF, U, j, xL, xR
     x1(i)  = 0.5d0*(xL + xR)
     dx1(i) = xR - xL
   ENDFOR

   READF, U, fn2
   NX2 = LONG64(fn2) &  x2  = DBLARR(NX2) &  dx2 = DBLARR(NX2)
   FOR i = 0L, NX2 - 1L DO BEGIN
     READF, U, j, xL, xR
     x2(i)  = 0.5d0*(xL + xR)
     dx2(i) = xR - xL
   ENDFOR

   READF, U, fn3
   NX3 = LONG64(fn3) & x3 = DBLARR(NX3) & dx3 = DBLARR(NX3)
   FOR i = 0L, NX3 - 1L DO BEGIN
     READF, U, j, xL, xR
     x3(i)  = 0.5d0*(xL + xR)
     dx3(i) = xR - xL
   ENDFOR

 ENDIF ELSE BEGIN; reading PLUTO 3 Grid data file

   PRINT,"> Using PLUTO 3.0 Grid file"
   fn1 = FLOAT(s)
   NX1 = LONG64(fn1) &  x1  = FLTARR(NX1) &  dx1 = FLTARR(NX1)
   FOR i = 0L,NX1-1L DO BEGIN
     READF, U, j,xl_0,xc_0,xr_0,dx_0
     x1(i)  = xc_0
     dx1(i) = dx_0
   ENDFOR

   READF, U, fn2
   NX2  = LONG64(fn2) & x2  = FLTARR(NX2) & dx2 = FLTARR(NX2)
   for i = 0L,NX2-1L DO BEGIN
     READF, U,j,xl_0,xc_0,xr_0,dx_0
     x2(i)  = xc_0
     dx2(i) = dx_0
   ENDFOR

   READF, U, fn3
   NX3  = LONG64(fn3) &  x3  = FLTARR(NX3) & dx3 = FLTARR(NX3)
   FOR i = 0L,NX3-1L DO BEGIN
     READF, U,j,xl_0,xc_0,xr_0,dx_0
     x3(i)  = xc_0
     dx3(i) = dx_0
   ENDFOR
 ENDELSE; read pluto grid

 IF (NOT KEYWORD_SET (silent)) THEN BEGIN
   xb = x1(0)-0.5*dx1(0) & xe = x1(NX1-1) + 0.5*dx1(NX1-1)
   yb = x2(0)-0.5*dx2(0) & ye = x2(NX2-1) + 0.5*dx2(NX2-1)
   zb = x3(0)-0.5*dx3(0) & ze = x3(NX3-1) + 0.5*dx3(NX3-1)

   domstring = "["+ARG2STR(xb)+", "+ARG2STR(xe)+"]"
   resstring = ARG2STR(NX1)

   IF (NX2 GT 1) THEN BEGIN
     domstring += " x ["+ARG2STR(yb)+", "+ARG2STR(ye)+"]"
     resstring += " x "+ARG2STR(NX2)
   ENDIF

   IF (NX3 GT 1) THEN BEGIN
     domstring += " x ["+ARG2STR(zb)+", "+ARG2STR(ze)+"]"
     resstring += " x "+ARG2STR(NX3)
   ENDIF

   print,"> Domain: ",domstring
   print,"> Size:   ",resstring
 ENDIF
 CLOSE,U
 FREE_LUN,U

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read definitions.h if present
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 OPENR, U, DIR+'definitions.h', /GET_LUN, ERROR=err

 IF (NOT (err EQ 0)) THEN BEGIN
   IF (KEYWORD_SET(verbose)) THEN print,"! Can not find definitions.h, skipping geometry definition..."
 ENDIF ELSE BEGIN
   str = ' '
   FOR i = 0,30 DO BEGIN
     readf, U,format='(A)',str
     str_el = strsplit(str,' ',/extract)
     IF (str_el(1) EQ 'GEOMETRY') THEN GOTO, DEF_GEOMETRY
   ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     define geometry
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 DEF_GEOMETRY:
   IF (str_el(2) EQ 'CARTESIAN')   THEN BEGIN
     geometry = 'cartesian'
     x = x1
     y = x2
     z = x3
     IF (NOT KEYWORD_SET (silent)) THEN BEGIN
       print,"> Geometry:    Cartesian (x = x1; y = x2; z = x3)"
     ENDIF
   ENDIF

   IF (str_el(2) EQ 'CYLINDRICAL') THEN BEGIN
     geometry = 'cylindrical'
     r   = x1
     z   = x2
     phi = x3
     IF (NOT KEYWORD_SET (silent)) THEN BEGIN
       print,"> Geometry:    Cylindrical (r = x1; z = x2; phi = x3)"
     ENDIF
   ENDIF

   IF (str_el(2) EQ 'SPHERICAL') THEN BEGIN
     geometry = 'spherical'
     r     = x1
     theta = x2

     xpos = fltarr(NX1,NX2)  ; define projection on cartesian coordinates
     ypos = fltarr(NX1,NX2)

     for i = 0,NX1-1 do begin
     for j = 0,NX2-1 do begin
       xpos(i,j) = (x1(i)-0.0*dx1(i))*sin(x2(j)-0.0*dx2(j))
       ypos(i,j) = (x1(i)-0.0*dx1(i))*cos(x2(j)-0.0*dx2(j))
     endfor
     endfor

     IF (NOT KEYWORD_SET (silent)) THEN BEGIN
       print,"> Geometry:    Spherical (xpos = r*sin(theta); ypos = r*cos(theta))"
     ENDIF
   ENDIF

   IF (str_el(2) EQ 'POLAR') THEN BEGIN
     geometry = 'polar'
     r   = x1
     phi = x2
     xpos = fltarr(NX1,NX2)  ; define projection on cartesian coordinates
     ypos = fltarr(NX1,NX2)

     for i = 0,NX1-1 do begin
     for j = 0,NX2-1 do begin
       xpos(i,j) = x1(i)*cos(x2(j))
       ypos(i,j) = x1(i)*sin(x2(j))
     endfor
     endfor

     IF (NOT KEYWORD_SET (silent)) THEN BEGIN
       print,"> Geometry:    Polar (xpos = r*cos(phi); ypos = r*sin(phi))"
     ENDIF
   ENDIF

 ENDELSE
 CLOSE,U
 FREE_LUN,U
 
END

;************************************************************
;  NAME:      REBIN3D
;
;  AUTHOR:    Andrea Mignone (mignone@ph.unito.it)
;
;  PURPOSE:   do fast shrinking of a 3D dataset
;
;  SYNTAX:    REBIN3D, q, r
;
;  ARGUMENTS:
;
;    q = a 3-D array.
;
;    r = an integer number > 1 telling how many zones have to
;        lumped together
;
;  LAST MODIFIED:  Nov 9, 2009
;
;************************************************************
PRO REBIN3D, q, r

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

 nxx = NX1/r
 nyy = NX2/r
 nzz = NX3/r
 qnew = fltarr(nxx,nyy,nzz)

 print,"> Rebinning to ["+ARG2STR(nxx)+","+ARG2STR(nyy)+","+ARG2STR(nzz)+$
       "]...please wait"

; **** create grid in z direction ****

 IF ( ( (NX1 MOD r)+(NX2 MOD r)+(NX3 MOD r) ) EQ 0) THEN BEGIN
   FOR k = 0, nzz - 1 DO qnew(*,*,k) = REBIN(q(k), nxx, nyy, 1)
 ENDIF ELSE BEGIN
   FOR k = 0, nzz - 1 DO BEGIN
     scrh = q(k)
     qnew(*,*,k) = CONGRID(scrh, nxx, nyy, 1, /CENTER)
   ENDFOR
 ENDELSE

 q  = qnew
END

;+
; NAME:      ARRAY_EXTRACT
;
; AUTHOR:    Andrea Mignone (mignone@ph.unito.it)
;
; PURPOSE:   extract a sub-array from a 2 or 3D data set
;
; SYNTAX:    ARRAY_EXTRACT, q, x1range, x2range, x3range, indx[, /SILENT][, /VECTOR]
;
; ARGUMENT:
;
;   q          an associated file variable with a 2D structure,
;              e.g., q[k] = f(i,j). 
;
;   x1range    a 2D vector [min,max] specifying the x1 bounds of the 
;              sub-domain to be extracted.
;
;   x2range    a 2D vector [min,max] specifying the x2 bounds of the 
;              sub-domain to be extracted.
;
;   x3range    a 2D vector [min,max] specifying the x1 bounds of the 
;              sub-domain to be extracted.
;
;   indx       a 6 element integer array containing the grid indices of
;              the extracted array [ibeg, iend, ..., kbeg, kend]
;
; KEYWORDS: 
;
;   /SILENT      suppress output
;
;   /VECTOR      When this keyword is supplied, the input associated file variable
;                is assumed to have a 3D structure of the type q[k] = f(c, i, j) where 
;                'c' specify the vector component (useful for VTK extraction).
;
; LAST MODIFIED:
;
;   August 1st, 2013 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO ARRAY_EXTRACT, q, x1range, x2range, x3range, indx, SILENT=SILENT, VECTOR=VECTOR

 COMMON PLUTO_GRID
 COMMON PLUTO_VAR
 COMMON PLUTO_RUN

; ---------------------------------------------------
; Transfrom from spatial coordinate to grid indexes
; ---------------------------------------------------

 xbeg = MAX([x1range(0),MIN(x1)]) & xend = MIN([x1range(1),MAX(x1)])
 ybeg = MAX([x2range(0),MIN(x2)]) & yend = MIN([x2range(1),MAX(x2)])
 zbeg = MAX([x3range(0),MIN(x3)]) & zend = MIN([x3range(1),MAX(x3)])

 ibeg = FIX(MAX(WHERE(x1 LE xbeg))) & iend = FIX(MIN(WHERE(x1 GE xend)))
 jbeg = FIX(MAX(WHERE(x2 LE ybeg))) & jend = FIX(MIN(WHERE(x2 GE yend)))
 kbeg = FIX(MAX(WHERE(x3 LE zbeg))) & kend = FIX(MIN(WHERE(x3 GE zend)))

 nxx = iend - ibeg + 1
 nyy = jend - jbeg + 1
 nzz = kend - kbeg + 1

 IF (NOT KEYWORD_SET(SILENT)) THEN BEGIN
   print,"> Extracting array: [", ARG2STR(ibeg), ":",ARG2STR(iend),$
                            "][", ARG2STR(jbeg), ":",ARG2STR(jend),$
                            "][", ARG2STR(kbeg), ":",ARG2STR(kend),"]",$
         "  (size: ",ARG2STR(nxx)," x ", ARG2STR(nyy)," x ", ARG2STR(nzz),")"
 ENDIF
 
 IF (KEYWORD_SET(VECTOR)) THEN BEGIN
   
   qnew = FLTARR(3, nxx, nyy, nzz) 
   FOR k = 0, nzz-1 DO BEGIN
     scrh = q(k+kbeg)
     qnew[0:2,*,*,k] = scrh[0:2,ibeg:iend,jbeg:jend]
   ENDFOR
   
 ENDIF ELSE BEGIN
 
   qnew = FLTARR(nxx,nyy,nzz)
   FOR k = 0, nzz-1 DO BEGIN
     scrh = q(k+kbeg)
     qnew[*,*,k] = scrh[ibeg:iend,jbeg:jend]
   ENDFOR
   
 ENDELSE
 
 q  = qnew
 indx = [ibeg, jbeg, kbeg, iend, jend, kend]

END


;+
;
; NAME:      PLOAD
;
; AUTHOR:    Andrea Mignone
;
; PURPOSE:  
;
;   Provides a convenient and flexible way to read PLUTO data files written
;   in different formats: binary (.dbl/.flt), hdf5 (.dbl.h5/.flt.h5), VTK (.vtk)
;   and  chombo-hdf5 (.hdf5).
;   Data is read and stored in the three common blocks PLUTO_GRID (grid-related
;   information), PLUTO_VAR (variable arrays) and PLUTO_RUN (time-step info).
;   PLOAD must be invoked prior to any other function for properly initializing 
;   the common blocks containing useful information specific to  the data being
;   read. 
;   The common blocks are:
;
;   PLUTO_GRID -> contains grid information: this is read using
;                 the READ_PLUTO_GRID function.
;
;      NX1  =  number of points in the x1 direction
;      NX2  =  number of points in the x2 direction
;      NX3  =  number of points in the x3 direction
;      x1  =  1-D array of NX1 points containing the x1 zone centers
;      x2  =  1-D array of NX2 points containing the x2 zone centers
;      x3  =  1-D array of NX3 points containing the x3 zone centers
;      dx1 =  1-D array of NX1 points containing mesh widths
;      dx2 =  1-D array of NX2 points containing mesh widths
;      dx3 =  1-D array of NX3 points containing mesh widths
;      xpos = a 2D array containing the x-coordinate values (useful for 
;             non-Cartesian geometry);
;      ypos = same as xpos but containing y-coordinate values.
;      AMRLevel = an array of structures, each containing a collection of
;                 rectangular AMR boxes (see HDF5LOAD) 
;      AMRBoxes = an array giving the level number as a function of the 
;                 coordinate indices.
;
;   PLUTO_VAR -> contains the number of variables (=NVAR) and the corresponding
;                file names:
;
;      rho           = density
;      vx1, vx2, vx3 = velocity components
;      Bx1, Bx2, Bx3 = magnetic field components
;      prs           = pressure
;      tr1, ...      = passive scalars
;      Bx1s, Bx2s, Bx3s = staggered components of magnetic field (only with 
;                         double precision binary data files) available with 
;                         constrained transport.
;                      
;   PLUTO_RUN -> contains time-related information:
;
;      t     = time at the corresponding output number
;      dt    = time step at the corresponding output number
;      nlast = the last written file
;
;   File and variable names are read from the corresponding .out file (e.g. 
;   flt.out, dbl.out, dbl.h5.out, etc...)
;
;  
; SYNTAX:
;
;   PLOAD, nout[, /ASSOC,][,DIR=string][, /FLOAT][,/H5][,/HDF5][,/NODATA]
;              [,SHRINK=integer][,/SILENT][, VAR=string][,/VTK][,/XYASSOC]
;              [{X1 | X2 | X3}RANGE=[min,max]][, _EXTRA = extra]
;
; ARGUMENTS:
;
;   nout = the sequential output number, e.g., 0,1,2, ...
;
; KEYWORDS:
;
;   /ASSOC   = use IDL file-association rather than storing files into memory,
;              This can be particularly useful for large data sets.
;
;   DIR      = the directory name where output files are located (e.g. 
;              DIR = 'Data1/'). Default is current directory.
;
;   /FLOAT   = set this keyword if you wish to read single precision
;              data rather than double.
;              The default will try to read dbl (if possible) and
;              will automatically switch to /float if dbl.out is not found.
;
;   /H5      = assume data is in (pixie) HDF5 format used by the static grid
;              version of PLUTO. Double precision data files (.dbl.h5) are 
;              searched for by default. 
;              Single precision files (.flt.h5) can be read by specifying 
;              also the /FLOAT keyword.
;
;   /HDF5    = assume data is in a Chombo HDF5 data file.
;              When this keyword is specified, PLOAD only acts as a wrapper 
;              for the HDF5LOAD  procedure which actually reads the HDF5
;              datafile(s).
;              In this case, different KEYWORDs are supported, see the 
;              documentation for HDF5LOAD.
;
;   /NODATA  = read grid information only. Do not attempt to read/store
;              data (useful for large data sets)
;
;   SHRINK   = an integer number > 1 specifying how much the data has to be 
;              shrinked in all directions (only for .dbl and .flt data formats). 
;              The original array is never stored into memory and resampling is 
;              done slice by slice.
;              Beware that the grid is redefined accordingly.
;              Useful for large 3-D data sets.
;              IMPORTANT: non-uniform mesh are not properly accounted for
;                         and the interpolation may be inaccurate.
;
;   /SILENT  = set this keyword to suppress output generated by pload.
;
;   VAR      = the name of a particular variable to be loaded, e.g., "prs" (only
;              for .dbl or .flt data formats). The remainig ones are ignored.
;
;   /VTK     = assume data is in VTK format (single precision).
;
;   /XYASSOC = associate an array of 2-D data out of a 3-D data file (only for 
;              .dbl and .flt data formats). Each index 'k' of the array 
;              corresponds to an XY plane sliced at different z=z(k).
;
;   X1RANGE   = a 2-element vector in the form [xbeg,xend]
;               giving the x1 axis range to be extracted from the data
;               (only for .dbl and .flt data formats).
;               Only the specified subset is stored into memory and the grid
;               is redefined accordingly. Useful for large datasets.
;
;   X2RANGE   = same as x1range but for the x2-axis.
;
;   X3RANGE   = same as x1range but for the x3-axis.
;
;
;   Please note that '/XYASSOC' and 'SHRINK' are mutually incompatible.
;
; EXAMPLES:
;
;  * Example #1: retrieves the 10th output data set in double precision 
;                from the "My_Data/" directory and store data in the 
;                variables {rho, vx1, vx2, ...}:
;
;    IDL> PLOAD, 10, DIR = 'My_Data/'
;
;
;  * Example #2:  open all files and associate the variables {rho,vx1,...}
;                 to of 2D XY slices at constant z planes.
;                 For instance, rho(k) will be a 2D slice containing density
;                 in the XY plane at z=z(k):
;
;    IDL> PLOAD, /FLOAT, /XYASSOC, 7
;
;
;  * Example #3: open the single-precision datafile containing 'rho'
;                and resize the whole array by making it a factor of 2 smaller.
;                The new grid will be halved respect to the original one.
;                Other variables will be discarded:
;
;    IDL> PLOAD, /FLOAT, 84, VAR="rho", SHRINK=2
;
;
;  * Example #4: extract a 3D portion of the domain from a VTK file using 
;                the range keywords. 
;                The content is stored into {rho,vx1,vx2,...}.
;
;    IDL> PLOAD, 1, /VTK, X1RANGE=[0, 40], X2RANGE=[-10,10], X3RANGE=[-10,10]
;
;
;  * Example #5: load "vx2" from a HDF5 Chombo-data file "data.0004.hdf5"
;                with the uniform resolution corresponding to level 3:
;
;    IDL> PLOAD,/hdf5,4,level=3, var="vx2"
;
;
; LAST MODIFIED:
;
;  May 25, 2018 by A. Mignone (mignone@ph.unito.it)
;
;-

PRO PLOAD, NOUT, ASSOC=ASSOC, DIR=DIR, FLOAT=FLOAT, H5=H5, HDF5=HDF5, $
           NODATA=NODATA, SHRINK=SHRINK, SILENT=SILENT, VAR=VAR, VTK=VTK,$
           XYASSOC=XYASSOC, X1RANGE=X1RANGE, X2RANGE=X2RANGE, X3RANGE=X3RANGE,$
           _EXTRA=extra

  COMMON PLUTO_GRID
  COMMON PLUTO_VAR
  COMMON PLUTO_RUN

; ----------------------------------------------------
;     Set environment on first call
; ----------------------------------------------------

  IF (NOT KEYWORD_SET(first_call)) THEN BEGIN
;   RESOLVE_ROUTINE,"regrid",/COMPILE_FULL_FILE,/either
;   RESOLVE_ALL
    print, "> PLOAD (Aug 2015) - by A. Mignone"
    CASE !D.NAME OF
       'WIN': DEVICE, retain=2, decomposed=0
       'X':   DEVICE, retain=2, decomposed=0, true_color = 24
    ELSE:
    ENDCASE
    first_call = 1
  ENDIF
 
  IF (NOT KEYWORD_SET (SILENT)) THEN silent = 0
  IF (NOT KEYWORD_SET(DIR))     THEN dir = '.'
  IF (NOT SILENT)               THEN PRINT,"> Directory: "+DIR

; -- add a trailing "/" to directory name (if not present) -- 

  IF (STRPOS(dir,'/',strlen(DIR)-1) EQ -1) THEN dir = dir+'/'

; -----------------------------------------------
;  If data is data is from HDF5 (Chombo),
;  go directly to HDF5LOAD and skip the rest
; -----------------------------------------------

  IF (KEYWORD_SET(HDF5)) THEN BEGIN
    filename = dir+"data."+string(nout,format='(I4.4)')+".hdf5"
    ifil = H5F_OPEN(filename)
    igen = H5G_OPEN(ifil,'/Chombo_global')
    dim  = H5A_READ(H5A_OPEN_NAME(igen,'SpaceDim'))
    H5G_CLOSE,igen
    H5F_CLOSE,ifil

    IF ( dim EQ 1) THEN BEGIN
      HDF5LOAD_ONED, nout, dir,  _EXTRA=extra
    ENDIF ELSE BEGIN
      HDF5LOAD, nout, dir, var=var, $
                x1range=x1range,x2range=x2range,x3range=x3range,$
                silent=silent, nodata=nodata,  _EXTRA=extra
    ENDELSE
   
    PRINT_MEMORY_USAGE,SILENT=SILENT
    RETURN
  ENDIF

; ----------------------------------
;  Read pluto grid file
; ----------------------------------

  READ_PLUTO_GRID, DIR, SILENT=SILENT

; -------------------------------------------------
;  Read the corresponding <ext>.out file, where
;  <ext> is the extension. Accepted extensions are
;
;  <ext> = flt, dbl, vtk, flt.h5, dbl.h5.
;
;  If no specific keyword is given, attempt to 
;  read dbl.out (default).
; -------------------------------------------------

  outname = (KEYWORD_SET(FLOAT) ? "flt.out":"dbl.out")
  IF (KEYWORD_SET(VTK)) THEN BEGIN
     outname = "vtk.out"
  ENDIF
  IF (KEYWORD_SET(H5)) THEN BEGIN
    outname = (KEYWORD_SET(FLOAT) ? "flt.h5.out":"dbl.h5.out")
  ENDIF 

  outname = DIR+outname

  fout = FILE_INFO(outname)
  IF (NOT fout.exists) THEN BEGIN
    print,"! cannot find "+outname
    RETURN
  ENDIF

  IF (NOT KEYWORD_SET(SILENT)) THEN PRINT, "> Reading ",outname
  scrh = READ_ASCII(outname,count = nlast)

  timestr = STRARR(nlast)
  t       = FLTARR(nlast)
  dt      = FLTARR(nlast)
  mode    = STRARR(nlast)
  endn    = STRARR(nlast)

  numvar = intarr(nlast)
  names  = strarr(nlast,64)

  nlast = nlast - 1
  IF (nout GT nlast) THEN BEGIN
    print,"! Index out of range in "+outname+$
          " (nlast = ",strcompress(string(nlast)),')'
    RETURN
  ENDIF

  OPENR,1, outname, ERROR = err
  READF, 1, timestr
  CLOSE,1

  FOR n = 0, nlast DO BEGIN
    q = STRSPLIT (timestr(n),' ',/extract)
    t(n)     = q[1]
    dt(n)    = q[2]
    mode(n)  = q[4]
    endn(n)  = q[5]

    sq6 = size(q(6:*))
    numvar(n) = sq6[1]

    names(n, 0:numvar(n)-1) = q(6:*)
  ENDFOR

  NVAR = numvar(nout)

; ---------------------------------------
;  Return if /NODATA has been given
; ---------------------------------------

  IF (KEYWORD_SET(NODATA)) THEN RETURN

; -------------------------------------------------------
;  If data is in H5 format (/H5),  go directly to H5LOAD
; -------------------------------------------------------

  IF (KEYWORD_SET(H5)) THEN BEGIN
    extnum = STRING(format='(I4.4)', nout)
    extnum = STRCOMPRESS(extnum,/remove_all)
    fname = dir+"data."+extnum+ (KEYWORD_SET(FLOAT) ? ".flt.h5":".dbl.h5")
    H5LOAD, fname, silent=silent
    RETURN
  ENDIF

; ---------------------------------------
;      Open files for reading
; ---------------------------------------

  extnum = string(format='(I4.4)', nout)
  extnum = strcompress(extnum,/remove_all)

  SINGLE_FILE = 0
  IF (mode(nout) EQ "single_file") THEN SINGLE_FILE = 1
  MULTIPLE_FILES = 1-SINGLE_FILE

; -----------------------------------------------------------
;  If data is in VTK format (/VTK), go directly to VTK_LOAD
; -----------------------------------------------------------

  IF (KEYWORD_SET(VTK)) THEN BEGIN
    extnum = STRING(format='(I4.4)', nout)
    extnum = STRCOMPRESS(extnum,/remove_all)
    ext    = "."+extnum+".vtk"
    IF (SINGLE_FILE) THEN BEGIN
      fname = dir+"data"+ext
      VTK_LOAD, fname, SILENT=SILENT, VAR=VAR, $
                X1RANGE=X1RANGE, X2RANGE=X2RANGE, X3RANGE=X3RANGE, domain_indx, _EXTRA=extra
    ENDIF
    IF (MULTIPLE_FILES) THEN BEGIN; with multiple files, one has to separately
                                  ; handle vector fields, since data file
                                  ; do not begin with the variable name.

     ; !! Important: numvar does not contain "vx2" and "vx3" and we add 
     ;    these variable names separately. 

      IF (NOT KEYWORD_SET(VAR)) THEN BEGIN
        FOR nv = 0, numvar[nout] - 1 DO BEGIN
          fname = dir+names[nout,nv]+ext
          IF (names[nout,nv] EQ "vx1") THEN fname = dir+"vfield"+ext
          IF (names[nout,nv] EQ "bx1") THEN fname = dir+"bfield"+ext
          VTK_LOAD, fname, SILENT=SILENT, $
                    X1RANGE=X1RANGE, X2RANGE=X2RANGE, X3RANGE=X3RANGE, domain_indx, _EXTRA=extra
        ENDFOR
      ENDIF ELSE BEGIN
        fname = dir+VAR+ext
        IF (VAR EQ "vx1" OR VAR EQ "vx2" OR VAR EQ "vx3") THEN fname = dir+"vfield"+ext
        IF (VAR EQ "bx1" OR VAR EQ "bx2" OR VAR EQ "bx3") THEN fname = dir+"bfield"+ext
        VTK_LOAD, fname, SILENT=SILENT, VAR=VAR, $
                  X1RANGE=X1RANGE, X2RANGE=X2RANGE, X3RANGE=X3RANGE, domain_indx, _EXTRA=extra
      ENDELSE
    ENDIF

    IF (KEYWORD_SET(X1RANGE) OR KEYWORD_SET(X2RANGE) OR KEYWORD_SET(X3RANGE)) THEN BEGIN
      ibeg = domain_indx(0) & iend = domain_indx(3)
      jbeg = domain_indx(1) & jend = domain_indx(4)
      kbeg = domain_indx(2) & kend = domain_indx(5)
      NX1 = iend - ibeg + 1
      NX2 = jend - jbeg + 1
      NX3 = kend - kbeg + 1
      x1 = x1(ibeg:iend) & dx1 = dx1(ibeg:iend)
      x2 = x2(jbeg:jend) & dx2 = dx2(jbeg:jend)
      x3 = x3(kbeg:kend) & dx3 = dx3(kbeg:kend)
    ENDIF
   
    PRINT_MEMORY_USAGE,SILENT=SILENT
    RETURN
  ENDIF

; --------------------------------------
;  Get precision and number of bytes
; --------------------------------------

  REAL   = 5; -- default is double --
  NBYTES = 8
  ext = "."+extnum+".dbl"

  IF (KEYWORD_SET(FLOAT)) THEN BEGIN
    REAL   = 4
    NBYTES = 4
    ext = "."+extnum+".flt"
  ENDIF

  IF (NOT KEYWORD_SET(SILENT)) THEN BEGIN
    IF (NBYTES EQ 4) THEN BEGIN
      print, "> PLOAD: using single precision binary data"
    ENDIF

    IF (NBYTES EQ 8) THEN BEGIN
      print, "> PLOAD: using double precision binary data"
    ENDIF
  ENDIF

; --------------------------------------------
;   Open .dbl or .flt using SINGE_FILE mode
; --------------------------------------------

  IF (SINGLE_FILE) THEN BEGIN
    fname = DIR+"data"+ext
    IF (endn[nout] EQ "little") THEN BEGIN
      OPENR, unit, fname, /GET_LUN, /SWAP_IF_BIG_ENDIAN, _EXTRA=extra, ERROR = err
    ENDIF ELSE BEGIN
      OPENR, unit, fname, /GET_LUN, /SWAP_IF_LITTLE_ENDIAN, _EXTRA=extra, ERROR = err
    ENDELSE
    IF (err NE 0) THEN BEGIN
      print,"! Cannot open file "+fname
      CLOSE,unit
      FREE_LUN,unit
      RETURN
    ENDIF
  ENDIF

; ----------------------------------------------------------
;  Begin main loop on variables:
;  for each variable:
;
;  * open the corresponding file (in MULTIPLE_FILES mode)
;  * associate a pointer to the the file
;  * store the arrays into memory (only if /ASSOC has not been given)
;  * copy the pointer content into a variable
;  * increment offset
; ----------------------------------------------------------

  offset = ulong64(0)
  FOR nv = 0, numvar[nout] - 1 DO BEGIN

  ; --- check if the VAR keyword has been given ---

    IF (KEYWORD_SET (VAR)) THEN BEGIN
      IF (names[nout,nv] NE var) THEN BEGIN
    ; -- skip to end of the loop for single file, continue otherwse --
        IF (SINGLE_FILE) THEN GOTO,OFFSET ELSE CONTINUE
      ENDIF
    ENDIF  

  ; ----------------------------------------------
  ;   Open .dbl or .flt using MULTIPLE_FILES mode
  ; ----------------------------------------------

    IF (MULTIPLE_FILES) THEN BEGIN
      unit = 48 + nv; ** Units begin at 48 **
      fname = dir+names(nout,nv)+ext
      IF (endn[nout] EQ "little") THEN BEGIN
        OPENR, unit, fname, /SWAP_IF_BIG_ENDIAN, _EXTRA=extra, ERROR = err
      ENDIF ELSE BEGIN
        OPENR, unit, fname, /SWAP_IF_LITTLE_ENDIAN, _EXTRA=extra, ERROR = err
      ENDELSE
      IF (err NE 0) THEN BEGIN
        print,"! Cannot open file "+fname
        CLOSE,unit
        FREE_LUN,unit
        CONTINUE
      ENDIF
    ENDIF

  ; --------------------------------------------
  ;  start by associating a pointer to the file
  ; --------------------------------------------

  ; **** set extra zones for staggered fields ****

    n1p = 0 & n2p = 0 & n3p = 0
    IF (names(nout,nv) EQ "Bx1s") THEN n1p = 1
    IF (names(nout,nv) EQ "Bx2s") THEN n2p = 1
    IF (names(nout,nv) EQ "Bx3s") THEN n3p = 1

  ; **** check if either one of shrink/xyassoc/xranges have been given ****

    CASE 1 OF

      (KEYWORD_SET (SHRINK)): BEGIN; SHRINK
        vpt = ASSOC(unit, MAKE_ARRAY(NX1+n1p,NX2+n2p,shrink, type = REAL),offset)
        REBIN3D, vpt,shrink
      END; SHRINK

      (KEYWORD_SET (XYASSOC)): BEGIN; xyassoc
        vpt = ASSOC(unit, MAKE_ARRAY(NX1+n1p,NX2+n2p, type = REAL),offset)
      END; xyassoc

      (KEYWORD_SET (X1RANGE) OR KEYWORD_SET (X2RANGE) OR KEYWORD_SET (X3RANGE)): BEGIN;

        IF (NOT KEYWORD_SET (X1RANGE)) THEN x1range=[min(x1),max(x1)]
        IF (NOT KEYWORD_SET (X2RANGE)) THEN x2range=[min(x2),max(x2)]
        IF (NOT KEYWORD_SET (X3RANGE)) THEN x3range=[min(x3),max(x3)]
        vpt = ASSOC(unit, MAKE_ARRAY(NX1+n1p,NX2+n2p, type = REAL),offset)
        ARRAY_EXTRACT, vpt, x1range, x2range, x3range,domain_indx,silent=silent
      END; DOMAIN

      ELSE: BEGIN
        vpt = ASSOC(unit, make_array(NX1+n1p,NX2+n2p,NX3+n3p, type = REAL),offset)

        ; ** store arrays into memory if the   **
        ; ** /assoc keyword has not been given **

        IF (NOT KEYWORD_SET(ASSOC)) THEN vpt = vpt(0)

      ENDELSE

    ENDCASE

  ; ---------------------------------------------
  ;    label pointer data with suitable names
  ; ---------------------------------------------

    MATCH_VARNAME,vpt, names(nout,nv), silent=silent

    vpt = 0; free memory

    OFFSET:
    IF (SINGLE_FILE) THEN BEGIN
      n1_off = ulong64(NX1)
      n2_off = ulong64(NX2)
      n3_off = ulong64(NX3)

      IF (names(nout,nv) EQ "Bx1s") THEN n1_off = ulong64(NX1 + 1)
      IF (names(nout,nv) EQ "Bx2s") THEN n2_off = ulong64(NX2 + 1)
      IF (names(nout,nv) EQ "Bx3s") THEN n3_off = ulong64(NX3 + 1)

      offset = offset + ulong64(n1_off*n2_off*n3_off*NBYTES)
      POINT_LUN,unit,offset; move file pointer 
    ENDIF

    IF (MULTIPLE_FILES) THEN BEGIN
      CLOSE,unit
      FREE_LUN,unit
    ENDIF
  ENDFOR ; nv = 0, NVAR-1
  
  IF (SINGLE_FILE) THEN BEGIN
    CLOSE,unit
    FREE_LUN,unit
  ENDIF

; **** print memory usage ****

  PRINT_MEMORY_USAGE,SILENT=SILENT

; **** change the grid if one of SHRINK/XRANGE's keywords has been given ****

  IF (KEYWORD_SET(SHRINK)) THEN BEGIN
    NX1 = NX1/shrink
    NX2 = NX2/shrink
    NX3 = NX3/shrink
    x1 = CONGRID(x1,NX1,/CENTER)
    x2 = CONGRID(x2,NX2,/CENTER)
    x3 = CONGRID(x3,NX3,/CENTER)

    dx1 = CONGRID(dx1,NX1,/CENTER)
    dx2 = CONGRID(dx2,NX2,/CENTER)
    dx3 = CONGRID(dx3,NX3,/CENTER)
  ENDIF

  IF (KEYWORD_SET(X1RANGE) OR KEYWORD_SET(X2RANGE) OR KEYWORD_SET(X3RANGE)) THEN BEGIN
    ibeg = domain_indx(0) & iend = domain_indx(3)
    jbeg = domain_indx(1) & jend = domain_indx(4)
    kbeg = domain_indx(2) & kend = domain_indx(5)
    NX1 = iend - ibeg + 1
    NX2 = jend - jbeg + 1
    NX3 = kend - kbeg + 1
    x1 = x1(ibeg:iend) & dx1 = dx1(ibeg:iend)
    x2 = x2(jbeg:jend) & dx2 = dx2(jbeg:jend)
    x3 = x3(kbeg:kend) & dx3 = dx3(kbeg:kend)
  ENDIF

END

