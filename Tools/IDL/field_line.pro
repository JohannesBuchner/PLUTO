;+
; NAME:     FIELD_INTERP
; 
; AUTHOR:   A. Mignone
; 
; PURPOSE:  Linearly interpolate the vector (Vx, Vy, Vz) defined on the 
;           grid (x,y,z) at the point with coordinate (xp, yp, zp) 
;
; SYNTAX    FIELD_INTERP, Vx, Vy, x, y, z, coord
; 
; ARGUMENTS:
;
;   Vx,Vy, Vz:    Vector components to be interpolated;
;   x,y,z:        1-D coordinate arrays; 
;   coord:        a three-element vector with the coordinates of 
;                 the point where interpolation is required.
; 
; LAST MODIFIED:  June 19, 2013 by A. Mignone (mignone@ph.unito.it)
;-
FUNCTION FIELD_INTERP, Vx, Vy, Vz, x, y, z, coord

; ----------------------------------
;  check the number of dimensions
;  and sizes
; ----------------------------------

 ndim = SIZE(Vx, /N_DIM)

; -----------------------------------------
; q is the two element vector with the
; interpolated values of U and V
; at xp and yp
; -----------------------------------------

 q  = fltarr(ndim)
 xp = coord[0]
 yp = coord[1]
 IF (ndim EQ 3) THEN  zp = coord[2] ELSE zp = 0.0 

; ----------------------------------------
;  Locate i0 and j0, indexes of
;  cell to which xp,yp belong to
; ----------------------------------------

 i0 = (WHERE(xp LT x))[0]
 j0 = (WHERE(yp LT y))[0]
 IF (ndim EQ 3) THEN k0 = (WHERE(zp LT z))[0]  ELSE k0 = 0

 IF i0 EQ -1 THEN i0 = 0
 IF j0 EQ -1 THEN j0 = 0
 IF k0 EQ -1 THEN k0 = 0

; --------------------------------------------
;   Interpolate
; --------------------------------------------

 IF (ndim EQ 2) THEN BEGIN

 ; -- interpolate Vx --

   Vx_x = INTERPOL(Vx(*,j0),x,xp)
   Vx_y = INTERPOL(Vx(i0,*),y,yp)
   q(0) = Vx_x + Vx_y - Vx(i0,j0)

 ; -- interpolate Vy --

   Vy_x = INTERPOL(Vy(*,j0),x,xp)
   Vy_y = INTERPOL(Vy(i0,*),y,yp)

   q(1) = Vy_x + Vy_y - Vy(i0,j0)

 ENDIF

 IF (ndim EQ 3) THEN BEGIN

 ; -- interpolate Vx --

   Vx_x = interpol(Vx(*,j0,k0),x,xp)
   Vx_y = interpol(Vx(i0,*,k0),y,yp)
   Vx_z = interpol(Vx(i0,j0,*),z,zp)
   q(0) = Vx_x + Vx_y + Vx_z - 2.0*Vx(i0,j0,k0)

 ; -- interpolate Vy --

   Vy_x = interpol(Vy(*,j0,k0),x,xp)
   Vy_y = interpol(Vy(i0,*,k0),y,yp)
   Vy_z = interpol(Vy(i0,j0,*),z,zp)
   q(1) = Vy_x + Vy_y + Vy_z - 2.0*Vy(i0,j0,k0)

 ; -- interpolate Vz --

   Vz_x = interpol(Vz(*,j0,k0),x,xp)
   Vz_y = interpol(Vz(i0,*,k0),y,yp)
   Vz_z = interpol(Vz(i0,j0,*),z,zp)
   q(2) = Vz_x + Vz_y + Vz_z - 2.0*Vz(i0,j0,k0)
 ENDIF

 RETURN, q
END


;+
;
; NAME:    FIELD_LINE
;
; AUTHOR:  A. Mignone (mignone@ph.unito.it)
;
; PURPOSE: Given a 2 or 3D vector field (Vx, Vy) or (Vx, Vy, Vz) computes the 
;          field line passing through the point (seed) [xseed, yseed, zseed].
;          The line is computed by solving a system of ODE of the form
;          
;            dx/dt = Vx(x,y,z)
;            dy/dt = Vy(x,y,z)
;            dz/dt = Vz(x,y,z)
;           
;          Integration stops when either the domain boundaries are reached or 
;          the max number of iteration is exceeded.
;
; SYNTAX: FIELD_LINE, Vx, Vy, Vz, x, y, z, seed = seed, pnt_list[,step=step]
;                     [,method=method][,maxstep=maxstep][,minstep=minstep][,tol=tol]
;
;
; ARGUMENTS:
;
;   Vx,Vy,Vz: 3D arrays giving the three vector components. In 2D, both Vz
;             and z must be scalars and equal to 0.
;
;   x,y,z:    1D coordinate arrays on top of which Vx, Vy and Vz are defined.
;             In 2D, set z to be 0.0
;
;   seed:     a 3-element array giving the point coordinates through which the
;             field line goes. 
;
;   pnt_list: on output, in contains 2D array giving the field line coordinates
;             {x,y,z} = {pnt_list[0,*], pnt_list[1,*], pnt_list[2,*]} (in 3D) or
;             {x,y }  = {pnt_list[0,*], pnt_list[1,*]} (in 3D) 
;             
;
; KEYWORDS:
;
;   step:   a scalar giving the initial step size. Default is (mean) grid spacing.
;
;   method: a string giving the integration method. The possible choices are:
;
;            "RK2"   explicit, fixed step, 2nd order Runge-Kutta methods.
;            "RK4"   explicit, fixed step, 4th order Runge-Kutta methods.
;            "BS23"  explicit, adaptive stepsize Runge-Kutta-Fehlberg of order 
;                    3 with local truncation error based on a 2nd-order accurate
;                    embedded solution.
;            "CK45"  explicit, adaptive stepsize Cask-Karp of order 
;                    5 with local truncation error based on a 4th-order accurate
;                    embedded solution.
;
;           The default is "RK2". Use an adaptive stepsize integrator
;           when the field line does not close on itself inside the domain.
; 
;
;   maxstep: a scalar giving the maximum allowed integration step.
;            Default is 100*step.
;
;   minstep: a scalar giving the minimum allowed integration step. 
;            Default is 0.05*step.
;
;   tol:   a scalar value giving the relative tolerance during adaptive stepsize
;          integration. It is ignored for fixed step-size integration (such as RK2, RK4)
;         
;     
;
; EXAMPLE:
;
;   * compute a field line tangent to the vector field (Bx1,Bx2) in 2D at the 
;     point with coordinate (-1,2) using the Bogacki-Shampine method with relative
;     tolerance 1.e-4:
;
;   IDL> field_line, Bx1, Bx2, 0.0 x1, x2, 0.0, seed=[-1,2], pl, method="BS23", tol=1.e-4
;   IDL> oplot, pl[0,*], pl[1,*]  ; overplot on current window
;
;   * Same as before but in 3D and at the point [-1,2,0.5]:
;
;   IDL> field_line, Bx1, Bx2, Bx3, x1, x2, x3, seed=[-1,2,0.5], pl, method="BS23", tol=1.e-4
;    
;
; LAST MODIFIED:    June 19, 2013 by A. Mignone (mignone@ph.unito.it)
;
;-
PRO FIELD_LINE, Vx, Vy, Vz, x, y, z, seed=seed, pnt_list, step=step,$
                method=method,maxstep=maxstep,minstep=minstep,tol=tol

 ON_ERROR,2

 sz = size(Vx)

; --------------------------------------------------
;  Check if we're running 2D or 3D
; --------------------------------------------------

 nx = sz[1]
 ny = sz[2]
 ndim = 2
 nz   = 0
 IF (sz[0] EQ 3) THEN BEGIN
   ndim = 3
   nz = sz[3]
 ENDIF

 npt  = INTARR(3)
 dbeg = FLTARR(3)
 dend = FLTARR(3)
 L    = FLTARR(3)

; ------------------------------------------
;  Get domain sizes. 
;  Take the initial and final coordinates
;  slightly larger to allow a seed to be 
;  specified on the boundary. 
; ------------------------------------------

 dbeg[0] = x(0)    - 0.51*(x(1) - x(0))
 dend[0] = x(nx-1) + 0.51*(x(nx-1) - x(nx-2))
 L[0]    = dend[0] - dbeg[0]
 
 dbeg[1] = y(0)    - 0.51*(y(1) - y(0))
 dend[1] = y(ny-1) + 0.51*(y(ny-1) - y(ny-2))
 L[1]    = dend[1] - dbeg[1]

 IF (ndim EQ 3) THEN BEGIN
   dbeg[2] = z(0)    - 0.51*(z(1) - z(0))
   dend[2] = z(nz-1) + 0.51*(z(nz-1) - z(nz-2))
   L[2]    = dend[2] - dbeg[2]
 ENDIF
 
; -------------------------------------
;  Normalize vector field to 1.
;  Only direction can change.
; -------------------------------------

 Vx = REFORM(Vx) 
 Vy = REFORM(Vy)
 Vz = REFORM(Vz) 
 
 IF (ndim EQ 2) THEN BEGIN
   norm = 1.0/sqrt(Vx*Vx + Vy*Vy + 1.e-18)
   Vx = Vx*norm
   Vy = Vy*norm
 ENDIF

 IF (ndim EQ 3) THEN BEGIN
   norm = 1.0/sqrt(Vx*Vx + Vy*Vy + Vz*Vz + 1.e-18)
   Vx = Vx*norm
   Vy = Vy*norm
   Vz = Vz*norm
 ENDIF

; ------------------------------------------
;  Make sure initial seed point falls 
;  inside the computational domain.
; ------------------------------------------

 inside_domain = 1
 FOR nd = 0, ndim - 1 DO BEGIN
   in = (seed[nd] LT dend[nd]) AND (seed[nd] GT dbeg[nd])
   inside_domain = inside_domain AND in
 ENDFOR

 IF (inside_domain EQ 0) THEN BEGIN
   print,"! Seed Point ",seed[0],seed[1]," falls outside grid range"
   RETURN
 ENDIF

; ----------------------------------------------------------
;  xfwd = coordinates for forward  integration
;  xbck = coordinates for backward integration
; ----------------------------------------------------------

 max_steps = 16384
 max_fail  = 1024

 xfwd = FLTARR(ndim, max_steps)
 xbck = FLTARR(ndim, max_steps)
 xk0  = FLTARR(ndim) & xk1 = FLTARR(ndim)
 xk2  = FLTARR(ndim) & xk3 = FLTARR(ndim)
 xk4  = FLTARR(ndim) & xk5 = FLTARR(ndim)

; ----------------------------------
;  Set initial conditions
; ----------------------------------

 FOR nd = 0, ndim-1 DO xfwd[nd,0] = seed[nd]
 FOR nd = 0, ndim-1 DO xbck[nd,0] = seed[nd]

; ---------------------------------------------
;  Check keywords: step, method and tolerance
; ---------------------------------------------

 IF (NOT KEYWORD_SET(method))  THEN method = "RK2"
 IF (NOT KEYWORD_SET(step))    THEN BEGIN
   step = MIN((dend[0:ndim-1] - dbeg[0:ndim-1])/sz[1:ndim])
 ENDIF
 IF (NOT KEYWORD_SET(tol))     THEN tol = 1.e-6
 IF (NOT KEYWORD_SET(maxstep)) THEN maxstep = 100*step
 IF (NOT KEYWORD_SET(minstep)) THEN minstep = 0.05*step

; -------------------------------------
;  tolerance factor should scale down 
;  actual to domain size.
; -------------------------------------

 tol = tol*MAX(L)

; --------------------------------------------------------------------
;  Set the coefficients for adaptive step size integration:
;  Cash-Karp 45 (CK45) and Bogacki-Shampine 23 (BS23).
;  Taken from:
;
;  http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
; --------------------------------------------------------------------

 IF (method EQ "CK45") THEN BEGIN
   b1 = 37.0/378.0  & b2 = 0.0 & b3 = 250.0/621.0
   b4 = 125.0/594.0 & b5 = 0.0 & b6 = 512.0/1771.0

   bs1 = 2825.0/27648.0  & bs2 = 0.0           & bs3 = 18575.0/48384.0
   bs4 = 13525.0/55296.0 & bs5 = 277.0/14336.0 & bs6 = 0.25

   a21 = 0.2
   a31 = 3.0/40.0       & a32 = 9.0/40.0
   a41 = 0.3            & a42 = -0.9        & a43 = 1.2
   a51 = -11.0/54.0     & a52 = 2.5         & a53 = -70.0/27.0    & a54 = 35.0/27.0
   a61 = 1631.0/55296.0 & a62 = 175.0/512.0 & a63 = 575.0/13824.0 & a64 = 44275.0/110592.0 & a65 = 253.0/4096.0
 ENDIF

 IF (method EQ "BS23") THEN BEGIN
   b1  = 2.0/9.0  & b2  = 1.0/3.0 & b3  = 4.0/9.0 & b4  = 0.0
   bs1 = 7.0/24.0 & bs2 = 1.0/4.0 & bs3 = 1.0/3.0 & bs4 = 1.0/8.0

   a21 = 0.5
   a31 = 0.0     & a32 = 0.75
   a41 = 2.0/9.0 & a42 = 1.0/3.0 & a43 = 4.0/9.0
 ENDIF

; ------------------------------------------------
;    Integrate Backward (s=-1) and Forward (s=1)
; ------------------------------------------------

 FOR s = -1,1,2 DO BEGIN
   dh = s*step
   inside_domain = 1
   k             = 0
   kfail         = 0
   WHILE (inside_domain AND (k LT MAX_STEPS-1)) DO BEGIN  ; attempt to integrate from k 
                                                          ; to k+1.

     ; --------------------------------
     ;  restrict dh to lie between 
     ;  minstep and maxstep
     ; --------------------------------

     dh = s*min([abs(dh),maxstep])
     dh = s*max([abs(dh),minstep])

;     IF (abs(dh)/minstep LE 1.0) THEN print,"Minimum step reached"
;     IF (abs(dh)/maxstep GE 1.0) THEN print,"Maximum step reached"

     ; -- set initial condition 

     xk0 = xfwd(*,k)

     ; ----------------------------------------------------------
     ;   Explicit Runge-Kutta method with 2nd order accuracy.
     ;   Fixed step size. Requires 2 function evaluations.
     ; ----------------------------------------------------------
     IF (method EQ "RK2") THEN BEGIN  
       k = k+1

       k1  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk0)
       xk1 = xk0 + 0.5*dh*k1

       k2 = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xfwd[*,k] = xk0 + dh*k2
     ENDIF

     ; ----------------------------------------------------------
     ;   Explicit Runge-Kutta method with 4th order accuracy.
     ;   Fixed step size. Requires 4 function evaluations.
     ; ----------------------------------------------------------
     IF (method EQ "RK4") THEN BEGIN  
       k = k+1

       k1  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk0)
       xk1 = xk0 + 0.5*dh*k1

       k2  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk1 = xk0 + 0.5*dh*k2

       k3 = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk1 = xk0 + dh*k3

       k4 = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xfwd[*,k] = xk0 + dh*(k1 + 2.0*(k2 + k3) + k4)/6.0
     ENDIF

     ; ---------------------------------------------------------------
     ;  Explicit Runge-Kutta-Fehlberg pair (2,3) with adaptive 
     ;  step size. It is also known as Bogacki-Shampine and provide
     ;  third-order accuracy using a total of 4 function evaluations.
     ; ---------------------------------------------------------------
     IF (method EQ "BS23") THEN BEGIN ; -- use BS23

       k1  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk0)
       xk1 = xk0 + dh*a21*k1

       k2  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk1 = xk0 + dh*(a31*k1 + a32*k2)

       k3  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk1 = xk0 + dh*(a41*k1 + a42*k2 + a43*k3)

       k4  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk3 = xk0 + dh*(b1*k1 + b2*k2 + b3*k3 + b4*k4)

       xk2 = xk0 + dh*(bs1*k1 + bs2*k2 + bs3*k3 + bs4*k4)

       ; ---- compute error ----

       err = MAX(ABS(xk3 - xk2))/tol
       IF ((err LT 1.0) OR (ABS(dh)/minstep LE 1.0)) THEN BEGIN  ; -- accept step
         k      = k + 1
         err    = MAX([err,1.e-12])
         dhnext = 0.9*abs(dh)*err^(-0.3333)
         dhnext = MIN([dhnext,3.0*ABS(dh)])
         dh     = s*dhnext
         xfwd[*,k] = xk3
       ENDIF ELSE BEGIN

         dh = 0.9*s*ABS(dh)*err^(-0.5)
         IF (kfail GT MAX_FAIL) THEN BEGIN
           print,"! Too many failures"
           RETURN
         ENDIF
         CONTINUE
        ENDELSE
      ENDIF ; method BS23

     ; ---------------------------------------------------------------
     ;  Cash-Karp fifth-order method using a (4,5) pair.
     ;  Provides adaptive step-size control with monitoring of local 
     ;  truncation error. It requires 6 function evaluations.
     ; ---------------------------------------------------------------
     IF (method EQ "CK45") THEN BEGIN 

       k1  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk0)
       xk1 = xk0 + dh*a21*k1

       k2  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk1 = xk0 + dh*(a31*k1 + a32*k2)

       k3  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk1 = xk0 + dh*(a41*k1 + a42*k2 + a43*k3)

       k4  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk1 = xk0 + dh*(a51*k1 + a52*k2 + a53*k3 + a54*k4)

       k5  = FIELD_INTERP(Vx, Vy, Vz, x, y, z,  xk1)
       xk1 = xk0 + dh*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5)

       k6  = FIELD_INTERP(Vx, Vy, Vz, x, y, z, xk1)
       xk5 = xk0 + dh*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6)

       xk4 = xk0 + dh*(bs1*k1 + bs2*k2 + bs3*k3 + bs4*k4 + bs5*k5 + bs6*k6)

       ; ---- compute error ----

       err = MAX(ABS(xk5 - xk4))/tol
       IF ((err LT 1.0) OR (ABS(dh)/minstep LE 1.0)) THEN BEGIN  ; -- accept step
         k      = k + 1
         err    = MAX([err,1.e-12])
         dhnext = 0.9*ABS(dh)*err^(-0.2)
         dhnext = MIN([dhnext,3.0*ABS(dh)])
         dh     = s*dhnext
         xfwd[*,k] = xk5
       ENDIF ELSE BEGIN

         dh = 0.9*s*ABS(dh)*err^(-0.25)
         IF (kfail GT MAX_FAIL) THEN BEGIN
           print,"! Too many failueres"
           RETURN
         ENDIF
         CONTINUE
        ENDELSE
      ENDIF ; method CK45

     ; ----------------------------------------------------
     ;  Check whether we're still inside the domain...
     ; ----------------------------------------------------

     inside_domain = 1
     FOR nd = 0, ndim - 1 DO BEGIN
       in = (xfwd[nd,k] GT dbeg[nd]) AND (xfwd[nd,k] LT dend[nd])
       inside_domain = inside_domain AND in
     ENDFOR

   ENDWHILE
   IF (S EQ -1) THEN BEGIN
     xbck  = xfwd
     k_bck = k
   ENDIF
 ENDFOR
 k_fwd = k
 if (k_fwd GE (MAX_STEPS-2)) THEN PRINT,"! Max number of iteration exceeded"
 if (k_bck GE (MAX_STEPS-2)) THEN PRINT,"! Max number of iteration exceeded"

 print,"Method: ",method,$      ; "; tol = ",strcompress(string(tol,format='(e10.2)')),$
       "; Forward steps: ",strcompress(string(k_fwd,format='(i5)')),$
       "; Bckward steps: ",strcompress(string(k_bck,format='(i5)'))

; --------------------------------------------
;         return arrays
; --------------------------------------------

  ;xfield = [reverse(REFORM(xbck(0,0:k_bck))),REFORM(xfwd(0,0:k_fwd))]
  ;yfield = [reverse(REFORM(xbck(1,0:k_bck))),REFORM(xfwd(1,0:k_fwd))]
  ;zfield = [reverse(REFORM(xbck(2,0:k_bck))),REFORM(xfwd(2,0:k_fwd))]

  npt = k_bck + k_fwd + 2
  
  pnt_list = FLTARR(ndim, npt)
  FOR nd = 0, ndim-1 DO BEGIN
    pnt_list[nd,*] = [reverse(REFORM(xbck(nd,0:k_bck))),REFORM(xfwd(nd,0:k_fwd))]
  ENDFOR
     
END

