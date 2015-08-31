;+
FUNCTION FIELD_INTERP, U, V, W, x, y, z, coord

;
;
;
;
; PURPOSE:  Linearly interpolate the vector
;           fields  U and V at point xp and yp
;
;
;
; ARGUMENTS:
;
;   U,V:    2-D vectors to be interpolated
;   nx,ny:  number of points in the two directions
;   x,y:    1-D abscissas and ordinates
;   xp,yp:  location where U and V should be interpolated
;
;
;  EXAMPLES
;
;    #1 Plot the logarithm of the density at the initial time and overplot
;       2 magnetic fieldlines rooted at the (x,y) = (1.0,0.1) and (2.0,0.1)
;    IDL> display, alog(rho(0)), x1 = x1, x2 = x2, /vbar
;    IDL> field_line, b1(0), b2(0), x1, x2, 1.0, 0.1, qx, qy
;    IDL> oplot, qx, qy
;    IDL> field_line, b1(0), b2(0), x1, x2, 2.0, 0.1, qx, qy
;    IDL> oplot, qx, qy
;
;
;-

; ----------------------------------
;  check the number of dimensions
;  and sizes
; ----------------------------------

 DIM = 3
 IF (size(w,/dim) EQ 0) THEN DIM = 2

 nx = SIZE(x,/DIM)
 ny = SIZE(y,/DIM)
 nz = SIZE(z,/DIM)

; -----------------------------------------
; q is the two element vector with the
; interpolated values of U and V
; at xp and yp
; -----------------------------------------

 q  = fltarr(DIM)
 xp = coord[0] & yp = coord[1]
 IF (DIM EQ 3) THEN zp = coord[2] ELSE zp = 0.0

; ----------------------------------------
;  Locate i0 and j0, indexes of
;  cell to which xp,yp belong to
; ----------------------------------------

 i0 = (WHERE(xp LT x))[0]
 j0 = (WHERE(yp LT y))[0]
 IF (DIM EQ 3) THEN k0 = (WHERE(zp LT z))[0]  ELSE k0 = 0

 IF i0 EQ -1 THEN i0 = 0
 IF j0 EQ -1 THEN j0 = 0
 IF k0 EQ -1 THEN k0 = 0

; interpolate U along x
;

 scrh1 = interpol(U(*,j0),x,xp)
 scrh2 = interpol(U(i0,*),y,yp)
 q(0)  = scrh1 + scrh2 - U(i0,j0)

;print,"q0 :",q(0)
;
; interpolate V along x
;

 scrh1 = interpol(V(*,j0),x,xp)
 scrh2 = interpol(V(i0,*),y,yp)
 q(1)  = scrh1 + scrh2 - V(i0,j0)
;print,"q1 :",q(1)

 RETURN, q
END


;+
;
; NAME:    FIELD_LINE
;
; AUTHOR:  A. Mignone (mignone@ph.unito.it)
;
; PURPOSE: Given a 2-D vector field (U,V), it computes the field line
;          passing through x0, y0. Integration stops when either 
;          the domain boundaries are reached or the max number of iteration
;          is exceeded.
;
; SYNTAX: FIELD_LINE, U, V, x, y, x0, y0, xfield, yfield[,step=step]
;                     [,method=method][,maxstep=maxstep][,minstep=minstep][,tol=tol]
;
;
; ARGUMENTS:
;
;   U,V: 2-D arrays giving the 1st (U) and 2nd (V) component of the 
;         vector field
;
;   x,y: 1-D arrays giving abscissas and ordinates
;
;   x0,y0: two scalars giving the point coordinates through which the
;          field line goes.
;
;   xfield,yfield: 1-D arrays containing the output abscissas and ordinates
;                  of the field line.
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
;   * compute a field line tangent to the vector field (B1,B2) in the point
;     with coordinate (-1,2) using the Bogacki-Shampine method with relative
;     tolerance 1.e-4:
;
;   IDL> field_line, b1, b2, x1, x2, -1, 2, xl, yl, method="BS23", tol=1.e-4
;   IDL> oplot,xl,yl  ; overplot on current window
;
;
; LAST MODIFIED
;
;   September 23, 2010 by A. Mignone (mignone@ph.unito.it)
;
;
;-
PRO FIELD_LINE, U0, V0, x, y, x0, y0, xfield, yfield,step=step,$
                method=method,maxstep=maxstep,minstep=minstep,tol=tol

; ---------------------------------
;    get domain sizes
; ---------------------------------

 sa = size(U0)
 nx = sa[1]
 ny = sa[2]

 xbeg = x(0)    - 0.51*(x(1) - x(0))
 xend = x(nx-1) + 0.51*(x(nx-1) - x(nx-2))

 ybeg = y(0)    - 0.51*(y(1) - y(0))
 yend = y(ny-1) + 0.51*(y(ny-1) - y(ny-2))

 Lx = xend - xbeg
 Ly = yend - ybeg

; ---------------------------------
;  Normalize vector field to 1.
;  Only direction can change.
; ---------------------------------

 U = REFORM(U0) & V = REFORM(V0)
 NORM = 1.0/sqrt(U*U + V*V + 1.e-18)
 U = U*NORM & V = V*NORM

; -------------------------------
;  check whether the point x0,y0
;  falls  outside grid ranges
; -------------------------------

 inside_domain = (x0 LT xend) AND (x0 GT xbeg) AND $
                 (y0 LT yend) AND (y0 GT ybeg)

 IF (inside_domain EQ 0) THEN BEGIN
   print,"! Point ",x0,y0," falls outside grid range"
   print,xbeg, xend
   print,ybeg, yend
   RETURN
 ENDIF

; ----------------------------------------------------------
;  xfwd = coordinates for forward  integration
;  xbck = coordinates for backward integration
; ----------------------------------------------------------

 DIM       = 2
 MAX_STEPS = 16384
 MAX_FAIL  = 1024

 xfwd = fltarr(DIM,MAX_STEPS)
 xbck = fltarr(DIM,MAX_STEPS)
 xk0  = fltarr(DIM) & xk1 = fltarr(DIM)
 xk2  = fltarr(DIM) & xk3 = fltarr(DIM)
 xk4  = fltarr(DIM) & xk5 = fltarr(DIM)

; ----------------------------------
;       Initial conditions
; ----------------------------------

 xfwd(0,0) = x0 & xfwd(1,0) = y0
 xbck(0,0) = x0 & xbck(1,0) = y0

; ---------------------------------------------
;  Check keywords: step, method and tolerance
; ---------------------------------------------

 IF (NOT KEYWORD_SET(method))  THEN method = "RK2"
 IF (NOT KEYWORD_SET(step))    THEN step = min([ (xend-xbeg)/nx, (yend-ybeg)/ny ])
 IF (NOT KEYWORD_SET(tol))     THEN tol = 1.e-6
 IF (NOT KEYWORD_SET(maxstep)) THEN maxstep = 100*step
 IF (NOT KEYWORD_SET(minstep)) THEN minstep = 0.05*step

; ---------------------------------
;  tol should scale to domain size
; ---------------------------------

 tol = tol*max([Lx, Ly])

; --------------------------------------------------------------------
;  set the coefficients for adaptive step size integration:
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

       k1  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk0)
       xk1 = xk0 + 0.5*dh*k1

       k2 = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xfwd[*,k] = xk0 + dh*k2
     ENDIF

     ; ----------------------------------------------------------
     ;   Explicit Runge-Kutta method with 4th order accuracy.
     ;   Fixed step size. Requires 4 function evaluations.
     ; ----------------------------------------------------------
     IF (method EQ "RK4") THEN BEGIN  
       k = k+1

       k1  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk0)
       xk1 = xk0 + 0.5*dh*k1

       k2  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk1 = xk0 + 0.5*dh*k2

       k3 = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk1 = xk0 + dh*k3

       k4 = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xfwd[*,k] = xk0 + dh*(k1 + 2.0*(k2 + k3) + k4)/6.0
     ENDIF

     ; ---------------------------------------------------------------
     ;  Explicit Runge-Kutta-Fehlberg pair (2,3) with adaptive 
     ;  step size. It is also known as Bogacki-Shampine and provide
     ;  third-order accuracy using a total of 4 function evaluations.
     ; ---------------------------------------------------------------
     IF (method EQ "BS23") THEN BEGIN ; -- use BS23

       k1  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk0)
       xk1 = xk0 + dh*a21*k1

       k2  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk1 = xk0 + dh*(a31*k1 + a32*k2)

       k3  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk1 = xk0 + dh*(a41*k1 + a42*k2 + a43*k3)

       k4  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk3 = xk0 + dh*(b1*k1 + b2*k2 + b3*k3 + b4*k4)

       xk2 = xk0 + dh*(bs1*k1 + bs2*k2 + bs3*k3 + bs4*k4)

       ; ---- compute error ----

       err = max(abs(xk3 - xk2))/tol
       IF ((err LT 1.0) OR (abs(dh)/minstep LE 1.0)) THEN BEGIN  ; -- accept step
         k      = k + 1
         err    = max([err,1.e-12])
         dhnext = 0.9*abs(dh)*err^(-0.3333)
         dhnext = min([dhnext,3.0*abs(dh)])
         dh     = s*dhnext
         xfwd[*,k] = xk3
       ENDIF ELSE BEGIN

         dh = 0.9*s*abs(dh)*err^(-0.5)
         IF (kfail GT MAX_FAIL) THEN BEGIN
           print,"! Too many fails"
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

       k1  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk0)
       xk1 = xk0 + dh*a21*k1

       k2  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk1 = xk0 + dh*(a31*k1 + a32*k2)

       k3  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk1 = xk0 + dh*(a41*k1 + a42*k2 + a43*k3)

       k4  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk1 = xk0 + dh*(a51*k1 + a52*k2 + a53*k3 + a54*k4)

       k5  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk1 = xk0 + dh*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5)

       k6  = FIELD_INTERP(U, V, 0.0, x, y, 0.0, xk1)
       xk5 = xk0 + dh*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6)

       xk4 = xk0 + dh*(bs1*k1 + bs2*k2 + bs3*k3 + bs4*k4 + bs5*k5 + bs6*k6)

       ; ---- compute error ----

       err = max(abs(xk5 - xk4))/tol
       IF ((err LT 1.0) OR (abs(dh)/minstep LE 1.0)) THEN BEGIN  ; -- accept step
         k      = k + 1
         err    = max([err,1.e-12])
         dhnext = 0.9*abs(dh)*err^(-0.2)
         dhnext = min([dhnext,3.0*abs(dh)])
         dh     = s*dhnext
         xfwd[*,k] = xk5
       ENDIF ELSE BEGIN

         dh = 0.9*s*abs(dh)*err^(-0.25)
         IF (kfail GT MAX_FAIL) THEN BEGIN
           print,"! Too many fails"
           RETURN
         ENDIF
         CONTINUE
        ENDELSE
      ENDIF ; method CK45

     ; ----------------------------------------------------
     ;  Check to see whether we're still inside the domain
     ; ----------------------------------------------------

     inside_domain = (xfwd(0,k) GT xbeg) AND (xfwd(0,k) LT xend) AND $
                     (xfwd(1,k) GT ybeg) AND (xfwd(1,k) LT yend)

   ENDWHILE
   IF (S EQ -1) THEN BEGIN
     xbck=xfwd
     k_bck = k
   ENDIF
 ENDFOR
 k_fwd = k
 if (k_fwd GE (MAX_STEPS-2)) THEN PRINT,"! Max number of iteration exceeded"
 if (k_bck GE (MAX_STEPS-2)) THEN PRINT,"! Max number of iteration exceeded"

 print,"Method: ",method,$      ; "; tol = ",strcompress(string(tol,format='(e10.2)')),$
       "; Fwd steps: ",strcompress(string(k_fwd,format='(i5)')),$
       "; bck steps: ",strcompress(string(k_bck,format='(i5)'))

; --------------------------------------------
;         return arrays
; --------------------------------------------

  xfield = [reverse(REFORM(xbck(0,0:k_bck))),REFORM(xfwd(0,0:k_fwd))]
  yfield = [reverse(REFORM(xbck(1,0:k_bck))),REFORM(xfwd(1,0:k_fwd))]

END

