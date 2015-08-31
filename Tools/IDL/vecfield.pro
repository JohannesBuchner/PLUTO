;+
;
;  NAME:      VECFIELD
;
;  AUTHOR:    A. Mignone (mignone@ph.unito.it)
;
;  PURPOSE:   Produces a two-dimensional velocity field plot. 
;
;  SYNTAX:    VECFIELD, U, V, X, Y[,SPACING=value][,/POLAR][,EXTRA]
;
;  KEYWORDS:
;
;       U = the X component of the two-dimensional field. 
;           u must be a two-dimensional array.
;
;       V = the Y component of the two-dimensional field. 
;           v must be a two-dimensional array.
;
;       X = the x-coordinate
;
;       Y = the y-coordinate
;
;       SPACING = the spacing between one arrow and the next
;                 Higher values will produce a more "rarefied" plot
;
;       POLAR = set this keyword if the input data are in polar coordinates
;               (i.e. U = vr, V = vphi, x = r, y = phi) and the vector plot
;               should be output in cartesian coordinates.
;               This is useful when combined with DISPLAY, e.g.,
;
;               IDL> DISPLAY,x1=x1,x2=x2,/polar,rho
;               IDL> VECFIELD,v1,v2,x1,x2,/POLAR,/OVERPLOT
;               
;               Note: for spherical input data (u,v)=(vr,vth) simply call 
;                     VECFIELD with x2 -> 0.5*!PI-x2 and v -> -v.
;
;       EXTRA = all the keyword accepted by velovect, such as
;
;               /overplot, length, etc...
;
;
;  LAST MODIFIED:    Aug 27, 2012
;
;-

PRO VECFIELD, uin, vin, xin, yin, _EXTRA=extra, spacing=spacing,$
              polar=polar

 u = reform(uin)
 v = reform(vin)
 
 x = xin
 y = yin

 IF (KEYWORD_SET(POLAR)) THEN BEGIN

 ; -- transform to cartesian components -- 

   r   = x
   phi = y
   nr   = SIZE(r,   /N_ELEM)
   nphi = SIZE(phi, /N_ELEM)

   phi2D = replicate(1,nr)#phi
   ux  = u*cos(phi2D) - v*sin(phi2D);
   uy  = u*sin(phi2D) + v*cos(phi2D);

   u = ux
   v = uy

 ; -- map to cartesian coordinates -- 

   x = r & y = phi
   POLAR, u, x, y,missing=1.e-12
   x = r & y = phi
   POLAR, v, x, y,missing=1.e-12
  
 ENDIF

 su = size(u)

 nx = su(1)
 ny = su(2)

 IF (NOT KEYWORD_SET(spacing)) THEN spacing = 10

 nxr = nx/spacing
 nyr = ny/spacing

 ur = congrid(u, nxr, nyr)
 vr = congrid(v, nxr, nyr)
 xr = congrid(x, nxr)
 yr = congrid(y, nyr)

 VELOVECT, ur, vr, xr, yr, $ 
           xstyle=1,ystyle=1, _extra=EXTRA
END

