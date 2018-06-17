;+
; NAME:       PARTICLES_SLICE
;
; AUTHOR:     DIPANJAN MUKHEJRE
;
; PURPoSE:    Returns a structure with particles located in a plane x (or y,z) = xs.
;	      Particles within xs-dx/s and xs+dx/2 are chosen.
;	      Useful for making 2D plots of a slice through a 3D data.
;
; SYNTAX:     structure = PARTICLES_SLICE(parts, xar=x, dx=dx, xs=xs, /XY [or ,/YZ, or /ZX])
;
; ARGUMENTS:  
;
;   parts   = Input structure of particles. 
;
;   xar     = Array of the coordinate axis where the planar slice xar=xs is defined.
;         
;   xs      = The coordinate of xar at which the planar slice is defined.
;
;   dx      = Array cell widths along the coordinate axis defined in xar. 
;	      Must be of same size as xar.
;
; KEYWORDS  
;
;  XY or YX       = The slice is in the X-Y plane. The input coordinate axis should be x3.
; 
;  YZ or ZY      = The slice is in the Y-Z plane. The input coordinate axis should be x1.
;
;  ZX or XZ      = The slice is in the Z-X plane. The input coordinate axis should be x2.
;
;  SILENT        = Supress messages.
;
; EXAMPLES:
;
;  part=particles_slice(particles,xar=x1,dx=dx1,/YZ)
; 
;  part=particles_slice(particles,xar=x2,dx=dx2,/ZX)
;
;  part=particles_slice(particles,xar=x3,dx=dx3,/XY)
;
;
;  LAST MODIFIED:
;
;  May 15, D. Mukherjee 
;
;-

FUNCTION PARTICLES_SLICE, parts, xar=xar, dx=dx, xs=xs, XY=XY, YX=YX, YZ=YZ, ZY=ZY, ZX=ZX, XZ=XZ, SILENT=SILENT
  IF (NOT KEYWORD_SET(xs)) THEN BEGIN
     PRINT,'> Particle slice coordinate xs not given. Assuming 0.'
	xs = 0.
  ENDIF

  IF (n_elements(xar) NE n_elements(dx)) THEN BEGIN
     PRINT, '! xar and dx arrays must have same length. Abort!'
     STOP
  ENDIF
  ns = VALUE_LOCATE(xar,xs) ;---locate index of xs in array

  IF (KEYWORD_SET(XY) OR KEYWORD_SET(YX)) THEN BEGIN
     IF (NOT KEYWORD_SET(SILENT)) THEN PRINT,'> XY Slice. Input arrays to be x3, dx3.'
     indx=WHERE((parts.x3 le xar[ns]+dx[ns]/2.) and (parts.x3 ge xar[ns]-dx[ns]/2.))
  ENDIF

  IF (KEYWORD_SET(YZ) OR KEYWORD_SET(ZY)) THEN BEGIN
     IF (NOT KEYWORD_SET(SILENT)) THEN PRINT,'> YZ Slice. Input arrays to be x1, dx1.'
     indx=WHERE((parts.x1 le xar[ns]+dx[ns]/2.) and (parts.x1 ge xar[ns]-dx[ns]/2.))
  ENDIF

  IF (KEYWORD_SET(XZ) OR KEYWORD_SET(ZX)) THEN BEGIN
     IF (NOT KEYWORD_SET(SILENT)) THEN PRINT,'> XZ Slice. Input arrays to be x2, dx2.'
     indx=WHERE((parts.x2 le xar[ns]+dx[ns]/2.) and (parts.x2 ge xar[ns]-dx[ns]/2.))
  ENDIF

RETURN, parts[indx]

END
